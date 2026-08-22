# ==============================================================================
# 4_synthetic_control.R
# Per-project synthetic control with tidysynth, parallelized one country per task.
#
# For each country (selected by array task_id), loop over its protected areas.
# For each PA, treat the PA *as a whole* (area-weighted mean of the outcome
# across all of the PA's grid cells) as the single treated unit, and use the
# country's never-protected grid cells as the donor pool. Every other project's
# cells are excluded (they are neither treated nor valid donors). Standard
# errors are obtained within tidysynth via placebo (in-space) inference.
#
# Outputs (per country, consolidated downstream by 4a_consolidate_synth.R):
#   1. effects_<country>.csv  - one row per project:
#        wdpa_pid, country_iso3, outcome, desig_year, treatment_effect (post ATT),
#        placebo_se, fishers_exact_pvalue, z_score, mspe_ratio,
#        pre_rmspe, post_rmspe, project_area_km2, n_treated_cells, n_donors
#   2. weights_<country>.csv  - one row per (project, donor control cell):
#        wdpa_pid, country_iso3, control_grid_id, weight, weight_rank, in_top10
#
# Usage: Rscript 4_synthetic_control.R <task_id>
#
# Donor pool (per PA): the K never-protected cells most similar to the PA on a
# set of baseline observables (travel times, population access/density, terrain,
# and maize/oil-palm/wheat/soybean observed yields), by standardized
# Euclidean distance. The synthetic control is then fit on those same
# observables plus 5 evenly-spaced pre-period forest-cover levels.
#
# Env vars (optional):
#   SYNTH_OUTCOME   outcome column (default "forest_cover")
#   SYNTH_K_DONORS  number of nearest donors per PA (default 100; fewer avail -> take all)
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(duckdb)
  library(DBI)
})

# ------------------------------------------------------------------------------
# Setup paths
# ------------------------------------------------------------------------------

if (Sys.getenv("SLURM_SUBMIT_DIR") != "") {
  project_root <- Sys.getenv("SLURM_SUBMIT_DIR")
} else {
  project_root <- here::here()
}

data_dir         <- file.path(project_root, "Data/build/final")
synth_dir        <- file.path(project_root, "Data/build/intermediate/synth")
effects_dir      <- file.path(synth_dir, "effects")    # per-country ATT tables
dynamic_dir      <- file.path(synth_dir, "dynamic")    # per-country dynamic effects
stacked_dir      <- file.path(synth_dir, "stacked")    # per-country stacked panels
fig_dir          <- file.path(project_root, "output/figures/analysis/synth_control")
panel_path       <- file.path(data_dir, "TMF_5km_panel.parquet")

for (d_ in c(effects_dir, dynamic_dir, stacked_dir, fig_dir)) {
  if (!dir.exists(d_)) dir.create(d_, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

OUTCOME  <- Sys.getenv("SYNTH_OUTCOME", unset = "forest_cover")
K_DONORS <- suppressWarnings(as.integer(Sys.getenv("SYNTH_K_DONORS", unset = "100")))
if (is.na(K_DONORS) || K_DONORS < 2) K_DONORS <- 100L

# Eligibility thresholds for running SC on a project
MIN_PRE_YEARS     <- 5   # need >=5 pre periods (also gives 5 evenly-spaced lags)
MIN_POST_YEARS    <- 3   # need enough post periods to define an ATT
MIN_TREATED_CELLS <- 1   # PA taken as a whole; >=1 cell
N_PRE_LAGS        <- 5   # number of evenly-spaced pre-period forest-cover lags

# Observable covariates used BOTH to pick the K nearest donors (stage 1) and as
# SC predictors (stage 2). All are time-invariant / baseline. A cell must be
# non-missing on every one of these to be matchable.
COV_COLS <- c("travel_time_cities", "travel_time_ports", "pop_access_50km",
              "pop_density_1990", "elevation_m", "slope_degrees",
              "terrain_ruggedness",
              # EarthStat observed yields. Coffee has no counterpart in the
              # extracted EarthStat set and is dropped; rice / sugarcane are
              # available if this pool should be widened.
              "yield_maize", "yield_oilpalm", "yield_wheat", "yield_soybean")

# ------------------------------------------------------------------------------
# Output schema + writers. Every country that runs writes its two files, even
# when there is nothing to estimate (empty file with headers), so a completed
# task always leaves a trace on disk.
#   effects_<country>.csv  - one row per PA: the post-period ATT + diagnostics
#   dynamic_<country>.csv  - one row per (PA, year): per-period effect tau_t
# ------------------------------------------------------------------------------

EFFECTS_COLS <- c("wdpa_pid", "country_iso3", "outcome", "desig_year",
                  "treatment_effect", "placebo_se", "fishers_exact_pvalue",
                  "z_score", "mspe_ratio", "pre_rmspe", "post_rmspe",
                  "project_area_km2", "n_treated_cells", "n_donors")
DYNAMIC_COLS <- c("wdpa_pid", "country_iso3", "outcome", "desig_year",
                  "year", "event_time", "real_y", "synth_y", "tau",
                  "placebo_sd")
# Stacked treated-and-control panel for a feols event study (one block per PA):
#   unit_id      - PA_<pid> for the treated aggregate, grid_id for each donor
#   wdpa_pid     - the PA this block belongs to (cohort / treatment ID)
#   treated      - 1 for the treated unit, 0 for donors
#   sc_weight    - 1 for the treated unit, synthetic-control weight for donors
#   reg_weight   - area_PA (treated) or area_PA * sc_weight (donor)
STACKED_COLS <- c("wdpa_pid", "country_iso3", "unit_id", "treated",
                  "year", "event_time", "Y", "sc_weight", "reg_weight",
                  "lon", "lat")

empty_effects <- function()
  setNames(data.table(matrix(nrow = 0, ncol = length(EFFECTS_COLS))), EFFECTS_COLS)
empty_dynamic <- function()
  setNames(data.table(matrix(nrow = 0, ncol = length(DYNAMIC_COLS))), DYNAMIC_COLS)
empty_stacked <- function()
  setNames(data.table(matrix(nrow = 0, ncol = length(STACKED_COLS))), STACKED_COLS)

write_country_outputs <- function(effects_dt, dynamic_dt, stacked_dt) {
  eff_file <- file.path(effects_dir, sprintf("effects_%s.csv", country))
  dyn_file <- file.path(dynamic_dir, sprintf("dynamic_%s.csv", country))
  stk_file <- file.path(stacked_dir, sprintf("stacked_%s.csv", country))
  fwrite(effects_dt, eff_file)
  fwrite(dynamic_dt, dyn_file)
  fwrite(stacked_dt, stk_file)
  cat("Effects ->", eff_file, "(", nrow(effects_dt), "rows )\n")
  cat("Dynamic ->", dyn_file, "(", nrow(dynamic_dt), "rows )\n")
  cat("Stacked ->", stk_file, "(", nrow(stacked_dt), "rows )\n")
}

# ------------------------------------------------------------------------------
# Task ID -> country
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Task ID required. Usage: Rscript 4_synthetic_control.R <task_id>")
}
task_id <- as.integer(args[1])

cat("==============================================\n")
cat("4_synthetic_control.R - Per-project synthetic control\n")
cat("Task ID:", task_id, "| Outcome:", OUTCOME,
    "| K donors:", K_DONORS, "\n")
cat("==============================================\n\n")

con <- dbConnect(duckdb())

# Countries that have at least one treated cell (protected, desig >= 1991)
country_query <- sprintf("
  SELECT DISTINCT country_iso3
  FROM read_parquet('%s')
  WHERE year = 1990
    AND is_protected = true
    AND desig_year >= 1991
  ORDER BY country_iso3
", panel_path)
countries <- dbGetQuery(con, country_query)$country_iso3
dbDisconnect(con, shutdown = TRUE)

cat("Found", length(countries), "countries with treated cells\n")
cat("Country -> task_id mapping (use as array index):\n")
cat(paste(sprintf("  [%3d] %s", seq_along(countries), countries),
          collapse = "\n"), "\n\n")

if (task_id > length(countries)) {
  cat("Task ID", task_id, "exceeds country count", length(countries),
      "- nothing to do.\n")
  quit(status = 0)
}

country <- countries[task_id]
cat("Processing country:", country, "(task", task_id, "of",
    length(countries), ")\n\n")
# NB: the per-country figure folder (country_fig_dir) is created later, only
# after the >=30 donor check passes, so skipped countries leave no empty folder.

# ------------------------------------------------------------------------------
# Load the full panel for this country (all cells, all years)
# ------------------------------------------------------------------------------

cat("Loading panel for", country, "...\n")

con <- dbConnect(duckdb())
query <- sprintf("
  SELECT
    grid_id, area_km2, wdpa_pid, country_iso3,
    centroid_lon, centroid_lat,
    year, desig_year, is_protected,
    is_interior, is_frontier,
    %s AS Y,
    %s
  FROM read_parquet('%s')
  WHERE country_iso3 = '%s'
", OUTCOME, paste(COV_COLS, collapse = ", "), panel_path, country)

dt <- as.data.table(dbGetQuery(con, query))
dbDisconnect(con, shutdown = TRUE)

cat("Loaded", format(nrow(dt), big.mark = ","), "cell-years |",
    format(uniqueN(dt$grid_id), big.mark = ","), "cells\n")

# Drop rows with a missing year (poisons the balanced-panel year count)
n_before <- nrow(dt)
dt <- dt[!is.na(year)]
cat("Dropped", format(n_before - nrow(dt), big.mark = ","),
    "rows with NA year\n")

gc()

year_min      <- min(dt$year)
year_max      <- max(dt$year)
n_years_total <- uniqueN(dt$year)
cat("Year range:", year_min, "-", year_max, "(", n_years_total, "years)\n")

# Country-wide cell coordinates (one row per cell) for full-extent weight maps.
# Captured before any cleaning subset so maps show the whole country as context.
coords <- unique(dt[, .(grid_id = as.character(grid_id),
                        centroid_lon, centroid_lat)])

# Load tidysynth after the DuckDB query is finished (mirrors matching scripts)
suppressPackageStartupMessages({
  library(tidysynth)
  library(ggplot2)
})

# ==============================================================================
# USER CLEANING RULES (applied FIRST, before any other eligibility checks)
#   (1) Keep only PAs and control cells with >= 1 cell EVER classified as
#       forest interior or frontier (is_interior | is_frontier in any year).
#   (2) Keep only PAs whose area-weighted average forest cover is > 0 in at
#       least one period.
# ==============================================================================

cat("\n=== Applying user cleaning rules ===\n")

n_cells_start    <- uniqueN(dt$grid_id)
n_prot_start     <- uniqueN(dt[is_protected == TRUE, grid_id])
n_ctrl_start     <- uniqueN(dt[is_protected == FALSE, grid_id])

# --- Rule 1: cell is EVER interior or frontier --------------------------------
ever_if <- dt[, .(is_if = any((is_interior | is_frontier) %in% TRUE)), by = grid_id]
keep_if_ids <- ever_if[is_if == TRUE, grid_id]

# Control cells: keep the cell itself if it is ever interior/frontier
ctrl_keep_ids <- intersect(dt[is_protected == FALSE, unique(grid_id)], keep_if_ids)

# PAs: keep the whole PA if >= 1 of its cells is ever interior/frontier
pa_if <- dt[is_protected == TRUE & !is.na(wdpa_pid),
            .(has_if = any(grid_id %in% keep_if_ids)), by = wdpa_pid]
pa_keep_rule1 <- pa_if[has_if == TRUE, wdpa_pid]

cat(sprintf("Rule 1 (interior/frontier): controls %d -> %d | PAs %d -> %d\n",
            n_ctrl_start, length(ctrl_keep_ids),
            uniqueN(dt[is_protected == TRUE & !is.na(wdpa_pid), wdpa_pid]),
            length(pa_keep_rule1)))

# --- Rule 2: PA average forest cover > 0 in at least one period ----------------
pa_fc <- dt[is_protected == TRUE & wdpa_pid %in% pa_keep_rule1,
            .(fc = weighted.mean(Y, area_km2, na.rm = TRUE)),
            by = .(wdpa_pid, year)][
            , .(max_fc = max(fc, na.rm = TRUE)), by = wdpa_pid]
pa_keep_rule2 <- pa_fc[is.finite(max_fc) & max_fc > 0, wdpa_pid]

cat(sprintf("Rule 2 (avg forest cover > 0): PAs %d -> %d\n",
            length(pa_keep_rule1), length(pa_keep_rule2)))

pa_keep_ids <- pa_keep_rule2

# Restrict the working panel: kept control cells + all cells of kept PAs
dt <- dt[(is_protected == FALSE & grid_id %in% ctrl_keep_ids) |
         (is_protected == TRUE  & wdpa_pid %in% pa_keep_ids)]

cat(sprintf("After user rules: %d cells (was %d) | %d PAs | %d control cells\n",
            uniqueN(dt$grid_id), n_cells_start,
            length(pa_keep_ids), uniqueN(dt[is_protected == FALSE, grid_id])))

if (length(pa_keep_ids) == 0) {
  cat("No PAs survive user cleaning rules in", country, "- skipping country.\n")
  write_country_outputs(empty_effects(), empty_dynamic(), empty_stacked())
  quit(status = 0)
}

# ------------------------------------------------------------------------------
# Build the donor pool once: never-protected cells with a complete record
# ------------------------------------------------------------------------------

cat("\nBuilding donor pool (never-protected cells)...\n")

donors_raw <- dt[is_protected == FALSE]

# --- Covariate rule: require non-missing values on ALL matching covariates ----
# A cell can only be distance-ranked if it has every matching covariate. Report
# how many control AND protected cells are dropped.
cov_bad <- function(d) Reduce(`|`, lapply(COV_COLS, function(v) is.na(d[[v]])))

ctrl_cov_bad_ids <- unique(donors_raw[cov_bad(donors_raw), grid_id])
prot_cells       <- dt[is_protected == TRUE]
prot_cov_bad_ids <- unique(prot_cells[cov_bad(prot_cells), grid_id])

cat(sprintf(paste0("Covariate rule (all %d matching covariates non-missing): ",
                   "dropped %d control cells, %d protected cells\n"),
            length(COV_COLS), length(ctrl_cov_bad_ids), length(prot_cov_bad_ids)))

donors_raw <- donors_raw[!grid_id %in% ctrl_cov_bad_ids]
dt <- dt[!grid_id %in% prot_cov_bad_ids]   # keep treated aggregates covariate-clean

# --- Balanced-panel rule (OUTCOME only) ---------------------------------------
# tidysynth needs every unit observed in every year with a non-NA outcome, else
# the units x years matrix it builds has holes and synth() fails. Enforced on
# the outcome only (covariates no longer required to be complete across years).
donor_ok <- donors_raw[, .(nyr = uniqueN(year), n_na_y = sum(is.na(Y))),
                       by = grid_id][nyr == n_years_total & n_na_y == 0, grid_id]
donors_dt <- donors_raw[grid_id %in% donor_ok]
n_donor_cells <- uniqueN(donors_dt$grid_id)
cat("Eligible candidate donor cells:", format(n_donor_cells, big.mark = ","),
    "(stage 1 picks up to", K_DONORS, "nearest per PA)\n")

if (n_donor_cells < 2) {
  cat("Fewer than 2 eligible donors in", country, "- skipping country.\n")
  write_country_outputs(empty_effects(), empty_dynamic(), empty_stacked())
  quit(status = 0)
}

# Country will run -> create the per-country figure folder
country_fig_dir <- file.path(fig_dir, country)
if (!dir.exists(country_fig_dir)) dir.create(country_fig_dir, recursive = TRUE)

# Candidate donor long table (subset to the K nearest per PA in run_one_pa)
donor_long <- donors_dt[, c("grid_id", "year", "Y", COV_COLS), with = FALSE]
setnames(donor_long, "grid_id", "unit_name")
donor_long[, unit_name := as.character(unit_name)]

# --- Stage-1 matching setup: standardized donor covariate matrix --------------
# Covariates are time-invariant, so one row per donor cell. Standardize on the
# candidate-pool mean/sd (shared by treated + donors). Drop covariates with no
# spread across donors (they cannot discriminate).
donor_cov <- unique(donors_dt[, c("grid_id", COV_COLS), with = FALSE])
donor_cov[, grid_id := as.character(grid_id)]

cov_center <- vapply(COV_COLS, function(v) mean(donor_cov[[v]], na.rm = TRUE), numeric(1))
cov_scale  <- vapply(COV_COLS, function(v) sd(donor_cov[[v]],   na.rm = TRUE), numeric(1))
match_cols <- COV_COLS[is.finite(cov_scale) & cov_scale > 0]
if (length(match_cols) < length(COV_COLS)) {
  cat("Note: dropping", length(COV_COLS) - length(match_cols),
      "zero-variance covariate(s) from matching:",
      paste(setdiff(COV_COLS, match_cols), collapse = ", "), "\n")
}

donor_z <- as.matrix(donor_cov[, ..match_cols])
donor_z <- sweep(donor_z, 2, cov_center[match_cols], "-")
donor_z <- sweep(donor_z, 2, cov_scale[match_cols],  "/")
rownames(donor_z) <- donor_cov$grid_id

# ------------------------------------------------------------------------------
# Enumerate projects (protected areas) to estimate
# ------------------------------------------------------------------------------

pa_meta <- dt[is_protected == TRUE & !is.na(wdpa_pid) &
              !is.na(desig_year) & desig_year >= 1991,
              .(desig_year = min(desig_year)), by = wdpa_pid][order(wdpa_pid)]

cat("\nProjects (PAs) to consider:", nrow(pa_meta), "\n")

# Accumulators
effects_list <- list()
dynamic_list <- list()
stacked_list <- list()
map_data     <- list()   # per estimated PA: weights + treated-cell ids, for maps

# ------------------------------------------------------------------------------
# Write one observed-vs-synthetic plot for a single project
# ------------------------------------------------------------------------------

plot_one_pa <- function(path, pid, d) {
  pl <- melt(
    path,
    id.vars = c("wdpa_pid", "desig_year", "time_unit"),
    measure.vars = c("real_y", "synth_y"),
    variable.name = "series", value.name = "forest_frac"
  )
  pl[, series := factor(
    series, levels = c("real_y", "synth_y"),
    labels = c("Observed project", "Synthetic control")
  )]

  p <- ggplot(pl, aes(x = time_unit, y = forest_frac,
                      color = series, linetype = series)) +
    geom_vline(xintercept = d, linetype = "dotted",
               color = "gray60", linewidth = 0.5) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = c("Observed project" = "steelblue",
                                  "Synthetic control" = "coral")) +
    scale_linetype_manual(values = c("Observed project" = "solid",
                                     "Synthetic control" = "dashed")) +
    labs(
      x = "Year", y = "Forest fraction", color = NULL, linetype = NULL,
      title = sprintf("Synthetic Control: PA %s (%s)", pid, country),
      subtitle = sprintf("Designation year %d (dotted line)", d)
    ) +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.text = element_text(color = "black", size = 9),
      axis.title = element_text(color = "black", size = 11),
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0),
      legend.position = "bottom",
      plot.margin = margin(15, 15, 15, 15)
    )

  fig_path <- file.path(country_fig_dir, sprintf("synth_PA%s.png", pid))
  ggsave(fig_path, p, width = 7, height = 5, dpi = 300, bg = "white")
  fig_path
}

# ------------------------------------------------------------------------------
# Map a PA's donor weights over the full country extent: all country cells in
# grey, donor cells colored by SC weight, the PA's own cells marked in red.
# Styling follows visualize_wdpa_tmf.R (viridis fill, theme_void).
# ------------------------------------------------------------------------------

plot_weight_map <- function(pid, d, weights, pa_cell_ids) {
  donor_xy <- merge(weights[, .(grid_id = control_grid_id, weight)],
                    coords, by = "grid_id")
  pa_xy    <- coords[grid_id %in% as.character(pa_cell_ids)]

  p <- ggplot() +
    # Whole-country context
    geom_point(data = coords, aes(centroid_lon, centroid_lat),
               color = "grey85", size = 0.3) +
    # Donor cells colored by weight
    geom_point(data = donor_xy,
               aes(centroid_lon, centroid_lat, color = weight), size = 1.6) +
    viridis::scale_color_viridis(name = "Donor\nweight", option = "viridis") +
    # The protected area itself
    geom_point(data = pa_xy, aes(centroid_lon, centroid_lat),
               color = "red", shape = 17, size = 1.4) +
    coord_quickmap() +
    labs(
      title = sprintf("Donor weights: PA %s (%s)", pid, country),
      subtitle = sprintf("Designation %d | red = PA cells, colored = donors", d),
      x = NULL, y = NULL
    ) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0),
      legend.position = "right",
      plot.margin = margin(10, 10, 10, 10)
    )

  map_dir <- file.path(country_fig_dir, "maps")
  if (!dir.exists(map_dir)) dir.create(map_dir, recursive = TRUE)
  fig_path <- file.path(map_dir, sprintf("weightmap_PA%s.png", pid))
  ggsave(fig_path, p, width = 7, height = 6, dpi = 300, bg = "white")
  fig_path
}

# ------------------------------------------------------------------------------
# Per-project synthetic control
# ------------------------------------------------------------------------------

run_one_pa <- function(pid, d) {
  treated_cells <- dt[is_protected == TRUE & wdpa_pid == pid]
  n_tr <- uniqueN(treated_cells$grid_id)
  if (n_tr < MIN_TREATED_CELLS) return(NULL)

  # Eligibility on pre/post window length
  pre_years  <- year_min:(d - 1)
  post_years <- d:year_max
  if (length(pre_years) < MIN_PRE_YEARS)  return(NULL)
  if (length(post_years) < MIN_POST_YEARS) return(NULL)

  # Matching window: the N years immediately before treatment, d-N .. d-1.
  # The synthetic control matches covariates and outcome lags on THESE years.
  match_years <- (d - N_PRE_LAGS):(d - 1)

  # Project area: per-cell area is time-invariant; sum across distinct cells
  proj_area <- treated_cells[, .(a = area_km2[1]), by = grid_id][, sum(a, na.rm = TRUE)]

  # Aggregate the PA into a single time series (area-weighted mean per year)
  treated_long <- treated_cells[
    , c(
        list(Y = if (all(is.na(Y))) NA_real_ else
                 weighted.mean(Y, area_km2, na.rm = TRUE)),
        lapply(COV_COLS, function(v)
                 weighted.mean(get(v), area_km2, na.rm = TRUE))
      ),
    by = year]
  setnames(treated_long, c("year", "Y", COV_COLS))

  # Treated unit must be complete over the full window
  if (uniqueN(treated_long$year) != n_years_total || any(is.na(treated_long$Y)))
    return(NULL)

  treated_name <- paste0("PA_", pid)
  treated_long[, unit_name := treated_name]
  setcolorder(treated_long, c("unit_name", "year", "Y", COV_COLS))

  # --- Stage 1: pick the K nearest donors by standardized covariate distance ---
  t_vec <- unlist(treated_long[1, ..COV_COLS])
  if (length(match_cols) > 0) {
    t_z <- (t_vec[match_cols] - cov_center[match_cols]) / cov_scale[match_cols]
    d2  <- rowSums(sweep(donor_z, 2, t_z, "-")^2)
    ord <- order(d2)
  } else {
    ord <- seq_len(nrow(donor_z))
  }
  k_sel   <- min(K_DONORS, length(ord))
  sel_ids <- rownames(donor_z)[ord[seq_len(k_sel)]]

  sel_donor_long <- donor_long[unit_name %in% sel_ids]
  panel <- rbindlist(list(treated_long, sel_donor_long), use.names = TRUE)

  # --- Stage 2: SC predictors = matching covariates + N pre-period lags --------
  # Keep only covariates that vary across the chosen units (constant predictors
  # make synth() fail with "no variation across control units").
  sel_cov   <- donor_cov[grid_id %in% sel_ids, ..COV_COLS]
  unit_vals <- rbind(as.data.table(as.list(t_vec))[, ..COV_COLS], sel_cov)
  use_cov   <- COV_COLS[vapply(COV_COLS, function(v) {
    x <- unit_vals[[v]]; length(unique(x[!is.na(x)])) > 1
  }, logical(1))]
  # Build covariate-mean predictor expressions for tidy-eval injection. Each is
  # mean(<col>, na.rm = TRUE); names become m_<col>. generate_predictor() is a
  # data-masking verb, so these are spliced in with !!! (not passed as values).
  cov_exprs <- setNames(
    lapply(use_cov, function(v) rlang::expr(mean(!!rlang::sym(v), na.rm = TRUE))),
    paste0("m_", use_cov)
  )

  # Outcome lags: the individual years d-N .. d-1 (one predictor per year).
  lag_years <- match_years

  sc <- panel %>%
    synthetic_control(
      outcome = Y, unit = unit_name, time = year,
      i_unit = treated_name, i_time = d,
      generate_placebos = TRUE
    )
  # Covariate-mean predictors over the matching window d-N..d-1 (spliced via !!!)
  sc <- generate_predictor(sc, time_window = match_years, !!!cov_exprs)
  # Pre-treatment outcome lags, one per year d-N..d-1 (fc1..fcN)
  for (li in seq_along(lag_years)) {
    nm <- paste0("fc", li)
    sc <- generate_predictor(sc, time_window = lag_years[li], !!nm := Y)
  }
  sc <- sc %>%
    generate_weights(optimization_window = match_years) %>%
    generate_control()

  # --- treatment effect (post-period ATT) -----------------------------------
  ts <- as.data.table(grab_synthetic_control(sc, placebo = FALSE))
  # columns: time_unit, real_y, synth_y
  att      <- ts[time_unit >= d, mean(real_y - synth_y, na.rm = TRUE)]
  pre_rmspe  <- sqrt(ts[time_unit < d,  mean((real_y - synth_y)^2, na.rm = TRUE)])
  post_rmspe <- sqrt(ts[time_unit >= d, mean((real_y - synth_y)^2, na.rm = TRUE)])

  # --- placebo (in-space) inference ------------------------------------------
  plc <- as.data.table(grab_synthetic_control(sc, placebo = TRUE))
  # columns: .id, .placebo, time_unit, real_y, synth_y  (.placebo == 1 for donors)
  if (".id" %in% names(plc))      setnames(plc, ".id", "plc_id")
  if (".placebo" %in% names(plc)) setnames(plc, ".placebo", "is_placebo")
  plc_att <- plc[is_placebo == 1 & time_unit >= d,
                 .(att = mean(real_y - synth_y, na.rm = TRUE)), by = plc_id]
  placebo_se <- if (nrow(plc_att) > 1) sd(plc_att$att, na.rm = TRUE) else NA_real_

  # Per-event-time placebo SE: spread of the placebo gap across donor units at
  # each year. This per-period SE is what gets area-weighted-averaged when the
  # PAs are pooled into an event study (4b).
  plc_sd <- plc[is_placebo == 1,
                .(placebo_sd = sd(real_y - synth_y, na.rm = TRUE)), by = time_unit]
  setnames(plc_sd, "time_unit", "year")

  # tidysynth significance table (Fisher exact p-value, z-score, MSPE ratio)
  sig <- as.data.table(grab_significance(sc))
  sig_tr <- sig[unit_name == treated_name]
  fisher_p   <- if ("fishers_exact_pvalue" %in% names(sig_tr) && nrow(sig_tr))
                  sig_tr$fishers_exact_pvalue[1] else NA_real_
  zscore     <- if ("z_score" %in% names(sig_tr) && nrow(sig_tr))
                  sig_tr$z_score[1] else NA_real_
  mspe_ratio <- if ("mspe_ratio" %in% names(sig_tr) && nrow(sig_tr))
                  sig_tr$mspe_ratio[1] else NA_real_

  # --- donor weights ---------------------------------------------------------
  w <- as.data.table(grab_unit_weights(sc))  # columns: unit, weight
  setnames(w, c("unit", "weight"), c("control_grid_id", "weight"))
  setorder(w, -weight)
  w[, weight_rank := .I]
  w[, in_top10 := weight_rank <= 10L]
  w[, `:=`(wdpa_pid = pid, country_iso3 = country)]
  setcolorder(w, c("wdpa_pid", "country_iso3", "control_grid_id",
                   "weight", "weight_rank", "in_top10"))

  eff <- data.table(
    wdpa_pid = pid,
    country_iso3 = country,
    outcome = OUTCOME,
    desig_year = d,
    treatment_effect = att,
    placebo_se = placebo_se,
    fishers_exact_pvalue = fisher_p,
    z_score = zscore,
    mspe_ratio = mspe_ratio,
    pre_rmspe = pre_rmspe,
    post_rmspe = post_rmspe,
    project_area_km2 = proj_area,
    n_treated_cells = n_tr,
    n_donors = k_sel
  )

  # Observed vs synthetic path for the per-project plot
  path <- ts[, .(time_unit, real_y, synth_y)]
  path[, `:=`(wdpa_pid = pid, desig_year = d)]

  # Dynamic effects: one row per year, tau_t = real_y - synth_y, with event time
  # and the per-event-time placebo SE.
  dyn <- ts[, .(year = time_unit, real_y, synth_y)]
  dyn[, `:=`(
    wdpa_pid = pid, country_iso3 = country, outcome = OUTCOME,
    desig_year = d, event_time = year - d, tau = real_y - synth_y
  )]
  dyn <- merge(dyn, plc_sd, by = "year", all.x = TRUE)
  setcolorder(dyn, DYNAMIC_COLS)

  # --- Stacked treated+control panel for a feols event study ------------------
  # Treated block: the PA aggregate, weight = project area, sc_weight = 1.
  # Coordinates: area-weighted centroid of the PA's cells (one point per PA).
  pa_lon <- treated_cells[, weighted.mean(centroid_lon, area_km2, na.rm = TRUE)]
  pa_lat <- treated_cells[, weighted.mean(centroid_lat, area_km2, na.rm = TRUE)]
  stk_tr <- treated_long[, .(unit_id = treated_name, year, Y)]
  stk_tr[, `:=`(treated = 1L, sc_weight = 1, reg_weight = proj_area,
                lon = pa_lon, lat = pa_lat)]

  # Control block: each selected donor's own series, weighted by its SC weight.
  # reg_weight = project area * sc_weight so each cohort's controls carry the
  # treated unit's mass distributed across donors. Coordinates: the donor cell's
  # own centroid (coords keyed by grid_id == unit_id).
  wmap <- w[, .(unit_id = as.character(control_grid_id), sc_weight = weight)]
  stk_ct <- merge(
    sel_donor_long[, .(unit_id = unit_name, year, Y)],
    wmap, by = "unit_id", all.x = TRUE
  )
  stk_ct[is.na(sc_weight), sc_weight := 0]
  stk_ct[, `:=`(treated = 0L, reg_weight = proj_area * sc_weight)]
  stk_ct <- merge(
    stk_ct,
    coords[, .(unit_id = grid_id, lon = centroid_lon, lat = centroid_lat)],
    by = "unit_id", all.x = TRUE
  )

  stacked <- rbindlist(list(stk_tr, stk_ct), use.names = TRUE)
  # Event time gated on treatment: controls get -9999 (absorbed into the ref).
  stacked[, `:=`(wdpa_pid = pid, country_iso3 = country,
                 event_time = fifelse(treated == 1L, year - d, -9999L))]
  setcolorder(stacked, STACKED_COLS)

  list(eff = eff, weights = w, path = path, dyn = dyn, stacked = stacked)
}

n_pa <- nrow(pa_meta)
for (i in seq_len(n_pa)) {
  pid <- pa_meta$wdpa_pid[i]
  d   <- pa_meta$desig_year[i]
  cat(sprintf("[%d/%d] PA %s (desig %d)... ", i, n_pa, pid, d))

  res <- tryCatch(
    run_one_pa(pid, d),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
  )

  if (is.null(res)) {
    cat("skipped/failed.\n")
    next
  }

  effects_list[[length(effects_list) + 1]] <- res$eff
  dynamic_list[[length(dynamic_list) + 1]] <- res$dyn
  stacked_list[[length(stacked_list) + 1]] <- res$stacked

  # Stash weights + PA cell ids for the (later) random sample of weight maps
  pa_cell_ids <- dt[is_protected == TRUE & wdpa_pid == pid, unique(grid_id)]
  map_data[[as.character(pid)]] <- list(pid = pid, d = d,
                                        weights = res$weights,
                                        pa_cell_ids = pa_cell_ids)

  # Write this project's observed-vs-synthetic plot immediately (real-time)
  fp <- tryCatch(plot_one_pa(res$path, pid, d),
                 error = function(e) { cat("(plot failed:",
                                           conditionMessage(e), ") "); NA_character_ })

  cat(sprintf("ATT=%.4f se=%.4f p=%.3f | plot: %s\n",
              res$eff$treatment_effect, res$eff$placebo_se,
              res$eff$fishers_exact_pvalue,
              if (is.na(fp)) "FAILED" else basename(fp)))
}

# ------------------------------------------------------------------------------
# Weight maps for up to 5 random estimated PAs (full country extent)
# ------------------------------------------------------------------------------

if (length(map_data) > 0) {
  set.seed(task_id)  # reproducible per task
  pick <- sample(seq_along(map_data), size = min(5L, length(map_data)))
  cat("\nDrawing", length(pick), "weight map(s)...\n")
  for (j in pick) {
    md <- map_data[[j]]
    mp <- tryCatch(
      plot_weight_map(md$pid, md$d, md$weights, md$pa_cell_ids),
      error = function(e) { cat("  PA", md$pid, "map failed:",
                                conditionMessage(e), "\n"); NA_character_ })
    if (!is.na(mp)) cat("  map ->", basename(mp), "\n")
  }
}

# ------------------------------------------------------------------------------
# Write per-country outputs (always, even if empty)
# ------------------------------------------------------------------------------

effects_dt <- if (length(effects_list)) rbindlist(effects_list) else empty_effects()
dynamic_dt <- if (length(dynamic_list)) rbindlist(dynamic_list) else empty_dynamic()
stacked_dt <- if (length(stacked_list)) rbindlist(stacked_list) else empty_stacked()

write_country_outputs(effects_dt, dynamic_dt, stacked_dt)

cat("\n==============================================\n")
cat("Completed:", country, "\n")
cat("Projects estimated:", nrow(effects_dt), "of", n_pa, "considered\n")
cat("Per-project plots ->", country_fig_dir, "\n")
cat("==============================================\n")
