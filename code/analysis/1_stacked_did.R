# ==============================================================================
# 1_stacked_did.R
# ==============================================================================
# Stacked difference-in-differences on PA designation, with clean never-treated
# controls drawn from the same country x biome as each treated cohort.
#
# Outputs (output/analysis/):
#   raw_trends_stacked.png        raw trends, ID x cohort intercept removed at t = -1
#   did_estimates.tex             static + dynamic DiD, two FE structures
#   event_study_stacked.png       dynamic coefficients
#   balance_treated_vs_never.tex  AER booktabs balance table
#   stacked_did_sample.txt        sample construction log (what was kept/dropped)
#
# REPLACES 1_balance_PA.R. See the "DROPPED" note at the bottom of this header.
#
# UNIT OF ANALYSIS: 5km grid cells, not municipalities. The panel has no
# municipal identifier; grid_id is the finest unit available and is what the
# PA boundaries were cut on.
#
# DROPPED relative to 1_balance_PA.R:
#   - Balance table by IUCN category (Table 2)
#   - Designation-year histogram (3-panel)
#   - Designation year by continent (stacked bars)
#   - Designation year by zone interior/frontier/other (stacked bars)
#   - Continent x zone cross-tabulation
#   Table 1 (protected vs unprotected, area-weighted, cross-section) is
#   superseded by the balance table here, which compares treated vs
#   NEVER-treated rather than protected vs unprotected in a single year.
#
# Usage: sbatch code/bash/analysis/1_stacked_did.sh
# ==============================================================================

suppressMessages({
  library(data.table)
  library(duckdb)
  library(DBI)
  library(fixest)
  library(ggplot2)
  library(sf)
  library(arrow)   # reading the Stage 1b cropland cross-section
})

# ------------------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------------------

# Both outcomes are run in ONE pass and written to outcome-suffixed filenames,
# so neither overwrites the other and every run produces the full set.
#   forest_cover     Undisturbed + Degraded + Regrowth  (forest EXTENT)
#   Undisturbed_TMF  undisturbed only                   (forest QUALITY)
# The pair is the substantive comparison: extent holding while undisturbed
# falls is degradation inside a protected boundary.
OUTCOMES       <- trimws(strsplit(Sys.getenv("DID_OUTCOMES",
                    "forest_cover,Undisturbed_TMF"), ",")[[1]])
PRETREND_VAR   <- Sys.getenv("DID_PRETREND_VAR", "Deforested")
EVENT_MIN      <- as.integer(Sys.getenv("DID_EVENT_MIN", "-10"))
EVENT_MAX      <- as.integer(Sys.getenv("DID_EVENT_MAX", "10"))
# Controls are capped per cohort: every never-treated cell in the country-biome
# would make the stack enormous (1.66M never-treated cells x ~28 cohorts x 21
# years). The cap is a random sample WITHIN each country-biome, so the control
# pool stays representative of the stratum. Raise via env var; the number kept
# and dropped is logged per cohort, never silently truncated.
CONTROL_CAP    <- as.integer(Sys.getenv("DID_CONTROL_CAP", "40000"))
MIN_PRE_YEARS  <- 5   # for the pre-trend column; units with fewer are dropped
SEED           <- 42

set.seed(SEED)

project_root <- if (Sys.getenv("SLURM_SUBMIT_DIR") != "") {
  Sys.getenv("SLURM_SUBMIT_DIR")
} else here::here()

panel_path  <- file.path(project_root, "Data/build/final/TMF_5km_panel.parquet")
# Stage 1b cross-section. Defined here rather than by sourcing BUILD_workspace.R,
# which pulls in sf/terra/exactextractr and the whole build-side path set that
# this analysis script does not otherwise need.
cropland_final_file <- file.path(project_root, "Data/build/final/TMF_5km_cropland.parquet")
biome_cache <- file.path(project_root, "Data/build/final/grid_biome.parquet")
eco_path    <- file.path(
  Sys.getenv("FRAG_SPATIAL",
             "/resnick/groups/MishraLab/GlobalForest/data/raw/spatial"),
  "ecoregions", "Ecoregions2017", "Ecoregions2017.shp")

out_dir <- file.path(project_root, "output", "analysis", "stacked_did")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_lines <- c()
say <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n"); log_lines <<- c(log_lines, msg)
}

say("Panel: %s", panel_path)
say("Outcomes: %s | event window [%d, %d] | control cap %s/cohort",
    paste(OUTCOMES, collapse = ", "), EVENT_MIN, EVENT_MAX,
    format(CONTROL_CAP, big.mark = ","))

# Time-invariant characteristics for the balance table.
# Stage 1b moved the crop/yield variables out of the panel into a cell-level
# cross-section (TMF_5km_cropland.parquet), joined on grid_id. The panel now
# carries 57 columns and none of the 172-crop set, so the two groups are read
# from different files and merged below. Asking DuckDB for yield_maize against
# the panel is what killed job 1447010:
#   Binder Error: Referenced column "yield_maize" not found in FROM clause!
BAL_VARS_PANEL <- c(
  "elevation_m", "slope_degrees", "terrain_ruggedness", "frac_slope_gt15",
  "dist_to_city_km", "travel_time_cities", "travel_time_ports",
  "pop_density_1990", "pop_access_50km",
  "aboveground_biomass_carbon_2010", "biomass_access_50km",
  "cattle_density_2010_da",
  "biodiv_combined_thr_sr_2022", "biodiv_rsr_crenvu",
  "area_km2")

# From the cropland cross-section. cropland_frac / cropshare_total summarise
# agricultural pressure better than any single crop now that 172 are available;
# yield_maize is kept as the one crop-productivity proxy.
BAL_VARS_CROP <- c("cropland_frac", "pasture_frac", "cropshare_total", "yield_maize")

BAL_VARS <- c(BAL_VARS_PANEL, BAL_VARS_CROP)

# Variables reported as percentiles of the pooled (treated + control)
# distribution rather than in native units. Range-size rarity is a sum of
# inverse species range areas, so its raw values sit around 1e-6 and carry no
# intuitive meaning; a percentile says directly where a cell sits in the
# rarity distribution. The transform is monotone, so it changes the scale of
# the row, not which side is rarer.
PCTL_VARS <- c("biodiv_rsr_crenvu")

BAL_LABELS <- c(
  elevation_m = "Elevation (m)", slope_degrees = "Slope (degrees)",
  terrain_ruggedness = "Terrain ruggedness", frac_slope_gt15 = "Share slope $>15^\\circ$",
  dist_to_city_km = "Distance to city (km)", travel_time_cities = "Travel time to cities",
  travel_time_ports = "Travel time to ports", pop_density_1990 = "Population density, 1990",
  pop_access_50km = "Population access, 50km",
  aboveground_biomass_carbon_2010 = "Aboveground biomass C, 2010",
  biomass_access_50km = "Biomass access, 50km",
  cropland_frac = "Cropland share of cell", pasture_frac = "Pasture share of cell",
  cropshare_total = "Harvested area / cell area", yield_maize = "Maize yield (t/ha)",
  cattle_density_2010_da = "Cattle density, 2010",
  biodiv_combined_thr_sr_2022 = "Threatened species richness",
  biodiv_rsr_crenvu = "Range-size rarity (CR/EN/VU), pctile",
  area_km2 = "Cell area (km$^2$)")

# ------------------------------------------------------------------------------
# BIOME LOOKUP (built once, cached)
#
# The panel has no biome column. Cell centroids are joined to RESOLVE
# Ecoregions2017 (847 polygons, WGS84) and the result cached, so this cost is
# paid once rather than on every run.
# ------------------------------------------------------------------------------

con <- dbConnect(duckdb(), shared_home = FALSE)
dbExecute(con, "SET memory_limit='48GB'")
dbExecute(con, "SET threads=4")

if (file.exists(biome_cache)) {
  say("Biome lookup: cached (%s)", basename(biome_cache))
  biome_dt <- as.data.table(arrow::read_parquet(biome_cache))
} else {
  say("Biome lookup: building from %s", basename(eco_path))
  stopifnot(file.exists(eco_path))

  cells <- as.data.table(dbGetQuery(con, sprintf(
    "SELECT DISTINCT grid_id, centroid_lon, centroid_lat
     FROM read_parquet('%s') WHERE centroid_lon IS NOT NULL", panel_path)))
  say("  %s cells to assign", format(nrow(cells), big.mark = ","))

  eco <- st_read(eco_path, quiet = TRUE)[, c("BIOME_NAME", "BIOME_NUM")]
  eco <- st_make_valid(eco)

  pts <- st_as_sf(cells, coords = c("centroid_lon", "centroid_lat"), crs = 4326)
  sf_use_s2(FALSE)
  # st_intersects, NOT st_join: a point inside overlapping ecoregion polygons
  # gets one row per match from st_join, so the result is longer than the input
  # and cbind-ing it back against grid_id silently shifts every assignment after
  # the first duplicate. Taking the first match per point keeps 1 row per cell
  # by construction.
  hits <- st_intersects(pts, eco)
  sf_use_s2(TRUE)

  first_hit <- vapply(hits, function(i) if (length(i)) i[[1]] else NA_integer_,
                      integer(1))
  n_multi <- sum(lengths(hits) > 1)
  if (n_multi > 0) {
    say("  %s cells fell in overlapping ecoregions; first match used",
        format(n_multi, big.mark = ","))
  }

  biome_dt <- data.table(grid_id = cells$grid_id,
                         biome = as.character(eco$BIOME_NAME)[first_hit])
  stopifnot(nrow(biome_dt) == nrow(cells))
  n_missing <- sum(is.na(biome_dt$biome))
  say("  unassigned (no ecoregion polygon): %s (%.2f%%)",
      format(n_missing, big.mark = ","), 100 * n_missing / nrow(biome_dt))

  arrow::write_parquet(biome_dt, biome_cache)
  rm(pts, joined, eco, cells); gc()
}

# ------------------------------------------------------------------------------
# CELL-LEVEL ATTRIBUTES (one row per grid_id)
# ------------------------------------------------------------------------------

say("Loading cell attributes...")

attr_cols <- paste(sprintf("any_value(%s) AS %s", BAL_VARS_PANEL, BAL_VARS_PANEL), collapse = ",\n    ")
cell_attr <- as.data.table(dbGetQuery(con, sprintf("
  SELECT grid_id,
         any_value(country_iso3)                       AS country_iso3,
         max(CASE WHEN is_protected THEN 1 ELSE 0 END) AS ever_protected,
         min(desig_year)                               AS desig_year,
         %s
  FROM read_parquet('%s')
  WHERE year IS NOT NULL
  GROUP BY grid_id", attr_cols, panel_path)))

# Attach the crop variables from the Stage 1b cross-section.
if (file.exists(cropland_final_file)) {
  crop_avail <- arrow::open_dataset(cropland_final_file, format = "parquet")$schema$names
  crop_want  <- intersect(BAL_VARS_CROP, crop_avail)
  if (length(crop_want) < length(BAL_VARS_CROP)) {
    say("  WARNING: absent from the cropland cross-section: %s",
        paste(setdiff(BAL_VARS_CROP, crop_avail), collapse = ", "))
  }
  crop_dt <- as.data.table(arrow::read_parquet(
    cropland_final_file, col_select = arrow::all_of(c("grid_id", crop_want))))
  say("  cropland cross-section: %s cells, %d crop variables",
      format(nrow(crop_dt), big.mark = ","), length(crop_want))
  cell_attr <- merge(cell_attr, crop_dt, by = "grid_id", all.x = TRUE)
  rm(crop_dt)
} else {
  say("  WARNING: %s not found; crop variables dropped from the balance table",
      basename(cropland_final_file))
  BAL_VARS <- BAL_VARS_PANEL
}

cell_attr <- merge(cell_attr, biome_dt, by = "grid_id", all.x = TRUE)
cell_attr <- cell_attr[!is.na(biome) & !is.na(country_iso3)]
cell_attr[, cb := paste(country_iso3, biome, sep = " | ")]

# Treated = ever protected with a designation year inside the panel window.
# Never-treated = never protected in any year, and no designation year.
cell_attr[, treated_cell := ever_protected == 1 & !is.na(desig_year)]
cell_attr[, never_cell   := ever_protected == 0 & is.na(desig_year)]

say("Cells: %s treated, %s never-treated, %s neither (protected but no desig year, or partial)",
    format(sum(cell_attr$treated_cell), big.mark = ","),
    format(sum(cell_attr$never_cell),   big.mark = ","),
    format(sum(!cell_attr$treated_cell & !cell_attr$never_cell), big.mark = ","))

# Cohorts must have a usable pre-period and post-period inside 1990-2023.
panel_yrs <- as.data.table(dbGetQuery(con, sprintf(
  "SELECT min(year) AS y0, max(year) AS y1 FROM read_parquet('%s')", panel_path)))
Y0 <- panel_yrs$y0; Y1 <- panel_yrs$y1

cohorts <- sort(unique(cell_attr[treated_cell == TRUE, desig_year]))
cohorts <- cohorts[cohorts > Y0 & cohorts <= Y1]
say("Cohorts (designation years) in %d-%d: %d", Y0, Y1, length(cohorts))

# ------------------------------------------------------------------------------
# BUILD THE STACK
#
# For each cohort c: treated cells designated in c, plus never-treated cells
# drawn from the SAME country x biome as those treated cells. Controls are
# clean - never treated at any point in the panel - so no already-treated or
# later-treated units enter any comparison.
# ------------------------------------------------------------------------------

say("Building stack...")
stack_members <- list()

for (cc in cohorts) {
  trt <- cell_attr[treated_cell == TRUE & desig_year == cc, .(grid_id, cb)]
  if (nrow(trt) == 0) next
  cbs <- unique(trt$cb)

  ctl_pool <- cell_attr[never_cell == TRUE & cb %in% cbs, .(grid_id, cb)]
  n_pool <- nrow(ctl_pool)
  if (n_pool == 0) {
    say("  cohort %d: %s treated, NO never-treated controls in its country-biomes - cohort dropped",
        cc, format(nrow(trt), big.mark = ","))
    next
  }
  if (n_pool > CONTROL_CAP) {
    ctl_pool <- ctl_pool[sample(.N, CONTROL_CAP)]
  }

  stack_members[[as.character(cc)]] <- rbind(
    data.table(grid_id = trt$grid_id,      cohort = cc, treated = 1L),
    data.table(grid_id = ctl_pool$grid_id, cohort = cc, treated = 0L))

  say("  cohort %d: %s treated, %s controls kept of %s in pool%s",
      cc, format(nrow(trt), big.mark = ","), format(nrow(ctl_pool), big.mark = ","),
      format(n_pool, big.mark = ","),
      if (n_pool > CONTROL_CAP) sprintf(" (capped, %s not used)",
        format(n_pool - CONTROL_CAP, big.mark = ",")) else "")
}

stack_map <- rbindlist(stack_members)
say("Stack membership rows: %s across %d cohorts",
    format(nrow(stack_map), big.mark = ","), uniqueN(stack_map$cohort))

# ------------------------------------------------------------------------------
# PULL THE OUTCOME PANEL FOR ONLY THE CELLS THAT APPEAR IN SOME STACK
# ------------------------------------------------------------------------------

needed <- unique(stack_map$grid_id)
say("Pulling outcome panel for %s distinct cells...", format(length(needed), big.mark = ","))

dbWriteTable(con, "needed_cells", data.frame(grid_id = needed), overwrite = TRUE)
panel <- as.data.table(dbGetQuery(con, sprintf("
  SELECT p.grid_id, p.year, %s, p.%s AS pretrend_y
  FROM read_parquet('%s') p
  INNER JOIN needed_cells n ON p.grid_id = n.grid_id
  WHERE p.year IS NOT NULL",
  paste(sprintf("p.%s AS %s", OUTCOMES, OUTCOMES), collapse = ", "),
  PRETREND_VAR, panel_path)))
say("  %s cell-year rows", format(nrow(panel), big.mark = ","))

dbDisconnect(con, shutdown = TRUE)

# ------------------------------------------------------------------------------
# ASSEMBLE: one row per (cell, cohort, year) within the event window
# ------------------------------------------------------------------------------

stk <- merge(stack_map, panel, by = "grid_id", allow.cartesian = TRUE)
stk[, rel_year := year - cohort]
stk <- stk[rel_year >= EVENT_MIN & rel_year <= EVENT_MAX]
stk[, id_cohort := paste(grid_id, cohort, sep = "_")]
say("Stacked rows in event window: %s", format(nrow(stk), big.mark = ","))

# Keep only units observed at the reference period, since every normalisation
# and the ref = -1 event study depend on it.
has_ref <- stk[rel_year == -1, .(id_cohort = unique(id_cohort))]
n_before <- uniqueN(stk$id_cohort)
stk <- stk[id_cohort %in% has_ref$id_cohort]
say("Units dropped for missing t = -1: %s of %s",
    format(n_before - uniqueN(stk$id_cohort), big.mark = ","),
    format(n_before, big.mark = ","))

# ------------------------------------------------------------------------------
# RAW TRENDS: remove the ID x cohort intercept at t = -1
# ------------------------------------------------------------------------------

stk <- merge(stk, cell_attr[, .(grid_id, biome, country_iso3, cb)], by = "grid_id", all.x = TRUE)
stk[, post := as.integer(rel_year >= 0)]
stk[, treat_post := treated * post]

for (OUTCOME in OUTCOMES) {
  say("--- outcome: %s ---", OUTCOME)
    stk[, y := get(OUTCOME)]
    say("Building raw trends...")
  stk[, y_ref := y[rel_year == -1][1], by = id_cohort]
  stk[, y_norm := y - y_ref]

  trends <- stk[, .(mean_y = mean(y_norm, na.rm = TRUE),
                    se = sd(y_norm, na.rm = TRUE) / sqrt(.N),
                    n = .N),
                by = .(rel_year, treated)]
  trends[, group := fifelse(treated == 1L, "Treated (PA designated)", "Never treated")]
  setorder(trends, treated, rel_year)

  p_raw <- ggplot(trends, aes(x = rel_year, y = mean_y, colour = group, fill = group)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    geom_vline(xintercept = -0.5, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
    geom_ribbon(aes(ymin = mean_y - 1.96 * se, ymax = mean_y + 1.96 * se),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
    scale_colour_manual(values = c("Treated (PA designated)" = "steelblue",
                                   "Never treated" = "grey35")) +
    scale_fill_manual(values = c("Treated (PA designated)" = "steelblue",
                                 "Never treated" = "grey35")) +
    labs(x = "Years since designation",
         y = sprintf("%s, normalised to t = -1", OUTCOME),
         title = "Raw trends, stacked by cohort",
         subtitle = sprintf("ID x cohort intercept removed at t = -1; controls never treated, same country x biome"),
         colour = NULL, fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())

  ggsave(file.path(out_dir, sprintf("raw_trends_stacked_%s.png", OUTCOME)), p_raw,
         width = 9, height = 5.5, dpi = 200, bg = "white")
  fwrite(trends, file.path(out_dir, sprintf("raw_trends_stacked_%s.csv", OUTCOME)))
  say("  wrote raw_trends_stacked_%s.png / .csv", OUTCOME)

  # ------------------------------------------------------------------------------
  # DIFFERENCE-IN-DIFFERENCES
  #   (1) ID x cohort FE + cohort x year FE
  #   (2) ID x cohort FE + cohort x biome x year FE
  # ------------------------------------------------------------------------------

  say("Estimating DiD...")
  # Logged AFTER the merge that creates cb. Reported before it, uniqueN(NULL)
  # printed "0 clusters" and looked like a broken clustering variable.
  say("  clustering on country x biome: %d clusters (%s rows with NA cb)",
      uniqueN(stk$cb), format(sum(is.na(stk$cb)), big.mark = ","))

  # Clustered on country x biome, the stratum controls are drawn from. NOT on
  # grid_id: treatment is assigned at the protected area, so every cell in a PA
  # flips together and cell-level clustering would treat thousands of mechanically
  # identical observations as independent (Moulton). Country x biome absorbs both
  # the within-PA correlation and the spatial correlation among neighbouring
  # control cells.
  m_static_1 <- feols(y ~ treat_post | id_cohort + cohort^year,
                      data = stk, cluster = ~cb, lean = TRUE)
  m_static_2 <- feols(y ~ treat_post | id_cohort + cohort^biome^year,
                      data = stk, cluster = ~cb, lean = TRUE)

  m_dyn_1 <- feols(y ~ i(rel_year, treated, ref = -1) | id_cohort + cohort^year,
                   data = stk, cluster = ~cb)
  m_dyn_2 <- feols(y ~ i(rel_year, treated, ref = -1) | id_cohort + cohort^biome^year,
                   data = stk, cluster = ~cb)

  etable(m_static_1, m_static_2,
         file = file.path(out_dir, sprintf("did_estimates_%s.tex", OUTCOME)), replace = TRUE,
         title = "Stacked difference-in-differences: effect of PA designation",
         style.tex = style.tex("aer"),
         dict = c(treat_post = "Treated $\\times$ Post", y = OUTCOME,
                  id_cohort = "ID $\\times$ cohort", `cohort^year` = "Cohort $\\times$ year",
                  `cohort^biome^year` = "Cohort $\\times$ biome $\\times$ year"))
  say("  wrote did_estimates_%s.tex", OUTCOME)

  dyn <- rbindlist(lapply(list(`Cohort x year` = m_dyn_1,
                               `Cohort x biome x year` = m_dyn_2), function(m) {
    ct <- as.data.table(coeftable(m), keep.rownames = "term")
    ct <- ct[grepl("rel_year::", term)]
    ct[, rel_year := as.integer(sub(".*rel_year::(-?\\d+).*", "\\1", term))]
    ct[, .(rel_year, est = Estimate, se = `Std. Error`)]
  }), idcol = "spec")
  dyn <- rbind(dyn, data.table(spec = unique(dyn$spec), rel_year = -1, est = 0, se = 0))

  p_es <- ggplot(dyn, aes(x = rel_year, y = est, colour = spec)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    geom_vline(xintercept = -0.5, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
    geom_linerange(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se),
                   position = position_dodge(width = 0.4)) +
    geom_point(position = position_dodge(width = 0.4), size = 1.8) +
    scale_colour_manual(values = c("Cohort x year" = "steelblue",
                                   "Cohort x biome x year" = "darkorange")) +
    labs(x = "Years since designation", y = sprintf("Effect on %s", OUTCOME),
         title = "Stacked event study", colour = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())

  ggsave(file.path(out_dir, sprintf("event_study_stacked_%s.png", OUTCOME)), p_es,
         width = 9, height = 5.5, dpi = 200, bg = "white")
  fwrite(dyn, file.path(out_dir, sprintf("event_study_stacked_%s.csv", OUTCOME)))
  say("  wrote event_study_stacked_%s.png / .csv", OUTCOME)

  # ------------------------------------------------------------------------------
  # BALANCE TABLE: treated vs never-treated
  #
  # Time-invariant characteristics, plus the pre-treatment 5-year trend in
  # deforestation. Units without at least MIN_PRE_YEARS pre-period observations
  # are dropped from the trend row only, not from the rest of the table.
  # ------------------------------------------------------------------------------
}

say("Building balance table...")

# The requirement is MIN_PRE_YEARS of USABLE pre-period observations: counting
# rows alone let a unit whose pretrend_y is entirely NA reach lm(), which fails
# with "0 (non-NA) cases" and killed the first run.
slope_safe <- function(yv, xv) {
  ok <- !is.na(yv) & !is.na(xv)
  if (sum(ok) < MIN_PRE_YEARS) return(NA_real_)
  if (length(unique(xv[ok])) < 2) return(NA_real_)
  unname(coef(lm(yv[ok] ~ xv[ok]))[2])
}

pre <- stk[rel_year >= -MIN_PRE_YEARS & rel_year <= -1,
           .(n_pre = sum(!is.na(pretrend_y)),
             slope = slope_safe(pretrend_y, rel_year)),
           by = .(id_cohort, grid_id, treated)]
n_trend_drop <- pre[is.na(slope), .N]
say("  pre-trend: %s of %s unit-cohorts dropped for < %d pre-period years",
    format(n_trend_drop, big.mark = ","), format(nrow(pre), big.mark = ","), MIN_PRE_YEARS)

bal_cells <- unique(stk[, .(grid_id, treated)])
bal <- merge(bal_cells, cell_attr[, c("grid_id", "cb", BAL_VARS), with = FALSE],
             by = "grid_id", all.x = TRUE)

# Percentile transform, applied to the pooled sample so treated and control
# are ranked on one common distribution.
for (v in intersect(PCTL_VARS, names(bal))) {
  x <- bal[[v]]
  ok <- !is.na(x)
  pct <- rep(NA_real_, length(x))
  pct[ok] <- 100 * ecdf(x[ok])(x[ok])
  bal[[v]] <- pct
  say("  %s reported as percentile of the pooled distribution", v)
}

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d, big.mark = ",")

# Pick a display scale per row. Range-size rarity is a sum of inverse species
# range areas, so its values are ~1e-6 and print as "0.000" at three decimals.
# Rows whose magnitude is below 0.001 are rescaled by a power of ten and the
# multiplier is stated in the row label, rather than padding the whole table out
# to eight decimals for one variable.
row_scale <- function(vals, label) {
  m <- suppressWarnings(max(abs(vals[is.finite(vals)]), na.rm = TRUE))
  if (!is.finite(m) || m == 0) return(list(mult = 1, label = label))
  if (m < 0.001) {
    k <- ceiling(-log10(m)) + 2
    return(list(mult = 10^k,
                label = sprintf("%s ($\\times 10^{-%d}$)", label, k)))
  }
  list(mult = 1, label = label)
}

# Differences and their standard errors come from a regression on the treated
# indicator clustered at country x biome, matching the DiD specifications. An
# unclustered t.test here would understate the SEs for the same Moulton reason.
diff_clustered <- function(dt, v) {
  d <- dt[!is.na(get(v))]
  if (!nrow(d) || uniqueN(d$treated) < 2) return(NULL)
  m <- feols(as.formula(sprintf("`%s` ~ treated", v)), data = d, cluster = ~cb)
  ct <- coeftable(m)
  list(diff = ct["treated", "Estimate"], se = ct["treated", "Std. Error"],
       p = ct["treated", "Pr(>|t|)"])
}

rows <- lapply(BAL_VARS, function(v) {
  t1 <- bal[treated == 1L][[v]]; t0 <- bal[treated == 0L][[v]]
  t1 <- t1[!is.na(t1)]; t0 <- t0[!is.na(t0)]
  if (!length(t1) || !length(t0)) return(NULL)
  cl <- diff_clustered(bal, v)
  if (is.null(cl)) return(NULL)
  sc <- row_scale(c(mean(t1), mean(t0), cl$diff), BAL_LABELS[[v]])
  # Normalised difference (Imbens & Rubin): scale-free, so it is comparable
  # across variables whose units are not. Computed on RAW values - it is
  # invariant to the display rescaling above.
  nd <- (mean(t1) - mean(t0)) / sqrt((var(t1) + var(t0)) / 2)
  data.table(label = sc$label,
             m1 = mean(t1) * sc$mult, s1 = sd(t1) * sc$mult,
             m0 = mean(t0) * sc$mult, s0 = sd(t0) * sc$mult,
             diff = cl$diff * sc$mult, se = cl$se * sc$mult, p = cl$p, nd = nd)
})
bal_tab <- rbindlist(rows)

pt <- merge(pre[!is.na(slope)], cell_attr[, .(grid_id, cb)], by = "grid_id", all.x = TRUE)
t1 <- pt[treated == 1L, slope]; t0 <- pt[treated == 0L, slope]
cl <- diff_clustered(pt, "slope")
sc <- row_scale(c(mean(t1), mean(t0), cl$diff),
                sprintf("Pre-treatment %d-year trend in %s", MIN_PRE_YEARS, PRETREND_VAR))
trend_row <- data.table(
  label = sc$label, m1 = mean(t1) * sc$mult, s1 = sd(t1) * sc$mult,
  m0 = mean(t0) * sc$mult, s0 = sd(t0) * sc$mult,
  diff = cl$diff * sc$mult, se = cl$se * sc$mult, p = cl$p,
  nd = (mean(t1) - mean(t0)) / sqrt((var(t1) + var(t0)) / 2))

stars <- function(p) fifelse(p < 0.01, "$^{***}$", fifelse(p < 0.05, "$^{**}$",
                     fifelse(p < 0.1, "$^{*}$", "")))

tex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Balance between treated and never-treated cells}",
  "\\label{tab:balance_treated_never}",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  " & Treated & Never treated & Difference & Norm. diff. \\\\",
  " & (1) & (2) & (3) & (4) \\\\",
  "\\midrule",
  "\\multicolumn{5}{l}{\\emph{Panel A: Time-invariant characteristics}} \\\\",
  # interleave: estimate row then its SE row, per variable
  c(rbind(
    sprintf("%s & %s & %s & %s%s & %s \\\\", bal_tab$label, fmt(bal_tab$m1),
            fmt(bal_tab$m0), fmt(bal_tab$diff), stars(bal_tab$p), fmt(bal_tab$nd, 2)),
    sprintf(" & (%s) & (%s) & (%s) & \\\\", fmt(bal_tab$s1), fmt(bal_tab$s0),
            fmt(bal_tab$se)))),
  "\\midrule",
  "\\multicolumn{5}{l}{\\emph{Panel B: Pre-treatment trend}} \\\\",
  sprintf("%s & %s & %s & %s%s & %s \\\\", trend_row$label, fmt(trend_row$m1),
          fmt(trend_row$m0), fmt(trend_row$diff), stars(trend_row$p), fmt(trend_row$nd, 2)),
  sprintf(" & (%s) & (%s) & (%s) & \\\\", fmt(trend_row$s1), fmt(trend_row$s0),
          fmt(trend_row$se)),
  "\\midrule",
  sprintf("Cells & %s & %s & & \\\\",
          format(bal[treated == 1L, .N], big.mark = ","),
          format(bal[treated == 0L, .N], big.mark = ",")),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{0.92\\textwidth}",
  "\\footnotesize",
  "\\vspace{0.5em}",
  sprintf("\\emph{Notes:} The unit of observation is a 5\\,km grid cell. Column (1) reports means for cells in a protected area with a designation year between %d and %d; column (2) reports means for cells that are never protected at any point in the panel. Standard deviations are in parentheses below each mean. Column (3) reports the difference in means with its standard error in parentheses. Column (4) reports the normalised difference, $(\\bar{x}_1-\\bar{x}_0)/\\sqrt{(s_1^2+s_0^2)/2}$, which is unit-free and therefore comparable across rows; values above 0.25 in absolute value are conventionally treated as substantial imbalance. Never-treated controls are drawn from the same country $\\times$ biome as the treated cells of each cohort, where biome is the RESOLVE Ecoregions 2017 classification. Control cells are randomly sampled within each country $\\times$ biome, capped at %s per cohort. Panel B reports the slope from a cell-level regression of %s on the year index over the %d years preceding designation; cells with fewer than %d pre-period observations are excluded from that row, which drops %s of %s unit-cohorts. Standard errors are clustered at the country $\\times$ biome level. Range-size rarity is reported as a percentile of the pooled treated-plus-control distribution; its native units are a sum of inverse species range areas and are not directly interpretable. Stars denote %s.",
          min(cohorts), max(cohorts), format(CONTROL_CAP, big.mark = ","),
          PRETREND_VAR, MIN_PRE_YEARS, MIN_PRE_YEARS,
          format(n_trend_drop, big.mark = ","), format(nrow(pre), big.mark = ","),
          "$^{*}p<0.1$, $^{**}p<0.05$, $^{***}p<0.01$"),
  "\\end{minipage}",
  "\\end{table}")

writeLines(tex, file.path(out_dir, "balance_treated_vs_never.tex"))
say("  wrote balance_treated_vs_never.tex")

writeLines(log_lines, file.path(out_dir, "stacked_did_sample.txt"))
say("Done. Outputs in %s", out_dir)
