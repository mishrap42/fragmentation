# ==============================================================================
# 2_pa_event_study.R
# ==============================================================================
# Per-PA dynamic event studies. Each protected area's cells are matched 1:1,
# nearest-neighbour WITHOUT replacement, to never-treated donor cells in the
# same country x biome; a dynamic event study is then estimated on that PA's
# matched sample alone.
#
# OUTPUT (one row per PA x event time):
#   wdpaid  year_treated  country  event  year  coefficient  std_error
#   plus n_treated / n_matched, kept because a coefficient from 8 matched pairs
#   should not be read the same way as one from 800.
#   The reference period (event = -1) is NOT emitted: it is normalised to zero
#   by construction and has no standard error.
#
# PARALLELISM: SLURM array over COUNTRIES. Donor pools never cross a country x
# biome boundary, so a country is self-contained and needs no cross-task state.
# Out-of-range task ids exit 0. Set PA_ES_CONSOLIDATE=1 to stitch the shards.
#
# DROP RULE - exactly one, as specified: a PA is dropped only when its country x
# biome cannot supply at least one donor per treated cell (matching is without
# replacement, so n_donors >= n_treated is required). Everything else is kept
# and reported, however thin.
#
# Usage:
#   sbatch --array=1-140 code/bash/analysis/2_pa_event_study.sh
#   sbatch --export=ALL,PA_ES_CONSOLIDATE=1 code/bash/analysis/2_pa_event_study.sh
# ==============================================================================

suppressMessages({
  library(data.table); library(duckdb); library(DBI)
  library(fixest); library(MatchIt); library(arrow)
})

project_root <- if (Sys.getenv("SLURM_SUBMIT_DIR") != "") Sys.getenv("SLURM_SUBMIT_DIR") else here::here()
panel_path   <- file.path(project_root, "Data/build/final/TMF_5km_panel.parquet")
crop_path    <- file.path(project_root, "Data/build/final/TMF_5km_cropland.parquet")
biome_cache  <- file.path(project_root, "Data/build/final/grid_biome.parquet")
shard_dir    <- file.path(project_root, "Data/build/intermediate/pa_event_study")
out_dir      <- file.path(project_root, "output", "analysis", "pa_event_study")
dir.create(shard_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir,   recursive = TRUE, showWarnings = FALSE)

OUTCOME    <- Sys.getenv("PA_ES_OUTCOME", "forest_cover")
EVENT_MIN  <- as.integer(Sys.getenv("PA_ES_EVENT_MIN", "-10"))
EVENT_MAX  <- as.integer(Sys.getenv("PA_ES_EVENT_MAX", "15"))
PRE_LO     <- -5L   # matching window on pre-period forest, inclusive
PRE_HI     <- -1L
CONLEY_KM  <- 50
CONSOLIDATE <- Sys.getenv("PA_ES_CONSOLIDATE", "0") == "1"

say <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")

# ==============================================================================
# CONSOLIDATION MODE
# ==============================================================================
if (CONSOLIDATE) {
  fs <- list.files(shard_dir, pattern = "^pa_es_.*\\.parquet$", full.names = TRUE)
  say("Consolidating %d shards", length(fs))
  stopifnot(length(fs) > 0)
  all <- rbindlist(lapply(fs, function(f) as.data.table(read_parquet(f))), fill = TRUE)
  setorder(all, country, wdpaid, event)
  write_parquet(all, file.path(out_dir, "pa_event_study.parquet"))
  fwrite(all, file.path(out_dir, "pa_event_study.csv"))
  say("PAs: %s | rows: %s | countries: %d",
      format(uniqueN(all$wdpaid), big.mark = ","), format(nrow(all), big.mark = ","),
      uniqueN(all$country))
  quit(save = "no", status = 0)
}

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

# ==============================================================================
# CELL ATTRIBUTES
# ==============================================================================
con <- dbConnect(duckdb(), shared_home = FALSE)
dbExecute(con, "SET memory_limit='60GB'"); dbExecute(con, "SET threads=4")

say("Loading cell attributes...")
cell <- as.data.table(dbGetQuery(con, sprintf("
  SELECT grid_id,
         any_value(country_iso3)                       AS country,
         any_value(centroid_lon)                       AS centroid_lon,
         any_value(centroid_lat)                       AS centroid_lat,
         max(CASE WHEN is_protected THEN 1 ELSE 0 END) AS ever_protected,
         min(desig_year)                               AS desig_year,
         any_value(wdpa_pid)                           AS wdpaid,
         any_value(travel_time_cities)                 AS travel_time_cities,
         any_value(pop_density_1990)                   AS pop_density_1990
  FROM read_parquet('%s') WHERE year IS NOT NULL GROUP BY grid_id", panel_path)))

# yields live in the Stage 1b cross-section, not the panel
cell <- merge(cell, as.data.table(read_parquet(crop_path,
                col_select = all_of(c("grid_id", "yield_maize")))),
              by = "grid_id", all.x = TRUE)
cell <- merge(cell, as.data.table(read_parquet(biome_cache)), by = "grid_id", all.x = TRUE)
cell <- cell[!is.na(biome) & !is.na(country) & !is.na(centroid_lon)]
cell[, cb := paste(country, biome, sep = " | ")]
cell[, treated_cell := ever_protected == 1 & !is.na(desig_year) & !is.na(wdpaid)]
cell[, never_cell   := ever_protected == 0 & is.na(desig_year)]

countries <- sort(unique(cell[treated_cell == TRUE, country]))
if (task_id > length(countries)) {
  say("task %d > %d countries with treated cells; nothing to do", task_id, length(countries))
  dbDisconnect(con, shutdown = TRUE); quit(save = "no", status = 0)
}
CTRY <- countries[task_id]
say("Country %d/%d: %s", task_id, length(countries), CTRY)

cell <- cell[country == CTRY]
pas  <- cell[treated_cell == TRUE, .N, by = .(wdpaid, desig_year)][order(wdpaid)]
say("  %s cells | %d PAs | %s never-treated",
    format(nrow(cell), big.mark = ","), nrow(pas),
    format(sum(cell$never_cell), big.mark = ","))
if (nrow(pas) == 0) { dbDisconnect(con, shutdown = TRUE); quit(save = "no", status = 0) }

# ==============================================================================
# OUTCOME PANEL for this country only
# ==============================================================================
dbWriteTable(con, "want", data.frame(grid_id = cell$grid_id), overwrite = TRUE)
pan <- as.data.table(dbGetQuery(con, sprintf("
  SELECT p.grid_id, p.year, p.%s AS y
  FROM read_parquet('%s') p INNER JOIN want w ON p.grid_id = w.grid_id
  WHERE p.year IS NOT NULL", OUTCOME, panel_path)))
Y0 <- min(pan$year); Y1 <- max(pan$year)
dbDisconnect(con, shutdown = TRUE)
say("  panel: %s rows, %d-%d", format(nrow(pan), big.mark = ","), Y0, Y1)

setkey(pan, grid_id, year)

# ==============================================================================
# PER-PA MATCH + EVENT STUDY
# ==============================================================================
# Matching covariates, per the specification:
#   forest share LEVELS at each of t-5 ... t-1        (5 targets)
#   forest share TREND over that same window          (1 target)
#   travel time, population, yield                    (3 targets)
# [-10, -6] is deliberately left untargeted so it can serve as an out-of-sample
# check on whether the match holds where it was not fitted.
MATCH_STATIC <- c("travel_time_cities", "pop_density_1990", "yield_maize")

pre_features <- function(ids, c_year) {
  w <- pan[.(ids), on = "grid_id"][year >= c_year + PRE_LO & year <= c_year + PRE_HI]
  w[, rel := year - c_year]
  lv <- dcast(w, grid_id ~ rel, value.var = "y")
  # lev_m5 ... lev_m1, NOT lev-5: a minus sign makes the name invalid inside a
  # model formula and matchit() fails with "invalid model formula in ExtractVars".
  setnames(lv, setdiff(names(lv), "grid_id"), sprintf("lev_m%d", abs(PRE_LO:PRE_HI)))
  tr <- w[, .(trend = if (sum(!is.na(y)) >= 2) coef(lm(y ~ rel))[2] else NA_real_), by = grid_id]
  merge(lv, tr, by = "grid_id")
}

results <- list()
for (i in seq_len(nrow(pas))) {
  pid <- pas$wdpaid[i]; c_year <- pas$desig_year[i]

  trt <- cell[treated_cell == TRUE & wdpaid == pid]
  don <- cell[never_cell == TRUE & cb %in% unique(trt$cb)]

  # THE ONLY DROP RULE
  if (nrow(don) < nrow(trt)) {
    say("  PA %s (%d): %d treated, only %d donors in its country x biome - DROPPED",
        pid, c_year, nrow(trt), nrow(don))
    next
  }

  e_lo <- max(EVENT_MIN, Y0 - c_year)
  e_hi <- min(EVENT_MAX, Y1 - c_year)
  if (e_lo > PRE_LO || e_hi < 1) next          # no usable pre-window or no post

  feat <- pre_features(c(trt$grid_id, don$grid_id), c_year)
  md <- merge(rbind(trt[, .(grid_id, treated = 1L)], don[, .(grid_id, treated = 0L)]),
              feat, by = "grid_id")
  md <- merge(md, cell[, c("grid_id", MATCH_STATIC), with = FALSE], by = "grid_id")
  md <- md[complete.cases(md)]
  if (md[treated == 1L, .N] == 0 || md[treated == 0L, .N] < md[treated == 1L, .N]) next

  mvars <- c(sprintf("lev_m%d", abs(PRE_LO:PRE_HI)), "trend", MATCH_STATIC)
  # standardise so no covariate dominates the Euclidean metric by its units
  for (v in mvars) { s <- sd(md[[v]]); md[, (v) := if (s > 0) (get(v) - mean(get(v))) / s else 0] }

  mo <- tryCatch(matchit(as.formula(paste("treated ~", paste(mvars, collapse = " + "))),
                         data = md, method = "nearest", distance = "euclidean",
                         replace = FALSE),
                 error = function(e) NULL)
  if (is.null(mo)) { say("  PA %s: matchit failed", pid); next }
  mm <- as.data.table(match.data(mo))[, .(grid_id, treated)]

  es <- merge(pan[.(mm$grid_id), on = "grid_id"], mm, by = "grid_id")
  es[, event := year - c_year]
  es <- es[event >= e_lo & event <= e_hi]
  es <- merge(es, cell[, .(grid_id, centroid_lon, centroid_lat)], by = "grid_id")
  if (uniqueN(es$event) < 3) next

  m <- tryCatch(feols(y ~ i(event, treated, ref = -1) | grid_id + year, data = es,
                      vcov = vcov_conley(lat = ~centroid_lat, lon = ~centroid_lon,
                                         cutoff = CONLEY_KM)),
                error = function(e) NULL)
  if (is.null(m)) { say("  PA %s: feols failed", pid); next }

  ct <- as.data.table(coeftable(m), keep.rownames = "term")[grepl("event::", term)]
  if (!nrow(ct)) next
  ct[, event := as.integer(sub(".*event::(-?\\d+).*", "\\1", term))]

  results[[length(results) + 1]] <- data.table(
    wdpaid = pid, year_treated = c_year, country = CTRY,
    event = ct$event, year = c_year + ct$event,
    coefficient = ct$Estimate, std_error = ct$`Std. Error`,
    n_treated = mm[treated == 1L, .N], n_matched = mm[treated == 0L, .N])

  if (i %% 25 == 0) say("  ... %d/%d PAs", i, nrow(pas))
}

out <- if (length(results)) rbindlist(results) else data.table(
  wdpaid = character(0), year_treated = integer(0), country = character(0),
  event = integer(0), year = integer(0), coefficient = numeric(0),
  std_error = numeric(0), n_treated = integer(0), n_matched = integer(0))

shard <- file.path(shard_dir, sprintf("pa_es_%s.parquet", CTRY))
write_parquet(out, shard)
say("Wrote %s: %d PAs, %s rows", basename(shard), uniqueN(out$wdpaid),
    format(nrow(out), big.mark = ","))
