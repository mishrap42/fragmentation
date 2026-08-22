# ==============================================================================
# STAGE 6: ASSEMBLE FINAL DATASET
# ==============================================================================
# This script merges all pipeline outputs into the final panel dataset.
#
# Input: Grid metadata (includes WDPA), TMF time series, classifications
# Output: Data/build/final/TMF_5km_panel.parquet
#
# Note: WDPA data is baked into grid cells at Stage 0 - each cell is either
# fully protected or not (grid is cut on WDPA boundaries).
# ==============================================================================

# Load configuration
here::i_am('code/build/6_assemble_final.R')
source("code/build/BUILD_workspace.R")

# Record start time
start_time <- Sys.time()

log_job_start("6_assemble_final.R", task_id = 1)

# Check if output already exists
output_file <- panel_final_file
skip_if_exists(output_file, "final assembly")

# ==============================================================================
# LOAD ALL GRID METADATA (from sub-tile files)
# ==============================================================================

log_message("Loading grid metadata...")

# Get all sub-tile grid files
grid_files <- sapply(1:N_SUB_TILES, function(s) get_grid_filename(s, "parquet"))
files_exist <- file.exists(grid_files)

log_message(sprintf("Found %d/%d sub-tile grid files", sum(files_exist), N_SUB_TILES))

grid_list <- lapply(grid_files[files_exist], function(f) {
  dt <- read_parquet(f)
  setDT(dt)
  return(dt)
})

grid_all <- rbindlist(grid_list, fill = TRUE)
log_message(sprintf("Total grid cells: %s from %d sub-tiles",
                    format(nrow(grid_all), big.mark = ","),
                    sum(files_exist)))

rm(grid_list)
gc()

# ==============================================================================
# LOAD TMF TIME SERIES (using Arrow dataset API for memory efficiency)
# ==============================================================================

log_message("Loading TMF time series from tile-level parquet files...")

# Use Arrow's dataset API to read directly from tile parquet files
# This avoids the 2.1 billion row limit of rbindlist
tile_files <- sapply(1:N_TMF_TILES, get_consolidated_tile_filename)
files_exist <- file.exists(tile_files)

log_message(sprintf("Found %d/%d consolidated tile files", sum(files_exist), N_TMF_TILES))

if (sum(files_exist) == 0) {
  stop("No consolidated tile files found. Run Stage 2a first.")
}

# Open as Arrow dataset (lazy evaluation - doesn't load into memory yet)
tmf_dataset <- arrow::open_dataset(tile_files[files_exist], format = "parquet")

log_message("Arrow dataset opened. Processing in chunks to reshape...")

# Process and reshape in chunks by tile to manage memory
# Each tile produces a wide-format data.table
tmf_wide_list <- lapply(which(files_exist), function(i) {
  tile_file <- tile_files[i]

  # Read single tile
  tile_data <- arrow::read_parquet(tile_file)
  setDT(tile_data)

  if (nrow(tile_data) == 0) return(NULL)

  # Reshape to wide format within the tile
  tile_wide <- dcast(tile_data, grid_id + year ~ tmf_class,
                     value.var = "fraction", fill = 0)

  rm(tile_data)
  return(tile_wide)
})

# Remove NULL entries
tmf_wide_list <- tmf_wide_list[!sapply(tmf_wide_list, is.null)]

log_message(sprintf("Processed %d tiles, combining...", length(tmf_wide_list)))

# Combine wide-format tiles (much fewer rows than long format)
tmf_wide <- rbindlist(tmf_wide_list, use.names = TRUE, fill = TRUE)

log_message(sprintf("TMF wide format: %s cell-years",
                    format(nrow(tmf_wide), big.mark = ",")))

rm(tmf_wide_list)
gc()

# ==============================================================================
# LOAD INTERIOR CLASSIFICATIONS
# ==============================================================================

log_message("Loading interior classifications...")

interior_files <- sapply(1:N_TMF_TILES, get_interior_filename)
files_exist <- file.exists(interior_files)

interior_list <- lapply(interior_files[files_exist], function(f) {
  dt <- read_parquet(f)
  setDT(dt)
  return(dt)
})

interior_all <- rbindlist(interior_list, fill = TRUE)
interior_all <- unique(interior_all, by = "grid_id")

log_message(sprintf("Interior classifications: %s cells (%s interior)",
                    format(nrow(interior_all), big.mark = ","),
                    format(sum(interior_all$is_interior), big.mark = ",")))

rm(interior_list)
gc()

# ==============================================================================
# LOAD FRONTIER CLASSIFICATIONS
# ==============================================================================

log_message("Loading frontier classifications...")

frontier_files <- sapply(1:N_TMF_TILES, get_frontier_filename)
files_exist <- file.exists(frontier_files)

frontier_list <- lapply(frontier_files[files_exist], function(f) {
  dt <- read_parquet(f)
  setDT(dt)
  return(dt)
})

frontier_all <- rbindlist(frontier_list, fill = TRUE)
frontier_all <- unique(frontier_all, by = "grid_id")

log_message(sprintf("Frontier classifications: %s cells (%s frontier)",
                    format(nrow(frontier_all), big.mark = ","),
                    format(sum(frontier_all$is_frontier, na.rm = TRUE), big.mark = ",")))

rm(frontier_list)
gc()

# Note: WDPA data is already in grid_all from Stage 0 (grid cells are cut on WDPA boundaries)
# Each cell has: is_protected, wdpa_pid, iucn_cat, desig_year, gov_type

# ==============================================================================
# YIELD / CROP-SHARE DATA IS *NOT* MERGED HERE
# ==============================================================================
# Stage 0c produces 345 EarthStat columns (172 yield_*, 172 cropshare_*,
# cropland_frac / pasture_frac / cropshare_total). All of it is a year-2000
# CROSS-SECTION. Merging it into this cell-year panel would replicate a
# time-invariant object across 34 years - roughly 66.8M x 345 x 8 bytes, ~184 GB
# uncompressed - to add exactly zero information per row.
#
# It is consolidated on the extraction side of the build, by Stage 1b
# (1b_assemble_cropland.R), into
#   Data/build/final/TMF_5km_cropland.parquet     (one row per grid_id)
# and joined on grid_id at analysis time. 1b and this stage are independent:
# neither waits on the other, and neither reads the other's output.
#
# The GAEZ-era yield_*200a_yld columns in any panel built before 2026-08-21 came
# through this merge under the pre-EarthStat config; they are gone by design,
# not by accident.

# ==============================================================================
# LOAD COVARIATE DATA (from Stage 5)
# ==============================================================================

log_message("Loading covariate data...")

covariate_files <- sapply(1:N_TMF_TILES, get_covariates_filename)
files_exist <- file.exists(covariate_files)

if (sum(files_exist) > 0) {
  log_message(sprintf("Found %d/%d covariate files", sum(files_exist), N_TMF_TILES))

  covariate_list <- lapply(covariate_files[files_exist], function(f) {
    dt <- read_parquet(f)
    setDT(dt)
    return(dt)
  })

  covariates_all <- rbindlist(covariate_list, fill = TRUE)
  log_message(sprintf("Covariate data: %s cells, %d variables",
                      format(nrow(covariates_all), big.mark = ","),
                      ncol(covariates_all) - 1))

  rm(covariate_list)
  gc()
} else {
  log_message("No covariate files found. Skipping covariate data.")
  covariates_all <- NULL
}

# ==============================================================================
# LOAD GFED BURNED AREA DATA (from Stage 2b)
# ==============================================================================

log_message("Loading GFED burned area data...")

gfed_files <- sapply(1:N_TMF_TILES, get_gfed_consolidated_filename)
files_exist <- file.exists(gfed_files)

if (sum(files_exist) > 0) {
  log_message(sprintf("Found %d/%d GFED tile files", sum(files_exist), N_TMF_TILES))

  gfed_list <- lapply(gfed_files[files_exist], function(f) {
    dt <- read_parquet(f)
    setDT(dt)
    return(dt)
  })

  gfed_all <- rbindlist(gfed_list, fill = TRUE)
  log_message(sprintf("GFED data: %s cell-years",
                      format(nrow(gfed_all), big.mark = ",")))

  rm(gfed_list)
  gc()
} else {
  log_message("No GFED files found. Skipping GFED data.")
  gfed_all <- NULL
}

# ==============================================================================
# MERGE CLASSIFICATIONS
# ==============================================================================

log_message("Merging grid metadata with classifications...")

# Start with grid metadata (includes WDPA and transport costs from Stage 0)
panel_base <- grid_all[, .(grid_id, tile_id, country_iso3, country_name,
                           centroid_lon, centroid_lat, area_km2,
                           is_protected, wdpa_pid, iucn_cat, desig_year, gov_type,
                           travel_time_cities, travel_time_ports)]

# Fill NA protection status (shouldn't happen but just in case)
panel_base[is.na(is_protected), is_protected := FALSE]

# Add interior classification
panel_base <- merge(
  panel_base,
  interior_all[, .(grid_id, is_interior, min_undisturbed_frac,
                   mean_undisturbed_frac, n_years_observed)],
  by = "grid_id",
  all.x = TRUE
)

# Fill NAs for cells not in interior data
panel_base[is.na(is_interior), is_interior := FALSE]

# Add frontier classification
panel_base <- merge(
  panel_base,
  frontier_all[, .(grid_id, is_frontier, dist_to_interior_km)],
  by = "grid_id",
  all.x = TRUE
)

# Fill NAs
panel_base[is.na(is_frontier), is_frontier := FALSE]

# Add covariate data if available
if (!is.null(covariates_all)) {
  log_message("Merging covariate data...")
  panel_base <- merge(
    panel_base,
    covariates_all,
    by = "grid_id",
    all.x = TRUE
  )
  rm(covariates_all)
  gc()
}

# Note: WDPA data already included from grid (no separate merge needed)

log_message(sprintf("Base panel: %s cells",
                    format(nrow(panel_base), big.mark = ",")))

# Clean up
rm(grid_all, interior_all, frontier_all)
gc()

# ==============================================================================
# CREATE PANEL DATASET (CROSS CELL x YEAR)
# ==============================================================================

log_message("Creating panel dataset (cell x year)...")

# Merge base with TMF time series
final_panel <- merge(
  panel_base,
  tmf_wide,
  by = "grid_id",
  all = TRUE  # Keep all cells and all TMF observations
)

log_message(sprintf("Panel after TMF merge: %s observations",
                    format(nrow(final_panel), big.mark = ",")))

# Clean up
rm(panel_base, tmf_wide)
gc()

# Add GFED burned area data (time-varying, merges on grid_id + year)
if (!is.null(gfed_all)) {
  log_message("Merging GFED burned area data...")

  final_panel <- merge(
    final_panel,
    gfed_all[, .(grid_id, year, burned_annual)],
    by = c("grid_id", "year"),
    all.x = TRUE
  )

  # Fill NA with 0 (no detected burning)
  final_panel[is.na(burned_annual), burned_annual := 0]

  log_message(sprintf("Cells with burns: %s",
                      format(sum(final_panel$burned_annual > 0), big.mark = ",")))

  rm(gfed_all)
  gc()
}

log_message(sprintf("Final panel: %s observations",
                    format(nrow(final_panel), big.mark = ",")))

gc_verbose()

# ==============================================================================
# CREATE DERIVED VARIABLES
# ==============================================================================

log_message("Creating derived variables...")

# Create zone variable
final_panel[, zone := fifelse(is_interior, "interior",
                               fifelse(is_frontier, "frontier", "other"))]

# Create forest cover variable (Undisturbed + Degraded + Regrowth)
forest_cols <- c("Undisturbed_TMF", "Degraded_TMF", "Regrowth")
existing_forest_cols <- intersect(forest_cols, names(final_panel))

if (length(existing_forest_cols) > 0) {
  final_panel[, forest_cover := rowSums(.SD, na.rm = TRUE),
              .SDcols = existing_forest_cols]
}

# Create deforestation indicator (any deforestation in cell-year)
if ("Deforested" %in% names(final_panel)) {
  final_panel[, has_deforestation := Deforested > 0]
}

# ==============================================================================
# QUALITY CONTROL CHECKS
# ==============================================================================

log_message("Running quality control checks...")

# Check 1: No duplicate grid_id-year combinations
n_dups <- nrow(final_panel) - nrow(unique(final_panel, by = c("grid_id", "year")))
if (n_dups > 0) {
  log_message(sprintf("WARNING: Found %d duplicate grid_id-year combinations", n_dups))
  final_panel <- unique(final_panel, by = c("grid_id", "year"))
}

# Check 2: Interior cells should have data in all years
if ("n_years_observed" %in% names(final_panel)) {
  interior_missing <- final_panel[is_interior == TRUE & n_years_observed < N_TMF_YEARS]
  if (nrow(interior_missing) > 0) {
    log_message(sprintf("WARNING: %d interior cells have < %d years of data",
                        length(unique(interior_missing$grid_id)), N_TMF_YEARS))
  }
}

# Check 3: Frontier should not overlap with interior
both_flags <- final_panel[is_interior == TRUE & is_frontier == TRUE]
if (nrow(both_flags) > 0) {
  log_message(sprintf("NOTE: %d observations are both interior and frontier (expected: 0)",
                      nrow(both_flags)))
}

# Check 4: Binary protection (cells are fully protected or not, cut on WDPA boundaries)
n_protected <- sum(final_panel$is_protected, na.rm = TRUE)
n_unique_protected <- length(unique(final_panel[is_protected == TRUE]$grid_id))
log_message(sprintf("Protected observations: %s (%s unique cells)",
                    format(n_protected, big.mark = ","),
                    format(n_unique_protected, big.mark = ",")))

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

log_message("Final dataset summary:")
log_message(sprintf("  Total observations: %s",
                    format(nrow(final_panel), big.mark = ",")))
log_message(sprintf("  Unique grid cells: %s",
                    format(length(unique(final_panel$grid_id)), big.mark = ",")))
log_message(sprintf("  Years: %d-%d (%d years)",
                    min(final_panel$year, na.rm = TRUE),
                    max(final_panel$year, na.rm = TRUE),
                    length(unique(final_panel$year))))
log_message(sprintf("  Countries: %d",
                    length(unique(final_panel$country_iso3))))

# Zone distribution
zone_dist <- final_panel[!is.na(year), .(
  n_cells = length(unique(grid_id)),
  n_observations = .N
), by = zone]

log_message("Zone distribution:")
for (i in seq_len(nrow(zone_dist))) {
  log_message(sprintf("  %s: %s cells, %s obs",
                      zone_dist$zone[i],
                      format(zone_dist$n_cells[i], big.mark = ","),
                      format(zone_dist$n_observations[i], big.mark = ",")))
}

# Protection status (cells are cut on WDPA boundaries, so each is fully protected or not)
prot_dist <- final_panel[!is.na(year), .(
  n_cells = length(unique(grid_id)),
  n_observations = .N
), by = is_protected]

log_message("Protection status:")
for (i in seq_len(nrow(prot_dist))) {
  log_message(sprintf("  Protected=%s: %s cells, %s obs",
                      prot_dist$is_protected[i],
                      format(prot_dist$n_cells[i], big.mark = ","),
                      format(prot_dist$n_observations[i], big.mark = ",")))
}

# Forest cover by zone
if ("forest_cover" %in% names(final_panel)) {
  fc_by_zone <- final_panel[!is.na(year), .(
    mean_forest_cover = mean(forest_cover, na.rm = TRUE)
  ), by = zone]

  log_message("Mean forest cover by zone:")
  for (i in seq_len(nrow(fc_by_zone))) {
    log_message(sprintf("  %s: %.3f", fc_by_zone$zone[i], fc_by_zone$mean_forest_cover[i]))
  }
}

# ==============================================================================
# SORT AND WRITE OUTPUT
# ==============================================================================

log_message("Sorting data...")

setorder(final_panel, grid_id, year)

log_message(sprintf("Writing output: %s", output_file))

write_atomic(final_panel, output_file)

# Verify write
if (file.exists(output_file)) {
  file_size <- file.info(output_file)$size / 1024 / 1024 / 1024  # GB
  log_message(sprintf("Output file size: %.2f GB", file_size))
}

# ==============================================================================
# WRITE SUMMARY FILE
# ==============================================================================

summary_file <- file.path(final_output_path, "summary_stats.txt")
log_message(sprintf("Writing summary: %s", summary_file))

sink(summary_file)
cat("TMF 5km Panel Dataset Summary\n")
cat("=============================\n")
cat(sprintf("Generated: %s\n\n", Sys.time()))

cat(sprintf("Total observations: %s\n", format(nrow(final_panel), big.mark = ",")))
cat(sprintf("Unique grid cells: %s\n", format(length(unique(final_panel$grid_id)), big.mark = ",")))
cat(sprintf("Years: %d-%d\n", min(final_panel$year, na.rm = TRUE), max(final_panel$year, na.rm = TRUE)))
cat(sprintf("Countries: %d\n\n", length(unique(final_panel$country_iso3))))

cat("Columns:\n")
for (col in names(final_panel)) {
  cat(sprintf("  - %s (%s)\n", col, class(final_panel[[col]])[1]))
}

sink()

gc_verbose()
log_job_end(start_time)

log_message("Pipeline complete!")
