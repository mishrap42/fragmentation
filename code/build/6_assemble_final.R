# ==============================================================================
# STAGE 6: ASSEMBLE FINAL DATASET
# ==============================================================================
# This script merges all pipeline outputs into the final panel dataset.
#
# Input: Grid metadata, TMF time series, classifications, WDPA data
# Output: Data/build/final/TMF_1km_panel.parquet
# ==============================================================================

# Load configuration
source("Code/BUILD_workspace.R")

# Record start time
start_time <- Sys.time()

log_job_start("6_assemble_final.R", task_id = 1)

# Check if output already exists
output_file <- file.path(final_output_path, "TMF_1km_panel.parquet")
skip_if_exists(output_file, "final assembly")

# ==============================================================================
# LOAD ALL GRID METADATA
# ==============================================================================

log_message("Loading grid metadata...")

grid_files <- sapply(1:N_TMF_TILES, function(t) get_grid_filename(t, "parquet"))
files_exist <- file.exists(grid_files)

log_message(sprintf("Found %d/%d grid files", sum(files_exist), length(grid_files)))

grid_list <- lapply(grid_files[files_exist], function(f) {
  dt <- read_parquet(f)
  setDT(dt)
  return(dt)
})

grid_all <- rbindlist(grid_list, fill = TRUE)
log_message(sprintf("Total grid cells: %s",
                    format(nrow(grid_all), big.mark = ",")))

rm(grid_list)
gc()

# ==============================================================================
# LOAD TMF TIME SERIES
# ==============================================================================

log_message("Loading TMF time series...")

tmf_file <- file.path(consolidated_path, "tmf_all.parquet")

if (!file.exists(tmf_file)) {
  stop(sprintf("TMF consolidated file not found: %s. Run Stage 2b first.", tmf_file))
}

tmf_all <- read_parquet(tmf_file)
setDT(tmf_all)

log_message(sprintf("TMF observations: %s",
                    format(nrow(tmf_all), big.mark = ",")))

# Reshape to wide format for easier analysis
log_message("Reshaping TMF to wide format...")

tmf_wide <- dcast(tmf_all, grid_id + year ~ tmf_class,
                  value.var = "fraction", fill = 0)

log_message(sprintf("TMF wide format: %s cell-years",
                    format(nrow(tmf_wide), big.mark = ",")))

rm(tmf_all)
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

# ==============================================================================
# LOAD WDPA DATA
# ==============================================================================

log_message("Loading WDPA data...")

wdpa_file <- get_wdpa_filename(NULL)  # Consolidated

if (!file.exists(wdpa_file)) {
  log_message("WDPA consolidated file not found, loading from tiles...")

  wdpa_files <- sapply(1:N_TMF_TILES, get_wdpa_filename)
  files_exist <- file.exists(wdpa_files)

  wdpa_list <- lapply(wdpa_files[files_exist], function(f) {
    dt <- read_parquet(f)
    setDT(dt)
    return(dt)
  })

  wdpa_all <- rbindlist(wdpa_list, fill = TRUE)
  wdpa_all <- unique(wdpa_all, by = "grid_id")

  rm(wdpa_list)
} else {
  wdpa_all <- read_parquet(wdpa_file)
  setDT(wdpa_all)
}

log_message(sprintf("WDPA data: %s cells (%s protected)",
                    format(nrow(wdpa_all), big.mark = ","),
                    format(sum(wdpa_all$is_protected, na.rm = TRUE), big.mark = ",")))

gc()

# ==============================================================================
# MERGE CLASSIFICATIONS
# ==============================================================================

log_message("Merging grid metadata with classifications...")

# Start with grid metadata
panel_base <- grid_all[, .(grid_id, tile_id, country_iso3, country_name,
                           centroid_lon, centroid_lat, area_km2)]

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

# Add WDPA data
panel_base <- merge(
  panel_base,
  wdpa_all[, .(grid_id, is_protected, protected_frac,
               wdpa_pid, iucn_cat, desig_year, gov_type)],
  by = "grid_id",
  all.x = TRUE
)

# Fill NAs
panel_base[is.na(is_protected), is_protected := FALSE]
panel_base[is.na(protected_frac), protected_frac := 0]

log_message(sprintf("Base panel: %s cells",
                    format(nrow(panel_base), big.mark = ",")))

# Clean up
rm(grid_all, interior_all, frontier_all, wdpa_all)
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

log_message(sprintf("Final panel: %s observations",
                    format(nrow(final_panel), big.mark = ",")))

# Clean up
rm(panel_base, tmf_wide)
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

# Check 4: Protected fraction should be <= 1
if ("protected_frac" %in% names(final_panel)) {
  invalid_frac <- final_panel[protected_frac > 1]
  if (nrow(invalid_frac) > 0) {
    log_message(sprintf("WARNING: %d cells have protected_frac > 1, capping at 1",
                        length(unique(invalid_frac$grid_id))))
    final_panel[protected_frac > 1, protected_frac := 1]
  }
}

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

# Protection status
prot_dist <- final_panel[!is.na(year), .(
  n_cells = length(unique(grid_id)),
  mean_protected_frac = mean(protected_frac, na.rm = TRUE)
), by = is_protected]

log_message("Protection status:")
for (i in seq_len(nrow(prot_dist))) {
  log_message(sprintf("  Protected=%s: %s cells, mean fraction=%.3f",
                      prot_dist$is_protected[i],
                      format(prot_dist$n_cells[i], big.mark = ","),
                      prot_dist$mean_protected_frac[i]))
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
cat("TMF 1km Panel Dataset Summary\n")
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
