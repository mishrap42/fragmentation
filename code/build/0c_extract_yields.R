# ==============================================================================
# STAGE 0c: EXTRACT AGRICULTURAL YIELD RASTERS
# ==============================================================================
# This script extracts mean values from agricultural yield rasters for each
# grid cell. Runs after grid creation (Stage 0) and before TMF extraction.
#
# Extracts three things, all EarthStat, all a year-2000 cross-section:
#   yield_<crop>     tons/ha on the crop's harvested area   (172 crops)
#   cropshare_<crop> harvested area as a fraction of CELL   (172 crops)
#   cropland_frac / pasture_frac   cell fractions, the denominators the
#                    per-crop layers cannot supply on their own
#
# cropshare_* is a fraction of the CELL, not of its cropland, and the 172 of
# them can sum past cropland_frac because multi-cropping counts a hectare once
# per harvest. Normalizing is an analysis decision, not an extraction one, so
# nothing is normalized here - both the numerators and the denominator are
# written raw and the choice is made downstream.
#
# Input: SLURM_ARRAY_TASK_ID (1 to N_SUB_TILES)
# Output: Data/build/yields/yields_sub_{sub_tile_id}.parquet
#
# PREREQUISITE: Run Stage 0 first to create grid files
# ==============================================================================

# Load configuration
here::i_am('code/build/0c_extract_yields.R')
source("code/build/BUILD_workspace.R")

# Record start time
start_time <- Sys.time()

# Get task ID (= sub_tile_id)
task_id <- get_slurm_task_id()
sub_tile_id <- task_id

log_job_start("0c_extract_yields.R", task_id)
log_message(sprintf("Processing sub-tile %d", sub_tile_id))

# Check if output already exists
output_file <- get_yields_filename(sub_tile_id)
skip_if_exists(output_file, sprintf("sub-tile %d", sub_tile_id))

# ==============================================================================
# CHECK PREREQUISITES
# ==============================================================================

# Check grid file exists
# A missing .gpkg is NOT necessarily a missing stage 0. Stage 0 writes an empty
# .parquet and NO .gpkg for sub-tiles with no land (40 of 344 - open ocean).
# stop()ing on those turns an expected state into a task failure, which is what
# left TMF_CROP permanently PD with DependencyNeverSatisfied: afterok can never
# be satisfied when 40 tasks "fail" by design.
#
# Distinguish the two cases by the .parquet: present means stage 0 ran and the
# sub-tile is genuinely empty; absent means stage 0 really has not run.
grid_gpkg    <- get_grid_filename(sub_tile_id, "gpkg")
grid_parquet <- get_grid_filename(sub_tile_id, "parquet")
if (!file.exists(grid_gpkg)) {
  if (file.exists(grid_parquet)) {
    log_message(sprintf(
      "Sub-tile %d has no grid geometry (ocean; stage 0 wrote an empty grid). Creating empty output.",
      sub_tile_id))
    write_atomic(data.table(grid_id = character(0)), output_file)
    log_job_end(start_time)
    quit(save = "no", status = 0)
  }
  stop(sprintf("Grid file not found: %s\nRun Stage 0 first.", grid_gpkg))
}

# ==============================================================================
# LIST YIELD RASTERS
# ==============================================================================

yield_rasters   <- list_yield_rasters()
landuse_rasters <- list_landuse_rasters()

# Bail only if BOTH sources are missing. Guarding on yield_rasters alone would
# throw away cropland_frac / pasture_frac whenever the 172-crop tree is
# unreachable, which is exactly the mount failure this stage has hit before.
if (length(yield_rasters) == 0 && length(landuse_rasters) == 0) {
  log_message("No yield or land-use rasters found. Creating empty output.")
  empty_output <- data.table(grid_id = character(0))
  write_atomic(empty_output, output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

log_message(sprintf("Found %d crop rasters and %d land-use rasters to extract",
                    length(yield_rasters), length(landuse_rasters)))

# ==============================================================================
# LOAD GRID GEOMETRIES
# ==============================================================================

log_message("Loading grid geometries...")

grid_sf <- st_read(grid_gpkg, quiet = TRUE)
log_message(sprintf("Loaded %d grid cells", nrow(grid_sf)))

if (nrow(grid_sf) == 0) {
  log_message("Empty grid. Creating empty output.")
  empty_output <- data.table(grid_id = character(0))
  write_atomic(empty_output, output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# Ensure WGS84 CRS for raster extraction
grid_sf <- st_transform(grid_sf, WGS84_CRS)

# ==============================================================================
# EXTRACT YIELD VALUES
# ==============================================================================

# Initialize result with grid_id
yield_data <- data.table(grid_id = grid_sf$grid_id)

#' Coverage-weighted mean of one raster onto the sub-tile grid.
#'
#' exact_extract's "mean" weights each source pixel by the share of it falling
#' inside the cell and skips NA. For EarthStat that is the behaviour we want:
#' ocean is NA, so a half-ocean coastal cell reports the mean over its LAND
#' pixels - i.e. a fraction of land area, not of land-plus-water. Cells that are
#' entirely NA come back NaN; coerced to NA here so the parquet does not carry
#' both sentinels for the same condition.
#'
#' A raster that cannot be read costs its own column, not the sub-tile.
extract_mean_layer <- function(raster_path, varname) {
  tryCatch({
    r <- terra::rast(raster_path)
    values <- exact_extract(r, grid_sf, "mean")
    values[is.nan(values)] <- NA_real_
    yield_data[, (varname) := values]
  }, error = function(e) {
    log_message(sprintf("  WARNING: Failed to extract %s: %s", varname, e$message))
    yield_data[, (varname) := NA_real_]
  })
  gc()
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# Per-crop layers: yield_<crop> and cropshare_<crop>
# ------------------------------------------------------------------------------

log_message("Extracting per-crop yield / harvested-area-fraction values...")

n_rasters <- length(yield_rasters)
for (i in seq_along(yield_rasters)) {
  raster_path <- yield_rasters[i]
  varname <- get_yield_varname(raster_path)

  # Progress message every 10 rasters
  if (i %% 10 == 1 || i == n_rasters) {
    log_message(sprintf("  Extracting raster %d/%d: %s", i, n_rasters, varname))
  }

  extract_mean_layer(raster_path, varname)
}

# ------------------------------------------------------------------------------
# Cell-level land-use fractions: cropland_frac, pasture_frac
# ------------------------------------------------------------------------------
# Extracted here rather than in Stage 5 so that every EarthStat layer shares one
# provenance, one CRS handling path, and one validator. cropland_frac is the
# denominator for a "share of cropland" reading of cropshare_*; pasture_frac is
# the only pasture signal in the build (GLW3 gives livestock HEADS, not land).

if (length(landuse_rasters) > 0) {
  log_message("Extracting EarthStat cropland / pasture fractions...")

  for (varname in names(landuse_rasters)) {
    log_message(sprintf("  Extracting %s", varname))
    extract_mean_layer(landuse_rasters[[varname]], varname)
  }
} else {
  log_message("No EarthStat land-use rasters available; skipping.")
}

# ------------------------------------------------------------------------------
# Total harvested fraction across all crops
# ------------------------------------------------------------------------------
# Written here because it is cheap on a sub-tile and expensive on the assembled
# cross-section. NOTE it can exceed cropland_frac: EarthStat counts a hectare
# once per harvest, so double-cropped land is double-counted. That ratio is the
# multi-cropping intensity, not an error - 6b reports it so it stays visible.

cropshare_cols <- grep("^cropshare_", names(yield_data), value = TRUE)
if (length(cropshare_cols) > 0) {
  yield_data[, cropshare_total := rowSums(.SD, na.rm = TRUE),
             .SDcols = cropshare_cols]
  # All-NA rows sum to 0 under na.rm; that is "no crops", which is correct for
  # a harvested-area fraction and NOT the same as "unobserved".
  log_message(sprintf("  cropshare_total: %d crops summed, mean %.5f",
                      length(cropshare_cols),
                      mean(yield_data$cropshare_total, na.rm = TRUE)))
}

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

log_message("Summary statistics:")
log_message(sprintf("  Grid cells: %d", nrow(yield_data)))
log_message(sprintf("  Total variables: %d", ncol(yield_data) - 1))
log_message(sprintf("    yield_*:     %d", sum(grepl("^yield_", names(yield_data)))))
log_message(sprintf("    cropshare_*: %d", sum(grepl("^cropshare_", names(yield_data)))))
log_message(sprintf("    land-use:    %d",
                    sum(c("cropland_frac", "pasture_frac") %in% names(yield_data))))

# Show summary of a few variables
yield_cols <- setdiff(names(yield_data), "grid_id")
for (col in head(yield_cols, 5)) {
  vals <- yield_data[[col]]
  non_na <- sum(!is.na(vals))
  if (non_na > 0) {
    log_message(sprintf("    %s: mean=%.2f, non-NA=%d/%d",
                        col, mean(vals, na.rm = TRUE), non_na, length(vals)))
  }
}
if (length(yield_cols) > 5) {
  log_message(sprintf("    ... and %d more variables", length(yield_cols) - 5))
}

# ==============================================================================
# WRITE OUTPUT
# ==============================================================================

log_message(sprintf("Writing output: %s", output_file))
write_atomic(yield_data, output_file)

if (file.exists(output_file)) {
  file_size <- file.info(output_file)$size / 1024
  log_message(sprintf("Output file size: %.1f KB", file_size))
}

gc_verbose()
log_job_end(start_time)
