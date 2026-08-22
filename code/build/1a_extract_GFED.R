# ==============================================================================
# STAGE 1a: EXTRACT GFED5 BURNED AREA DATA
# ==============================================================================
# This script extracts GFED5 burned area for each grid cell using exact_extract.
# Runs at SUB-TILE level, processing all available GFED years per job.
#
# Input: SLURM_ARRAY_TASK_ID (1 to N_SUB_TILES = 344)
# Output: Data/build/gfed/gfed_sub_{sub_tile_id}.parquet
#
# PREREQUISITE: Run Stage 0 first to create grid files
# ==============================================================================

# Load configuration
here::i_am('code/build/1a_extract_GFED.R')
source("code/build/BUILD_workspace.R")

# Record start time
start_time <- Sys.time()

# Get task ID (= sub_tile_id)
task_id <- get_slurm_task_id()
sub_tile_id <- task_id

log_job_start("1a_extract_GFED.R", task_id)
log_message(sprintf("Processing sub-tile %d", sub_tile_id))

# Validate
if (sub_tile_id < 1 || sub_tile_id > N_SUB_TILES) {
  stop(sprintf("Invalid sub_tile_id: %d (must be 1-%d)", sub_tile_id, N_SUB_TILES))
}

# Check if output already exists
output_file <- get_gfed_subtile_filename(sub_tile_id)
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
    write_atomic(data.table(grid_id = character(0), year = integer(0),
                            month = integer(0), burned_area = numeric(0)),
                 output_file)
    log_job_end(start_time)
    quit(save = "no", status = 0)
  }
  stop(sprintf("Grid file not found: %s\nRun Stage 0 first.", grid_gpkg))
}

# Get available GFED years
GFED_YEARS <- get_gfed_years()
N_GFED_YEARS <- length(GFED_YEARS)
log_message(sprintf("GFED years available: %d-%d (%d years)",
                    min(GFED_YEARS), max(GFED_YEARS), N_GFED_YEARS))

# ==============================================================================
# LOAD GRID GEOMETRIES
# ==============================================================================

log_message("Loading grid geometries...")

grid_sf <- st_read(grid_gpkg, quiet = TRUE)
log_message(sprintf("Loaded %d grid cells", nrow(grid_sf)))

if (nrow(grid_sf) == 0) {
  log_message("Empty grid. Creating empty output.")
  empty_output <- data.table(
    grid_id = character(0),
    year = integer(0),
    month = integer(0),
    burned_area = numeric(0)
  )
  write_atomic(empty_output, output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# ==============================================================================
# NORMALIZE GEOMETRY (following TMF pattern from 1_extract_TMF.R lines 62-90)
# ==============================================================================

# Check and fix geometry types
geom_types <- as.character(unique(st_geometry_type(grid_sf)))
geom_counts <- table(st_geometry_type(grid_sf))
log_message(sprintf("Geometry types: %s",
                    paste(names(geom_counts), geom_counts, sep = "=", collapse = ", ")))

# Track split grid_ids for area-weighted aggregation
split_grid_ids <- character(0)

if (length(geom_types) > 1 || any(grepl("COLLECTION", geom_types))) {
  log_message("Normalizing to MULTIPOLYGON...")
  n_before <- nrow(grid_sf)
  grid_sf <- st_collection_extract(grid_sf, "POLYGON")
  grid_sf <- st_cast(grid_sf, "MULTIPOLYGON")
  n_after <- nrow(grid_sf)
  log_message(sprintf("After normalization: %d features (was %d)", n_after, n_before))

  # Identify split grid_ids
  id_counts <- table(grid_sf$grid_id)
  split_grid_ids <- names(id_counts[id_counts > 1])

  if (length(split_grid_ids) > 0) {
    log_message(sprintf("WARNING: %d grid_ids split into multiple polygons",
                        length(split_grid_ids)))
    log_message("Will compute area-weighted means for these cells")

    # Calculate area for split polygons
    split_idx <- grid_sf$grid_id %in% split_grid_ids
    grid_sf$split_area <- NA_real_
    grid_sf$split_area[split_idx] <- as.numeric(st_area(grid_sf[split_idx, ]))
  }
}

# Ensure WGS84 CRS (GFED is in WGS84)
grid_sf <- st_transform(grid_sf, WGS84_CRS)

# ==============================================================================
# EXTRACT GFED DATA FOR ALL YEARS AND MONTHS
# ==============================================================================

log_message("Extracting GFED burned area...")

# Determine append columns
append_cols <- "grid_id"
has_split_area <- "split_area" %in% names(grid_sf)
if (has_split_area) {
  append_cols <- c("grid_id", "split_area")
}

results_list <- list()

for (year in GFED_YEARS) {
  gfed_file <- get_gfed_file(year)

  if (!file.exists(gfed_file)) {
    stop(sprintf("GFED file not found for year %d: %s", year, gfed_file))
  }

  log_message(sprintf("  Processing year %d...", year))

  # Process each month
  for (month in 1:12) {
    # Convert NetCDF to raster for this month
    gfed_rast <- gfed_to_raster(gfed_file, month = month)

    # Extract mean burned fraction using exact_extract
    extraction <- exact_extract(
      gfed_rast,
      grid_sf,
      fun = "mean",
      append_cols = append_cols,
      progress = FALSE
    )

    setDT(extraction)
    extraction[, year := year]
    extraction[, month := month]

    # Rename mean column
    setnames(extraction, "mean", "burned_area")

    results_list[[length(results_list) + 1]] <- extraction

    rm(gfed_rast)
  }

  gc()
}

# Combine all results
output_data <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
rm(results_list)

log_message(sprintf("Raw extraction: %d cell-year-month combinations", nrow(output_data)))

# ==============================================================================
# HANDLE SPLIT POLYGONS (area-weighted aggregation)
# ==============================================================================

if (has_split_area && length(split_grid_ids) > 0) {
  log_message("Computing area-weighted means for split grid_ids...")

  split_rows <- output_data[grid_id %in% split_grid_ids]
  nonsplit_rows <- output_data[!(grid_id %in% split_grid_ids)]

  if (nrow(split_rows) > 0) {
    split_agg <- split_rows[, .(
      burned_area = weighted.mean(burned_area, split_area, na.rm = TRUE)
    ), by = .(grid_id, year, month)]

    log_message(sprintf("  Split rows: %d -> %d after aggregation",
                        nrow(split_rows), nrow(split_agg)))
  } else {
    split_agg <- data.table(grid_id = character(0), year = integer(0),
                            month = integer(0), burned_area = numeric(0))
  }

  nonsplit_out <- nonsplit_rows[, .(grid_id, year, month, burned_area)]
  output_data <- rbind(nonsplit_out, split_agg)
}

# Select final columns
output_data <- output_data[, .(grid_id, year, month, burned_area)]

# ==============================================================================
# FILTER AND CLEAN
# ==============================================================================

# Check for NA values
n_na <- sum(is.na(output_data$burned_area))
if (n_na > 0) {
  na_grid_ids <- unique(output_data[is.na(burned_area)]$grid_id)
  log_message(sprintf("Removing %d NA values (%d grid_ids outside GFED coverage)",
                      n_na, length(na_grid_ids)))
}
output_data <- output_data[!is.na(burned_area)]

# Remove exact zeros to save space (no burning detected)
n_zeros <- sum(output_data$burned_area == 0)
output_data <- output_data[burned_area > 0]
log_message(sprintf("Removed %d zero values", n_zeros))

# Sort for consistent output
setorder(output_data, grid_id, year, month)

log_message(sprintf("Final output: %d observations with burned_area > 0", nrow(output_data)))

# ==============================================================================
# VERIFY DATA QUALITY
# ==============================================================================

if (nrow(output_data) > 0) {
  # Check for duplicates (should not exist after aggregation)
  dup_check <- output_data[, .N, by = .(grid_id, year, month)][N > 1]
  if (nrow(dup_check) > 0) {
    stop(sprintf("Found %d unexpected duplicates after aggregation!", nrow(dup_check)))
  }

  # Check burned_area is in valid range [0, 1]
  invalid_range <- output_data[burned_area < 0 | burned_area > 1]
  if (nrow(invalid_range) > 0) {
    log_message(sprintf("WARNING: %d values outside [0,1] range (max: %.4f)",
                        nrow(invalid_range), max(output_data$burned_area)))
  }
}

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

if (nrow(output_data) > 0) {
  log_message("Summary statistics:")
  log_message(sprintf("  Unique grid cells with burns: %d",
                      length(unique(output_data$grid_id))))
  log_message(sprintf("  Years with data: %d-%d",
                      min(output_data$year), max(output_data$year)))
  log_message(sprintf("  Mean burned fraction: %.6f", mean(output_data$burned_area)))
  log_message(sprintf("  Max burned fraction: %.6f", max(output_data$burned_area)))
} else {
  log_message("No burned area detected in this sub-tile")
}

# ==============================================================================
# WRITE OUTPUT
# ==============================================================================

log_message(sprintf("Writing output: %s", output_file))
write_atomic(output_data, output_file)

if (file.exists(output_file)) {
  file_size <- file.info(output_file)$size / 1024
  log_message(sprintf("Output file size: %.1f KB", file_size))
}

gc_verbose()
log_job_end(start_time)
