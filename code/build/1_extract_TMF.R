# ==============================================================================
# STAGE 1: EXTRACT TMF DATA
# ==============================================================================
# This script extracts TMF land cover fractions for each grid cell
# for a specific tile-year combination.
#
# Input: SLURM_ARRAY_TASK_ID (1 to N_EXTRACTION_JOBS = tiles * years)
# Output: Data/build/tmf/tmf_{tile_id}_{year}.csv.gz
# ==============================================================================

# Load configuration
here::i_am('code/build/1_extract_TMF.R')
source("code/build/BUILD_workspace.R")

# Record start time
start_time <- Sys.time()

# Get task ID and convert to tile/year
task_id <- get_slurm_task_id()
job_info <- job_to_year_tile(task_id)
tile_id <- job_info$tile_id
year <- job_info$year

log_job_start("1_extract_TMF.R", task_id)
log_message(sprintf("Tile: %d, Year: %d", tile_id, year))

# Validate
if (tile_id < 1 || tile_id > N_TMF_TILES) {
  stop(sprintf("Invalid tile_id: %d", tile_id))
}
if (!(year %in% TMF_YEARS)) {
  stop(sprintf("Invalid year: %d", year))
}

# Check if output already exists
output_file <- get_tmf_filename(tile_id, year)
skip_if_exists(output_file, sprintf("tile %d year %d", tile_id, year))

# ==============================================================================
# LOAD GRID (from sub-tile files within this TMF tile)
# ==============================================================================

log_message("Loading grid cells for tile...")

# Get all sub-tile grid files for this TMF tile
grid_files <- get_grid_files_for_tmf_tile(tile_id, "gpkg")
grid_files <- grid_files[file.exists(grid_files)]

if (length(grid_files) == 0) {
  stop(sprintf("No grid files found for TMF tile %d. Run Stage 0 first.", tile_id))
}

log_message(sprintf("Found %d sub-tile grid files", length(grid_files)))

# Load and combine all sub-tile grids
grid_list <- lapply(grid_files, function(f) {
  st_read(f, quiet = TRUE)
})
grid_sf <- do.call(rbind, grid_list)
rm(grid_list)

# Check and fix geometry types
geom_types <- as.character(unique(st_geometry_type(grid_sf)))
geom_counts <- table(st_geometry_type(grid_sf))
log_message(sprintf("Geometry types: %s", paste(names(geom_counts), geom_counts, sep = "=", collapse = ", ")))

# Normalize if mixed types or geometry collections present
split_grid_ids <- character(0)  # Track which grid_ids got split
if (length(geom_types) > 1 || any(grepl("COLLECTION", geom_types))) {
  log_message("Normalizing to MULTIPOLYGON...")
  n_before <- nrow(grid_sf)
  grid_sf <- st_collection_extract(grid_sf, "POLYGON")
  grid_sf <- st_cast(grid_sf, "MULTIPOLYGON")
  n_after <- nrow(grid_sf)
  log_message(sprintf("After normalization: %d features (was %d)", n_after, n_before))

  # Identify which grid_ids got split (have duplicates)
  id_counts <- table(grid_sf$grid_id)
  split_grid_ids <- names(id_counts[id_counts > 1])

  if (length(split_grid_ids) > 0) {
    log_message(sprintf("WARNING: %d grid_ids split into multiple polygons", length(split_grid_ids)))
    log_message("Will compute area-weighted means for these cells")

    # Calculate area for split polygons (needed for weighted aggregation)
    split_idx <- grid_sf$grid_id %in% split_grid_ids
    grid_sf$split_area <- NA_real_
    grid_sf$split_area[split_idx] <- as.numeric(st_area(grid_sf[split_idx, ]))
  }
}

log_message(sprintf("Loaded %d grid cells from %d sub-tiles",
                    nrow(grid_sf), length(grid_files)))

if (nrow(grid_sf) == 0) {
  log_message("Empty grid (likely ocean tile). Creating empty output.")
  empty_output <- data.table(
    grid_id = character(0),
    year = integer(0),
    tmf_class = character(0),
    fraction = numeric(0)
  )
  fwrite(empty_output, output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# ==============================================================================
# LOAD TMF RASTER
# ==============================================================================

log_message("Loading TMF raster...")

tmf_raster_file <- get_tmf_raster_path(tile_id, year)

if (!file.exists(tmf_raster_file)) {
  log_message(sprintf("TMF raster not found: %s", tmf_raster_file))
  log_message("This tile may not have TMF coverage for this year.")
  log_message("Creating empty output.")

  empty_output <- data.table(
    grid_id = character(0),
    year = integer(0),
    tmf_class = character(0),
    fraction = numeric(0)
  )
  fwrite(empty_output, output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

tmf_raster <- terra::rast(tmf_raster_file)
log_message(sprintf("TMF raster loaded: %d x %d pixels",
                    ncol(tmf_raster), nrow(tmf_raster)))
log_message(sprintf("Resolution: %.1f m x %.1f m",
                    terra::res(tmf_raster)[1],
                    terra::res(tmf_raster)[2]))

# ==============================================================================
# REPROJECT GRID TO MATCH RASTER
# ==============================================================================

log_message("Reprojecting grid to match raster CRS...")

raster_crs <- terra::crs(tmf_raster, proj = TRUE)
grid_reproj <- st_transform(grid_sf, raster_crs)
rm(grid_sf)

# ==============================================================================
# EXTRACT TMF DATA (using built-in frac function for speed)
# ==============================================================================

log_message(sprintf("Extracting TMF data for %d cells...", nrow(grid_reproj)))

# Determine which columns to append (include split_area if present)
append_cols <- "grid_id"
if ("split_area" %in% names(grid_reproj)) {
  append_cols <- c("grid_id", "split_area")
}

# Use built-in "frac" function - much faster than custom R function
# Returns wide format: grid_id, frac_0, frac_1, frac_2, ...
extraction_wide <- exact_extract(
  tmf_raster,
  grid_reproj,
  fun = "frac",
  append_cols = append_cols,
  progress = TRUE
)

rm(grid_reproj)
log_message(sprintf("Extraction complete: %d rows, %d columns",
                    nrow(extraction_wide), ncol(extraction_wide)))

# ==============================================================================
# PROCESS RESULTS (reshape wide to long)
# ==============================================================================

log_message("Processing results...")

setDT(extraction_wide)

# Get the frac columns (frac_0, frac_1, etc.)
frac_cols <- grep("^frac_", names(extraction_wide), value = TRUE)
log_message(sprintf("Found %d TMF classes in data", length(frac_cols)))

# Determine id.vars for melt (include split_area if present)
has_split_area <- "split_area" %in% names(extraction_wide)
id_vars <- if (has_split_area) c("grid_id", "split_area") else "grid_id"

# Reshape from wide to long
extraction_long <- melt(
  extraction_wide,
  id.vars = id_vars,
  measure.vars = frac_cols,
  variable.name = "frac_col",
  value.name = "fraction"
)

rm(extraction_wide)

# Extract the TMF value from column name (frac_1 -> 1)
extraction_long[, value := as.integer(sub("frac_", "", frac_col))]
extraction_long[, frac_col := NULL]

# Filter out zero/NA fractions
extraction_long <- extraction_long[!is.na(fraction) & fraction > 0]

log_message(sprintf("After filtering: %d non-zero cell-class combinations",
                    nrow(extraction_long)))

# Add year column
extraction_long[, year := year]

# Map TMF values to class names
extraction_long <- merge(extraction_long, TMF_LEGEND,
                         by = "value", all.x = TRUE)

# Handle any unmapped values
extraction_long[is.na(tmf_class), tmf_class := paste0("Unknown_", value)]

# ==============================================================================
# AGGREGATE SPLIT POLYGONS (area-weighted means)
# ==============================================================================

# For grid_ids that were split by st_collection_extract, we need to compute
# area-weighted means. Non-split grid_ids pass through unchanged.

n_before_agg <- nrow(extraction_long)

if (has_split_area && length(split_grid_ids) > 0) {
  log_message(sprintf("Computing area-weighted means for %d split grid_ids...",
                      length(split_grid_ids)))

  # Separate split and non-split rows
  split_rows <- extraction_long[grid_id %in% split_grid_ids]
  nonsplit_rows <- extraction_long[!(grid_id %in% split_grid_ids)]

  # For split rows: compute area-weighted mean
  if (nrow(split_rows) > 0) {
    split_agg <- split_rows[, .(
      fraction = weighted.mean(fraction, split_area, na.rm = TRUE)
    ), by = .(grid_id, year, tmf_class)]

    log_message(sprintf("  Split rows: %d -> %d after aggregation",
                        nrow(split_rows), nrow(split_agg)))
  } else {
    split_agg <- data.table(grid_id = character(0), year = integer(0),
                            tmf_class = character(0), fraction = numeric(0))
  }

  # For non-split rows: just select columns (no aggregation needed)
  nonsplit_out <- nonsplit_rows[, .(grid_id, year, tmf_class, fraction)]

  # Combine
  output_data <- rbind(nonsplit_out, split_agg)

} else {
  # No split polygons - just select columns
  output_data <- extraction_long[, .(grid_id, year, tmf_class, fraction)]

  # Still check for any unexpected duplicates and aggregate with simple mean
  dup_check <- output_data[, .N, by = .(grid_id, year, tmf_class)][N > 1]
  if (nrow(dup_check) > 0) {
    log_message(sprintf("WARNING: Found %d unexpected duplicates, using simple mean",
                        nrow(dup_check)))
    output_data <- output_data[, .(fraction = mean(fraction)), by = .(grid_id, year, tmf_class)]
  }
}

n_after_agg <- nrow(output_data)

if (n_before_agg != n_after_agg) {
  log_message(sprintf("Aggregated %d duplicate (grid_id, tmf_class) combinations",
                      n_before_agg - n_after_agg))
}

# Sort for consistent output
setorder(output_data, grid_id, tmf_class)

log_message(sprintf("Output data: %d rows for %d unique cells",
                    nrow(output_data),
                    length(unique(output_data$grid_id))))

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

log_message("Summary statistics:")

# Summarize by TMF class
class_summary <- output_data[, .(
  n_cells = length(unique(grid_id)),
  mean_fraction = mean(fraction),
  total_fraction = sum(fraction)
), by = tmf_class]

for (i in seq_len(nrow(class_summary))) {
  log_message(sprintf("  %s: %d cells, mean fraction %.3f",
                      class_summary$tmf_class[i],
                      class_summary$n_cells[i],
                      class_summary$mean_fraction[i]))
}

# ==============================================================================
# WRITE OUTPUT
# ==============================================================================

log_message(sprintf("Writing output: %s", output_file))

# Write as compressed CSV (gzip)
fwrite(output_data, output_file)

# Verify write
if (file.exists(output_file)) {
  file_size <- file.info(output_file)$size / 1024  # KB
  log_message(sprintf("Output file size: %.1f KB", file_size))
} else {
  stop("Failed to write output file!")
}

# Clean up
gc_verbose()
log_job_end(start_time)
