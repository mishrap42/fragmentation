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
source("Code/BUILD_workspace.R")

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
# LOAD GRID
# ==============================================================================

log_message("Loading grid cells...")

grid_file <- get_grid_filename(tile_id, "gpkg")

if (!file.exists(grid_file)) {
  stop(sprintf("Grid file not found: %s. Run Stage 0 first.", grid_file))
}

grid_sf <- st_read(grid_file, quiet = TRUE)
log_message(sprintf("Loaded %d grid cells", nrow(grid_sf)))

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

# TMF rasters are in WGS84
raster_crs <- terra::crs(tmf_raster, proj = TRUE)
grid_reproj <- st_transform(grid_sf, raster_crs)

# ==============================================================================
# EXTRACT TMF DATA
# ==============================================================================

log_message("Extracting TMF data using exactextractr...")

# Use exactextractr for coverage-weighted extraction
# This handles partial cell coverage accurately
extraction_result <- exact_extract(
  tmf_raster,
  grid_reproj,
  fun = function(df) {
    # df has columns: value, coverage_fraction
    setDT(df)

    # Filter out NA values
    df <- df[!is.na(value)]

    if (nrow(df) == 0) {
      return(data.table(value = integer(0), fraction = numeric(0)))
    }

    # Calculate total coverage fraction per value
    # Normalize by total coverage to get proportions within the cell
    total_coverage <- sum(df$coverage_fraction, na.rm = TRUE)

    if (total_coverage == 0) {
      return(data.table(value = integer(0), fraction = numeric(0)))
    }

    result <- df[, .(
      fraction = sum(coverage_fraction, na.rm = TRUE) / total_coverage
    ), by = value]

    return(result)
  },
  summarize_df = TRUE,
  include_cols = "grid_id",
  progress = TRUE
)

log_message(sprintf("Extraction complete: %d rows", nrow(extraction_result)))

# ==============================================================================
# PROCESS RESULTS
# ==============================================================================

log_message("Processing results...")

setDT(extraction_result)

# Filter out empty results
extraction_result <- extraction_result[!is.na(value) & fraction > 0]

# Add year column
extraction_result[, year := year]

# Map TMF values to class names
extraction_result <- merge(extraction_result, TMF_LEGEND,
                           by = "value", all.x = TRUE)

# Handle any unmapped values
extraction_result[is.na(tmf_class), tmf_class := paste0("Unknown_", value)]

# Select and order columns
output_data <- extraction_result[, .(grid_id, year, tmf_class, fraction)]

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
