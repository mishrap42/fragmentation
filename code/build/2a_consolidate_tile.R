# ==============================================================================
# STAGE 2A: CONSOLIDATE TMF DATA BY TILE
# ==============================================================================
# This script consolidates all years of TMF data for a single tile
# into a single parquet file.
#
# Input: SLURM_ARRAY_TASK_ID (1 to N_TMF_TILES)
# Output: Data/build/tmf_consolidated/tile_{tile_id}.parquet
# ==============================================================================

# Load configuration
source("Code/BUILD_workspace.R")

# Record start time
start_time <- Sys.time()

# Get task ID
task_id <- get_slurm_task_id()
tile_id <- task_id

log_job_start("2a_consolidate_tile.R", task_id)
log_message(sprintf("Consolidating tile: %d", tile_id))

# Validate
if (tile_id < 1 || tile_id > N_TMF_TILES) {
  stop(sprintf("Invalid tile_id: %d", tile_id))
}

# Check if output already exists
output_file <- get_consolidated_tile_filename(tile_id)
skip_if_exists(output_file, sprintf("tile %d", tile_id))

# ==============================================================================
# FIND INPUT FILES
# ==============================================================================

log_message("Finding input files...")

input_files <- sapply(TMF_YEARS, function(year) {
  get_tmf_filename(tile_id, year)
})

# Check which files exist
files_exist <- file.exists(input_files)
log_message(sprintf("Found %d/%d year files",
                    sum(files_exist), length(input_files)))

if (sum(files_exist) == 0) {
  log_message("No input files found for this tile. Creating empty output.")
  empty_output <- data.table(
    grid_id = character(0),
    year = integer(0),
    tmf_class = character(0),
    fraction = numeric(0)
  )
  write_atomic(empty_output, output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# ==============================================================================
# READ AND COMBINE DATA
# ==============================================================================

log_message("Reading and combining data...")

# Read all files into a list
data_list <- lapply(which(files_exist), function(i) {
  year <- TMF_YEARS[i]
  file_path <- input_files[i]

  # Check file size
  file_size <- file.info(file_path)$size

  if (file_size < 10) {
    # Empty or near-empty file
    return(NULL)
  }

  # Read data
  dt <- fread(file_path)

  if (nrow(dt) == 0) {
    return(NULL)
  }

  return(dt)
})

# Remove NULL entries
data_list <- data_list[!sapply(data_list, is.null)]

if (length(data_list) == 0) {
  log_message("All input files are empty. Creating empty output.")
  empty_output <- data.table(
    grid_id = character(0),
    year = integer(0),
    tmf_class = character(0),
    fraction = numeric(0)
  )
  write_atomic(empty_output, output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# Combine into single data.table
combined_data <- rbindlist(data_list, use.names = TRUE, fill = TRUE)

log_message(sprintf("Combined data: %d rows", nrow(combined_data)))

# Clean up to free memory
rm(data_list)
gc_verbose()

# ==============================================================================
# VALIDATE AND CLEAN DATA
# ==============================================================================

log_message("Validating and cleaning data...")

# Check for required columns
required_cols <- c("grid_id", "year", "tmf_class", "fraction")
missing_cols <- setdiff(required_cols, names(combined_data))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

# Remove any rows with NA in key columns
combined_data <- combined_data[!is.na(grid_id) & !is.na(year) & !is.na(tmf_class)]

# Remove zero or negative fractions
combined_data <- combined_data[fraction > 0]

# Sort for consistent output
setorder(combined_data, grid_id, year, tmf_class)

log_message(sprintf("After cleaning: %d rows", nrow(combined_data)))

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

log_message("Summary statistics:")
log_message(sprintf("  Unique grid cells: %d", length(unique(combined_data$grid_id))))
log_message(sprintf("  Years covered: %d-%d (%d years)",
                    min(combined_data$year),
                    max(combined_data$year),
                    length(unique(combined_data$year))))

# Class distribution
class_summary <- combined_data[, .(
  n_observations = .N,
  mean_fraction = mean(fraction)
), by = tmf_class]
setorder(class_summary, -n_observations)

for (i in seq_len(nrow(class_summary))) {
  log_message(sprintf("  %s: %d obs, mean=%.3f",
                      class_summary$tmf_class[i],
                      class_summary$n_observations[i],
                      class_summary$mean_fraction[i]))
}

# ==============================================================================
# WRITE OUTPUT
# ==============================================================================

log_message(sprintf("Writing output: %s", output_file))

write_atomic(combined_data, output_file)

# Verify write
if (file.exists(output_file)) {
  file_size <- file.info(output_file)$size / 1024 / 1024  # MB
  log_message(sprintf("Output file size: %.1f MB", file_size))
} else {
  stop("Failed to write output file!")
}

# Clean up
gc_verbose()
log_job_end(start_time)
