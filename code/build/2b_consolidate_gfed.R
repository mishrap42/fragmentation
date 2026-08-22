# ==============================================================================
# STAGE 2B: CONSOLIDATE GFED DATA BY TILE
# ==============================================================================
# This script consolidates all sub-tile GFED extractions for a single TMF tile.
# Aggregates monthly burned area to annual sums.
#
# Input: SLURM_ARRAY_TASK_ID (1 to N_TMF_TILES = 86)
# Output: Data/build/gfed_consolidated/gfed_tile_{tile_id}.parquet
#
# PREREQUISITE: Run Stage 1a first to extract GFED data
# ==============================================================================

# Load configuration
here::i_am("code/build/2b_consolidate_gfed.R")
source("code/build/BUILD_workspace.R")

# Record start time
start_time <- Sys.time()

# Get task ID
task_id <- get_slurm_task_id()
tile_id <- task_id

log_job_start("2b_consolidate_gfed.R", task_id)
log_message(sprintf("Consolidating GFED for tile: %d", tile_id))

# Validate
if (tile_id < 1 || tile_id > N_TMF_TILES) {

  stop(sprintf("Invalid tile_id: %d (must be 1-%d)", tile_id, N_TMF_TILES))
}

# Check if output already exists
output_file <- get_gfed_consolidated_filename(tile_id)
skip_if_exists(output_file, sprintf("tile %d", tile_id))

# ==============================================================================
# FIND INPUT FILES (sub-tile files for this TMF tile)
# ==============================================================================

log_message("Finding GFED sub-tile files...")

# Get sub-tiles for this TMF tile
sub_tile_ids <- get_sub_tiles_for_tmf_tile(tile_id)
log_message(sprintf("TMF tile %d contains %d sub-tiles", tile_id, length(sub_tile_ids)))

# Get file paths
input_files <- sapply(sub_tile_ids, get_gfed_subtile_filename)
files_exist <- file.exists(input_files)

log_message(sprintf("Found %d/%d sub-tile GFED files",
                    sum(files_exist), length(input_files)))

if (sum(files_exist) == 0) {
  log_message("No GFED files found for this tile. Creating empty output.")
  empty_output <- data.table(
    grid_id = character(0),
    year = integer(0),
    burned_annual = numeric(0)
  )
  write_atomic(empty_output, output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# ==============================================================================
# READ AND COMBINE SUB-TILE DATA
# ==============================================================================

log_message("Reading and combining data...")

data_list <- lapply(input_files[files_exist], function(f) {
  # Check file size
  file_size <- file.info(f)$size

  if (file_size < 100) {
    # Empty or near-empty file
    return(NULL)
  }

  dt <- read_parquet(f)
  setDT(dt)
  return(dt)
})

# Remove NULL entries
data_list <- data_list[!sapply(data_list, is.null)]
data_list <- data_list[sapply(data_list, function(x) nrow(x) > 0)]

if (length(data_list) == 0) {
  log_message("All input files are empty. Creating empty output.")
  empty_output <- data.table(
    grid_id = character(0),
    year = integer(0),
    burned_annual = numeric(0)
  )
  write_atomic(empty_output, output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

combined_data <- rbindlist(data_list, use.names = TRUE, fill = TRUE)
log_message(sprintf("Combined monthly data: %d rows", nrow(combined_data)))

rm(data_list)
gc()

# ==============================================================================
# VALIDATE INPUT DATA
# ==============================================================================

log_message("Validating input data...")

# Check for required columns
required_cols <- c("grid_id", "year", "month", "burned_area")
missing_cols <- setdiff(required_cols, names(combined_data))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

# Remove any rows with NA in key columns
combined_data <- combined_data[!is.na(grid_id) & !is.na(year) & !is.na(burned_area)]

# ==============================================================================
# AGGREGATE MONTHLY TO ANNUAL
# ==============================================================================

log_message("Aggregating monthly to annual...")

# Sum burned_area across months for each grid_id-year
# Note: burned_annual can exceed 1 if same area burns multiple times per year
annual_data <- combined_data[, .(
  burned_annual = sum(burned_area, na.rm = TRUE)
), by = .(grid_id, year)]

log_message(sprintf("Annual data: %d grid-year observations", nrow(annual_data)))

rm(combined_data)
gc()

# ==============================================================================
# DEDUPLICATE (safety check - should not happen)
# ==============================================================================

n_before <- nrow(annual_data)
dup_check <- annual_data[, .N, by = .(grid_id, year)][N > 1]

if (nrow(dup_check) > 0) {
  log_message(sprintf("WARNING: Found %d unexpected duplicates, using mean",
                      nrow(dup_check)))
  annual_data <- annual_data[, .(burned_annual = mean(burned_annual)),
                              by = .(grid_id, year)]
}

n_after <- nrow(annual_data)
if (n_before != n_after) {
  log_message(sprintf("Deduplicated: %d -> %d rows", n_before, n_after))
}

# ==============================================================================
# FILTER AND CLEAN
# ==============================================================================

# Remove zeros to save space (cells with no burning)
n_zeros <- sum(annual_data$burned_annual == 0)
annual_data <- annual_data[burned_annual > 0]
log_message(sprintf("Removed %d zero values", n_zeros))

# Sort for consistent output
setorder(annual_data, grid_id, year)

log_message(sprintf("Final output: %d observations", nrow(annual_data)))

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

if (nrow(annual_data) > 0) {
  log_message("Summary statistics:")
  log_message(sprintf("  Unique grid cells: %d", length(unique(annual_data$grid_id))))
  log_message(sprintf("  Years: %d-%d (%d years)",
                      min(annual_data$year), max(annual_data$year),
                      length(unique(annual_data$year))))
  log_message(sprintf("  Mean annual burned: %.6f", mean(annual_data$burned_annual)))
  log_message(sprintf("  Max annual burned: %.6f", max(annual_data$burned_annual)))

  # Year distribution
  year_summary <- annual_data[, .(n_cells = .N), by = year]
  setorder(year_summary, year)
  log_message("  Cells with burns by year (sample):")
  log_message(sprintf("    %d: %d cells", year_summary$year[1], year_summary$n_cells[1]))
  if (nrow(year_summary) > 1) {
    log_message(sprintf("    %d: %d cells",
                        year_summary$year[nrow(year_summary)],
                        year_summary$n_cells[nrow(year_summary)]))
  }
} else {
  log_message("No burned area detected in this tile")
}

# ==============================================================================
# WRITE OUTPUT
# ==============================================================================

log_message(sprintf("Writing output: %s", output_file))

write_atomic(annual_data, output_file)

# Verify write
if (file.exists(output_file)) {
  file_size <- file.info(output_file)$size / 1024
  log_message(sprintf("Output file size: %.1f KB", file_size))
} else {
  stop("Failed to write output file!")
}

# Clean up
gc_verbose()
log_job_end(start_time)
