# ==============================================================================
# STAGE 2B: CONSOLIDATE TMF DATA GLOBALLY
# ==============================================================================
# This script merges all tile-level consolidated files into a single
# global TMF dataset.
#
# Input: Data/build/tmf_consolidated/tile_*.parquet (86 files)
# Output: Data/build/tmf_consolidated/tmf_all.parquet
# ==============================================================================

# Load configuration
source("Code/BUILD_workspace.R")

# Record start time
start_time <- Sys.time()

log_job_start("2b_consolidate_global.R", task_id = 1)

# Check if output already exists
output_file <- file.path(consolidated_path, "tmf_all.parquet")
skip_if_exists(output_file, "global TMF consolidation")

# ==============================================================================
# FIND INPUT FILES
# ==============================================================================

log_message("Finding consolidated tile files...")

input_files <- sapply(1:N_TMF_TILES, get_consolidated_tile_filename)

# Check which files exist
files_exist <- file.exists(input_files)
log_message(sprintf("Found %d/%d tile files",
                    sum(files_exist), length(input_files)))

if (sum(files_exist) == 0) {
  stop("No consolidated tile files found. Run Stage 2a first.")
}

# Check for missing tiles
missing_tiles <- which(!files_exist)
if (length(missing_tiles) > 0) {
  log_message(sprintf("Warning: Missing tiles: %s",
                      paste(missing_tiles, collapse = ", ")))
}

# ==============================================================================
# READ AND COMBINE DATA (Memory-efficient approach)
# ==============================================================================

log_message("Reading and combining data...")

# Process in batches to manage memory
batch_size <- 10
n_batches <- ceiling(sum(files_exist) / batch_size)

# Create temporary storage for intermediate results
temp_dir <- file.path(consolidated_path, "temp_global")
if (!dir.exists(temp_dir)) dir.create(temp_dir)

existing_files <- input_files[files_exist]

for (batch in 1:n_batches) {
  batch_start <- (batch - 1) * batch_size + 1
  batch_end <- min(batch * batch_size, length(existing_files))

  log_message(sprintf("Processing batch %d/%d (files %d-%d)...",
                      batch, n_batches, batch_start, batch_end))

  batch_files <- existing_files[batch_start:batch_end]

  # Read batch files
  batch_data <- lapply(batch_files, function(f) {
    dt <- read_parquet(f)
    setDT(dt)
    return(dt)
  })

  # Combine batch
  batch_combined <- rbindlist(batch_data, use.names = TRUE, fill = TRUE)

  # Write batch to temp file
  batch_file <- file.path(temp_dir, sprintf("batch_%03d.parquet", batch))
  write_parquet(batch_combined, batch_file)

  log_message(sprintf("  Batch %d: %d rows", batch, nrow(batch_combined)))

  # Clean up
  rm(batch_data, batch_combined)
  gc()
}

# ==============================================================================
# MERGE BATCHES
# ==============================================================================

log_message("Merging batches...")

batch_files <- list.files(temp_dir, pattern = "batch_.*\\.parquet",
                          full.names = TRUE)
batch_files <- sort(batch_files)

log_message(sprintf("Found %d batch files to merge", length(batch_files)))

# Read all batches
all_batches <- lapply(batch_files, function(f) {
  dt <- read_parquet(f)
  setDT(dt)
  return(dt)
})

# Combine
combined_data <- rbindlist(all_batches, use.names = TRUE, fill = TRUE)

log_message(sprintf("Total combined rows: %d", nrow(combined_data)))

# Clean up batches
rm(all_batches)
gc_verbose()

# ==============================================================================
# FINAL PROCESSING
# ==============================================================================

log_message("Final processing...")

# Ensure correct column types
combined_data[, grid_id := as.character(grid_id)]
combined_data[, year := as.integer(year)]
combined_data[, tmf_class := as.character(tmf_class)]
combined_data[, fraction := as.numeric(fraction)]

# Remove duplicates if any (shouldn't happen, but safety check)
combined_data <- unique(combined_data)

# Sort for consistent output
setorder(combined_data, grid_id, year, tmf_class)

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

log_message("Summary statistics:")
log_message(sprintf("  Total rows: %s", format(nrow(combined_data), big.mark = ",")))
log_message(sprintf("  Unique grid cells: %s",
                    format(length(unique(combined_data$grid_id)), big.mark = ",")))
log_message(sprintf("  Years: %d-%d",
                    min(combined_data$year),
                    max(combined_data$year)))
log_message(sprintf("  Unique years: %d", length(unique(combined_data$year))))

# Class distribution
class_summary <- combined_data[, .(
  n_observations = .N,
  n_cells = length(unique(grid_id)),
  mean_fraction = mean(fraction)
), by = tmf_class]
setorder(class_summary, -n_observations)

log_message("TMF class distribution:")
for (i in seq_len(nrow(class_summary))) {
  log_message(sprintf("  %s: %s obs in %s cells, mean=%.3f",
                      class_summary$tmf_class[i],
                      format(class_summary$n_observations[i], big.mark = ","),
                      format(class_summary$n_cells[i], big.mark = ","),
                      class_summary$mean_fraction[i]))
}

# Year coverage
year_summary <- combined_data[, .(
  n_cells = length(unique(grid_id)),
  n_observations = .N
), by = year]
setorder(year_summary, year)

log_message("Year coverage (sample):")
log_message(sprintf("  First year (%d): %s cells",
                    year_summary$year[1],
                    format(year_summary$n_cells[1], big.mark = ",")))
log_message(sprintf("  Last year (%d): %s cells",
                    year_summary$year[nrow(year_summary)],
                    format(year_summary$n_cells[nrow(year_summary)], big.mark = ",")))

# ==============================================================================
# WRITE OUTPUT
# ==============================================================================

log_message(sprintf("Writing output: %s", output_file))

write_atomic(combined_data, output_file)

# Verify write
if (file.exists(output_file)) {
  file_size <- file.info(output_file)$size / 1024 / 1024 / 1024  # GB
  log_message(sprintf("Output file size: %.2f GB", file_size))
} else {
  stop("Failed to write output file!")
}

# ==============================================================================
# CLEANUP
# ==============================================================================

log_message("Cleaning up temporary files...")

# Remove temp batch files
unlink(temp_dir, recursive = TRUE)

gc_verbose()
log_job_end(start_time)
