# ==============================================================================
# STAGE 3: CLASSIFY FOREST INTERIOR
# ==============================================================================
# This script identifies grid cells that are "forest interior" - cells that
# had >0% undisturbed forest in EVERY year of the TMF data (1990-2024).
#
# Input: SLURM_ARRAY_TASK_ID (1 to N_TMF_TILES)
# Output: Data/build/classifications/interior_tile_{tile_id}.parquet
# ==============================================================================

# Load configuration
source("Code/BUILD_workspace.R")

# Record start time
start_time <- Sys.time()

# Get task ID
task_id <- get_slurm_task_id()
tile_id <- task_id

log_job_start("3_classify_interior.R", task_id)
log_message(sprintf("Classifying interior for tile: %d", tile_id))

# Validate
if (tile_id < 1 || tile_id > N_TMF_TILES) {
  stop(sprintf("Invalid tile_id: %d", tile_id))
}

# Check if output already exists
output_file <- get_interior_filename(tile_id)
skip_if_exists(output_file, sprintf("tile %d", tile_id))

# ==============================================================================
# LOAD TMF DATA
# ==============================================================================

log_message("Loading consolidated TMF data...")

input_file <- get_consolidated_tile_filename(tile_id)

if (!file.exists(input_file)) {
  stop(sprintf("Consolidated file not found: %s. Run Stage 2a first.", input_file))
}

tmf_data <- read_parquet(input_file)
setDT(tmf_data)

log_message(sprintf("Loaded %s rows", format(nrow(tmf_data), big.mark = ",")))

if (nrow(tmf_data) == 0) {
  log_message("Empty input data. Creating empty output.")
  empty_output <- data.table(
    grid_id = character(0),
    tile_id = integer(0),
    is_interior = logical(0),
    min_undisturbed_frac = numeric(0),
    n_years_observed = integer(0)
  )
  write_atomic(empty_output, output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# ==============================================================================
# FILTER TO UNDISTURBED TMF
# ==============================================================================

log_message("Filtering to Undisturbed TMF class...")

# Filter to only undisturbed forest
undisturbed <- tmf_data[tmf_class == "Undisturbed_TMF"]

log_message(sprintf("Undisturbed TMF observations: %s",
                    format(nrow(undisturbed), big.mark = ",")))

if (nrow(undisturbed) == 0) {
  log_message("No undisturbed forest in this tile. All cells are non-interior.")

  # Get all unique grid cells from the original data
  all_cells <- unique(tmf_data$grid_id)

  output_data <- data.table(
    grid_id = all_cells,
    tile_id = tile_id,
    is_interior = FALSE,
    min_undisturbed_frac = 0,
    n_years_observed = 0L
  )

  write_atomic(output_data, output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# ==============================================================================
# CALCULATE INTERIOR STATUS
# ==============================================================================

log_message("Calculating interior status...")

# For each grid cell, calculate:
# 1. Number of years with undisturbed forest data
# 2. Minimum undisturbed fraction across all years

cell_stats <- undisturbed[, .(
  min_undisturbed_frac = min(fraction, na.rm = TRUE),
  max_undisturbed_frac = max(fraction, na.rm = TRUE),
  mean_undisturbed_frac = mean(fraction, na.rm = TRUE),
  n_years_observed = length(unique(year))
), by = grid_id]

log_message(sprintf("Calculated stats for %s cells",
                    format(nrow(cell_stats), big.mark = ",")))

# ==============================================================================
# CLASSIFY INTERIOR
# ==============================================================================

log_message("Classifying interior cells...")

# Interior definition: >0% undisturbed in ALL years
# Note: We check for presence in all TMF years
expected_years <- N_TMF_YEARS

cell_stats[, is_interior := (min_undisturbed_frac > 0) &
             (n_years_observed >= expected_years)]

# Count
n_interior <- sum(cell_stats$is_interior)
n_total <- nrow(cell_stats)
pct_interior <- 100 * n_interior / n_total

log_message(sprintf("Interior cells: %s / %s (%.1f%%)",
                    format(n_interior, big.mark = ","),
                    format(n_total, big.mark = ","),
                    pct_interior))

# ==============================================================================
# ADD NON-FORESTED CELLS
# ==============================================================================

log_message("Adding cells without any undisturbed forest...")

# Get all unique grid cells
all_cells <- unique(tmf_data$grid_id)

# Find cells not in our stats (never had undisturbed forest)
missing_cells <- setdiff(all_cells, cell_stats$grid_id)

if (length(missing_cells) > 0) {
  log_message(sprintf("Adding %s cells without undisturbed forest",
                      format(length(missing_cells), big.mark = ",")))

  missing_data <- data.table(
    grid_id = missing_cells,
    min_undisturbed_frac = 0,
    max_undisturbed_frac = 0,
    mean_undisturbed_frac = 0,
    n_years_observed = 0L,
    is_interior = FALSE
  )

  cell_stats <- rbind(cell_stats, missing_data, fill = TRUE)
}

# ==============================================================================
# PREPARE OUTPUT
# ==============================================================================

log_message("Preparing output...")

# Add tile_id
cell_stats[, tile_id := tile_id]

# Select and order columns
output_data <- cell_stats[, .(
  grid_id,
  tile_id,
  is_interior,
  min_undisturbed_frac,
  max_undisturbed_frac,
  mean_undisturbed_frac,
  n_years_observed
)]

# Sort
setorder(output_data, grid_id)

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

log_message("Summary statistics:")
log_message(sprintf("  Total cells: %s", format(nrow(output_data), big.mark = ",")))
log_message(sprintf("  Interior cells: %s (%.1f%%)",
                    format(sum(output_data$is_interior), big.mark = ","),
                    100 * mean(output_data$is_interior)))
log_message(sprintf("  Non-interior cells: %s",
                    format(sum(!output_data$is_interior), big.mark = ",")))

# Distribution of years observed
log_message("Years observed distribution:")
year_dist <- output_data[, .N, by = n_years_observed]
setorder(year_dist, n_years_observed)
for (i in seq_len(nrow(year_dist))) {
  log_message(sprintf("  %d years: %s cells",
                      year_dist$n_years_observed[i],
                      format(year_dist$N[i], big.mark = ",")))
}

# Undisturbed fraction distribution (for interior cells)
if (sum(output_data$is_interior) > 0) {
  interior_cells <- output_data[is_interior == TRUE]
  log_message("Undisturbed fraction for interior cells:")
  log_message(sprintf("  Min: %.3f", min(interior_cells$min_undisturbed_frac)))
  log_message(sprintf("  Mean: %.3f", mean(interior_cells$mean_undisturbed_frac)))
  log_message(sprintf("  Max: %.3f", max(interior_cells$max_undisturbed_frac)))
}

# ==============================================================================
# WRITE OUTPUT
# ==============================================================================

log_message(sprintf("Writing output: %s", output_file))

write_atomic(output_data, output_file)

# Verify write
if (file.exists(output_file)) {
  file_size <- file.info(output_file)$size / 1024  # KB
  log_message(sprintf("Output file size: %.1f KB", file_size))
}

gc_verbose()
log_job_end(start_time)
