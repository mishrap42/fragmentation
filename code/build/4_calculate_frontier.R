# ==============================================================================
# STAGE 4: CALCULATE FOREST FRONTIER
# ==============================================================================
# This script identifies grid cells that are "forest frontier" - cells within
# 100km of any forest interior cell.
#
# Uses efficient spatial indexing to avoid O(n^2) distance calculations.
# Each tile job processes cells in that tile plus a 100km buffer of interior
# cells from neighboring tiles.
#
# Input: SLURM_ARRAY_TASK_ID (1 to N_TMF_TILES)
# Output: Data/build/classifications/frontier_tile_{tile_id}.parquet
# ==============================================================================

# Load configuration
here::i_am('code/build/4_calculate_frontier.R')
source("code/build/BUILD_workspace.R")

# Record start time
start_time <- Sys.time()

# Get task ID
task_id <- get_slurm_task_id()
current_tile_id <- task_id

log_job_start("4_calculate_frontier.R", task_id)
log_message(sprintf("Calculating frontier for tile: %d", current_tile_id))

# Validate
if (current_tile_id < 1 || current_tile_id > N_TMF_TILES) {
  stop(sprintf("Invalid tile_id: %d", current_tile_id))
}

# Check if output already exists
output_file <- get_frontier_filename(current_tile_id)
skip_if_exists(output_file, sprintf("tile %d", current_tile_id))

# ==============================================================================
# LOAD GRID CELLS FOR THIS TILE (from sub-tile files)
# ==============================================================================

log_message("Loading grid cells for this tile...")

# Get all sub-tile grid files for this TMF tile
grid_files <- get_grid_files_for_tmf_tile(current_tile_id, "parquet")
grid_files <- grid_files[file.exists(grid_files)]

if (length(grid_files) == 0) {
  stop(sprintf("No grid files found for tile %d. Run Stage 0 first.", current_tile_id))
}

log_message(sprintf("Found %d sub-tile grid files", length(grid_files)))

# Load and combine all sub-tile grids
grid_list <- lapply(grid_files, function(f) {
  dt <- read_parquet(f)
  setDT(dt)
  dt
})
grid_data <- rbindlist(grid_list, fill = TRUE)

log_message(sprintf("Loaded %s grid cells from %d sub-tiles",
                    format(nrow(grid_data), big.mark = ","), length(grid_files)))

if (nrow(grid_data) == 0) {
  log_message("Empty grid. Creating empty output.")
  empty_output <- data.table(
    grid_id = character(0),
    tile_id = integer(0),
    is_frontier = logical(0),
    dist_to_interior_km = numeric(0)
  )
  write_atomic(empty_output, output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# ==============================================================================
# LOAD INTERIOR CLASSIFICATIONS
# ==============================================================================

log_message("Loading interior classifications...")

# Load interior data for this tile
interior_file <- get_interior_filename(current_tile_id)

if (!file.exists(interior_file)) {
  stop(sprintf("Interior file not found: %s. Run Stage 3 first.", interior_file))
}

interior_this_tile <- read_parquet(interior_file)
setDT(interior_this_tile)

log_message(sprintf("Interior cells in this tile: %s",
                    format(sum(interior_this_tile$is_interior), big.mark = ",")))

# ==============================================================================
# DETERMINE NEIGHBORING TILES
# ==============================================================================

log_message("Finding neighboring tiles for 100km buffer...")

# Get tile bounds
tile_info <- TMF_TILE_INDEX[tile_id == current_tile_id]

# Calculate which tiles could have interior cells within 100km
# 100km buffer in degrees (approximate: 1 degree ~ 111km at equator)
buffer_degrees <- BUFFER_DISTANCE_M / 111000 * 1.5  # Add 50% margin for projection

# Find tiles that overlap with buffered extent
neighboring_tiles <- TMF_TILE_INDEX[
  lat_north >= (tile_info$lat_south - buffer_degrees) &
    lat_south <= (tile_info$lat_north + buffer_degrees) &
    lon_east >= (tile_info$lon_west - buffer_degrees) &
    lon_west <= (tile_info$lon_east + buffer_degrees)
]$tile_id

# Remove current tile
neighboring_tiles <- setdiff(neighboring_tiles, current_tile_id)

log_message(sprintf("Neighboring tiles: %s", paste(neighboring_tiles, collapse = ", ")))

# ==============================================================================
# LOAD INTERIOR CELLS FROM NEIGHBORING TILES
# ==============================================================================

log_message("Loading interior cells from neighboring tiles...")

all_interior <- list()
all_interior[[1]] <- interior_this_tile[is_interior == TRUE]

for (neighbor_id in neighboring_tiles) {
  neighbor_file <- get_interior_filename(neighbor_id)

  if (file.exists(neighbor_file)) {
    neighbor_interior <- read_parquet(neighbor_file)
    setDT(neighbor_interior)

    # Only keep interior cells
    neighbor_interior <- neighbor_interior[is_interior == TRUE]

    if (nrow(neighbor_interior) > 0) {
      all_interior[[length(all_interior) + 1]] <- neighbor_interior
    }
  }
}

# Combine all interior cells
interior_combined <- rbindlist(all_interior, fill = TRUE)
log_message(sprintf("Total interior cells (including neighbors): %s",
                    format(nrow(interior_combined), big.mark = ",")))

# Clean up
rm(all_interior)
gc()

# ==============================================================================
# HANDLE CASE: NO INTERIOR CELLS NEARBY
# ==============================================================================

if (nrow(interior_combined) == 0) {
  log_message("No interior cells in or near this tile. All cells are non-frontier.")

  output_data <- data.table(
    grid_id = grid_data$grid_id,
    tile_id = current_tile_id,
    is_frontier = FALSE,
    dist_to_interior_km = NA_real_
  )

  write_atomic(output_data, output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# ==============================================================================
# GET COORDINATES FOR DISTANCE CALCULATION
# ==============================================================================

log_message("Preparing coordinates for distance calculation...")

# Get coordinates for grid cells
grid_coords <- grid_data[, .(grid_id, centroid_lon, centroid_lat)]

# Get coordinates for interior cells
# Need to load grid parquet files (from sub-tiles) to get centroid coordinates
interior_coords_list <- list()

for (tid in unique(interior_combined$tile_id)) {
  # Get all sub-tile grid files for this TMF tile
  grid_files_tid <- get_grid_files_for_tmf_tile(tid, "parquet")
  grid_files_tid <- grid_files_tid[file.exists(grid_files_tid)]

  for (grid_file_tid in grid_files_tid) {
    grid_tid <- read_parquet(grid_file_tid)
    setDT(grid_tid)

    # Filter to interior cells
    interior_ids <- interior_combined[tile_id == tid]$grid_id
    grid_tid <- grid_tid[grid_id %in% interior_ids]

    if (nrow(grid_tid) > 0) {
      interior_coords_list[[length(interior_coords_list) + 1]] <- grid_tid[, .(
        grid_id,
        centroid_lon,
        centroid_lat
      )]
    }
  }
}

interior_coords <- rbindlist(interior_coords_list)
log_message(sprintf("Interior cells with coordinates: %s",
                    format(nrow(interior_coords), big.mark = ",")))

rm(interior_coords_list)
gc()

# ==============================================================================
# EFFICIENT FRONTIER CALCULATION USING SPATIAL INDEXING
# ==============================================================================

log_message("Calculating frontier using spatial indexing...")

# Convert to sf objects
sf_use_s2(FALSE)

grid_sf <- st_as_sf(grid_coords, coords = c("centroid_lon", "centroid_lat"),
                    crs = WGS84_CRS)
interior_sf <- st_as_sf(interior_coords, coords = c("centroid_lon", "centroid_lat"),
                        crs = WGS84_CRS)

# Project to Mollweide for distance calculations
grid_sf_moll <- st_transform(grid_sf, MOLLWEIDE_CRS)
interior_sf_moll <- st_transform(interior_sf, MOLLWEIDE_CRS)

# ==============================================================================
# STRATEGY: K-NN WITH RANN FOR FAST NEAREST NEIGHBOR SEARCH
# ==============================================================================
# Use k-d tree based nearest neighbor search (RANN package).
# Complexity: O(n log m) instead of O(n * m) - orders of magnitude faster.
# ==============================================================================

log_message("Using k-NN strategy for fast distance calculation...")

# Extract projected coordinates as matrices for RANN
grid_coords_matrix <- st_coordinates(grid_sf_moll)
interior_coords_matrix <- st_coordinates(interior_sf_moll)

log_message(sprintf("Building k-d tree for %s interior cells...",
                    format(nrow(interior_coords_matrix), big.mark = ",")))

# Find nearest interior cell for each grid cell using k-d tree
# nn2() returns distances in the same units as input (meters for Mollweide)
nn_result <- RANN::nn2(
  data = interior_coords_matrix,   # Interior cells (build tree on these)
  query = grid_coords_matrix,      # Grid cells (query these)
  k = 1,                           # Only need nearest neighbor
  searchtype = "priority"          # Optimized for k=1
)

# Extract minimum distances (in km)
min_dists_km <- nn_result$nn.dists[, 1] / 1000

log_message(sprintf("Computed distances for %s grid cells",
                    format(length(min_dists_km), big.mark = ",")))

# Classify as frontier (within 100km of interior)
grid_coords[, dist_to_interior_km := min_dists_km]
grid_coords[, is_frontier := (dist_to_interior_km <= (BUFFER_DISTANCE_M / 1000))]

# Clean up
rm(nn_result, grid_coords_matrix, interior_coords_matrix)
gc()

sf_use_s2(TRUE)

# ==============================================================================
# MARK INTERIOR CELLS (they are NOT frontier, they ARE interior)
# ==============================================================================

log_message("Marking interior cells...")

interior_ids <- interior_this_tile[is_interior == TRUE]$grid_id
grid_coords[grid_id %in% interior_ids, is_frontier := FALSE]
grid_coords[grid_id %in% interior_ids, dist_to_interior_km := 0]

# ==============================================================================
# PREPARE OUTPUT
# ==============================================================================

log_message("Preparing output...")

output_data <- data.table(
  grid_id = grid_coords$grid_id,
  tile_id = current_tile_id,
  is_frontier = grid_coords$is_frontier,
  dist_to_interior_km = grid_coords$dist_to_interior_km
)

setorder(output_data, grid_id)

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

log_message("Summary statistics:")
log_message(sprintf("  Total cells: %s", format(nrow(output_data), big.mark = ",")))
log_message(sprintf("  Frontier cells: %s (%.1f%%)",
                    format(sum(output_data$is_frontier, na.rm = TRUE), big.mark = ","),
                    100 * mean(output_data$is_frontier, na.rm = TRUE)))
log_message(sprintf("  Interior cells: %s",
                    format(sum(output_data$dist_to_interior_km == 0, na.rm = TRUE), big.mark = ",")))
log_message(sprintf("  Cells beyond 100km: %s",
                    format(sum(!output_data$is_frontier & is.na(output_data$dist_to_interior_km), na.rm = TRUE), big.mark = ",")))

# Distance distribution for frontier cells
frontier_cells <- output_data[is_frontier == TRUE & dist_to_interior_km > 0]
if (nrow(frontier_cells) > 0) {
  log_message("Distance distribution for frontier cells (km):")
  log_message(sprintf("  Min: %.1f", min(frontier_cells$dist_to_interior_km, na.rm = TRUE)))
  log_message(sprintf("  Mean: %.1f", mean(frontier_cells$dist_to_interior_km, na.rm = TRUE)))
  log_message(sprintf("  Median: %.1f", median(frontier_cells$dist_to_interior_km, na.rm = TRUE)))
  log_message(sprintf("  Max: %.1f", max(frontier_cells$dist_to_interior_km, na.rm = TRUE)))
}

# ==============================================================================
# WRITE OUTPUT
# ==============================================================================

log_message(sprintf("Writing output: %s", output_file))

write_atomic(output_data, output_file)

if (file.exists(output_file)) {
  file_size <- file.info(output_file)$size / 1024  # KB
  log_message(sprintf("Output file size: %.1f KB", file_size))
}

gc_verbose()
log_job_end(start_time)
