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
source("Code/BUILD_workspace.R")

# Record start time
start_time <- Sys.time()

# Get task ID
task_id <- get_slurm_task_id()
tile_id <- task_id

log_job_start("4_calculate_frontier.R", task_id)
log_message(sprintf("Calculating frontier for tile: %d", tile_id))

# Validate
if (tile_id < 1 || tile_id > N_TMF_TILES) {
  stop(sprintf("Invalid tile_id: %d", tile_id))
}

# Check if output already exists
output_file <- get_frontier_filename(tile_id)
skip_if_exists(output_file, sprintf("tile %d", tile_id))

# ==============================================================================
# LOAD GRID CELLS FOR THIS TILE
# ==============================================================================

log_message("Loading grid cells for this tile...")

grid_file <- get_grid_filename(tile_id, "parquet")

if (!file.exists(grid_file)) {
  stop(sprintf("Grid file not found: %s. Run Stage 0 first.", grid_file))
}

grid_data <- read_parquet(grid_file)
setDT(grid_data)

log_message(sprintf("Loaded %s grid cells", format(nrow(grid_data), big.mark = ",")))

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
interior_file <- get_interior_filename(tile_id)

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
tile_info <- TMF_TILE_INDEX[tile_id == tile_id]

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
neighboring_tiles <- setdiff(neighboring_tiles, tile_id)

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
    tile_id = tile_id,
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
# Need to load grid parquet files to get centroid coordinates
interior_coords_list <- list()

for (tid in unique(interior_combined$tile_id)) {
  grid_file_tid <- get_grid_filename(tid, "parquet")
  if (file.exists(grid_file_tid)) {
    grid_tid <- read_parquet(grid_file_tid)
    setDT(grid_tid)

    # Filter to interior cells
    interior_ids <- interior_combined[tile_id == tid]$grid_id
    grid_tid <- grid_tid[grid_id %in% interior_ids]

    interior_coords_list[[length(interior_coords_list) + 1]] <- grid_tid[, .(
      grid_id,
      centroid_lon,
      centroid_lat
    )]
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
# STRATEGY: CREATE BUFFER AROUND INTERIOR, THEN INTERSECT
# ==============================================================================

log_message("Creating 100km buffer around interior cells...")

# Union all interior points (more efficient than buffering each)
interior_union <- st_union(interior_sf_moll)

# Create 100km buffer
interior_buffer <- st_buffer(interior_union, dist = BUFFER_DISTANCE_M)

log_message("Identifying cells within buffer...")

# Check which grid cells intersect the buffer
in_buffer <- st_intersects(grid_sf_moll, interior_buffer, sparse = FALSE)[, 1]

n_in_buffer <- sum(in_buffer)
log_message(sprintf("Cells in 100km buffer: %s (%.1f%%)",
                    format(n_in_buffer, big.mark = ","),
                    100 * n_in_buffer / nrow(grid_coords)))

# ==============================================================================
# CALCULATE EXACT DISTANCES FOR CELLS IN BUFFER
# ==============================================================================

log_message("Calculating exact distances for cells in buffer...")

# Initialize results
grid_coords[, is_frontier := FALSE]
grid_coords[, dist_to_interior_km := NA_real_]

if (n_in_buffer > 0) {
  # Get indices of cells in buffer
  buffer_idx <- which(in_buffer)

  # Process in chunks to manage memory
  chunk_size <- 10000
  n_chunks <- ceiling(length(buffer_idx) / chunk_size)

  for (chunk in 1:n_chunks) {
    chunk_start <- (chunk - 1) * chunk_size + 1
    chunk_end <- min(chunk * chunk_size, length(buffer_idx))
    chunk_idx <- buffer_idx[chunk_start:chunk_end]

    if (chunk %% 10 == 1) {
      log_message(sprintf("Processing chunk %d/%d...", chunk, n_chunks))
    }

    # Calculate distances from chunk cells to all interior cells
    dist_matrix <- st_distance(grid_sf_moll[chunk_idx, ], interior_sf_moll)

    # Find minimum distance for each cell
    min_dists <- apply(dist_matrix, 1, min)
    min_dists_km <- as.numeric(min_dists) / 1000

    # Update results
    grid_coords[chunk_idx, dist_to_interior_km := min_dists_km]
    grid_coords[chunk_idx, is_frontier := (min_dists_km <= (BUFFER_DISTANCE_M / 1000))]
  }
}

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
  tile_id = tile_id,
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
