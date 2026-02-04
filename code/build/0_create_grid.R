# ==============================================================================
# STAGE 0: CREATE 1KM EQUAL-AREA GRID
# ==============================================================================
# This script creates a 1km x 1km equal-area grid (Mollweide projection)
# for a single TMF tile, with grid cells trimmed at country borders.
#
# Input: SLURM_ARRAY_TASK_ID (1 to N_TMF_TILES)
# Output: Data/build/grids/grid_tile_{tile_id}.parquet
#         Data/build/grids/grid_tile_{tile_id}.gpkg
# ==============================================================================

# Load configuration
source("Code/BUILD_workspace.R")

# Record start time
start_time <- Sys.time()

# Get task ID
task_id <- get_slurm_task_id()
log_job_start("0_create_grid.R", task_id)

# Validate task ID
if (task_id < 1 || task_id > N_TMF_TILES) {
  stop(sprintf("Invalid task_id: %d (must be 1-%d)", task_id, N_TMF_TILES))
}

tile_id <- task_id

# Check if output already exists
output_parquet <- get_grid_filename(tile_id, "parquet")
output_gpkg <- get_grid_filename(tile_id, "gpkg")
skip_if_exists(output_parquet, sprintf("tile %d", tile_id))

log_message(sprintf("Processing tile %d", tile_id))

# ==============================================================================
# LOAD COUNTRY BOUNDARIES
# ==============================================================================

log_message("Loading GADM country boundaries...")

# Read GADM - only country level (admin0)
gadm <- st_read(gadm_path, layer = "ADM_ADM_0", quiet = TRUE)
log_message(sprintf("Loaded %d countries", nrow(gadm)))

# Make geometries valid
gadm <- make_valid_safe(gadm)

# ==============================================================================
# GET TILE EXTENT
# ==============================================================================

tile_info <- TMF_TILE_INDEX[tile_id == tile_id]
log_message(sprintf("Tile bounds: lat [%d, %d], lon [%d, %d]",
                    tile_info$lat_south, tile_info$lat_north,
                    tile_info$lon_west, tile_info$lon_east))

# Create tile extent polygon in WGS84
tile_bbox <- st_bbox(c(
  xmin = tile_info$lon_west,
  ymin = tile_info$lat_south,
  xmax = tile_info$lon_east,
  ymax = tile_info$lat_north
), crs = st_crs(WGS84_CRS))

tile_extent <- st_as_sfc(tile_bbox)

# ==============================================================================
# CLIP COUNTRIES TO TILE
# ==============================================================================

log_message("Clipping countries to tile extent...")

# Disable s2 for planar operations (faster, avoids some edge cases)
sf_use_s2(FALSE)

# Find countries that intersect the tile
countries_intersecting <- st_intersects(gadm, tile_extent, sparse = FALSE)
gadm_tile <- gadm[countries_intersecting[, 1], ]

log_message(sprintf("Found %d countries intersecting tile", nrow(gadm_tile)))

if (nrow(gadm_tile) == 0) {
  log_message("No countries in this tile (likely ocean). Creating empty output.")
  # Create empty output
  empty_grid <- data.table(
    grid_id = character(0),
    tile_id = integer(0),
    country_iso3 = character(0),
    country_name = character(0),
    centroid_lon = numeric(0),
    centroid_lat = numeric(0),
    area_km2 = numeric(0)
  )
  write_atomic(empty_grid, output_parquet)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# Clip countries to tile extent
gadm_clipped <- st_intersection(gadm_tile, tile_extent)
gadm_clipped <- make_valid_safe(gadm_clipped)

# Re-enable s2
sf_use_s2(TRUE)

# ==============================================================================
# PROJECT TO MOLLWEIDE
# ==============================================================================

log_message("Projecting to Mollweide equal-area...")

gadm_moll <- st_transform(gadm_clipped, MOLLWEIDE_CRS)
tile_extent_moll <- st_transform(tile_extent, MOLLWEIDE_CRS)

# ==============================================================================
# CREATE GRID
# ==============================================================================

log_message(sprintf("Creating %dm x %dm grid...", GRID_CELLSIZE_M, GRID_CELLSIZE_M))

# Create grid over the tile extent
grid_raw <- st_make_grid(
  tile_extent_moll,
  cellsize = c(GRID_CELLSIZE_M, GRID_CELLSIZE_M),
  what = "polygons"
)

log_message(sprintf("Created %d raw grid cells", length(grid_raw)))

# Convert to sf object
grid_sf <- st_sf(
  cell_idx = seq_along(grid_raw),
  geometry = grid_raw
)

# ==============================================================================
# INTERSECT WITH COUNTRY BOUNDARIES
# ==============================================================================

log_message("Intersecting grid with country boundaries...")

sf_use_s2(FALSE)

# Perform intersection - this trims cells at country borders
grid_country <- st_intersection(grid_sf, gadm_moll)

sf_use_s2(TRUE)

log_message(sprintf("Created %d grid-country intersections", nrow(grid_country)))

# ==============================================================================
# FILTER SLIVERS AND CALCULATE ATTRIBUTES
# ==============================================================================

log_message("Filtering slivers and calculating attributes...")

# Calculate area of each cell
grid_country$area_m2 <- as.numeric(st_area(grid_country))

# Target area for a full cell
target_area_m2 <- GRID_CELLSIZE_M^2

# Calculate fraction of target area
grid_country$area_frac <- grid_country$area_m2 / target_area_m2

# Filter out slivers (cells less than MIN_CELL_AREA_FRAC of target)
n_before <- nrow(grid_country)
grid_country <- grid_country[grid_country$area_frac >= MIN_CELL_AREA_FRAC, ]
n_after <- nrow(grid_country)

log_message(sprintf("Filtered %d slivers (%.1f%% of target area threshold)",
                    n_before - n_after, MIN_CELL_AREA_FRAC * 100))

# ==============================================================================
# CALCULATE CENTROIDS
# ==============================================================================

log_message("Calculating centroids...")

sf_use_s2(FALSE)

# Get centroids in Mollweide
centroids_moll <- st_centroid(grid_country)

# Transform centroids to WGS84 for lat/lon
centroids_wgs84 <- st_transform(centroids_moll, WGS84_CRS)
coords <- st_coordinates(centroids_wgs84)

sf_use_s2(TRUE)

grid_country$centroid_lon <- coords[, "X"]
grid_country$centroid_lat <- coords[, "Y"]

# ==============================================================================
# CREATE UNIQUE GRID IDs
# ==============================================================================

log_message("Creating unique grid IDs...")

# Grid ID format: tile_XXX_XXXXXXX (tile_id + sequential within tile)
grid_country$grid_id <- sprintf("tile_%03d_%07d", tile_id, seq_len(nrow(grid_country)))
grid_country$tile_id <- tile_id

# Extract country info
grid_country$country_iso3 <- grid_country$GID_0
grid_country$country_name <- grid_country$COUNTRY

# ==============================================================================
# PREPARE OUTPUT
# ==============================================================================

log_message("Preparing output...")

# Select columns for parquet (attributes only)
output_cols <- c("grid_id", "tile_id", "country_iso3", "country_name",
                 "centroid_lon", "centroid_lat", "area_m2", "area_frac")

grid_attrs <- as.data.table(st_drop_geometry(grid_country))
grid_attrs <- grid_attrs[, ..output_cols]

# Rename area column
setnames(grid_attrs, "area_m2", "area_km2")
grid_attrs[, area_km2 := area_km2 / 1e6]

log_message(sprintf("Final grid: %d cells in %d countries",
                    nrow(grid_attrs),
                    length(unique(grid_attrs$country_iso3))))

# ==============================================================================
# WRITE OUTPUT
# ==============================================================================

# Write parquet (attributes)
write_atomic(grid_attrs, output_parquet)

# Write geopackage (with geometry) for spatial operations
grid_output <- grid_country[, c("grid_id", "tile_id", "country_iso3",
                                 "centroid_lon", "centroid_lat", "area_frac")]

# Transform back to WGS84 for storage
grid_output_wgs84 <- st_transform(grid_output, WGS84_CRS)

log_message(sprintf("Writing geopackage: %s", output_gpkg))
st_write(grid_output_wgs84, output_gpkg, delete_dsn = TRUE, quiet = TRUE)

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

log_message("Summary statistics:")
log_message(sprintf("  Total cells: %d", nrow(grid_attrs)))
log_message(sprintf("  Countries: %s",
                    paste(unique(grid_attrs$country_iso3), collapse = ", ")))
log_message(sprintf("  Total area: %.2f km^2", sum(grid_attrs$area_km2)))
log_message(sprintf("  Mean cell area: %.4f km^2", mean(grid_attrs$area_km2)))
log_message(sprintf("  Min cell area: %.4f km^2 (%.1f%% of target)",
                    min(grid_attrs$area_km2),
                    min(grid_attrs$area_frac) * 100))

# Clean up
gc_verbose()
log_job_end(start_time)
