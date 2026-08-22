# ==============================================================================
# STAGE 5: EXTRACT STATIC COVARIATES
# ==============================================================================
# This script extracts static covariate data for each tile.
#
# Input: Grid files (stage 0), frontier classifications (stage 4)
# Output: Data/build/covariates/covariates_tile_{tile_id}.parquet
#
# Variables extracted:
#   1. Distance to major cities (km)
#   2. Population density (multiple years from GHS folders)
#   3. Population access (mean pop density within 50km, excl. own cell)
#   4. Above-ground biomass
#   5. Below-ground biomass
#   6. Biomass access (mean biomass within 50km, excl. own cell)
#   7. Signed distance to forest frontier
#   8. Elevation, slope, terrain ruggedness (ASTER)
#
# Usage: Rscript code/build/5_extract_covariates.R [tile_id]
#        sbatch code/bash/5_covariates.sh
#
# ==============================================================================
# CRITICAL: LESSONS LEARNED FROM exact_extract FAILURES
# ==============================================================================
#
# PROBLEM: "Mixed-type geometries not supported" errors during exact_extract
#
# ROOT CAUSE: Double CRS transformation corrupts polygon geometries.
#
# WHAT WENT WRONG (DO NOT DO THIS):
#   1. Load grid from gpkg (original CRS, e.g., WGS84)
#   2. Transform to WGS84: st_transform(grid_sf, WGS84_CRS)  <-- FIRST TRANSFORM
#   3. Normalize geometries (st_collection_extract, st_cast)
#   4. For each raster:
#      - Transform to raster CRS: st_transform(grid_sf, raster_crs)  <-- SECOND
#      - exact_extract() --> FAILS with "Mixed-type geometries"
#
#   The second transform (WGS84 -> Mollweide or other projection) corrupts
#   the polygon geometries, creating mixed types that exact_extract rejects.
#
# CORRECT APPROACH (follows 1_extract_TMF.R which has ZERO extraction errors):
#   1. Load grid from gpkg (keep original CRS - do NOT transform)
#   2. Normalize geometries IF mixed types detected (on original CRS)
#   3. For each raster:
#      - Transform ONCE to raster CRS: st_transform(grid_sf, raster_crs)
#      - exact_extract() --> WORKS
#
# WHY 1_extract_TMF.R WORKS:
#   - Grid loaded in WGS84, TMF rasters in WGS84
#   - st_transform is essentially a no-op (same CRS)
#   - No geometry corruption occurs
#
# KEY PRINCIPLES FOR exact_extract:
#   1. SINGLE TRANSFORM: Grid should only be transformed ONCE per raster
#   2. KEEP ORIGINAL CRS: Do not pre-transform grid to a "standard" CRS
#   3. NORMALIZE EARLY: Fix geometry types BEFORE any reprojection
#   4. NO BAND-AIDS: Do not use st_cast after st_transform as a workaround;
#      fix the root cause instead
#   5. DISABLE s2: Use sf_use_s2(FALSE) to avoid spherical geometry errors
#      with grids created using planar (GEOS) operations
#
# DEBUGGING TIPS:
#   - If exact_extract warns "Polygons transformed to raster CRS (EPSG:NA)",
#     it means YOUR transform didn't match exactly - check CRS handling
#   - Remove all tryCatch blocks during debugging to see actual errors
#   - Test with a single tile before running full array
#
# ==============================================================================

# Load configuration
here::i_am('code/build/5_extract_covariates.R')
source("code/build/BUILD_workspace.R")

# Disable s2 spherical geometry to avoid "Edge X crosses edge Y" errors
sf_use_s2(FALSE)

# Record start time
start_time <- Sys.time()

# Get tile ID from SLURM or command line
task_id <- get_slurm_task_id()
tile_id <- task_id

log_job_start("5_extract_covariates.R", task_id = tile_id)

# Check if output already exists
output_file <- get_covariates_filename(tile_id)
skip_if_exists(output_file, sprintf("tile %d covariates", tile_id))

# ==============================================================================
# LOAD GRID CELLS FOR THIS TILE
# ==============================================================================

log_message(sprintf("Loading grid cells for tile %d...", tile_id))

grid_sf <- load_grid_for_tile(tile_id, with_geometry = TRUE)
log_message(sprintf("Loaded %d grid cells", nrow(grid_sf)))

# Fix any invalid geometries (keep original CRS - do NOT transform to WGS84)
grid_sf <- st_make_valid(grid_sf)

# Normalize geometry types BEFORE any reprojection (same as 1_extract_TMF.R)
geom_types <- as.character(unique(st_geometry_type(grid_sf)))
geom_counts <- table(st_geometry_type(grid_sf))
log_message(sprintf("Geometry types: %s", paste(names(geom_counts), geom_counts, sep = "=", collapse = ", ")))

if (length(geom_types) > 1 || any(grepl("COLLECTION", geom_types))) {
  log_message("Normalizing geometries to MULTIPOLYGON...")
  n_before <- nrow(grid_sf)
  grid_sf <- st_collection_extract(grid_sf, "POLYGON")
  grid_sf <- st_cast(grid_sf, "MULTIPOLYGON")
  log_message(sprintf("  After normalization: %d features (was %d)", nrow(grid_sf), n_before))
}

# Store original CRS for reference
original_crs <- st_crs(grid_sf)
log_message(sprintf("Grid CRS: %s", original_crs$input))

# Compute centroids for distance calculations (need WGS84 for lat/lon)
grid_centroids_wgs <- st_transform(st_centroid(grid_sf), WGS84_CRS)
centroid_coords <- st_coordinates(grid_centroids_wgs)

# Initialize output data.table
covariate_data <- data.table(grid_id = grid_sf$grid_id)

# ==============================================================================
# 1. DISTANCE TO MAJOR CITIES
# ==============================================================================

log_message("Extracting distance to cities...")

stopifnot(file.exists(cities_path))
cities <- fread(cities_path)

stopifnot(all(c("lat", "lon") %in% names(cities)))

# Filter out rows with NA/NaN/Inf coordinates
cities <- cities[is.finite(lat) & is.finite(lon)]
log_message(sprintf("  Cities with valid coordinates: %d", nrow(cities)))

cities_coords <- as.matrix(cities[, .(lon, lat)])

# Use RANN for nearest neighbor search
nn <- RANN::nn2(cities_coords, centroid_coords, k = 1, searchtype = "priority")

# Convert degree distance to approximate km
lat_factor <- cos(centroid_coords[, 2] * pi / 180)
dist_deg <- nn$nn.dists[, 1]
covariate_data[, dist_to_city_km := dist_deg * 111 * sqrt((1 + lat_factor^2) / 2)]

log_message(sprintf("  Mean distance to city: %.1f km",
                    mean(covariate_data$dist_to_city_km, na.rm = TRUE)))

gc()

# ==============================================================================
# 2. POPULATION DENSITY (multiple years)
# ==============================================================================

log_message("Extracting population density...")

stopifnot(dir.exists(population_path))

pop_folders <- list.dirs(population_path, recursive = FALSE, full.names = TRUE)
pop_folders <- pop_folders[grepl("GHS_POP_E\\d{4}", basename(pop_folders))]

log_message(sprintf("  Found %d population folders", length(pop_folders)))
stopifnot(length(pop_folders) > 0)

latest_pop_year <- NULL
latest_pop_rast <- NULL

for (folder in pop_folders) {
  # Extract year from folder name
  year <- str_extract(basename(folder), "E(\\d{4})", group = 1)
  stopifnot(!is.na(year))

  # Find TIF file in folder
  tif_files <- list.files(folder, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
  stopifnot(length(tif_files) > 0)

  tif_file <- tif_files[1]
  varname <- sprintf("pop_density_%s", year)

  log_message(sprintf("  Processing %s...", varname))

  pop_rast <- terra::rast(tif_file)

  # Reproject grid to match raster CRS (single transform, like 1_extract_TMF.R)
  raster_crs <- terra::crs(pop_rast, proj = TRUE)
  grid_reproj <- st_transform(grid_sf, raster_crs)

  values <- exact_extract(pop_rast, grid_reproj, "mean")
  covariate_data[, (varname) := values]

  n_valid <- sum(!is.na(values))
  log_message(sprintf("    %d/%d valid (%.1f%%), mean=%.2f",
                      n_valid, length(values), 100 * n_valid / length(values),
                      mean(values, na.rm = TRUE)))

  # Track latest year for population access calculation
  year_int <- as.integer(year)
  if (is.null(latest_pop_year) || year_int > latest_pop_year) {
    latest_pop_year <- year_int
    latest_pop_rast <- pop_rast
  }

  gc()
}

# ==============================================================================
# 3. POPULATION ACCESS (50km buffer, excl. own cell)
# ==============================================================================

# exactextractr preloads the whole working area only if it fits in
# max_cells_in_memory (default 3e7 cells, ~240 MB); otherwise it re-reads the
# raster for EVERY feature. The 50km buffers overlap heavily, so that fallback
# is what made pop access take 58 minutes on 2026-08-19:
#     "Cannot preload entire working area of 152354220 cells
#      with max_cells_in_memory = 3e+07"
# Those working areas are ~128-152M cells, i.e. ~1.0-1.2 GB at 8 bytes/cell -
# trivial against this job's 120G. 3e8 covers them with headroom to spare.
#
# NOTE: this is deliberately NOT applied to the tile-sized extracts in
# 1_extract_TMF.R, whose working area is 1.38e9 cells (~11 GB). There the
# per-feature path is the right trade: it keeps --mem=8G and 200-way
# concurrency, which finished 1,000 tasks in ~10 minutes of wall clock.
BUFFER_MAX_CELLS <- 3e8

log_message("Extracting population access (50km buffer)...")

# Create 50km buffers around centroids in Mollweide (accurate distance)
grid_centroids_moll <- st_transform(grid_centroids_wgs, MOLLWEIDE_CRS)
buffers_50km <- st_buffer(grid_centroids_moll, 50000)  # 50km in meters

# Reproject buffers to match population raster CRS
raster_crs <- terra::crs(latest_pop_rast, proj = TRUE)
buffers_reproj <- st_transform(buffers_50km, raster_crs)

buffer_pop <- exact_extract(latest_pop_rast, buffers_reproj, "mean",
                            max_cells_in_memory = BUFFER_MAX_CELLS)
covariate_data[, pop_access_50km := buffer_pop]

log_message(sprintf("  Mean pop access (50km): %.1f", mean(buffer_pop, na.rm = TRUE)))

gc()

# ==============================================================================
# 4-5. BIOMASS (above/below ground)
# ==============================================================================

log_message("Extracting biomass data...")

stopifnot(dir.exists(biomass_path))

biomass_files <- list.files(biomass_path, pattern = "\\.tif$",
                             full.names = TRUE, ignore.case = TRUE)
biomass_files <- biomass_files[!grepl("uncertainty", biomass_files, ignore.case = TRUE)]

log_message(sprintf("  Found %d biomass files (excl. uncertainty)", length(biomass_files)))
stopifnot(length(biomass_files) > 0)

agb_rast <- NULL  # Track above-ground biomass raster for access calculation

for (f in biomass_files) {
  basename_clean <- tools::file_path_sans_ext(basename(f))
  varname <- tolower(gsub("[- ]+", "_", basename_clean))

  log_message(sprintf("  Processing %s...", varname))

  bio_rast <- terra::rast(f)

  # Identify above-ground biomass for access calculation
  if (grepl("above|agb", varname, ignore.case = TRUE)) {
    agb_rast <- bio_rast
  }

  # Reproject grid to match raster CRS (single transform, like 1_extract_TMF.R)
  raster_crs <- terra::crs(bio_rast, proj = TRUE)
  grid_reproj <- st_transform(grid_sf, raster_crs)

  values <- exact_extract(bio_rast, grid_reproj, "mean")
  covariate_data[, (varname) := values]

  n_valid <- sum(!is.na(values))
  log_message(sprintf("    %d/%d valid (%.1f%%), mean=%.2f",
                      n_valid, length(values), 100 * n_valid / length(values),
                      mean(values, na.rm = TRUE)))

  gc()
}

# ==============================================================================
# 6. BIOMASS ACCESS (50km buffer)
# ==============================================================================

log_message("Extracting biomass access (50km buffer)...")

stopifnot(!is.null(agb_rast))

# Reproject buffers to match biomass raster CRS
raster_crs <- terra::crs(agb_rast, proj = TRUE)
buffers_reproj <- st_transform(buffers_50km, raster_crs)

buffer_biomass <- exact_extract(agb_rast, buffers_reproj, "mean",
                                max_cells_in_memory = BUFFER_MAX_CELLS)
covariate_data[, biomass_access_50km := buffer_biomass]

log_message(sprintf("  Mean biomass access (50km): %.1f", mean(buffer_biomass, na.rm = TRUE)))

gc()

# ==============================================================================
# 7. SIGNED DISTANCE TO FRONTIER
# ==============================================================================

log_message("Computing signed distance to frontier...")

frontier_file <- get_frontier_filename(tile_id)
interior_file <- get_interior_filename(tile_id)

stopifnot(file.exists(frontier_file))
stopifnot(file.exists(interior_file))

frontier <- arrow::read_parquet(frontier_file)
interior <- arrow::read_parquet(interior_file)
setDT(frontier)
setDT(interior)

# Merge frontier data
covariate_data <- merge(
  covariate_data,
  frontier[, .(grid_id, is_frontier, dist_to_interior_km)],
  by = "grid_id",
  all.x = TRUE
)

# Merge interior data
covariate_data <- merge(
  covariate_data,
  interior[, .(grid_id, is_interior)],
  by = "grid_id",
  all.x = TRUE
)

# Fill NAs
covariate_data[is.na(is_frontier), is_frontier := FALSE]
covariate_data[is.na(is_interior), is_interior := FALSE]

# Create signed distance:
# - Interior: negative (-1 as placeholder)
# - Frontier: 0
# - Other: positive dist_to_interior_km
covariate_data[, signed_dist_to_frontier_km := fifelse(
  is_interior == TRUE, -1,
  fifelse(is_frontier == TRUE, 0, dist_to_interior_km)
)]

# Clean up intermediate columns
covariate_data[, c("is_frontier", "is_interior", "dist_to_interior_km") := NULL]

interior_count <- sum(covariate_data$signed_dist_to_frontier_km == -1, na.rm = TRUE)
frontier_count <- sum(covariate_data$signed_dist_to_frontier_km == 0, na.rm = TRUE)
log_message(sprintf("  Interior cells: %d, Frontier cells: %d", interior_count, frontier_count))

gc()

# ------------------------------------------------------------------------------
# Safe per-raster extraction.
#
# Task 1287248 lost 83 of 86 tasks to a single bad raster: roadDistance.tif has
# NO CRS, so st_transform() died with "crs not found: is it missing?" - AFTER
# ~1.5h of buffer extraction, and before livestock and elevation ran. Same
# late-death shape as the ASTER 403. A raster that cannot be used should cost
# its own column, not the whole job.
# ------------------------------------------------------------------------------
extract_raster_mean <- function(f, varname, grid_sf, covariate_data) {
  if (!file.exists(f)) {
    log_message(sprintf("    SKIPPED %s: file not found", varname)); return(FALSE)
  }
  r <- tryCatch(terra::rast(f), error = function(e) NULL)
  if (is.null(r)) {
    log_message(sprintf("    SKIPPED %s: unreadable", varname)); return(FALSE)
  }
  cs <- tryCatch(terra::crs(r, proj = TRUE), error = function(e) "")
  if (is.null(cs) || is.na(cs) || !nzchar(cs)) {
    log_message(sprintf("    SKIPPED %s: raster has no CRS - cannot be georeferenced", varname))
    return(FALSE)
  }
  vals <- tryCatch({
    grid_reproj <- st_transform(grid_sf, cs)
    exact_extract(r, grid_reproj, "mean")
  }, error = function(e) {
    log_message(sprintf("    SKIPPED %s: extraction failed (%s)", varname,
                        conditionMessage(e)))
    NULL
  })
  if (is.null(vals)) return(FALSE)

  covariate_data[, (varname) := vals]
  log_message(sprintf("    %d/%d valid (%.1f%%), mean=%.3f",
                      sum(!is.na(vals)), length(vals),
                      100 * sum(!is.na(vals)) / length(vals),
                      mean(vals, na.rm = TRUE)))
  TRUE
}

# ==============================================================================
# 7b. BIODIVERSITY (species richness + range-size rarity)
# ==============================================================================
# Proxies for conservation value - what made a site worth protecting. These are
# the covariates a PA-boundary RD leans on for its continuity check: richness
# should not jump at the boundary, and if it does, siting was endogenous to it.

log_message("Extracting biodiversity layers...")

if (!dir.exists(biodiversity_path)) {
  log_message(sprintf("  WARNING: biodiversity path not found: %s", biodiversity_path))
} else {
  biodiv_files <- c(
    list.files(biodiversity_path, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE),
    list.files(file.path(biodiversity_path, "IUCN"), pattern = "\\.tif$",
               full.names = TRUE, ignore.case = TRUE)
  )
  # .ovr pyramids and .aux.xml sidecars share the .tif stem - drop them.
  biodiv_files <- biodiv_files[!grepl("\\.(ovr|aux|vat)", biodiv_files, ignore.case = TRUE)]

  log_message(sprintf("  Found %d biodiversity rasters", length(biodiv_files)))

  for (f in biodiv_files) {
    varname <- paste0("biodiv_", tolower(gsub("[- ]+", "_",
                      tools::file_path_sans_ext(basename(f)))))
    log_message(sprintf("  Processing %s...", varname))

    extract_raster_mean(f, varname, grid_sf, covariate_data)
  }
  gc()
}

# ==============================================================================
# 7c. TRANSPORT: ROAD AND PORT DISTANCE
# ==============================================================================
# Distinct from the Nelson travel-time surfaces already extracted: those are
# cost-weighted accessibility, these are straight-line distance to the nearest
# road / port. Roads are the proximate deforestation driver, so the two are
# complementary rather than redundant.

log_message("Extracting road and port distance...")

# roadDistance.tif and portDistance.tif carry NO CRS and their extents are in
# metres with no recoverable projection (road: x [-16.9M, 17.1M]; port: origin
# at 0,0). There is no .prj/.aux.xml sidecar and the Global Forest Repo never
# used them either, so the projection cannot be recovered - assuming one would
# silently sample the wrong locations. They are attempted here so that a future
# fix to the source files is picked up automatically; until then
# extract_raster_mean() skips them with a logged reason.
#
# Accessibility is already covered by the Nelson travel-time surfaces
# (travel_time_cities, travel_time_ports), which are properly georeferenced.

for (nm in c("roadDistance", "portDistance")) {
  f <- file.path(dirname(transport_path), paste0(nm, ".tif"))
  varname <- tolower(gsub("([a-z])([A-Z])", "\\1_\\2", nm))
  log_message(sprintf("  Processing %s...", varname))
  extract_raster_mean(f, varname, grid_sf, covariate_data)
}
gc()

# ==============================================================================
# 7d. LIVESTOCK DENSITY (GLW3, 2010)
# ==============================================================================
# Cattle density under three allocation methods (Da = dasymetric, Aw = areal
# weighted, Ps = ...). Static 2010 snapshot, so a baseline pressure measure
# rather than a time-varying one. 8_Areakm is the cell-area denominator and is
# not itself a covariate.

log_message("Extracting livestock density (GLW3)...")

if (!dir.exists(livestock_path)) {
  log_message(sprintf("  WARNING: livestock path not found: %s", livestock_path))
} else {
  glw_files <- list.files(livestock_path, pattern = "_Ct_2010_[A-Za-z]+\\.tif$",
                          full.names = TRUE)
  glw_files <- glw_files[!grepl("\\.(ovr|aux|vat)", glw_files, ignore.case = TRUE)]
  log_message(sprintf("  Found %d GLW3 cattle rasters", length(glw_files)))

  for (f in glw_files) {
    stem <- tools::file_path_sans_ext(basename(f))          # e.g. 5_Ct_2010_Da
    suffix <- tolower(sub(".*_Ct_2010_", "", stem))          # da / aw / ps
    varname <- paste0("cattle_density_2010_", suffix)
    log_message(sprintf("  Processing %s...", varname))

    extract_raster_mean(f, varname, grid_sf, covariate_data)
  }
  gc()
}

# ==============================================================================
# 8. ELEVATION, SLOPE, TERRAIN RUGGEDNESS
# ==============================================================================

# ------------------------------------------------------------------------------
# Copernicus DEM 30m, not ASTER.
#
# The previous source was a single global GeoTIFF streamed from
# storage.googleapis.com/natcap-data-cache. That bucket revoked anonymous read:
# every one of the 86 tasks on 2026-08-20 died here with
#     GDAL Error 11: HTTP response code: 403  (AccessDenied)
# after ~1.5h of upstream extraction. It is dead for everyone, not blocked by
# the cluster - compute nodes do have outbound HTTPS.
#
# get_copernicus_dem() is ported from the Global Forest Repo's
# code/1a_add_elevation.R, which has been using this source successfully. The
# DEM is served as 1-degree COG tiles from AWS; tiles over open ocean simply do
# not exist, so missing ones are skipped rather than treated as errors.
# ------------------------------------------------------------------------------

#' Fetch and merge the Copernicus 30m DEM tiles covering a bbox
#' @param extent_bbox an st_bbox in WGS84
#' @return a SpatRaster, or NULL if no tiles exist (all-ocean extent)
get_copernicus_dem <- function(extent_bbox) {
  # Expand by 1 degree: slope/TRI at the tile edge need neighbouring pixels.
  xmin <- floor(extent_bbox[["xmin"]]) - 1
  xmax <- ceiling(extent_bbox[["xmax"]]) + 1
  ymin <- floor(extent_bbox[["ymin"]]) - 1
  ymax <- ceiling(extent_bbox[["ymax"]]) + 1

  urls <- c()
  for (lon in xmin:(xmax - 1)) {
    for (lat in ymin:(ymax - 1)) {
      lat_str <- ifelse(lat >= 0, sprintf("N%02d_00", lat), sprintf("S%02d_00", abs(lat)))
      lon_str <- ifelse(lon >= 0, sprintf("E%03d_00", lon), sprintf("W%03d_00", abs(lon)))
      urls <- c(urls, sprintf(
        "/vsicurl/https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_%s_%s_DEM/Copernicus_DSM_COG_10_%s_%s_DEM.tif",
        lat_str, lon_str, lat_str, lon_str))
    }
  }

  log_message(sprintf("  Attempting %d Copernicus DEM tiles...", length(urls)))

  # terra::rast() is lazy - headers only - so holding the list costs no pixels.
  valid_rasters <- list()
  for (url in urls) {
    r <- tryCatch(terra::rast(url),
                  error   = function(e) NULL,   # tile does not exist (ocean)
                  warning = function(w) NULL)
    if (!is.null(r)) valid_rasters[[length(valid_rasters) + 1]] <- r
  }
  log_message(sprintf("  Located %d of %d tiles", length(valid_rasters), length(urls)))

  if (length(valid_rasters) == 0) return(NULL)
  if (length(valid_rasters) == 1) return(valid_rasters[[1]])

  # Merge to a temp GeoTIFF rather than into RAM: block-by-block, so memory
  # stays bounded however many tiles a grid tile spans.
  merged_tmp <- tempfile(fileext = ".tif")
  do.call(terra::merge,
          c(valid_rasters, list(filename = merged_tmp, overwrite = TRUE)))
}

log_message("Extracting elevation and terrain data from Copernicus DEM 30m...")

# Get tile extent with buffer for slope edge effects
tile_extent <- get_tile_extent_sf(tile_id, crs = WGS84_CRS)
tile_bbox_buf <- st_buffer(tile_extent, 0.1)  # 0.1 degree buffer (~10km)

dem_full <- get_copernicus_dem(st_bbox(tile_bbox_buf))
if (is.null(dem_full)) {
  stop("No Copernicus DEM tiles cover tile ", tile_id,
       " - unexpected for a tile that has land cells.")
}

# Crop to tile extent
log_message("  Cropping DEM to tile extent...")
aster_tile <- terra::crop(dem_full, terra::vect(tile_bbox_buf))

stopifnot(terra::ncell(aster_tile) > 0)
log_message(sprintf("  DEM tile: %d x %d pixels", ncol(aster_tile), nrow(aster_tile)))

# Reproject grid to match DEM CRS (single transform, like 1_extract_TMF.R)
raster_crs <- terra::crs(aster_tile, proj = TRUE)
grid_reproj <- st_transform(grid_sf, raster_crs)
log_message(sprintf("  Reprojected grid to DEM CRS: %s", raster_crs))

# Mean elevation
log_message("  Extracting mean elevation...")
elev_values <- exact_extract(aster_tile, grid_reproj, "mean")
covariate_data[, elevation_m := elev_values]
log_message(sprintf("    %d/%d valid, mean=%.1fm",
                    sum(!is.na(elev_values)), length(elev_values),
                    mean(elev_values, na.rm = TRUE)))

# Calculate slope (degrees)
log_message("  Calculating slope...")
slope_rast <- terra::terrain(aster_tile, v = "slope", unit = "degrees")
slope_values <- exact_extract(slope_rast, grid_reproj, "mean")
covariate_data[, slope_degrees := slope_values]
log_message(sprintf("    %d/%d valid, mean=%.1f degrees",
                    sum(!is.na(slope_values)), length(slope_values),
                    mean(slope_values, na.rm = TRUE)))

# Calculate TRI (Terrain Ruggedness Index)
log_message("  Calculating terrain ruggedness...")
tri_rast <- terra::terrain(aster_tile, v = "TRI")
tri_values <- exact_extract(tri_rast, grid_reproj, "mean")
covariate_data[, terrain_ruggedness := tri_values]
log_message(sprintf("    %d/%d valid, mean=%.1f",
                    sum(!is.na(tri_values)), length(tri_values),
                    mean(tri_values, na.rm = TRUE)))

# Fraction steep (>15 degrees)
log_message("  Calculating steep slope fraction...")
steep_mask <- slope_rast > 15
steep_values <- exact_extract(steep_mask, grid_reproj, "mean")
covariate_data[, frac_slope_gt15 := steep_values]
log_message(sprintf("    %d/%d valid, mean=%.2f",
                    sum(!is.na(steep_values)), length(steep_values),
                    mean(steep_values, na.rm = TRUE)))

gc()

# ==============================================================================
# FINAL VALIDATION REPORT
# ==============================================================================

log_message("")
log_message("========================================")
log_message("EXTRACTION VALIDATION REPORT")
log_message("========================================")

numeric_cols <- names(covariate_data)[sapply(covariate_data, is.numeric)]

for (col in numeric_cols) {
  vals <- covariate_data[[col]]
  n_total <- length(vals)
  n_na <- sum(is.na(vals))
  n_valid <- n_total - n_na
  pct_valid <- 100 * n_valid / n_total

  status <- if (n_valid == 0) {
    "FAILED"
  } else if (pct_valid < 50) {
    "PARTIAL"
  } else {
    "OK"
  }

  log_message(sprintf("  %-30s: %6d/%6d valid (%5.1f%%) - %s",
                      col, n_valid, n_total, pct_valid, status))
}

log_message("========================================")

# ==============================================================================
# WRITE OUTPUT
# ==============================================================================

log_message("Writing output...")

# Ensure output directory exists
if (!dir.exists(dirname(output_file))) {
  dir.create(dirname(output_file), recursive = TRUE)
}

# Handle duplicate grid_ids from geometry normalization (aggregate with mean)
if (anyDuplicated(covariate_data$grid_id)) {
  n_before <- nrow(covariate_data)
  numeric_cols <- names(covariate_data)[sapply(covariate_data, is.numeric)]
  covariate_data <- covariate_data[, lapply(.SD, mean, na.rm = TRUE),
                                   by = grid_id, .SDcols = numeric_cols]
  log_message(sprintf("Aggregated duplicates: %d -> %d rows", n_before, nrow(covariate_data)))
}

log_message(sprintf("Output columns: %s", paste(names(covariate_data), collapse = ", ")))
log_message(sprintf("Output rows: %d", nrow(covariate_data)))

write_atomic(covariate_data, output_file)

if (file.exists(output_file)) {
  file_size <- file.info(output_file)$size / 1024
  log_message(sprintf("Output file size: %.1f KB", file_size))
}

gc_verbose()
log_job_end(start_time)
