# ==============================================================================
# STAGE 0: CREATE 5KM EQUAL-AREA GRID (5x5 DEGREE SUB-TILES)
# ==============================================================================
# This script creates a 5km x 5km equal-area grid (Mollweide projection)
# for a single 5x5 degree sub-tile, with grid cells cut at:
#   - Country borders (each cell fully in one country)
#   - WDPA protected area borders (each cell fully protected or not)
#
# Input: SLURM_ARRAY_TASK_ID (1 to N_SUB_TILES)
# Output: Data/build/grids/grid_sub_{sub_tile_id}.parquet
#         Data/build/grids/grid_sub_{sub_tile_id}.gpkg
#
# PREREQUISITE: Run 0a_clean_wdpa.R first to create wdpa_cleaned.gpkg
# ==============================================================================
#
# GEOMETRY CLEANING STEPS (applied to handle problematic WDPA polygons):
# ------------------------------------------------------------------------------
# 1. BBOX FILTERING (lines ~222-244): Drop WDPA geometries with bounding boxes
#    that span >90 degrees or don't overlap the tile extent. These are typically
#    corrupt/invalid geometries from the source data.
#
# 2. ST_SIMPLIFY (lines ~395-398): Simplify WDPA geometries with dTolerance=1m
#    before st_union(). This removes degenerate vertices (zero-length line
#    segments) that cause TopologyException in GEOS operations.
#
# 3. GEOMETRYCOLLECTION HANDLING (lines ~544-570): After all intersections/
#    differences, extract only POLYGON/MULTIPOLYGON components from any
#    GEOMETRYCOLLECTIONs, then cast to MULTIPOLYGON. This ensures exact_extract
#    can process the final grid (it doesn't support mixed geometry types).
#
# To quantify the impact of these cleaning steps on WDPA coverage, run:
#   Rscript code/build/diagnose_wdpa_cleaning.R
# ==============================================================================

# Load configuration
here::i_am('code/build/0_create_grid.R')
source("code/build/BUILD_workspace.R")

log_message("TRACE: BUILD_workspace.R loaded successfully")

# Record start time
start_time <- Sys.time()

# Get task ID (= sub_tile_id)
task_id <- get_slurm_task_id()
sub_tile_id <- task_id

log_job_start("0_create_grid.R", task_id)

# Get sub-tile info
log_message(sprintf("TRACE: Getting sub-tile info for sub_tile_id=%d", sub_tile_id))
sub_tile <- get_sub_tile_info(sub_tile_id)
log_message(sprintf("TRACE: Sub-tile %d: lat [%d, %d], lon [%d, %d], tmf_tile_id=%d",
                    sub_tile_id,
                    sub_tile$lat_south, sub_tile$lat_north,
                    sub_tile$lon_west, sub_tile$lon_east,
                    sub_tile$tmf_tile_id))

# Check if output already exists
output_parquet <- get_grid_filename(sub_tile_id, "parquet")
output_gpkg <- get_grid_filename(sub_tile_id, "gpkg")
log_message(sprintf("TRACE: Output files: %s, %s", output_parquet, output_gpkg))
skip_if_exists(output_parquet, sprintf("sub-tile %d", sub_tile_id))

# ==============================================================================
# CHECK WDPA PREREQUISITE
# ==============================================================================

if (!file.exists(wdpa_clean_path)) {
  stop(sprintf("Cleaned WDPA not found: %s\nRun 0a_clean_wdpa.R first.", wdpa_clean_path))
}

# ==============================================================================
# CREATE SUB-TILE EXTENT
# ==============================================================================

log_message("TRACE: Creating tile_bbox...")
tile_bbox <- st_bbox(c(
  xmin = sub_tile$lon_west,
  ymin = sub_tile$lat_south,
  xmax = sub_tile$lon_east,
  ymax = sub_tile$lat_north
), crs = st_crs(WGS84_CRS))

log_message("TRACE: Creating tile_extent from tile_bbox...")
tile_extent <- st_as_sfc(tile_bbox)
log_message("TRACE: Creating tile_extent_wkt...")
tile_extent_wkt <- st_as_text(tile_extent)
log_message("TRACE: Sub-tile extent created successfully")

# ==============================================================================
# LOAD COUNTRY BOUNDARIES (for this sub-tile extent only)
# ==============================================================================

log_message("Loading GADM country boundaries...")

log_message(sprintf("TRACE: Reading GADM from: %s", gadm_path))
gadm <- tryCatch({
  log_message("TRACE: Attempting st_read with wkt_filter...")
  result <- st_read(gadm_path, layer = "ADM_0", quiet = TRUE, wkt_filter = tile_extent_wkt)
  log_message("TRACE: st_read with wkt_filter succeeded")
  result
}, error = function(e) {
  log_message(sprintf("TRACE: Spatial filter failed with error: %s", e$message))
  log_message("TRACE: Using manual filter instead...")
  gadm_full <- st_read(gadm_path, layer = "ADM_0", quiet = TRUE)
  log_message(sprintf("TRACE: Loaded full GADM with %d rows", nrow(gadm_full)))
  sf_use_s2(FALSE)
  intersects <- st_intersects(gadm_full, tile_extent, sparse = FALSE)[, 1]
  sf_use_s2(TRUE)
  log_message(sprintf("TRACE: Manual filter found %d intersecting countries", sum(intersects)))
  gadm_full[intersects, ]
})

log_message(sprintf("TRACE: Loaded %d countries intersecting sub-tile", nrow(gadm)))

if (nrow(gadm) == 0) {
  log_message("No countries in this sub-tile (ocean). Creating empty output.")
  empty_grid <- data.table(
    grid_id = character(0), sub_tile_id = integer(0), tile_id = integer(0),
    country_iso3 = character(0), country_name = character(0),
    centroid_lon = numeric(0), centroid_lat = numeric(0), area_km2 = numeric(0),
    is_protected = logical(0), wdpa_pid = integer(0), iucn_cat = character(0),
    desig_year = integer(0), gov_type = character(0),
    travel_time_cities = numeric(0), travel_time_ports = numeric(0)
  )
  write_atomic(empty_grid, output_parquet)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

log_message("TRACE: Calling make_valid_safe on GADM (before transform)...")
gadm <- make_valid_safe(gadm)
log_message("TRACE: Transforming GADM to WGS84...")
gadm <- st_transform(gadm, WGS84_CRS)
log_message("TRACE: Calling make_valid_safe on GADM (after transform)...")
gadm <- make_valid_safe(gadm)

# Drop any geometries still invalid after repair
still_invalid <- !st_is_valid(gadm)
if (any(still_invalid)) {
  log_message(sprintf("TRACE: Dropping %d still-invalid GADM geometries", sum(still_invalid)))
  gadm <- gadm[!still_invalid, ]
}

# ------------------------------------------------------------------------------
# Drop distant island parts, THEN test for antimeridian crossing.
#
# GADM ADM_0 stores each country as a single multipolygon covering every
# offshore territory it holds. The USA's therefore includes the Aleutians,
# which straddle 180 deg, so the whole-country bbox spans -180..180 - and the
# old test discarded the ENTIRE country, the lower 48 included, for every tile
# the USA touched. Observed in job 1109967: sub-tiles 2 and 3 silently lost
# their US half, and sub-tile 4 (US only) produced an empty grid and exited 0.
# Russia, Fiji and Kiribati would fail the same way.
#
# Fix: explode to individual polygons, keep only the parts that are actually
# near this sub-tile (the distant islands are dropped here), then apply the
# antimeridian test per part. A reassembled country therefore keeps its
# mainland, and the test only fires on geometry genuinely sitting on 180 deg.
#
# Parts are reassembled to one row per country before returning: the later
# st_intersection(grid_sf, gadm_moll) at line ~406 would otherwise emit one
# duplicate cell per island part.
# ------------------------------------------------------------------------------

n_countries_before <- nrow(gadm)
gadm_parts <- suppressWarnings(st_cast(st_cast(gadm, "MULTIPOLYGON"), "POLYGON"))
log_message(sprintf("TRACE: Exploded %d countries into %d polygon parts",
                    n_countries_before, nrow(gadm_parts)))

# Keep parts intersecting the sub-tile, with a 1 degree margin so geometry
# just off the edge still informs boundary cells.
part_zone <- st_as_sfc(st_bbox(c(
  xmin = as.numeric(st_bbox(tile_extent)["xmin"]) - 1,
  ymin = as.numeric(st_bbox(tile_extent)["ymin"]) - 1,
  xmax = as.numeric(st_bbox(tile_extent)["xmax"]) + 1,
  ymax = as.numeric(st_bbox(tile_extent)["ymax"]) + 1
), crs = st_crs(gadm_parts)))

sf_use_s2(FALSE)
near_part <- lengths(st_intersects(gadm_parts, part_zone)) > 0
sf_use_s2(TRUE)
log_message(sprintf("TRACE: Keeping %d of %d parts near this sub-tile (dropped %d distant island parts)",
                    sum(near_part), length(near_part), sum(!near_part)))
gadm_parts <- gadm_parts[near_part, ]

# Now the antimeridian test, per surviving part
if (nrow(gadm_parts) > 0) {
  am_part <- vapply(seq_len(nrow(gadm_parts)), function(i) {
    bb <- st_bbox(gadm_parts[i, ])
    isTRUE(bb[["xmin"]] <= -179 && bb[["xmax"]] >= 179)
  }, logical(1))
  if (any(am_part)) {
    log_message(sprintf("TRACE: Dropping %d parts genuinely spanning the antimeridian",
                        sum(am_part)))
    gadm_parts <- gadm_parts[!am_part, ]
  }
}

# Reassemble: one row per country, attributes taken from its first part
if (nrow(gadm_parts) == 0) {
  gadm <- gadm[0, ]
} else {
  key <- if ("GID_0" %in% names(gadm_parts)) as.character(gadm_parts$GID_0) else
         rep("ALL", nrow(gadm_parts))
  gadm <- do.call(rbind, lapply(unique(key), function(k) {
    sub_parts <- gadm_parts[key == k, ]
    out <- sub_parts[1, ]
    st_geometry(out) <- st_union(st_geometry(sub_parts))
    out
  }))
  gadm <- make_valid_safe(gadm)
}

log_message(sprintf("TRACE: %d countries retained after part filtering (was %d)",
                    nrow(gadm), n_countries_before))

if (nrow(gadm) == 0) {
  log_message("No valid GADM geometries after filtering. Creating empty output.")
  empty_grid <- data.table(
    grid_id = character(0), sub_tile_id = integer(0), tile_id = integer(0),
    country_iso3 = character(0), country_name = character(0),
    centroid_lon = numeric(0), centroid_lat = numeric(0), area_km2 = numeric(0),
    is_protected = logical(0), wdpa_pid = integer(0), iucn_cat = character(0),
    desig_year = integer(0), gov_type = character(0),
    travel_time_cities = numeric(0), travel_time_ports = numeric(0)
  )
  write_atomic(empty_grid, output_parquet)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

log_message(sprintf("TRACE: GADM bbox after transform: xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f",
                    st_bbox(gadm)["xmin"], st_bbox(gadm)["xmax"],
                    st_bbox(gadm)["ymin"], st_bbox(gadm)["ymax"]))

# Filter out GADM geometries with bounding boxes far outside tile extent
# (catches corrupt geometries or distant territories like Pacific islands)
log_message("TRACE: Filtering GADM geometries by bbox...")
tile_bbox <- st_bbox(tile_extent)
bbox_buffer <- 10  # degrees buffer

gadm_valid_idx <- sapply(seq_len(nrow(gadm)), function(i) {
  bb <- st_bbox(gadm[i, ])
  # Check if bbox is remotely near tile (within buffer)
  bb["xmin"] < (tile_bbox["xmax"] + bbox_buffer) &&
  bb["xmax"] > (tile_bbox["xmin"] - bbox_buffer) &&
  bb["ymin"] < (tile_bbox["ymax"] + bbox_buffer) &&
  bb["ymax"] > (tile_bbox["ymin"] - bbox_buffer)
})

n_filtered <- sum(!gadm_valid_idx)
if (n_filtered > 0) {
  log_message(sprintf("TRACE: Filtered %d GADM geometries with bounding boxes outside tile extent", n_filtered))
  gadm <- gadm[gadm_valid_idx, ]
} else {
  log_message("TRACE: All GADM geometries passed bbox filter")
}

if (nrow(gadm) == 0) {
  log_message("No valid GADM geometries after bbox filter. Creating empty output.")
  empty_grid <- data.table(
    grid_id = character(0), sub_tile_id = integer(0), tile_id = integer(0),
    country_iso3 = character(0), country_name = character(0),
    centroid_lon = numeric(0), centroid_lat = numeric(0), area_km2 = numeric(0),
    is_protected = logical(0), wdpa_pid = integer(0), iucn_cat = character(0),
    desig_year = integer(0), gov_type = character(0),
    travel_time_cities = numeric(0), travel_time_ports = numeric(0)
  )
  write_atomic(empty_grid, output_parquet)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# ==============================================================================
# LOAD WDPA (for this sub-tile extent only)
# ==============================================================================

log_message("Loading WDPA protected areas...")

log_message(sprintf("TRACE: Reading WDPA from: %s", wdpa_clean_path))
wdpa <- tryCatch({
  log_message("TRACE: Attempting st_read with wkt_filter...")
  result <- st_read(wdpa_clean_path, quiet = TRUE, wkt_filter = tile_extent_wkt)
  log_message(sprintf("TRACE: st_read with wkt_filter succeeded, got %d rows", nrow(result)))
  result
}, error = function(e) {
  log_message(sprintf("TRACE: Spatial filter failed with error: %s", e$message))
  log_message("TRACE: Using manual filter instead...")
  wdpa_full <- st_read(wdpa_clean_path, quiet = TRUE)
  log_message(sprintf("TRACE: Loaded full WDPA with %d rows", nrow(wdpa_full)))
  sf_use_s2(FALSE)
  intersects <- st_intersects(wdpa_full, tile_extent, sparse = FALSE)[, 1]
  sf_use_s2(TRUE)
  log_message(sprintf("TRACE: Manual filter found %d intersecting PAs", sum(intersects)))
  wdpa_full[intersects, ]
})

log_message(sprintf("TRACE: Loaded %d protected areas intersecting sub-tile", nrow(wdpa)))

# Standardize WDPA column names
log_message("TRACE: Standardizing WDPA column names...")
if ("WDPAID" %in% names(wdpa)) {
  wdpa$wdpa_pid <- as.integer(wdpa$WDPAID)
} else if ("WDPA_PID" %in% names(wdpa)) {
  wdpa$wdpa_pid <- as.integer(wdpa$WDPA_PID)
} else {
  wdpa$wdpa_pid <- seq_len(nrow(wdpa))
}

if ("IUCN_CAT" %in% names(wdpa)) wdpa$iucn_cat <- wdpa$IUCN_CAT
if ("STATUS_YR" %in% names(wdpa)) wdpa$desig_year <- as.integer(wdpa$STATUS_YR)
if ("GOV_TYPE" %in% names(wdpa)) wdpa$gov_type <- wdpa$GOV_TYPE

# Keep only needed columns
wdpa_cols <- intersect(c("wdpa_pid", "iucn_cat", "desig_year", "gov_type"), names(wdpa))
log_message(sprintf("TRACE: Keeping WDPA columns: %s", paste(wdpa_cols, collapse = ", ")))
wdpa <- wdpa[, wdpa_cols]
log_message(sprintf("TRACE: WDPA after column filter: %d rows", nrow(wdpa)))

if (nrow(wdpa) > 0) {
  log_message("TRACE: Calling make_valid_safe on WDPA...")
  wdpa <- make_valid_safe(wdpa)

  # Drop any geometries still invalid after repair (make_valid_safe warns but doesn't drop)
  still_invalid <- !st_is_valid(wdpa)
  if (any(still_invalid)) {
    log_message(sprintf("TRACE: Dropping %d still-invalid WDPA geometries", sum(still_invalid)))
    wdpa <- wdpa[!still_invalid, ]
  }

  log_message("TRACE: Transforming WDPA to WGS84...")
  wdpa <- st_transform(wdpa, WGS84_CRS)
  log_message(sprintf("TRACE: WDPA bbox after transform: xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f",
                      st_bbox(wdpa)["xmin"], st_bbox(wdpa)["xmax"],
                      st_bbox(wdpa)["ymin"], st_bbox(wdpa)["ymax"]))

  # Filter out corrupt geometries:
  # 1. Geometries with absurdly large bboxes (spanning > 90 degrees)
  # 2. Geometries with bboxes far outside tile extent
  tile_bbox <- st_bbox(tile_extent)
  max_bbox_span <- 90  # degrees - no valid PA should span more than this

  wdpa_valid_idx <- sapply(seq_len(nrow(wdpa)), function(i) {
    bb <- st_bbox(wdpa[i, ])
    lon_span <- bb["xmax"] - bb["xmin"]
    lat_span <- bb["ymax"] - bb["ymin"]

    # Reject if bbox spans too much (corrupt geometry)
    if (lon_span > max_bbox_span || lat_span > max_bbox_span) {
      return(FALSE)
    }

    # Reject if bbox doesn't overlap tile (with buffer)
    bbox_buffer <- 10
    bb["xmin"] < (tile_bbox["xmax"] + bbox_buffer) &&
    bb["xmax"] > (tile_bbox["xmin"] - bbox_buffer) &&
    bb["ymin"] < (tile_bbox["ymax"] + bbox_buffer) &&
    bb["ymax"] > (tile_bbox["ymin"] - bbox_buffer)
  })

  n_filtered <- sum(!wdpa_valid_idx)
  if (n_filtered > 0) {
    log_message(sprintf("Filtered %d WDPA geometries with invalid/oversized bounding boxes", n_filtered))
    wdpa <- wdpa[wdpa_valid_idx, ]
  }
}

# ==============================================================================
# CLIP TO SUB-TILE EXTENT
# ==============================================================================

log_message("Clipping to sub-tile extent...")

sf_use_s2(FALSE)

log_message("TRACE: About to st_crop GADM...")
log_message(sprintf("TRACE: GADM has %d rows, tile_extent bbox: [%.2f, %.2f] x [%.2f, %.2f]",
                    nrow(gadm),
                    st_bbox(tile_extent)["xmin"], st_bbox(tile_extent)["xmax"],
                    st_bbox(tile_extent)["ymin"], st_bbox(tile_extent)["ymax"]))
gadm_clipped <- st_crop(gadm, tile_extent)
log_message("TRACE: st_crop GADM succeeded")
gadm_clipped <- make_valid_safe(gadm_clipped)
log_message("TRACE: make_valid_safe on gadm_clipped succeeded")

if (nrow(wdpa) > 0) {
  log_message("TRACE: About to st_crop WDPA...")
  log_message(sprintf("TRACE: WDPA has %d rows", nrow(wdpa)))
  wdpa_clipped <- st_crop(wdpa, tile_extent)
  log_message("TRACE: st_crop WDPA succeeded")
  wdpa_clipped <- make_valid_safe(wdpa_clipped)
  log_message("TRACE: make_valid_safe on wdpa_clipped succeeded")
} else {
  wdpa_clipped <- wdpa
}

sf_use_s2(TRUE)

if (nrow(gadm_clipped) == 0) {
  log_message("No land after clipping. Creating empty output.")
  empty_grid <- data.table(
    grid_id = character(0), sub_tile_id = integer(0), tile_id = integer(0),
    country_iso3 = character(0), country_name = character(0),
    centroid_lon = numeric(0), centroid_lat = numeric(0), area_km2 = numeric(0),
    is_protected = logical(0), wdpa_pid = integer(0), iucn_cat = character(0),
    desig_year = integer(0), gov_type = character(0),
    travel_time_cities = numeric(0), travel_time_ports = numeric(0)
  )
  write_atomic(empty_grid, output_parquet)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# ==============================================================================
# PROJECT TO MOLLWEIDE
# ==============================================================================

log_message("Projecting to Mollweide equal-area...")

log_message("TRACE: Transforming gadm_clipped to Mollweide...")
gadm_moll <- st_transform(gadm_clipped, MOLLWEIDE_CRS)
log_message("TRACE: st_transform gadm_clipped succeeded")
gadm_moll <- make_valid_safe(gadm_moll)
log_message("TRACE: make_valid_safe gadm_moll succeeded")
tile_extent_moll <- st_transform(tile_extent, MOLLWEIDE_CRS)
log_message("TRACE: tile_extent_moll created")

if (nrow(wdpa_clipped) > 0) {
  log_message("TRACE: Transforming wdpa_clipped to Mollweide...")
  wdpa_moll <- st_transform(wdpa_clipped, MOLLWEIDE_CRS)
  log_message("TRACE: st_transform wdpa_clipped succeeded")
  wdpa_moll <- make_valid_safe(wdpa_moll)
  log_message("TRACE: make_valid_safe wdpa_moll succeeded")
}

# ==============================================================================
# CREATE GRID
# ==============================================================================

log_message(sprintf("Creating %dm x %dm grid...", GRID_CELLSIZE_M, GRID_CELLSIZE_M))

log_message("TRACE: Calling st_make_grid...")
grid_raw <- st_make_grid(
  tile_extent_moll,
  cellsize = c(GRID_CELLSIZE_M, GRID_CELLSIZE_M),
  what = "polygons"
)

log_message(sprintf("TRACE: Created %d raw grid cells", length(grid_raw)))

log_message("TRACE: Creating grid_sf from grid_raw...")
grid_sf <- st_sf(cell_idx = seq_along(grid_raw), geometry = grid_raw)
log_message(sprintf("TRACE: grid_sf created with %d rows", nrow(grid_sf)))

# ==============================================================================
# INTERSECT WITH COUNTRY BOUNDARIES
# ==============================================================================

log_message("Intersecting grid with country boundaries...")

sf_use_s2(FALSE)
log_message("TRACE: sf_use_s2 set to FALSE for country intersection")

log_message(sprintf("TRACE: About to st_intersection grid_sf (%d rows) with gadm_moll (%d rows)...",
                    nrow(grid_sf), nrow(gadm_moll)))
grid_country <- st_intersection(grid_sf, gadm_moll)
log_message("TRACE: st_intersection grid_sf with gadm_moll succeeded")
grid_country <- make_valid_safe(grid_country)
log_message("TRACE: make_valid_safe grid_country succeeded")

sf_use_s2(TRUE)
log_message("TRACE: sf_use_s2 set back to TRUE")

log_message(sprintf("TRACE: Created %d grid-country cells", nrow(grid_country)))

if (nrow(grid_country) == 0) {
  log_message("No grid cells after country intersection. Creating empty output.")
  empty_grid <- data.table(
    grid_id = character(0), sub_tile_id = integer(0), tile_id = integer(0),
    country_iso3 = character(0), country_name = character(0),
    centroid_lon = numeric(0), centroid_lat = numeric(0), area_km2 = numeric(0),
    is_protected = logical(0), wdpa_pid = integer(0), iucn_cat = character(0),
    desig_year = integer(0), gov_type = character(0),
    travel_time_cities = numeric(0), travel_time_ports = numeric(0)
  )
  write_atomic(empty_grid, output_parquet)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# ==============================================================================
# INTERSECT WITH WDPA BOUNDARIES
# ==============================================================================

if (nrow(wdpa_clipped) > 0) {
  log_message("Intersecting grid with WDPA boundaries...")

  sf_use_s2(FALSE)
  log_message("TRACE: sf_use_s2 set to FALSE for WDPA intersection")

  # Get grid cells that intersect with protected areas
  log_message(sprintf("TRACE: About to st_intersection grid_country (%d rows) with wdpa_moll (%d rows)...",
                      nrow(grid_country), nrow(wdpa_moll)))
  grid_in_pa <- st_intersection(grid_country, wdpa_moll)
  log_message(sprintf("TRACE: st_intersection succeeded, grid_in_pa has %d rows", nrow(grid_in_pa)))
  grid_in_pa <- make_valid_safe(grid_in_pa)
  log_message("TRACE: make_valid_safe grid_in_pa succeeded")

  # Only assign values if there are rows (avoid 0-row assignment error)
  if (nrow(grid_in_pa) > 0) {
    grid_in_pa$is_protected <- TRUE
  }

  # Get grid cells outside protected areas (difference)
  # This is the tricky part - we need the parts of grid_country NOT in any PA
  log_message("TRACE: Creating wdpa_union with st_union...")

  # Simplify WDPA geometries to remove degenerate vertices (within 1m tolerance)
  # This fixes TopologyException errors from zero-length segments
  log_message("TRACE: Simplifying wdpa_moll to remove degenerate vertices...")
  wdpa_moll_clean <- st_simplify(wdpa_moll, dTolerance = 1, preserveTopology = TRUE)
  wdpa_moll_clean <- wdpa_moll_clean[!st_is_empty(wdpa_moll_clean), ]
  log_message(sprintf("TRACE: After st_simplify: %d rows (was %d)", nrow(wdpa_moll_clean), nrow(wdpa_moll)))

  # Re-validate after simplify (can create holes outside shells)
  wdpa_moll_clean <- make_valid_safe(wdpa_moll_clean)
  log_message("TRACE: make_valid_safe after st_simplify succeeded")

  wdpa_union <- st_union(wdpa_moll_clean)
  log_message("TRACE: st_union succeeded")
  wdpa_union <- make_valid_safe(st_sf(geometry = wdpa_union))
  log_message("TRACE: make_valid_safe wdpa_union succeeded")

  log_message(sprintf("TRACE: About to st_difference grid_country (%d rows) with wdpa_union...",
                      nrow(grid_country)))
  grid_outside_pa <- st_difference(grid_country, wdpa_union)
  log_message(sprintf("TRACE: st_difference succeeded, grid_outside_pa has %d rows", nrow(grid_outside_pa)))
  grid_outside_pa <- make_valid_safe(grid_outside_pa)
  log_message("TRACE: make_valid_safe grid_outside_pa succeeded")

  # Only assign values if there are rows (avoid 0-row assignment error)
  if (nrow(grid_outside_pa) > 0) {
    grid_outside_pa$is_protected <- FALSE
    grid_outside_pa$wdpa_pid <- NA_integer_
    grid_outside_pa$iucn_cat <- NA_character_
    grid_outside_pa$desig_year <- NA_integer_
    grid_outside_pa$gov_type <- NA_character_
  }

  sf_use_s2(TRUE)
  log_message("TRACE: sf_use_s2 set back to TRUE")

  # Combine protected and unprotected cells
  # Handle edge cases where one is empty
  log_message(sprintf("TRACE: About to combine protected (%d rows) and unprotected (%d rows) cells",
                      nrow(grid_in_pa), nrow(grid_outside_pa)))

  if (nrow(grid_in_pa) == 0 && nrow(grid_outside_pa) == 0) {
    # Both empty - use grid_country with is_protected = FALSE
    grid_combined <- grid_country
    grid_combined$is_protected <- FALSE
    grid_combined$wdpa_pid <- NA_integer_
    grid_combined$iucn_cat <- NA_character_
    grid_combined$desig_year <- NA_integer_
    grid_combined$gov_type <- NA_character_
  } else if (nrow(grid_in_pa) == 0) {
    # Only unprotected cells
    grid_combined <- grid_outside_pa
  } else if (nrow(grid_outside_pa) == 0) {
    # Only protected cells
    grid_combined <- grid_in_pa
  } else {
    # Both have rows - rbind them
    log_message(sprintf("TRACE: grid_in_pa geom types: %s",
                        paste(unique(st_geometry_type(grid_in_pa)), collapse = ", ")))
    log_message(sprintf("TRACE: grid_outside_pa geom types: %s",
                        paste(unique(st_geometry_type(grid_outside_pa)), collapse = ", ")))
    common_cols <- intersect(names(grid_in_pa), names(grid_outside_pa))
    grid_combined <- rbind(
      grid_in_pa[, common_cols],
      grid_outside_pa[, common_cols]
    )
  }

  log_message(sprintf("TRACE: After combine: %d rows, class: %s",
                      nrow(grid_combined), paste(class(grid_combined), collapse = ", ")))
  log_message(sprintf("TRACE: grid_combined geom types: %s",
                      paste(unique(st_geometry_type(grid_combined)), collapse = ", ")))

  log_message(sprintf("After WDPA split: %d protected cells, %d unprotected cells",
                      sum(grid_combined$is_protected),
                      sum(!grid_combined$is_protected)))

  # Handle cells that overlap multiple PAs - keep oldest (earliest desig_year)
  if ("desig_year" %in% names(grid_combined) && any(grid_combined$is_protected)) {
    # Create a unique key for each original cell + country combo
    grid_combined$temp_key <- paste(grid_combined$cell_idx, grid_combined$GID_0, sep = "_")

    # For protected cells, if multiple PAs, keep the one with earliest desig_year
    grid_dt <- as.data.table(grid_combined)

    # Find cells with multiple PA overlaps
    protected_counts <- grid_dt[is_protected == TRUE, .N, by = temp_key]
    multi_pa <- protected_counts[N > 1]$temp_key

    if (length(multi_pa) > 0) {
      log_message(sprintf("Resolving %d cells with multiple PA overlaps (using oldest PA)",
                          length(multi_pa)))

      # Use base R operations directly on sf (preserves geometry correctly)
      grid_combined$desig_year[is.na(grid_combined$desig_year)] <- 1900L
      grid_combined <- grid_combined[order(grid_combined$temp_key, grid_combined$desig_year), ]
      grid_combined <- grid_combined[!duplicated(grid_combined$temp_key), ]
      grid_combined$desig_year[grid_combined$desig_year == 1900L] <- NA_integer_

      log_message(sprintf("TRACE: After deduplication: %d rows", nrow(grid_combined)))
    }

    grid_combined$temp_key <- NULL
  }

  log_message(sprintf("TRACE: grid_combined before assignment to grid_final: %d rows", nrow(grid_combined)))
  grid_final <- grid_combined
  log_message(sprintf("TRACE: grid_final assigned: %d rows, geometry types: %s",
                      nrow(grid_final),
                      paste(unique(st_geometry_type(grid_final)), collapse = ", ")))

} else {
  # No protected areas in this sub-tile
  log_message("No protected areas in this sub-tile")
  grid_final <- grid_country
  grid_final$is_protected <- FALSE
  grid_final$wdpa_pid <- NA_integer_
  grid_final$iucn_cat <- NA_character_
  grid_final$desig_year <- NA_integer_
  grid_final$gov_type <- NA_character_
}

# ==============================================================================
# CALCULATE ATTRIBUTES
# ==============================================================================

log_message("Calculating cell attributes...")

log_message(sprintf("TRACE: grid_final has %d rows, %d columns", nrow(grid_final), ncol(grid_final)))
log_message(sprintf("TRACE: grid_final geometry type: %s", paste(unique(st_geometry_type(grid_final)), collapse = ", ")))
log_message(sprintf("TRACE: grid_final class: %s", paste(class(grid_final), collapse = ", ")))

# Verify grid_final is an sf object
if (!inherits(grid_final, "sf")) {
  log_message("TRACE: grid_final is not sf, converting...")
  grid_final <- st_as_sf(grid_final)
}

# Check for geometry column attributes
log_message(sprintf("TRACE: geometry column name: %s", attr(grid_final, "sf_column")))

sf_use_s2(FALSE)
log_message("TRACE: sf_use_s2 set to FALSE")

log_message("TRACE: About to call st_area...")
grid_final$area_m2 <- as.numeric(st_area(grid_final))
log_message("TRACE: st_area succeeded")

grid_final$area_km2 <- grid_final$area_m2 / 1e6
log_message("TRACE: area_km2 calculated")

# Calculate centroids
log_message("TRACE: About to call st_centroid...")
centroids_moll <- st_centroid(grid_final)
log_message("TRACE: st_centroid succeeded")

log_message("TRACE: About to transform centroids to WGS84...")
centroids_wgs84 <- st_transform(centroids_moll, WGS84_CRS)
log_message("TRACE: st_transform centroids succeeded")

log_message("TRACE: About to call st_coordinates...")
coords <- st_coordinates(centroids_wgs84)
log_message("TRACE: st_coordinates succeeded")

sf_use_s2(TRUE)
log_message("TRACE: sf_use_s2 set back to TRUE")

grid_final$centroid_lon <- coords[, "X"]
grid_final$centroid_lat <- coords[, "Y"]

# ==============================================================================
# EXTRACT TRANSPORT COSTS
# ==============================================================================

log_message("Extracting transport cost rasters...")

# Need grid in WGS84 for raster extraction (rasters are typically WGS84)
log_message("TRACE: Transforming grid_final to WGS84 for raster extraction...")
grid_wgs84 <- st_transform(grid_final, WGS84_CRS)
log_message(sprintf("TRACE: grid_wgs84 created with %d rows", nrow(grid_wgs84)))

# GEOMETRYCOLLECTION HANDLING: exact_extract doesn't support mixed geometry types
# Extract polygon components from GEOMETRYCOLLECTIONs (don't drop them!)
geom_types <- st_geometry_type(grid_wgs84)
geom_type_table <- table(geom_types)
log_message("TRACE: Geometry type distribution before normalization:")
for (gtype in names(geom_type_table)) {
  log_message(sprintf("  %s: %d", gtype, geom_type_table[[gtype]]))
}

# Extract polygon components from GEOMETRYCOLLECTIONs
gc_idx <- which(geom_types == "GEOMETRYCOLLECTION")
if (length(gc_idx) > 0) {
  log_message(sprintf("TRACE: Extracting polygons from %d GEOMETRYCOLLECTIONs", length(gc_idx)))
  # Extract polygon components (keeps POLYGON and MULTIPOLYGON parts)
  grid_wgs84 <- st_collection_extract(grid_wgs84, "POLYGON")
  grid_final <- st_collection_extract(grid_final, "POLYGON")
  log_message(sprintf("TRACE: After st_collection_extract: %d rows", nrow(grid_wgs84)))
}

# Now drop any remaining non-polygon geometries (points, lines from edge artifacts)
geom_types <- st_geometry_type(grid_wgs84)
non_polygon_idx <- which(!geom_types %in% c("POLYGON", "MULTIPOLYGON"))
if (length(non_polygon_idx) > 0) {
  log_message(sprintf("TRACE: Dropping %d non-polygon geometries (points/lines)", length(non_polygon_idx)))
  grid_wgs84 <- grid_wgs84[-non_polygon_idx, ]
  grid_final <- grid_final[-non_polygon_idx, ]
}

# Cast all to MULTIPOLYGON for consistency
grid_wgs84 <- st_cast(grid_wgs84, "MULTIPOLYGON")
log_message(sprintf("TRACE: After st_cast to MULTIPOLYGON: %d rows", nrow(grid_wgs84)))

# Verify geometry types are now uniform
geom_types_after <- unique(st_geometry_type(grid_wgs84))
log_message(sprintf("TRACE: Geometry types after normalization: %s", paste(geom_types_after, collapse = ", ")))

# Travel time to cities
if (file.exists(travel_cities_path)) {
  log_message("TRACE: Loading cities raster...")
  cities_rast <- terra::rast(travel_cities_path)
  log_message("TRACE: Extracting travel time to cities...")
  grid_final$travel_time_cities <- exact_extract(cities_rast, grid_wgs84, "mean")
  log_message(sprintf("TRACE: Travel time to cities extracted (mean = %.1f min)",
                      mean(grid_final$travel_time_cities, na.rm = TRUE)))
} else {
  log_message(sprintf("TRACE: WARNING: Cities raster not found: %s", travel_cities_path))
  grid_final$travel_time_cities <- NA_real_
}

# Travel time to ports
if (file.exists(travel_ports_path)) {
  log_message("TRACE: Loading ports raster...")
  ports_rast <- terra::rast(travel_ports_path)
  log_message("TRACE: Extracting travel time to ports...")
  grid_final$travel_time_ports <- exact_extract(ports_rast, grid_wgs84, "mean")
  log_message(sprintf("TRACE: Travel time to ports extracted (mean = %.1f min)",
                      mean(grid_final$travel_time_ports, na.rm = TRUE)))
} else {
  log_message(sprintf("TRACE: WARNING: Ports raster not found: %s", travel_ports_path))
  grid_final$travel_time_ports <- NA_real_
}

log_message("TRACE: Cleaning up grid_wgs84...")
rm(grid_wgs84)
gc()
log_message("TRACE: Transport cost extraction complete")

# ==============================================================================
# CREATE UNIQUE GRID IDs
# ==============================================================================

log_message("Creating unique grid IDs...")

# Grid ID format: t{tile}_s{subtile}_{P/U}_{seq}
# P = protected, U = unprotected
log_message("TRACE: Creating prot_flag...")
prot_flag <- ifelse(grid_final$is_protected, "P", "U")
log_message("TRACE: Creating grid_id column...")
grid_final$grid_id <- sprintf("t%03d_s%04d_%s_%07d",
                               sub_tile$tmf_tile_id,
                               sub_tile_id,
                               prot_flag,
                               seq_len(nrow(grid_final)))
log_message(sprintf("TRACE: grid_id created, sample: %s", grid_final$grid_id[1]))

grid_final$sub_tile_id <- sub_tile_id
grid_final$tile_id <- sub_tile$tmf_tile_id

# Standardize country columns
log_message("TRACE: Standardizing country columns...")
if ("GID_0" %in% names(grid_final)) grid_final$country_iso3 <- grid_final$GID_0
if ("COUNTRY" %in% names(grid_final)) grid_final$country_name <- grid_final$COUNTRY
log_message("TRACE: Country columns standardized")

# ==============================================================================
# PREPARE OUTPUT
# ==============================================================================

log_message("Preparing output...")

output_cols <- c("grid_id", "sub_tile_id", "tile_id",
                 "country_iso3", "country_name",
                 "centroid_lon", "centroid_lat", "area_km2",
                 "is_protected", "wdpa_pid", "iucn_cat", "desig_year", "gov_type",
                 "travel_time_cities", "travel_time_ports")

log_message(sprintf("TRACE: grid_final has columns: %s", paste(names(grid_final), collapse = ", ")))

# Ensure all columns exist
log_message("TRACE: Ensuring all output columns exist...")
for (col in output_cols) {
  if (!col %in% names(grid_final)) {
    log_message(sprintf("TRACE: Adding missing column: %s", col))
    if (col == "is_protected") {
      grid_final[[col]] <- FALSE
    } else if (col %in% c("wdpa_pid", "desig_year")) {
      grid_final[[col]] <- NA_integer_
    } else if (col %in% c("travel_time_cities", "travel_time_ports")) {
      grid_final[[col]] <- NA_real_
    } else {
      grid_final[[col]] <- NA_character_
    }
  }
}
log_message("TRACE: All output columns ensured")

log_message("TRACE: Calling st_drop_geometry and converting to data.table...")
grid_attrs <- as.data.table(st_drop_geometry(grid_final))
log_message(sprintf("TRACE: grid_attrs created with %d rows, %d columns", nrow(grid_attrs), ncol(grid_attrs)))
grid_attrs <- grid_attrs[, ..output_cols]
log_message(sprintf("TRACE: grid_attrs filtered to %d output columns", ncol(grid_attrs)))

log_message(sprintf("TRACE: Final grid: %d cells (%d protected, %d unprotected)",
                    nrow(grid_attrs),
                    sum(grid_attrs$is_protected),
                    sum(!grid_attrs$is_protected)))

# ==============================================================================
# WRITE OUTPUT
# ==============================================================================

log_message(sprintf("TRACE: Writing parquet: %s", output_parquet))
write_atomic(grid_attrs, output_parquet)
log_message("TRACE: Parquet write succeeded")

# Write geopackage with geometry
log_message("TRACE: Selecting output columns from grid_final...")
grid_output <- grid_final[, output_cols]
log_message(sprintf("TRACE: grid_output has %d rows", nrow(grid_output)))
log_message("TRACE: Transforming grid_output to WGS84...")
grid_output_wgs84 <- st_transform(grid_output, WGS84_CRS)
log_message("TRACE: st_transform grid_output succeeded")

log_message(sprintf("TRACE: Writing geopackage: %s", output_gpkg))
st_write(grid_output_wgs84, output_gpkg, delete_dsn = TRUE, quiet = TRUE)
log_message("TRACE: Geopackage write succeeded")

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

log_message("TRACE: Computing summary statistics...")
log_message("Summary statistics:")
log_message(sprintf("  Total cells: %d", nrow(grid_attrs)))
log_message(sprintf("  Protected cells: %d (%.1f%%)",
                    sum(grid_attrs$is_protected),
                    100 * mean(grid_attrs$is_protected)))
log_message(sprintf("  Countries: %s",
                    paste(unique(grid_attrs$country_iso3), collapse = ", ")))
log_message(sprintf("  Total area: %.2f km^2", sum(grid_attrs$area_km2)))
log_message(sprintf("  Mean cell area: %.4f km^2", mean(grid_attrs$area_km2)))

if (sum(grid_attrs$is_protected) > 0) {
  iucn_dist <- grid_attrs[is_protected == TRUE, .N, by = iucn_cat]
  setorder(iucn_dist, -N)
  log_message("IUCN category distribution:")
  for (i in seq_len(min(nrow(iucn_dist), 5))) {
    log_message(sprintf("    %s: %d", iucn_dist$iucn_cat[i], iucn_dist$N[i]))
  }
}

log_message("TRACE: Final cleanup...")
gc_verbose()
log_message("TRACE: Job complete!")
log_job_end(start_time)
