# ==============================================================================
# TMF 5KM PIPELINE - SHARED CONFIGURATION
# ==============================================================================
# This file contains all shared configuration, paths, and utility functions
# for the Tropical Moist Forest data pipeline.
#
# Grid cells are 5km x 5km and cut on country + WDPA boundaries
# so each cell is fully within one country and either fully protected or not.
#
# Source this file at the beginning of every pipeline script:
#   source("code/build/BUILD_workspace.R")
# ==============================================================================

library(here)
here::i_am("code/build/BUILD_workspace.R")

# Load required packages
suppressMessages({
  library(sf)
  library(terra)
  library(exactextractr)
  library(data.table)
  library(arrow)
  library(dplyr)
  library(stringr)
  library(units)
  library(RANN)
  library(foreign)  # For reading DBF files
})

# ==============================================================================
# DATA PATHS - CONFIGURE THESE FOR YOUR SYSTEM
# ==============================================================================

# Project root (auto-detected via here package)
project_root <- here()

# ------------------------------------------------------------------------------
# INPUT DATA
#
# The big shared rasters live in the lab's Global Forest Repo tree, not in this
# project. Point at them where they actually are rather than copying ~49 GB.
# Use the space-free GlobalForest symlink: GDAL handles the space in
# "Global Forest Repo" badly. Override FRAG_SPATIAL if that tree moves.
# ------------------------------------------------------------------------------
spatial_path <- Sys.getenv(
  "FRAG_SPATIAL",
  "/resnick/groups/MishraLab/GlobalForest/data/raw/spatial"
)

# Project-local paths
raw_data_path <- file.path(project_root, "Data", "raw")
build_data_path <- file.path(project_root, "Data", "build")

# TMF lives in this project's own tree (120 GB, staged 2026-08-18), not in
# the shared spatial/ tree - same convention as GFED5 below.
tmf_path  <- file.path(raw_data_path, 'raster', 'TMF', 'AnnualChange')  # TMF GeoTIFF tiles
gadm_path <- file.path(spatial_path, 'GADM', 'gadm_410-levels.gpkg')    # GADM country boundaries
wdpa_path <- file.path(spatial_path, 'WDPA', 'WDPA_Oct2021',            # WDPA GDB
                       'WDPA_WDOECM_Oct2021_Public_all.gdb')

# Transport cost rasters (static metadata)
# The Nelson (2015) travel-time rasters are in the nelson/ subdirectory; the
# _12 / _5 variants are the ones this pipeline uses.
transport_path <- file.path(spatial_path, 'transpocosts', 'nelson')
travel_cities_path <- file.path(transport_path, 'travel_time_to_cities_12.tif')
travel_ports_path <- file.path(transport_path, 'travel_time_to_ports_5.tif')

# EarthStat yields + harvested-area (crop share), one directory per crop.
# Pointed at the real location rather than a symlink: a symlink created over
# the SMB mount from macOS is unreadable on the cluster ("cannot read symbolic
# link: Input/output error"), which silently produced 303 empty yield outputs.
earthstat_path <- file.path(spatial_path, 'landuse', 'EarthStat')
yields_path <- file.path(earthstat_path,
                         'HarvestedAreaYield175Crops_Geotiff', 'GeoTiff')

# EarthStat 5-arcmin cropland / pasture fractions (Ramankutty et al. 2008) - the
# layer ABOVE the 172 per-crop grids. Needed because HarvestedAreaFraction is a
# fraction of CELL area and the 172 of them can sum past the cropland fraction:
# multi-cropping counts one hectare once per harvest. cropland_frac is therefore
# the denominator that makes "share of the cell's cropland in crop c" a
# well-defined object, and pasture_frac carries the cattle margin that the crop
# layers miss entirely - the dominant tropical-forest land use in Latin America.
cropland_raster_path <- file.path(earthstat_path, 'Cropland2000_5m.tif')
pasture_raster_path  <- file.path(earthstat_path, 'Pasture2000_5m.tif')

# ------------------------------------------------------------------------------
# Non-spatial economic inputs, shared with the Global Forest Repo. spatial_path
# points at .../data/raw/spatial, so non-spatial is its sibling. Read-only here:
# this project consumes those files, the companion owns them.
# ------------------------------------------------------------------------------
gf_nonspatial_path <- file.path(dirname(spatial_path), 'non-spatial')

# FAOSTAT producer prices. 2022 vintage: annual, Y1991-Y2021, and it uses the
# PRE-2021-CPC item names ("Rice, paddy", "Oil, palm"), which is why the crop
# crosswalk carries fao_item_legacy alongside fao_item_current. Join on legacy.
fao_prices_path <- file.path(gf_nonspatial_path, 'prices',
                             'Prices_E_All_Data_NOFLAG.csv')

# Country trucking price, USD per ton-km, read off Figures 1.6/1.7 of the World
# Bank's "Shrinking Economic Distance". Coarse by construction - it is a
# country-level scalar, so it moves no within-country variation.
trucking_cost_path <- file.path(gf_nonspatial_path, 'trade costs',
                                'SED Trucking Cost', 'worldbank_chatgpt_extract.csv')

# Hand-maintained lookups owned by THIS repo (version-controlled, not derived).
lookup_path <- file.path(project_root, "Data", "lookup")
crop_crosswalk_file <- file.path(lookup_path, "crop_crosswalk_earthstat_fao.csv")

# Derived economic tables written by 0e_prepare_crop_prices.R
crop_prices_file    <- file.path(lookup_path, "fao_crop_prices.parquet")
crop_price_pre_file <- file.path(lookup_path, "crop_price_preperiod.parquet")
crop_price_glb_file <- file.path(lookup_path, "crop_price_global.parquet")
trucking_file       <- file.path(lookup_path, "trucking_cost.parquet")

# Covariate data paths (static covariates for Stage 5)
cities_path <- file.path(spatial_path, "cities", "cityLocations.csv")
# The GHS_POP_E#### folders sit inside GHSL/, and the extractor scans exactly
# one level deep for them - so this points at GHSL, not its parent.
population_path <- file.path(spatial_path, "population", "GHSL")
biomass_path <- file.path(spatial_path, "biomass")

# Biodiversity: species richness + range-size rarity. Combined_THR_SR_2022
# (threatened richness) and RSR_CRENVU (range-size rarity for CR/EN/VU species)
# proxy conservation value - i.e. why a PA was sited where it was.
biodiversity_path <- file.path(spatial_path, "biodiversity")

# Livestock density, GLW3 2010. Cattle under three allocation methods; the
# dasymetric (_Da) layer is the one to prefer.
livestock_path <- file.path(spatial_path, "pasture", "GLW3")
# DEAD as of 2026-08-20: the natcap bucket revoked anonymous read and this URL
# returns HTTP 403 for everyone. Elevation now comes from Copernicus DEM 30m,
# fetched per-tile by get_copernicus_dem() in 5_extract_covariates.R. Kept only
# so the failure is documented rather than silently rediscovered.
# aster_url <- "/vsicurl/https://storage.googleapis.com/natcap-data-cache/global/aster-v3-1s/aster-v3-1s.tif"

# Specific output directories
grid_output_path <- file.path(build_data_path, "grids")
tmf_output_path <- file.path(build_data_path, "tmf")
consolidated_path <- file.path(build_data_path, "tmf_consolidated")
classification_path <- file.path(build_data_path, "classifications")
wdpa_output_path <- file.path(build_data_path, "wdpa")
yields_output_path <- file.path(build_data_path, "yields")
covariates_output_path <- file.path(build_data_path, "covariates")
final_output_path <- file.path(build_data_path, "final")
gfed_output_path <- file.path(build_data_path, "gfed")
gfed_consolidated_path <- file.path(build_data_path, "gfed_consolidated")

# GFED data path
gfed_data_path <- file.path(raw_data_path, "raster", "GFED5")

# Final outputs. TWO files, at two different grains, produced by two different
# stages - this is deliberate, not an accident of the build order:
#
#   panel_final_file     cell-YEAR   66.8M rows   Stage 6  (6_assemble_final.R)
#   cropland_final_file  cell        ~2M rows     Stage 1b (1b_assemble_cropland.R)
#
# EarthStat is a year-2000 snapshot, so folding its 345 crop columns into the
# cell-year panel would replicate a time-invariant object 34 times over (~184 GB
# uncompressed) to say nothing new. Analysis joins the two on grid_id.
panel_final_file    <- file.path(final_output_path, "TMF_5km_panel.parquet")
cropland_final_file <- file.path(final_output_path, "TMF_5km_cropland.parquet")

# Create directories if they don't exist
sapply(c(grid_output_path, tmf_output_path, consolidated_path,
         classification_path, wdpa_output_path, yields_output_path,
         covariates_output_path, final_output_path, gfed_output_path,
         gfed_consolidated_path, lookup_path),
       function(p) if(!dir.exists(p)) dir.create(p, recursive = TRUE))

# ==============================================================================
# PROJECTION DEFINITIONS
# ==============================================================================

# Mollweide equal-area projection (for grid creation and area calculations)
MOLLWEIDE_CRS <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# WGS84 geographic coordinates (for data storage and TMF rasters)
WGS84_CRS <- "+proj=longlat +datum=WGS84 +no_defs"

# ==============================================================================
# GRID PARAMETERS
# ==============================================================================

# Grid cell size in meters (5 km x 5 km)
GRID_CELLSIZE_M <- 5000
GRID_CELLSIZE_KM <- GRID_CELLSIZE_M / 1000

# Buffer distance for forest frontier calculation (100 km)
BUFFER_DISTANCE_M <- 100000

# Minimum cell area fraction - set to 0 to keep all fragments
# (grid is cut on country and WDPA boundaries)
MIN_CELL_AREA_FRAC <- 0

# ==============================================================================
# WDPA PARAMETERS (Protected Areas)
# ==============================================================================

# WDPA attributes to extract for each grid cell
# When a cell overlaps multiple PAs, use the oldest one (earliest STATUS_YR)
WDPA_COLUMNS <- c("WDPAID", "WDPA_PID", "NAME", "IUCN_CAT", "STATUS_YR", "GOV_TYPE")

# Grid output schema - columns written to grid files
GRID_SCHEMA <- c(
  # Identifiers
  "grid_id", "tile_id", "subtile_id",
  # Location
  "country_iso3", "country_name", "centroid_lon", "centroid_lat", "area_km2",
  # WDPA (baked into grid at creation)
  "is_protected", "wdpa_pid", "iucn_cat", "desig_year", "gov_type",
  # Transport costs (baked into grid at creation)
  "travel_time_cities", "travel_time_ports"
)

# ==============================================================================
# TMF DATA PARAMETERS
# ==============================================================================

# Years of TMF data availability
# Crop-price window. FAO prices are averaged over these years and held fixed,
# so this is also the risk-set cutoff for anything that treats price as
# pre-determined: 6_protection_logit.R sets DESIG_CUTOFF <- max(PRICE_PRE_YEARS).
# Defined HERE rather than in 0e_prepare_crop_prices.R because both the build
# and the analysis side need it, and both source this file - 6_protection_logit.R
# does not source 0e, so referencing it there failed with
# "object 'PRICE_PRE_YEARS' not found" after ~90s of setup work.
PRICE_PRE_YEARS <- 1991:2000

TMF_YEARS <- 1990:2023   # JRC TMF v1 AnnualChange: 34 rasters per tile
N_TMF_YEARS <- length(TMF_YEARS)

# Number of TMF tiles (computed from actual files after index is created)
# N_TMF_TILES is set after TMF_TILE_INDEX is created below

# TMF land cover classification legend
TMF_LEGEND <- data.table(
  value = c(1, 2, 3, 4, 5, 6),
  tmf_class = c("Undisturbed_TMF", "Degraded_TMF", "Deforested",
                "Regrowth", "Water", "Other_landcover")
)

# ==============================================================================
# TMF TILE INDEX
# ==============================================================================

# TMF tile definitions extracted from actual TMF files in the data directory
# Format: tile_id, lat_north, lat_south, lon_west, lon_east

create_tmf_tile_index <- function() {
  # List TMF tile folders
  # Format: JRC_TMF_AnnualChange_v1_{REGION}_ID{N}_{LAT}{VAL}_{LON}{VAL}
  # e.g., JRC_TMF_AnnualChange_v1_AFR_ID6_S20_E40

  tmf_dirs <- list.dirs(tmf_path, full.names = FALSE, recursive = FALSE)
  tmf_dirs <- tmf_dirs[grepl("^JRC_TMF_AnnualChange", tmf_dirs)]

  if (length(tmf_dirs) == 0) {
    stop(sprintf("No TMF tile folders found in: %s", tmf_path))
  }

  # Parse folder names to extract lat/lon
  # Pattern: ..._S20_E40 or ..._N10_W080
  parsed <- str_match(tmf_dirs, "_([NS])(\\d+)_([EW])(\\d+)$")

  if (any(is.na(parsed[, 1]))) {
    warning("Some TMF folders did not match expected pattern")
    valid <- !is.na(parsed[, 1])
    parsed <- parsed[valid, , drop = FALSE]
    tmf_dirs <- tmf_dirs[valid]
  }

  # Convert to numeric coordinates
  lat_sign <- ifelse(parsed[, 2] == "N", 1, -1)
  lon_sign <- ifelse(parsed[, 4] == "E", 1, -1)

  lat_north <- lat_sign * as.integer(parsed[, 3])
  lon_west <- lon_sign * as.integer(parsed[, 5])

  # Create data.table with unique tiles
  tiles <- data.table(
    lat_north = lat_north,
    lon_west = lon_west,
    folder = tmf_dirs
  )
  tiles <- unique(tiles, by = c("lat_north", "lon_west"))

  # Add derived columns
  tiles[, lat_south := lat_north - 10]
  tiles[, lon_east := lon_west + 10]

  # Sort by lat (descending) then lon (ascending) for consistent ordering
  setorder(tiles, -lat_north, lon_west)
  tiles[, tile_id := .I]

  # Reorder columns
  setcolorder(tiles,
              c("tile_id", "lat_north", "lat_south", "lon_west", "lon_east", "folder"))

  message(sprintf("Found %d unique TMF tiles", nrow(tiles)))

  tiles
}

TMF_TILE_INDEX <- create_tmf_tile_index()
N_TMF_TILES <- nrow(TMF_TILE_INDEX)

# ==============================================================================
# SUB-TILE INDEX (5° x 5° sub-tiles for parallelized grid creation)
# With 5km cells, each sub-tile has ~100x100 = 10,000 cells max
# ==============================================================================

SUB_TILE_SIZE <- 5  # degrees

#' Create index of 5x5 degree sub-tiles from 10x10 TMF tiles
create_sub_tile_index <- function() {
  sub_tiles_list <- lapply(seq_len(nrow(TMF_TILE_INDEX)), function(i) {
    tmf_tile <- TMF_TILE_INDEX[i]

    # Generate 5x5 sub-tiles within this 10x10 tile
    lat_starts <- seq(tmf_tile$lat_south, tmf_tile$lat_north - SUB_TILE_SIZE,
                      by = SUB_TILE_SIZE)
    lon_starts <- seq(tmf_tile$lon_west, tmf_tile$lon_east - SUB_TILE_SIZE,
                      by = SUB_TILE_SIZE)

    # Create all combinations
    grid <- expand.grid(lat_south = lat_starts, lon_west = lon_starts)

    data.table(
      tmf_tile_id = tmf_tile$tile_id,
      tmf_folder = tmf_tile$folder,
      lat_south = grid$lat_south,
      lat_north = grid$lat_south + SUB_TILE_SIZE,
      lon_west = grid$lon_west,
      lon_east = grid$lon_west + SUB_TILE_SIZE
    )
  })

  sub_tiles <- rbindlist(sub_tiles_list)
  setorder(sub_tiles, -lat_north, lon_west)
  sub_tiles[, sub_tile_id := .I]

  setcolorder(sub_tiles, c("sub_tile_id", "tmf_tile_id", "tmf_folder",
                            "lat_north", "lat_south", "lon_west", "lon_east"))

  message(sprintf("Created %d sub-tiles (%d° x %d°) from %d TMF tiles",
                  nrow(sub_tiles), SUB_TILE_SIZE, SUB_TILE_SIZE, N_TMF_TILES))

  sub_tiles
}

SUB_TILE_INDEX <- create_sub_tile_index()
N_SUB_TILES <- nrow(SUB_TILE_INDEX)
N_EXTRACTION_JOBS <- N_TMF_TILES * N_TMF_YEARS

#' Get sub-tiles for a given TMF tile
#' @param tmf_tile_id Integer TMF tile ID
#' @return Integer vector of sub_tile_ids
get_sub_tiles_for_tmf_tile <- function(target_tile_id) {
  SUB_TILE_INDEX[tmf_tile_id == target_tile_id]$sub_tile_id
}

# ==============================================================================
# SLURM UTILITIES
# ==============================================================================

#' Get SLURM array task ID
#' @param default Default value if not running under SLURM
#' @return Integer task ID
#' This cluster has MaxArraySize = 1001, so no array index above 1000 can be
#' submitted - yet stage 1 has 2,924 tile-year tasks. TASK_ID_OFFSET lets the
#' work go out as several arrays that each number themselves 1-1000, with the
#' offset added here to recover the true task id. Defaults to 0, so every other
#' stage is unaffected. See code/bash/1_extract_TMF.sh for the submission plan.
get_slurm_task_id <- function(default = 1L) {
  offset <- suppressWarnings(as.integer(Sys.getenv("TASK_ID_OFFSET", "0")))
  if (is.na(offset)) offset <- 0L

  task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
  if (task_id == "") {
    # Also check command line arguments (for local testing)
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) >= 1) {
      return(as.integer(args[1]) + offset)
    }
    return(as.integer(default) + offset)
  }
  return(as.integer(task_id) + offset)
}

#' Get SLURM job ID
#' @return String job ID or "local" if not on SLURM
get_slurm_job_id <- function() {
  job_id <- Sys.getenv("SLURM_JOB_ID", "local")
  return(job_id)
}

#' Convert job ID to (year, tile) combination for Stage 1
#' @param job_id Integer job ID (1 to N_EXTRACTION_JOBS)
#' @return List with year and tile_id
job_to_year_tile <- function(job_id) {
  # Jobs are organized as: tile_1_year_1, tile_1_year_2, ..., tile_1_year_N,
  #                        tile_2_year_1, tile_2_year_2, ...
  year_idx <- ((job_id - 1) %% N_TMF_YEARS) + 1
  tile_idx <- floor((job_id - 1) / N_TMF_YEARS) + 1

  return(list(
    year = TMF_YEARS[year_idx],
    tile_id = tile_idx
  ))
}

#' Convert (year, tile) combination to job ID
#' @param year Integer year
#' @param tile_id Integer tile ID
#' @return Integer job ID
year_tile_to_job <- function(year, tile_id) {
  year_idx <- which(TMF_YEARS == year)
  job_id <- (tile_id - 1) * N_TMF_YEARS + year_idx
  return(job_id)
}

#' Get sub-tile info by ID
#' @param sub_tile_id Integer sub-tile ID (1 to N_SUB_TILES)
#' @return List with sub-tile info
get_sub_tile_info <- function(sub_tile_id) {
  if (sub_tile_id < 1 || sub_tile_id > N_SUB_TILES) {
    stop(sprintf("Invalid sub_tile_id: %d (must be 1-%d)", sub_tile_id, N_SUB_TILES))
  }
  row <- SUB_TILE_INDEX[sub_tile_id]
  as.list(row)
}

# ==============================================================================
# FILE NAMING CONVENTIONS
# ==============================================================================

#' Get grid output filename for a sub-tile
#' @param sub_tile_id Integer sub-tile ID
#' @param ext File extension ("parquet" or "gpkg")
#' @return Full path to output file
get_grid_filename <- function(sub_tile_id, ext = "parquet") {
  file.path(grid_output_path,
            sprintf("grid_sub_%04d.%s", sub_tile_id, ext))
}

#' Get all grid files for a TMF tile (all sub-tiles within it)
#' @param tmf_tile_id Integer TMF tile ID
#' @param ext File extension
#' @return Character vector of file paths
get_grid_files_for_tmf_tile <- function(target_tile_id, ext = "parquet") {
  sub_ids <- SUB_TILE_INDEX[tmf_tile_id == target_tile_id]$sub_tile_id
  sapply(sub_ids, function(sid) get_grid_filename(sid, ext))
}

#' Get TMF extraction output filename
#' @param tile_id Integer tile ID
#' @param year Integer year
#' @return Full path to output file
get_tmf_filename <- function(tile_id, year) {
  file.path(tmf_output_path,
            sprintf("tmf_%03d_%d.csv.gz", tile_id, year))
}

#' Get consolidated tile filename
#' @param tile_id Integer tile ID
#' @return Full path to output file
get_consolidated_tile_filename <- function(tile_id) {
  file.path(consolidated_path,
            sprintf("tile_%03d.parquet", tile_id))
}

#' Get interior classification filename
#' @param tile_id Integer tile ID (or NULL for global)
#' @return Full path to output file
get_interior_filename <- function(tile_id = NULL) {
  if (is.null(tile_id)) {
    file.path(classification_path, "interior_all.parquet")
  } else {
    file.path(classification_path,
              sprintf("interior_tile_%03d.parquet", tile_id))
  }
}

#' Get frontier classification filename
#' @param tile_id Integer tile ID (or NULL for global)
#' @return Full path to output file
get_frontier_filename <- function(tile_id = NULL) {
  if (is.null(tile_id)) {
    file.path(classification_path, "frontier_all.parquet")
  } else {
    file.path(classification_path,
              sprintf("frontier_tile_%03d.parquet", tile_id))
  }
}

#' Get WDPA extraction filename
#' @param tile_id Integer tile ID (or NULL for consolidated)
#' @return Full path to output file
get_wdpa_filename <- function(tile_id = NULL) {
  if (is.null(tile_id)) {
    file.path(wdpa_output_path, "wdpa_consolidated.parquet")
  } else {
    file.path(wdpa_output_path,
              sprintf("wdpa_tile_%03d.parquet", tile_id))
  }
}

#' Get yields extraction filename
#' @param sub_tile_id Integer sub-tile ID
#' @return Full path to output file
get_yields_filename <- function(sub_tile_id) {
  file.path(yields_output_path,
            sprintf("yields_sub_%04d.parquet", sub_tile_id))
}

#' List all yield / crop-share raster files
#'
#' EarthStat ships one directory per crop, each holding several layers
#' (Production, HarvestedAreaHectares, DataQuality_*, ...). We want exactly
#' two of them: YieldPerHectare and HarvestedAreaFraction.
#'
#' Verified 2026-08-21 against the source tree: all 172 crop directories carry
#' all three of HarvestedAreaFraction / YieldPerHectare / HarvestedAreaHectares,
#' so this returns 344 rasters. An earlier comment here claimed the abaca..areca
#' directories were partially extracted and got skipped by the suffix match -
#' they are not, and they are not skipped. Do not use this function's return
#' length as evidence that some crops were dropped.
#'
#' @return Character vector of full paths to the wanted TIF files
list_yield_rasters <- function() {
  if (!dir.exists(yields_path)) {
    warning(sprintf("Yields directory not found: %s", yields_path))
    return(character(0))
  }

  tif_files <- list.files(
    yields_path,
    pattern = "_(YieldPerHectare|HarvestedAreaFraction)\\.tif$",
    full.names = TRUE, recursive = TRUE, ignore.case = TRUE
  )

  if (length(tif_files) == 0) {
    warning(sprintf("No EarthStat yield/crop-share rasters found in: %s",
                    yields_path))
  }

  return(sort(tif_files))
}

#' Get yield variable name from filename
#' @param filepath Path to yield raster
#' @return Clean variable name (without .tif extension)
get_yield_varname <- function(filepath) {
  stem <- tools::file_path_sans_ext(basename(filepath))

  # EarthStat: "<crop>_YieldPerHectare"      -> yield_<crop>
  #            "<crop>_HarvestedAreaFraction" -> cropshare_<crop>
  if (grepl("_YieldPerHectare$", stem, ignore.case = TRUE)) {
    crop <- sub("_YieldPerHectare$", "", stem, ignore.case = TRUE)
    return(paste0("yield_", tolower(crop)))
  }
  if (grepl("_HarvestedAreaFraction$", stem, ignore.case = TRUE)) {
    crop <- sub("_HarvestedAreaFraction$", "", stem, ignore.case = TRUE)
    return(paste0("cropshare_", tolower(crop)))
  }

  # Fallback for any other raster dropped into yields_path.
  varname <- tolower(gsub("[- ]+", "_", stem))
  if (!grepl("^yield", varname)) {
    varname <- paste0("yield_", varname)
  }
  return(varname)
}

#' List the EarthStat cropland / pasture fraction rasters
#'
#' The two 5-arcmin land-use layers that sit one level above the per-crop
#' directories. Returned as a NAMED vector (names become output column names)
#' so Stage 0c can push them through the same extraction loop as the crop
#' layers. A missing file is dropped with a warning rather than aborting the
#' sub-tile, matching list_yield_rasters()'s tolerance - these are covariates,
#' not a reason to lose 344 tasks.
#'
#' @return Named character vector of existing raster paths (possibly empty)
list_landuse_rasters <- function() {
  candidates <- c(cropland_frac = cropland_raster_path,
                  pasture_frac  = pasture_raster_path)

  ok <- file.exists(candidates)
  if (!all(ok)) {
    warning(sprintf("EarthStat land-use raster(s) not found: %s",
                    paste(candidates[!ok], collapse = ", ")))
  }

  return(candidates[ok])
}


# ==============================================================================
# GFED FILE UTILITIES
# ==============================================================================

#' List available GFED ecosystem files
#' @return Character vector of full paths to GFED NetCDF files
list_gfed_files <- function() {
  if (!dir.exists(gfed_data_path)) {
    stop(sprintf("GFED directory not found: %s", gfed_data_path))
  }

  nc_files <- list.files(gfed_data_path,
                         pattern = "^GFED5\\.1_ecosystem_\\d{4}\\.nc$",
                         full.names = TRUE)

  if (length(nc_files) == 0) {
    stop(sprintf("No GFED ecosystem files found in: %s", gfed_data_path))
  }

  return(sort(nc_files))
}

#' Get GFED years from available files
#' @return Integer vector of available years
get_gfed_years <- function() {
  files <- list_gfed_files()
  years <- as.integer(gsub(".*_(\\d{4})\\.nc$", "\\1", basename(files)))
  return(sort(years))
}

#' Get GFED file path for a specific year
#' @param year Integer year
#' @return Full path to GFED NetCDF file
get_gfed_file <- function(year) {
  file.path(gfed_data_path, sprintf("GFED5.1_ecosystem_%d.nc", year))
}

#' Get GFED sub-tile output filename
#' @param sub_tile_id Integer sub-tile ID
#' @return Full path to output file
get_gfed_subtile_filename <- function(sub_tile_id) {
  file.path(gfed_output_path,
            sprintf("gfed_sub_%04d.parquet", sub_tile_id))
}

#' Get consolidated GFED filename by tile
#' @param tile_id Integer tile ID
#' @return Full path to output file
get_gfed_consolidated_filename <- function(tile_id) {
  file.path(gfed_consolidated_path,
            sprintf("gfed_tile_%03d.parquet", tile_id))
}

#' Convert GFED NetCDF burned area to terra raster for a specific month
#' @param nc_file Path to GFED NetCDF file
#' @param month Integer month (1-12)
#' @return terra SpatRaster with burned fraction (burned_area / grid_area)
gfed_to_raster <- function(nc_file, month) {
  # terra::rast() reads NetCDF directly with correct orientation
  # GFED has 12 layers (one per month) for burned_area
  ba_rast <- terra::rast(nc_file, subds = "burned_area", lyrs = month)
  ga_rast <- terra::rast(nc_file, subds = "grid_area")

  # Compute burned fraction
  burned_frac <- ba_rast / ga_rast

  # Replace NA/Inf with 0
  burned_frac[is.na(burned_frac) | is.infinite(burned_frac)] <- 0

  return(burned_frac)
}

#' Get covariates filename
#' @param tile_id Integer tile ID
#' @return Full path to output file
get_covariates_filename <- function(tile_id) {
  file.path(covariates_output_path,
            sprintf("covariates_tile_%03d.parquet", tile_id))
}

#' Load all grid cells for a TMF tile
#' @param target_tile_id Integer TMF tile ID
#' @param with_geometry If TRUE, return sf; if FALSE, return data.table
#' @return sf object or data.table with all grid cells for the tile
load_grid_for_tile <- function(target_tile_id, with_geometry = TRUE) {
  sub_ids <- get_sub_tiles_for_tmf_tile(target_tile_id)

  if (with_geometry) {
    grid_files <- sapply(sub_ids, function(sid) get_grid_filename(sid, "gpkg"))
    grid_files <- grid_files[file.exists(grid_files)]

    if (length(grid_files) == 0) {
      stop(sprintf("No grid files found for tile %d", target_tile_id))
    }

    grid_list <- lapply(grid_files, function(f) st_read(f, quiet = TRUE))
    grid_sf <- do.call(rbind, grid_list)
    return(grid_sf)
  } else {
    grid_files <- sapply(sub_ids, function(sid) get_grid_filename(sid, "parquet"))
    grid_files <- grid_files[file.exists(grid_files)]

    if (length(grid_files) == 0) {
      stop(sprintf("No grid files found for tile %d", target_tile_id))
    }

    grid_list <- lapply(grid_files, function(f) {
      dt <- arrow::read_parquet(f)
      setDT(dt)
      return(dt)
    })
    return(rbindlist(grid_list, fill = TRUE))
  }
}

# ==============================================================================
# RESUMABILITY UTILITIES
# ==============================================================================

#' Check if output file exists and is complete
#' @param filepath Path to output file
#' @param min_size Minimum file size in bytes (default 100)
#' @return TRUE if file exists and is valid
output_exists <- function(filepath, min_size = 100) {
  if (!file.exists(filepath)) return(FALSE)
  file_info <- file.info(filepath)
  return(!is.na(file_info$size) && file_info$size >= min_size)
}

#' Write data atomically (write to temp, then rename)
#' @param data Data to write
#' @param filepath Target output path
#' @param write_fn Function to use for writing (default: write_parquet)
write_atomic <- function(data, filepath, write_fn = arrow::write_parquet) {
  temp_file <- paste0(filepath, ".tmp.", get_slurm_job_id())
  write_fn(data, temp_file)
  file.rename(temp_file, filepath)
  message(sprintf("Wrote: %s", filepath))
}

#' Skip job if output already exists
#' @param filepath Path to check
#' @param job_desc Description for logging
skip_if_exists <- function(filepath, job_desc = "this job") {
  if (output_exists(filepath)) {
    message(sprintf("Output already exists for %s: %s", job_desc, filepath))
    message("Skipping.")
    quit(save = "no", status = 0)
  }
}

# ==============================================================================
# LOGGING UTILITIES
# ==============================================================================

#' Print timestamped message
#' @param msg Message to print
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, msg))
}

#' Log job start
#' @param script_name Name of the current script
#' @param task_id SLURM task ID
log_job_start <- function(script_name, task_id = NULL) {
  log_message("========================================")
  log_message(sprintf("Script: %s", script_name))
  log_message(sprintf("Job ID: %s", get_slurm_job_id()))
  if (!is.null(task_id)) {
    log_message(sprintf("Task ID: %d", task_id))
  }
  log_message(sprintf("Host: %s", Sys.info()["nodename"]))
  log_message("========================================")
}

#' Log job completion
#' @param start_time Start time from Sys.time()
log_job_end <- function(start_time) {
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  log_message("========================================")
  log_message(sprintf("Completed in %.1f minutes", as.numeric(elapsed)))
  log_message("========================================")
}

# ==============================================================================
# GEOMETRY UTILITIES
# ==============================================================================

#' Make geometries valid with fallback strategies
#' @param sf_obj sf object with potentially invalid geometries
#' @return sf object with valid geometries
make_valid_safe <- function(sf_obj) {
  # First pass: standard st_make_valid
  sf_obj <- st_make_valid(sf_obj)

  # Check for remaining invalid geometries
  invalid_idx <- which(!st_is_valid(sf_obj))

  if (length(invalid_idx) > 0) {
    log_message(sprintf("Found %d invalid geometries, attempting buffer(0) fix",
                        length(invalid_idx)))

    # Try buffer(0) trick
    sf_use_s2(FALSE)
    sf_obj[invalid_idx, ] <- st_buffer(sf_obj[invalid_idx, ], 0)
    sf_use_s2(TRUE)

    # Revalidate
    sf_obj <- st_make_valid(sf_obj)

    # Final check
    still_invalid <- sum(!st_is_valid(sf_obj))
    if (still_invalid > 0) {
      warning(sprintf("%d geometries still invalid after repair attempts",
                      still_invalid))
    }
  }

  return(sf_obj)
}

# ==============================================================================
# WDPA READING (handles corrupt geometries)
# ==============================================================================

# Path to cleaned WDPA (created once, reused by all jobs)
wdpa_clean_path <- file.path(build_data_path, "wdpa_cleaned.gpkg")

#' Create cleaned WDPA file using gdal_utils (ogr2ogr in R)
#' Drops corrupt geometries and reports statistics
#' @return Path to cleaned file
create_clean_wdpa <- function() {
  if (file.exists(wdpa_clean_path)) {
    log_message(sprintf("Cleaned WDPA exists: %s", wdpa_clean_path))
    return(wdpa_clean_path)
  }

  log_message("Creating cleaned WDPA file...")
  log_message(sprintf("Source: %s", wdpa_path))

  # Use gdal_utils (R binding for ogr2ogr) with -skipfailures
  # Also reproject to WGS84 for consistency with grid data
  # Note: -makevalid requires GDAL 3.1+, so we handle validation in R
  sf::gdal_utils(
    util = "vectortranslate",
    source = wdpa_path,
    destination = wdpa_clean_path,
    options = c(
      "-f", "GPKG",
      "-skipfailures",
      "-nlt", "PROMOTE_TO_MULTI",
      "-t_srs", "EPSG:4326"  # Reproject to WGS84
    )
  )

  # Validate geometries in R (since -makevalid may not be available)
  log_message("Validating geometries in R...")
  wdpa_temp <- st_read(wdpa_clean_path, quiet = TRUE)

  # Check and fix invalid geometries
  invalid_idx <- !st_is_valid(wdpa_temp)
  n_invalid <- sum(invalid_idx, na.rm = TRUE)

  if (n_invalid > 0) {
    log_message(sprintf("Found %d invalid geometries, attempting repair...", n_invalid))
    wdpa_temp <- st_make_valid(wdpa_temp)

    # Remove any empty geometries that result from repair
    empty_idx <- st_is_empty(wdpa_temp)
    if (any(empty_idx)) {
      log_message(sprintf("Removing %d empty geometries", sum(empty_idx)))
      wdpa_temp <- wdpa_temp[!empty_idx, ]
    }

    # Rewrite the cleaned file
    st_write(wdpa_temp, wdpa_clean_path, delete_dsn = TRUE, quiet = TRUE)
    log_message("Repaired geometries saved")
  }

  if (!file.exists(wdpa_clean_path)) {
    stop("Failed to create cleaned WDPA file")
  }

  # Check file size
  file_size <- file.info(wdpa_clean_path)$size
  log_message(sprintf("Cleaned file size: %.1f MB", file_size / 1024 / 1024))

  if (file_size < 1000) {
    # File is suspiciously small - delete and fail
    unlink(wdpa_clean_path)
    stop("Cleaned WDPA file is empty or too small. gdal_utils may have failed.")
  }

  # Report statistics
  log_message("Comparing original vs cleaned...")

  # Get original feature count from st_layers (works for GDB, GPKG, shapefiles)
  layers_info <- st_layers(wdpa_path)
  log_message(sprintf("Source layer: %s", layers_info$name[1]))
  n_orig <- layers_info$features[1]

  # Read cleaned file to get actual count
  wdpa_clean <- st_read(wdpa_clean_path, quiet = TRUE)
  n_clean <- nrow(wdpa_clean)

  if (!is.na(n_orig) && n_orig > 0) {
    n_dropped <- n_orig - n_clean
    pct_dropped <- 100 * n_dropped / n_orig

    log_message(sprintf("Original features: %d", n_orig))
    log_message(sprintf("Cleaned features: %d", n_clean))
    log_message(sprintf("Dropped: %d (%.2f%%)", n_dropped, pct_dropped))

    # Check for area column
    area_cols <- grep("area|AREA|Area|REP_AREA|GIS_AREA", names(wdpa_clean), value = TRUE)
    if (length(area_cols) > 0) {
      area_col <- area_cols[1]
      total_area_clean <- sum(as.numeric(wdpa_clean[[area_col]]), na.rm = TRUE)
      log_message(sprintf("Area column: %s", area_col))
      log_message(sprintf("Total area in cleaned data: %.0f sq km", total_area_clean))
    }
  } else {
    log_message(sprintf("Cleaned features: %d", n_clean))
  }

  log_message(sprintf("Cleaned WDPA saved: %s", wdpa_clean_path))
  return(wdpa_clean_path)
}

#' Read WDPA data with optional spatial filter
#' @param bbox Optional bounding box (xmin, ymin, xmax, ymax)
#' @return sf object with WDPA polygons
read_wdpa <- function(bbox = NULL) {
  # Ensure cleaned file exists
  wdpa_file <- create_clean_wdpa()

  log_message(sprintf("Reading WDPA from: %s", wdpa_file))

  # Read with spatial filter if provided
  if (!is.null(bbox)) {
    bbox_wkt <- sprintf(
      "POLYGON((%f %f, %f %f, %f %f, %f %f, %f %f))",
      bbox[1], bbox[2],  # xmin, ymin
      bbox[3], bbox[2],  # xmax, ymin
      bbox[3], bbox[4],  # xmax, ymax
      bbox[1], bbox[4],  # xmin, ymax
      bbox[1], bbox[2]   # close
    )
    wdpa_sf <- st_read(wdpa_file, quiet = TRUE, wkt_filter = bbox_wkt)
  } else {
    wdpa_sf <- st_read(wdpa_file, quiet = TRUE)
  }

  log_message(sprintf("Read %d features", nrow(wdpa_sf)))
  return(wdpa_sf)
}

#' Get tile extent as sf polygon
#' @param tile_id Integer tile ID
#' @param crs Target CRS (default WGS84)
#' @return sf polygon of tile extent
get_tile_extent_sf <- function(target_tile_id, crs = WGS84_CRS) {
  tile <- TMF_TILE_INDEX[tile_id == target_tile_id]

  bbox <- st_bbox(c(
    xmin = tile$lon_west,
    ymin = tile$lat_south,
    xmax = tile$lon_east,
    ymax = tile$lat_north
  ), crs = st_crs(WGS84_CRS))

  extent_sf <- st_as_sfc(bbox)

  if (!is.null(crs) && crs != WGS84_CRS) {
    extent_sf <- st_transform(extent_sf, crs)
  }

  return(extent_sf)
}

# ==============================================================================
# TMF FILE UTILITIES
# ==============================================================================

#' Get TMF raster filename for a tile and year
#' @param target_tile_id Integer tile ID
#' @param year Integer year
#' @return Full path to TMF raster file
get_tmf_raster_path <- function(target_tile_id, year) {
  tile <- TMF_TILE_INDEX[tile_id == target_tile_id]

  # Folder: JRC_TMF_AnnualChange_v1_AFR_ID6_S20_E40
  # File:   JRC_TMF_AnnualChange_v1_1990_AFR_ID6_S20_E40.tif
  # Insert year after "v1_"
  filename <- sub("_v1_", sprintf("_v1_%d_", year), tile$folder)
  filename <- paste0(filename, ".tif")

  file.path(tmf_path, tile$folder, filename)
}

#' Check if TMF raster exists for a tile and year
#' @param tile_id Integer tile ID
#' @param year Integer year
#' @return TRUE if raster file exists
tmf_raster_exists <- function(tile_id, year) {
  filepath <- get_tmf_raster_path(tile_id, year)
  return(file.exists(filepath))
}

# ==============================================================================
# MEMORY MANAGEMENT
# ==============================================================================

#' Force garbage collection and report memory
gc_verbose <- function() {
  gc_result <- gc(verbose = FALSE, reset = TRUE)
  used_mb <- sum(gc_result[, 2])
  log_message(sprintf("Memory used: %.1f MB", used_mb))
  invisible(gc_result)
}

# ==============================================================================
# STARTUP MESSAGE
# ==============================================================================

log_message(sprintf("BUILD_workspace.R loaded from: %s", project_root))
log_message(sprintf("TMF years: %d-%d (%d years)",
                    min(TMF_YEARS), max(TMF_YEARS), N_TMF_YEARS))
log_message(sprintf("TMF tiles (10x10): %d", N_TMF_TILES))
log_message(sprintf("Sub-tiles (%dx%d): %d", SUB_TILE_SIZE, SUB_TILE_SIZE, N_SUB_TILES))
log_message(sprintf("Total extraction jobs: %d", N_EXTRACTION_JOBS))
