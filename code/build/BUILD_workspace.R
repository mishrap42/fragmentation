# ==============================================================================
# TMF 1KM PIPELINE - SHARED CONFIGURATION
# ==============================================================================
# This file contains all shared configuration, paths, and utility functions
# for the Tropical Moist Forest data pipeline.
#
# Source this file at the beginning of every pipeline script:
#   source("Code/BUILD_workspace.R")
# ==============================================================================

library(here)
here::i_am("Code/BUILD_workspace.R")

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
})

# ==============================================================================
# DATA PATHS - CONFIGURE THESE FOR YOUR SYSTEM
# ==============================================================================

# Project root (auto-detected via here package)
project_root <- here::i_am('code/build/BUILD_workspace.R')

# Input data paths - UPDATE THESE TO YOUR ACTUAL DATA LOCATIONS
tmf_path <- file.path(here(), 'Data', 'raw', 'TMF', 'AnnualChange')        # TMF GeoTIFF tiles
gadm_path <- file.path(here(), 'Data', 'raw', 'GADM', 'gadm_410-levels.gpkg')  # GADM country boundaries
wdpa_path <- file.path(here(), 'Data', 'raw', 'TMF', 'protected sites', 'protectedsites.shp')    # WDPA protected areas

# Output paths (relative to project root)
raw_data_path <- file.path(project_root, "Data", "raw")
build_data_path <- file.path(project_root, "Data", "build")

# Specific output directories
grid_output_path <- file.path(build_data_path, "grids")
tmf_output_path <- file.path(build_data_path, "tmf")
consolidated_path <- file.path(build_data_path, "tmf_consolidated")
classification_path <- file.path(build_data_path, "classifications")
wdpa_output_path <- file.path(build_data_path, "wdpa")
final_output_path <- file.path(build_data_path, "final")

# Create directories if they don't exist
sapply(c(grid_output_path, tmf_output_path, consolidated_path,
         classification_path, wdpa_output_path, final_output_path),
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

# Grid cell size in meters (1 km x 1 km)
GRID_CELLSIZE_M <- 1000

# Buffer distance for forest frontier calculation (100 km)
BUFFER_DISTANCE_M <- 100000

# Minimum cell area fraction (filter slivers < 1% of target)
MIN_CELL_AREA_FRAC <- 0.01

# ==============================================================================
# TMF DATA PARAMETERS
# ==============================================================================

# Years of TMF data availability
TMF_YEARS <- 1990:2024
N_TMF_YEARS <- length(TMF_YEARS)

# Number of TMF tiles (10x10 degree tiles covering tropical belt)
N_TMF_TILES <- 86

# TMF land cover classification legend
TMF_LEGEND <- data.table(
  value = c(1, 2, 3, 4, 5, 6),
  tmf_class = c("Undisturbed_TMF", "Degraded_TMF", "Deforested",
                "Regrowth", "Water", "Other_landcover")
)

# Total number of extraction jobs (tiles x years)
N_EXTRACTION_JOBS <- N_TMF_TILES * N_TMF_YEARS

# ==============================================================================
# TMF TILE INDEX
# ==============================================================================

# TMF tile definitions (lat/lon of upper-left corner for each tile)
# This covers the pan-tropical belt from ~23.5N to ~23.5S
# Format: tile_id, lat_north, lon_west

create_tmf_tile_index <- function() {
  # Latitude bands: 20N to 30S in 10-degree increments
  lats <- c(20, 10, 0, -10, -20, -30)

  # Longitude bands: -180 to 180 in 10-degree increments
  lons <- seq(-180, 170, by = 10)

  # Create all combinations
  tiles <- expand.grid(lat_north = lats, lon_west = lons)
  tiles <- as.data.table(tiles)
  tiles[, tile_id := .I]
  tiles[, lat_south := lat_north - 10]
  tiles[, lon_east := lon_west + 10]

  # Filter to only tropical moist forest regions
  # (This is a simplified filter - actual TMF tiles cover specific regions)
  # In practice, you should use the actual TMF tile availability
  tiles <- tiles[lat_north >= -30 & lat_south <= 30]
  tiles[, tile_id := .I]  # Re-index

  return(tiles)
}

TMF_TILE_INDEX <- create_tmf_tile_index()

# ==============================================================================
# SLURM UTILITIES
# ==============================================================================

#' Get SLURM array task ID
#' @param default Default value if not running under SLURM
#' @return Integer task ID
get_slurm_task_id <- function(default = 1L) {
  task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
  if (task_id == "") {
    # Also check command line arguments (for local testing)
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) >= 1) {
      return(as.integer(args[1]))
    }
    return(as.integer(default))
  }
  return(as.integer(task_id))
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

# ==============================================================================
# FILE NAMING CONVENTIONS
# ==============================================================================

#' Get grid output filename for a tile
#' @param tile_id Integer tile ID
#' @param ext File extension ("parquet" or "gpkg")
#' @return Full path to output file
get_grid_filename <- function(tile_id, ext = "parquet") {
  file.path(grid_output_path,
            sprintf("grid_tile_%03d.%s", tile_id, ext))
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

#' Get tile extent as sf polygon
#' @param tile_id Integer tile ID
#' @param crs Target CRS (default WGS84)
#' @return sf polygon of tile extent
get_tile_extent_sf <- function(tile_id, crs = WGS84_CRS) {
  tile <- TMF_TILE_INDEX[tile_id == tile_id]

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
#' @param tile_id Integer tile ID
#' @param year Integer year
#' @return Full path to TMF raster file
get_tmf_raster_path <- function(tile_id, year) {
  tile <- TMF_TILE_INDEX[tile_id == tile_id]

  # Construct filename based on TMF naming convention
  # Format: AnnualChange_YYYY_latN/SXX_lonE/WXXX.tif
  lat_str <- ifelse(tile$lat_north >= 0,
                    sprintf("N%02d", abs(tile$lat_north)),
                    sprintf("S%02d", abs(tile$lat_north)))
  lon_str <- ifelse(tile$lon_west >= 0,
                    sprintf("E%03d", abs(tile$lon_west)),
                    sprintf("W%03d", abs(tile$lon_west)))

  filename <- sprintf("AnnualChange_%d_%s_%s.tif", year, lat_str, lon_str)

  return(file.path(tmf_path, filename))
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
log_message(sprintf("TMF tiles: %d", N_TMF_TILES))
log_message(sprintf("Total extraction jobs: %d", N_EXTRACTION_JOBS))
