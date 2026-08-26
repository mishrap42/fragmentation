# ==============================================================================
# STAGE 0d: EXTRACT GAEZ v4 ATTAINABLE YIELD RASTERS
# ==============================================================================
# Coverage-weighted mean of each GAEZ attainable-yield raster onto the grid.
# Structurally identical to Stage 0c; separate because the source, the units and
# the interpretation are all different, and because 0c has already run.
#
# Input:  SLURM_ARRAY_TASK_ID (1 to N_SUB_TILES)
# Output: Data/build/gaez/gaez_sub_{sub_tile_id}.parquet
#           gaez_<crop>200a   agro-climatic attainable yield,   kg/ha DM
#           gaez_<crop>200b   agro-ecological attainable yield, kg/ha DM
#
# WHY THIS EXISTS
# ---------------
# EarthStat yield is measured on the crop's HARVESTED area and is zero wherever
# the crop is not grown, so on undisturbed forest it is identically zero. A
# profitability measure built on it is therefore close to an indicator for "this
# land has already been cleared", which is why beta_p in the siting logit fell
# from -1.61 (t = 4.5) to -0.27 (t = 1.75) the moment 1990 forest cover or
# above-ground carbon entered. GAEZ is an agronomic simulation defined on all
# land irrespective of what is planted, so it carries land-quality variation
# that survives conditioning on forest state.
#
# UNITS: kg/ha dry matter. Divide by 1000 before meeting USD/tonne prices.
#
# PREREQUISITE: Stage 0 (grid files)
# ==============================================================================

here::i_am('code/build/0d_extract_gaez.R')
source("code/build/BUILD_workspace.R")

start_time <- Sys.time()

task_id <- get_slurm_task_id()
sub_tile_id <- task_id

log_job_start("0d_extract_gaez.R", task_id)
log_message(sprintf("Processing sub-tile %d", sub_tile_id))

output_file <- get_gaez_filename(sub_tile_id)
skip_if_exists(output_file, sprintf("sub-tile %d", sub_tile_id))

# ==============================================================================
# PREREQUISITES
# ==============================================================================

grid_gpkg <- get_grid_filename(sub_tile_id, "gpkg")
if (!file.exists(grid_gpkg)) {
  stop(sprintf("Grid file not found: %s\nRun Stage 0 first.", grid_gpkg))
}

gaez_rasters <- list_gaez_rasters()

if (length(gaez_rasters) == 0) {
  log_message("No GAEZ rasters found. Creating empty output.")
  write_atomic(data.table(grid_id = character(0)), output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

log_message(sprintf("Found %d GAEZ rasters to extract", length(gaez_rasters)))

# ==============================================================================
# GRID
# ==============================================================================

log_message("Loading grid geometries...")
grid_sf <- st_read(grid_gpkg, quiet = TRUE)
log_message(sprintf("Loaded %d grid cells", nrow(grid_sf)))

if (nrow(grid_sf) == 0) {
  log_message("Empty grid. Creating empty output.")
  write_atomic(data.table(grid_id = character(0)), output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# GAEZ is WGS84 at 5 arcmin, the same grid as EarthStat, so a single transform
# up front is enough and the per-raster loop needs no reprojection.
grid_sf <- st_transform(grid_sf, WGS84_CRS)

# ==============================================================================
# EXTRACT
# ==============================================================================

gaez_data <- data.table(grid_id = grid_sf$grid_id)

#' As 0c: coverage-weighted mean, NA-skipping, NaN normalised to NA, and a
#' raster that cannot be read costs its own column rather than the sub-tile.
extract_mean_layer <- function(raster_path, varname) {
  tryCatch({
    r <- terra::rast(raster_path)
    values <- exact_extract(r, grid_sf, "mean")
    values[is.nan(values)] <- NA_real_
    gaez_data[, (varname) := values]
  }, error = function(e) {
    log_message(sprintf("  WARNING: Failed to extract %s: %s", varname, e$message))
    gaez_data[, (varname) := NA_real_]
  })
  gc()
  invisible(NULL)
}

log_message("Extracting GAEZ attainable yields...")

n_rasters <- length(gaez_rasters)
for (i in seq_along(gaez_rasters)) {
  varname <- get_gaez_varname(gaez_rasters[i])
  if (i %% 25 == 1 || i == n_rasters) {
    log_message(sprintf("  Extracting raster %d/%d: %s", i, n_rasters, varname))
  }
  extract_mean_layer(gaez_rasters[i], varname)
}

# ==============================================================================
# SUMMARY
# ==============================================================================

gaez_cols <- setdiff(names(gaez_data), "grid_id")
log_message("Summary statistics:")
log_message(sprintf("  Grid cells: %d", nrow(gaez_data)))
log_message(sprintf("  GAEZ variables: %d (a: %d, b: %d)",
                    length(gaez_cols),
                    sum(grepl("200a$", gaez_cols)),
                    sum(grepl("200b$", gaez_cols))))

for (col in head(gaez_cols, 3)) {
  v <- gaez_data[[col]]
  if (sum(!is.na(v)) > 0) {
    log_message(sprintf("    %s: mean=%.1f kg/ha, non-NA=%d/%d",
                        col, mean(v, na.rm = TRUE), sum(!is.na(v)), length(v)))
  }
}

# ==============================================================================
# WRITE
# ==============================================================================

log_message(sprintf("Writing output: %s", output_file))
write_atomic(gaez_data, output_file)

if (file.exists(output_file)) {
  log_message(sprintf("Output file size: %.1f KB",
                      file.info(output_file)$size / 1024))
}

gc_verbose()
log_job_end(start_time)
