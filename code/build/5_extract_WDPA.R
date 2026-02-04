# ==============================================================================
# STAGE 5: EXTRACT WDPA PROTECTED AREA DATA
# ==============================================================================
# This script extracts WDPA (World Database on Protected Areas) metadata
# for grid cells classified as forest interior or frontier.
#
# Input: SLURM_ARRAY_TASK_ID (1-86 for tiles, 87 for consolidation)
# Output: Data/build/wdpa/wdpa_tile_{tile_id}.parquet
#         Data/build/wdpa/wdpa_consolidated.parquet (task 87)
# ==============================================================================

# Load configuration
source("Code/BUILD_workspace.R")

# Record start time
start_time <- Sys.time()

# Get task ID
task_id <- get_slurm_task_id()

log_job_start("5_extract_WDPA.R", task_id)

# ==============================================================================
# HANDLE CONSOLIDATION JOB (TASK 87)
# ==============================================================================

if (task_id == N_TMF_TILES + 1) {
  log_message("Running consolidation job...")

  output_file <- get_wdpa_filename(NULL)  # Consolidated file
  skip_if_exists(output_file, "WDPA consolidation")

  # Find all tile files
  tile_files <- sapply(1:N_TMF_TILES, get_wdpa_filename)
  files_exist <- file.exists(tile_files)

  log_message(sprintf("Found %d/%d tile files",
                      sum(files_exist), length(tile_files)))

  if (sum(files_exist) == 0) {
    stop("No WDPA tile files found. Run extraction jobs first.")
  }

  # Read and combine
  all_data <- lapply(tile_files[files_exist], function(f) {
    dt <- read_parquet(f)
    setDT(dt)
    return(dt)
  })

  combined <- rbindlist(all_data, fill = TRUE)
  log_message(sprintf("Combined: %s rows", format(nrow(combined), big.mark = ",")))

  # Remove duplicates (cells at tile boundaries)
  combined <- unique(combined, by = "grid_id")
  log_message(sprintf("After deduplication: %s rows",
                      format(nrow(combined), big.mark = ",")))

  write_atomic(combined, output_file)

  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# ==============================================================================
# REGULAR TILE EXTRACTION
# ==============================================================================

tile_id <- task_id
log_message(sprintf("Processing tile: %d", tile_id))

# Validate
if (tile_id < 1 || tile_id > N_TMF_TILES) {
  stop(sprintf("Invalid tile_id: %d", tile_id))
}

# Check if output already exists
output_file <- get_wdpa_filename(tile_id)
skip_if_exists(output_file, sprintf("tile %d", tile_id))

# ==============================================================================
# LOAD GRID CELLS FOR THIS TILE
# ==============================================================================

log_message("Loading grid cells...")

grid_file <- get_grid_filename(tile_id, "gpkg")

if (!file.exists(grid_file)) {
  stop(sprintf("Grid file not found: %s", grid_file))
}

grid_sf <- st_read(grid_file, quiet = TRUE)
log_message(sprintf("Loaded %d grid cells", nrow(grid_sf)))

if (nrow(grid_sf) == 0) {
  log_message("Empty grid. Creating empty output.")
  empty_output <- data.table(
    grid_id = character(0),
    tile_id = integer(0),
    is_protected = logical(0),
    protected_frac = numeric(0),
    wdpa_pid = integer(0),
    iucn_cat = character(0),
    desig_year = integer(0),
    gov_type = character(0)
  )
  write_atomic(empty_output, output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# ==============================================================================
# LOAD INTERIOR AND FRONTIER CLASSIFICATIONS
# ==============================================================================

log_message("Loading classifications...")

interior_file <- get_interior_filename(tile_id)
frontier_file <- get_frontier_filename(tile_id)

if (!file.exists(interior_file) || !file.exists(frontier_file)) {
  stop("Classification files not found. Run Stages 3 and 4 first.")
}

interior <- read_parquet(interior_file)
setDT(interior)

frontier <- read_parquet(frontier_file)
setDT(frontier)

# Identify cells of interest (interior OR frontier)
cells_of_interest <- unique(c(
  interior[is_interior == TRUE]$grid_id,
  frontier[is_frontier == TRUE]$grid_id
))

log_message(sprintf("Cells of interest (interior or frontier): %s",
                    format(length(cells_of_interest), big.mark = ",")))

if (length(cells_of_interest) == 0) {
  log_message("No interior or frontier cells. Creating output with all cells marked as non-protected.")

  output_data <- data.table(
    grid_id = grid_sf$grid_id,
    tile_id = tile_id,
    is_protected = FALSE,
    protected_frac = 0,
    wdpa_pid = NA_integer_,
    iucn_cat = NA_character_,
    desig_year = NA_integer_,
    gov_type = NA_character_
  )

  write_atomic(output_data, output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# Filter grid to cells of interest
grid_sf <- grid_sf[grid_sf$grid_id %in% cells_of_interest, ]
log_message(sprintf("Filtered to %d cells of interest", nrow(grid_sf)))

# ==============================================================================
# LOAD WDPA DATA
# ==============================================================================

log_message("Loading WDPA data...")

# Get tile bounding box with buffer
tile_bbox <- st_bbox(grid_sf)

# Add 10km buffer for edge cases
buffer_deg <- 0.1  # ~10km
tile_bbox_buffered <- c(
  xmin = tile_bbox["xmin"] - buffer_deg,
  ymin = tile_bbox["ymin"] - buffer_deg,
  xmax = tile_bbox["xmax"] + buffer_deg,
  ymax = tile_bbox["ymax"] + buffer_deg
)

# Read WDPA with bounding box filter
log_message(sprintf("Reading WDPA for bbox: [%.2f, %.2f] to [%.2f, %.2f]",
                    tile_bbox_buffered[1], tile_bbox_buffered[2],
                    tile_bbox_buffered[3], tile_bbox_buffered[4]))

# Note: This assumes WDPA file supports spatial filtering
# If not, we read all and filter
wdpa <- tryCatch({
  st_read(wdpa_path, quiet = TRUE,
          wkt_filter = st_as_text(st_as_sfc(st_bbox(tile_bbox_buffered))))
}, error = function(e) {
  log_message("Spatial filter failed, reading full WDPA and filtering...")
  wdpa_full <- st_read(wdpa_path, quiet = TRUE)
  bbox_poly <- st_as_sfc(st_bbox(tile_bbox_buffered), crs = st_crs(wdpa_full))
  wdpa_full[st_intersects(wdpa_full, bbox_poly, sparse = FALSE)[, 1], ]
})

log_message(sprintf("Loaded %d protected areas", nrow(wdpa)))

if (nrow(wdpa) == 0) {
  log_message("No protected areas in this tile.")

  output_data <- data.table(
    grid_id = grid_sf$grid_id,
    tile_id = tile_id,
    is_protected = FALSE,
    protected_frac = 0,
    wdpa_pid = NA_integer_,
    iucn_cat = NA_character_,
    desig_year = NA_integer_,
    gov_type = NA_character_
  )

  write_atomic(output_data, output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# ==============================================================================
# CLEAN WDPA GEOMETRIES
# ==============================================================================

log_message("Cleaning WDPA geometries...")

# Make valid
wdpa <- make_valid_safe(wdpa)

# Remove empty geometries
wdpa <- wdpa[!st_is_empty(wdpa), ]

log_message(sprintf("After cleaning: %d protected areas", nrow(wdpa)))

# ==============================================================================
# EXTRACT WDPA METADATA FOR EACH CELL
# ==============================================================================

log_message("Extracting WDPA metadata for each cell...")

sf_use_s2(FALSE)

# Calculate intersection between grid cells and protected areas
# This gives us the area of each PA in each cell
intersections <- st_intersection(grid_sf, wdpa)

sf_use_s2(TRUE)

log_message(sprintf("Found %d cell-PA intersections", nrow(intersections)))

if (nrow(intersections) == 0) {
  log_message("No intersections with protected areas.")

  output_data <- data.table(
    grid_id = grid_sf$grid_id,
    tile_id = tile_id,
    is_protected = FALSE,
    protected_frac = 0,
    wdpa_pid = NA_integer_,
    iucn_cat = NA_character_,
    desig_year = NA_integer_,
    gov_type = NA_character_
  )

  write_atomic(output_data, output_file)
  log_job_end(start_time)
  quit(save = "no", status = 0)
}

# Calculate area of each intersection
intersections$intersection_area <- as.numeric(st_area(intersections))

# Convert to data.table for processing
int_dt <- as.data.table(st_drop_geometry(intersections))

# ==============================================================================
# CALCULATE PROTECTED FRACTION AND SELECT DOMINANT PA
# ==============================================================================

log_message("Calculating protected fractions...")

# Calculate total cell area (from original grid)
grid_areas <- data.table(
  grid_id = grid_sf$grid_id,
  cell_area = as.numeric(st_area(grid_sf))
)

# Sum intersection areas per cell
cell_protection <- int_dt[, .(
  total_protected_area = sum(intersection_area, na.rm = TRUE)
), by = grid_id]

# Merge with cell areas
cell_protection <- merge(cell_protection, grid_areas, by = "grid_id")
cell_protection[, protected_frac := pmin(total_protected_area / cell_area, 1)]

# For each cell, get the dominant PA (largest intersection)
dominant_pa <- int_dt[, .SD[which.max(intersection_area)], by = grid_id]

# Select relevant PA attributes
# Standard WDPA columns (may vary by version)
pa_cols <- intersect(names(dominant_pa),
                     c("WDPAID", "WDPA_PID", "NAME", "IUCN_CAT",
                       "STATUS_YR", "GOV_TYPE", "DESIG_TYPE"))

if (length(pa_cols) > 0) {
  dominant_pa <- dominant_pa[, c("grid_id", pa_cols), with = FALSE]
} else {
  log_message("Warning: No standard WDPA columns found")
  dominant_pa <- dominant_pa[, .(grid_id)]
}

# Standardize column names
if ("WDPAID" %in% names(dominant_pa)) {
  setnames(dominant_pa, "WDPAID", "wdpa_pid")
} else if ("WDPA_PID" %in% names(dominant_pa)) {
  setnames(dominant_pa, "WDPA_PID", "wdpa_pid")
}

if ("IUCN_CAT" %in% names(dominant_pa)) {
  setnames(dominant_pa, "IUCN_CAT", "iucn_cat")
}

if ("STATUS_YR" %in% names(dominant_pa)) {
  setnames(dominant_pa, "STATUS_YR", "desig_year")
}

if ("GOV_TYPE" %in% names(dominant_pa)) {
  setnames(dominant_pa, "GOV_TYPE", "gov_type")
}

# Merge protection info with dominant PA
protected_cells <- merge(cell_protection, dominant_pa, by = "grid_id", all.x = TRUE)
protected_cells[, is_protected := protected_frac > 0]

# ==============================================================================
# ADD NON-PROTECTED CELLS
# ==============================================================================

log_message("Adding non-protected cells...")

all_grid_ids <- grid_sf$grid_id
protected_ids <- protected_cells$grid_id
non_protected_ids <- setdiff(all_grid_ids, protected_ids)

if (length(non_protected_ids) > 0) {
  non_protected <- data.table(
    grid_id = non_protected_ids,
    is_protected = FALSE,
    protected_frac = 0,
    wdpa_pid = NA_integer_,
    iucn_cat = NA_character_,
    desig_year = NA_integer_,
    gov_type = NA_character_
  )

  protected_cells <- rbind(protected_cells, non_protected, fill = TRUE)
}

# ==============================================================================
# PREPARE OUTPUT
# ==============================================================================

log_message("Preparing output...")

# Add tile_id
protected_cells[, tile_id := tile_id]

# Ensure all columns exist
if (!"wdpa_pid" %in% names(protected_cells)) protected_cells[, wdpa_pid := NA_integer_]
if (!"iucn_cat" %in% names(protected_cells)) protected_cells[, iucn_cat := NA_character_]
if (!"desig_year" %in% names(protected_cells)) protected_cells[, desig_year := NA_integer_]
if (!"gov_type" %in% names(protected_cells)) protected_cells[, gov_type := NA_character_]

# Select final columns
output_data <- protected_cells[, .(
  grid_id,
  tile_id,
  is_protected,
  protected_frac,
  wdpa_pid,
  iucn_cat,
  desig_year,
  gov_type
)]

setorder(output_data, grid_id)

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

log_message("Summary statistics:")
log_message(sprintf("  Total cells: %s", format(nrow(output_data), big.mark = ",")))
log_message(sprintf("  Protected cells: %s (%.1f%%)",
                    format(sum(output_data$is_protected), big.mark = ","),
                    100 * mean(output_data$is_protected)))

if (sum(output_data$is_protected) > 0) {
  protected <- output_data[is_protected == TRUE]
  log_message(sprintf("  Mean protected fraction: %.3f",
                      mean(protected$protected_frac)))
  log_message(sprintf("  Fully protected cells (>99%%): %s",
                      format(sum(protected$protected_frac > 0.99), big.mark = ",")))

  # IUCN category distribution
  iucn_dist <- protected[!is.na(iucn_cat), .N, by = iucn_cat]
  setorder(iucn_dist, -N)
  log_message("IUCN category distribution:")
  for (i in seq_len(min(nrow(iucn_dist), 6))) {
    log_message(sprintf("    %s: %s", iucn_dist$iucn_cat[i],
                        format(iucn_dist$N[i], big.mark = ",")))
  }
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
