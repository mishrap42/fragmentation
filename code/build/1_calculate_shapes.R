library(here)        # Manage file paths relative to project root
here::i_am("code/build/1_calculate_shapes.R")
library(sf)          # Spatial data handling (reading shapefiles, geometric operations)
library(dplyr)       # Data manipulation (filtering, grouping, summarizing)

set.seed(42)
# This stores your repository path as a function "here()"
if (grepl("ahaanj", here())) {
  dropbox_dir <- "/Users/ahaanj/Dartmouth College Dropbox/Ahaan Jindal/Protected Area Fragmentation/"
  i = 1
} else if (grepl("mishrap", here())) {
  dropbox_dir <- "/Users/mishrap/Dropbox (Personal)/Protected Area Fragmentation/"
  i = 1
} else {
  # Suppress summarise info
  options(dplyr.summarise.inform = FALSE)

  # Set paths
  path <- here()
  dropbox_dir <- file.path(here(), 'Protected Area Fragmentation')

  # Get task ID from SLURM array
  i <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID', '1'))

  print(paste0('Processing chunk ', i))

}

intermediate_data_dir = file.path(dropbox_dir, 'Data/build')
figure_dir = file.path(dropbox_dir, 'Figures', '1_calculate_shapes')
if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

# Read in the shapefile of protected areas
pa <- st_read(file.path(dropbox_dir, 'Data/raw/protected sites/protectedsites.shp'))

# print('Country counts: validation of sampling')
# pa %>% st_drop_geometry() %>% group_by(cddaCountr) %>% summarise(n = n()) %>% print(n = Inf)

# print('Foundation year counts:')
# pa %>% st_drop_geometry() %>% group_by(legalFound) %>% summarise(n = n()) %>% print(n = Inf)

print('Checking shape validity:')
pa$valid <- st_is_valid(pa)
sprintf('Number of valid geometries: %d out of %d\n', sum( (pa %>% st_drop_geometry())$valid), nrow(pa))

# ========== NINA'S SHAPE METRICS - CORRECTED IMPLEMENTATION ==========
# All require sampling points and computing pairwise distances

# Function to sample interior points from a polygon
sample_interior_points <- function(geom, n_points = 100) {
  # Generate random points inside the polygon
  points <- st_sample(geom, size = n_points, type = "random")
  return(points)
}

# Function to process one chunk
chunk_size <- 100
n_total <- nrow(pa)
n_chunks <- ceiling(n_total / chunk_size)
pa_shapes <- list()

start_idx <- (i - 1) * chunk_size + 1
end_idx <- min(i * chunk_size, n_total)

cat(sprintf("Processing chunk %d of %d (rows %d to %d)...\n", i, n_chunks, start_idx, end_idx))

chunk = pa[start_idx:end_idx, ]
chunk_result <- chunk %>%
  rowwise() %>%
  mutate(
    interior_pts = list(sample_interior_points(geometry, n_points = 100))
  ) %>%
  mutate(
    dist_matrix = list(units::set_units(st_distance(interior_pts, interior_pts), "km")),
    dist_to_centroid = list(units::set_units(st_distance(interior_pts, st_centroid(geometry)), "km"))
  ) %>%
  mutate(
    disconnection_index = as.numeric(sum(dist_matrix, na.rm = TRUE) / (nrow(dist_matrix) * (nrow(dist_matrix) - 1))),
    remoteness = as.numeric(mean(dist_to_centroid, na.rm = TRUE)),
    spin = as.numeric(mean(dist_to_centroid^2, na.rm = TRUE)),
    range_km = as.numeric(units::set_units(st_length(st_cast(st_convex_hull(geometry), "LINESTRING")), "km"))
  ) %>%
  select(-interior_pts, -dist_matrix, -dist_to_centroid) %>%
  ungroup()

pa_shapes[[i]] <- chunk_result


pa_shapes <- bind_rows(pa_shapes)

st_write(pa_shapes, file.path(intermediate_data_dir, paste0('pa_', i, '.shp')), delete_dsn = TRUE)
