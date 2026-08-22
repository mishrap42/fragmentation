# ==============================================================================
# VALIDATE TMF DATA EXTRACTION PIPELINE
# ==============================================================================
# Produces validation outputs:
# [1] Summary stats on land use shares by country (single cross-section)
# [2] Global time series of land use shares
# [3] Map of undisturbed forest share overlaid on country boundaries
# [4] Maps of deforestation, degradation, and regrowth (final year)
# [5] Protected area maps (overlay + undisturbed forest in protected areas)
#
# IMPORTANT: Processes data TILE-BY-TILE to avoid R's row limit (~2.1B rows)
# ==============================================================================

.libPaths(c('/dartfs-hpc/rc/lab/M/MishraP/Rlibs/4.4.0/', .libPaths()))

library(here)
here::i_am("code/build/validate_pipeline.R")
source("code/build/BUILD_workspace.R")

library(ggplot2)
library(sf)
library(viridis)
library(ggnewscale)  # For new_scale_color() in protected area overlay map

cat(rep("=", 71), "\n", sep = "")
cat("TMF PIPELINE VALIDATION (Chunked Processing)\n")
cat(rep("=", 71), "\n", sep = "")

# Create output directory
output_dir <- file.path(project_root, "output", "validation")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# IDENTIFY AVAILABLE TILES
# ==============================================================================

cat("\nIdentifying available data...\n")

tile_files <- sapply(1:N_TMF_TILES, get_consolidated_tile_filename)
files_exist <- file.exists(tile_files)
available_tiles <- which(files_exist)
n_tiles <- length(available_tiles)

cat(sprintf("Found %d/%d consolidated tile files\n", n_tiles, N_TMF_TILES))

if (n_tiles == 0) {
  stop("No consolidated tile files found. Run Stage 2a first.")
}

# ==============================================================================
# [0] TILE INDEX MAP - Visual reference for tile numbering
# ==============================================================================

cat("\nGenerating tile index map...\n")

# Create tile polygons from TMF_TILE_INDEX
tile_polys <- lapply(1:N_TMF_TILES, function(i) {

  tile <- TMF_TILE_INDEX[i]
  coords <- matrix(c(
    tile$lon_west, tile$lat_south,
    tile$lon_east, tile$lat_south,
    tile$lon_east, tile$lat_north,
    tile$lon_west, tile$lat_north,
    tile$lon_west, tile$lat_south
  ), ncol = 2, byrow = TRUE)
  st_polygon(list(coords))
})
tiles_sf <- st_sf(
  tile_id = 1:N_TMF_TILES,
  available = 1:N_TMF_TILES %in% available_tiles,
  geometry = st_sfc(tile_polys, crs = 4326)
)

# Calculate tile centers for labels
tile_centers <- st_centroid(tiles_sf)
tile_centers_coords <- st_coordinates(tile_centers)
tiles_sf$label_x <- tile_centers_coords[, 1]
tiles_sf$label_y <- tile_centers_coords[, 2]

# Load world boundaries
world <- st_read(gadm_path, layer = "ADM_0", quiet = TRUE)
world <- st_transform(world, 4326)

# Create tile index map
p_tiles <- ggplot() +

  geom_sf(data = world, fill = "gray90", color = "gray70", linewidth = 0.1) +
  geom_sf(data = tiles_sf, aes(fill = available), color = "black", linewidth = 0.3, alpha = 0.5) +
  scale_fill_manual(values = c("TRUE" = "forestgreen", "FALSE" = "red"),
                    labels = c("TRUE" = "Available", "FALSE" = "Missing"),
                    name = "Status") +
  geom_text(data = tiles_sf, aes(x = label_x, y = label_y, label = tile_id),
            size = 2.5, fontface = "bold") +
  coord_sf(xlim = c(-180, 180), ylim = c(-35, 35)) +
  labs(title = "TMF Tile Index Map",
       subtitle = sprintf("%d tiles total, %d available", N_TMF_TILES, n_tiles),
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "00_tile_index_map.png"), p_tiles,
       width = 16, height = 8, dpi = 150)
cat(sprintf("Saved: %s\n", file.path(output_dir, "00_tile_index_map.png")))

rm(tile_polys, tiles_sf, tile_centers, world)
gc()

# ==============================================================================
# LOAD GRID METADATA (smaller - can load all at once)
# ==============================================================================

cat("\nLoading grid metadata...\n")

grid_files <- sapply(1:N_SUB_TILES, function(s) get_grid_filename(s, "parquet"))
files_exist_g <- file.exists(grid_files)
cat(sprintf("Found %d/%d grid files\n", sum(files_exist_g), N_SUB_TILES))

grid_list <- lapply(grid_files[files_exist_g], function(f) {
  dt <- arrow::read_parquet(f)
  setDT(dt)
  dt[, .(grid_id, centroid_lon, centroid_lat, country_iso3, country_name, area_km2)]
})
grid_all <- rbindlist(grid_list, fill = TRUE)
rm(grid_list)
gc()

cat(sprintf("Loaded %s grid cells with metadata\n", format(nrow(grid_all), big.mark = ",")))

# Create lookup tables for efficiency
grid_area <- grid_all[, .(grid_id, area_km2)]
setkey(grid_area, grid_id)

grid_country <- grid_all[, .(grid_id, country_iso3, country_name)]
setkey(grid_country, grid_id)

grid_coords <- grid_all[, .(grid_id, centroid_lon, centroid_lat)]

# Check for duplicate grid_ids (critical bug if present)
dup_ids <- grid_coords[, .N, by = grid_id][N > 1]
if (nrow(dup_ids) > 0) {
  cat(sprintf("\n*** WARNING: %d duplicate grid_ids found! ***\n", nrow(dup_ids)))
  cat("Sample duplicates:\n")
  print(head(dup_ids, 10))
  cat("\nDe-duplicating by taking first occurrence...\n")
  grid_coords <- unique(grid_coords, by = "grid_id")
}

setkey(grid_coords, grid_id)

# ==============================================================================
# PASS 1: AGGREGATE BY COUNTRY AND YEAR (TILE BY TILE)
# ==============================================================================

cat("\n", rep("=", 71), "\n", sep = "")
cat("PASS 1: Computing aggregations tile-by-tile...\n")
cat(rep("=", 71), "\n", sep = "")

# Initialize accumulators
country_year_class_agg <- list()
global_year_class_agg <- list()
unique_grids <- c()
year_range <- c(Inf, -Inf)
all_classes <- c()

for (i in seq_along(available_tiles)) {
  tile_id <- available_tiles[i]
  tile_file <- tile_files[tile_id]

  if (i %% 10 == 1 || i == length(available_tiles)) {
    cat(sprintf("  Processing tile %d/%d (tile_id=%d)...\n", i, n_tiles, tile_id))
  }

  # Load this tile's TMF data
  tmf_tile <- arrow::read_parquet(tile_file)
  setDT(tmf_tile)

  # Track metadata
  unique_grids <- union(unique_grids, unique(tmf_tile$grid_id))
  year_range[1] <- min(year_range[1], min(tmf_tile$year))
  year_range[2] <- max(year_range[2], max(tmf_tile$year))
  all_classes <- union(all_classes, unique(tmf_tile$tmf_class))

  # Merge with area and country
  tmf_tile <- merge(tmf_tile, grid_area, by = "grid_id", all.x = TRUE)
  tmf_tile <- merge(tmf_tile, grid_country, by = "grid_id", all.x = TRUE)

  # Aggregate by country, year, class
  tile_country_agg <- tmf_tile[!is.na(country_iso3), .(
    area_km2 = sum(fraction * area_km2, na.rm = TRUE),
    n_cells = length(unique(grid_id))
  ), by = .(country_iso3, country_name, year, tmf_class)]

  country_year_class_agg[[i]] <- tile_country_agg

  # Aggregate by year, class only (for global time series)
  tile_global_agg <- tmf_tile[, .(
    area_km2 = sum(fraction * area_km2, na.rm = TRUE),
    n_cells = length(unique(grid_id))
  ), by = .(year, tmf_class)]

  global_year_class_agg[[i]] <- tile_global_agg

  rm(tmf_tile, tile_country_agg, tile_global_agg)
}

gc()

# Combine aggregations
cat("\nCombining aggregations...\n")

country_agg_all <- rbindlist(country_year_class_agg)
country_summary <- country_agg_all[, .(
  total_area_km2 = sum(area_km2, na.rm = TRUE),
  n_cells = sum(n_cells)
), by = .(country_iso3, country_name, year, tmf_class)]
rm(country_agg_all)

global_agg_all <- rbindlist(global_year_class_agg)
global_ts <- global_agg_all[, .(
  total_area_km2 = sum(area_km2, na.rm = TRUE),
  n_cells = sum(n_cells)
), by = .(year, tmf_class)]
rm(global_agg_all)

cat(sprintf("Unique grid cells: %s\n", format(length(unique_grids), big.mark = ",")))
cat(sprintf("Years: %d to %d\n", year_range[1], year_range[2]))
cat(sprintf("TMF classes: %s\n", paste(sort(all_classes), collapse = ", ")))

# ==============================================================================
# [1] COUNTRY-LEVEL SUMMARY STATS (LATEST YEAR CROSS-SECTION)
# ==============================================================================

cat("\n", rep("=", 71), "\n", sep = "")
cat("[OUTPUT 1] COUNTRY-LEVEL LAND USE SUMMARY (Latest Year)\n")
cat(rep("=", 71), "\n", sep = "")

latest_year <- year_range[2]
cat(sprintf("Cross-section year: %d\n\n", latest_year))

# Filter to latest year
country_latest <- country_summary[year == latest_year]

# Reshape to wide format
country_wide <- dcast(country_latest, country_iso3 + country_name ~ tmf_class,
                      value.var = "total_area_km2", fill = 0)

# Add cell counts
cell_counts <- country_latest[, .(n_cells = max(n_cells)), by = country_iso3]
country_wide <- merge(country_wide, cell_counts, by = "country_iso3", all.x = TRUE)

# Calculate totals and percentages
total_cols <- setdiff(names(country_wide), c("country_iso3", "country_name", "n_cells"))
country_wide[, total_area := rowSums(.SD, na.rm = TRUE), .SDcols = total_cols]

if ("Undisturbed_TMF" %in% names(country_wide)) {
  country_wide[, undisturbed_pct := 100 * Undisturbed_TMF / total_area]
}

setorder(country_wide, -total_area)

# Print top 20 countries
cat("Top 20 countries by total TMF coverage:\n")
cat(rep("-", 100), "\n", sep = "")
print(head(country_wide[, .(country_iso3, country_name, n_cells,
                            total_area = round(total_area, 0),
                            undisturbed_pct = round(undisturbed_pct, 1))], 20))

# Save full table
country_file <- file.path(output_dir, sprintf("country_landuse_%d.csv", latest_year))
fwrite(country_wide, country_file)
cat(sprintf("\nFull table saved to: %s\n", country_file))

# ==============================================================================
# [2] GLOBAL TIME SERIES OF LAND USE SHARES
# ==============================================================================

cat("\n", rep("=", 71), "\n", sep = "")
cat("[OUTPUT 2] GLOBAL TIME SERIES OF LAND USE SHARES\n")
cat(rep("=", 71), "\n", sep = "")

# Calculate percentage of total for each year
year_totals <- global_ts[, .(year_total = sum(total_area_km2)), by = year]
global_ts <- merge(global_ts, year_totals, by = "year")
global_ts[, pct_of_total := 100 * total_area_km2 / year_total]

setorder(global_ts, year, tmf_class)

# Print summary
cat("\nGlobal land use time series (area in million km2):\n")
cat(rep("-", 80), "\n", sep = "")

ts_wide <- dcast(global_ts, year ~ tmf_class, value.var = "total_area_km2")
ts_wide_print <- copy(ts_wide)
ts_wide_print[, (setdiff(names(ts_wide_print), "year")) := lapply(.SD, function(x) round(x/1e6, 3)),
              .SDcols = setdiff(names(ts_wide_print), "year")]
print(ts_wide_print)

# Save time series
ts_file <- file.path(output_dir, "global_timeseries.csv")
fwrite(global_ts, ts_file)
cat(sprintf("\nTime series saved to: %s\n", ts_file))

# Create time series plot
cat("\nCreating time series plot...\n")

class_colors <- c(
  "Undisturbed_TMF" = "#1a9641",
  "Degraded_TMF" = "#a6d96a",
  "Deforested" = "#d73027",
  "Regrowth" = "#4575b4",
  "Water" = "#74add1",
  "Other_land" = "#d9d9d9"
)

main_classes <- c("Undisturbed_TMF", "Degraded_TMF", "Deforested", "Regrowth")
plot_data_ts <- global_ts[tmf_class %in% main_classes]

p_ts <- ggplot(plot_data_ts,
               aes(x = year, y = total_area_km2 / 1e6, color = tmf_class, group = tmf_class)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(values = class_colors, name = "Land Use Class") +
  scale_x_continuous(breaks = seq(1990, 2025, 5)) +
  labs(
    title = "Global Tropical Moist Forest Land Use Time Series (1990-2023)",
    subtitle = sprintf("Data from %s grid cells across %d TMF tiles",
                       format(length(unique_grids), big.mark = ","),
                       n_tiles),
    x = "Year",
    y = "Area (million km²)",
    caption = "Source: JRC Tropical Moist Forest dataset"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ts_plot_file <- file.path(output_dir, "global_timeseries.png")
ggsave(ts_plot_file, p_ts, width = 12, height = 7, dpi = 300)
cat(sprintf("Time series plot saved to: %s\n", ts_plot_file))

# ==============================================================================
# PASS 2: SAMPLE POINTS FOR MAPS (TILE BY TILE)
# ==============================================================================

cat("\n", rep("=", 71), "\n", sep = "")
cat("PASS 2: Sampling points for maps (tile-by-tile)...\n")
cat(rep("=", 71), "\n", sep = "")

# Target sample sizes per tile
max_total_points <- 500000
points_per_tile <- ceiling(max_total_points / n_tiles)

# Expected TMF columns (some tiles may not have all classes)
expected_cols <- c("Undisturbed_TMF", "Degraded_TMF", "Deforested", "Regrowth")

map_samples <- list()

for (i in seq_along(available_tiles)) {
  tile_id <- available_tiles[i]
  tile_file <- tile_files[tile_id]

  # Load this tile's TMF data for latest year only
  tmf_tile <- arrow::read_parquet(tile_file)
  setDT(tmf_tile)
  tmf_latest <- tmf_tile[year == latest_year]
  rm(tmf_tile)

  if (nrow(tmf_latest) == 0) next

  # DEBUG: Print class distribution for this tile
  tile_class_dist <- tmf_latest[, .(mean_frac = mean(fraction, na.rm = TRUE)), by = tmf_class]
  regrowth_mean <- tile_class_dist[tmf_class == "Regrowth"]$mean_frac
  undist_mean <- tile_class_dist[tmf_class == "Undisturbed_TMF"]$mean_frac
  if (length(regrowth_mean) == 0) regrowth_mean <- 0
  if (length(undist_mean) == 0) undist_mean <- 0

  # Flag tiles with suspicious data
  if (regrowth_mean > 0.5 || undist_mean < 0.01) {
    cat(sprintf("  ** TILE %d: Regrowth=%.1f%%, Undisturbed=%.1f%% **\n",
                tile_id, 100*regrowth_mean, 100*undist_mean))
  } else if (i %% 10 == 1 || i == length(available_tiles)) {
    cat(sprintf("  Tile %d/%d (id=%d): Regrowth=%.1f%%, Undisturbed=%.1f%%\n",
                i, n_tiles, tile_id, 100*regrowth_mean, 100*undist_mean))
  }

  # Reshape to wide format for this tile (average duplicates from geometry splitting)
  tile_wide <- dcast(tmf_latest, grid_id ~ tmf_class, value.var = "fraction",
                     fill = 0, fun.aggregate = mean)

  # Merge with coordinates
  tile_wide <- merge(tile_wide, grid_coords, by = "grid_id", all.x = TRUE)

  # Ensure expected columns exist (some tiles may not have all classes)
  for (col in expected_cols) {
    if (!col %in% names(tile_wide)) {
      tile_wide[, (col) := 0]
    }
  }

  # Sample from this tile
  n_sample <- min(nrow(tile_wide), points_per_tile)
  if (nrow(tile_wide) > n_sample) {
    # Prioritize cells with forest
    has_forest <- tile_wide[Undisturbed_TMF > 0.1 | Degraded_TMF > 0.1 |
                            Deforested > 0.1 | Regrowth > 0.1]
    no_forest <- tile_wide[Undisturbed_TMF <= 0.1 & Degraded_TMF <= 0.1 &
                           Deforested <= 0.1 & Regrowth <= 0.1]

    n_forest <- min(nrow(has_forest), ceiling(n_sample * 0.8))
    n_other <- min(nrow(no_forest), n_sample - n_forest)

    tile_sample <- rbind(
      if (n_forest > 0) has_forest[sample(.N, n_forest)] else NULL,
      if (n_other > 0) no_forest[sample(.N, n_other)] else NULL
    )
  } else {
    tile_sample <- tile_wide
  }

  map_samples[[i]] <- tile_sample
  rm(tmf_latest, tile_wide, tile_sample)
}

gc()

map_data <- rbindlist(map_samples, fill = TRUE)
rm(map_samples)

cat(sprintf("Total sampled points: %s\n", format(nrow(map_data), big.mark = ",")))

# Ensure all expected columns exist in combined data
for (col in expected_cols) {
  if (!col %in% names(map_data)) {
    map_data[, (col) := 0]
  }
}

# ==============================================================================
# [3] MAP OF UNDISTURBED FOREST SHARE (LATEST YEAR)
# ==============================================================================

cat("\n", rep("=", 71), "\n", sep = "")
cat("[OUTPUT 3] MAP OF UNDISTURBED FOREST SHARE\n")
cat(rep("=", 71), "\n", sep = "")

# Stats
n_with_forest <- sum(map_data$Undisturbed_TMF > 0, na.rm = TRUE)
mean_forest <- mean(map_data[Undisturbed_TMF > 0]$Undisturbed_TMF, na.rm = TRUE)

cat(sprintf("Sampled cells with undisturbed forest: %s\n", format(n_with_forest, big.mark = ",")))
cat(sprintf("Mean undisturbed fraction (where >0): %.1f%%\n", 100 * mean_forest))

# Load country boundaries from GADM
cat("Loading country boundaries...\n")
gadm <- st_read(gadm_path, layer = "ADM_0", quiet = TRUE)
gadm <- st_transform(gadm, WGS84_CRS)

# Filter to tropical countries
sf_use_s2(FALSE)
gadm_centroids <- st_coordinates(st_centroid(gadm))
sf_use_s2(TRUE)
gadm_tropical <- gadm[gadm_centroids[,2] > -35 & gadm_centroids[,2] < 35, ]

cat(sprintf("Creating undisturbed forest map (%s points)...\n",
            format(sum(map_data$Undisturbed_TMF > 0), big.mark = ",")))

p_undisturbed <- ggplot() +
  geom_sf(data = gadm_tropical, fill = "gray95", color = "gray60", linewidth = 0.2) +
  geom_point(data = map_data[Undisturbed_TMF > 0],
             aes(x = centroid_lon, y = centroid_lat, color = Undisturbed_TMF),
             size = 0.1, alpha = 0.7) +
  scale_color_viridis(option = "G", direction = -1,
                      limits = c(0, 1),
                      labels = scales::percent,
                      name = "Undisturbed\nForest") +
  coord_sf(xlim = c(-180, 180), ylim = c(-30, 30), expand = FALSE) +
  labs(
    title = sprintf("Undisturbed Tropical Moist Forest Cover (%d)", latest_year),
    subtitle = sprintf("%s sampled cells | Mean coverage: %.1f%% (where present)",
                       format(n_with_forest, big.mark = ","),
                       100 * mean_forest),
    x = "Longitude", y = "Latitude",
    caption = "Source: JRC Tropical Moist Forest dataset | 1km grid cells"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "gray40"),
    legend.position = "right",
    panel.background = element_rect(fill = "aliceblue", color = NA),
    panel.grid = element_line(color = "gray90")
  )

undisturbed_file <- file.path(output_dir, sprintf("map_undisturbed_%d.png", latest_year))
ggsave(undisturbed_file, p_undisturbed, width = 16, height = 8, dpi = 300)
cat(sprintf("Undisturbed forest map saved to: %s\n", undisturbed_file))

# ==============================================================================
# [4] MAPS OF DEFORESTATION, DEGRADATION, REGROWTH (LATEST YEAR)
# ==============================================================================

cat("\n", rep("=", 71), "\n", sep = "")
cat("[OUTPUT 4] MAPS OF DEFORESTATION, DEGRADATION, REGROWTH\n")
cat(rep("=", 71), "\n", sep = "")

# Print summary stats from sampled data
cat(sprintf("\nDeforestation (%d) - from sample:\n", latest_year))
cat(sprintf("  Cells with deforestation: %s\n",
            format(sum(map_data$Deforested > 0, na.rm = TRUE), big.mark = ",")))
cat(sprintf("  Mean fraction (where >0): %.1f%%\n",
            100 * mean(map_data[Deforested > 0]$Deforested, na.rm = TRUE)))

cat(sprintf("\nDegradation (%d) - from sample:\n", latest_year))
cat(sprintf("  Cells with degradation: %s\n",
            format(sum(map_data$Degraded_TMF > 0, na.rm = TRUE), big.mark = ",")))
cat(sprintf("  Mean fraction (where >0): %.1f%%\n",
            100 * mean(map_data[Degraded_TMF > 0]$Degraded_TMF, na.rm = TRUE)))

cat(sprintf("\nRegrowth (%d) - from sample:\n", latest_year))
cat(sprintf("  Cells with regrowth: %s\n",
            format(sum(map_data$Regrowth > 0, na.rm = TRUE), big.mark = ",")))
cat(sprintf("  Mean fraction (where >0): %.1f%%\n",
            100 * mean(map_data[Regrowth > 0]$Regrowth, na.rm = TRUE)))

# Function to create change map
create_change_map <- function(data, var_name, color_scale, title, subtitle) {
  data_filtered <- data[get(var_name) > 0]

  ggplot() +
    geom_sf(data = gadm_tropical, fill = "gray95", color = "gray60", linewidth = 0.2) +
    geom_point(data = data_filtered,
               aes(x = centroid_lon, y = centroid_lat, color = get(var_name)),
               size = 0.15, alpha = 0.8) +
    scale_color_viridis(option = color_scale, direction = -1,
                        limits = c(0, 1),
                        labels = scales::percent,
                        name = "Cell\nFraction") +
    coord_sf(xlim = c(-180, 180), ylim = c(-30, 30), expand = FALSE) +
    labs(title = title, subtitle = subtitle,
         x = "Longitude", y = "Latitude",
         caption = "Source: JRC Tropical Moist Forest dataset") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "gray40"),
      legend.position = "right",
      panel.background = element_rect(fill = "aliceblue", color = NA),
      panel.grid = element_line(color = "gray90")
    )
}

# Deforestation map
cat("\nCreating deforestation map...\n")
p_deforest <- create_change_map(
  map_data, "Deforested", "B",
  sprintf("Deforested Land (%d)", latest_year),
  sprintf("%s sampled cells with deforestation",
          format(sum(map_data$Deforested > 0, na.rm = TRUE), big.mark = ","))
)
ggsave(file.path(output_dir, sprintf("map_deforested_%d.png", latest_year)),
       p_deforest, width = 16, height = 8, dpi = 300)

# Degradation map
cat("Creating degradation map...\n")
p_degraded <- create_change_map(
  map_data, "Degraded_TMF", "E",
  sprintf("Degraded Forest (%d)", latest_year),
  sprintf("%s sampled cells with degradation",
          format(sum(map_data$Degraded_TMF > 0, na.rm = TRUE), big.mark = ","))
)
ggsave(file.path(output_dir, sprintf("map_degraded_%d.png", latest_year)),
       p_degraded, width = 16, height = 8, dpi = 300)

# Regrowth map
cat("Creating regrowth map...\n")
p_regrowth <- create_change_map(
  map_data, "Regrowth", "D",
  sprintf("Forest Regrowth (%d)", latest_year),
  sprintf("%s sampled cells with regrowth",
          format(sum(map_data$Regrowth > 0, na.rm = TRUE), big.mark = ","))
)
ggsave(file.path(output_dir, sprintf("map_regrowth_%d.png", latest_year)),
       p_regrowth, width = 16, height = 8, dpi = 300)

cat(sprintf("\nChange maps saved to: %s/\n", output_dir))

# ==============================================================================
# [5] PROTECTED AREA MAPS
# ==============================================================================

cat("\n", rep("=", 71), "\n", sep = "")
cat("[OUTPUT 5] PROTECTED AREA MAPS\n")
cat(rep("=", 71), "\n", sep = "")

# Load protection status from grid files
cat("Loading protection status from grid files...\n")
grid_protection_list <- lapply(grid_files[files_exist_g], function(f) {
  dt <- arrow::read_parquet(f)
  setDT(dt)
  if ("is_protected" %in% names(dt)) {
    dt[, .(grid_id, is_protected)]
  } else {
    NULL
  }
})
grid_protection <- rbindlist(grid_protection_list[!sapply(grid_protection_list, is.null)], fill = TRUE)
rm(grid_protection_list)
setkey(grid_protection, grid_id)

cat(sprintf("Loaded protection status for %s grid cells\n",
            format(nrow(grid_protection), big.mark = ",")))
cat(sprintf("Protected cells: %s (%.1f%%)\n",
            format(sum(grid_protection$is_protected, na.rm = TRUE), big.mark = ","),
            100 * mean(grid_protection$is_protected, na.rm = TRUE)))

# Merge protection status with map data
map_data <- merge(map_data, grid_protection, by = "grid_id", all.x = TRUE)
map_data[is.na(is_protected), is_protected := FALSE]

# Stats
n_protected <- sum(map_data$is_protected, na.rm = TRUE)
n_unprotected <- sum(!map_data$is_protected, na.rm = TRUE)
cat(sprintf("\nIn sampled map data:\n"))
cat(sprintf("  Protected cells: %s\n", format(n_protected, big.mark = ",")))
cat(sprintf("  Unprotected cells: %s\n", format(n_unprotected, big.mark = ",")))

# [5a] Map showing protected areas highlighted in color
cat("\nCreating protected areas overlay map...\n")

p_protected_overlay <- ggplot() +
  geom_sf(data = gadm_tropical, fill = "gray95", color = "gray60", linewidth = 0.2) +
  # Plot unprotected first (background)
  geom_point(data = map_data[is_protected == FALSE & Undisturbed_TMF > 0],
             aes(x = centroid_lon, y = centroid_lat, color = Undisturbed_TMF),
             size = 0.1, alpha = 0.5) +
  scale_color_gradient(low = "lightgray", high = "darkgray",
                       limits = c(0, 1), guide = "none") +
  # New scale for protected areas
  ggnewscale::new_scale_color() +
  geom_point(data = map_data[is_protected == TRUE],
             aes(x = centroid_lon, y = centroid_lat, color = Undisturbed_TMF),
             size = 0.2, alpha = 0.9) +
  scale_color_viridis(option = "G", direction = -1,
                      limits = c(0, 1),
                      labels = scales::percent,
                      name = "Protected Area\nUndisturbed Forest") +
  coord_sf(xlim = c(-180, 180), ylim = c(-30, 30), expand = FALSE) +
  labs(
    title = sprintf("Protected Areas Highlighted (%d)", latest_year),
    subtitle = sprintf("%s protected cells (colored) | %s unprotected cells (gray)",
                       format(n_protected, big.mark = ","),
                       format(n_unprotected, big.mark = ",")),
    x = "Longitude", y = "Latitude",
    caption = "Protected areas from WDPA | Unprotected shown in gray"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "gray40"),
    legend.position = "right",
    panel.background = element_rect(fill = "aliceblue", color = NA),
    panel.grid = element_line(color = "gray90")
  )

protected_overlay_file <- file.path(output_dir, sprintf("map_protected_overlay_%d.png", latest_year))
ggsave(protected_overlay_file, p_protected_overlay, width = 16, height = 8, dpi = 300)
cat(sprintf("Protected areas overlay map saved to: %s\n", protected_overlay_file))

# [5b] Map showing ONLY protected areas, shaded by undisturbed forest fraction
cat("\nCreating protected areas only map (by undisturbed forest fraction)...\n")

# Calculate mean undisturbed forest in protected vs unprotected
mean_undist_protected <- mean(map_data[is_protected == TRUE]$Undisturbed_TMF, na.rm = TRUE)
mean_undist_unprotected <- mean(map_data[is_protected == FALSE]$Undisturbed_TMF, na.rm = TRUE)
cat(sprintf("Mean undisturbed forest fraction:\n"))
cat(sprintf("  Protected areas: %.1f%%\n", 100 * mean_undist_protected))
cat(sprintf("  Unprotected areas: %.1f%%\n", 100 * mean_undist_unprotected))

p_protected_only <- ggplot() +
  # Gray background for unprotected areas
  geom_sf(data = gadm_tropical, fill = "gray90", color = "gray70", linewidth = 0.2) +
  geom_point(data = map_data[is_protected == FALSE],
             aes(x = centroid_lon, y = centroid_lat),
             color = "gray80", size = 0.05, alpha = 0.3) +
  # Protected areas colored by undisturbed forest
  geom_point(data = map_data[is_protected == TRUE],
             aes(x = centroid_lon, y = centroid_lat, color = Undisturbed_TMF),
             size = 0.25, alpha = 0.9) +
  scale_color_viridis(option = "G", direction = -1,
                      limits = c(0, 1),
                      labels = scales::percent,
                      name = "Undisturbed\nForest") +
  coord_sf(xlim = c(-180, 180), ylim = c(-30, 30), expand = FALSE) +
  labs(
    title = sprintf("Protected Areas: Undisturbed Forest Coverage (%d)", latest_year),
    subtitle = sprintf("Protected: %.1f%% undisturbed | Unprotected: %.1f%% undisturbed",
                       100 * mean_undist_protected, 100 * mean_undist_unprotected),
    x = "Longitude", y = "Latitude",
    caption = "Gray regions = unprotected | Colored = protected areas (WDPA)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "gray40"),
    legend.position = "right",
    panel.background = element_rect(fill = "aliceblue", color = NA),
    panel.grid = element_line(color = "gray90")
  )

protected_only_file <- file.path(output_dir, sprintf("map_protected_undisturbed_%d.png", latest_year))
ggsave(protected_only_file, p_protected_only, width = 16, height = 8, dpi = 300)
cat(sprintf("Protected areas undisturbed map saved to: %s\n", protected_only_file))

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n", rep("=", 71), "\n", sep = "")
cat("VALIDATION COMPLETE\n")
cat(rep("=", 71), "\n", sep = "")

cat(sprintf("\nOutput files in: %s/\n", output_dir))
cat("  - country_landuse_YYYY.csv        : Country-level land use summary\n")
cat("  - global_timeseries.csv           : Global time series data\n")
cat("  - global_timeseries.png           : Time series plot\n")
cat("  - map_undisturbed_YYYY.png        : Undisturbed forest map\n")
cat("  - map_deforested_YYYY.png         : Deforestation map\n")
cat("  - map_degraded_YYYY.png           : Degradation map\n")
cat("  - map_regrowth_YYYY.png           : Regrowth map\n")
cat("  - map_protected_overlay_YYYY.png  : Protected areas highlighted (unprotected gray)\n")
cat("  - map_protected_undisturbed_YYYY.png : Protected areas by undisturbed fraction\n")

cat("\n", rep("=", 71), "\n", sep = "")
cat("DATA QUALITY SUMMARY\n")
cat(rep("=", 71), "\n", sep = "")
cat(sprintf("Unique grid cells: %s\n", format(length(unique_grids), big.mark = ",")))
cat(sprintf("Countries covered: %d\n", length(unique(grid_all$country_iso3))))
cat(sprintf("Years covered: %d-%d (%d years)\n",
            year_range[1], year_range[2],
            year_range[2] - year_range[1] + 1))
cat(sprintf("TMF classes found: %s\n", paste(sort(all_classes), collapse = ", ")))

cat("\nDone.\n")
