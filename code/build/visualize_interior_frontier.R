# ==============================================================================
# VISUALIZE FOREST INTERIOR AND FRONTIER DISTRIBUTION
# ==============================================================================
# Creates a global map showing interior (core forest) and frontier zones
# ==============================================================================

# Add lab library path (before loading any packages)
.libPaths(c('/dartfs-hpc/rc/lab/M/MishraP/Rlibs/4.4.0/', .libPaths()))

library(here)
here::i_am("code/build/visualize_interior_frontier.R")
source("code/build/BUILD_workspace.R")

library(ggplot2)
library(sf)

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading interior classifications...\n")

interior_files <- sapply(1:N_TMF_TILES, get_interior_filename)
files_exist <- file.exists(interior_files)
cat(sprintf("Found %d/%d interior files\n", sum(files_exist), N_TMF_TILES))

if (sum(files_exist) == 0) {
  stop("No interior classification files found. Run Stage 3 first.")
}

interior_list <- lapply(interior_files[files_exist], function(f) {
  dt <- read_parquet(f)
  setDT(dt)
  dt
})
interior_all <- rbindlist(interior_list, fill = TRUE)
rm(interior_list)

cat(sprintf("Loaded %s interior classifications\n",
            format(nrow(interior_all), big.mark = ",")))

# Load frontier data
cat("Loading frontier classifications...\n")

frontier_files <- sapply(1:N_TMF_TILES, get_frontier_filename)
files_exist_f <- file.exists(frontier_files)
cat(sprintf("Found %d/%d frontier files\n", sum(files_exist_f), N_TMF_TILES))

if (sum(files_exist_f) > 0) {
  frontier_list <- lapply(frontier_files[files_exist_f], function(f) {
    dt <- read_parquet(f)
    setDT(dt)
    dt
  })
  frontier_all <- rbindlist(frontier_list, fill = TRUE)
  rm(frontier_list)
  cat(sprintf("Loaded %s frontier classifications\n",
              format(nrow(frontier_all), big.mark = ",")))
} else {
  frontier_all <- NULL
  cat("No frontier files found - will show interior only\n")
}

# ==============================================================================
# LOAD GRID COORDINATES
# ==============================================================================

cat("Loading grid coordinates...\n")

grid_files <- sapply(1:N_SUB_TILES, function(s) get_grid_filename(s, "parquet"))
files_exist_g <- file.exists(grid_files)
cat(sprintf("Found %d/%d grid files\n", sum(files_exist_g), N_SUB_TILES))

grid_list <- lapply(grid_files[files_exist_g], function(f) {
  dt <- read_parquet(f)
  setDT(dt)
  dt[, .(grid_id, centroid_lon, centroid_lat)]
})
grid_coords <- rbindlist(grid_list, fill = TRUE)
rm(grid_list)

cat(sprintf("Loaded coordinates for %s cells\n",
            format(nrow(grid_coords), big.mark = ",")))

# ==============================================================================
# MERGE AND CLASSIFY
# ==============================================================================

cat("Merging data...\n")

# Merge interior status with coordinates
plot_data <- merge(grid_coords,
                   interior_all[, .(grid_id, is_interior)],
                   by = "grid_id", all.x = TRUE)

# Add frontier status if available
if (!is.null(frontier_all)) {
  plot_data <- merge(plot_data,
                     frontier_all[, .(grid_id, is_frontier)],
                     by = "grid_id", all.x = TRUE)
  plot_data[is.na(is_frontier), is_frontier := FALSE]
} else {
  plot_data[, is_frontier := FALSE]
}

plot_data[is.na(is_interior), is_interior := FALSE]

# Create zone classification
plot_data[, zone := fifelse(is_interior, "Interior",
                            fifelse(is_frontier, "Frontier", "Other"))]

cat("\nZone distribution:\n")
zone_counts <- plot_data[, .N, by = zone]
print(zone_counts)

# ==============================================================================
# SAMPLE FOR PLOTTING (if too many points)
# ==============================================================================

n_points <- nrow(plot_data)
max_points <- 500000  # Limit for reasonable plotting

if (n_points > max_points) {
  cat(sprintf("\nSampling %s points for visualization (from %s total)...\n",
              format(max_points, big.mark = ","),
              format(n_points, big.mark = ",")))

  # Stratified sample - keep more interior/frontier points
  interior_sample <- plot_data[zone == "Interior"][sample(.N, min(.N, max_points * 0.4))]
  frontier_sample <- plot_data[zone == "Frontier"][sample(.N, min(.N, max_points * 0.4))]
  other_sample <- plot_data[zone == "Other"][sample(.N, min(.N, max_points * 0.2))]

  plot_sample <- rbind(interior_sample, frontier_sample, other_sample)
} else {
  plot_sample <- plot_data
}

cat(sprintf("Plotting %s points\n", format(nrow(plot_sample), big.mark = ",")))

# ==============================================================================
# CREATE MAP
# ==============================================================================

cat("Creating map...\n")

# Define colors
zone_colors <- c(
  "Interior" = "#1a9641",   # Dark green - core forest
  "Frontier" = "#fdae61",   # Orange - forest edge
  "Other" = "#d7d7d7"       # Light gray - non-forest/cleared
)

# Create the plot
p <- ggplot() +
  # Plot points (other first, then frontier, then interior on top)
  geom_point(data = plot_sample[zone == "Other"],
             aes(x = centroid_lon, y = centroid_lat, color = zone),
             size = 0.1, alpha = 0.3) +
  geom_point(data = plot_sample[zone == "Frontier"],
             aes(x = centroid_lon, y = centroid_lat, color = zone),
             size = 0.2, alpha = 0.5) +
  geom_point(data = plot_sample[zone == "Interior"],
             aes(x = centroid_lon, y = centroid_lat, color = zone),
             size = 0.2, alpha = 0.7) +

  # Colors
  scale_color_manual(values = zone_colors,
                     name = "Forest Zone",
                     breaks = c("Interior", "Frontier", "Other")) +

  # Limit to tropical regions (TMF coverage)
  coord_cartesian(xlim = c(-180, 180), ylim = c(-30, 30), expand = FALSE) +

  # Theme
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"),
    panel.grid = element_line(color = "gray90"),
    panel.background = element_rect(fill = "aliceblue", color = NA)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +

  # Labels
  labs(
    title = "Global Distribution of Forest Interior and Frontier Zones",
    subtitle = sprintf("Interior: %s cells | Frontier: %s cells | Total: %s cells",
                       format(sum(plot_data$is_interior), big.mark = ","),
                       format(sum(plot_data$is_frontier), big.mark = ","),
                       format(nrow(plot_data), big.mark = ",")),
    x = "Longitude",
    y = "Latitude",
    caption = "Interior = >0% undisturbed forest in all years (1990-2023)\nFrontier = within 100km of interior forest"
  )

# ==============================================================================
# SAVE OUTPUT
# ==============================================================================

output_dir <- file.path(project_root, "output", "validation", "interior_frontier")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_file <- file.path(output_dir, "interior_frontier_map.png")

cat(sprintf("Saving map to: %s\n", output_file))

ggsave(output_file, p, width = 16, height = 8, dpi = 300)

cat("Done!\n")

# Also display summary statistics
cat("\n========================================\n")
cat("SUMMARY STATISTICS\n")
cat("========================================\n")

cat(sprintf("\nTotal grid cells: %s\n", format(nrow(plot_data), big.mark = ",")))
cat(sprintf("Interior cells: %s (%.1f%%)\n",
            format(sum(plot_data$is_interior), big.mark = ","),
            100 * mean(plot_data$is_interior)))
cat(sprintf("Frontier cells: %s (%.1f%%)\n",
            format(sum(plot_data$is_frontier), big.mark = ","),
            100 * mean(plot_data$is_frontier)))
cat(sprintf("Other cells: %s (%.1f%%)\n",
            format(sum(plot_data$zone == "Other"), big.mark = ","),
            100 * mean(plot_data$zone == "Other")))

# Interior statistics
if (sum(plot_data$is_interior) > 0) {
  interior_stats <- interior_all[is_interior == TRUE]
  cat(sprintf("\nInterior cell statistics:\n"))
  cat(sprintf("  Mean undisturbed fraction: %.1f%%\n",
              100 * mean(interior_stats$mean_undisturbed_frac, na.rm = TRUE)))
  cat(sprintf("  Min undisturbed fraction: %.1f%%\n",
              100 * min(interior_stats$min_undisturbed_frac, na.rm = TRUE)))
}

cat("\n========================================\n")
cat(sprintf("Map saved to: %s\n", output_file))
cat("========================================\n")
