# ==============================================================================
# VISUALIZATION: Protected Area + TMF Undisturbed Forest Share
# ==============================================================================
# Creates a visualization of one randomly selected protected area
# overlaid with TMF undisturbed forest share for grid cells within it,
# faceted by year.
#
# Note: WDPA data is now baked into grid cells (each cell is fully
# protected or not, since grid is cut on WDPA boundaries).
#
# Usage: Rscript code/analysis/visualize_wdpa_tmf.R [seed]
#        Optional seed argument for reproducibility
# ==============================================================================

here::i_am('code/analysis/visualize_wdpa_tmf.R')
source("code/build/BUILD_workspace.R")

library(ggplot2)
library(viridis)
library(cowplot)

# Set seed for reproducibility (optional argument)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  set.seed(as.integer(args[1]))
} else {
  set.seed(42)
}

log_message("=== WDPA + TMF Visualization ===")

# ==============================================================================
# LOAD DATA
# ==============================================================================

log_message("Loading final panel data...")
panel_file <- file.path(final_output_path, "TMF_5km_panel.parquet")
if (!file.exists(panel_file)) {
  stop("Final panel file not found. Run the full pipeline first.")
}
panel_data <- read_parquet(panel_file)
setDT(panel_data)

log_message(sprintf("Panel data: %s rows, years %d-%d",
                    format(nrow(panel_data), big.mark = ","),
                    min(panel_data$year, na.rm = TRUE),
                    max(panel_data$year, na.rm = TRUE)))

# Get protected cells from panel data (WDPA is baked into grid)
log_message("Extracting protected cells from panel data...")

# Get unique grid cells that are protected
protected_cells <- unique(panel_data[is_protected == TRUE, .(
  grid_id, tile_id, wdpa_pid, iucn_cat, desig_year, gov_type,
  country_iso3, centroid_lon, centroid_lat
)])

log_message(sprintf("Protected cells: %s",
                    format(nrow(protected_cells), big.mark = ",")))

# ==============================================================================
# SELECT A RANDOM PROTECTED AREA
# ==============================================================================

# Get unique protected area IDs
unique_pa_ids <- unique(protected_cells[!is.na(wdpa_pid)]$wdpa_pid)
log_message(sprintf("Unique protected areas: %d", length(unique_pa_ids)))

if (length(unique_pa_ids) == 0) {
  stop("No protected areas found in the data.")
}

# Select random PA
selected_pa <- sample(unique_pa_ids, 1)
log_message(sprintf("Selected PA ID: %d", selected_pa))

# Get cells in this PA
pa_cells <- protected_cells[wdpa_pid == selected_pa]
log_message(sprintf("Cells in selected PA: %d", nrow(pa_cells)))

# Get PA metadata
pa_info <- pa_cells[1]
log_message(sprintf("PA Name: %s", ifelse(!is.na(pa_info$iucn_cat),
                                          pa_info$iucn_cat, "Unknown IUCN")))

# ==============================================================================
# LOAD GRID GEOMETRY FOR SELECTED CELLS
# ==============================================================================

log_message("Loading grid geometry for selected cells...")

# tile_id is already in the data
tile_ids <- unique(pa_cells$tile_id)
log_message(sprintf("Tiles containing PA cells: %s", paste(tile_ids, collapse = ", ")))

# Load grid files for these tiles
grid_list <- list()
for (tid in tile_ids) {
  # Get all sub-tile files for this TMF tile
  grid_files <- get_grid_files_for_tmf_tile(tid, "gpkg")
  grid_files <- grid_files[file.exists(grid_files)]

  if (length(grid_files) > 0) {
    for (f in grid_files) {
      grid_sf <- st_read(f, quiet = TRUE)
      # Filter to cells in our PA
      grid_sf <- grid_sf[grid_sf$grid_id %in% pa_cells$grid_id, ]
      if (nrow(grid_sf) > 0) {
        grid_list[[length(grid_list) + 1]] <- grid_sf
      }
    }
  }
}

if (length(grid_list) == 0) {
  stop("Could not load grid geometries for selected PA cells.")
}

grid_sf <- do.call(rbind, grid_list)
log_message(sprintf("Loaded %d grid cell geometries", nrow(grid_sf)))

# ==============================================================================
# EXTRACT TMF DATA FOR PA CELLS
# ==============================================================================

log_message("Extracting TMF data for PA cells...")

# Get panel data for cells in this PA
tmf_pa <- panel_data[grid_id %in% pa_cells$grid_id]
log_message(sprintf("TMF observations for PA cells: %s", format(nrow(tmf_pa), big.mark = ",")))

# Undisturbed fraction is Undisturbed_TMF column from wide format
# Check column name (might be Undisturbed_TMF or similar)
if ("Undisturbed_TMF" %in% names(tmf_pa)) {
  tmf_pa[, undisturbed_frac := Undisturbed_TMF]
} else if ("frac_1" %in% names(tmf_pa)) {
  tmf_pa[, undisturbed_frac := frac_1]
} else {
  # List available columns for debugging
  log_message(sprintf("Available columns: %s", paste(names(tmf_pa), collapse = ", ")))
  stop("Could not find undisturbed forest column in data.")
}

# ==============================================================================
# LOAD WDPA BOUNDARY FOR OVERLAY
# ==============================================================================

log_message("Loading WDPA boundary for overlay...")

# Read the full WDPA and filter to selected PA
wdpa_sf <- read_wdpa()
pa_boundary <- wdpa_sf[wdpa_sf$WDPAID == selected_pa | wdpa_sf$WDPA_PID == selected_pa, ]

if (nrow(pa_boundary) == 0) {
  # Try alternative column names
  if ("wdpa_pid" %in% names(wdpa_sf)) {
    pa_boundary <- wdpa_sf[wdpa_sf$wdpa_pid == selected_pa, ]
  }
}

has_boundary <- nrow(pa_boundary) > 0
if (has_boundary) {
  log_message(sprintf("Loaded PA boundary with %d features", nrow(pa_boundary)))
} else {
  log_message("Warning: Could not load PA boundary polygon. Using cell extent.")
}

# ==============================================================================
# CREATE VISUALIZATION
# ==============================================================================

log_message("Creating visualization...")

# Get available years
years <- sort(unique(tmf_pa$year))
log_message(sprintf("Years available: %d-%d", min(years), max(years)))

# Select subset of years for visualization (every 5 years + first and last)
year_subset <- c(min(years), seq(1995, 2020, by = 5), max(years))
year_subset <- sort(unique(year_subset[year_subset %in% years]))
log_message(sprintf("Years for visualization: %s", paste(year_subset, collapse = ", ")))

# Filter to selected years
tmf_pa_subset <- tmf_pa[year %in% year_subset]

# Create plot list
plot_list <- list()

for (yr in year_subset) {
  # Get data for this year
  yr_data <- tmf_pa_subset[year == yr]

  # Merge with grid geometry
  grid_yr <- merge(grid_sf, yr_data, by = "grid_id", all.x = TRUE)

  # Create plot
  p <- ggplot() +
    geom_sf(data = grid_yr, aes(fill = undisturbed_frac), color = NA) +
    scale_fill_viridis(
      name = "Undisturbed\nForest",
      limits = c(0, 1),
      labels = scales::percent,
      na.value = "grey80"
    ) +
    labs(title = as.character(yr)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )

  # Add PA boundary if available
  if (has_boundary) {
    p <- p + geom_sf(data = pa_boundary, fill = NA, color = "red", linewidth = 0.5)
  }

  plot_list[[as.character(yr)]] <- p
}

# Create combined plot with shared legend
combined <- plot_grid(
  plotlist = plot_list,
  ncol = min(4, length(plot_list)),
  align = "hv"
)

# Add title
pa_name <- ifelse(!is.na(pa_info$iucn_cat),
                  sprintf("Protected Area (IUCN: %s)", pa_info$iucn_cat),
                  sprintf("Protected Area (ID: %d)", selected_pa))

title <- ggdraw() +
  draw_label(
    sprintf("%s\nUndisturbed Forest Share by Year", pa_name),
    fontface = 'bold',
    size = 14
  )

# Create legend
legend_plot <- ggplot(data.frame(x = 1, y = 1, fill = 0.5)) +
  geom_tile(aes(x, y, fill = fill)) +
  scale_fill_viridis(
    name = "Undisturbed\nForest Share",
    limits = c(0, 1),
    labels = scales::percent
  ) +
  theme_void() +
  theme(legend.position = "bottom")

legend <- get_legend(legend_plot)

# Combine all
final_plot <- plot_grid(
  title,
  combined,
  legend,
  ncol = 1,
  rel_heights = c(0.1, 0.8, 0.1)
)

# ==============================================================================
# SAVE OUTPUT
# ==============================================================================

output_dir <- file.path(project_root, "output", "validation", "wdpa_tmf")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_file <- file.path(output_dir, sprintf("wdpa_tmf_pa_%d.png", selected_pa))

ggsave(
  output_file,
  plot = final_plot,
  width = 12,
  height = 8,
  dpi = 150
)

log_message(sprintf("Saved visualization: %s", output_file))

# Also save as PDF for publication quality
pdf_file <- file.path(output_dir, sprintf("wdpa_tmf_pa_%d.pdf", selected_pa))
ggsave(pdf_file, plot = final_plot, width = 12, height = 8)
log_message(sprintf("Saved PDF: %s", pdf_file))

log_message("=== Visualization complete ===")
