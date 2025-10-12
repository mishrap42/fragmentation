if (!require('binsreg')) install.packages("binsreg")
library(here)        # Manage file paths relative to project root
library(sf)          # Spatial data handling (reading shapefiles, geometric operations)
library(ggplot2)     # Data visualization and plotting
library(dplyr)       # Data manipulation (filtering, grouping, summarizing)
library(gridExtra)   # Arrange multiple plots side-by-side
library(readxl)      # Read Excel files
library(binsreg)     # Binned regression plots
library(lubridate)   # Date/time manipulation
library(Hmisc)       # Miscellaneous data analysis functions
library(fixest)      # Fast fixed-effects regressions
library(countrycode) # Convert country codes/names

# This stores your repository path as a function "here()"
if (grepl("ahaanj", here())) {
  dropbox_dir <- "/Users/ahaanj/Dartmouth College Dropbox/Ahaan Jindal/Protected Area Fragmentation/"
} else if (grepl("mishrap", here())) {
  dropbox_dir <- "/Users/mishrap/Dropbox (Personal)/Protected Area Fragmentation/"
} else {
  stop("Username not recognized. Please update the script with your Dropbox path.")
}

intermediate_data_dir = file.path(dropbox_dir, 'Data/build')
figure_dir = file.path(dropbox_dir, 'Figures', '1_calculate_shapes')
if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

# Read in the shapefile of protected areas
pa <- st_read(paste0(dropbox_dir, 'Data/raw/protected sites/protectedsites.shp'))

# Get a stratified sample - 50 areas from each country
set.seed(123)
pa_sample <- pa %>%
  group_by(legalFound, cddaCountr) %>%
  slice_sample(n = 5) %>%
  ungroup()

cat("Total sampled areas:", nrow(pa_sample), "\n")

print('Country counts: validation of sampling')
pa_sample %>% st_drop_geometry() %>% group_by(cddaCountr) %>% summarize(n = n()) %>% print(n = Inf)

print('Foundation year counts:')
pa_sample %>% st_drop_geometry() %>% group_by(legalFound) %>% summarize(n = n()) %>% print(n = Inf)

print('Checking shape validity:')
pa_sample$valid <- st_is_valid(pa_sample)
sprintf('Number of valid geometries: %d out of %d\n', sum( (pa_sample %>% st_drop_geometry())$valid), nrow(pa_sample))

# ========== NINA'S SHAPE METRICS - CORRECTED IMPLEMENTATION ==========
# All require sampling points and computing pairwise distances

# Function to sample interior points from a polygon
sample_interior_points <- function(geom, n_points = 100) {
  # Generate random points inside the polygon
  points <- st_sample(geom, size = n_points, type = "random")
  return(points)
}

# Function to sample perimeter points from a polygon
sample_perimeter_points <- function(geom, n_points = 100) {
  # Cast to linestring and sample along the boundary
  boundary <- st_cast(st_boundary(geom), "LINESTRING")
  points <- st_line_sample(boundary, n = n_points)
  # st_line_sample returns multipoint, convert to individual points
  points <- st_cast(points, "POINT")
  return(points)
}

# Calculate shape metrics for all protected areas
pa_sample <- pa_sample %>%
  mutate(
    area_km2 = as.numeric(units::set_units(st_area(geometry), "km^2")),
    perimeter_km = as.numeric(units::set_units(st_length(st_cast(geometry, "MULTILINESTRING")), "km"))
  ) 
  
pa_centroid = pa_sample %>% st_centroid() %>% st_coordinates() %>% as.data.frame()

point_sample = sample_interior_points(pa_sample$geometry[1], n_points = 100)

# Visualize sampled points and centroid for the first polygon
p <- ggplot() + geom_sf(data = pa_sample[1,], fill = "lightgreen", color = "darkgreen") +
  geom_sf(data = point_sample, color = "red", size = 2) +
  geom_point(data = as.data.frame(pa_centroid[1,]), aes(x = X, y = Y), color = "blue", size = 3) +
  theme_minimal()
ggsave(file.path(figure_dir, "example_point_sampling.png"), p, width = 8, height = 6)

# P. 3 of Nina's appendix:
dist_matrix = units::set_units(st_distance(point_sample, point_sample), "km")
example_disconnection_index = sum(dist_matrix, na.rm = TRUE) / (nrow(dist_matrix) * (nrow(dist_matrix) - 1))

# Process in chunks of 100. TODO: parallelize this.
chunk_size <- 100
n_total <- nrow(pa_sample)
n_chunks <- ceiling(n_total / chunk_size)

pa_sample_disconnection <- list()

for (i in 1:n_chunks) {
  start_idx <- (i - 1) * chunk_size + 1
  end_idx <- min(i * chunk_size, n_total)

  cat(sprintf("Processing chunk %d of %d (rows %d to %d)...\n", i, n_chunks, start_idx, end_idx))

  chunk_result <- pa_sample[start_idx:end_idx, ] %>%
    rowwise() %>%
    mutate(
      # Disconnection: avg distance between interior point pairs
      disconnection_index = {
        interior_pts <- sample_interior_points(geometry, n_points = 100)
        if (length(interior_pts) < 2) {
          NA_real_
        } else {
          dist_matrix <- units::set_units(st_distance(interior_pts, interior_pts), "km")
          sum(dist_matrix, na.rm = TRUE) / (nrow(dist_matrix) * (nrow(dist_matrix) - 1))
        }
      }
    ) %>%
    ungroup()

  pa_sample_disconnection[[i]] <- chunk_result
}

# Combine all chunks
pa_sample_disconnection <- bind_rows(pa_sample_disconnection) %>%
  mutate(disconnection_index = as.numeric(disconnection_index))
st_write(pa_sample_disconnection, paste0(intermediate_data_dir, '/pa_sample_disconnection.shp'), delete_dsn = TRUE)

# Histogram: distribution of disconnection index
p <- pa_sample_disconnection %>%
  st_drop_geometry() %>%
  ggplot() +
  geom_histogram(aes(x = disconnection_index), bins = 30, fill = "lightblue", color = "black") +
  theme_minimal() +
  scale_x_log10() +
  labs(title = "Distribution of Disconnection Index", x = "Disconnection Index", y = "Count")
ggsave(file.path(figure_dir, "disconnection_histogram.png"), p, width = 8, height = 6)

# How does disconnection index vary by legal found year?
sample_summary = pa_sample_disconnection %>%
  filter(legalFound > 1980 & legalFound < 2030) %>%
  st_drop_geometry() %>%
  group_by(legalFound) %>%
  summarise(mean_disconnection = mean(disconnection_index, na.rm = TRUE),
            area_wt_mean_disconnection = wtd.mean(disconnection_index, area_km2, na.rm = TRUE),
            n = n())

p <- sample_summary %>%
  ggplot() +
  geom_col(aes(x = legalFound, y = mean_disconnection), fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Disconnection Index", x = "Legal Found Year", y = "Mean Disconnection Index")
ggsave(file.path(figure_dir, "disconnection_by_year.png"), p, width = 10, height = 6)

# What is the trend in disconnection index over time, residualizing on area?
trend_test = feols(sw(log(disconnection_index), disconnection_rescaled) ~ legalFound + log(area_km2) | csw0(cddaCountr), data = pa_sample_disconnection %>% 
   st_drop_geometry() %>% filter(legalFound > 1980 & legalFound < 2030) %>%
   mutate(disconnection_rescaled = disconnection_index / sd(disconnection_index, na.rm = TRUE))) 

etable(
  trend_test,
  style.tex = style.tex("aer"),
  tex = TRUE,
  dict = c("legalFound" = "Legal Found Year",
  "log(area_km2)" = "Log Area (km^2)",
  "cddaCountr" = "Country",
  "disconnection_rescaled" = "Disconnection Index Std. Dev.",
  "log(disconnection_index)" = "Log Disconnection Index"),
  fitstat = c("n", "r2", "my"),
  digits = 4,
  file = file.path(figure_dir, "trend_test_table.tex")
)

# Country-specific trends: allow legalFound effect to vary by country
country_trends = feols(disconnection_index ~ i(countryname, legalFound) + log(area_km2) | countryname,
                       data = pa_sample_disconnection %>%
                         st_drop_geometry() %>%
                         filter(legalFound > 1980 & legalFound < 2030))

# Extract country-specific trend coefficients
country_trend_coefs <- data.frame(
  countryname = gsub("countryname::", "", gsub(":legalFound", "",
    names(coef(country_trends))[grepl("countryname.*legalFound", names(coef(country_trends)))])),
  trend = coef(country_trends)[grepl("countryname.*legalFound", names(coef(country_trends)))]
)

# Plot country-specific trends
p <- country_trend_coefs %>%
  ggplot() +
  geom_col(aes(x = reorder(countryname, trend), y = trend), fill = "lightblue", color = "black") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Country-Specific Trends in Disconnection Over Time",
       x = "Country",
       y = "Trend Coefficient (Change in Disconnection per Year)")
ggsave(file.path(figure_dir, "country_trends.png"), p, width = 10, height = 8)

# Which countries have the most disconnected protected areas on average?
pa_sample_disconnection$countryname = countrycode(pa_sample_disconnection$cddaCountr, origin = 'eurostat', destination = 'country.name')
country_summary = pa_sample_disconnection %>%
  st_drop_geometry() %>%
  group_by(countryname) %>%
  summarise(mean_disconnection = mean(disconnection_index, na.rm = TRUE),
            area_wt_mean_disconnection = wtd.mean(disconnection_index, area_km2, na.rm = TRUE),
            n = n()) %>%
  arrange(desc(mean_disconnection))

p <- country_summary %>%
  ggplot() +
  geom_col(aes(x = reorder(countryname, mean_disconnection), y = mean_disconnection), fill = "lightblue", color = "black") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Mean Disconnection Index by Country", x = "Country", y = "Mean Disconnection Index")
ggsave(file.path(figure_dir, "disconnection_by_country.png"), p, width = 10, height = 8)


area_residualized_country_ranking = feols(disconnection_index ~ i(countryname) + log(area_km2) - 1, data = pa_sample_disconnection %>%
   st_drop_geometry() %>% filter(legalFound > 1980 & legalFound < 2030) )

# Extract country coefficients (residualized disconnection)
country_coefs <- data.frame(
  countryname = gsub("countryname::", "", names(coef(area_residualized_country_ranking))[grepl("countryname", names(coef(area_residualized_country_ranking)))]),
  residualized_disconnection = coef(area_residualized_country_ranking)[grepl("countryname", names(coef(area_residualized_country_ranking)))]
)

# Plot residualized disconnection by country
p <- country_coefs %>%
  ggplot() +
  geom_col(aes(x = reorder(countryname, residualized_disconnection), y = residualized_disconnection), fill = "lightblue", color = "black") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Area-Residualized Disconnection Index by Country", x = "Country", y = "Residualized Disconnection Index")
ggsave(file.path(figure_dir, "residualized_disconnection_by_country.png"), p, width = 10, height = 8)

# Decomposition: merge mean disconnection with residualized disconnection
decomposition_data <- country_summary %>%
  left_join(country_coefs, by = "countryname") %>%
  mutate(
    area_effect = mean_disconnection - residualized_disconnection,
    true_disconnection = residualized_disconnection
  ) %>%
  select(countryname, mean_disconnection, true_disconnection, area_effect)

decomposition_data <- rbind(
  decomposition_data %>% mutate(component = "true_disconnection", value = true_disconnection),
  decomposition_data %>% mutate(component = "area_effect", value = area_effect)
) %>% select(countryname, mean_disconnection, component, value)

# Stacked bar chart showing decomposition
p <- decomposition_data %>%
  ggplot(aes(x = reorder(countryname, mean_disconnection), y = value, fill = component)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("true_disconnection" = "steelblue", "area_effect" = "coral"),
                    labels = c("Area Effect", "Residualized Disconnection")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Disconnection Index Decomposition by Country",
       x = "Country",
       y = "Disconnection Index",
       fill = "Component")
ggsave(file.path(figure_dir, "disconnection_decomposition_by_country.png"), p, width = 10, height = 8) 

#   %>%
#   rowwise() %>%
#   mutate(
#     centroid = list(st_centroid(geometry)),
#
#     # Disconnection: avg distance between interior point pairs
#     disconnection_index = {
#       interior_pts <- sample_interior_points(geometry, n_points = 100)
#       if (length(interior_pts) < 2) {
#         NA_real_
#       } else {
#         dist_matrix <- st_distance(interior_pts, interior_pts)
#         dist_values <- as.numeric(dist_matrix[upper.tri(dist_matrix)])
#         round(mean(dist_values, na.rm = TRUE) / 1000, 3)
#       }
#     },
#
#     # Remoteness: avg distance from interior points to center
#     remoteness = {
#       interior_pts <- sample_interior_points(geometry, n_points = 100)
#       if (length(interior_pts) < 1) {
#         NA_real_
#       } else {
#         distances <- st_distance(interior_pts, centroid[[1]])
#         round(mean(as.numeric(distances), na.rm = TRUE) / 1000, 3)
#       }
#     },
#
#     # Spin: avg distance from boundary points to center
#     spin = {
#       perimeter_pts <- sample_perimeter_points(geometry, n_points = 100)
#       if (length(perimeter_pts) < 1) {
#         NA_real_
#       } else {
#         distances <- st_distance(perimeter_pts, centroid[[1]])
#         round(mean(as.numeric(distances), na.rm = TRUE) / 1000, 3)
#       }
#     },
#
#     # Range: avg distance between boundary point pairs
#     range_km = {
#       perimeter_pts <- sample_perimeter_points(geometry, n_points = 100)
#       if (length(perimeter_pts) < 2) {
#         NA_real_
#       } else {
#         dist_matrix <- st_distance(perimeter_pts, perimeter_pts)
#         dist_values <- as.numeric(dist_matrix[upper.tri(dist_matrix)])
#         round(mean(dist_values, na.rm = TRUE) / 1000, 3)
#       }
#     }
#   ) %>%
#   ungroup() %>%
#   select(-centroid)

# print(pa_shapes$disconnection_index[1], digits = 10)

# # Summary by country (medians)
# shape_summary <- pa_shapes %>%
#   st_drop_geometry() %>%
#   group_by(cddaCountr) %>%
#   summarize(
#     n_areas = n(),
#     med_area = median(area_km2, na.rm = TRUE),
#     med_disconnection = median(disconnection_index, na.rm = TRUE),
#     med_remoteness = median(remoteness, na.rm = TRUE),
#     med_spin = median(spin, na.rm = TRUE),
#     med_range = median(range_km, na.rm = TRUE),
#     .groups = 'drop'
#   ) %>%
#   arrange(desc(med_area))

# View(shape_summary)

# pa_shapes$disconnection_index[1:10]


# # ========== FIND EXTREMES FOR EACH METRIC ==========

# # 1. DISCONNECTION INDEX (higher = more spread out internally)
# most_disconnected <- pa_shapes %>%
#   st_drop_geometry() %>%
#   arrange(desc(disconnection_index)) %>%
#   select(cddaId, siteName, cddaCountr, disconnection_index, area_km2) %>%
#   head(2)

# least_disconnected <- pa_shapes %>%
#   st_drop_geometry() %>%
#   arrange(disconnection_index) %>%
#   select(cddaId, siteName, cddaCountr, disconnection_index, area_km2) %>%
#   head(2)

# # 2. REMOTENESS (higher = farther from center on average)
# most_remote <- pa_shapes %>%
#   st_drop_geometry() %>%
#   arrange(desc(remoteness)) %>%
#   select(cddaId, siteName, cddaCountr, remoteness, area_km2) %>%
#   head(2)

# least_remote <- pa_shapes %>%
#   st_drop_geometry() %>%
#   arrange(remoteness) %>%
#   select(cddaId, siteName, cddaCountr, remoteness, area_km2) %>%
#   head(2)

# # 3. SPIN (higher = perimeter is farther from center)
# highest_spin <- pa_shapes %>%
#   st_drop_geometry() %>%
#   arrange(desc(spin)) %>%
#   select(cddaId, siteName, cddaCountr, spin, area_km2) %>%
#   head(2)

# lowest_spin <- pa_shapes %>%
#   st_drop_geometry() %>%
#   arrange(spin) %>%
#   select(cddaId, siteName, cddaCountr, spin, area_km2) %>%
#   head(2)

# # 4. RANGE (higher = larger perimeter extent)
# largest_range <- pa_shapes %>%
#   st_drop_geometry() %>%
#   arrange(desc(range_km)) %>%
#   select(cddaId, siteName, cddaCountr, range_km, area_km2) %>%
#   head(2)

# smallest_range <- pa_shapes %>%
#   st_drop_geometry() %>%
#   arrange(range_km) %>%
#   select(cddaId, siteName, cddaCountr, range_km, area_km2) %>%
#   head(2)

# # PRINT ALL RESULTS
# cat("\n========== DISCONNECTION INDEX ==========\n")
# cat("Economic interpretation: Average internal travel distance\n")
# cat("\nMOST DISCONNECTED (largest internal distances):\n")
# print(most_disconnected)
# cat("\nLEAST DISCONNECTED (smallest internal distances):\n")
# print(least_disconnected)

# cat("\n========== REMOTENESS ==========\n")
# cat("Economic interpretation: Average distance to center\n")
# cat("\nMOST REMOTE (farthest from center):\n")
# print(most_remote)
# cat("\nLEAST REMOTE (closest to center):\n")
# print(least_remote)

# cat("\n========== SPIN ==========\n")
# cat("Economic interpretation: Average boundary-to-center distance\n")
# cat("\nHIGHEST SPIN (boundary far from center):\n")
# print(highest_spin)
# cat("\nLOWEST SPIN (boundary close to center):\n")
# print(lowest_spin)

# cat("\n========== RANGE ==========\n")
# cat("Economic interpretation: Average boundary extent\n")
# cat("\nLARGEST RANGE (biggest boundary extent):\n")
# print(largest_range)
# cat("\nSMALLEST RANGE (smallest boundary extent):\n")
# print(smallest_range)

# # ========== VISUALIZE EXAMPLE ==========
# specific_id <- most_disconnected$cddaId[1]
# pa_single <- pa_shapes %>% filter(cddaId == specific_id)

# # First plot
# p1 <- ggplot(pa_single) +
#   geom_sf(fill = "lightgreen", color = "darkgreen", size = 1) +
#   theme_minimal() +
#   labs(
#     title = paste("Example:", pa_single$siteName),
#     subtitle = paste("Country:", pa_single$cddaCountr)
#   )

# # Example: another polygon to compare
# specific_id2 <- least_disconnected$cddaId[2]
# pa_single2 <- pa_shapes %>% filter(cddaId == specific_id2)

# # Second plot
# p2 <- ggplot(pa_single2) +
#   geom_sf(fill = "lightblue", color = "navy", size = 1) +
#   theme_minimal() +
#   labs(
#     title = paste("Example:", pa_single2$siteName),
#     subtitle = paste("Country:", pa_single2$cddaCountr)
#   )

# # Display both side by side
# grid.arrange(p1, p2, ncol = 2)


# cat("\n=== EXAMPLE AREA METRICS ===\n")
# cat("Site Name:", pa_single$siteName, "\n")
# cat("Country:", pa_single$cddaCountr, "\n")
# cat("Area:", round(pa_single$area_km2, 2), "km²\n")
# cat("Disconnection Index:", round(pa_single$disconnection_index, 3), "km (avg internal distance)\n")
# cat("Remoteness:", round(pa_single$remoteness, 3), "km (avg distance to center)\n")
# cat("Spin:", round(pa_single$spin, 3), "km (avg boundary-to-center distance)\n")
# cat("Range:", round(pa_single$range_km, 3), "km (avg boundary extent)\n")


# # ------------ PLOT INDEX AVERAGES BY TIME 

# all_data <- pa_shapes %>%
#   st_drop_geometry() %>%
#   mutate(legal_year = as.numeric(substr(as.character(legalFound), 1, 4))) %>%
#   filter(!is.na(legal_year))

# by_year <- all_data %>%
#   group_by(legal_year) %>%
#   summarize(
#     range = median(range_km, na.rm = TRUE),
#     disconnection = median(disconnection_index, na.rm = TRUE),
#     spin = median(spin, na.rm = TRUE),
#     remoteness = median(remoteness, na.rm = TRUE),
#     .groups = 'drop'
#   ) %>%
#   arrange(legal_year)

# p1 <- ggplot(by_year, aes(x = legal_year, y = range)) +
#   geom_point(color = "#4DAF4A", size = 3) +
#   geom_line(color = "#4DAF4A", linewidth = 1) +
#   scale_x_continuous(limits = c(1975, 2025)) +
#   scale_y_continuous(limits = c(0, 15)) +
#   labs(title = "Range", x = "Time", y = "Median (km)") +
#   theme_minimal()

# p2 <- ggplot(by_year, aes(x = legal_year, y = disconnection)) +
#   geom_point(color = "#E41A1C", size = 3) +
#   geom_line(color = "#E41A1C", linewidth = 1) +
#   scale_x_continuous(limits = c(1975, 2025)) +
#   scale_y_continuous(limits = c(0.8, 1.2)) +
#   labs(title = "Disconnection", x = "Time", y = "Median (km)") +
#   theme_minimal()

# p3 <- ggplot(by_year, aes(x = legal_year, y = spin)) +
#   geom_point(color = "#377EB8", size = 3) +
#   geom_line(color = "#377EB8", linewidth = 1) +
#   scale_x_continuous(limits = c(1975, 2025)) +
#   scale_y_continuous(limits = c(0.15, 0.6)) +
#   labs(title = "Spin", x = "Time", y = "Median (km)") +
#   theme_minimal()

# p4 <- ggplot(by_year, aes(x = legal_year, y = remoteness)) +
#   geom_point(color = "#984EA3", size = 3) +
#   geom_line(color = "#984EA3", linewidth = 1) +
#   scale_x_continuous(limits = c(1975, 2025)) +
#   scale_y_continuous(limits = c(0, 1)) +
#   labs(title = "Remoteness", x = "Time", y = "Median (km)") +
#   theme_minimal()

# gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)

# # -----------

# # Prepare data
# binned_data <- pa_shapes %>%
#   st_drop_geometry() %>%
#   filter(!is.na(disconnection_index) & !is.na(spin))

# # Simple regression plot
# ggplot(binned_data, aes(x = disconnection_index, y = spin)) +
#   geom_point(alpha = 0.4, color = "#377EB8", size = 2) +
#   geom_smooth(method = "lm", color = "#E41A1C", se = TRUE, linewidth = 1.5) +
#   scale_x_continuous(limits = c(0, 2)) +
#   scale_y_continuous(limits = c(0.0, 1.0)) +
#   labs(
#     title = "Spin vs Disconnection Index",
#     x = "Disconnection Index (km)",
#     y = "Spin Index (km)"
#   ) +
#   theme_minimal()

