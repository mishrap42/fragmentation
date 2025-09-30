library(here)
library(sf)
library(ggplot2)
library(dplyr)
library(gridExtra)
# This stores your repository path as a function "here()"
here::i_am('code/build/1_calculate_shapes.R')
# Depending on your computer, you may need to change the path to your Dropbox folder
# For example, on Prakash's computer:
if(grep('mishrap', here()) == 1) {
  dropbox_dir = '/Users/mishrap/Dropbox (Personal)/Protected Area Fragmentation/'
} else if(grep('ahaan', here()) == 1){
    # fill in the blank
    dropbox_dir = '...'
}

# Read in the shapefile of protected areas
pa <- st_read(paste0(dropbox_dir, 'Data/protected sites/protectedsites.shp'))

# Basic data checks
cat("Total protected areas:", nrow(pa), "\n")
cat("Invalid geometries:", sum(!st_is_valid(pa)), "\n")
cat("Countries represented:", length(unique(pa$cddaCountr)), "\n")
cat("IUCN categories:", paste(sort(unique(pa$iucnCatego)), collapse = ", "),
    "\n")
cat("Year range:", min(pa$legalFound, na.rm = TRUE), "-",
    max(pa$legalFound, na.rm = TRUE), "\n")
cat("Duplicate IDs:", sum(duplicated(pa$cddaId)), "\n")

# Visualization: Timeline of protected area establishment
pa_timeline <- pa %>%
  mutate(year = as.numeric(legalFound),
         area_km2 = as.numeric(st_area(geometry)) / 1e6) %>%
  st_set_geometry(NULL) %>%
  filter(!is.na(year))

# Aggregate by year for area plot
pa_area_by_year <- pa_timeline %>%
  st_drop_geometry(NULL) %>%
  group_by(year) %>%
  summarize(total_area = sum(area_km2, na.rm = TRUE))

# Panel A: Number of protected areas
p1 <- ggplot(pa_timeline, aes(x = year)) +
  geom_histogram(binwidth = 1, fill = "gray30", color = "white", size = 0.2) +
  scale_x_continuous(breaks = seq(1900, 2020, 20), limits = c(1900, 2020)) +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  labs(x = "Year of Establishment",
       y = "Number of Protected Areas",
       title = "Panel A: Number of Protected Areas Established") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold", size = 11, hjust = 0),
        axis.title = element_text(face = "plain"),
        axis.text = element_text(color = "black"),
        plot.margin = margin(10, 10, 10, 10))

# Panel B: Total area protected
p2 <- ggplot(pa_area_by_year, aes(x = year, y = total_area)) +
  geom_bar(stat = "identity", fill = "gray30", color = "white", size = 0.2) +
  scale_x_continuous(breaks = seq(1900, 2020, 20), limits = c(1900, 2020)) +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  labs(x = "Year of Establishment",
       y = expression(paste("Total Area Protected (", km^2, ")")),
       title = "Panel B: Total Area of Protected Areas Established") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold", size = 11, hjust = 0),
        axis.title = element_text(face = "plain"),
        axis.text = element_text(color = "black"),
        plot.margin = margin(10, 10, 10, 10))

# Combine panels
grid.arrange(p1, p2, ncol = 1)
