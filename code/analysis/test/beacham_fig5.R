# ==============================================================================
# REPLICATION: Beacham (2022), "Conserving What's Left: The Political Economy
# of Protected Area Location", Figure 5
# ==============================================================================
# Figure 5 is a map of the Uatuma-Trombetas Moist Forests ecoregion (Brazil)
# showing the growth of the protected-area network over time. Each PA is shaded
# by the earliest decadal snapshot in which it was already part of the network
# (1980, 1990, 2000, 2010, 2020), so light green = long-standing protection and
# dark green = protection added in the 2010s. Equivalent to drawing the
# cumulative 2020 network, then the 2010 network on top of it, and so on.
#
# Sources (as in the paper): WDPA (PA polygons + STATUS_YR), Ecoregions 2017
# (Dinerstein et al. 2017), GADM/Natural Earth for the Brazil inset. EPSG:4326.
#
# One-off script - not part of the pipeline.
#
# Usage: Rscript code/analysis/test/beacham_fig5.R
# ==============================================================================

# Deliberately does not source BUILD_workspace.R: this needs none of the
# pipeline's heavy dependencies, and staying self-contained keeps it runnable
# from a plain R install.
here::i_am('code/analysis/test/beacham_fig5.R')
project_root <- here::here()
log_message <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(cowplot)
  library(rnaturalearth)
})

sf_use_s2(FALSE)  # planar ops on lon/lat; fine at this scale and avoids
                  # s2 edge failures on WDPA's self-touching rings

ECO_NAME  <- "Uatumã-Trombetas moist forests"
SNAPSHOTS <- c(1980, 1990, 2000, 2010, 2020)

# Greens, light (oldest protection) -> dark (newest), as in the paper
PAL <- c("1980" = "#EDF8E9", "1990" = "#BAE4B3", "2000" = "#74C476",
         "2010" = "#238B45", "2020" = "#00441B")

# s2 is off, so st_area on lon/lat would need lwgeom. Measure in a local
# equal-area projection instead.
EQAREA <- "+proj=laea +lat_0=0 +lon_0=-57 +datum=WGS84 +units=m +no_defs"

# The lab mount holds the built data when running off-cluster; fall back to it
# whenever the project-local copy is absent.
MOUNT <- "/Volumes/MishraLab/Fragmentation"
resolve <- function(...) {
  local <- file.path(project_root, ...)
  if (file.exists(local)) return(local)
  file.path(MOUNT, ...)
}

out_dir <- file.path(project_root, "output", "analysis", "beacham_fig5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# ECOREGION
# ==============================================================================

eco_path <- resolve("Data", "raw", "spatial", "ecoregions2017",
                    "Ecoregions2017.shp")
if (!file.exists(eco_path)) {
  stop("Ecoregions 2017 not found at ", eco_path,
       "\nDownload: https://storage.googleapis.com/teow2016/Ecoregions2017.zip")
}

log_message("Loading ecoregion polygon...")
eco <- st_read(
  eco_path,
  query = sprintf("SELECT ECO_NAME, ECO_ID FROM Ecoregions2017 WHERE ECO_NAME = '%s'",
                  ECO_NAME),
  quiet = TRUE
)
stopifnot(nrow(eco) == 1)
eco <- st_make_valid(st_transform(eco, 4326))
log_message(sprintf("  %s (ECO_ID %s)", eco$ECO_NAME, eco$ECO_ID))

# ==============================================================================
# PROTECTED AREAS
# ==============================================================================
# Spatial-filter the read on the ecoregion bbox so we never pull the full
# 255k-feature WDPA layer off the network mount.

# The read + clip takes ~1 min against the lab mount; cache it so the plotting
# can be iterated on cheaply. Delete the file (or set REFRESH) to rebuild.
REFRESH   <- FALSE
cache_dir <- file.path(project_root, "Data", "temp")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
pa_cache  <- file.path(cache_dir, "beacham_fig5_pa.rds")

if (!REFRESH && file.exists(pa_cache)) {

log_message("Loading clipped PAs from cache...")
pa <- readRDS(pa_cache)

} else {

wdpa_file  <- resolve("Data", "build", "wdpa_cleaned.gpkg")
wdpa_layer <- st_layers(wdpa_file)$name[1]

# Read the layer directly rather than through a SQL query: only a plain layer
# read lets the GPKG driver use its RTree index, and a full scan of the 255k
# features across the lab mount takes tens of minutes.
log_message("Reading WDPA within the ecoregion bounding box...")
pa <- st_read(
  wdpa_file, layer = wdpa_layer,
  wkt_filter = st_as_text(st_as_sfc(st_bbox(eco))),
  quiet = TRUE
)
pa <- pa[grepl("BRA", pa$ISO3), ]
log_message(sprintf("  %d Brazilian PAs in the bounding box", nrow(pa)))

# Designated terrestrial PAs in force by the end of 2020, as in the paper.
# STATUS_YR == 0 is "not reported" and cannot be placed on the timeline.
n_noyear <- sum(pa$STATUS_YR == 0)
pa <- pa[pa$PA_DEF == 1 &
         pa$STATUS %in% c("Designated", "Inscribed", "Established") &
         pa$MARINE != "2" &
         pa$STATUS_YR > 0 & pa$STATUS_YR <= 2020, ]
log_message(sprintf("  %d after filtering (%d dropped for unreported STATUS_YR)",
                    nrow(pa), n_noyear))

log_message("Clipping to the ecoregion...")
pa <- st_make_valid(st_transform(pa, 4326))
pa <- suppressWarnings(st_intersection(pa, st_geometry(eco)))
pa <- pa[!st_is_empty(pa), ]
pa <- pa[st_dimension(pa) == 2, ]  # drop line/point slivers from the clip
saveRDS(pa, pa_cache)

}

log_message(sprintf("  %d PAs intersect the ecoregion, earliest %d",
                    nrow(pa), min(pa$STATUS_YR)))

# ==============================================================================
# CUMULATIVE NETWORKS
# ==============================================================================
# One copy of the network per snapshot year, stacked darkest-first so that the
# lighter (older) networks are painted over the darker (later) ones. The visible
# colour of any polygon is therefore its earliest snapshot.

net <- do.call(rbind, lapply(rev(SNAPSHOTS), function(y) {
  s <- pa[pa$STATUS_YR <= y, ]
  s$snapshot <- factor(y, levels = SNAPSHOTS)
  s
}))

log_message("Cumulative PA share of the ecoregion:")
eco_ea   <- st_transform(eco, EQAREA)
pa_ea    <- st_transform(pa, EQAREA)
eco_area <- as.numeric(sum(st_area(eco_ea)))
for (y in SNAPSHOTS) {
  sel <- pa_ea[pa_ea$STATUS_YR <= y, ]
  cov <- if (nrow(sel) == 0) 0 else
    as.numeric(st_area(st_union(st_make_valid(sel)))) / eco_area
  log_message(sprintf("  %d: %5.1f%% of the ecoregion (%d PAs)",
                      y, 100 * cov, nrow(sel)))
}

# ==============================================================================
# MAIN MAP
# ==============================================================================

bb <- st_bbox(eco)

main <- ggplot() +
  geom_sf(data = eco, fill = "white", colour = "grey25", linewidth = 0.3) +
  geom_sf(data = net, aes(fill = snapshot), colour = "grey15", linewidth = 0.12) +
  scale_fill_manual(values = PAL, name = "PA Network in:", drop = FALSE) +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
           ylim = c(bb["ymin"], bb["ymax"]), expand = TRUE) +
  labs(title = "Growth in Protected Area Coverage Over Time",
       subtitle = "Uatumã-Trombetas Moist Forests Ecoregion, Brazil") +
  theme_void(base_family = "serif") +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 17),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position      = c(0.03, 0.20),
    legend.justification = c(0, 0.5),
    legend.title  = element_text(size = 10),
    legend.text   = element_text(size = 9),
    legend.key.size = unit(0.42, "cm"),
    # the ecoregion reaches under the legend here, unlike in the original
    legend.background = element_rect(fill = "white", colour = NA),
    legend.margin = margin(4, 6, 4, 4),
    plot.margin   = margin(8, 8, 8, 8)
  )

# ==============================================================================
# BRAZIL INSET
# ==============================================================================

brazil <- st_transform(
  ne_countries(country = "Brazil", scale = "medium", returnclass = "sf"), 4326)
eco_box <- st_as_sfc(st_bbox(eco), crs = 4326)

inset <- ggplot() +
  geom_sf(data = brazil, fill = "grey96", colour = "grey30", linewidth = 0.3) +
  geom_sf(data = eco_box, fill = "red", alpha = 0.25,
          colour = "red", linewidth = 0.3) +
  geom_sf(data = eco, fill = "red", colour = "red", linewidth = 0.1) +
  theme_void() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        plot.background = element_rect(fill = "white", colour = NA))

# ==============================================================================
# ASSEMBLE
# ==============================================================================

credit <- paste("Replication of Beacham (2022), Figure 5",
                "Data: WDPA, Ecoregions 2017, Natural Earth",
                "Projection: EPSG:4326",
                sprintf("Date created: %s", format(Sys.Date(), "%d %B %Y")),
                sep = "\n")

fig <- ggdraw(main) +
  draw_plot(inset, x = 0.68, y = 0.06, width = 0.28, height = 0.30) +
  draw_label(credit, x = 0.035, y = 0.055, hjust = 0, vjust = 0,
             size = 6, fontfamily = "serif", colour = "grey20") +
  theme(plot.background = element_rect(fill = "white", colour = "black"))

out_file <- file.path(out_dir, "beacham_fig5.pdf")
ggsave(out_file, fig, width = 9, height = 6.2, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", out_file), fig, width = 9, height = 6.2, dpi = 300)
log_message(sprintf("Wrote %s", out_file))
