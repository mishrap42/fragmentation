# ==============================================================================
# DIAGNOSE WDPA CLEANING: Tabulate areas dropped at each cleaning step
# ==============================================================================
# Quantifies the impact of geometry cleaning steps on WDPA coverage:
#   1. BBOX filtering (oversized/corrupt geometries)
#   2. st_simplify (degenerate vertices)
#   3. GEOMETRYCOLLECTION extraction (polygon components only)
#
# Run: Rscript code/build/diagnose_wdpa_cleaning.R
# ==============================================================================

.libPaths(c('/dartfs-hpc/rc/lab/M/MishraP/Rlibs/4.4.0/', .libPaths()))

library(here)
here::i_am("code/build/diagnose_wdpa_cleaning.R")
source("code/build/BUILD_workspace.R")

cat(rep("=", 71), "\n", sep = "")
cat("WDPA CLEANING IMPACT ANALYSIS\n")
cat(rep("=", 71), "\n", sep = "")
cat(sprintf("Date: %s\n\n", Sys.time()))

# ==============================================================================
# LOAD RAW WDPA
# ==============================================================================

cat("Loading cleaned WDPA data...\n")
if (!file.exists(wdpa_clean_path)) {
  stop(sprintf("WDPA file not found: %s", wdpa_clean_path))
}

wdpa_raw <- st_read(wdpa_clean_path, quiet = TRUE)
cat(sprintf("Loaded %d WDPA polygons\n", nrow(wdpa_raw)))

# Calculate areas in original CRS (should be WGS84)
wdpa_raw <- st_transform(wdpa_raw, MOLLWEIDE_CRS)
wdpa_raw$area_km2 <- as.numeric(st_area(wdpa_raw)) / 1e6

total_count_raw <- nrow(wdpa_raw)
total_area_raw <- sum(wdpa_raw$area_km2)

cat(sprintf("\nOriginal WDPA:\n"))
cat(sprintf("  Count: %s polygons\n", format(total_count_raw, big.mark = ",")))
cat(sprintf("  Total area: %s km^2\n", format(round(total_area_raw), big.mark = ",")))

# Store results
results <- data.frame(
  step = character(),
  count_before = integer(),
  count_after = integer(),
  count_dropped = integer(),
  pct_count_dropped = numeric(),
  area_before_km2 = numeric(),
  area_after_km2 = numeric(),
  area_dropped_km2 = numeric(),
  pct_area_dropped = numeric(),
  stringsAsFactors = FALSE
)

# ==============================================================================
# STEP 1: BBOX FILTERING
# ==============================================================================

cat("\n", rep("-", 40), "\n", sep = "")
cat("STEP 1: BBOX Filtering\n")
cat(rep("-", 40), "\n", sep = "")
cat("Filters geometries with bounding boxes spanning >90 degrees\n")
cat("(typically corrupt/invalid source geometries)\n\n")

max_bbox_span <- 90  # degrees

# Need to check in WGS84
wdpa_wgs84 <- st_transform(wdpa_raw, WGS84_CRS)

valid_bbox_idx <- sapply(seq_len(nrow(wdpa_wgs84)), function(i) {
  bb <- st_bbox(wdpa_wgs84[i, ])
  lon_span <- bb["xmax"] - bb["xmin"]
  lat_span <- bb["ymax"] - bb["ymin"]
  lon_span <= max_bbox_span && lat_span <= max_bbox_span
})

n_bbox_filtered <- sum(!valid_bbox_idx)
area_bbox_filtered <- sum(wdpa_raw$area_km2[!valid_bbox_idx])

cat(sprintf("Dropped: %d polygons (%.2f%% of count)\n",
            n_bbox_filtered, 100 * n_bbox_filtered / total_count_raw))
cat(sprintf("Dropped: %.2f km^2 (%.4f%% of area)\n",
            area_bbox_filtered, 100 * area_bbox_filtered / total_area_raw))

if (n_bbox_filtered > 0 && n_bbox_filtered <= 20) {
  cat("\nDropped polygons:\n")
  for (i in which(!valid_bbox_idx)) {
    bb <- st_bbox(wdpa_wgs84[i, ])
    cat(sprintf("  Row %d: bbox [%.2f,%.2f] x [%.2f,%.2f], area=%.2f km^2\n",
                i, bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"],
                wdpa_raw$area_km2[i]))
  }
}

wdpa_step1 <- wdpa_raw[valid_bbox_idx, ]

results <- rbind(results, data.frame(
  step = "1_bbox_filter",
  count_before = total_count_raw,
  count_after = nrow(wdpa_step1),
  count_dropped = n_bbox_filtered,
  pct_count_dropped = 100 * n_bbox_filtered / total_count_raw,
  area_before_km2 = total_area_raw,
  area_after_km2 = sum(wdpa_step1$area_km2),
  area_dropped_km2 = area_bbox_filtered,
  pct_area_dropped = 100 * area_bbox_filtered / total_area_raw
))

rm(wdpa_wgs84)
gc()

# ==============================================================================
# STEP 2: ST_SIMPLIFY (degenerate vertices)
# ==============================================================================

cat("\n", rep("-", 40), "\n", sep = "")
cat("STEP 2: st_simplify (dTolerance=1m)\n")
cat(rep("-", 40), "\n", sep = "")
cat("Removes degenerate vertices (zero-length segments) that cause\n")
cat("TopologyException in GEOS st_union operations\n\n")

count_before_simplify <- nrow(wdpa_step1)
area_before_simplify <- sum(wdpa_step1$area_km2)

wdpa_step2 <- st_simplify(wdpa_step1, dTolerance = 1, preserveTopology = TRUE)

# Recalculate areas after simplification
wdpa_step2$area_km2 <- as.numeric(st_area(wdpa_step2)) / 1e6

# Check for empty geometries created by simplification
empty_after_simplify <- st_is_empty(wdpa_step2)
n_empty <- sum(empty_after_simplify)

if (n_empty > 0) {
  area_empty <- sum(wdpa_step1$area_km2[empty_after_simplify])
  cat(sprintf("Geometries became empty: %d (%.2f%% of count)\n",
              n_empty, 100 * n_empty / count_before_simplify))
  cat(sprintf("Area lost to empty: %.2f km^2 (%.4f%% of area)\n",
              area_empty, 100 * area_empty / area_before_simplify))
  wdpa_step2 <- wdpa_step2[!empty_after_simplify, ]
}

# Calculate area change due to simplification (not counting empties)
area_change <- area_before_simplify - sum(wdpa_step2$area_km2)
cat(sprintf("Area change from simplification: %.2f km^2 (%.4f%% of area)\n",
            area_change, 100 * area_change / area_before_simplify))

cat(sprintf("Remaining: %d polygons, %.2f km^2\n",
            nrow(wdpa_step2), sum(wdpa_step2$area_km2)))

results <- rbind(results, data.frame(
  step = "2_st_simplify",
  count_before = count_before_simplify,
  count_after = nrow(wdpa_step2),
  count_dropped = count_before_simplify - nrow(wdpa_step2),
  pct_count_dropped = 100 * (count_before_simplify - nrow(wdpa_step2)) / count_before_simplify,
  area_before_km2 = area_before_simplify,
  area_after_km2 = sum(wdpa_step2$area_km2),
  area_dropped_km2 = area_before_simplify - sum(wdpa_step2$area_km2),
  pct_area_dropped = 100 * (area_before_simplify - sum(wdpa_step2$area_km2)) / area_before_simplify
))

# ==============================================================================
# STEP 3: GEOMETRYCOLLECTION EXTRACTION
# ==============================================================================

cat("\n", rep("-", 40), "\n", sep = "")
cat("STEP 3: GEOMETRYCOLLECTION Handling\n")
cat(rep("-", 40), "\n", sep = "")
cat("Extracts only POLYGON components from GEOMETRYCOLLECTIONs\n")
cat("(drops points, lines, and other non-polygon components)\n\n")

# First, let's see geometry type distribution
geom_types <- st_geometry_type(wdpa_step2)
geom_type_table <- table(geom_types)

cat("Geometry type distribution:\n")
for (gtype in names(geom_type_table)) {
  cat(sprintf("  %s: %d\n", gtype, geom_type_table[[gtype]]))
}

count_before_gc <- nrow(wdpa_step2)
area_before_gc <- sum(wdpa_step2$area_km2)

# Process GEOMETRYCOLLECTIONs
gc_idx <- which(geom_types == "GEOMETRYCOLLECTION")

if (length(gc_idx) > 0) {
  cat(sprintf("\nProcessing %d GEOMETRYCOLLECTIONs...\n", length(gc_idx)))

  wdpa_geom <- st_geometry(wdpa_step2)
  area_before_gc_extraction <- sum(wdpa_step2$area_km2[gc_idx])

  for (i in gc_idx) {
    gc_parts <- st_collection_extract(wdpa_geom[i], "POLYGON")
    if (length(gc_parts) > 0 && !st_is_empty(gc_parts)) {
      wdpa_geom[i] <- st_union(gc_parts)
    } else {
      wdpa_geom[i] <- st_polygon()
    }
  }
  st_geometry(wdpa_step2) <- wdpa_geom

  # Recalculate areas
  wdpa_step2$area_km2 <- as.numeric(st_area(wdpa_step2)) / 1e6
  area_after_gc_extraction <- sum(wdpa_step2$area_km2[gc_idx])

  cat(sprintf("GC polygons area before: %.2f km^2\n", area_before_gc_extraction))
  cat(sprintf("GC polygons area after: %.2f km^2\n", area_after_gc_extraction))
  cat(sprintf("Area change: %.2f km^2\n", area_before_gc_extraction - area_after_gc_extraction))

  # Remove empty geometries
  empty_after_gc <- st_is_empty(wdpa_step2)
  if (any(empty_after_gc)) {
    cat(sprintf("Dropping %d geometries with no polygon components\n", sum(empty_after_gc)))
    wdpa_step2 <- wdpa_step2[!empty_after_gc, ]
  }
} else {
  cat("\nNo GEOMETRYCOLLECTIONs found.\n")
}

wdpa_step3 <- wdpa_step2

results <- rbind(results, data.frame(
  step = "3_gc_extraction",
  count_before = count_before_gc,
  count_after = nrow(wdpa_step3),
  count_dropped = count_before_gc - nrow(wdpa_step3),
  pct_count_dropped = 100 * (count_before_gc - nrow(wdpa_step3)) / count_before_gc,
  area_before_km2 = area_before_gc,
  area_after_km2 = sum(wdpa_step3$area_km2),
  area_dropped_km2 = area_before_gc - sum(wdpa_step3$area_km2),
  pct_area_dropped = 100 * (area_before_gc - sum(wdpa_step3$area_km2)) / area_before_gc
))

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n", rep("=", 71), "\n", sep = "")
cat("SUMMARY\n")
cat(rep("=", 71), "\n\n", sep = "")

cat(sprintf("%-20s %12s %12s %12s %10s\n",
            "Step", "Count Drop", "% Count", "Area Drop", "% Area"))
cat(rep("-", 70), "\n", sep = "")

for (i in seq_len(nrow(results))) {
  cat(sprintf("%-20s %12s %11.4f%% %10.0f %9.4f%%\n",
              results$step[i],
              format(results$count_dropped[i], big.mark = ","),
              results$pct_count_dropped[i],
              results$area_dropped_km2[i],
              results$pct_area_dropped[i]))
}

cat(rep("-", 70), "\n", sep = "")
total_count_dropped <- total_count_raw - nrow(wdpa_step3)
total_area_dropped <- total_area_raw - sum(wdpa_step3$area_km2)
cat(sprintf("%-20s %12s %11.4f%% %10.0f %9.4f%%\n",
            "TOTAL",
            format(total_count_dropped, big.mark = ","),
            100 * total_count_dropped / total_count_raw,
            total_area_dropped,
            100 * total_area_dropped / total_area_raw))

cat("\n")
cat(sprintf("Original: %s polygons, %s km^2\n",
            format(total_count_raw, big.mark = ","),
            format(round(total_area_raw), big.mark = ",")))
cat(sprintf("Final:    %s polygons, %s km^2\n",
            format(nrow(wdpa_step3), big.mark = ","),
            format(round(sum(wdpa_step3$area_km2)), big.mark = ",")))

# ==============================================================================
# WRITE RESULTS
# ==============================================================================

output_file <- file.path(project_root, "LOGS", "wdpa_cleaning_impact.csv")
write.csv(results, output_file, row.names = FALSE)
cat(sprintf("\nDetailed results written to: %s\n", output_file))
