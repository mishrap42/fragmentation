# ==============================================================================
# VALIDATE YIELDS (Stage 0c) AND GFED BURNED AREA (Stage 1a)
# ==============================================================================
# Neither dataset had any value-level validation. diagnose_pipeline.R only calls
# file.exists() on the expected outputs - a check that passed happily on
# 2026-08-19 while all 303 yield files were 964-byte EMPTY parquets, which is
# exactly the failure it should have caught.
#
# This matters most for yields, which were rewired most recently: new source
# (EarthStat, replacing GAEZ), new column naming (yield_* / cropshare_*), new
# units. Three analysis scripts already consume those columns.
#
# Three sections: the Stage 0c sub-tile extraction, the Stage 1b assembled
# cross-section (skipped if 1b has not run), and Stage 1a GFED.
#
# Numeric checks cover ALL columns. Plots cover a curated subset - the land-use
# layers plus the top PLOT_TOP_N_CROPS crops by mean share - because the full
# 172-crop extraction would otherwise emit ~690 PNGs and overrun the wall clock.
# The skipped set is logged, never silently dropped.
#
# Outputs -> output/figures/yields_gfed/
#   summary_yields.txt / summary_gfed.txt   per-column checks with PASS/FAIL
#   yield_hist_<col>.png / yield_map_<col>.png   (curated subset)
#   gfed_annual_timeseries.png / gfed_map_mean_annual.png / gfed_hist.png
#
# Usage: sbatch code/bash/validate_yields_gfed.sh
# ==============================================================================

here::i_am("code/build/validate_yields_gfed.R")
source("code/build/BUILD_workspace.R")

library(ggplot2)

output_dir <- file.path(project_root, "output", "validation", "yields_gfed")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

start_time <- Sys.time()
log_message("========================================")
log_message("YIELD + GFED VALIDATION")
log_message("========================================")

# ==============================================================================
# GRID COORDINATES (both datasets are keyed on grid_id only)
# ==============================================================================

log_message("Loading grid centroids for mapping...")
grid_files <- sapply(1:N_SUB_TILES, function(i) get_grid_filename(i, "parquet"))
grid_files <- grid_files[file.exists(grid_files)]

grid_coords <- rbindlist(lapply(grid_files, function(f) {
  dt <- tryCatch(
    arrow::read_parquet(f, col_select = c("grid_id", "centroid_lon",
                                          "centroid_lat", "country_iso3")),
    error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  setDT(dt)
  dt
}), fill = TRUE)

log_message(sprintf("  %s grid cells with coordinates",
                    format(nrow(grid_coords), big.mark = ",")))

# ==============================================================================
# SHARED PLOT HELPERS (mirroring validate_covariates.R)
# ==============================================================================

create_histogram <- function(vals, col, output_path, unit = "") {
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) {
    log_message(sprintf("  Skipping histogram for %s (no valid values)", col))
    return(invisible(NULL))
  }
  p <- ggplot(data.frame(x = vals), aes(x = x)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
    labs(title = sprintf("Distribution of %s", col),
         subtitle = sprintf("N = %s, Mean = %.3f, Median = %.3f, Max = %.3f",
                            format(length(vals), big.mark = ","),
                            mean(vals), median(vals), max(vals)),
         x = if (nzchar(unit)) sprintf("%s (%s)", col, unit) else col,
         y = "Count") +
    theme_minimal() + theme(plot.title = element_text(face = "bold"))
  ggsave(output_path, p, width = 8, height = 5, dpi = 150)
}

create_map <- function(dt, col, output_path, title_prefix = "Spatial distribution") {
  plot_data <- dt[!is.na(get(col)) & !is.na(centroid_lon)]
  if (nrow(plot_data) == 0) {
    log_message(sprintf("  Skipping map for %s (no valid values)", col))
    return(invisible(NULL))
  }
  if (nrow(plot_data) > 100000) {
    set.seed(42)
    plot_data <- plot_data[sample(.N, 100000)]
  }
  vals <- plot_data[[col]]
  q01 <- quantile(vals, 0.01, na.rm = TRUE)
  q99 <- quantile(vals, 0.99, na.rm = TRUE)
  if (!is.finite(q01) || !is.finite(q99) || q99 <= q01) { q01 <- min(vals); q99 <- max(vals) }

  p <- ggplot(plot_data, aes(x = centroid_lon, y = centroid_lat,
                             color = .data[[col]])) +
    geom_point(size = 0.1, alpha = 0.5) +
    scale_color_viridis_c(name = col, limits = c(q01, q99), oob = scales::squish) +
    coord_fixed(ratio = 1.3) +
    labs(title = sprintf("%s: %s", title_prefix, col),
         x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
  ggsave(output_path, p, width = 12, height = 8, dpi = 150)
}

# Accumulates check rows; every check records PASS/FAIL rather than stopping, so
# one bad column does not hide the state of the other 24.
checks <- list()
add_check <- function(dataset, item, test, pass, detail = "") {
  checks[[length(checks) + 1]] <<- data.table(
    dataset = dataset, item = item, test = test,
    result = if (isTRUE(pass)) "PASS" else "FAIL", detail = detail)
}

# ==============================================================================
# SECTION 1: YIELDS
# ==============================================================================

log_message("----------------------------------------")
log_message("SECTION 1: YIELDS (EarthStat)")
log_message("----------------------------------------")

yield_files <- sapply(1:N_SUB_TILES, get_yields_filename)
yield_exists <- file.exists(yield_files)
log_message(sprintf("Found %d/%d yield files", sum(yield_exists), N_SUB_TILES))

yields <- rbindlist(lapply(yield_files[yield_exists], function(f) {
  dt <- tryCatch(arrow::read_parquet(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  setDT(dt); dt
}), fill = TRUE)

log_message(sprintf("Loaded %s rows, %d columns",
                    format(nrow(yields), big.mark = ","), ncol(yields)))

yield_cols <- setdiff(names(yields), "grid_id")
y_cols <- grep("^yield_", yield_cols, value = TRUE)
s_cols <- setdiff(grep("^cropshare_", yield_cols, value = TRUE), "cropshare_total")

# Cell-area fractions that must lie in [0,1]: the 172 per-crop harvested-area
# fractions plus the two EarthStat land-use layers. cropshare_total is
# deliberately EXCLUDED - it is a sum over crops and legitimately exceeds 1
# wherever land is double-cropped, because EarthStat counts a hectare once per
# harvest. Bounds-checking it would fail on correct data.
frac_cols <- c(s_cols, intersect(c("cropland_frac", "pasture_frac"), yield_cols))

log_message(sprintf("  %d yield_* columns, %d cropshare_* columns, %d bounded fractions",
                    length(y_cols), length(s_cols), length(frac_cols)))

# The empty-output failure mode: files present but carrying no value columns.
add_check("yields", "columns", "value columns present",
          length(yield_cols) > 0,
          sprintf("%d value columns", length(yield_cols)))
add_check("yields", "rows", "non-empty",
          nrow(yields) > 0, sprintf("%s rows", format(nrow(yields), big.mark = ",")))

# ------------------------------------------------------------------------------
# WHICH COLUMNS GET PLOTTED
# ------------------------------------------------------------------------------
# Numeric checks run on ALL columns. Plots do not. At 25 columns (the 14-crop
# partial extraction this script was written against) plotting everything cost
# 50 PNGs; at the full 172-crop extraction it is ~690, which overruns the 3h
# wall and buries the informative figures among 600 near-empty ones.
#
# So: the land-use/summary layers always, plus the top N crops by mean cell
# share and their matching yield layers. What is dropped is logged rather than
# silently truncated - a validator that quietly narrows its own coverage is the
# failure mode this whole script exists to catch.

PLOT_TOP_N_CROPS <- 12

# Upper bound for a harvested-area fraction; see the bounds check below.
MULTICROP_CEILING <- 4

# Crop columns that are constant (almost always all-zero: crop absent from the
# tropics). Collected during the loop and reported once.
constant_cols <- character(0)

always_plot <- intersect(c("cropland_frac", "pasture_frac", "cropshare_total"),
                         yield_cols)

if (length(s_cols) > 0) {
  crop_rank  <- sort(sapply(s_cols, function(c) mean(yields[[c]], na.rm = TRUE)),
                     decreasing = TRUE)
  top_crops  <- sub("^cropshare_", "", head(names(crop_rank), PLOT_TOP_N_CROPS))
} else {
  top_crops <- character(0)
}

plot_cols <- intersect(unique(c(always_plot,
                                paste0("cropshare_", top_crops),
                                paste0("yield_", top_crops))),
                       yield_cols)

log_message(sprintf("Plotting %d of %d columns (top %d crops by mean share + %d land-use layers)",
                    length(plot_cols), length(yield_cols),
                    length(top_crops), length(always_plot)))
log_message(sprintf("  plotted: %s", paste(plot_cols, collapse = ", ")))
log_message(sprintf("  NOT plotted (checked numerically only): %d columns",
                    length(setdiff(yield_cols, plot_cols))))

# Only the plotted columns are joined to coordinates. Merging the full 345-column
# table would duplicate ~5.5 GB for the sake of a dozen scatterplots.
yields_geo <- merge(yields[, c("grid_id", plot_cols), with = FALSE],
                    grid_coords, by = "grid_id", all.x = TRUE)

summary_rows <- list()
for (col in yield_cols) {
  v <- yields[[col]]
  n_valid <- sum(!is.na(v))
  pct <- 100 * n_valid / length(v)
  vv <- v[!is.na(v)]

  if (n_valid == 0) {
    add_check("yields", col, "has non-NA values", FALSE, "all NA")
  } else {
    add_check("yields", col, "has non-NA values", TRUE, sprintf("%.1f%% non-NA", pct))
    # A constant column is a genuine bug for the land-use layers, which must vary
    # across a global tropical grid. For an individual CROP it is expected: this
    # grid is tropical moist forest and the 172-crop set is global, so carob,
    # kiwi, chicory, peppermint and the fodder crops are legitimately absent
    # everywhere. Counted and listed below rather than failed 25 times.
    if (col %in% c("cropland_frac", "pasture_frac", "cropshare_total")) {
      add_check("yields", col, "non-degenerate (varies)",
                length(unique(vv)) > 1,
                sprintf("%d distinct values", length(unique(vv))))
    } else if (length(unique(vv)) <= 1) {
      constant_cols <<- c(constant_cols, col)
    }
    # Harvested-area fractions of cell area. The upper bound is NOT 1: EarthStat
    # counts a hectare once per harvest, so a double- or triple-cropped cell
    # exceeds 1 for the individual crop too, not just for the total. The 2026-08-22
    # run hit this on exactly the crops you would predict - rice 1.93, cotton 1.82,
    # soybean 1.62, maize 1.58, sugarcane 1.44 - all real multi-cropping.
    # MULTICROP_CEILING sits above any plausible harvest count and far below what
    # a HarvestedAreaHectares mix-up would produce (a 5 km cell is ~2500 ha), which
    # is the error this check actually exists to catch.
    if (col %in% frac_cols) {
      add_check("yields", col,
                sprintf("bounded in [0,%g]", MULTICROP_CEILING),
                min(vv) >= -1e-9 && max(vv) <= MULTICROP_CEILING + 1e-9,
                sprintf("range [%.4f, %.4f]", min(vv), max(vv)))
    }

    # cropland_frac / pasture_frac are true land-cover shares with no harvest
    # dimension, so they ARE bounded by 1.
    if (col %in% c("cropland_frac", "pasture_frac")) {
      add_check("yields", col, "land-cover share bounded in [0,1]",
                max(vv) <= 1 + 1e-9, sprintf("max = %.4f", max(vv)))
    }
    # EarthStat YieldPerHectare is tons/ha. Single digits for most grains, but
    # the 172-crop set includes fresh-weight sugar and fodder crops (sugarcane,
    # beetfor, cabbagefor, alfalfa, grass) that legitimately reach the low
    # hundreds. The failure this catches is a UNIT error - a kg/ha layer read as
    # t/ha is off by 1000x - so the ceiling only has to sit above real fresh
    # weights and far below 1000x. 200 was calibrated against the 14-crop
    # partial extraction and would now FAIL on correct fodder data.
    if (startsWith(col, "yield_")) {
      add_check("yields", col, "non-negative", min(vv) >= -1e-9,
                sprintf("min = %.4f", min(vv)))
      add_check("yields", col, "plausible magnitude (max < 500 t/ha)",
                max(vv) < 500, sprintf("max = %.2f", max(vv)))
    }
  }

  summary_rows[[col]] <- data.table(
    column = col, n_valid = n_valid, pct_valid = round(pct, 1),
    min = if (n_valid) round(min(vv), 4) else NA_real_,
    median = if (n_valid) round(median(vv), 4) else NA_real_,
    mean = if (n_valid) round(mean(vv), 4) else NA_real_,
    max = if (n_valid) round(max(vv), 4) else NA_real_)

  log_message(sprintf("  %-26s %6.1f%% valid  median=%9.4f  max=%9.4f",
                      col, pct,
                      if (n_valid) median(vv) else NA_real_,
                      if (n_valid) max(vv) else NA_real_))

  if (col %in% plot_cols) {
    create_histogram(v, col, file.path(output_dir, sprintf("yield_hist_%s.png", col)),
                     unit = if (startsWith(col, "yield_")) "t/ha" else "fraction")
    create_map(yields_geo, col, file.path(output_dir, sprintf("yield_map_%s.png", col)))
  }
}

yield_summary <- rbindlist(summary_rows)

log_message(sprintf("Constant (all-zero) crop columns: %d of %d",
                    length(constant_cols), length(yield_cols)))
if (length(constant_cols) > 0) {
  log_message(sprintf("  %s", paste(constant_cols, collapse = ", ")))
}
add_check("yields", "constant columns", "fewer than half the crop layers are all-zero",
          length(constant_cols) < 0.5 * length(yield_cols),
          sprintf("%d of %d constant", length(constant_cols), length(yield_cols)))

# ==============================================================================
# SECTION 1b: ASSEMBLED CROPLAND CROSS-SECTION (Stage 1b)
# ==============================================================================
# Section 1 checks the per-sub-tile extraction. This checks what 1b made of it:
# the consolidation can lose cells (a sub-tile that never ran) or gain NAs (a
# grid cell with no extraction row) without either side looking wrong alone.
# Skipped rather than failed when 1b has not run yet - 0c is validatable on its
# own and the two stages are submitted separately.

log_message("----------------------------------------")
log_message("SECTION 1b: CROPLAND CROSS-SECTION")
log_message("----------------------------------------")

if (!file.exists(cropland_final_file)) {
  log_message(sprintf("Not found, skipping: %s", cropland_final_file))
  log_message("  (run code/bash/1b_assemble_cropland.sh to produce it)")
} else {
  meta_cols <- c("grid_id", "cropland_frac", "pasture_frac", "cropshare_total",
                 "n_crops_present", "multicrop_index")
  # col_select goes through tidyselect, which deprecates bare external vectors -
  # a literal c(...) is fine but a variable warns. arrow re-exports all_of(), so
  # this is warning-free and needs no extra dependency.
  avail <- arrow::open_dataset(cropland_final_file, format = "parquet")$schema$names
  have  <- intersect(meta_cols, avail)
  cs <- setDT(arrow::read_parquet(cropland_final_file,
                                  col_select = arrow::all_of(have)))

  log_message(sprintf("Loaded %s cells, %d summary columns",
                      format(nrow(cs), big.mark = ","), ncol(cs)))

  add_check("cropland", "rows", "non-empty", nrow(cs) > 0,
            sprintf("%s rows", format(nrow(cs), big.mark = ",")))
  add_check("cropland", "grid_id", "unique (one row per cell)",
            !any(duplicated(cs$grid_id)),
            sprintf("%d duplicated", sum(duplicated(cs$grid_id))))

  # Every cell that Section 1 extracted must survive into the cross-section.
  n_extracted <- length(unique(yields$grid_id))
  n_covered   <- if ("cropshare_total" %in% names(cs)) cs[!is.na(cropshare_total), .N] else NA_integer_
  add_check("cropland", "coverage", "all extracted cells carried through",
            !is.na(n_covered) && n_covered >= n_extracted,
            sprintf("%s of %s extracted cells have data",
                    format(n_covered, big.mark = ","),
                    format(n_extracted, big.mark = ",")))

  for (col in intersect(c("cropland_frac", "pasture_frac"), names(cs))) {
    v <- cs[[col]]; vv <- v[!is.na(v)]
    add_check("cropland", col, "bounded in [0,1]",
              length(vv) > 0 && min(vv) >= -1e-9 && max(vv) <= 1 + 1e-9,
              sprintf("range [%.4f, %.4f]", min(vv), max(vv)))
  }

  if ("cropshare_total" %in% names(cs)) {
    v <- cs$cropshare_total; vv <- v[!is.na(v)]
    add_check("cropland", "cropshare_total", "non-negative",
              min(vv) >= -1e-9, sprintf("min = %.6f", min(vv)))
    # NOT bounded above by 1: EarthStat counts a hectare once per harvest, so
    # double-cropped cells legitimately exceed it. Reported, not failed.
    log_message(sprintf("  cropshare_total: median %.5f, max %.4f, %.2f%% of cells > 1",
                        median(vv), max(vv), 100 * mean(vv > 1)))
  }

  if ("multicrop_index" %in% names(cs)) {
    v <- cs$multicrop_index; vv <- v[is.finite(v)]
    add_check("cropland", "multicrop_index", "finite where defined",
              length(vv) > 0, sprintf("%s finite values",
                                      format(length(vv), big.mark = ",")))
    if (length(vv) > 0) {
      log_message(sprintf("  multicrop_index: median %.3f, p99 %.3f (>1 = double-cropping)",
                          median(vv), quantile(vv, 0.99)))
    }
  }

  if ("n_crops_present" %in% names(cs)) {
    log_message(sprintf("  n_crops_present: median %g, max %g",
                        median(cs$n_crops_present, na.rm = TRUE),
                        max(cs$n_crops_present, na.rm = TRUE)))
  }
  rm(cs); gc()
}

# ==============================================================================
# SECTION 2: GFED BURNED AREA
# ==============================================================================

log_message("----------------------------------------")
log_message("SECTION 2: GFED BURNED AREA")
log_message("----------------------------------------")

gfed_files <- sapply(1:N_SUB_TILES, get_gfed_subtile_filename)
gfed_exists <- file.exists(gfed_files)
log_message(sprintf("Found %d/%d GFED files", sum(gfed_exists), N_SUB_TILES))

# ~100M rows across all sub-tiles: aggregate per file rather than rbinding the
# lot, so peak memory is one sub-tile regardless of total size.
per_cell <- list(); per_year <- list()
n_rows_total <- 0; ba_min <- Inf; ba_max <- -Inf; n_na <- 0
yr_min <- Inf; yr_max <- -Inf; mo_bad <- 0

for (f in gfed_files[gfed_exists]) {
  dt <- tryCatch(arrow::read_parquet(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) next
  setDT(dt)

  n_rows_total <- n_rows_total + nrow(dt)
  n_na <- n_na + sum(is.na(dt$burned_area))
  ba <- dt$burned_area[!is.na(dt$burned_area)]
  if (length(ba)) { ba_min <- min(ba_min, min(ba)); ba_max <- max(ba_max, max(ba)) }
  if ("year" %in% names(dt)) {
    yr_min <- min(yr_min, min(dt$year, na.rm = TRUE))
    yr_max <- max(yr_max, max(dt$year, na.rm = TRUE))
  }
  if ("month" %in% names(dt)) {
    mo_bad <- mo_bad + sum(dt$month < 1 | dt$month > 12, na.rm = TRUE)
  }

  n_years <- length(unique(dt$year))
  per_cell[[length(per_cell) + 1]] <-
    dt[, .(mean_annual_burned = sum(burned_area, na.rm = TRUE) / max(n_years, 1)),
       by = grid_id]
  per_year[[length(per_year) + 1]] <-
    dt[, .(total_burned = sum(burned_area, na.rm = TRUE)), by = year]
}

cell_burn <- rbindlist(per_cell)[, .(mean_annual_burned = sum(mean_annual_burned)),
                                 by = grid_id]
year_burn <- rbindlist(per_year)[, .(total_burned = sum(total_burned)), by = year]
setorder(year_burn, year)

log_message(sprintf("  %s rows, %s cells with burns, years %d-%d",
                    format(n_rows_total, big.mark = ","),
                    format(nrow(cell_burn), big.mark = ","), yr_min, yr_max))

add_check("gfed", "rows", "non-empty", n_rows_total > 0,
          sprintf("%s rows", format(n_rows_total, big.mark = ",")))
add_check("gfed", "burned_area", "no NA values", n_na == 0, sprintf("%d NA", n_na))
add_check("gfed", "burned_area", "non-negative", is.finite(ba_min) && ba_min >= 0,
          sprintf("min = %s", format(ba_min)))
add_check("gfed", "year", "within GFED5 range 2002-2022",
          yr_min >= 2002 && yr_max <= 2022, sprintf("observed %d-%d", yr_min, yr_max))
add_check("gfed", "month", "within 1-12", mo_bad == 0,
          sprintf("%d out-of-range", mo_bad))

cell_burn_geo <- merge(cell_burn, grid_coords, by = "grid_id", all.x = TRUE)

create_histogram(cell_burn$mean_annual_burned, "mean_annual_burned",
                 file.path(output_dir, "gfed_hist.png"), unit = "per cell")
create_map(cell_burn_geo, "mean_annual_burned",
           file.path(output_dir, "gfed_map_mean_annual.png"),
           title_prefix = "GFED mean annual burned area")

p <- ggplot(year_burn, aes(x = year, y = total_burned)) +
  geom_line(color = "steelblue", linewidth = 0.9) +
  geom_point(color = "steelblue", size = 1.6) +
  labs(title = "GFED total burned area by year",
       subtitle = sprintf("%d-%d, summed across all grid cells", yr_min, yr_max),
       x = "Year", y = "Total burned area") +
  theme_minimal() + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(output_dir, "gfed_annual_timeseries.png"), p,
       width = 9, height = 5, dpi = 150)

# ==============================================================================
# SUMMARY
# ==============================================================================

check_dt <- rbindlist(checks)
n_fail <- sum(check_dt$result == "FAIL")

y_path <- file.path(output_dir, "summary_yields.txt")
sink(y_path)
cat("YIELD VALIDATION (EarthStat)\n")
cat(sprintf("Generated: %s\n", format(Sys.time())))
cat(sprintf("Files: %d/%d   Rows: %s   Columns: %d\n\n",
            sum(yield_exists), N_SUB_TILES,
            format(nrow(yields), big.mark = ","), length(yield_cols)))
# print.data.table truncates to head+tail past 100 rows. With 345 value columns
# and ~1000 checks that would silently hide most of the report - the default was
# harmless only while this stage produced 25 columns.
print(yield_summary, nrows = Inf)
cat("\nCHECKS\n")
print(check_dt[dataset == "yields"], nrows = Inf)
if (any(check_dt$dataset == "cropland")) {
  cat("\nCROPLAND CROSS-SECTION (Stage 1b)\n")
  print(check_dt[dataset == "cropland"], nrows = Inf)
}
sink()

g_path <- file.path(output_dir, "summary_gfed.txt")
sink(g_path)
cat("GFED BURNED AREA VALIDATION\n")
cat(sprintf("Generated: %s\n", format(Sys.time())))
cat(sprintf("Files: %d/%d   Rows: %s   Cells with burns: %s\n\n",
            sum(gfed_exists), N_SUB_TILES,
            format(n_rows_total, big.mark = ","),
            format(nrow(cell_burn), big.mark = ",")))
print(year_burn)
cat("\nCHECKS\n")
print(check_dt[dataset == "gfed"])
sink()

log_message("========================================")
if (n_fail == 0) {
  log_message(sprintf("ALL %d CHECKS PASSED", nrow(check_dt)))
} else {
  log_message(sprintf("%d of %d CHECKS FAILED:", n_fail, nrow(check_dt)))
  print(check_dt[result == "FAIL"])
}
log_message(sprintf("Figures + summaries: %s", output_dir))
log_message(sprintf("Completed in %.1f minutes",
                    as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
log_message("========================================")
