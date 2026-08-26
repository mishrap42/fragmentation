# ==============================================================================
# STAGE 1b: ASSEMBLE EARTHSTAT CROPLAND CROSS-SECTION
# ==============================================================================
# Consolidates the per-sub-tile Stage 0c outputs into one cell-level file.
#
# Input:  Data/build/yields/yields_sub_####.parquet   (N_SUB_TILES files)
#         Data/build/grids/*.parquet                  (for cell metadata)
# Output: Data/build/final/TMF_5km_cropland.parquet   (one row per grid_id)
#
# WHERE THIS SITS IN THE PIPELINE
# -------------------------------
# Numbered 1b, immediately after the 0c extraction it consolidates, because it
# belongs to the EXTRACTION side of the build - not the assembly side. Its only
# prerequisites are Stage 0 (grids) and Stage 0c (EarthStat rasters). It touches
# no TMF output, so it does not wait on Stages 1-4 and nothing downstream waits
# on it. Compare 2b_consolidate_gfed.R, which stands in the same relation to the
# 1a GFED extraction.
#
# It is emphatically NOT part of Stage 6. Stage 6 builds the cell-YEAR panel;
# this builds a cell-level CROSS-SECTION and the two are joined on grid_id at
# analysis time. EarthStat is a year-2000 snapshot: merging its 345 columns into
# 66.8M cell-year rows would replicate the same numbers 34 times (~184 GB
# uncompressed) and add no information. Keeping them apart is the whole point of
# splitting this out, so do not "helpfully" fold it back into 6_assemble_final.R.
#
# WHAT COMES OUT
# --------------
#   yield_<crop>       tons/ha ON THE CROP'S HARVESTED AREA (172 crops).
#                      This is a conditional-on-being-grown yield, so it is
#                      missing/zero exactly where the crop is not grown. It is
#                      NOT a potential yield and must not be read as one - for
#                      an opportunity cost on land that is currently forest,
#                      GAEZ potential yields are the right object.
#   cropshare_<crop>   harvested area as a fraction of CELL area (172 crops).
#   cropland_frac      cropland as a fraction of cell area (Ramankutty 2008).
#   pasture_frac       pasture as a fraction of cell area.
#   cropshare_total    sum of the 172 cropshare_* (computed in 0c).
#   n_crops_present    count of crops with cropshare_* > 0.
#   multicrop_index    cropshare_total / cropland_frac.
#
# NORMALIZATION IS LEFT TO THE ANALYST - deliberately. The two readings of a
# "crop land use share" are different objects and the choice is a modelling
# decision, not a build decision:
#
#   s_c = cropshare_c                        -> USD per hectare OF CELL
#   s_c = cropshare_c / cropshare_total      -> USD per hectare OF CROPLAND
#
# The second sums to 1 by construction and is what the companion repo's
# crop/6_log_relative_share.R uses; the first keeps the extensive margin.
#
# Usage: sbatch code/bash/1b_assemble_cropland.sh
# ==============================================================================

here::i_am('code/build/1b_assemble_cropland.R')
source("code/build/BUILD_workspace.R")

start_time <- Sys.time()
log_job_start("1b_assemble_cropland.R", task_id = 1)

output_file <- cropland_final_file
skip_if_exists(output_file, "cropland assembly")

# ==============================================================================
# LOAD STAGE 0c OUTPUTS
# ==============================================================================

log_message("Loading yield / crop-share sub-tile files...")

yield_files <- sapply(1:N_SUB_TILES, get_yields_filename)
files_exist <- file.exists(yield_files)

if (sum(files_exist) == 0) {
  stop(sprintf("No Stage 0c outputs found in %s. Run 0c_yields.sh first.",
               yields_output_path))
}

log_message(sprintf("Found %d/%d sub-tile files", sum(files_exist), N_SUB_TILES))

if (sum(files_exist) < N_SUB_TILES) {
  missing <- which(!files_exist)
  log_message(sprintf("  WARNING: %d sub-tiles MISSING: %s",
                      length(missing),
                      paste(head(missing, 25), collapse = ", ")))
  log_message("  Output will be incomplete. Use rerun_missing.sh before trusting it.")
}

# Track empty files separately: a 964-byte parquet with zero rows is the exact
# signature of the 2026-08-19 failure where the EarthStat mount was unreadable
# and every task "succeeded". file.exists() cannot see that; row counts can.
empty_subtiles <- integer(0)

cropland_list <- lapply(which(files_exist), function(i) {
  dt <- read_parquet(yield_files[i])
  setDT(dt)
  if (nrow(dt) == 0) empty_subtiles <<- c(empty_subtiles, i)
  dt
})

cropland_all <- rbindlist(cropland_list, fill = TRUE)
rm(cropland_list); gc()

if (length(empty_subtiles) > 0) {
  log_message(sprintf("  WARNING: %d sub-tile files had ZERO rows: %s",
                      length(empty_subtiles),
                      paste(head(empty_subtiles, 25), collapse = ", ")))
}

log_message(sprintf("Cropland data: %s cells, %d variables",
                    format(nrow(cropland_all), big.mark = ","),
                    ncol(cropland_all) - 1))

if (nrow(cropland_all) == 0) {
  stop("All Stage 0c outputs are empty. Check that yields_path is readable ",
       "from the compute nodes - this is the SMB-symlink failure mode.")
}

# ==============================================================================
# ATTACH GAEZ ATTAINABLE YIELDS (Stage 0d)
# ==============================================================================
# Joined here rather than in a third file because it shares the grain exactly -
# one row per grid_id, a time-invariant agronomic cross-section - and every
# consumer wants it alongside the EarthStat layers. Optional: 0d is independent
# of 0c and may not have run.

gaez_files <- sapply(1:N_SUB_TILES, get_gaez_filename)
gaez_exists <- file.exists(gaez_files)

if (sum(gaez_exists) == 0) {
  log_message("No Stage 0d (GAEZ) outputs found; cross-section will carry EarthStat only.")
} else {
  log_message(sprintf("Loading GAEZ sub-tile files (%d/%d present)...",
                      sum(gaez_exists), N_SUB_TILES))
  if (sum(gaez_exists) < N_SUB_TILES) {
    log_message(sprintf("  WARNING: %d GAEZ sub-tiles MISSING",
                        N_SUB_TILES - sum(gaez_exists)))
  }

  gaez_all <- rbindlist(lapply(gaez_files[gaez_exists], function(f) {
    dt <- read_parquet(f); setDT(dt); dt
  }), fill = TRUE)

  log_message(sprintf("  GAEZ: %s cells, %d variables",
                      format(nrow(gaez_all), big.mark = ","),
                      ncol(gaez_all) - 1))

  if (nrow(gaez_all) > 0) {
    # An inner-join-shaped left join: cells absent from 0d keep NA GAEZ columns
    # rather than being dropped, so a partial 0d cannot silently shrink the
    # cross-section.
    cropland_all <- merge(cropland_all, gaez_all, by = "grid_id", all.x = TRUE)
    n_nogaez <- cropland_all[, sum(is.na(get(setdiff(names(gaez_all), "grid_id")[1])))]
    log_message(sprintf("  cells with no GAEZ row: %s (%.2f%%)",
                        format(n_nogaez, big.mark = ","),
                        100 * n_nogaez / nrow(cropland_all)))
  }
  rm(gaez_all); gc()

  # ----------------------------------------------------------------------------
  # Composite GAEZ layers where one crop maps to several GAEZ codes.
  # ----------------------------------------------------------------------------
  # The crosswalk's gaez field is pipe-separated when GAEZ ships no aggregate
  # layer for the crop. Rice is the case that forces this: GAEZ has ricd
  # (dryland) and ricw (wetland) and no combined `rice`, and the attainable
  # yield on a cell is the MAX of the two - a grower adopts the better system
  # available there, so neither picking one nor averaging is right.
  #
  # pmax(na.rm = TRUE) returns NA only where every constituent is NA, so a cell
  # suited to just one system keeps that system's yield rather than going NA.
  xw_g <- fread(crop_crosswalk_file)
  multi <- xw_g[nzchar(gaez) & grepl("|", gaez, fixed = TRUE)]

  for (i in seq_len(nrow(multi))) {
    crop  <- multi$earthstat[i]
    parts <- strsplit(multi$gaez[i], "|", fixed = TRUE)[[1]]
    for (v in c("200a", "200b")) {
      src <- intersect(paste0("gaez_", parts, v), names(cropland_all))
      if (length(src) < 2) {
        log_message(sprintf("  WARNING: %s %s - only %d of %d source layers present, skipping composite",
                            crop, v, length(src), length(parts)))
        next
      }
      tgt <- paste0("gaez_", crop, v)
      cropland_all[, (tgt) := do.call(pmax, c(.SD, list(na.rm = TRUE))),
                   .SDcols = src]
      log_message(sprintf("  composite %s = pmax(%s)  [mean %.1f kg/ha]",
                          tgt, paste(src, collapse = ", "),
                          mean(cropland_all[[tgt]], na.rm = TRUE)))
    }
  }
}

# ==============================================================================
# ATTACH CELL METADATA
# ==============================================================================
# Carried so the cross-section stands alone for descriptives and mapping without
# having to touch the 4.6 GB panel. The panel remains the source of truth for
# anything time-varying.

log_message("Loading grid metadata...")

grid_cols <- c("grid_id", "tile_id", "country_iso3", "country_name",
               "centroid_lon", "centroid_lat", "area_km2")

grid_files <- sapply(1:N_SUB_TILES, function(s) get_grid_filename(s, "parquet"))
grid_files <- grid_files[file.exists(grid_files)]

grid_all <- rbindlist(lapply(grid_files, function(f) {
  dt <- tryCatch(read_parquet(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  setDT(dt)
  dt[, intersect(grid_cols, names(dt)), with = FALSE]
}), fill = TRUE)

log_message(sprintf("  %s grid cells with metadata",
                    format(nrow(grid_all), big.mark = ",")))

cropland_all <- merge(grid_all, cropland_all, by = "grid_id", all.x = TRUE)
rm(grid_all); gc()

log_message(sprintf("After metadata join: %s cells",
                    format(nrow(cropland_all), big.mark = ",")))

# ==============================================================================
# DERIVED COLUMNS
# ==============================================================================

gaez_cols <- grep("^gaez_", names(cropland_all), value = TRUE)
if (length(gaez_cols) > 0) {
  log_message(sprintf("GAEZ layers present: %d (a: %d, b: %d)",
                      length(gaez_cols),
                      sum(grepl("200a$", gaez_cols)),
                      sum(grepl("200b$", gaez_cols))))
}

cropshare_cols <- grep("^cropshare_", names(cropland_all), value = TRUE)
cropshare_cols <- setdiff(cropshare_cols, "cropshare_total")
yield_cols     <- grep("^yield_", names(cropland_all), value = TRUE)

log_message(sprintf("Crop layers present: %d cropshare_*, %d yield_*",
                    length(cropshare_cols), length(yield_cols)))

if (length(cropshare_cols) == 0) {
  stop("No cropshare_* columns found. Stage 0c ran against the wrong source ",
       "tree - yields_path must point at the EarthStat 175-crop GeoTiff dir, ",
       "not at GAEZ.")
}

# cropshare_total is written by 0c; recompute only if an older run predates it.
if (!"cropshare_total" %in% names(cropland_all)) {
  log_message("cropshare_total absent (pre-2026-08-21 Stage 0c); recomputing...")
  cropland_all[, cropshare_total := rowSums(.SD, na.rm = TRUE),
               .SDcols = cropshare_cols]
}

log_message("Computing derived columns...")

cropland_all[, n_crops_present := rowSums(.SD > 0, na.rm = TRUE),
             .SDcols = cropshare_cols]

# Multi-cropping intensity. > 1 is expected and correct: EarthStat counts a
# hectare once per harvest, so double-cropped land shows up twice in the
# numerator and once in the denominator. Guarded against the divide-by-zero that
# covers most of the tropical-forest sample, where cropland_frac is 0.
if ("cropland_frac" %in% names(cropland_all)) {
  cropland_all[, multicrop_index := fifelse(
    is.na(cropland_frac) | cropland_frac <= 0,
    NA_real_,
    cropshare_total / cropland_frac
  )]
} else {
  log_message("  WARNING: cropland_frac absent; multicrop_index not computed.")
  cropland_all[, multicrop_index := NA_real_]
}

# ==============================================================================
# DIAGNOSTICS
# ==============================================================================
# The failure this guards against is silent: a shifted CRS or a bad mount gives
# all-NA or all-zero columns that merge cleanly and estimate cleanly, and only
# look wrong in a coefficient nobody can attribute. Check coverage here.

log_message("========================================")
log_message("COVERAGE DIAGNOSTICS")
log_message("========================================")

n_cells <- nrow(cropland_all)

report_coverage <- function(col) {
  if (!col %in% names(cropland_all)) {
    log_message(sprintf("  %-18s ABSENT", col))
    return(invisible(NULL))
  }
  v <- cropland_all[[col]]
  log_message(sprintf("  %-18s non-NA %5.1f%%  >0 %5.1f%%  mean %.5f  max %.4f",
                      col,
                      100 * sum(!is.na(v)) / n_cells,
                      100 * sum(v > 0, na.rm = TRUE) / n_cells,
                      mean(v, na.rm = TRUE),
                      suppressWarnings(max(v, na.rm = TRUE))))
}

for (col in c("cropland_frac", "pasture_frac", "cropshare_total",
              "n_crops_present", "multicrop_index")) {
  report_coverage(col)
}

log_message("Top 15 crops by mean cell share:")
crop_means <- sort(sapply(cropshare_cols,
                          function(c) mean(cropland_all[[c]], na.rm = TRUE)),
                   decreasing = TRUE)
for (nm in head(names(crop_means), 15)) {
  log_message(sprintf("  %-28s %.6f", nm, crop_means[[nm]]))
}

# All-NA columns are the loud version of the silent failure above.
all_na <- names(which(sapply(c(cropshare_cols, yield_cols),
                            function(c) all(is.na(cropland_all[[c]])))))
if (length(all_na) > 0) {
  log_message(sprintf("  WARNING: %d columns are entirely NA: %s",
                      length(all_na), paste(head(all_na, 20), collapse = ", ")))
} else {
  log_message("  No all-NA crop columns.")
}

# Cells the extraction never reached.
n_unmatched <- cropland_all[is.na(cropshare_total), .N]
log_message(sprintf("Cells with no Stage 0c row: %s (%.2f%%)",
                    format(n_unmatched, big.mark = ","),
                    100 * n_unmatched / n_cells))

# ==============================================================================
# WRITE OUTPUT
# ==============================================================================

setcolorder(cropland_all,
            c(intersect(grid_cols, names(cropland_all)),
              intersect(c("cropland_frac", "pasture_frac", "cropshare_total",
                          "n_crops_present", "multicrop_index"),
                        names(cropland_all))))
setorder(cropland_all, grid_id)

log_message(sprintf("Writing output: %s", output_file))
write_atomic(cropland_all, output_file)

if (file.exists(output_file)) {
  log_message(sprintf("Output file size: %.2f GB",
                      file.info(output_file)$size / 1024^3))
}

# ==============================================================================
# SUMMARY FILE
# ==============================================================================

summary_file <- file.path(final_output_path, "summary_cropland.txt")
log_message(sprintf("Writing summary: %s", summary_file))

sink(summary_file)
cat("EarthStat Cropland Cross-Section Summary\n")
cat("=======================================\n")
cat(sprintf("Generated: %s\n", Sys.time()))
cat(sprintf("Source: %s\n\n", yields_path))

cat(sprintf("Grid cells: %s\n", format(nrow(cropland_all), big.mark = ",")))
cat(sprintf("Countries: %d\n", length(unique(cropland_all$country_iso3))))
cat(sprintf("Crops (cropshare_*): %d\n", length(cropshare_cols)))
cat(sprintf("Crops (yield_*): %d\n", length(yield_cols)))
cat(sprintf("GAEZ layers (gaez_*): %d\n\n", length(gaez_cols)))

cat("Reference year: 2000 (EarthStat is a cross-section; join on grid_id)\n\n")

cat("Columns:\n")
for (col in names(cropland_all)) {
  cat(sprintf("  - %s (%s)\n", col, class(cropland_all[[col]])[1]))
}
sink()

gc_verbose()
log_job_end(start_time)

log_message("Cropland cross-section complete.")
