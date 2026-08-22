# ==============================================================================
# STAGE 0a: WDPA PREPROCESSING
# ==============================================================================
# This script cleans the WDPA (World Database on Protected Areas) file.
# Run ONCE before the main pipeline to create wdpa_cleaned.gpkg.
#
# Input: Raw WDPA GDB file
# Output: Data/build/wdpa_cleaned.gpkg
#
# Usage: Rscript code/build/0a_clean_wdpa.R
# ==============================================================================

here::i_am('code/build/0a_clean_wdpa.R')
source("code/build/BUILD_workspace.R")

# Record start time
start_time <- Sys.time()

log_message("========================================")
log_message("WDPA Preprocessing")
log_message("========================================")

# ==============================================================================
# CRASH SAFETY DESIGN
# ==============================================================================
# Job 1067292 was killed by the 4h wall clock partway through the GDB->GPKG
# conversion, leaving a 4.55 GB truncated wdpa_cleaned.gpkg that the old
# file.exists() guard would have accepted as finished. This version is built so
# that no failure - timeout, OOM, node loss - can either lose completed work or
# leave something a later run mistakes for a valid output.
#
# Three properties:
#
#   1. TWO CHECKPOINTED PHASES. Phase A (convert, the expensive part) writes
#      wdpa_raw.gpkg. Phase B (validate/repair, the memory-heavy part) reads it
#      and writes wdpa_cleaned.gpkg. If B dies - which is where an OOM would hit
#      - a rerun skips A entirely and resumes at B. The 4h+ conversion is paid
#      at most once.
#
#   2. ATOMIC PUBLICATION. Every phase writes to <target>.tmp.<jobid> and only
#      file.rename()s into place after the result verifies. A killed job leaves
#      a .tmp file, never a target file, so a partial output cannot be mistaken
#      for a complete one. rename() within one filesystem is atomic.
#
#   3. VERIFICATION BY FEATURE COUNT, NOT FILE SIZE. output_exists() only checks
#      size >= 100 bytes, which the 4.55 GB truncated file would have passed.
#      gpkg_ok() opens the GPKG and compares its feature count against what we
#      expect - cheap (SQLite reads the count, it does not load geometries) and
#      it actually catches truncation.
# ==============================================================================

wdpa_raw_path <- file.path(build_data_path, "wdpa_raw.gpkg")

#' Feature count of a GPKG layer, or NA if unreadable/corrupt
gpkg_feature_count <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  info <- tryCatch(sf::st_layers(path), error = function(e) NULL)
  if (is.null(info) || length(info$name) == 0) return(NA_integer_)

  n <- suppressWarnings(as.integer(info$features[1]))
  # GDAL reports -1 when the count is not cached; fall back to a COUNT(*),
  # which is still cheap in SQLite.
  if (is.na(n) || n < 0) {
    n <- tryCatch({
      q <- sprintf("SELECT COUNT(*) AS n FROM \"%s\"", info$name[1])
      as.integer(sf::st_read(path, query = q, quiet = TRUE)$n[1])
    }, error = function(e) NA_integer_)
  }
  return(n)
}

#' TRUE only if the GPKG opens AND holds at least min_features features
gpkg_ok <- function(path, min_features) {
  n <- gpkg_feature_count(path)
  !is.na(n) && n >= min_features
}

#' Temp path for an atomic write.
#' The .gpkg extension MUST come last: sf::st_write() and GDAL infer the driver
#' from the extension, so a name ending in ".tmp.<jobid>" fails with
#' "Could not guess driver" - and it would fail in Phase B, i.e. AFTER the
#' expensive load and geometry repair. Verified by test.
tmp_path_for <- function(final_path) {
  file.path(dirname(final_path),
            sprintf("%s.tmp.%s.gpkg",
                    tools::file_path_sans_ext(basename(final_path)),
                    get_slurm_job_id()))
}

#' Publish a verified temp file into its final location, atomically
publish <- function(tmp_path, final_path, min_features, what) {
  n <- gpkg_feature_count(tmp_path)
  if (is.na(n)) {
    unlink(tmp_path)
    stop(sprintf("%s: output is unreadable, discarded: %s", what, tmp_path))
  }
  if (n < min_features) {
    unlink(tmp_path)
    stop(sprintf("%s: only %s features, expected >= %s. Truncated - discarded.",
                 what, format(n, big.mark = ","), format(min_features, big.mark = ",")))
  }
  if (!file.rename(tmp_path, final_path)) {
    stop(sprintf("%s: could not rename %s -> %s", what, tmp_path, final_path))
  }
  log_message(sprintf("%s: published %s features -> %s",
                      what, format(n, big.mark = ","), basename(final_path)))
}

# Clear .tmp files left by earlier killed jobs. They are job-id suffixed, so
# this cannot touch a temp file belonging to a concurrently running job.
stale <- list.files(build_data_path,
                    pattern = "^wdpa_(raw|cleaned)\\.tmp\\..*\\.gpkg$",
                    full.names = TRUE)
stale <- stale[!grepl(sprintf("\\.tmp\\.%s\\.gpkg$", get_slurm_job_id()), stale)]
if (length(stale) > 0) {
  log_message(sprintf("Removing %d stale temp file(s) from earlier runs:", length(stale)))
  for (f in stale) log_message(sprintf("  %s (%.1f GB)", basename(f),
                                       file.info(f)$size / 1024^3))
  unlink(stale)
}

# ==============================================================================
# SOURCE INSPECTION - establishes the expected feature count
# ==============================================================================

log_message(sprintf("Source: %s", wdpa_path))

if (!file.exists(wdpa_path)) {
  stop(sprintf("WDPA file not found: %s", wdpa_path))
}

layers_info <- st_layers(wdpa_path)
log_message(sprintf("Available layers: %s", paste(layers_info$name, collapse = ", ")))
log_message(sprintf("Using layer: %s", layers_info$name[1]))

n_source <- as.integer(layers_info$features[1])
log_message(sprintf("Source features: %s", format(n_source, big.mark = ",")))

# -skipfailures may legitimately drop a few unconvertible geometries, so the
# conversion floor is 99% of source rather than an exact match.
MIN_RAW <- as.integer(floor(0.99 * n_source))

# ==============================================================================
# ALREADY COMPLETE?
# ==============================================================================

if (gpkg_ok(wdpa_clean_path, MIN_RAW * 0.95)) {
  log_message(sprintf("Cleaned WDPA already exists and verifies: %s", wdpa_clean_path))
  log_message(sprintf("Features: %s", format(gpkg_feature_count(wdpa_clean_path),
                                             big.mark = ",")))
  log_message("Delete this file to regenerate.")
  quit(save = "no", status = 0)
}

if (file.exists(wdpa_clean_path)) {
  log_message("WARNING: wdpa_cleaned.gpkg exists but FAILED verification.")
  log_message("         Treating it as incomplete and rebuilding.")
  unlink(wdpa_clean_path)
}

# ==============================================================================
# PHASE A: CONVERT GDB -> GPKG  (expensive; checkpointed to wdpa_raw.gpkg)
# ==============================================================================

if (gpkg_ok(wdpa_raw_path, MIN_RAW)) {
  log_message(sprintf("PHASE A: skipped - %s already verifies (%s features)",
                      basename(wdpa_raw_path),
                      format(gpkg_feature_count(wdpa_raw_path), big.mark = ",")))
} else {
  if (file.exists(wdpa_raw_path)) {
    log_message("PHASE A: existing wdpa_raw.gpkg failed verification, rebuilding.")
    unlink(wdpa_raw_path)
  }

  log_message("PHASE A: converting GDB to GPKG, reprojecting to WGS84...")
  log_message("  The previous attempt exceeded 4h at this step. Tuning applied:")
  log_message("    OGR_SQLITE_SYNCHRONOUS=OFF, 512MB page cache, -gt 65536,")
  log_message("    SPATIAL_INDEX=NO (Phase B reads sequentially and st_write")
  log_message("    rebuilds the index on the final output anyway).")

  # GDAL reads these from the environment, which works across sf versions -
  # unlike the config_options argument, which is not present in all of them.
  Sys.setenv(OGR_SQLITE_SYNCHRONOUS = "OFF")
  Sys.setenv(OGR_SQLITE_CACHE = "512")

  raw_tmp <- tmp_path_for(wdpa_raw_path)
  unlink(raw_tmp)

  phase_a_start <- Sys.time()
  sf::gdal_utils(
    util = "vectortranslate",
    source = wdpa_path,
    destination = raw_tmp,
    options = c(
      "-f", "GPKG",
      "-skipfailures",
      "-nlt", "PROMOTE_TO_MULTI",
      "-t_srs", "EPSG:4326",
      "-lco", "SPATIAL_INDEX=NO",
      "-gt", "65536"
    )
  )
  log_message(sprintf("PHASE A: conversion took %.1f minutes",
                      as.numeric(difftime(Sys.time(), phase_a_start, units = "mins"))))

  publish(raw_tmp, wdpa_raw_path, MIN_RAW, "PHASE A")
}

# ==============================================================================
# PHASE B: VALIDATE AND REPAIR  (memory-heavy; resumable without redoing A)
# ==============================================================================

log_message("PHASE B: loading converted file for geometry validation...")
log_message("  This is the peak-memory step. If the job dies here, rerunning")
log_message("  resumes from this point - Phase A will not be repeated.")

wdpa <- st_read(wdpa_raw_path, quiet = TRUE)
log_message(sprintf("PHASE B: loaded %d features", nrow(wdpa)))

log_message("PHASE B: checking geometry validity...")
sf_use_s2(FALSE)  # Use planar geometry for validation
invalid_idx <- !st_is_valid(wdpa)
n_invalid <- sum(invalid_idx, na.rm = TRUE)

if (n_invalid > 0) {
  log_message(sprintf("Found %d invalid geometries (%.2f%%), attempting repair...",
                      n_invalid, 100 * n_invalid / nrow(wdpa)))

  wdpa <- st_make_valid(wdpa)

  still_invalid <- sum(!st_is_valid(wdpa), na.rm = TRUE)
  if (still_invalid > 0) {
    log_message(sprintf("WARNING: %d geometries still invalid after repair", still_invalid))
  } else {
    log_message("All geometries repaired successfully")
  }
}

# Remove empty geometries
empty_idx <- st_is_empty(wdpa)
n_empty <- sum(empty_idx)
if (n_empty > 0) {
  log_message(sprintf("Removing %d empty geometries", n_empty))
  wdpa <- wdpa[!empty_idx, ]
}

sf_use_s2(TRUE)

# ==============================================================================
# WRITE CLEANED FILE (atomically - never write directly to the target path)
# ==============================================================================

log_message("PHASE B: writing cleaned WDPA file...")
clean_tmp <- tmp_path_for(wdpa_clean_path)
unlink(clean_tmp)

# layer = "wdpa_cleaned" explicitly: st_write would otherwise name the layer
# after the temp FILE (wdpa_cleaned.tmp.<jobid>), and renaming the file on
# publish does not rename the layer inside it.
st_write(wdpa, clean_tmp, layer = "wdpa_cleaned", driver = "GPKG", quiet = TRUE)

# The in-memory row count is the exact expectation for the written file.
publish(clean_tmp, wdpa_clean_path, nrow(wdpa), "PHASE B")

log_message(sprintf("Checkpoint retained for future reruns: %s (%.1f GB)",
                    basename(wdpa_raw_path),
                    file.info(wdpa_raw_path)$size / 1024^3))

# ==============================================================================
# REPORT STATISTICS
# ==============================================================================

log_message("========================================")
log_message("WDPA Preprocessing Complete")
log_message("========================================")

file_size_mb <- file.info(wdpa_clean_path)$size / 1024 / 1024
log_message(sprintf("Output file: %s", wdpa_clean_path))
log_message(sprintf("File size: %.1f MB", file_size_mb))
log_message(sprintf("Total features: %s", format(nrow(wdpa), big.mark = ",")))

# Column summary
log_message("Columns:")
for (col in names(wdpa)) {
  if (col != "geom" && col != "geometry") {
    log_message(sprintf("  - %s", col))
  }
}

# IUCN category distribution
if ("IUCN_CAT" %in% names(wdpa)) {
  iucn_dist <- as.data.table(wdpa)[, .N, by = IUCN_CAT]
  setorder(iucn_dist, -N)
  log_message("IUCN category distribution:")
  for (i in seq_len(min(nrow(iucn_dist), 10))) {
    log_message(sprintf("  %s: %s",
                        iucn_dist$IUCN_CAT[i],
                        format(iucn_dist$N[i], big.mark = ",")))
  }
}

# Report area
area_cols <- grep("area|AREA|Area|REP_AREA|GIS_AREA|REP_M_AREA", names(wdpa), value = TRUE)
if (length(area_cols) > 0) {
  area_col <- area_cols[1]
  total_area <- sum(as.numeric(wdpa[[area_col]]), na.rm = TRUE)
  log_message(sprintf("Total area (%s): %.0f sq km", area_col, total_area))
}

# Timing
elapsed <- difftime(Sys.time(), start_time, units = "mins")
log_message(sprintf("Completed in %.1f minutes", as.numeric(elapsed)))
