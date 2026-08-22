# ==============================================================================
# DIAGNOSE PIPELINE: Missing files + data quality + blank tile analysis
# ==============================================================================
# Consolidated diagnostic that checks:
#   1. Missing output files at each stage
#   2. Data quality issues (empty grids, bad class distributions)
#   3. Blank tiles (no undisturbed forest when expected)
#
# Output: LOGS/rerun/stage*_missing.txt (job IDs to rerun)
#
# Usage: sbatch code/bash/diagnose_pipeline.sh
#        bash code/bash/rerun_missing.sh <stage>
# ==============================================================================

.libPaths(c('/dartfs-hpc/rc/lab/M/MishraP/Rlibs/4.4.0/', .libPaths()))

library(here)
here::i_am("code/build/diagnose_pipeline.R")
source("code/build/BUILD_workspace.R")

cat(rep("=", 71), "\n", sep = "")
cat("PIPELINE DIAGNOSTIC\n")
cat(rep("=", 71), "\n", sep = "")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("Project: %s\n\n", project_root))

# ==============================================================================
# STAGE 0: GRID FILES (missing + quality + geometry issues)
# ==============================================================================

cat("STAGE 0: Grid Creation\n")
cat(rep("-", 40), "\n", sep = "")

stage0_missing <- c()
stage0_empty <- c()
stage0_low_cells <- c()  # Potential GEOMETRYCOLLECTION drop victims
grid_cell_ratios <- list()

for (sid in 1:N_SUB_TILES) {
  grid_pq <- get_grid_filename(sid, "parquet")
  grid_gpkg <- get_grid_filename(sid, "gpkg")

  if (!file.exists(grid_pq) && !file.exists(grid_gpkg)) {
    stage0_missing <- c(stage0_missing, sid)
  } else {
    grid_file <- if (file.exists(grid_pq)) grid_pq else grid_gpkg
    tryCatch({
      if (grepl("\\.parquet$", grid_file)) {
        g <- arrow::read_parquet(grid_file)
      } else {
        g <- st_read(grid_file, quiet = TRUE)
      }

      n_cells <- nrow(g)
      info <- get_sub_tile_info(sid)

      # Calculate expected cells based on tile area (5km grid = 25km² per cell)
      lat_mid <- (info$lat_south + info$lat_north) / 2
      lon_width <- (info$lon_east - info$lon_west) * cos(lat_mid * pi / 180) * 111
      lat_height <- (info$lat_north - info$lat_south) * 111
      tile_area_km2 <- lon_width * lat_height
      expected_cells <- as.integer(tile_area_km2 / 25)

      cell_ratio <- if (expected_cells > 0) n_cells / expected_cells else 0

      if (n_cells == 0) {
        stage0_empty <- c(stage0_empty, sid)
      } else if (cell_ratio < 0.5 && abs(lat_mid) < 30 && expected_cells > 50) {
        # Less than 50% of expected cells in tropical region = GEOMETRYCOLLECTION issue
        stage0_low_cells <- c(stage0_low_cells, sid)
        grid_cell_ratios[[as.character(sid)]] <- list(
          actual = n_cells,
          expected = expected_cells,
          ratio = cell_ratio,
          lat = lat_mid,
          lon = (info$lon_west + info$lon_east) / 2
        )
      }
      rm(g)
    }, error = function(e) {
      stage0_missing <- c(stage0_missing, sid)
    })
  }
}

stage0_issues <- unique(c(stage0_missing, stage0_empty, stage0_low_cells))

cat(sprintf("  Total sub-tiles: %d\n", N_SUB_TILES))
cat(sprintf("  Missing files: %d\n", length(stage0_missing)))
cat(sprintf("  Empty grids: %d\n", length(stage0_empty)))
cat(sprintf("  Low cell count (<50%% expected): %d\n", length(stage0_low_cells)))

# Show details for low cell count tiles (GEOMETRYCOLLECTION victims)
if (length(stage0_low_cells) > 0) {
  cat("\n  LOW CELL COUNT TILES (potential GEOMETRYCOLLECTION issue):\n")
  for (sid in stage0_low_cells[1:min(10, length(stage0_low_cells))]) {
    r <- grid_cell_ratios[[as.character(sid)]]
    cat(sprintf("    Sub-tile %d: %d cells (expected %d, ratio=%.1f%%) at lat=%.1f, lon=%.1f\n",
                sid, r$actual, r$expected, r$ratio * 100, r$lat, r$lon))
  }
  if (length(stage0_low_cells) > 10) {
    cat(sprintf("    ... and %d more\n", length(stage0_low_cells) - 10))
  }
}

# ==============================================================================
# STAGE 0c: YIELD EXTRACTION (optional)
# ==============================================================================

cat("\nSTAGE 0c: Yield Extraction\n")
cat(rep("-", 40), "\n", sep = "")

yield_files <- sapply(1:N_SUB_TILES, get_yields_filename)
yield_exists <- file.exists(yield_files)
yield_missing <- which(!yield_exists)

cat(sprintf("  Total: %d | Found: %d | Missing: %d\n",
            N_SUB_TILES, sum(yield_exists), length(yield_missing)))

# ==============================================================================
# STAGE 1: TMF EXTRACTION
# ==============================================================================

cat("\nSTAGE 1: TMF Extraction\n")
cat(rep("-", 40), "\n", sep = "")

stage1_issues <- c()

for (tile_id in 1:N_TMF_TILES) {
  for (year in TMF_YEARS) {
    tmf_file <- get_tmf_filename(tile_id, year)
    if (!file.exists(tmf_file)) {
      job_id <- year_tile_to_job(year, tile_id)
      stage1_issues <- c(stage1_issues, job_id)
    }
  }
}

n_total_tmf <- N_TMF_TILES * N_TMF_YEARS
cat(sprintf("  Total tile-years: %d\n", n_total_tmf))
cat(sprintf("  Missing: %d\n", length(stage1_issues)))

# ==============================================================================
# STAGE 2a: CONSOLIDATED TILES (missing + quality)
# ==============================================================================

cat("\nSTAGE 2a: Tile Consolidation\n")
cat(rep("-", 40), "\n", sep = "")

stage2a_issues <- c()
blank_tiles <- c()

for (tile_id in 1:N_TMF_TILES) {
  cons_file <- get_consolidated_tile_filename(tile_id)

  if (!file.exists(cons_file)) {
    stage2a_issues <- c(stage2a_issues, tile_id)
  } else {
    tryCatch({
      cons <- arrow::read_parquet(cons_file)
      setDT(cons)

      if (nrow(cons) == 0) {
        stage2a_issues <- c(stage2a_issues, tile_id)
        blank_tiles <- c(blank_tiles, tile_id)
      } else if ("tmf_class" %in% names(cons) && "fraction" %in% names(cons)) {
        total_frac <- sum(cons$fraction, na.rm = TRUE)
        if (total_frac > 0) {
          und_frac <- sum(cons[tmf_class == "Undisturbed_TMF"]$fraction, na.rm = TRUE)
          unknown_frac <- sum(cons[grepl("^Unknown", tmf_class)]$fraction, na.rm = TRUE)

          # No undisturbed OR >99% unknown in tropical band = blank tile
          if (und_frac == 0 || unknown_frac / total_frac > 0.99) {
            tile_info <- TMF_TILE_INDEX[tile_id]
            if (nrow(tile_info) > 0) {
              lat_center <- (tile_info$lat_south + tile_info$lat_north) / 2
              if (abs(lat_center) < 30) {
                stage2a_issues <- c(stage2a_issues, tile_id)
                blank_tiles <- c(blank_tiles, tile_id)
              }
            }
          }
        }
      }
      rm(cons)
    }, error = function(e) {
      stage2a_issues <- c(stage2a_issues, tile_id)
    })
  }
}

stage2a_issues <- unique(stage2a_issues)
blank_tiles <- unique(blank_tiles)

cat(sprintf("  Total tiles: %d\n", N_TMF_TILES))
cat(sprintf("  Issues (missing/empty/blank): %d\n", length(stage2a_issues)))
cat(sprintf("  Blank tiles (no undisturbed): %d\n", length(blank_tiles)))

if (length(blank_tiles) > 0 && length(blank_tiles) <= 20) {
  cat(sprintf("  Blank tile IDs: %s\n", paste(blank_tiles, collapse = ", ")))
}

# ==============================================================================
# STAGE 3: INTERIOR CLASSIFICATION
# ==============================================================================

cat("\nSTAGE 3: Interior Classification\n")
cat(rep("-", 40), "\n", sep = "")

interior_files <- sapply(1:N_TMF_TILES, get_interior_filename)
stage3_issues <- which(!file.exists(interior_files))

cat(sprintf("  Total: %d | Missing: %d\n", N_TMF_TILES, length(stage3_issues)))

# ==============================================================================
# STAGE 4: FRONTIER CALCULATION
# ==============================================================================

cat("\nSTAGE 4: Frontier Calculation\n")
cat(rep("-", 40), "\n", sep = "")

frontier_files <- sapply(1:N_TMF_TILES, get_frontier_filename)
stage4_issues <- which(!file.exists(frontier_files))

cat(sprintf("  Total: %d | Missing: %d\n", N_TMF_TILES, length(stage4_issues)))

# ==============================================================================
# STAGE 5: COVARIATE EXTRACTION
# ==============================================================================

cat("\nSTAGE 5: Covariate Extraction\n")
cat(rep("-", 40), "\n", sep = "")

covariate_files <- sapply(1:N_TMF_TILES, get_covariates_filename)
stage5_issues <- which(!file.exists(covariate_files))

cat(sprintf("  Total: %d | Missing: %d\n", N_TMF_TILES, length(stage5_issues)))

# ==============================================================================
# STAGE 6: FINAL ASSEMBLY
# ==============================================================================

cat("\nSTAGE 6: Final Assembly\n")
cat(rep("-", 40), "\n", sep = "")

final_file <- file.path(final_output_path, "TMF_5km_panel.parquet")
final_exists <- file.exists(final_file)

cat(sprintf("  Status: %s\n", ifelse(final_exists, "COMPLETE", "MISSING")))
if (final_exists) {
  file_size <- file.info(final_file)$size / 1024 / 1024 / 1024
  cat(sprintf("  Size: %.2f GB\n", file_size))
}

# ==============================================================================
# DEPENDENCY CHECK
# ==============================================================================

cat("\n", rep("=", 71), "\n", sep = "")
cat("DEPENDENCY ANALYSIS\n")
cat(rep("=", 71), "\n", sep = "")

# Which TMF tiles have missing grids?
tmf_tiles_blocked <- c()
for (tile_id in 1:N_TMF_TILES) {
  sub_ids <- SUB_TILE_INDEX[tmf_tile_id == tile_id]$sub_tile_id
  if (any(sub_ids %in% stage0_issues)) {
    tmf_tiles_blocked <- c(tmf_tiles_blocked, tile_id)
  }
}

if (length(tmf_tiles_blocked) > 0) {
  cat(sprintf("TMF tiles blocked by Stage 0 issues: %d\n", length(tmf_tiles_blocked)))
} else {
  cat("All TMF tiles have complete grid coverage.\n")
}

# ==============================================================================
# WRITE RERUN FILES
# ==============================================================================

cat("\n", rep("=", 71), "\n", sep = "")
cat("RERUN JOB LISTS\n")
cat(rep("=", 71), "\n", sep = "")

rerun_dir <- file.path(project_root, "LOGS", "rerun")
if (!dir.exists(rerun_dir)) dir.create(rerun_dir, recursive = TRUE)

write_jobs <- function(stage, ids, label) {
  if (length(ids) > 0) {
    ids <- sort(unique(ids))
    f <- file.path(rerun_dir, sprintf("stage%s_missing.txt", stage))
    writeLines(as.character(ids), f)
    cat(sprintf("Stage %s: %d %s -> %s\n", stage, length(ids), label, basename(f)))
  }
}

write_jobs("0", stage0_issues, "jobs")
if (sum(yield_exists) > 0 && length(yield_missing) > 0) {
  write_jobs("0c", yield_missing, "jobs")
}
write_jobs("1", stage1_issues, "jobs")
write_jobs("2a", stage2a_issues, "jobs")
write_jobs("3", stage3_issues, "jobs")
write_jobs("4", stage4_issues, "jobs")
write_jobs("5", stage5_issues, "jobs")

if (!final_exists) {
  f <- file.path(rerun_dir, "stage6_missing.txt")
  writeLines("1", f)
  cat(sprintf("Stage 6: 1 job -> %s\n", basename(f)))
}

# Write detailed cell ratio report for GEOMETRYCOLLECTION tracking
if (length(grid_cell_ratios) > 0) {
  ratio_df <- data.table(
    sub_tile_id = as.integer(names(grid_cell_ratios)),
    actual_cells = sapply(grid_cell_ratios, function(x) x$actual),
    expected_cells = sapply(grid_cell_ratios, function(x) x$expected),
    cell_ratio = sapply(grid_cell_ratios, function(x) x$ratio),
    lat = sapply(grid_cell_ratios, function(x) x$lat),
    lon = sapply(grid_cell_ratios, function(x) x$lon)
  )
  f <- file.path(rerun_dir, "stage0_low_cell_ratio.csv")
  fwrite(ratio_df, f)
  cat(sprintf("\nGEOMETRYCOLLECTION tracking: %d tiles -> %s\n", nrow(ratio_df), basename(f)))
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n", rep("=", 71), "\n", sep = "")
cat("SUMMARY\n")
cat(rep("=", 71), "\n", sep = "")

total_issues <- length(stage0_issues) + length(stage1_issues) +
                length(stage2a_issues) + length(stage3_issues) +
                length(stage4_issues) + length(stage5_issues) +
                as.integer(!final_exists)

if (total_issues == 0 && length(stage0_low_cells) == 0) {
  cat("\nAll stages complete. No geometry issues detected.\n")
} else {
  cat(sprintf("\nTotal issues: %d\n", total_issues))
  if (length(stage0_low_cells) > 0) {
    cat(sprintf("GEOMETRYCOLLECTION fix needed: %d sub-tiles with <50%% expected cells\n",
                length(stage0_low_cells)))
    cat("  These tiles likely lost cells due to geometry type dropping.\n")
    cat("  Check: LOGS/rerun/stage0_low_cell_ratio.csv\n")
  }
  cat("\nTo rerun: bash code/bash/rerun_missing.sh <stage>\n")
}

cat("\n")
