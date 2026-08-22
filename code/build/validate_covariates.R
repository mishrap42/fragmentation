# ==============================================================================
# VALIDATE COVARIATES: Visualize all extracted covariate columns
# ==============================================================================
# Creates summary statistics and visualizations for Stage 5 covariate data.
#
# Output: LOGS/validation/covariates/
#   - summary_stats.txt (text summary)
#   - *.png (visualizations for each variable)
#
# Usage: Rscript code/build/validate_covariates.R [tile_id]
#        If tile_id omitted, validates all available tiles
# ==============================================================================

here::i_am("code/build/validate_covariates.R")
source("code/build/BUILD_workspace.R")

library(ggplot2)

# Get optional tile_id argument
args <- commandArgs(trailingOnly = TRUE)
single_tile <- if (length(args) > 0) as.integer(args[1]) else NULL

# Output directory (same as other visualization scripts)
output_dir <- file.path(project_root, "output", "figures", "covariates")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

log_message("========================================")
log_message("COVARIATE VALIDATION")
log_message("========================================")

# ==============================================================================
# LOAD COVARIATE DATA
# ==============================================================================

if (!is.null(single_tile)) {
  log_message(sprintf("Loading covariate data for tile %d...", single_tile))
  cov_file <- get_covariates_filename(single_tile)
  if (!file.exists(cov_file)) {
    stop(sprintf("Covariate file not found: %s", cov_file))
  }
  cov_data <- arrow::read_parquet(cov_file)
  setDT(cov_data)
  output_prefix <- sprintf("tile_%03d_", single_tile)
} else {
  log_message("Loading all covariate files...")
  cov_files <- sapply(1:N_TMF_TILES, get_covariates_filename)
  files_exist <- file.exists(cov_files)

  if (sum(files_exist) == 0) {
    stop("No covariate files found. Run Stage 5 first.")
  }

  log_message(sprintf("Found %d/%d covariate files", sum(files_exist), N_TMF_TILES))

  cov_list <- lapply(cov_files[files_exist], function(f) {
    dt <- arrow::read_parquet(f)
    setDT(dt)
    return(dt)
  })

  cov_data <- rbindlist(cov_list, fill = TRUE)
  rm(cov_list)
  gc()
  output_prefix <- "all_tiles_"
}

log_message(sprintf("Loaded %s cells", format(nrow(cov_data), big.mark = ",")))

# ==============================================================================
# IDENTIFY NUMERIC COLUMNS (skip grid_id)
# ==============================================================================

numeric_cols <- names(cov_data)[sapply(cov_data, is.numeric)]
log_message(sprintf("Numeric columns to visualize: %d", length(numeric_cols)))
log_message(sprintf("  %s", paste(numeric_cols, collapse = ", ")))

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

log_message("Computing summary statistics...")

summary_file <- file.path(output_dir, paste0(output_prefix, "summary_stats.txt"))

sink(summary_file)
cat("COVARIATE VALIDATION SUMMARY\n")
cat("============================\n")
cat(sprintf("Generated: %s\n", Sys.time()))
cat(sprintf("Total cells: %s\n\n", format(nrow(cov_data), big.mark = ",")))

cat("COLUMN SUMMARY:\n")
cat(rep("-", 60), "\n", sep = "")

for (col in numeric_cols) {
  vals <- cov_data[[col]]
  n_na <- sum(is.na(vals))
  n_valid <- sum(!is.na(vals))

  cat(sprintf("\n%s:\n", col))
  cat(sprintf("  N valid: %s (%.1f%%)\n",
              format(n_valid, big.mark = ","),
              100 * n_valid / nrow(cov_data)))
  cat(sprintf("  N missing: %s (%.1f%%)\n",
              format(n_na, big.mark = ","),
              100 * n_na / nrow(cov_data)))

  if (n_valid > 0) {
    cat(sprintf("  Min: %.4f\n", min(vals, na.rm = TRUE)))
    cat(sprintf("  Q1: %.4f\n", quantile(vals, 0.25, na.rm = TRUE)))
    cat(sprintf("  Median: %.4f\n", median(vals, na.rm = TRUE)))
    cat(sprintf("  Mean: %.4f\n", mean(vals, na.rm = TRUE)))
    cat(sprintf("  Q3: %.4f\n", quantile(vals, 0.75, na.rm = TRUE)))
    cat(sprintf("  Max: %.4f\n", max(vals, na.rm = TRUE)))
    cat(sprintf("  SD: %.4f\n", sd(vals, na.rm = TRUE)))
  }
}

sink()
log_message(sprintf("Summary written: %s", basename(summary_file)))

# ==============================================================================
# LOAD GRID COORDINATES FOR MAPPING
# ==============================================================================

log_message("Loading grid coordinates for mapping...")

# Load grid metadata (parquet files, no geometry needed)
if (!is.null(single_tile)) {
  grid_data <- load_grid_for_tile(single_tile, with_geometry = FALSE)
} else {
  grid_files <- sapply(1:N_SUB_TILES, function(s) get_grid_filename(s, "parquet"))
  files_exist <- file.exists(grid_files)
  grid_list <- lapply(grid_files[files_exist], function(f) {
    dt <- arrow::read_parquet(f)
    setDT(dt)
    return(dt[, .(grid_id, centroid_lon, centroid_lat)])
  })
  grid_data <- rbindlist(grid_list, fill = TRUE)
  rm(grid_list)
}

# Merge coordinates with covariate data
cov_data <- merge(cov_data,
                  grid_data[, .(grid_id, centroid_lon, centroid_lat)],
                  by = "grid_id", all.x = TRUE)

log_message(sprintf("Merged coordinates for %s cells",
                    format(sum(!is.na(cov_data$centroid_lon)), big.mark = ",")))

# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

create_histogram <- function(data, col, output_path) {
  vals <- data[[col]]
  vals <- vals[!is.na(vals)]

  if (length(vals) == 0) {
    log_message(sprintf("  Skipping %s (no valid values)", col))
    return(NULL)
  }

  p <- ggplot(data.frame(x = vals), aes(x = x)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
    labs(title = sprintf("Distribution of %s", col),
         subtitle = sprintf("N = %s, Mean = %.2f, Median = %.2f",
                            format(length(vals), big.mark = ","),
                            mean(vals), median(vals)),
         x = col, y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))

  ggsave(output_path, p, width = 8, height = 5, dpi = 150)
}

create_map <- function(data, col, output_path) {
  # Sample if too many points
  plot_data <- data[!is.na(get(col)) & !is.na(centroid_lon)]

  if (nrow(plot_data) == 0) {
    log_message(sprintf("  Skipping map for %s (no valid values)", col))
    return(NULL)
  }

  if (nrow(plot_data) > 100000) {
    set.seed(42)
    plot_data <- plot_data[sample(.N, 100000)]
  }

  # Determine color scale limits (clip outliers)
  vals <- plot_data[[col]]
  q01 <- quantile(vals, 0.01, na.rm = TRUE)
  q99 <- quantile(vals, 0.99, na.rm = TRUE)

  p <- ggplot(plot_data, aes(x = centroid_lon, y = centroid_lat, color = get(col))) +
    geom_point(size = 0.1, alpha = 0.5) +
    scale_color_viridis_c(name = col, limits = c(q01, q99), oob = scales::squish) +
    coord_fixed(ratio = 1.3) +
    labs(title = sprintf("Spatial distribution: %s", col),
         x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right")

  ggsave(output_path, p, width = 12, height = 8, dpi = 150)
}

# ==============================================================================
# CREATE VISUALIZATIONS
# ==============================================================================

log_message("Creating visualizations...")

for (col in numeric_cols) {
  log_message(sprintf("  Processing: %s", col))

  # Histogram
  hist_path <- file.path(output_dir, paste0(output_prefix, "hist_", col, ".png"))
  tryCatch({
    create_histogram(cov_data, col, hist_path)
  }, error = function(e) {
    log_message(sprintf("    ERROR creating histogram: %s", e$message))
  })

  # Map
  map_path <- file.path(output_dir, paste0(output_prefix, "map_", col, ".png"))
  tryCatch({
    create_map(cov_data, col, map_path)
  }, error = function(e) {
    log_message(sprintf("    ERROR creating map: %s", e$message))
  })
}

# ==============================================================================
# CORRELATION MATRIX
# ==============================================================================

log_message("Creating correlation matrix...")

# Select numeric columns with enough valid data
valid_cols <- character(0)
for (col in numeric_cols) {
  if (sum(!is.na(cov_data[[col]])) > 100) {
    valid_cols <- c(valid_cols, col)
  }
}

if (length(valid_cols) >= 2) {
  cor_data <- cov_data[, ..valid_cols]
  cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")

  # Convert to long format for ggplot
  cor_df <- as.data.frame(as.table(cor_matrix))
  names(cor_df) <- c("Var1", "Var2", "Correlation")

  p <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 2) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = "Covariate Correlation Matrix") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold"))

  cor_path <- file.path(output_dir, paste0(output_prefix, "correlation_matrix.png"))
  ggsave(cor_path, p, width = 10, height = 8, dpi = 150)
  log_message(sprintf("Correlation matrix saved: %s", basename(cor_path)))
}

# ==============================================================================
# MISSING DATA SUMMARY
# ==============================================================================

log_message("Creating missing data summary...")

missing_df <- data.frame(
  variable = numeric_cols,
  n_missing = sapply(numeric_cols, function(col) sum(is.na(cov_data[[col]]))),
  pct_missing = sapply(numeric_cols, function(col)
    100 * sum(is.na(cov_data[[col]])) / nrow(cov_data))
)
missing_df <- missing_df[order(-missing_df$pct_missing), ]

p <- ggplot(missing_df, aes(x = reorder(variable, pct_missing), y = pct_missing)) +
  geom_bar(stat = "identity", fill = "coral") +
  geom_text(aes(label = sprintf("%.1f%%", pct_missing)), hjust = -0.1, size = 3) +
  coord_flip() +
  labs(title = "Missing Data by Variable",
       x = "", y = "Percent Missing") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold")) +
  ylim(0, max(missing_df$pct_missing) * 1.15)

missing_path <- file.path(output_dir, paste0(output_prefix, "missing_data.png"))
ggsave(missing_path, p, width = 8, height = 6, dpi = 150)

# ==============================================================================
# SUMMARY
# ==============================================================================

log_message("========================================")
log_message("VALIDATION COMPLETE")
log_message("========================================")
log_message(sprintf("Output directory: %s", output_dir))
log_message(sprintf("Files created: %d",
                    length(list.files(output_dir, pattern = output_prefix))))
log_message("")
log_message("Key files:")
log_message(sprintf("  Summary stats: %ssummary_stats.txt", output_prefix))
log_message(sprintf("  Correlation: %scorrelation_matrix.png", output_prefix))
log_message(sprintf("  Missing data: %smissing_data.png", output_prefix))
