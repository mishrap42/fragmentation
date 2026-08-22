# ==============================================================================
# 3_matching.R
# MatchIt nearest-neighbor 1:1 matching for protected area analysis
# Then re-run TWFE event study on matched sample
# ==============================================================================

library(data.table)
library(duckdb)
library(DBI)
# Note: MatchIt, fixest, cobalt loaded AFTER DuckDB query to avoid dependency conflicts

# ------------------------------------------------------------------------------
# Setup paths
# ------------------------------------------------------------------------------

if (Sys.getenv("SLURM_SUBMIT_DIR") != "") {
  project_root <- Sys.getenv("SLURM_SUBMIT_DIR")
} else {
  project_root <- here::here()
}

data_dir <- file.path(project_root, "Data/build/final")
intermediate_dir <- file.path(project_root, "Data/build/intermediate")
fig_dir <- file.path(project_root, "output/figures/analysis")

panel_path <- file.path(data_dir, "TMF_5km_panel.parquet")

# Create directories if needed
if (!dir.exists(intermediate_dir)) dir.create(intermediate_dir, recursive = TRUE)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

cat("Loading matching data from:", panel_path, "\n")

# ==============================================================================
# PART 1: Load matching covariates (year = 1990 cross-section)
# ==============================================================================

cat("\nPART 1: Loading 1990 cross-section for matching...\n")

con <- dbConnect(duckdb())

# Test query: same pattern as 1_balance_PA.R but year = 1990
cat("Testing: SELECT * WHERE year = 1990...\n")
test_query <- sprintf("SELECT * FROM read_parquet('%s') WHERE year = 1990", panel_path)
dt_match <- as.data.table(dbGetQuery(con, test_query))
cat("Test passed! Loaded", nrow(dt_match), "rows\n")

# Now select only the columns we need
dt_match <- dt_match[, .(
  grid_id, country_iso3, is_protected, desig_year,
  pop_density_1990, travel_time_cities, dist_to_city_km,
  signed_dist_to_frontier_km, elevation_m,
  maize_yield = yield_maize,
  Undisturbed_TMF, Degraded_TMF, Deforested, forest_cover
)]

dbDisconnect(con, shutdown = TRUE)

cat("Loaded", format(nrow(dt_match), big.mark = ","), "cells for year 1990\n")

gc()

# Load remaining libraries now that DuckDB is done
library(MatchIt)
library(fixest)
library(ggplot2)
library(cobalt)

# ==============================================================================
# PART 2: Prepare treatment and exact-match variables
# ==============================================================================

cat("\nPART 2: Preparing treatment and exact-match variables...\n")

# Treatment: protected with designation year in panel period (1991+)
# Controls: never protected
dt_match[, treated := is_protected == TRUE & !is.na(desig_year) & desig_year >= 1991]
dt_match[, control := is_protected == FALSE]

cat("Treated cells (protected, desig >= 1991):", format(sum(dt_match$treated), big.mark = ","), "\n")
cat("Control cells (unprotected):", format(sum(dt_match$control), big.mark = ","), "\n")

# Keep only treated and control cells (already filtered in SQL, but verify)
dt_match <- dt_match[treated | control]
cat("Cells for matching:", format(nrow(dt_match), big.mark = ","), "\n")
cat("  Treated:", sum(dt_match$treated), "\n")
cat("  Control:", sum(!dt_match$treated), "\n")

# Create 1990 land use code (categorical based on dominant class)
dt_match[, luc_1990 := fcase(
  Undisturbed_TMF >= 0.5, "undisturbed",
  Degraded_TMF >= 0.5, "degraded",
  Deforested >= 0.5, "deforested",
  forest_cover >= 0.5, "mixed_forest",
  default = "other"
)]

cat("\n1990 Land Use Code distribution:\n")
print(dt_match[, .N, by = luc_1990][order(-N)])

cat("\nLand use by treatment status:\n")
print(dt_match[, .(N = .N), by = .(treated, luc_1990)][order(treated, -N)])

# ==============================================================================
# PART 3: Standardize continuous matching variables
# ==============================================================================

cat("\nPART 3: Standardizing continuous matching variables...\n")

# Define continuous variables to standardize
std_vars <- c("pop_density_1990", "travel_time_cities", "dist_to_city_km",
              "signed_dist_to_frontier_km", "elevation_m", "maize_yield")

# Check for missing values
cat("\nMissing values in matching variables:\n")
for (v in std_vars) {
  n_missing <- sum(is.na(dt_match[[v]]))
  pct_missing <- 100 * n_missing / nrow(dt_match)
  cat(sprintf("  %s: %s (%.1f%%)\n", v, format(n_missing, big.mark = ","), pct_missing))
}

# Remove rows with any missing matching variables
dt_match_complete <- dt_match[complete.cases(dt_match[, ..std_vars])]
cat("\nCells with complete matching data:", format(nrow(dt_match_complete), big.mark = ","), "\n")
cat("Treated:", format(sum(dt_match_complete$treated), big.mark = ","), "\n")
cat("Control:", format(sum(!dt_match_complete$treated), big.mark = ","), "\n")

# Standardize variables (z-score)
for (v in std_vars) {
  z_var <- paste0(v, "_z")
  dt_match_complete[, (z_var) := scale(get(v))]
}

cat("\nStandardized variable means and SDs:\n")
for (v in std_vars) {
  z_var <- paste0(v, "_z")
  cat(sprintf("  %s: mean=%.3f, sd=%.3f\n", z_var,
              mean(dt_match_complete[[z_var]], na.rm = TRUE),
              sd(dt_match_complete[[z_var]], na.rm = TRUE)))
}

# Convert to data.frame for MatchIt (required)
df_match <- as.data.frame(dt_match_complete)

# ==============================================================================
# PART 4: Run MatchIt (nearest neighbor, Euclidean, exact on country + LUC)
# ==============================================================================

cat("\nPART 4: Running MatchIt...\n")
cat("Method: nearest neighbor 1:1, Euclidean distance, exact on country + LUC\n")

# Timing
start_time <- Sys.time()

m.out <- matchit(
  treated ~ pop_density_1990_z + travel_time_cities_z + dist_to_city_km_z +
            signed_dist_to_frontier_km_z + elevation_m_z + maize_yield_z,
  data = df_match,
  method = "nearest",
  distance = "euclidean",
  exact = ~ country_iso3 + luc_1990,
  ratio = 1,
  replace = TRUE
)

end_time <- Sys.time()
cat("Matching completed in:", round(difftime(end_time, start_time, units = "mins"), 1), "minutes\n")

# Print summary
cat("\nMatchIt Summary:\n")
print(summary(m.out))

# ==============================================================================
# PART 5: Save matched grid_ids to intermediate parquet
# ==============================================================================

cat("\nPART 5: Saving matched grid_ids...\n")

# Extract matched data
matched_data <- match.data(m.out)
matched_dt <- as.data.table(matched_data)

cat("Matched sample size:", format(nrow(matched_dt), big.mark = ","), "\n")
cat("Matched treated:", format(sum(matched_dt$treated), big.mark = ","), "\n")
cat("Matched controls:", format(sum(!matched_dt$treated), big.mark = ","), "\n")

# Save matched grid_ids (using fwrite - no arrow dependency)
matched_ids <- matched_dt[, .(grid_id)]
fwrite(matched_ids, file.path(intermediate_dir, "matched_grid_ids.csv"))
cat("Saved matched grid_ids to:", file.path(intermediate_dir, "matched_grid_ids.csv"), "\n")

# Save full matched data with weights for analysis
cat("Columns in matched data:", paste(names(matched_dt), collapse = ", "), "\n")
matched_pairs <- matched_dt[, .(grid_id, treated, weights)]
fwrite(matched_pairs, file.path(intermediate_dir, "matched_pairs.csv"))
cat("Saved matched pairs to:", file.path(intermediate_dir, "matched_pairs.csv"), "\n")

# ==============================================================================
# PART 6: Generate balance diagnostics (Love plot)
# ==============================================================================

cat("\nPART 6: Generating balance diagnostics...\n")

# Get balance table and filter to continuous vars only
bt <- bal.tab(m.out, stats = c("mean.diffs"), un = TRUE)

# Extract SMDs for continuous matching vars only
cont_vars <- c("pop_density_1990_z", "travel_time_cities_z", "dist_to_city_km_z",
               "signed_dist_to_frontier_km_z", "elevation_m_z", "maize_yield_z")
var_labels <- c("Population Density (1990)", "Travel Time to Cities",
                "Distance to City (km)", "Signed Dist. to Frontier (km)",
                "Elevation (m)", "Maize Yield")

bal_df <- as.data.frame(bt$Balance)
bal_df$var <- rownames(bal_df)
bal_df <- bal_df[bal_df$var %in% cont_vars, ]

# Reshape for plotting
bal_long <- data.table(
  variable = factor(var_labels, levels = rev(var_labels)),
  Unadjusted = abs(bal_df$Diff.Un),
  Adjusted = abs(bal_df$Diff.Adj)
)
bal_melt <- melt(bal_long, id.vars = "variable", variable.name = "Sample",
                 value.name = "SMD")

# Create Love plot manually
p_love <- ggplot(bal_melt, aes(x = SMD, y = variable, color = Sample, shape = Sample)) +

  geom_point(size = 3) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Unadjusted" = "coral", "Adjusted" = "steelblue")) +
  scale_shape_manual(values = c("Unadjusted" = 17, "Adjusted" = 16)) +
  labs(
    x = "Absolute Standardized Mean Difference",
    y = NULL,
    title = "Covariate Balance: Before vs After Matching"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    axis.line = element_line(color = "black", linewidth = 0.3),
    plot.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "matching_balance.png"), p_love,
       width = 8, height = 6, dpi = 300, bg = "white")
cat("Saved balance plot to:", file.path(fig_dir, "matching_balance.png"), "\n")

# Print balance table
cat("\nBalance Table:\n")
print(bal.tab(m.out, stats = c("mean.diffs", "variance.ratios")))

# ==============================================================================
# PART 7: Load full panel for matched cells only
# ==============================================================================

cat("\nPART 7: Loading full panel for matched cells...\n")

gc()

# Get list of matched grid_ids
matched_grid_ids <- matched_dt$grid_id

con <- dbConnect(duckdb())

# Register the list of IDs as a DuckDB table
dbExecute(con, "CREATE OR REPLACE TABLE matched_ids (grid_id VARCHAR)")
dbAppendTable(con, "matched_ids", data.frame(grid_id = matched_grid_ids))

# Query panel filtering to matched IDs
query_panel <- sprintf("
  SELECT
    p.grid_id,
    p.area_km2,
    p.wdpa_pid,
    p.country_iso3,
    p.year,
    p.desig_year,
    p.is_protected,
    p.Undisturbed_TMF,
    p.Degraded_TMF,
    p.Deforested,
    p.Regrowth,
    p.forest_cover
  FROM read_parquet('%s') p
  INNER JOIN matched_ids m ON p.grid_id = m.grid_id
", panel_path)

dt_panel <- as.data.table(dbGetQuery(con, query_panel))
dbDisconnect(con, shutdown = TRUE)

cat("Loaded", format(nrow(dt_panel), big.mark = ","), "observations for matched cells\n")
cat("Unique matched cells:", format(uniqueN(dt_panel$grid_id), big.mark = ","), "\n")

gc()

# ==============================================================================
# PART 8: Run TWFE event study on matched sample
# ==============================================================================

cat("\nPART 8: Running TWFE event study on matched sample...\n")

# Create event time variable
dt_panel[, event_time := fifelse(is_protected, year - desig_year, -9999)]
dt_panel[, event_time := fifelse(event_time <= -10 & event_time > -9999, -10,
                                  fifelse(event_time >= 15, 15, event_time))]

cat("\nEvent time distribution (matched protected cells):\n")
print(dt_panel[is_protected == TRUE, .N, by = event_time][order(event_time)])

# Run TWFE model
cat("\nEstimating TWFE model...\n")

twfe.model <- feols(
  sw(Undisturbed_TMF, Degraded_TMF, Deforested, Regrowth) ~
    i(event_time, ref = c(-1, -9999)) | grid_id + country_iso3^year + year,
  data = dt_panel,
  weights = ~area_km2,
  cluster = ~wdpa_pid,
  lean = TRUE
)

cat("TWFE estimation complete.\n")

# ==============================================================================
# PART 9: Generate 4-panel event study figure
# ==============================================================================

cat("\nPART 9: Generating event study figure...\n")

outcome_names <- c("Undisturbed_TMF", "Degraded_TMF", "Deforested", "Regrowth")
outcome_labels <- c("Undisturbed TMF", "Degraded TMF", "Deforested", "Regrowth")

# Extract coefficients for all outcomes
coef_list <- lapply(seq_along(outcome_names), function(i) {
  model_i <- twfe.model[[i]]

  ct <- coeftable(model_i, keep = "event_time")
  coef_df <- as.data.table(ct)
  setnames(coef_df, c("estimate", "se", "t", "p"))

  coef_df[, event_time := as.integer(gsub("event_time::(-?[0-9]+)", "\\1", rownames(ct)))]
  coef_df[, outcome := outcome_labels[i]]

  coef_df
})

coef_all <- rbindlist(coef_list)

# Add reference period
ref_rows <- data.table(
  estimate = 0, se = 0, t = NA_real_, p = NA_real_, event_time = -1L,
  outcome = outcome_labels
)
coef_all <- rbind(coef_all, ref_rows)

# Calculate 95% CIs
coef_all[, ci_lower := estimate - 1.96 * se]
coef_all[, ci_upper := estimate + 1.96 * se]

setorder(coef_all, outcome, event_time)
coef_all[, outcome := factor(outcome, levels = outcome_labels)]

# X-axis breaks
x_breaks <- c(-10, -5, -1, 0, 5, 10, 15)
x_labels <- c("\u2264-10", "-5", "-1", "0", "5", "10", "\u226515")

# Create plot
p_event <- ggplot(coef_all, aes(x = event_time, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "gray60", linewidth = 0.5) +
  geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), color = "steelblue", linewidth = 0.7) +
  geom_point(color = "steelblue", size = 2) +
  facet_wrap(~ outcome, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  labs(
    x = "Years Relative to Protected Area Designation",
    y = "Effect (pp)",
    title = "Effect of Protected Area Designation on Forest Cover (Matched Sample)",
    subtitle = "TWFE with grid cell + country\u00D7year FE, SEs clustered at PA level"
  ) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 11),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0),
    plot.margin = margin(15, 15, 15, 15),
    panel.spacing = unit(1, "lines")
  )

ggsave(file.path(fig_dir, "event_study_matched.png"), p_event,
       width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved event study figure to:", file.path(fig_dir, "event_study_matched.png"), "\n")

# ==============================================================================
# Summary
# ==============================================================================

cat("\n")
cat("=======================================================\n")
cat("MATCHING AND EVENT STUDY COMPLETE\n")
cat("=======================================================\n")
cat("\nMatched sample:\n")
cat("  Treated cells:", format(sum(matched_dt$treated), big.mark = ","), "\n")
cat("  Control cells:", format(sum(!matched_dt$treated), big.mark = ","), "\n")
cat("  Total cell-years:", format(nrow(dt_panel), big.mark = ","), "\n")
cat("\nOutput files:\n")
cat("  ", file.path(intermediate_dir, "matched_grid_ids.parquet"), "\n")
cat("  ", file.path(intermediate_dir, "matched_pairs.parquet"), "\n")
cat("  ", file.path(fig_dir, "matching_balance.png"), "\n")
cat("  ", file.path(fig_dir, "event_study_matched.png"), "\n")
cat("=======================================================\n")
