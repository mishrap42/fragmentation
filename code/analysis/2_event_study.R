# ==============================================================================
# 2_event_study.R
# Event study analysis of protected area designation on forest cover
# ==============================================================================

library(data.table)
library(duckdb)
library(DBI)
library(fixest)
library(ggplot2)

# ------------------------------------------------------------------------------
# Setup paths
# ------------------------------------------------------------------------------

if (Sys.getenv("SLURM_SUBMIT_DIR") != "") {
  project_root <- Sys.getenv("SLURM_SUBMIT_DIR")
} else {
  project_root <- here::here()
}

data_dir <- file.path(project_root, "Data/build/final")
temp_dir <- file.path(project_root, "Data/temp")
panel_path <- file.path(data_dir, "TMF_5km_panel.parquet")

cat("Loading data from:", panel_path, "\n")

# ------------------------------------------------------------------------------
# Load only necessary variables for event study
# ------------------------------------------------------------------------------

cat("Connecting to DuckDB and querying parquet...\n")
con <- dbConnect(duckdb())

# Select only the columns needed for event study:
# - Identifiers: grid_id, wdpa_pid, country_iso3
# - Time: year, desig_year
# - Outcomes: TMF land cover fractions
query <- sprintf("
  SELECT
    grid_id,
    area_km2,
    centroid_lon,
    centroid_lat,
    wdpa_pid,
    country_iso3,
    year,
    desig_year,
    is_protected,
    Undisturbed_TMF,
    Degraded_TMF,
    Deforested,
    Regrowth,
    forest_cover
  FROM read_parquet('%s')
", panel_path)

dt <- as.data.table(dbGetQuery(con, query))

dbDisconnect(con, shutdown = TRUE)
cat("Loaded", format(nrow(dt), big.mark = ","), "observations\n")
cat("Unique grid cells:", format(uniqueN(dt$grid_id), big.mark = ","), "\n")
cat("Year range:", min(dt$year), "-", max(dt$year), "\n")

gc()

fwrite(dt[country_iso3 == 'PER',], file.path(temp_dir, "peru_panel_sample.csv"))

# ------------------------------------------------------------------------------
# Create event time variable
# ------------------------------------------------------------------------------

# Event time: years relative to designation
# Negative = before designation, 0 = year of designation, positive = after
dt[, event_time := fifelse(is_protected, year - desig_year, -9999)]
dt[, event_time := fifelse(event_time <= -10 & event_time > -9999, -10, fifelse(event_time >= 15, 15, event_time))]

cat("\nEvent time distribution (protected cells only):\n")
print(dt[!is.na(event_time), .N, by = event_time][order(event_time)])

# ------------------------------------------------------------------------------
# Summary statistics
# ------------------------------------------------------------------------------

cat("\n=== Event Study Data Summary ===\n")
cat("Protected cells:", format(dt[is_protected == TRUE, uniqueN(grid_id)], big.mark = ","), "\n")
cat("Unprotected cells:", format(dt[is_protected == FALSE, uniqueN(grid_id)], big.mark = ","), "\n")
cat("Unique protected areas:", format(dt[!is.na(wdpa_pid), uniqueN(wdpa_pid)], big.mark = ","), "\n")
cat("Countries:", uniqueN(dt$country_iso3), "\n")

cat("\nDesignation year distribution:\n")
print(dt[!is.na(desig_year), .(n_cells = uniqueN(grid_id)), by = desig_year][order(desig_year)])

# ------------------------------------------------------------------------------
# TWFE Model
# ------------------------------------------------------------------------------

twfe.model <- feols(
  sw(Undisturbed_TMF, Degraded_TMF,
     Deforested, Regrowth) ~ i(event_time, ref = c(-1, -9999)) | grid_id + country_iso3^year + year,
  data = dt,
  weights = ~area_km2,
  cluster = ~wdpa_pid,
  lean = TRUE
)

# ------------------------------------------------------------------------------
# Event Study Plot (4-panel: one per outcome)
# ------------------------------------------------------------------------------

cat("\nGenerating event study plots...\n")

# Outcome names and labels
outcome_names <- c("Undisturbed_TMF", "Degraded_TMF", "Deforested", "Regrowth")
outcome_labels <- c("Undisturbed TMF", "Degraded TMF", "Deforested", "Regrowth")

# Extract coefficients for all outcomes
coef_list <- lapply(seq_along(outcome_names), function(i) {
  model_i <- twfe.model[[i]]

  # Get coefficient table
  ct <- coeftable(model_i, keep = "event_time")
  coef_df <- as.data.table(ct)
  setnames(coef_df, c("estimate", "se", "t", "p"))

  # Extract event_time from row names
  coef_df[, event_time := as.integer(gsub("event_time::(-?[0-9]+)", "\\1", rownames(ct)))]

  # Add outcome identifier
 coef_df[, outcome := outcome_labels[i]]

  coef_df
})

# Combine all outcomes
coef_all <- rbindlist(coef_list)

# Add reference period (event_time = -1) for each outcome
ref_rows <- data.table(
  estimate = 0, se = 0, t = NA_real_, p = NA_real_, event_time = -1L,
  outcome = outcome_labels
)
coef_all <- rbind(coef_all, ref_rows)

# Calculate 95% confidence intervals
coef_all[, ci_lower := estimate - 1.96 * se]
coef_all[, ci_upper := estimate + 1.96 * se]

# Sort
setorder(coef_all, outcome, event_time)

# Factor for panel ordering
coef_all[, outcome := factor(outcome, levels = outcome_labels)]

# Define breaks for x-axis
x_breaks <- c(-10, -5, -1, 0, 5, 10, 15)
x_labels <- c("\u2264-10", "-5", "-1", "0", "5", "10", "\u226515")

# Create the 4-panel plot
p_event <- ggplot(coef_all, aes(x = event_time, y = estimate)) +
  # Zero reference line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  # Vertical line at treatment (t = 0)
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "gray60", linewidth = 0.5) +
  # Confidence intervals (no caps)
  geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), color = "steelblue", linewidth = 0.7) +
  # Point estimates
  geom_point(color = "steelblue", size = 2) +
  # Facet by outcome (2x2 grid)
  facet_wrap(~ outcome, ncol = 2, scales = "free_y") +
  # Scales
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  # Labels
  labs(
    x = "Years Relative to Protected Area Designation",
    y = "Effect (pp)",
    title = "Effect of Protected Area Designation on Forest Cover",
    subtitle = "TWFE with grid cell + country\u00D7year FE, SEs clustered at PA level"
  ) +
  # Theme: minimal/void-like with clear axes
  theme(
    # Remove all background elements
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    # Keep axes
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 11),
    # Facet styling
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_blank(),
    # Title styling
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0),
    # Margins
    plot.margin = margin(15, 15, 15, 15),
    panel.spacing = unit(1, "lines")
  )

# Create output directory if it doesn't exist
fig_dir <- file.path(project_root, "output/figures/analysis")
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
  cat("Created directory:", fig_dir, "\n")
}

# Save figure
fig_path <- file.path(fig_dir, "event_study_all_outcomes.png")
ggsave(fig_path, p_event, width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved figure to:", fig_path, "\n")
