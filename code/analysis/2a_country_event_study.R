# ==============================================================================
# 2a_country_event_study.R
# Event study analysis of protected area designation on forest cover: parallel country-specific execution
# ==============================================================================

library(data.table)
library(duckdb)
library(DBI)
library(fixest)
library(ggplot2)
library(patchwork)

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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Task ID required.")
}
task_id <- as.integer(args[1])


# ------------------------------------------------------------------------------
# Get list of countries with treated cells
# ------------------------------------------------------------------------------

cat("Getting country list...\n")

con <- dbConnect(duckdb())

# Get countries that have treated cells (protected with desig_year >= 1991)
country_query <- sprintf("
  SELECT DISTINCT country_iso3
  FROM read_parquet('%s')
  WHERE year = 1990
    AND is_protected = true
    AND desig_year >= 1991
  ORDER BY country_iso3
", panel_path)

countries <- dbGetQuery(con, country_query)$country_iso3
dbDisconnect(con, shutdown = TRUE)

cat("Found", length(countries), "countries with treated cells\n")

# Check if task_id is valid
if (task_id > length(countries)) {
  cat("Task ID", task_id, "exceeds country count", length(countries), "- skipping\n")
  quit(status = 0)
}

country <- countries[task_id]
cat("Processing country:", country, "(task", task_id, "of", length(countries), ")\n\n")

# ------------------------------------------------------------------------------
# Load data for this country only
# ------------------------------------------------------------------------------

cat("Loading data for country:", country, "...\n")
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
  WHERE country_iso3 = '%s'
", panel_path, country)

dt <- as.data.table(dbGetQuery(con, query))

dbDisconnect(con, shutdown = TRUE)
cat("Loaded", format(nrow(dt), big.mark = ","), "observations\n")
cat("Unique grid cells:", format(uniqueN(dt$grid_id), big.mark = ","), "\n")
cat("Year range:", min(dt$year), "-", max(dt$year), "\n")

gc()

# ------------------------------------------------------------------------------
# Create event time variable
# ------------------------------------------------------------------------------

dt <- dt[(is.na(desig_year) | desig_year == 0) | (is_protected == TRUE & desig_year >= 1991),]

# Event time: years relative to designation
# Negative = before designation, 0 = year of designation, positive = after
dt[, event_time := fifelse(is_protected, year - desig_year, -9999)]

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
     Deforested, Regrowth) ~ i(event_time, ref = c(-1, -9999)) | grid_id + year,
  data = dt,
  weights = ~area_km2,
  cluster = ~wdpa_pid,
  lean = TRUE
)

# ------------------------------------------------------------------------------
# Event Study Plot (single outcome with obs counts)
# ------------------------------------------------------------------------------

cat("\nGenerating event study plots...\n")

# Extract coefficients for Undisturbed_TMF (primary outcome)
ct <- coeftable(twfe.model[[1]], keep = "event_time")
coef_df <- as.data.table(ct)
setnames(coef_df, c("estimate", "se", "t", "p"))
coef_df[, event_time := as.integer(gsub("event_time::(-?[0-9]+)", "\\1", rownames(ct)))]

# Add reference period
ref_row <- data.table(estimate = 0, se = 0, t = NA_real_, p = NA_real_, event_time = -1L)
coef_df <- rbind(coef_df, ref_row)

# Calculate 95% CIs
coef_df[, ci_lower := estimate - 1.96 * se]
coef_df[, ci_upper := estimate + 1.96 * se]
setorder(coef_df, event_time)

# Get observation counts by event_time (treated obs only)
obs_counts <- dt[event_time != -9999, .N, by = event_time][order(event_time)]

# ------------------------------------------------------------------------------
# Find "preferred" region based on 50% drop from t=-1
# ------------------------------------------------------------------------------

n_ref <- obs_counts[event_time == -1, N]
if (length(n_ref) == 0) n_ref <- max(obs_counts$N)
threshold <- 0.5 * n_ref

# Find left bound
left_times <- obs_counts[event_time < -1][order(-event_time)]
left_bound <- min(obs_counts$event_time)
for (i in seq_len(nrow(left_times))) {
  if (left_times[i, N] < threshold) {
    if (i == 1) {
      left_bound <- -1
    } else {
      left_bound <- left_times[i - 1, event_time]
    }
    break
  }
}

# Find right bound
right_times <- obs_counts[event_time > -1][order(event_time)]
right_bound <- max(obs_counts$event_time)
for (i in seq_len(nrow(right_times))) {
  if (right_times[i, N] < threshold) {
    if (i == 1) {
      right_bound <- -1
    } else {
      right_bound <- right_times[i - 1, event_time]
    }
    break
  }
}

cat("Preferred region: [", left_bound, ", ", right_bound, "]\n")

# X-axis settings
x_breaks <- sort(unique(coef_df$event_time))
x_min <- min(x_breaks)
x_max <- max(x_breaks)

# Top panel: Event study coefficients
p_coef <- ggplot(coef_df, aes(x = event_time, y = estimate)) +
  annotate("rect", xmin = left_bound - 0.5, xmax = right_bound + 0.5,
           ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "gray60", linewidth = 0.5) +
  geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 2.5) +
  scale_x_continuous(breaks = x_breaks, limits = c(x_min - 0.5, x_max + 0.5)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  labs(x = NULL, y = "Effect on Undisturbed TMF (pp)") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 11),
    plot.margin = margin(10, 10, 0, 10)
  )

# Bottom panel: Observation counts
p_obs <- ggplot(obs_counts, aes(x = event_time, y = N)) +
  annotate("rect", xmin = left_bound - 0.5, xmax = right_bound + 0.5,
           ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.3) +
  geom_col(fill = "gray70", width = 0.7) +
  geom_text(aes(label = scales::comma(N)), vjust = -0.3, size = 2.5) +
  scale_x_continuous(breaks = x_breaks, limits = c(x_min - 0.5, x_max + 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Years Relative to PA Designation", y = "N obs") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 11),
    plot.margin = margin(0, 10, 10, 10)
  )

# Combine panels
p_combined <- p_coef / p_obs + plot_layout(heights = c(3, 1)) +
  plot_annotation(
    title = sprintf("TWFE Event Study: %s", country),
    subtitle = sprintf("Grid + year FE, clustered SEs | Shaded: preferred region [%d, %d]",
                       left_bound, right_bound),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0)
    )
  )

# Create output directory if it doesn't exist
fig_dir <- file.path(project_root, "output/figures/analysis/country_event_studies/")
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
  cat("Created directory:", fig_dir, "\n")
}

# Save figure
fig_path <- file.path(fig_dir, paste0("event_study_twfe_", country, ".png"))
ggsave(fig_path, p_combined, width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved TWFE figure to:", fig_path, "\n")

# ------------------------------------------------------------------------------
# Re-run as a stacked diff-in-diff
# ------------------------------------------------------------------------------

designation_years <- dt[!is.na(desig_year),.N, by = desig_year][order(desig_year),]
did_stack <- list()
for(dyr in designation_years$desig_year){
  cat("Designation year:", dyr, "| N protected cells:", designation_years[desig_year == dyr, N], "\n")
  stacked_element <- dt[(is.na(desig_year) | desig_year == 0) | desig_year == dyr,][, event_time := fifelse(is_protected == TRUE & !is.na(desig_year) & desig_year > 0, year - desig_year, -9999)][,cohort := dyr][,.(cohort, grid_id, year, event_time, is_protected, area_km2, wdpa_pid, country_iso3, Undisturbed_TMF, Degraded_TMF,
                        Deforested, Regrowth)]
  did_stack <- c(did_stack, list(stacked_element))
}
did_stack <- rbindlist(did_stack)
did_stack[,cohort_id := .GRP, by = .(cohort, grid_id)]
did_stack[,cohort_year := .GRP, by = .(cohort, year)]
did_stack[,cluster_id := fifelse(is.na(wdpa_pid), -9999, wdpa_pid)]

stacked_model <- feols(sw(Undisturbed_TMF, Degraded_TMF,
                        Deforested, Regrowth) ~ i(event_time, ref = c(-1, -9999))| cohort_id + cohort_year,
                        data = did_stack, vcov = ~ cluster_id + grid_id,
                        weights = ~ area_km2, lean = TRUE)

cat("Stacked DID estimation complete.\n")

# ------------------------------------------------------------------------------
# Stacked DID Event Study Plot
# ------------------------------------------------------------------------------

cat("\nGenerating stacked DID event study plot...\n")

# Extract coefficients for Undisturbed_TMF (primary outcome)
ct_stack <- coeftable(stacked_model[[1]], keep = "event_time")
coef_stack <- as.data.table(ct_stack)
setnames(coef_stack, c("estimate", "se", "t", "p"))
coef_stack[, event_time := as.integer(gsub("event_time::(-?[0-9]+)", "\\1", rownames(ct_stack)))]

# Add reference period
ref_row_stack <- data.table(estimate = 0, se = 0, t = NA_real_, p = NA_real_, event_time = -1L)
coef_stack <- rbind(coef_stack, ref_row_stack)

# Calculate 95% CIs
coef_stack[, ci_lower := estimate - 1.96 * se]
coef_stack[, ci_upper := estimate + 1.96 * se]
setorder(coef_stack, event_time)

# Get observation counts by event_time (treated obs only, from stacked data)
obs_counts_stack <- did_stack[event_time != -9999, .N, by = event_time][order(event_time)]

# Find "preferred" region based on 50% drop from t=-1
n_ref_stack <- obs_counts_stack[event_time == -1, N]
if (length(n_ref_stack) == 0) n_ref_stack <- max(obs_counts_stack$N)
threshold_stack <- 0.5 * n_ref_stack

# Find left bound
left_times_stack <- obs_counts_stack[event_time < -1][order(-event_time)]
left_bound_stack <- min(obs_counts_stack$event_time)
for (i in seq_len(nrow(left_times_stack))) {
  if (left_times_stack[i, N] < threshold_stack) {
    if (i == 1) {
      left_bound_stack <- -1
    } else {
      left_bound_stack <- left_times_stack[i - 1, event_time]
    }
    break
  }
}

# Find right bound
right_times_stack <- obs_counts_stack[event_time > -1][order(event_time)]
right_bound_stack <- max(obs_counts_stack$event_time)
for (i in seq_len(nrow(right_times_stack))) {
  if (right_times_stack[i, N] < threshold_stack) {
    if (i == 1) {
      right_bound_stack <- -1
    } else {
      right_bound_stack <- right_times_stack[i - 1, event_time]
    }
    break
  }
}

cat("Stacked preferred region: [", left_bound_stack, ", ", right_bound_stack, "]\n")

# X-axis settings
x_breaks_stack <- sort(unique(coef_stack$event_time))
x_min_stack <- min(x_breaks_stack)
x_max_stack <- max(x_breaks_stack)

# Top panel: Event study coefficients
p_coef_stack <- ggplot(coef_stack, aes(x = event_time, y = estimate)) +
  annotate("rect", xmin = left_bound_stack - 0.5, xmax = right_bound_stack + 0.5,
           ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "gray60", linewidth = 0.5) +
  geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), color = "coral", linewidth = 0.8) +
  geom_point(color = "coral", size = 2.5) +
  scale_x_continuous(breaks = x_breaks_stack, limits = c(x_min_stack - 0.5, x_max_stack + 0.5)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  labs(x = NULL, y = "Effect on Undisturbed TMF (pp)") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 11),
    plot.margin = margin(10, 10, 0, 10)
  )

# Bottom panel: Observation counts
p_obs_stack <- ggplot(obs_counts_stack, aes(x = event_time, y = N)) +
  annotate("rect", xmin = left_bound_stack - 0.5, xmax = right_bound_stack + 0.5,
           ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.3) +
  geom_col(fill = "gray70", width = 0.7) +
  geom_text(aes(label = scales::comma(N)), vjust = -0.3, size = 2.5) +
  scale_x_continuous(breaks = x_breaks_stack, limits = c(x_min_stack - 0.5, x_max_stack + 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Years Relative to PA Designation", y = "N obs") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 11),
    plot.margin = margin(0, 10, 10, 10)
  )

# Combine panels
p_combined_stack <- p_coef_stack / p_obs_stack + plot_layout(heights = c(3, 1)) +
  plot_annotation(
    title = sprintf("Stacked DID Event Study: %s", country),
    subtitle = sprintf("Cohort + cohort-year FE, two-way clustered SEs | Shaded: preferred region [%d, %d]",
                       left_bound_stack, right_bound_stack),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0)
    )
  )

# Create output directory for stacked figures
stacked_fig_dir <- file.path(project_root, "output/figures/analysis/stacked_country_event_studies/")
if (!dir.exists(stacked_fig_dir)) {
  dir.create(stacked_fig_dir, recursive = TRUE)
  cat("Created directory:", stacked_fig_dir, "\n")
}

# Save stacked DID figure
stacked_fig_path <- file.path(stacked_fig_dir, paste0("event_study_stacked_", country, ".png"))
ggsave(stacked_fig_path, p_combined_stack, width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved stacked DID figure to:", stacked_fig_path, "\n")

cat("\n==============================================\n")
cat("Completed:", country, "\n")
cat("==============================================")
