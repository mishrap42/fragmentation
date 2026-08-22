# ==============================================================================
# 3a_matching.R
# Parallelized MatchIt matching - one country per array task
# Usage: Rscript 3a_matching.R <task_id>
# ==============================================================================

library(data.table)
library(duckdb)
library(DBI)
# Note: MatchIt, cobalt loaded AFTER DuckDB query to avoid dependency conflicts

# ------------------------------------------------------------------------------
# Setup paths
# ------------------------------------------------------------------------------

if (Sys.getenv("SLURM_SUBMIT_DIR") != "") {
  project_root <- Sys.getenv("SLURM_SUBMIT_DIR")
} else {
  project_root <- here::here()
}

data_dir <- file.path(project_root, "Data/build/final")
intermediate_dir <- file.path(project_root, "Data/build/intermediate/matching")
fig_dir <- file.path(project_root, "output/figures/analysis")

panel_path <- file.path(data_dir, "TMF_5km_panel.parquet")

# Create directories if needed
if (!dir.exists(intermediate_dir)) dir.create(intermediate_dir, recursive = TRUE)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# Get task ID from command line
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Task ID required. Usage: Rscript 3a_matching.R <task_id>")
}
task_id <- as.integer(args[1])

cat("==============================================\n")
cat("3a_matching.R - Parallel Matching by Country\n")
cat("Task ID:", task_id, "\n")
cat("==============================================\n\n")

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
# Load data for this country only (using SELECT * pattern that works)
# ------------------------------------------------------------------------------

cat("Loading data for country:", country, "...\n")

con <- dbConnect(duckdb())

# Use SELECT * pattern (same as 1_balance_PA.R) then subset in R
query <- sprintf("
  SELECT *
  FROM read_parquet('%s')
  WHERE year BETWEEN 1990 AND 1995 AND country_iso3 = '%s'
", panel_path, country)

dt_match <- as.data.table(dbGetQuery(con, query))
dbDisconnect(con, shutdown = TRUE)

cat("Loaded", format(nrow(dt_match), big.mark = ","), "cells\n")

# Subset to needed columns
dt_match <- dt_match[, .(
  grid_id, country_iso3, year, is_protected, desig_year,
  Undisturbed_TMF, Degraded_TMF, Deforested, forest_cover
)]

gc()

# Load remaining libraries now that DuckDB is done
library(MatchIt)
library(cobalt)
library(fixest)
library(ggplot2)

# ------------------------------------------------------------------------------
# Pivot land use shares wide (1990-1995)
# ------------------------------------------------------------------------------

cat("Pivoting land use shares wide...\n")

share_cols <- c("Undisturbed_TMF", "Degraded_TMF", "Deforested")

dt_wide <- dcast(
  dt_match,
  grid_id + country_iso3 ~ year,
  value.var = share_cols
)
# Resulting columns: Undisturbed_TMF_1990, ..., Deforested_1995

# Bring in cell-invariant protection attributes from the 1990 slice
dt_meta <- unique(dt_match[year == 1990, .(grid_id, is_protected, desig_year)])
dt_wide <- merge(dt_wide, dt_meta, by = "grid_id")

# ------------------------------------------------------------------------------
# Define treatment: protected with desig_year >= 1996 (drop 1991-1995 designations)
# ------------------------------------------------------------------------------

dt_wide[, treated := is_protected == TRUE & !is.na(desig_year) & desig_year >= 1996]
dt_wide[, control := is_protected == FALSE]

# Drop cells designated 1991-1995 (neither treated nor control)
dt_wide <- dt_wide[treated | control]

n_treated <- sum(dt_wide$treated)
n_control <- sum(dt_wide$control)
cat("  Treated (desig_year >= 1996):", n_treated, "\n")
cat("  Control (never protected):  ", n_control, "\n")

if (n_treated == 0) {
  cat("No treated cells in", country, "- skipping\n")
  quit(status = 0)
}
if (n_control == 0) {
  cat("No control cells in", country, "- skipping\n")
  quit(status = 0)
}

# Require complete share trajectory 1990-1995
share_vars <- as.vector(outer(share_cols, 1990:1995, paste, sep = "_"))
missing_shares <- setdiff(share_vars, names(dt_wide))
if (length(missing_shares) > 0) {
  cat("Missing share columns in", country, ":", paste(missing_shares, collapse = ", "), "- skipping\n")
  quit(status = 0)
}

dt_wide_complete <- dt_wide[complete.cases(dt_wide[, ..share_vars])]

if (sum(dt_wide_complete$treated) == 0 || sum(!dt_wide_complete$treated) == 0) {
  cat("Missing treated or control after requiring full share trajectory in", country, "- skipping\n")
  quit(status = 0)
}

cat("Complete cases:", nrow(dt_wide_complete), "\n")
cat("  1990 share means by group:\n")
print(dt_wide_complete[, .(
  undisturbed = mean(Undisturbed_TMF_1990),
  degraded    = mean(Degraded_TMF_1990),
  deforested  = mean(Deforested_1990)
), by = treated])

df_match <- as.data.frame(dt_wide_complete)

# ------------------------------------------------------------------------------
# Run MatchIt: nearest-neighbor on logit propensity score over 1990-1995 shares
# ------------------------------------------------------------------------------

cat("\nRunning MatchIt (logit pscore on 1990-1995 LUC shares) for", country, "...\n")

start_time <- Sys.time()

ps_formula <- as.formula(paste("treated ~", paste(share_vars, collapse = " + ")))

m.out <- matchit(
  ps_formula,
  data = df_match,
  method = "nearest",
  distance = "glm",
  ratio = 1,
  replace = TRUE
)

end_time <- Sys.time()
cat("Matching completed in:", round(difftime(end_time, start_time, units = "secs"), 1), "seconds\n")

# ------------------------------------------------------------------------------
# Extract and save matched data
# ------------------------------------------------------------------------------

matched_data <- match.data(m.out)
matched_dt <- as.data.table(matched_data)

if (nrow(matched_dt) == 0 || sum(matched_dt$treated) == 0 || sum(!matched_dt$treated) == 0) {
  cat("No matched units in", country, "- skipping\n")
  quit(status = 0)
}

cat("\nMatched sample:\n")
cat("  Total:", nrow(matched_dt), "\n")
cat("  Treated:", sum(matched_dt$treated), "\n")
cat("  Controls:", sum(!matched_dt$treated), "\n")

# Save matched pairs for this country
output_file <- file.path(intermediate_dir, sprintf("matched_pairs_%s.csv", country))
matched_pairs <- matched_dt[, .(grid_id, country_iso3, treated, weights)]
fwrite(matched_pairs, output_file)
cat("Saved to:", output_file, "\n")

# ------------------------------------------------------------------------------
# Print balance summary
# ------------------------------------------------------------------------------

cat("\nBalance summary:\n")
print(bal.tab(m.out, stats = c("mean.diffs")))

# ------------------------------------------------------------------------------
# Part 7: Load full panel for this country's matched cells
# ------------------------------------------------------------------------------

cat("\nPart 7: Loading full panel for matched cells...\n")

matched_grid_ids <- matched_dt$grid_id

con <- dbConnect(duckdb())
dbExecute(con, "CREATE OR REPLACE TABLE matched_ids (grid_id VARCHAR)")
dbAppendTable(con, "matched_ids", data.frame(grid_id = matched_grid_ids))

query_panel <- sprintf("
  SELECT p.grid_id, p.area_km2, p.wdpa_pid, p.country_iso3, p.year,
         p.desig_year, p.is_protected,
         p.Undisturbed_TMF, p.Degraded_TMF, p.Deforested, p.Regrowth,
         p.forest_cover
  FROM read_parquet('%s') p
  INNER JOIN matched_ids m ON p.grid_id = m.grid_id
  WHERE p.country_iso3 = '%s'
", panel_path, country)

dt_panel <- as.data.table(dbGetQuery(con, query_panel))
dbDisconnect(con, shutdown = TRUE)

# Merge in MatchIt weights so only matchit-returned treat-control pairs are used
dt_panel <- merge(
  dt_panel,
  matched_dt[, .(grid_id, matchit_weight = weights)],
  by = "grid_id"
)

cat("Loaded", format(nrow(dt_panel), big.mark = ","), "cell-years\n")
cat("Unique cells:", uniqueN(dt_panel$grid_id), "\n")

gc()

# ------------------------------------------------------------------------------
# Part 8: Run TWFE event study (within-country)
# ------------------------------------------------------------------------------

cat("\nPart 8: Running TWFE event study...\n")

# Create event time
dt_panel[, event_time := fifelse(is_protected, year - desig_year, -9999)]
dt_panel[, event_time := fifelse(
  event_time <= -10 & event_time > -9999, -10,
  fifelse(event_time >= 15, 15, event_time)
)]

cat("Event time distribution:\n")
print(dt_panel[is_protected == TRUE, .N, by = event_time][order(event_time)])

# Check if we have variation in event_time
n_event_times <- uniqueN(dt_panel[is_protected == TRUE, event_time])
if (n_event_times < 3) {
  cat("Insufficient event time variation in", country, "- skipping event study\n")
  cat("\n==============================================\n")
  cat("Completed matching only:", country, "\n")
  cat("==============================================\n")
  quit(status = 0)
}

# Within-country TWFE: grid_id + year FE
# Weights: MatchIt stratum weights * area_km2
twfe_model <- feols(
  sw(Undisturbed_TMF, Degraded_TMF, Deforested, Regrowth) ~
    i(event_time, ref = c(-1, -9999)) | grid_id + year,
  data = dt_panel,
  weights = ~ I(matchit_weight * area_km2),
  cluster = ~wdpa_pid,
  lean = TRUE
)

cat("TWFE estimation complete.\n")

# ------------------------------------------------------------------------------
# Part 9: Generate country event study figure
# ------------------------------------------------------------------------------

cat("\nPart 9: Generating event study figure...\n")

outcome_names <- c("Undisturbed_TMF", "Degraded_TMF", "Deforested", "Regrowth")
outcome_labels <- c("Undisturbed TMF", "Degraded TMF", "Deforested", "Regrowth")

# Extract coefficients
coef_list <- lapply(seq_along(outcome_names), function(i) {
  model_i <- twfe_model[[i]]
  ct <- coeftable(model_i, keep = "event_time")
  coef_df <- as.data.table(ct)
  setnames(coef_df, c("estimate", "se", "t", "p"))
  coef_df[, event_time := as.integer(
    gsub("event_time::(-?[0-9]+)", "\\1", rownames(ct))
  )]
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
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40",
             linewidth = 0.5) +
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "gray60",
             linewidth = 0.5) +
  geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), color = "steelblue",
                 linewidth = 0.7) +
  geom_point(color = "steelblue", size = 2) +
  facet_wrap(~ outcome, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  labs(
    x = "Years Relative to PA Designation",
    y = "Effect (pp)",
    title = sprintf("Event Study: %s (Matched Sample)", country),
    subtitle = "TWFE with grid cell + year FE, SEs clustered at PA level"
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

output_fig <- file.path(fig_dir, sprintf("event_study_%s.png", country))
ggsave(output_fig, p_event, width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved:", output_fig, "\n")

cat("\n==============================================\n")
cat("Completed:", country, "\n")
cat("==============================================\n")
