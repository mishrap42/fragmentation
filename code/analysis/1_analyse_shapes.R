library(here)        # Manage file paths relative to project root
here::i_am("code/analysis/1_analyse_shapes.R")
library(sf)          # Spatial data handling (reading shapefiles, geometric operations)
library(ggplot2)     # Data visualization and plotting
library(dplyr)       # Data manipulation (filtering, grouping, summarizing)
library(gridExtra)   # Arrange multiple plots side-by-side
library(readxl)      # Read Excel files
library(binsreg)     # Binned regression plots
library(lubridate)   # Date/time manipulation
library(Hmisc)       # Miscellaneous data analysis functions
library(fixest)      # Fast fixed-effects regressions
library(countrycode) # Convert country codes/names
library(data.table)
library(janitor)      # Data cleaning functions
library(quantreg)      # Quantile regression

set.seed(42)

# This stores your repository path as a function "here()"
if (grepl("ahaanj", here())) {
  dropbox_dir <- "/Users/ahaanj/Dartmouth College Dropbox/Ahaan Jindal/Protected Area Fragmentation/"
  i = 1
} else if (grepl("mishrap", here())) {
  dropbox_dir <- "/Users/mishrap/Dropbox (Personal)/Protected Area Fragmentation/"
  i = 1
} 

intermediate_data_dir = file.path(dropbox_dir, 'Data/build')
figure_dir = file.path(dropbox_dir, 'Figures', '1_analyse_shapes')
if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

# ========== LOAD CONSOLIDATED SHAPE METRICS DATASET ==========
pa_shapes_complete <- fread(file.path(intermediate_data_dir, 'pa_shapes_complete.gz'))
pa_shapes_complete <- pa_shapes_complete %>%
  clean_names()
pa_shapefile = st_read(file.path(dropbox_dir, 'Data/raw/protected sites/protectedsites.shp'))
pa_shapefile$area = units::drop_units(units::set_units(st_area(pa_shapefile), 'km^2'))
pa_shapes_complete = merge(pa_shapes_complete, pa_shapefile %>% st_drop_geometry() %>% select(cddaId, area), by.x = 'cdda_id', by.y = c('cddaId'))

pa_shapes_complete = pa_shapes_complete[legl_fnd < 2025 & legl_fnd > 1950]

max_year_prior_to_natura = pa_shapes_complete[legl_fnd < 1992,.(max(legl_fnd, na.rm = TRUE))]
m.resid = feols(dscnnc ~ i(legl_fnd, ref = .[max_year_prior_to_natura]) + splines::bs(area, 5) | cdd_cntr^iucn_ctg, data = pa_shapes_complete)

# Extract year coefficients (excluding area)
coef_names <- names(coef(m.resid))
year_coefs <- coef_names[grepl("legl_fnd::", coef_names)]

coef_data <- data.frame(
  year = as.numeric(gsub("legl_fnd::", "", year_coefs)),
  coef = coef(m.resid)[year_coefs],
  se = se(m.resid)[year_coefs]
)

# Add reference year (coefficient = 0) - dynamically get it from the model
ref_year <- as.numeric(max_year_prior_to_natura)
coef_data <- rbind(
  data.frame(year = ref_year, coef = 0, se = 0),
  coef_data
) %>% arrange(year)

# 95% confidence intervals
coef_data$ci_lower <- coef_data$coef - 1.96 * coef_data$se
coef_data$ci_upper <- coef_data$coef + 1.96 * coef_data$se

# AER-style plot
p_aer <- ggplot(coef_data, aes(x = year, y = coef)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "gray50") +
  geom_line(linewidth = 0.8, color = "black") +
  geom_point(size = 1.5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.5) +
  geom_vline(xintercept = 1992, linetype = "dashed", color = "red", linewidth = 0.5) +
  annotate("text", x = 1992, y = 0.25, label = "Natura 2000",
           color = "red", hjust = -0.1, size = 3.5) +
  labs(
    x = "Legal Foundation Year",
    y = "Disconnection Index Coefficient",
    title = "Evolution of Protected Area Fragmentation Over Time",
    subtitle = "Controlling for protected area size"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )

ggsave(file.path(figure_dir, "disconnection_by_year_aer.pdf"),
       p_aer, width = 8, height = 5)
ggsave(file.path(figure_dir, "disconnection_by_year_aer.png"),
       p_aer, width = 8, height = 5, dpi = 300)

print(p_aer)

m.iucn = feols(dscnnc ~ i(legl_fnd, ref = .[max_year_prior_to_natura]) + splines::bs(area, 5) | cdd_cntr, data = pa_shapes_complete, split = ~ iucn_ctg)

# Extract split-specific coefficients
split_coef_data <- list()
split_names <- names(m.iucn)

for (i in seq_along(m.iucn)) {
  coef_names <- names(coef(m.iucn[[i]]))
  year_coefs <- coef_names[grepl("legl_fnd::", coef_names)]

  if (length(year_coefs) > 0) {
    split_df <- data.frame(
      year = as.numeric(gsub("legl_fnd::", "", year_coefs)),
      coef = coef(m.iucn[[i]])[year_coefs],
      se = se(m.iucn[[i]])[year_coefs],
      iucn_category = split_names[i]
    )

    # Add reference year
    split_df <- rbind(
      data.frame(year = as.numeric(max_year_prior_to_natura),
                 coef = 0, se = 0, iucn_category = split_names[i]),
      split_df
    )

    split_coef_data[[i]] <- split_df
  }
}

# Combine all splits
split_coef_data <- bind_rows(split_coef_data) %>%
  mutate(
    ci_lower = coef - 1.96 * se,
    ci_upper = coef + 1.96 * se
  )

# AER-style plot with facets by IUCN category
p_iucn <- ggplot(split_coef_data, aes(x = year, y = coef)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "gray50") +
  geom_line(linewidth = 0.6, color = "black") +
  geom_point(size = 1, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.4) +
  geom_vline(xintercept = 1992, linetype = "dashed", color = "red", linewidth = 0.4) +
  facet_wrap(~ iucn_category, scales = "free_y") +
  labs(
    x = "Legal Foundation Year",
    y = "Disconnection Index Coefficient",
    title = "Evolution of Protected Area Fragmentation by IUCN Category",
    subtitle = "Controlling for area and country fixed effects"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 9)
  )

ggsave(file.path(figure_dir, "disconnection_by_year_iucn.pdf"),
       p_iucn, width = 10, height = 8)
ggsave(file.path(figure_dir, "disconnection_by_year_iucn.png"),
       p_iucn, width = 10, height = 8, dpi = 300)

print(p_iucn)



m.country = feols(dscnnc ~ i(legl_fnd, ref = .[max_year_prior_to_natura]) + splines::bs(area, 5) | iucn_ctg, data = pa_shapes_complete, split = ~ cdd_cntr)

# Extract country-specific coefficients
country_coef_data <- list()
split_labels <- names(m.country)

for (i in seq_along(m.country)) {
  # Extract country code from split name like "sample.var: cdd_cntr; sample: XK"
  country_code <- gsub(".*sample: ([A-Z]{2}).*", "\\1", split_labels[i])

  coef_names <- names(coef(m.country[[i]]))
  year_coefs <- coef_names[grepl("legl_fnd::", coef_names)]

  if (length(year_coefs) > 0) {
    country_df <- data.frame(
      year = as.numeric(gsub("legl_fnd::", "", year_coefs)),
      coef = coef(m.country[[i]])[year_coefs],
      se = se(m.country[[i]])[year_coefs],
      country_code = country_code
    )

    # Add reference year
    country_df <- rbind(
      data.frame(year = as.numeric(max_year_prior_to_natura),
                 coef = 0, se = 0, country_code = country_code),
      country_df
    )

    country_coef_data[[i]] <- country_df
  }
}

# Combine all countries
country_coef_data <- rbindlist(country_coef_data)
country_coef_data[, ci_lower := coef - 1.96 * se]
country_coef_data[, ci_upper := coef + 1.96 * se]
country_coef_data[, country_name := countrycode(country_code, origin = "eurostat", destination = "country.name.en")]

# Drop countries with too few observations
country_n_obs <- country_coef_data[, .N, by = country_code]
countries_to_keep <- country_n_obs[N >= 10, country_code]
country_coef_data <- country_coef_data[country_code %in% countries_to_keep]

# Calculate average trend for each country (mean coefficient across years)
country_trends <- country_coef_data[, .(trend = mean(coef[year >= 1992], na.rm = TRUE) -
mean(coef[year < 1992 & year > 1980], na.rm = TRUE)), by = .(country_code, country_name)]
setorder(country_trends, -trend)

# Get top 5 most positive trends
top5_positive <- country_trends[1:5, country_code]

# Get 5 countries with trends closest to zero (least clear trends)
country_trends[, abs_coef := abs(trend)]
setorder(country_trends, abs_coef)
top5_unclear <- country_trends[1:5, country_code]

# Filter data for each group
country_coef_positive <- country_coef_data[country_code %in% top5_positive]
country_coef_unclear <- country_coef_data[country_code %in% top5_unclear]

# Plot: Most positive trends
p_country_positive <- ggplot(country_coef_positive, aes(x = year, y = coef)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "gray50") +
  geom_line(linewidth = 0.6, color = "black") +
  geom_point(size = 1, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.4) +
  geom_vline(xintercept = 1992, linetype = "dashed", color = "red", linewidth = 0.4) +
  facet_wrap(~ country_name, scales = "free_y") +
  labs(
    x = "Legal Foundation Year",
    y = "Disconnection Index Coefficient",
    title = "Countries with Most Positive Fragmentation Trends",
    subtitle = "Controlling for area and IUCN category"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 9)
  )

ggsave(file.path(figure_dir, "disconnection_country_positive.pdf"),
       p_country_positive, width = 10, height = 8)
ggsave(file.path(figure_dir, "disconnection_country_positive.png"),
       p_country_positive, width = 10, height = 8, dpi = 300)

# Plot: Least clear trends
p_country_unclear <- ggplot(country_coef_unclear, aes(x = year, y = coef)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "gray50") +
  geom_line(linewidth = 0.6, color = "black") +
  geom_point(size = 1, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.4) +
  geom_vline(xintercept = 1992, linetype = "dashed", color = "red", linewidth = 0.4) +
  facet_wrap(~ country_name, scales = "free_y") +
  labs(
    x = "Legal Foundation Year",
    y = "Disconnection Index Coefficient",
    title = "Countries with Least Clear Fragmentation Trends",
    subtitle = "Controlling for area and IUCN category"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 9)
  )

ggsave(file.path(figure_dir, "disconnection_country_unclear.pdf"),
       p_country_unclear, width = 10, height = 8)
ggsave(file.path(figure_dir, "disconnection_country_unclear.png"),
       p_country_unclear, width = 10, height = 8, dpi = 300)

print(p_country_positive)
print(p_country_unclear)


# Quantile regression analysis
pa_shapes_complete[,post_natura := fifelse(legl_fnd >= 1992, 1, 0)]
pa_shapes_complete[,pre_natura := fifelse(legl_fnd < 1992, 1, 0)]

quantile_model = rq(log(dscnnc) ~ legl_fnd * post_natura + log(area),
                    tau = c(0.1, 0.25, 0.5, 0.75, 0.9),
                    data = pa_shapes_complete)

# Create regression table manually
summary_qr <- summary(quantile_model)
taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)

# Extract coefficients and standard errors
coef_list <- lapply(summary_qr, function(x) x$coefficients[, 1])
se_list <- lapply(summary_qr, function(x) x$coefficients[, 2])

var_names <- c("Year", "Post-Natura 2000", "Log(Area)", "Year $\\times$ Post-Natura")

# Build LaTeX table
latex_lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Quantile Regression: Protected Area Fragmentation}",
  "\\label{tab:quantile_reg}",
  "\\begin{tabular}{lccccc}",
  "\\hline\\hline",
  " & Q10 & Q25 & Q50 & Q75 & Q90 \\\\",
  "\\hline"
)

# Add rows for each variable
for (i in seq_along(var_names)) {
  coef_row <- paste0(var_names[i], " & ",
                     paste(sprintf("%.4f", sapply(coef_list, `[`, i)), collapse = " & "),
                     " \\\\")
  se_row <- paste0(" & ",
                   paste(sprintf("(%.4f)", sapply(se_list, `[`, i)), collapse = " & "),
                   " \\\\")
  latex_lines <- c(latex_lines, coef_row, se_row)
}

# Calculate % change in trend from interaction
# % change = (Year x Post-Natura coefficient / Year coefficient) * 100
pct_change <- sapply(coef_list, function(x) (x[4] / x[1]) * 100)
pct_change_row <- paste0("\\% Change in Trend & ",
                         paste(sprintf("%.1f\\%%", pct_change), collapse = " & "),
                         " \\\\")

latex_lines <- c(latex_lines, "\\hline", pct_change_row)

# Add footer
latex_lines <- c(latex_lines,
                 "\\hline",
                 sprintf("N & %d & %d & %d & %d & %d \\\\",
                         nrow(pa_shapes_complete), nrow(pa_shapes_complete),
                         nrow(pa_shapes_complete), nrow(pa_shapes_complete),
                         nrow(pa_shapes_complete)),
                 "\\hline\\hline",
                 "\\end{tabular}",
                 "\\begin{tablenotes}[flushleft]",
                 "\\small",
                 "\\item \\textit{Notes:} Standard errors in parentheses. Dependent variable is log(Disconnection Index).",
                 "\\end{tablenotes}",
                 "\\end{table}")

# Write to file
writeLines(latex_lines, file.path(figure_dir, "quantile_regression_table.tex"))
cat(paste(latex_lines, collapse = "\n"))