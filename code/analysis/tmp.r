# library(data.table)
# library(ggplot2)
# library(magrittr)
# library(Hmisc)
# library(fixest)
# library(MatchIt)
# library(did2s)

# tmp = fread('/Users/mishrap/Documents/fragmentation/Data/temp/peru_panel_sample.csv')

# tmp = tmp[desig_year > 1992 | is.na(desig_year) | is_protected == FALSE,]
# tmp[,forest_1990 := Undisturbed_TMF[year == 1990], by = grid_id]
# tmp[,forest_1991 := Undisturbed_TMF[year == 1991], by = grid_id]
# tmp[,forest_1992 := Undisturbed_TMF[year == 1992], by = grid_id]

# tmp[,init_deforest_trend_1991 :=  forest_1991 - forest_1990, by = grid_id]
# tmp[,init_deforest_trend_1992 :=  forest_1992 - forest_1991, by = grid_id]

# matching_output = matchit(
#     is_protected ~ forest_1990 + forest_1991 + forest_1992 + init_deforest_trend_1991 + init_deforest_trend_1992,
#     data = tmp[year == 1990, .SD[1], by = grid_id],
#     method = "nearest",
#     distance = "euclidean",
#     replace = TRUE,
#     ratio = 3
# )

# plt = tmp[is_protected == TRUE,.(area = sum(area_km2)), by = desig_year][!is.na(desig_year),] %>%
# ggplot() + geom_col(aes(x = desig_year, y = area)) + theme_minimal() +
#   labs(title = "PA Designation Year Distribution (Peru Sample)",
#        x = "Designation Year", y = "Area of Protected Cells, km2")

# plt = tmp[,.(
#   mean_undisturbed = wtd.mean(Undisturbed_TMF, area_km2, na.rm = TRUE),
#   mean_degraded = wtd.mean(Degraded_TMF, area_km2, na.rm = TRUE),
#   mean_deforested = wtd.mean(Deforested, area_km2, na.rm = TRUE),
#   mean_regrowth = wtd.mean(Regrowth, area_km2, na.rm = TRUE)
# ), by = .(year)][order(year)] %>% ggplot() +
#   geom_line(aes(x = year, y = mean_undisturbed), color = 'darkgreen') +
#   geom_line(aes(x = year, y = mean_degraded), color = 'orange') +
#   geom_line(aes(x = year, y = mean_deforested), color = 'red') +
#   geom_line(aes(x = year, y = mean_regrowth), color = 'blue') +
#   labs(title = "Average Land Cover Fractions in Peru (5km grid cells)",
#        x = "Year",
#        y = "Area-Weighted Mean Fraction") +
#   theme_minimal()

# plt = melt(tmp[,.(
#   mean_undisturbed = wtd.mean(Undisturbed_TMF, area_km2, na.rm = TRUE),
#   mean_degraded = wtd.mean(Degraded_TMF, area_km2, na.rm = TRUE),
#   mean_deforested = wtd.mean(Deforested, area_km2, na.rm = TRUE),
#   mean_regrowth = wtd.mean(Regrowth, area_km2, na.rm = TRUE)
# ), by = .(year, is_protected)], id.vars = c('year', 'is_protected')) %>% ggplot() +
#   geom_line(aes(x = year, y = value, color = variable)) +
#   facet_wrap(~ is_protected) +
#   labs(title = "Average Land Cover Fractions in Peru (5km grid cells)",
#        x = "Year",
#        y = "Area-Weighted Mean Fraction") +
#   theme_minimal()

# tmp[,Forest_Share := Undisturbed_TMF] #  + Degraded_TMF + Regrowth
# tmp[, event_time := fifelse(is_protected == TRUE & !is.na(desig_year) & desig_year > 0, year - desig_year, -9999)]


designation_years <- matched_forest[!is.na(desig_year),.N, by = desig_year][order(desig_year),]
did_stack <- list()
for(dyr in designation_years$desig_year){
  cat("Designation year:", dyr, "| N protected cells:", designation_years[desig_year == dyr, N], "\n")
  stacked_element <- matched_forest[(is.na(desig_year) | desig_year == 0) | desig_year == dyr,][, event_time := fifelse(is_protected == TRUE & !is.na(desig_year) & desig_year > 0, year - desig_year, -9999)][,cohort := dyr][,.(cohort, grid_id, year, event_time, is_protected, area_km2, wdpa_pid, country_iso3, Forest_Share, forest_1990, forest_1991, forest_1992)]
  did_stack <- c(did_stack, list(stacked_element))
}
did_stack <- rbindlist(did_stack)
did_stack[,cohort_id := .GRP, by = .(cohort, grid_id)]
did_stack[,cohort_year := .GRP, by = .(cohort, year)]
did_stack[,cluster_id := fifelse(is.na(wdpa_pid), -9999, wdpa_pid)]

evt_model = feols(Forest_Share ~ i(event_time, ref = c(-1, -9999))| cohort_id + cohort_year, 
                    data = did_stack, vcov = ~ cluster_id + grid_id, 
                    weights = ~ area_km2, lean = TRUE)

# match_key <- match.data(matching_output)
# matched_forest <- merge(tmp, match_key[,.(grid_id, weights)], by = "grid_id")
# matched_forest[,weights := weights * area_km2] # reweight by area
# evt_match_model = feols(Forest_Share ~ i(event_time, ref = c(-1, -9999)) | grid_id + year, data = matched_forest, vcov = ~ wdpa_pid, weights = ~ weights)

# tmp[,DID := fifelse(is_protected & event_time >= 0, 1, 0)]
# matched_forest[,DID := fifelse(is_protected & event_time >= 0, 1, 0)]
# evt_did2s = did2s(
#             data = tmp,
#             yname = "Forest_Share",
#             first_stage = ~ 0 | wdpa_pid + year,
#             second_stage = ~ i(event_time, ref = c(-1, -9999)),
#             treatment = 'DID',
#             cluster_var = 'wdpa_pid', 
#             weights = 'area_km2',
#             verbose = TRUE)

# tmp[,gcs := fifelse(!is.na(desig_year) & desig_year > 0, desig_year, 0 )]
# tmp[,grid_id_num := .GRP, by = grid_id]
# tmp[,wdpa_pid_recode := fifelse(is.na(wdpa_pid), -9999, wdpa_pid)]
# forest_att <- att_gt(
#     yname = 'Forest_Share',
#     tname = 'year',
#     idname = 'grid_id_num',
#     gname = 'gcs',
#     data = tmp,
#     weightsname = 'area_km2',
#     base_period = 'universal',
#     cluster = c('grid_id_num', 'wdpa_pid_recode')
# )

# forest_att_dyn <- aggte(forest_att, type = "dynamic", na.rm = TRUE, min_e = -15)

# ggdid(forest_att_dyn)

# forest_att_match <- att_gt(
#     yname = 'Forest_Share',
#     tname = 'year',
#     idname = 'grid_id_num',
#     gname = 'gcs',
#     data = tmp,
#     xformla = ~ forest_1990 + init_deforest_trend_1991 + init_deforest_trend_1992,
#     weightsname = 'area_km2',
#     base_period = 'universal',
#     cluster = c('grid_id_num', 'wdpa_pid_recode')
# )

# forest_att_dyn_match <- aggte(forest_att_match, type = "dynamic", na.rm = TRUE, min_e = -15)

# ggdid(forest_att_dyn_match)

# evt_match_did2s = did2s(
#             data = matched_forest,
#             yname = "Forest_Share",
#             first_stage = ~ 0 | wdpa_pid + year,
#             second_stage = ~ i(event_time, ref = c(-1, -9999)),
#             treatment = 'DID',
#             cluster_var = 'wdpa_pid', 
#             weights = 'weights',
#             verbose = TRUE)

# ------------------------------------------------------------------------------
# Plot event study coefficients + observation counts
# ------------------------------------------------------------------------------

# Extract coefficients
ct <- coeftable(evt_model, keep = "event_time")
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
obs_counts <- tmp[event_time != -9999, .N, by = event_time][order(event_time)]

# ------------------------------------------------------------------------------
# Find "preferred" region based on 50% drop from t=-1
# ------------------------------------------------------------------------------

# Get obs count at reference period (t = -1)
n_ref <- obs_counts[event_time == -1, N]
threshold <- 0.5 * n_ref

# Find left bound: move left from t=-1 until N < threshold
left_times <- obs_counts[event_time < -1][order(-event_time)]  # descending (closest to -1 first)
left_bound <- min(obs_counts$event_time)  # default to min if all pass
for (i in seq_len(nrow(left_times))) {
  if (left_times[i, N] < threshold) {
    # This bin fails; bound is the previous bin (or -1 if first fails)
    if (i == 1) {
      left_bound <- -1
    } else {
      left_bound <- left_times[i - 1, event_time]
    }
    break
  }
}

# Find right bound: move right from t=-1 until N < threshold
right_times <- obs_counts[event_time > -1][order(event_time)]  # ascending
right_bound <- max(obs_counts$event_time)  # default to max if all pass
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
cat("Reference N (t=-1):", n_ref, "| Threshold (50%):", threshold, "\n")

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
  labs(x = NULL, y = "Effect on Forest Share (pp)") +
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
library(patchwork)
p_combined <- p_coef / p_obs + plot_layout(heights = c(3, 1)) +
  plot_annotation(
    title = "Event Study: Effect of PA Designation on Forest Share (Peru)",
    subtitle = sprintf("TWFE with grid + year FE, clustered SEs | Shaded: preferred region [%d, %d]",
                       left_bound, right_bound),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0)
    )
  )

print(p_combined)

# ------------------------------------------------------------------------------
# Overlay plot: Original vs Matched event study
# ------------------------------------------------------------------------------

# Extract coefficients from matched model
ct_match <- coeftable(evt_match_model, keep = "event_time")
coef_match <- as.data.table(ct_match)
setnames(coef_match, c("estimate", "se", "t", "p"))
coef_match[, event_time := as.integer(gsub("event_time::(-?[0-9]+)", "\\1", rownames(ct_match)))]

# Add reference period
ref_row_match <- data.table(estimate = 0, se = 0, t = NA_real_, p = NA_real_, event_time = -1L)
coef_match <- rbind(coef_match, ref_row_match)

# Calculate 95% CIs
coef_match[, ci_lower := estimate - 1.96 * se]
coef_match[, ci_upper := estimate + 1.96 * se]
setorder(coef_match, event_time)

# Add model labels
coef_df[, model := "Full Sample"]
coef_match[, model := "Matched Sample"]

# Combine
coef_overlay <- rbind(coef_df, coef_match)

# Offset for visibility
coef_overlay[, x_pos := event_time + fifelse(model == "Full Sample", -0.15, 0.15)]

# Overlay plot
p_overlay <- ggplot(coef_overlay, aes(x = x_pos, y = estimate, color = model)) +
  annotate("rect", xmin = left_bound - 0.5, xmax = right_bound + 0.5,
           ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "gray60", linewidth = 0.5) +
  geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Full Sample" = "steelblue", "Matched Sample" = "coral")) +
  scale_x_continuous(breaks = x_breaks, limits = c(x_min - 0.5, x_max + 0.5)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  labs(
    x = "Years Relative to PA Designation",
    y = "Effect on Forest Share (pp)",
    title = "Event Study: Full Sample vs Matched Sample",
    subtitle = "TWFE with grid + year FE, clustered SEs | Shaded: preferred region",
    color = NULL
  ) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 11),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0),
    legend.position = "bottom",
    plot.margin = margin(15, 15, 15, 15)
  )

print(p_overlay)

# # ==============================================================================
# # EXPLORATORY MAP: Undisturbed TMF 1990 vs 2023 with PA Boundaries
# # ==============================================================================

# cat("\n\n========== EXPLORATORY MAP: TMF 1990 vs 2023 ==========\n\n")

# # ------------------------------------------------------------------------------
# # Prepare data for mapping
# # ------------------------------------------------------------------------------

# # Get 1990 and 2023 data
# map_1990 <- tmp[year == 1990, .(
#   grid_id, centroid_lon, centroid_lat, is_protected, wdpa_pid,
#   Undisturbed_TMF_1990 = Undisturbed_TMF
# )]

# map_2023 <- tmp[year == 2023, .(
#   grid_id,
#   Undisturbed_TMF_2023 = Undisturbed_TMF
# )]

# # Merge
# map_df <- merge(map_1990, map_2023, by = "grid_id", all.x = TRUE)

# # Calculate change
# map_df[, tmf_change := Undisturbed_TMF_2023 - Undisturbed_TMF_1990]

# cat("Map data prepared:\n")
# cat("  Grid cells:", nrow(map_df), "\n")
# cat("  Protected cells:", sum(map_df$is_protected, na.rm = TRUE), "\n")
# cat("  Mean TMF 1990:", round(mean(map_df$Undisturbed_TMF_1990, na.rm = TRUE), 3), "\n")
# cat("  Mean TMF 2023:", round(mean(map_df$Undisturbed_TMF_2023, na.rm = TRUE), 3), "\n")
# cat("  Mean change:", round(mean(map_df$tmf_change, na.rm = TRUE), 3), "\n")

# # ------------------------------------------------------------------------------
# # Create PA boundary polygons (convex hulls per PA)
# # ------------------------------------------------------------------------------

# # Get PA boundaries by finding the convex hull of each protected area
# library(sf)

# # Create sf object from protected cells only
# pa_cells <- map_df[is_protected == TRUE & !is.na(centroid_lon) & !is.na(centroid_lat)]

# if (nrow(pa_cells) > 0) {
#   pa_sf <- st_as_sf(pa_cells, coords = c("centroid_lon", "centroid_lat"), crs = 4326)

#   # Create convex hull for each wdpa_pid (simplified PA boundaries)
#   pa_boundaries <- pa_sf %>%
#     group_by(wdpa_pid) %>%
#     summarise(geometry = st_convex_hull(st_combine(geometry)), .groups = "drop")

#   cat("  PA boundaries created:", nrow(pa_boundaries), "protected areas\n")
# } else {
#   pa_boundaries <- NULL
#   cat("  Warning: No protected cells found for boundaries\n")
# }

# # # ------------------------------------------------------------------------------
# # # Panel 1: Undisturbed TMF 1990
# # # ------------------------------------------------------------------------------

# p_1990 <- ggplot() +
#   geom_point(data = map_df[!is.na(centroid_lon)],
#              aes(x = centroid_lon, y = centroid_lat, color = Undisturbed_TMF_1990),
#              size = 0.8, alpha = 0.8) +
#   scale_color_viridis_c(limits = c(0, 1), option = "G", direction = -1,
#                         name = "Undisturbed\nTMF") +
#   coord_sf() +
#   labs(title = "1990", x = NULL, y = NULL) +
#   theme(
#     panel.background = element_rect(fill = "gray95"),
#     panel.grid = element_blank(),
#     axis.text = element_text(size = 8),
#     plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
#     legend.position = "right"
#   )

# # ------------------------------------------------------------------------------
# # Panel 2: Undisturbed TMF 2023
# # ------------------------------------------------------------------------------

# p_2023 <- ggplot() +
#   geom_point(data = map_df[!is.na(centroid_lon)],
#              aes(x = centroid_lon, y = centroid_lat, color = Undisturbed_TMF_2023),
#              size = 0.8, alpha = 0.8) +
#   scale_color_viridis_c(limits = c(0, 1), option = "G", direction = -1,
#                         name = "Undisturbed\nTMF") +
#   coord_sf() +
#   labs(title = "2023", x = NULL, y = NULL)+ theme_void() +
#   theme(
#     panel.background = element_rect(fill = "gray95"),
#     panel.grid = element_blank(),
#     axis.text = element_text(size = 8),
#     plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
#     legend.position = "right"
#   )

# # ------------------------------------------------------------------------------
# # Panel 3: Change (2023 - 1990)
# # ------------------------------------------------------------------------------

# p_change <- ggplot() +
#   geom_point(data = map_df[!is.na(centroid_lon)],
#              aes(x = centroid_lon, y = centroid_lat, color = tmf_change),
#              size = 0.8, alpha = 0.8) +
#   scale_color_gradient2(low = "darkred", mid = "white", high = "darkgreen",
#                         midpoint = 0, limits = c(-1, 1),
#                         name = "Change\n(2023-1990)") +
#   coord_sf() +
#   labs(title = "Change (2023 - 1990)", x = NULL, y = NULL) + theme_void() +
#   theme(
#     panel.background = element_rect(fill = "gray95"),
#     panel.grid = element_blank(),
#     axis.text = element_text(size = 8),
#     plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
#     legend.position = "right"
#   ) 

# # ------------------------------------------------------------------------------
# # Combine panels
# # ------------------------------------------------------------------------------

# p_map_combined <- (p_1990 | p_2023) / p_change +
#   plot_annotation(
#     title = "Undisturbed TMF: 1990 vs 2023 (Peru Sample)",
#     subtitle = "Red/black outlines show protected area boundaries (convex hulls)",
#     theme = theme_void() + theme(
#       plot.title = element_text(face = "bold", size = 14, hjust = 0),
#       plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0)
#     )
#   )

# print(p_map_combined)

# map_df = map_df[map_df$is_protected == TRUE,]

# # # ------------------------------------------------------------------------------
# # # Panel 1: Undisturbed TMF 1990
# # # ------------------------------------------------------------------------------

# p_1990 <- ggplot() +
#   geom_point(data = map_df[!is.na(centroid_lon)],
#              aes(x = centroid_lon, y = centroid_lat, color = Undisturbed_TMF_1990),
#              size = 0.8, alpha = 0.8) +
#   scale_color_viridis_c(limits = c(0, 1), option = "G", direction = -1,
#                         name = "Undisturbed\nTMF") +
#   coord_sf() +
#   labs(title = "1990", x = NULL, y = NULL) +
#   theme(
#     panel.background = element_rect(fill = "gray95"),
#     panel.grid = element_blank(),
#     axis.text = element_text(size = 8),
#     plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
#     legend.position = "right"
#   )

# # ------------------------------------------------------------------------------
# # Panel 2: Undisturbed TMF 2023
# # ------------------------------------------------------------------------------

# p_2023 <- ggplot() +
#   geom_point(data = map_df[!is.na(centroid_lon)],
#              aes(x = centroid_lon, y = centroid_lat, color = Undisturbed_TMF_2023),
#              size = 0.8, alpha = 0.8) +
#   scale_color_viridis_c(limits = c(0, 1), option = "G", direction = -1,
#                         name = "Undisturbed\nTMF") +
#   coord_sf() +
#   labs(title = "2023", x = NULL, y = NULL)+ theme_void() +
#   theme(
#     panel.background = element_rect(fill = "gray95"),
#     panel.grid = element_blank(),
#     axis.text = element_text(size = 8),
#     plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
#     legend.position = "right"
#   )

# # ------------------------------------------------------------------------------
# # Panel 3: Change (2023 - 1990)
# # ------------------------------------------------------------------------------

# p_change <- ggplot() +
#   geom_point(data = map_df[!is.na(centroid_lon)],
#              aes(x = centroid_lon, y = centroid_lat, color = tmf_change),
#              size = 0.8, alpha = 0.8) +
#   scale_color_gradient2(low = "darkred", mid = "white", high = "darkgreen",
#                         midpoint = 0, limits = c(-1, 1),
#                         name = "Change\n(2023-1990)") +
#   coord_sf() +
#   labs(title = "Change (2023 - 1990)", x = NULL, y = NULL) + theme_void() +
#   theme(
#     panel.background = element_rect(fill = "gray95"),
#     panel.grid = element_blank(),
#     axis.text = element_text(size = 8),
#     plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
#     legend.position = "right"
#   ) 

# # ------------------------------------------------------------------------------
# # Combine panels
# # ------------------------------------------------------------------------------

# p_map_combined <- (p_1990 | p_2023) / p_change +
#   plot_annotation(
#     title = "Undisturbed TMF: 1990 vs 2023 (Peru Sample)",
#     subtitle = "Red/black outlines show protected area boundaries (convex hulls)",
#     theme = theme_void() + theme(
#       plot.title = element_text(face = "bold", size = 14, hjust = 0),
#       plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0)
#     )
#   )

# print(p_map_combined)


# cat("\nSaved: output/figures/analysis/map_tmf_1990_vs_2023.png\n")
# stopifnot(0)
# # ------------------------------------------------------------------------------
# # Summary by protection status
# # ------------------------------------------------------------------------------

# cat("\nSummary by protection status:\n")
# summary_by_protection <- map_df[, .(
#   n_cells = .N,
#   mean_1990 = round(mean(Undisturbed_TMF_1990, na.rm = TRUE), 3),
#   mean_2023 = round(mean(Undisturbed_TMF_2023, na.rm = TRUE), 3),
#   mean_change = round(mean(tmf_change, na.rm = TRUE), 3),
#   pct_loss = round(100 * mean(tmf_change < -0.1, na.rm = TRUE), 1),
#   pct_stable = round(100 * mean(abs(tmf_change) <= 0.1, na.rm = TRUE), 1),
#   pct_gain = round(100 * mean(tmf_change > 0.1, na.rm = TRUE), 1)
# ), by = is_protected]

# print(summary_by_protection)

# cat("\n========== EXPLORATORY MAP COMPLETE ==========\n")