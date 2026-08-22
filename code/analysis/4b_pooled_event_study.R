# ==============================================================================
# 4b_pooled_event_study.R
# Treatment-area-weighted pooled graphs built from the per-PA synthetic control
# output (4_synthetic_control.R -> 4a_consolidate_synth.R).
#
# Reads the consolidated files:
#   Data/build/intermediate/synth/synth_project_effects.csv  (PA areas)
#   Data/build/intermediate/synth/synth_dynamic_effects.csv  (per-PA paths + SE)
#
# Produces, both GLOBALLY (all PAs) and PER COUNTRY:
#   (1) pooled_levels   - area-weighted average observed PA line and area-weighted
#                         average synthetic-control line, on one graph (levels,
#                         no standard errors).
#   (2) pooled_event_study - area-weighted average treatment effect (gap) by event
#                         time, normalized to event_time = -1, as points + 95% CI
#                         error bars; area-weighted average of the per-PA placebo
#                         SEs. Event times -5..-1 (the SC matching window) are
#                         shaded. A lower panel shows the number of contributing
#                         PAs at each event time (a la 2a_country_event_study.R).
#
# Usage: Rscript 4b_pooled_event_study.R
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(fixest)
  library(rnaturalearth)
  library(sf)
  library(classInt)
})

if (Sys.getenv("SLURM_SUBMIT_DIR") != "") {
  project_root <- Sys.getenv("SLURM_SUBMIT_DIR")
} else {
  project_root <- here::here()
}

synth_dir <- file.path(project_root, "Data/build/intermediate/synth")
fig_dir   <- file.path(project_root, "output/figures/analysis/synth_control")

eff_path <- file.path(synth_dir, "synth_project_effects.csv")
dyn_path <- file.path(synth_dir, "synth_dynamic_effects.csv")

if (!file.exists(eff_path) || !file.exists(dyn_path)) {
  stop("Consolidated files not found. Run 4a_consolidate_synth.R first.")
}

effects <- fread(eff_path)
dynamic <- fread(dyn_path)

if (nrow(dynamic) == 0) stop("No dynamic effects to pool.")

# ------------------------------------------------------------------------------
# Attach each PA's treated area as the pooling weight
# ------------------------------------------------------------------------------

dt <- merge(
  dynamic,
  unique(effects[, .(wdpa_pid, country_iso3, project_area_km2)]),
  by = c("wdpa_pid", "country_iso3"), all.x = TRUE
)

# Drop PA-years with no usable area weight
dt <- dt[is.finite(project_area_km2) & project_area_km2 > 0]

# Normalize each PA's gap to event_time = -1 (remove the level at -1)
ref <- dt[event_time == -1L, .(wdpa_pid, country_iso3, tau_ref = tau)]
dt  <- merge(dt, ref, by = c("wdpa_pid", "country_iso3"), all.x = TRUE)
dt[, tau_norm := tau - tau_ref]

# ------------------------------------------------------------------------------
# Plot builders
# ------------------------------------------------------------------------------

# Theme shared by both graphs (matches the repo's void-like event-study style)
base_theme <- theme(
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

# (1) Pooled levels: area-weighted observed vs synthetic, no SEs
make_levels_plot <- function(d, label) {
  lv <- d[, .(
    Observed  = weighted.mean(real_y,  project_area_km2, na.rm = TRUE),
    Synthetic = weighted.mean(synth_y, project_area_km2, na.rm = TRUE)
  ), by = event_time][order(event_time)]

  lvl <- melt(lv, id.vars = "event_time", variable.name = "series",
              value.name = "forest_frac")
  lvl[, series := factor(series, levels = c("Observed", "Synthetic"),
                         labels = c("Observed (area-weighted)",
                                    "Synthetic control (area-weighted)"))]

  ggplot(lvl, aes(event_time, forest_frac, color = series, linetype = series)) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray60",
               linewidth = 0.5) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = c("Observed (area-weighted)" = "steelblue",
                                  "Synthetic control (area-weighted)" = "coral")) +
    scale_linetype_manual(values = c("Observed (area-weighted)" = "solid",
                                     "Synthetic control (area-weighted)" = "dashed")) +
    labs(x = "Years relative to PA designation", y = "Forest fraction",
         color = NULL, linetype = NULL,
         title = sprintf("Area-weighted pooled forest fraction (%s)", label),
         subtitle = "Treated-area-weighted average of observed and synthetic paths") +
    base_theme
}

# (2) Pooled event study: area-weighted gap + area-weighted placebo SE, with a
#     contributing-PA count panel below.
make_event_study_plot <- function(d, label) {
  # Area-weighted mean effect, and the standard error OF that weighted mean:
  # sqrt(weighted variance of tau_norm across PAs / number of PAs).
  est <- d[, .(
    estimate = Hmisc::wtd.mean(tau_norm, project_area_km2, na.rm = TRUE),
    se       = sqrt(Hmisc::wtd.var(tau_norm, project_area_km2, na.rm = TRUE,
                                   normwt = TRUE) / uniqueN(wdpa_pid)),
    n_pa     = uniqueN(wdpa_pid)
  ), by = event_time][order(event_time)]

  # event_time == -1 is the reference period: 0 by construction, no CI
  est[event_time == -1L, `:=`(estimate = 0, se = 0)]
  est[, `:=`(ci_lo = estimate - 1.96 * se, ci_hi = estimate + 1.96 * se)]
  est[event_time == -1L, `:=`(ci_lo = 0, ci_hi = 0)]

  x_breaks <- sort(unique(est$event_time))
  x_min <- min(x_breaks); x_max <- max(x_breaks)

  # Shade the SC matching window, event_time -5..-1 (targeted periods)
  shade <- annotate("rect", xmin = -5.5, xmax = -0.5, ymin = -Inf, ymax = Inf,
                    fill = "gray80", alpha = 0.4)

  p_coef <- ggplot(est, aes(event_time, estimate)) +
    shade +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40",
               linewidth = 0.5) +
    geom_vline(xintercept = -0.5, linetype = "dotted", color = "gray60",
               linewidth = 0.5) +
    geom_linerange(aes(ymin = ci_lo, ymax = ci_hi), color = "steelblue",
                   linewidth = 0.7) +
    geom_point(color = "steelblue", size = 2) +
    scale_x_continuous(breaks = x_breaks, limits = c(x_min - 0.5, x_max + 0.5)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
    labs(x = NULL, y = "Effect on forest fraction (pp)",
         title = sprintf("Pooled synthetic-control event study (%s)", label),
         subtitle = paste0("Area-weighted gap, normalized to t = -1 | ",
                           "shaded = matching window (t = -5..-1)")) +
    base_theme +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(15, 15, 0, 15))

  p_obs <- ggplot(est, aes(event_time, n_pa)) +
    shade +
    geom_col(fill = "gray70", width = 0.7) +
    scale_x_continuous(breaks = x_breaks, limits = c(x_min - 0.5, x_max + 0.5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Years relative to PA designation", y = "N PAs") +
    base_theme +
    theme(plot.margin = margin(0, 15, 15, 15))

  p_coef / p_obs + plot_layout(heights = c(3, 1))
}

# ------------------------------------------------------------------------------
# Write one (levels, event-study) pair for a given subset
# ------------------------------------------------------------------------------

write_pair <- function(d, label, out_dir, suffix) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  lv <- make_levels_plot(d, label)
  ggsave(file.path(out_dir, sprintf("pooled_levels_%s.png", suffix)),
         lv, width = 9, height = 6, dpi = 300, bg = "white")

  es <- make_event_study_plot(d, label)
  ggsave(file.path(out_dir, sprintf("pooled_event_study_%s.png", suffix)),
         es, width = 9, height = 8, dpi = 300, bg = "white")

  cat(sprintf("  %-14s %d PAs, %d event-times -> %s\n",
              label, uniqueN(d$wdpa_pid), uniqueN(d$event_time), out_dir))
}

# ------------------------------------------------------------------------------
# Global pooled graphs
# ------------------------------------------------------------------------------

cat("Building global pooled graphs...\n")
write_pair(dt, "GLOBAL", file.path(fig_dir, "GLOBAL"), "GLOBAL")

# ------------------------------------------------------------------------------
# Per-country pooled graphs
# ------------------------------------------------------------------------------

cat("Building per-country pooled graphs...\n")
for (cc in sort(unique(dt$country_iso3))) {
  d_cc <- dt[country_iso3 == cc]
  if (uniqueN(d_cc$wdpa_pid) < 1) next
  write_pair(d_cc, cc, file.path(fig_dir, cc), cc)
}

cat("\nDone. Pooled graphs written under", fig_dir, "\n")

# ==============================================================================
# Stacked event study with Conley (50 km) standard errors
# ------------------------------------------------------------------------------
# Reads the stacked treated+control panel (synth_stacked_panel.csv) and runs a
# single feols event study, recovering SEs internally:
#   Y ~ i(event_time, ref = c(-1, -9999)) | unit_id^wdpa_pid + wdpa_pid^year
#   weights = reg_weight ; vcov = Conley(lat, lon, cutoff = 50 km)
# Controls carry event_time = -9999 and are absorbed into the reference.
# Produced GLOBALLY and per country (countries with >= MIN_PA_COUNTRY PAs).
# ==============================================================================

CONLEY_CUTOFF  <- 50   # km
MIN_PA_COUNTRY <- 1    # min estimated PAs to attempt a per-country event study

stk_path <- file.path(synth_dir, "synth_stacked_panel.csv")

if (!file.exists(stk_path)) {
  cat("\n[stacked] synth_stacked_panel.csv not found - skipping Conley ES.\n")
} else {
  cat("\n=== Stacked event study (Conley", CONLEY_CUTOFF, "km SEs) ===\n")
  sp <- fread(stk_path)

  # Conley needs coordinates: drop rows with missing lat/lon
  n0 <- nrow(sp)
  sp <- sp[!is.na(lon) & !is.na(lat)]
  cat("Dropped", format(n0 - nrow(sp), big.mark = ","),
      "rows with missing coordinates\n")

  # Floor negative SC weights (solver noise) at 0, then drop zero-weight rows;
  # they contribute nothing to the weighted fit.
  sp[reg_weight < 0, reg_weight := 0]
  n1 <- nrow(sp)
  sp <- sp[reg_weight > 0]
  cat("Dropped", format(n1 - nrow(sp), big.mark = ","),
      "zero/neg-weight rows\n")

  # Composite unit identifier (donor cell reused across cohorts stays distinct)
  sp[, unit_cohort := paste0(unit_id, "_", wdpa_pid)]

  # --- estimator + plotter for one subset ------------------------------------
  make_stacked_es <- function(d, label, suffix, out_dir) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    m <- tryCatch(
      feols(Y ~ i(event_time, ref = c(-1, -9999)) | unit_cohort + wdpa_pid^year,
            data = d, weights = ~reg_weight,
            vcov = vcov_conley(lat = ~lat, lon = ~lon, cutoff = CONLEY_CUTOFF)),
      error = function(e) { cat("  [", label, "] feols failed:",
                                conditionMessage(e), "\n"); NULL })
    if (is.null(m)) return(invisible(NULL))

    ct <- coeftable(m, keep = "event_time")
    cf <- as.data.table(ct)
    setnames(cf, c("estimate", "se", "t", "p"))
    cf[, event_time := as.integer(
      gsub("event_time::(-?[0-9]+)", "\\1", rownames(ct)))]
    cf <- rbind(cf, data.table(estimate = 0, se = 0, t = NA_real_,
                               p = NA_real_, event_time = -1L))
    cf[, ci_lo := estimate - 1.96 * se]
    cf[, ci_hi := estimate + 1.96 * se]
    cf[event_time == -1L, `:=`(ci_lo = 0, ci_hi = 0)]
    setorder(cf, event_time)
    fwrite(cf, file.path(out_dir, sprintf("stacked_event_study_%s_coefs.csv",
                                          suffix)))

    pa_counts <- d[treated == 1, .(n_pa = uniqueN(wdpa_pid)), by = event_time]

    x_breaks <- sort(unique(cf$event_time))
    x_min <- min(x_breaks); x_max <- max(x_breaks)
    shade <- annotate("rect", xmin = -5.5, xmax = -0.5, ymin = -Inf, ymax = Inf,
                      fill = "gray80", alpha = 0.4)

    p_coef <- ggplot(cf, aes(event_time, estimate)) +
      shade +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40",
                 linewidth = 0.5) +
      geom_vline(xintercept = -0.5, linetype = "dotted", color = "gray60",
                 linewidth = 0.5) +
      geom_linerange(aes(ymin = ci_lo, ymax = ci_hi), color = "steelblue",
                     linewidth = 0.7) +
      geom_point(color = "steelblue", size = 2) +
      scale_x_continuous(breaks = x_breaks, limits = c(x_min - 0.5, x_max + 0.5)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
      labs(x = NULL, y = "Effect on forest fraction (pp)",
           title = sprintf("Stacked synthetic-control event study (%s)", label),
           subtitle = sprintf(paste0("feols, area*SC weights, unit + cohort-year",
                                     " FE | Conley %d km SEs | shaded = matching",
                                     " window (t = -5..-1)"), CONLEY_CUTOFF)) +
      base_theme +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            plot.margin = margin(15, 15, 0, 15))

    p_obs <- ggplot(pa_counts, aes(event_time, n_pa)) +
      shade +
      geom_col(fill = "gray70", width = 0.7) +
      scale_x_continuous(breaks = x_breaks, limits = c(x_min - 0.5, x_max + 0.5)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(x = "Years relative to PA designation", y = "N PAs") +
      base_theme +
      theme(plot.margin = margin(0, 15, 15, 15))

    p <- p_coef / p_obs + plot_layout(heights = c(3, 1))
    ggsave(file.path(out_dir, sprintf("stacked_event_study_%s.png", suffix)),
           p, width = 10, height = 8, dpi = 300, bg = "white")
    cat(sprintf("  %-14s %d PAs -> %s\n", label, uniqueN(d$wdpa_pid), out_dir))
  }

  # Global
  make_stacked_es(sp, "GLOBAL", "GLOBAL", file.path(fig_dir, "GLOBAL"))

  # Per country (>= MIN_PA_COUNTRY estimated PAs)
  for (cc in sort(unique(sp$country_iso3))) {
    d_cc <- sp[country_iso3 == cc]
    if (uniqueN(d_cc[treated == 1, wdpa_pid]) < MIN_PA_COUNTRY) next
    make_stacked_es(d_cc, cc, cc, file.path(fig_dir, cc))
  }

  cat("Stacked event studies complete.\n")
}

# ==============================================================================
# Country ATT at horizon g: choropleth + coefplots
# ------------------------------------------------------------------------------
# g (per country) = min(max event_time in that country's stacked ES, G_CAP).
# The ATT at g and its Conley SE are read from the per-country stacked ES coef
# tables written above. Map colors are Jenks (natural-breaks) bins of the
# winsorized ATT; coefplots show the raw point estimate + Conley CIs.
# ==============================================================================

G_CAP    <- 15     # horizon cap for g
WINS_P   <- 0.05   # winsorize ATT at the 5th / 95th percentiles for map colors
N_JENKS  <- 6      # number of Jenks color bins

coef_files <- list.files(fig_dir, recursive = TRUE,
                         pattern = "^stacked_event_study_.*_coefs\\.csv$",
                         full.names = TRUE)
coef_files <- coef_files[!grepl("GLOBAL", coef_files)]

if (length(coef_files) == 0) {
  cat("\n[att-g] No per-country stacked coef tables - skipping ATT-g figures.\n")
} else {
  cat("\n=== Country ATT at g: map + coefplots ===\n")
  cat("Found", length(coef_files), "per-country stacked coef tables\n")

  att_g <- rbindlist(lapply(coef_files, function(f) {
    cc <- sub(".*stacked_event_study_(.*)_coefs\\.csv$", "\\1", basename(f))
    ct <- fread(f)
    g <- min(max(ct$event_time, na.rm = TRUE), G_CAP)
    row <- ct[event_time == g]
    if (nrow(row) == 0) row <- ct[event_time == max(event_time[event_time <= g])]
    data.table(country_iso3 = cc, g = g, att = row$estimate[1],
               se = row$se[1], ci_lo = row$ci_lo[1], ci_hi = row$ci_hi[1])
  }))

  area_by_country <- effects[, .(area = sum(project_area_km2, na.rm = TRUE),
                                 n_pa = uniqueN(wdpa_pid)), by = country_iso3]
  att_g <- merge(att_g, area_by_country, by = "country_iso3", all.x = TRUE)

  gdir <- file.path(fig_dir, "GLOBAL")
  if (!dir.exists(gdir)) dir.create(gdir, recursive = TRUE)
  fwrite(att_g, file.path(gdir, "att_at_g_by_country.csv"))
  cat("Countries with ATT at g:", nrow(att_g), "\n")

  # --- (1) Choropleth: Jenks bins of winsorized ATT, hard break at 0 ---------
  # Negatives and positives are binned separately with 0 as a forced boundary,
  # so no bin straddles zero. Negative bins are reds, positive bins are blues.
  cat("(1) Choropleth...\n")
  qlo <- quantile(att_g$att, WINS_P,     na.rm = TRUE)
  qhi <- quantile(att_g$att, 1 - WINS_P, na.rm = TRUE)
  att_g[, att_w := pmin(pmax(att, qlo), qhi)]

  jenks_breaks <- function(x, n) {
    x <- x[is.finite(x)]
    if (length(unique(x)) <= 1) return(range(x))
    n <- min(n, length(unique(x)))
    unique(classInt::classIntervals(x, n = n, style = "jenks")$brks)
  }

  neg <- att_g[att_w < 0, att_w]
  pos <- att_g[att_w > 0, att_w]
  half <- max(2, floor(N_JENKS / 2))
  neg_brks <- if (length(neg)) jenks_breaks(neg, half) else NULL
  pos_brks <- if (length(pos)) jenks_breaks(pos, half) else NULL

  # Assemble a single increasing break vector with 0 as an exact boundary
  brks <- sort(unique(c(neg_brks[neg_brks < 0], 0, pos_brks[pos_brks > 0])))
  att_g[, att_bin := cut(att_w, breaks = brks, include.lowest = TRUE,
                         dig.lab = 2)]

  world <- ne_countries(scale = "medium", returnclass = "sf")
  world_att <- merge(world, att_g, by.x = "iso_a3", by.y = "country_iso3",
                     all.x = TRUE)

  # Color: negative bins red, positive bins blue, diverging outward from 0 so
  # the strongest effects are darkest and bins fade toward white at 0. Bins are
  # in increasing order, so negatives run most-negative -> near-zero (dark->light
  # red) and positives near-zero -> most-positive (light->dark blue).
  lev      <- levels(att_g$att_bin)
  bin_hi   <- as.numeric(sub(".*,([-0-9.eE+]+)\\]$", "\\1", lev))
  n_neg    <- sum(bin_hi <= 0)
  n_pos    <- length(lev) - n_neg
  reds  <- if (n_neg > 0)
    colorRampPalette(c("#67000d", "#fee5d9"))(n_neg) else character(0)
  blues <- if (n_pos > 0)
    colorRampPalette(c("#deebf7", "#08306b"))(n_pos) else character(0)
  pal <- c(reds, blues)

  p_map <- ggplot(world_att) +
    geom_sf(aes(fill = att_bin), color = "gray70", linewidth = 0.1) +
    scale_fill_manual(values = pal, na.value = "gray92",
                      name = "ATT at g\n(forest frac,\nJenks bins)",
                      drop = FALSE) +
    labs(title = "Protected-area effect on forest fraction at horizon g",
         subtitle = sprintf(paste0("Stacked SC ATT at g = min(max event time,",
                                   " %d) | winsorized %d%% | %d countries"),
                            G_CAP, round(100 * WINS_P), nrow(att_g))) +
    coord_sf(ylim = c(-40, 40), expand = FALSE) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
          plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0),
          legend.position = "right")

  ggsave(file.path(gdir, "att_choropleth_g.png"), p_map,
         width = 12, height = 6, dpi = 300, bg = "white")
  cat("  -> att_choropleth_g.png\n")

  # --- coefplot helper -------------------------------------------------------
  make_coefplot <- function(d, title, subtitle, fname) {
    d <- copy(d)
    d[, country_iso3 := factor(country_iso3, levels = d$country_iso3)]
    p <- ggplot(d, aes(att, country_iso3)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40",
                 linewidth = 0.5) +
      geom_linerange(aes(xmin = ci_lo, xmax = ci_hi), color = "steelblue",
                     linewidth = 0.7) +
      geom_point(color = "steelblue", size = 2) +
      labs(x = "ATT at g (forest fraction)", y = NULL,
           title = title, subtitle = subtitle) +
      base_theme
    ggsave(file.path(gdir, fname), p, width = 8,
           height = max(4, 0.28 * nrow(d) + 1.5), dpi = 300, bg = "white")
    cat("  ->", fname, "\n")
  }

  sub_g <- sprintf("Stacked SC ATT at g = min(max event time, %d) | Conley 50 km CIs",
                   G_CAP)

  # --- (2) Top 20 by area protected ------------------------------------------
  cat("(2) Coefplot: top 20 by area...\n")
  top_area <- att_g[order(-area)][seq_len(min(20, nrow(att_g)))][order(att)]
  make_coefplot(top_area, "ATT at g: 20 countries with the most protected area",
                sub_g, "att_coefplot_top20area.png")

  # --- (3) Top 10 + bottom 10 by ATT -----------------------------------------
  cat("(3) Coefplot: top 10 + bottom 10 by ATT...\n")
  ord <- att_g[order(att)]
  n <- nrow(ord)
  sel <- if (n <= 20) ord else rbind(ord[seq_len(10)], ord[(n - 9):n])
  sel <- unique(sel)[order(att)]
  make_coefplot(sel, "ATT at g: 10 largest and 10 smallest country effects",
                sub_g, "att_coefplot_top10bottom10.png")

  # --- (4) Density of per-PA treatment effects at g --------------------------
  # Each PA's effect = its tau at g_PA = min(that PA's max event time, G_CAP).
  cat("(4) Density of per-PA ATT at g...\n")
  pa_att <- dynamic[, {
    g_pa <- min(max(event_time), G_CAP)
    .(att = tau[event_time == g_pa][1])
  }, by = .(country_iso3, wdpa_pid)]
  pa_att <- pa_att[is.finite(att)]

  p_dens <- ggplot(pa_att, aes(att)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60",
               linewidth = 0.5) +
    geom_density(color = "steelblue", linewidth = 1.1) +
    labs(x = "Treatment effect at g (forest fraction)", y = "Density",
         title = "Distribution of protected-area treatment effects at g",
         subtitle = sprintf("Per-PA synthetic-control effect | %d PAs",
                            nrow(pa_att))) +
    base_theme

  ggsave(file.path(gdir, "att_density.png"), p_dens,
         width = 8, height = 5, dpi = 300, bg = "white")
  cat("  -> att_density.png\n")

  cat("ATT-at-g figures complete.\n")
}
