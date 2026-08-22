# ==============================================================================
# 5_protection_lpm.R
# Protected-area SELECTION: how parcel accessibility predicts forest USE, the
# implied "inclusive value" of a parcel, and a balance table comparing parcels
# that are ever vs. never protected.
#
# STEP 1 -- Accessibility regression (cross-section of cells, year 1995):
#   feols( log((1 - forest_share_1995) / forest_share_1995)
#            ~ log(pop_density_1990 + 1)
#            + log(min(travel_time_cities, travel_time_ports) + 1)
#            + log(pop_access_50km + 1)
#            + log(biomass_access_50km + 1)
#          | country_iso3 + is_interior,
#          weights = ~area_km2,
#          vcov    = Conley(50 km) )
#   LHS = log-odds that a parcel is in NON-forest ("forest use") in 1995.
#
# STEP 2 -- Inclusive value in minute-equivalents (linear log-minute metric):
#   V_i        = b_pop*l_pop + b_tt*l_tt + b_popacc*l_popacc + b_bio*l_bio
#   IV_minutes = V_i / (-b_tt)
#   (the four estimated slope coefficients applied to the parcel's accessibility
#    terms, excluding the country / interior FE; expressed on the travel-time
#    scale, i.e. log-minute equivalents.)
#   Berry parcel unobservable (from OBSERVABLES ONLY, not the FE):
#     xi = y_use - V_i           (so xi absorbs the FE + the true residual)
#   IV with unobservable = (V_i + xi) / (-b_tt) = y_use / (-b_tt), which matches
#   each parcel's observed log-odds exactly.
#
# STEP 3 -- Implied land-use tax: for each protected parcel, the additive logit
#   tax that reproduces its PA's last-observed synthetic-control treatment effect
#   (forest predicted = fc_1995 + TE exactly). Control parcels are taxed at their
#   country's area-weighted mean tax, yielding a counterfactual forest share and
#   an IMPLIED treatment effect. Treated-parcel reconstruction recovers TE up to
#   forest-share clamping (checked in-script).
#
# STEP 4 -- Balance table (AER booktabs), ever-protected vs. never-protected,
#   area-weighted, on each individual sub-variable, the 1990-95 forest trend, the
#   1995 and 1990 forest levels, the inclusive value, and the treatment effect
#   (treated actual / control implied). Difference + clustered (country) SE +
#   stars. Written to a .tex file.
#
# STEP 5 -- Country-level optimal tax in minute-equivalents (.tex/.csv) and a
#   density plot overlaying actual treated vs. implied control treatment effects.
#
# STEP 6 -- Decision model: LPM (feols) and probit (feglm) for the probability of
#   protection as a tradeoff between the treatment-effect benefit and the
#   (with-unobservable) inclusive-value opportunity cost, with country and
#   frontier FE added cumulatively. Coefficients logged to console and .tex.
#   Also heterogeneity by continent, zone, and zone x continent.
#
# STEP 7 -- Carbon-benefit spec: committed carbon density (tCO2/ha) = forest-share
#   gained x aboveground carbon density x 44/12. Re-runs the decision model and
#   reports WTP for carbon = -b_co2 / b_oc (log-min-equiv. per tCO2/ha).
#
# STEP 8 -- Country WTP for carbon: a noisy country-interacted LPM gives free
#   per-country (b_co2, b_oc); a normal-normal model shrinks each toward the
#   POOLED area-weighted prior (Step 7's lpm_c, which is also the reported global
#   average WTP). Outputs raw vs. shrunk per-country WTP.
#
# Inputs:
#   Data/build/final/TMF_5km_panel.parquet
# Outputs (output/figures/analysis/protection_lpm/):
#   accessibility_regression.tex   first-stage regression (basis for the IV)
#   balance_ever_vs_never.tex      publication-ready balance table
#
# Usage: Rscript code/analysis/5_protection_lpm.R
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(duckdb)
  library(DBI)
  library(fixest)
  library(ggplot2)
  library(patchwork)
  library(countrycode)
})

if (Sys.getenv("SLURM_SUBMIT_DIR") != "") {
  project_root <- Sys.getenv("SLURM_SUBMIT_DIR")
} else {
  project_root <- here::here()
}

panel_path <- file.path(project_root, "Data/build/final/TMF_5km_panel.parquet")
out_dir    <- file.path(project_root, "output/figures/analysis/protection_lpm")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

stopifnot(file.exists(panel_path))

CONLEY_CUTOFF <- 50  # km, matches 4b_pooled_event_study.R

# ------------------------------------------------------------------------------
# 1. Build the cell-level cross-section from the panel (DuckDB collapse)
# ------------------------------------------------------------------------------
# Static cell attributes are year-invariant (ANY_VALUE is safe); ever_protected
# is MAX over the is_protected flag; forest shares come from the 1990 and 1995
# rows. One streaming pass over ~67M rows -> ~2M cells.

cat("Collapsing panel to one row per cell...\n")
con <- dbConnect(duckdb())

query <- sprintf("
  SELECT
    grid_id,
    MAX(CASE WHEN is_protected THEN 1 ELSE 0 END)    AS ever_protected,
    ANY_VALUE(wdpa_pid)                              AS wdpa_pid,
    ANY_VALUE(country_iso3)                          AS country_iso3,
    ANY_VALUE(is_interior)                           AS is_interior,
    ANY_VALUE(is_frontier)                           AS is_frontier,
    ANY_VALUE(zone)                                  AS zone,
    ANY_VALUE(desig_year)                            AS desig_year,
    ANY_VALUE(area_km2)                              AS area_km2,
    ANY_VALUE(centroid_lat)                          AS centroid_lat,
    ANY_VALUE(centroid_lon)                          AS centroid_lon,
    ANY_VALUE(pop_density_1990)                      AS pop_density_1990,
    ANY_VALUE(travel_time_cities)                    AS travel_time_cities,
    ANY_VALUE(travel_time_ports)                     AS travel_time_ports,
    ANY_VALUE(pop_access_50km)                       AS pop_access_50km,
    ANY_VALUE(biomass_access_50km)                   AS biomass_access_50km,
    ANY_VALUE(aboveground_biomass_carbon_2010)       AS agc_2010,
    MAX(CASE WHEN year = 1990 THEN forest_cover END) AS fc_1990,
    MAX(CASE WHEN year = 1995 THEN forest_cover END) AS fc_1995
  FROM read_parquet('%s')
  GROUP BY grid_id", panel_path)

dt <- as.data.table(dbGetQuery(con, query))
dbDisconnect(con, shutdown = TRUE)
cat("Collapsed to", format(nrow(dt), big.mark = ","), "cells\n")

# Clamp forest share to (0.001, 0.999) for ALL analysis: keeps the log-odds
# finite at the 0/1 boundary (no cells dropped) and bounds the forest levels.
FS_LO <- 0.001
FS_HI <- 0.999
dt[, fc_1990 := pmin(pmax(fc_1990, FS_LO), FS_HI)]
dt[, fc_1995 := pmin(pmax(fc_1995, FS_LO), FS_HI)]

# The raw aboveground biomass-carbon raster is stored in 10 t/ha units; rescale
# to t/ha (tC/ha) so all downstream carbon math is in tonnes.
dt[, agc_2010 := agc_2010 * 10]

# ------------------------------------------------------------------------------
# 2. Construct regressors, outcome, and the forest trend
# ------------------------------------------------------------------------------

# Nearest market = min of the two travel times (NA only if both missing)
dt[, min_travel_time := pmin(travel_time_cities, travel_time_ports, na.rm = TRUE)]
dt[!is.finite(min_travel_time), min_travel_time := NA_real_]

# log(x + 1) transforms used by the regression and the inclusive value
dt[, l_pop    := log(pop_density_1990    + 1)]
dt[, l_tt     := log(min_travel_time     + 1)]
dt[, l_popacc := log(pop_access_50km     + 1)]
dt[, l_bio    := log(biomass_access_50km + 1)]

# Forest-share trend (1990 -> 1995) and the log-odds of NON-forest ("use") in 1995
dt[, forest_share_diff := fc_1995 - fc_1990]
dt[, y_use := log((1 - fc_1995) / fc_1995)]
dt[!is.finite(y_use), y_use := NA_real_]  # NA only where fc_1995 itself is missing

# ------------------------------------------------------------------------------
# 3. Accessibility regression: forest use ~ log accessibility | country + interior
# ------------------------------------------------------------------------------

cat("\n=== Accessibility regression (log-odds of 1995 forest use) ===\n")

est_sample <- dt[is.finite(y_use) & is.finite(area_km2) & area_km2 > 0 &
                 !is.na(centroid_lat) & !is.na(centroid_lon)]
cat("Estimation sample:", format(nrow(est_sample), big.mark = ","), "cells\n")

m_use <- feols(
  y_use ~ l_pop + l_tt + l_popacc + l_bio | country_iso3 + is_interior,
  data    = est_sample,
  weights = ~area_km2,
  vcov    = vcov_conley(lat = ~centroid_lat, lon = ~centroid_lon,
                        cutoff = CONLEY_CUTOFF)
)

reg_labels <- c(
  l_pop    = "log(Pop. density 1990 + 1)",
  l_tt     = "log(Min. travel time + 1)",
  l_popacc = "log(Pop. access 50km + 1)",
  l_bio    = "log(Biomass access 50km + 1)"
)
etable(m_use, dict = reg_labels)  # echo to console
etable(m_use, dict = reg_labels,
       file = file.path(out_dir, "accessibility_regression.tex"), replace = TRUE,
       title = paste0("Log-odds of 1995 forest use on parcel accessibility ",
                      "(country + interior-forest FE, Conley 50 km SEs)"),
       fitstat = ~ n + r2 + war2, digits = 4)

# ------------------------------------------------------------------------------
# 4. Inclusive value in minute-equivalents (linear log-minute metric)
# ------------------------------------------------------------------------------
# V_i = b_pop*l_pop + b_tt*l_tt + b_popacc*l_popacc + b_bio*l_bio  (slopes only,
# FE excluded); IV in minute-equivalents = V_i / (-b_tt). Dividing by the
# negative travel-time coefficient (b_tt < 0, so this is |b_tt|) orients the
# metric so that more-accessible parcels score higher.

b <- coef(m_use)
cat(sprintf("\nCoefficients: pop=%.4f  travel=%.4f  popacc=%.4f  bio=%.4f\n",
            b["l_pop"], b["l_tt"], b["l_popacc"], b["l_bio"]))

dt[, iv_index := b["l_pop"] * l_pop + b["l_tt"] * l_tt +
                 b["l_popacc"] * l_popacc + b["l_bio"] * l_bio]
dt[, iv_minutes := iv_index / (-b["l_tt"])]

# Berry inversion: the parcel unobservable is the observed log-odds minus the
# prediction from OBSERVABLES ONLY (not the FE), so xi captures the FE + the true
# residual. Adding it back matches each parcel exactly, hence the with-xi
# inclusive value equals y_use / (-b_tt).
dt[, xi := y_use - iv_index]
dt[, iv_minutes_xi := (iv_index + xi) / (-b["l_tt"])]

# ------------------------------------------------------------------------------
# 5. Implied land-use tax reproducing the synthetic-control treatment effect
# ------------------------------------------------------------------------------
# The logit maps an additive index shift to a forest-share change:
#   forest(idx) = 1 / (1 + exp(idx)),  with idx = y_use the log-odds of USE.
# A tax t > 0 lowers the use index (idx -> y_use - t), raising forest share.
# Using the parcel's OWN observed log-odds (= the with-unobservable inclusive
# value) guarantees an exact baseline match. For each protected parcel we solve
# for the t that moves baseline forest fc_1995 to fc_1995 + TE, where TE is the
# PA's LAST-OBSERVED synthetic-control effect (tau at the final event time):
#   t_parcel = y_use - logit_use( clamp(fc_1995 + TE) ),
# so the logit then predicts exactly TE + the 1995 land use.
# Controls are taxed at their country's area-weighted mean of t_parcel, giving a
# counterfactual forest share and an IMPLIED treatment effect.

forest_from_index <- function(idx) 1 / (1 + exp(idx))   # inverse of y_use
logit_use         <- function(f)   log((1 - f) / f)     # = y_use at f

dyn_path <- file.path(project_root,
                      "Data/build/intermediate/synth/synth_dynamic_effects.csv")
stopifnot(file.exists(dyn_path))
dyn <- fread(dyn_path)
dyn <- dyn[outcome == "forest_cover"]

# Last-observed treatment effect per PA = tau at the maximum (final) year
setorder(dyn, wdpa_pid, year)
te_pa <- dyn[, .(te_last = tau[.N], te_year = year[.N]),
             by = .(wdpa_pid, country_iso3)]
te_pa <- unique(te_pa, by = "wdpa_pid")
cat(sprintf("\nLast-observed SC treatment effects: %d PAs (mean tau %.4f)\n",
            nrow(te_pa), mean(te_pa$te_last, na.rm = TRUE)))

# Attach each PA's last treatment effect to its (protected) parcels
dt <- merge(dt, te_pa[, .(wdpa_pid, te_last)], by = "wdpa_pid", all.x = TRUE)

# Per-parcel implied tax on protected parcels that have an SC effect
dt[, f_star := NA_real_]
dt[ever_protected == 1 & !is.na(te_last) & is.finite(y_use),
   f_star := pmin(pmax(fc_1995 + te_last, FS_LO), FS_HI)]
dt[, tax_parcel := NA_real_]
dt[!is.na(f_star), tax_parcel := y_use - logit_use(f_star)]

n_clamped <- dt[!is.na(f_star) & (fc_1995 + te_last > FS_HI |
                                  fc_1995 + te_last < FS_LO), .N]
cat(sprintf("Implied tax computed on %s protected parcels (%s clamped at bounds)\n",
            format(dt[!is.na(tax_parcel), .N], big.mark = ","),
            format(n_clamped, big.mark = ",")))

# Country area-weighted average implied tax over treated parcels
tau_country <- dt[!is.na(tax_parcel),
                  .(tax_bar  = weighted.mean(tax_parcel, area_km2, na.rm = TRUE),
                    n_parcel = .N,
                    n_pa     = uniqueN(wdpa_pid)),
                  by = country_iso3]

# Apply the country tax to control (never-protected) parcels -> counterfactual
dt <- merge(dt, tau_country[, .(country_iso3, tax_bar)],
            by = "country_iso3", all.x = TRUE)
dt[, te_control := NA_real_]
dt[ever_protected == 0 & is.finite(y_use) & !is.na(tax_bar),
   te_control := forest_from_index(y_use - tax_bar) - fc_1995]

# Treated reconstruction (recovers te_last exactly except where f_star clamped)
dt[, te_treated_implied := NA_real_]
dt[!is.na(f_star), te_treated_implied := f_star - fc_1995]

# Combined effect for the balance table: treated -> ACTUAL SC effect (te_last);
# control -> IMPLIED effect under the country tax (te_control).
dt[, te_combined := fifelse(ever_protected == 1L, te_last, te_control)]

# Carbon-denominated benefit, PER HECTARE (tCO2/ha): forest-share gained x
# aboveground carbon density x CO2/C. Intensive (density) basis so cell area is
# NOT double-counted against the area weights. agc_2010 is Spawn (2020) carbon
# density (tC/ha), so NO 0.47 biomass->carbon factor.
dt[, co2_benefit_ha := te_combined * agc_2010 * (44 / 12)]

# --- Internal consistency check ----------------------------------------------
chk_actual  <- dt[ever_protected == 1 & !is.na(te_last),
                  weighted.mean(te_last, area_km2)]
chk_implied <- dt[!is.na(te_treated_implied),
                  weighted.mean(te_treated_implied, area_km2)]
chk_ctrl    <- dt[!is.na(te_control), weighted.mean(te_control, area_km2)]
max_recon_gap <- dt[!is.na(te_treated_implied),
                    max(abs(te_treated_implied - te_last))]
cat(sprintf(paste0("Consistency check (area-weighted):\n",
            "  treated ACTUAL last tau      = %.5f\n",
            "  treated IMPLIED (own tax)    = %.5f  (max parcel gap %.5f from clamping)\n",
            "  control IMPLIED (country tax)= %.5f\n"),
            chk_actual, chk_implied, max_recon_gap, chk_ctrl))

# ------------------------------------------------------------------------------
# 6. Balance table: ever-protected vs. never-protected (area-weighted)
# ------------------------------------------------------------------------------
# For each variable, a weighted OLS on the ever_protected dummy gives the
# never-protected mean (intercept), the difference (slope), and a country-
# clustered SE in one shot -- so means and the tested difference are consistent.

# (variable, display label, decimals)
bal_vars <- list(
  list("pop_density_1990",    "Population density (1990)",        2),
  list("travel_time_cities",  "Travel time to cities (min)",      1),
  list("travel_time_ports",   "Travel time to ports (min)",       1),
  list("pop_access_50km",     "Population access, 50km",          2),
  list("biomass_access_50km", "Biomass access, 50km",             2),
  list("forest_share_diff",   "Forest-share trend (1990--1995)",  4),
  list("fc_1995",             "Forest share (1995)",              3),
  list("fc_1990",             "Forest share (1990)",              3),
  list("iv_minutes",     "Inclusive value, obs. only (log-min-equiv.)",     3),
  list("iv_minutes_xi",  "Inclusive value, w/ unobservable (log-min-equiv.)", 3),
  list("te_combined",    "Treatment effect (treated actual / control implied)", 4)
)

stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) return("$^{***}$")
  if (p < 0.05) return("$^{**}$")
  if (p < 0.10) return("$^{*}$")
  ""
}

dbal <- dt[is.finite(area_km2) & area_km2 > 0 & !is.na(ever_protected)]

rows <- list()
for (v in bal_vars) {
  col <- v[[1]]; lab <- v[[2]]; dec <- v[[3]]
  fmt <- paste0("%.", dec, "f")

  m <- feols(as.formula(paste(col, "~ ever_protected")), data = dbal,
             weights = ~area_km2, cluster = ~country_iso3)
  cf        <- coef(m)
  mean_never <- cf["(Intercept)"]
  diff       <- cf["ever_protected"]
  mean_ever  <- mean_never + diff
  se         <- se(m)["ever_protected"]
  p          <- pvalue(m)["ever_protected"]

  rows[[length(rows) + 1]] <- sprintf("%s & %s & %s & %s%s \\\\",
    lab, sprintf(fmt, mean_ever), sprintf(fmt, mean_never),
    sprintf(fmt, diff), stars(p))
  rows[[length(rows) + 1]] <- sprintf(" & & & (%s) \\\\", sprintf(fmt, se))
}

n_ever  <- dbal[ever_protected == 1 & !is.na(area_km2), .N]
n_never <- dbal[ever_protected == 0 & !is.na(area_km2), .N]

tex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Parcel accessibility and forest by protection status}",
  "\\label{tab:balance_protection}",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  " & \\multicolumn{2}{c}{Area-weighted mean} & \\\\",
  "\\cmidrule(lr){2-3}",
  " & Ever protected & Never protected & Difference \\\\",
  "\\midrule",
  unlist(rows),
  "\\midrule",
  sprintf("Observations & %s & %s & \\\\",
          format(n_ever, big.mark = ","), format(n_never, big.mark = ",")),
  "\\bottomrule",
  "\\end{tabular}",
  "\\par\\vspace{0.5em}",
  "\\begin{minipage}{0.92\\linewidth}",
  paste0("\\footnotesize \\textit{Notes:} Cell-level cross-section, ",
         "area-weighted. The difference is the coefficient on an ever-protected ",
         "indicator from a weighted regression; standard errors (in parentheses) ",
         "are clustered by country. The inclusive value is the fitted ",
         "accessibility index $V_i$ divided by the negative of the travel-time ",
         "coefficient, in log-minute equivalents; the w/-unobservable version ",
         "adds the Berry parcel residual $\\xi_i = y_i - V_i$ (observables only, ",
         "FE excluded), matching each parcel's observed log-odds. ",
         "$^{*}$, $^{**}$, $^{***}$ denote $p<0.10$, $0.05$, $0.01$."),
  "\\end{minipage}",
  "\\end{table}"
)

bal_path <- file.path(out_dir, "balance_ever_vs_never.tex")
writeLines(tex, bal_path)

cat("\n=== Balance table (ever vs. never protected) ===\n")
cat(paste(tex, collapse = "\n"), "\n")

# ------------------------------------------------------------------------------
# 7. Country-level optimal tax (minute-equivalents) + TE distribution plot
# ------------------------------------------------------------------------------
# Express the country area-weighted implied tax on the travel-time scale, the
# same metric as the inclusive value: tax_minutes = tax_bar / (-b_tt).

tau_country[, tax_minutes := tax_bar / (-b["l_tt"])]
setorder(tau_country, -tax_minutes)
fwrite(tau_country, file.path(out_dir, "country_implied_tax.csv"))

tax_rows <- tau_country[, sprintf("%s & %d & %.4f & %.2f \\\\",
                                  country_iso3, n_pa, tax_bar, tax_minutes)]
tax_tex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Country-level implied optimal land-use tax}",
  "\\label{tab:country_tax}",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  "Country & PAs & Tax (log-odds) & Tax (log-min-equiv.) \\\\",
  "\\midrule",
  tax_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\par\\vspace{0.5em}",
  "\\begin{minipage}{0.92\\linewidth}",
  paste0("\\footnotesize \\textit{Notes:} Area-weighted average over a country's ",
         "protected parcels of the additive logit tax that reproduces each ",
         "parcel's last-observed synthetic-control treatment effect. The minute-",
         "equivalent divides by the negative travel-time coefficient."),
  "\\end{minipage}",
  "\\end{table}"
)
tax_path <- file.path(out_dir, "country_implied_tax.tex")
writeLines(tax_tex, tax_path)

# Density overlay: actual treated effects vs. implied control effects.
# Treated parcels are collapsed to one unit per PA (te_last is constant within a
# PA; weight = total PA area); controls stay per-parcel. Both area-weighted.
dens_treated <- dt[ever_protected == 1 & !is.na(te_last),
                   .(te = te_last[1], w = sum(area_km2)), by = wdpa_pid]
dens <- rbindlist(list(
  dens_treated[, .(te, w, grp = "Treated (actual SC effect)")],
  dt[!is.na(te_control),
     .(te = te_control, w = area_km2, grp = "Control (implied under tax)")]
))
dens[, w := w / sum(w), by = grp]
# The two groups have wildly different concentration (controls pile near 0,
# treated effects are spread), so area-normalized densities overlaid on a shared
# y-axis are unreadable -- the control spike (~5000) flattens everything. Peak-
# scale each density to 1 (after_stat(scaled)) to compare SHAPES, and set the
# x-window to the UNION of both groups' central 1-99% (a pooled quantile is
# swamped by the 1.8M controls and would clip the treated spread).
qs <- dens[, .(lo = quantile(te, 0.01, na.rm = TRUE),
               hi = quantile(te, 0.99, na.rm = TRUE)), by = grp]
xr_te <- c(min(qs$lo), max(qs$hi))

p_dens <- ggplot(dens, aes(te, y = after_stat(scaled), weight = w,
                           color = grp, fill = grp)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40",
             linewidth = 0.5) +
  geom_density(alpha = 0.25, linewidth = 1) +
  coord_cartesian(xlim = xr_te) +
  scale_color_manual(values = c("Treated (actual SC effect)" = "steelblue",
                                "Control (implied under tax)" = "coral")) +
  scale_fill_manual(values = c("Treated (actual SC effect)" = "steelblue",
                               "Control (implied under tax)" = "coral")) +
  labs(x = "Treatment effect on forest share", y = "Relative density (peak = 1)",
       color = NULL, fill = NULL,
       title = "Treatment effects: treated (actual) vs. controls (implied under tax)",
       subtitle = paste0("Each density peak-scaled to 1 (controls concentrate ",
                          "near 0). Treated: per PA; controls: per parcel; ",
                          "area-weighted.")) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"), legend.position = "bottom")

dens_path <- file.path(out_dir, "te_distribution_treated_vs_control.png")
ggsave(dens_path, p_dens, width = 9, height = 5.5, dpi = 300, bg = "white")

# --- Path of treatment effect & inclusive value by designation year ----------
# Area-weighted mean within each designation year, on the SC-estimated treated
# parcels. Two stacked column panels share the year axis; each keeps its own
# y-units (forest share vs. minute-equivalents).
dpath <- dt[ever_protected == 1 & !is.na(desig_year) & !is.na(te_last) &
            is.finite(iv_minutes_xi) & is.finite(area_km2) & area_km2 > 0]
yr_agg <- dpath[, .(te_mean = weighted.mean(te_last, area_km2),
                    iv_mean = weighted.mean(iv_minutes_xi, area_km2)),
                by = desig_year][order(desig_year)]

yr_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"))

p_te_yr <- ggplot(yr_agg, aes(desig_year, te_mean)) +
  geom_col(fill = "steelblue", width = 0.8) +
  labs(x = NULL, y = "Mean treatment effect\n(forest share)",
       title = "Treatment effect and opportunity cost by designation year") +
  yr_theme + theme(axis.text.x = element_blank())

p_iv_yr <- ggplot(yr_agg, aes(desig_year, iv_mean)) +
  geom_col(fill = "coral", width = 0.8) +
  labs(x = "Designation year", y = "Mean inclusive value\n(log-min-equiv.)") +
  yr_theme

yr_path <- file.path(out_dir, "te_iv_by_desig_year.png")
ggsave(yr_path, p_te_yr / p_iv_yr, width = 9, height = 6.5, dpi = 300,
       bg = "white")

# ------------------------------------------------------------------------------
# 8. Decision model: probit + LPM for the probability of protection
# ------------------------------------------------------------------------------
# The decision maker trades off the conservation BENEFIT (the treatment effect,
# now filled in for every parcel: actual for treated, implied for controls)
# against the OPPORTUNITY COST of foregone land use (the logit inclusive value).
#   P(protected) = f( b_te * treatment_effect + b_oc * inclusive_value )
# Expect b_te > 0 (protect where protection works) and b_oc < 0 (avoid high
# opportunity-cost parcels). Estimated as an LPM (feols) and a probit (feglm),
# each with the obs.-only and the with-unobservable inclusive value, area-
# weighted, country-clustered.

cat("\n=== Decision model: P(protection) ~ treatment effect + opportunity cost ===\n")

ddec <- dt[!is.na(ever_protected) & is.finite(area_km2) & area_km2 > 0 &
           !is.na(te_combined) & !is.na(centroid_lat) & !is.na(centroid_lon)]
cat(sprintf("Decision-model sample: %s parcels (%s protected)\n",
            format(nrow(ddec), big.mark = ","),
            format(ddec[ever_protected == 1, .N], big.mark = ",")))

# With-unobservable inclusive value only; FE added cumulatively
# (none -> country -> country + frontier), for both LPM and probit.
lpm1 <- feols(ever_protected ~ te_combined + iv_minutes_xi, ddec,
              weights = ~area_km2,
              vcov = vcov_conley(lat = ~centroid_lat, lon = ~centroid_lon,
                                 cutoff = CONLEY_CUTOFF))
lpm2 <- feols(ever_protected ~ te_combined + iv_minutes_xi | country_iso3, ddec,
              weights = ~area_km2,
              vcov = vcov_conley(lat = ~centroid_lat, lon = ~centroid_lon,
                                 cutoff = CONLEY_CUTOFF))
lpm3 <- feols(ever_protected ~ te_combined + iv_minutes_xi |
              country_iso3 + is_frontier, ddec,
              weights = ~area_km2,
              vcov = vcov_conley(lat = ~centroid_lat, lon = ~centroid_lon,
                                 cutoff = CONLEY_CUTOFF))
prb1 <- feglm(ever_protected ~ te_combined + iv_minutes_xi, ddec,
              family = binomial("probit"),
              weights = ~area_km2,
              vcov = vcov_conley(lat = ~centroid_lat, lon = ~centroid_lon,
                                 cutoff = CONLEY_CUTOFF))
prb2 <- feglm(ever_protected ~ te_combined + iv_minutes_xi | country_iso3, ddec,
              family = binomial("probit"),
              weights = ~area_km2,
              vcov = vcov_conley(lat = ~centroid_lat, lon = ~centroid_lon,
                                 cutoff = CONLEY_CUTOFF))
prb3 <- feglm(ever_protected ~ te_combined + iv_minutes_xi |
              country_iso3 + is_frontier, ddec,
              family = binomial("probit"),
              weights = ~area_km2,
              vcov = vcov_conley(lat = ~centroid_lat, lon = ~centroid_lon,
                                 cutoff = CONLEY_CUTOFF))

dec_labels <- c(
  te_combined    = "Treatment effect",
  co2_benefit_ha = "Carbon benefit (tCO2/ha)",
  iv_minutes_xi  = "Opportunity cost (inclusive value)"
)
dec_hdr <- c("LPM", "LPM", "LPM", "Probit", "Probit", "Probit")

# Console log + LaTeX table
etable(lpm1, lpm2, lpm3, prb1, prb2, prb3, dict = dec_labels, headers = dec_hdr)
choice_path <- file.path(out_dir, "protection_choice_probit_lpm.tex")
etable(lpm1, lpm2, lpm3, prb1, prb2, prb3, dict = dec_labels, headers = dec_hdr,
       file = choice_path, replace = TRUE,
       title = paste0("Probability of protection: treatment-effect benefit vs. ",
                      "inclusive-value opportunity cost (area-weighted, ",
                      "Conley 50 km SEs; cumulative country/frontier FE)"),
       digits = 4)

# --- Carbon-benefit spec + implied WTP for carbon -----------------------------
# Re-run the full-FE decision model with the carbon benefit (tCO2/ha) replacing
# the forest-share benefit. Implied WTP for carbon = -b_co2 / b_oc, the
# opportunity cost (log-min-equiv.) the planner trades per tCO2/ha (delta SE).
lpm_c <- feols(ever_protected ~ co2_benefit_ha + iv_minutes_xi |
               country_iso3 + is_frontier, ddec, weights = ~area_km2,
               vcov = vcov_conley(lat = ~centroid_lat, lon = ~centroid_lon,
                                  cutoff = CONLEY_CUTOFF))
prb_c <- feglm(ever_protected ~ co2_benefit_ha + iv_minutes_xi |
               country_iso3 + is_frontier, ddec, family = binomial("probit"),
               weights = ~area_km2,
               vcov = vcov_conley(lat = ~centroid_lat, lon = ~centroid_lon,
                                  cutoff = CONLEY_CUTOFF))

cat("\n=== Decision model: forest-share vs. carbon benefit ===\n")
etable(lpm3, lpm_c, prb3, prb_c, dict = dec_labels,
       headers = c("LPM forest", "LPM carbon", "Probit forest", "Probit carbon"))
carbon_path <- file.path(out_dir, "protection_choice_carbon.tex")
etable(lpm3, lpm_c, prb3, prb_c, dict = dec_labels,
       headers = c("LPM forest", "LPM carbon", "Probit forest", "Probit carbon"),
       file = carbon_path, replace = TRUE,
       title = paste0("Probability of protection: forest-share vs. carbon ",
                      "benefit (country + frontier FE, Conley 50 km SEs)"),
       digits = 4)

wtp_delta <- function(m) {
  v <- c("co2_benefit_ha", "iv_minutes_xi")
  b <- coef(m)[v]; Vd <- vcov(m)[v, v]
  grad <- c(-1 / b[2], b[1] / b[2]^2)        # d(-b1/b2)
  c(wtp = unname(-b[1] / b[2]),
    se  = sqrt(as.numeric(t(grad) %*% Vd %*% grad)))
}
wtp_lpm <- wtp_delta(lpm_c); wtp_prb <- wtp_delta(prb_c)
cat(sprintf(paste0("Implied WTP for carbon (log-min-equiv. per tCO2/ha):\n",
            "  LPM    = %.4f (se %.4f)\n  Probit = %.4f (se %.4f)\n"),
            wtp_lpm["wtp"], wtp_lpm["se"], wtp_prb["wtp"], wtp_prb["se"]))
wtp_csv <- file.path(out_dir, "carbon_wtp_global.csv")
fwrite(data.table(model = c("LPM", "Probit"),
                  wtp_logmin_per_tco2ha = c(wtp_lpm["wtp"], wtp_prb["wtp"]),
                  se                    = c(wtp_lpm["se"],  wtp_prb["se"])),
       wtp_csv)

# --- Warm glow: probit residual expressed in travel-time (minute) terms -------
# The residual of the full probit (prb3) is the protection propensity NOT
# explained by the treatment-effect benefit, the inclusive-value cost, or the
# country/frontier FE -- an intrinsic "warm glow" of protection. Dividing by the
# IV coefficient (b_oc, utility per log-minute-equivalent) re-expresses it in the
# same travel-time metric as the opportunity cost.
b_oc <- coef(prb3)["iv_minutes_xi"]
stopifnot(is.finite(b_oc))
# predict() returns NA for rows whose country/frontier FE level was dropped by
# the probit (perfect-separation groups); keep only the finite residuals.
ddec[, prob_resid := ever_protected - predict(prb3, newdata = ddec,
                                              type = "response")]
ddec[, warm_glow_min := prob_resid / b_oc]
dwg <- ddec[is.finite(warm_glow_min) & is.finite(area_km2) & area_km2 > 0]
cat(sprintf("Warm glow computed on %s of %s parcels (%s dropped: FE separation)\n",
            format(nrow(dwg), big.mark = ","), format(nrow(ddec), big.mark = ","),
            format(nrow(ddec) - nrow(dwg), big.mark = ",")))

wg_mean <- weighted.mean(dwg$warm_glow_min, dwg$area_km2)
wg_sd   <- sqrt(weighted.mean((dwg$warm_glow_min - wg_mean)^2, dwg$area_km2))
cat(sprintf(paste0("Warm glow (probit residual / IV coef), travel-time terms:",
            "\n  mean = %.2f log-min-equiv,  sd = %.2f log-min-equiv\n"),
            wg_mean, wg_sd))

xr_wg <- quantile(dwg$warm_glow_min, c(0.01, 0.99), na.rm = TRUE)
p_warm <- ggplot(dwg, aes(warm_glow_min,
                          weight = area_km2 / sum(area_km2))) +
  geom_density(aes(y = after_stat(density)), fill = "steelblue", alpha = 0.3,
               color = "steelblue", linewidth = 1) +
  geom_vline(xintercept = wg_mean, linetype = "dashed", color = "coral",
             linewidth = 0.8) +
  coord_cartesian(xlim = xr_wg) +
  labs(x = "Warm glow (log-min-equiv.)", y = "Density (area-weighted)",
       title = "Warm glow: unexplained protection propensity in travel-time terms",
       subtitle = sprintf(paste0("Probit (prb3) residual / IV coef | ",
                          "mean %.1f, sd %.1f log-min-equiv"), wg_mean, wg_sd)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
warm_path <- file.path(out_dir, "warm_glow_distribution.png")
ggsave(warm_path, p_warm, width = 9, height = 5, dpi = 300, bg = "white")

# --- Probit coefficients broken out by continent (3 groups) -------------------
ddec[, continent := countrycode(country_iso3, "iso3c", "continent", warn = FALSE)]
ddec[, continent := fcase(
  continent == "Americas",              "Americas",
  continent == "Africa",                "Africa",
  continent %in% c("Asia", "Oceania"),  "Asia/Oceania",
  default = NA_character_)]
cont_levels <- c("Americas", "Africa", "Asia/Oceania")

prb_cont <- lapply(cont_levels, function(cc) {
  feglm(ever_protected ~ te_combined + iv_minutes_xi | country_iso3 + is_frontier,
        data = ddec[continent == cc], family = binomial("probit"),
        weights = ~area_km2,
        vcov = vcov_conley(lat = ~centroid_lat, lon = ~centroid_lon,
                           cutoff = CONLEY_CUTOFF))
})

cat("\n=== Probit by continent (Americas / Africa / Asia-Oceania) ===\n")
etable(prb_cont, dict = dec_labels, headers = cont_levels)
cont_path <- file.path(out_dir, "protection_choice_probit_by_continent.tex")
etable(prb_cont, dict = dec_labels, headers = cont_levels,
       file = cont_path, replace = TRUE,
       title = paste0("Probability of protection by continent (probit, ",
                      "country + frontier FE, Conley 50 km SEs)"),
       digits = 4)

# --- Probit by zone, and by zone x continent ---------------------------------
# Conditioning on zone makes the frontier FE collinear, so these drop it and use
# country FE only (te + IV otherwise unchanged, area-weighted, Conley 50 km).
zone_levels <- c("frontier", "interior", "other")

fit_zone <- function(d) {
  if (nrow(d) < 50L || d[, uniqueN(ever_protected)] < 2L) return(NULL)
  tryCatch(
    feglm(ever_protected ~ te_combined + iv_minutes_xi | country_iso3,
          data = d, family = binomial("probit"), weights = ~area_km2,
          vcov = vcov_conley(lat = ~centroid_lat, lon = ~centroid_lon,
                             cutoff = CONLEY_CUTOFF)),
    error = function(e) NULL)
}

# (1) Zone only
zmods <- lapply(zone_levels, function(z) fit_zone(ddec[zone == z]))
zok   <- !vapply(zmods, is.null, logical(1))
cat("\n=== Probit by zone (frontier / interior / other) ===\n")
etable(zmods[zok], dict = dec_labels, headers = zone_levels[zok])
zone_path <- file.path(out_dir, "protection_choice_probit_by_zone.tex")
etable(zmods[zok], dict = dec_labels, headers = zone_levels[zok],
       file = zone_path, replace = TRUE,
       title = paste0("Probability of protection by zone (probit, country FE, ",
                      "Conley 50 km SEs)"), digits = 4)

# (2) Zone x continent: one table per zone, continents as columns
zonecont_path <- file.path(out_dir,
                           "protection_choice_probit_by_zone_continent.tex")
first_written <- FALSE
for (z in zone_levels) {
  cmods <- lapply(cont_levels, function(cc) fit_zone(ddec[zone == z &
                                                          continent == cc]))
  cok <- !vapply(cmods, is.null, logical(1))
  if (!any(cok)) next
  cat(sprintf("\n=== Probit zone = %s by continent ===\n", z))
  etable(cmods[cok], dict = dec_labels, headers = cont_levels[cok])
  etable(cmods[cok], dict = dec_labels, headers = cont_levels[cok],
         file = zonecont_path, replace = !first_written,
         title = sprintf(paste0("Probability of protection: zone = %s by ",
                         "continent (probit, country FE, Conley 50 km SEs)"), z),
         digits = 4)
  first_written <- TRUE
}

# ------------------------------------------------------------------------------
# 11. Country WTP for carbon: noisy interacted LPM + multivariate EB shrinkage
# ------------------------------------------------------------------------------
# One area-weighted LPM with free per-country slopes on the carbon benefit and
# the opportunity cost (country FE for intercepts, HC1 SEs). Extract each
# country's (b_co2, b_oc) and its 2x2 sampling covariance, then shrink toward the
# POOLED prior (mean and precision = the pooled area-weighted model, lpm_c).
# The global average WTP that we report is the pooled model's WTP, not anything
# re-aggregated here. Per-country WTP = -b_co2 / b_oc (log-min-equiv. per tCO2/ha).

cat("\n=== Country WTP for carbon: interacted LPM + EB shrinkage ===\n")
elig <- ddec[is.finite(co2_benefit_ha),
             .(n = .N, npro = sum(ever_protected == 1),
               ncon = sum(ever_protected == 0)), by = country_iso3
            ][n >= 100 & npro >= 5 & ncon >= 5, country_iso3]
ddec_b <- ddec[country_iso3 %in% elig & is.finite(co2_benefit_ha)]
cat(sprintf("Countries estimated: %d (of %d); parcels: %s\n",
            length(elig), uniqueN(ddec$country_iso3),
            format(nrow(ddec_b), big.mark = ",")))

mi <- feols(ever_protected ~ i(country_iso3, co2_benefit_ha) +
            i(country_iso3, iv_minutes_xi) | country_iso3,
            data = ddec_b, weights = ~area_km2, vcov = "hetero")

cf <- coef(mi); Vall <- vcov(mi); nm <- names(cf)
is_co2 <- grepl("co2_benefit_ha", nm) & grepl("country_iso3", nm)
is_iv  <- grepl("iv_minutes_xi",  nm) & grepl("country_iso3", nm)
cc_of  <- function(x) sub(".*country_iso3::([^:]+).*", "\\1", x)
co2_map <- setNames(which(is_co2), cc_of(nm[is_co2]))
iv_map  <- setNames(which(is_iv),  cc_of(nm[is_iv]))
ctrys <- intersect(names(co2_map), names(iv_map))

Theta <- matrix(NA_real_, length(ctrys), 2,
                dimnames = list(ctrys, c("b_co2", "b_oc")))
Vlist <- vector("list", length(ctrys))
keep  <- logical(length(ctrys))
for (i in seq_along(ctrys)) {
  idx <- c(co2_map[ctrys[i]], iv_map[ctrys[i]])
  th  <- cf[idx]; Vc <- Vall[idx, idx, drop = FALSE]
  if (all(is.finite(th)) && all(is.finite(Vc))) {
    Theta[i, ] <- th; Vlist[[i]] <- Vc + diag(1e-12, 2); keep[i] <- TRUE
  }
}
Theta <- Theta[keep, , drop = FALSE]; Vlist <- Vlist[keep]; ctrys <- ctrys[keep]

# Shrink each country toward the POOLED prior. The prior mean is the pooled,
# area-weighted (Conley) coefficients; the prior precision is that model's
# coefficient precision. The prior is held FIXED (not learned from the country
# estimates). Closed-form normal-normal posterior:
#   theta~_c = (P0 + Vc^-1)^-1 (P0 mu0 + Vc^-1 theta^_c).
prior_v <- c("co2_benefit_ha", "iv_minutes_xi")
mu0 <- coef(lpm_c)[prior_v]
P0  <- solve(vcov(lpm_c)[prior_v, prior_v])

shrunk <- t(vapply(seq_along(ctrys), function(i) {
  Pc <- solve(Vlist[[i]])
  as.numeric(solve(P0 + Pc) %*% (P0 %*% mu0 + Pc %*% Theta[i, ]))
}, numeric(2)))
dimnames(shrunk) <- list(ctrys, c("b_co2", "b_oc"))

wtp_raw    <- -Theta[, "b_co2"]  / Theta[, "b_oc"]
wtp_shrunk <- -shrunk[, "b_co2"] / shrunk[, "b_oc"]

# Post-shrinkage SE of each country's WTP: delta method on the posterior
# covariance Sig_c = (P0 + Vc^-1)^-1 of the shrunk (b_co2, b_oc).
wtp_shrunk_se <- vapply(seq_along(ctrys), function(i) {
  Sig <- solve(P0 + solve(Vlist[[i]]))
  b1 <- shrunk[i, "b_co2"]; b2 <- shrunk[i, "b_oc"]
  g  <- c(-1 / b2, b1 / b2^2)
  sqrt(as.numeric(t(g) %*% Sig %*% g))
}, numeric(1))

# The GLOBAL average WTP IS the pooled (area-weighted) model's WTP -- the prior
# itself -- computed in the carbon-spec step above (wtp_lpm). Not derived here.
wtp_global <- unname(wtp_lpm["wtp"])
cat(sprintf(paste0("Global WTP for carbon = pooled area-weighted prior: ",
            "%.4f log-min-equiv per tCO2/ha\n"), wtp_global))
cat(sprintf("Per-country b_oc < 0: %d/%d (shrunk), %d/%d (raw)\n",
            sum(shrunk[, "b_oc"] < 0), nrow(shrunk),
            sum(Theta[, "b_oc"] < 0), nrow(Theta)))

wtp_dt <- data.table(country_iso3 = ctrys,
                     b_co2_raw = Theta[, "b_co2"], b_oc_raw = Theta[, "b_oc"],
                     b_co2_shrunk = shrunk[, "b_co2"],
                     b_oc_shrunk  = shrunk[, "b_oc"],
                     wtp_raw = wtp_raw, wtp_shrunk = wtp_shrunk,
                     wtp_shrunk_se = wtp_shrunk_se)
shrink_csv <- file.path(out_dir, "country_carbon_wtp.csv")
fwrite(wtp_dt, shrink_csv)

# Caterpillar: raw -> shrunk WTP per country (winsorized for display)
wq <- quantile(c(wtp_raw, wtp_shrunk), c(0.02, 0.98), na.rm = TRUE)
pl <- copy(wtp_dt)
pl[, `:=`(wr = pmin(pmax(wtp_raw, wq[1]), wq[2]),
          ws = pmin(pmax(wtp_shrunk, wq[1]), wq[2]))]
pl[, country_iso3 := factor(country_iso3,
                            levels = pl[order(ws), country_iso3])]
p_shrink <- ggplot(pl) +
  geom_vline(xintercept = wtp_global, linetype = "dashed", color = "gray40") +
  geom_segment(aes(x = wr, xend = ws, y = country_iso3, yend = country_iso3),
               color = "gray70") +
  geom_point(aes(wr, country_iso3, color = "Raw"), size = 1.6) +
  geom_point(aes(ws, country_iso3, color = "Shrunk"), size = 1.6) +
  scale_color_manual(values = c(Raw = "coral", Shrunk = "steelblue"),
                     name = NULL) +
  labs(x = "WTP for carbon (log-min-equiv. per tCO2/ha)", y = NULL,
       title = "Country WTP for carbon: raw vs. shrunk to pooled prior",
       subtitle = sprintf(paste0("Normal-normal shrinkage toward the pooled ",
                          "area-weighted prior | %d countries | dashed = pooled"),
                          nrow(pl))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"), legend.position = "bottom")
shrink_plot <- file.path(out_dir, "country_carbon_wtp_caterpillar.png")
ggsave(shrink_plot, p_shrink, width = 8,
       height = max(5, 0.16 * nrow(pl) + 1.5), dpi = 300, bg = "white")

# --- Bar figure: global + per-country shrunk WTP, 95% CI (per tCO2) -----------
# Units: log-min-equiv. per tCO2/ha. Point + CI ends winsorized to the 5/95% of
# country WTPs for display; per-country CIs are the post-shrinkage 95% intervals,
# the global bar uses the pooled delta-method CI.
bar_dt <- rbind(
  data.table(country = "Global (pooled prior)",
             wtp = wtp_global, se = unname(wtp_lpm["se"]),
             is_global = TRUE),
  data.table(country = ctrys,
             wtp = wtp_shrunk, se = wtp_shrunk_se,
             is_global = FALSE)
)
bar_dt[, `:=`(ci_lo = wtp - 1.96 * se, ci_hi = wtp + 1.96 * se)]
wq2 <- quantile(bar_dt[is_global == FALSE, wtp], c(0.05, 0.95), na.rm = TRUE)
clip <- function(x) pmin(pmax(x, wq2[1]), wq2[2])
bar_dt[, `:=`(wtp_w = clip(wtp), lo_w = clip(ci_lo), hi_w = clip(ci_hi))]
ord <- bar_dt[is_global == FALSE][order(wtp), country]
bar_dt[, country := factor(country, levels = c(ord, "Global (pooled prior)"))]
bar_dt[, fill_cat := fifelse(is_global, "Global",
                             fifelse(wtp < 0, "Negative", "Positive"))]

p_bar <- ggplot(bar_dt, aes(wtp_w, country, fill = fill_cat)) +
  geom_col(width = 0.75) +
  geom_errorbarh(aes(xmin = lo_w, xmax = hi_w), height = 0.3,
                 color = "gray30", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c(Negative = "#bf7b86", Positive = "#2e6e8e",
                               Global = "#222222"), guide = "none") +
  labs(x = "WTP for CO2 (log-min-equiv. per tCO2/ha)", y = NULL,
       title = "Willingness to pay for carbon by country",
       subtitle = paste0("Shrunk to pooled area-weighted prior | 95% CI | ",
                         "winsorized 5/95% for display")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 8))
wtp_bar_path <- file.path(out_dir, "country_carbon_wtp_bars.png")
ggsave(wtp_bar_path, p_bar, width = 8,
       height = max(5, 0.18 * nrow(bar_dt) + 1.5), dpi = 300, bg = "white")

cat("\nWrote:\n  ", bal_path,
    "\n  ", file.path(out_dir, "accessibility_regression.tex"),
    "\n  ", tax_path,
    "\n  ", file.path(out_dir, "country_implied_tax.csv"),
    "\n  ", dens_path,
    "\n  ", yr_path,
    "\n  ", choice_path,
    "\n  ", warm_path,
    "\n  ", cont_path,
    "\n  ", zone_path,
    "\n  ", zonecont_path,
    "\n  ", carbon_path,
    "\n  ", wtp_csv,
    "\n  ", shrink_csv,
    "\n  ", shrink_plot,
    "\n  ", wtp_bar_path, "\n")
