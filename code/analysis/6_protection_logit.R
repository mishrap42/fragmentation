# ==============================================================================
# 6_protection_logit.R
# Static pooled logit: ever-protected on a $-denominated agricultural pressure
# measure built from EarthStat crop shares and yields, FAO producer prices, and
# a cell-specific transport wedge.
#
# THE MEASURE
# -----------
#   pressure_r = sum_c  s_cr * (p_c - tau_r) * y_cr        [USD per ha of cell]
#     s_cr  EarthStat cropshare_<crop>   harvested area / CELL area
#     y_cr  EarthStat yield_<crop>       t/ha on that crop's harvested area
#     p_c   FAO producer price, country x crop, 1991-2000 mean, USD/t
#     tau_r = tau_country * travel_time_cities_r / 60      [USD/t]
#
#   tau_r is a per-TON wedge common to all crops in the cell; its per-HECTARE
#   bite scales with each crop's own yield. Factoring, the bracket is
#   (p_c - tau_r) * y_cr, so a cell far from market can carry NEGATIVE pressure -
#   the haul costs more than the crop fetches. That is meaningful, not an error,
#   and it is much of the deep-forest sample.
#
# THREE WEIGHTING VARIANTS
#   pressure_cell   s_cr = cropshare_cr            USD per ha of CELL   [headline]
#   pressure_crop   s_cr = cropshare_cr / total    USD per ha of CROPLAND
#   pressure_max    max_c (p_c - tau_r) * y_cr     best-use frontier value
#
# WHAT THIS SPECIFICATION CAN AND CANNOT SAY
# ------------------------------------------
# y_cr is EarthStat ACTUAL yield, which is ZERO wherever the crop is not grown
# (EarthStat fills 0, not NA - confirmed by the Stage 0c validator: 100% non-NA,
# min 0). So pressure_r is ~0 on undisturbed forest and positive on the cleared
# frontier. It measures DEMONSTRATED agricultural value, not agronomic
# potential, and it conflates "this land is unsuitable" with "this land has not
# been reached yet". The share of the estimation sample sitting at exactly zero
# is reported below and should be read BEFORE the coefficient is.
#
# The clean fix is GAEZ potential yield, defined on land regardless of what is
# planted. Deliberately out of scope here.
#
# Inputs:  Data/build/final/TMF_5km_panel.parquet
#          Data/build/final/TMF_5km_cropland.parquet
#          Data/lookup/{crop_price_preperiod,trucking_cost}.parquet
#          Data/lookup/crop_crosswalk_earthstat_fao.csv
# Outputs: output/figures/analysis/protection_logit/
#            protection_logit.tex        spec ladder
#            pressure_summary.txt        measure diagnostics
#
# Usage: Rscript code/analysis/6_protection_logit.R
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(duckdb)
  library(DBI)
  library(arrow)
  library(fixest)
})

here::i_am('code/analysis/6_protection_logit.R')
source("code/build/BUILD_workspace.R")

output_dir <- file.path(project_root, "output", "figures", "analysis",
                        "protection_logit")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

FOREST_BASELINE_YEAR <- 1990
FOREST_THRESHOLD     <- 0.5     # undisturbed share to count as a candidate
CONLEY_KM            <- 50      # matches 5_protection_lpm.R

start_time <- Sys.time()
log_message("========================================")
log_message("PROTECTION LOGIT: $-denominated pressure")
log_message("========================================")

# ==============================================================================
# 1. CELL-LEVEL PANEL COLLAPSE
# ==============================================================================
# Arrow cannot open this file - the footer trips thrift's size limit - so the
# collapse goes through DuckDB, same as 5_protection_lpm.R.

panel_path <- file.path(final_output_path, "TMF_5km_panel.parquet")
stopifnot(file.exists(panel_path))

log_message("Collapsing panel to one row per cell...")
con <- dbConnect(duckdb())

query <- sprintf("
  SELECT
    grid_id,
    MAX(CASE WHEN is_protected THEN 1 ELSE 0 END)        AS ever_protected,
    ANY_VALUE(country_iso3)                              AS country_iso3,
    ANY_VALUE(area_km2)                                  AS area_km2,
    ANY_VALUE(centroid_lat)                              AS centroid_lat,
    ANY_VALUE(centroid_lon)                              AS centroid_lon,
    ANY_VALUE(desig_year)                                AS desig_year,
    ANY_VALUE(is_frontier)                               AS is_frontier,
    ANY_VALUE(travel_time_cities)                        AS travel_time_cities,
    ANY_VALUE(slope_degrees)                             AS slope_degrees,
    ANY_VALUE(elevation_m)                               AS elevation_m,
    ANY_VALUE(terrain_ruggedness)                        AS terrain_ruggedness,
    ANY_VALUE(pop_density_1990)                          AS pop_density_1990,
    ANY_VALUE(aboveground_biomass_carbon_2010)           AS agc_2010,
    ANY_VALUE(biodiv_combined_thr_sr_2022)               AS biodiv_thr,
    ANY_VALUE(cattle_density_2010_da)                    AS cattle_density,
    MAX(CASE WHEN year = %d THEN Undisturbed_TMF END)    AS undist_base
  FROM read_parquet('%s')
  GROUP BY grid_id", FOREST_BASELINE_YEAR, panel_path)

dt <- as.data.table(dbGetQuery(con, query))
dbDisconnect(con, shutdown = TRUE)
log_message(sprintf("  %s cells, %d countries",
                    format(nrow(dt), big.mark = ","), uniqueN(dt$country_iso3)))

# ==============================================================================
# 2. CROPLAND CROSS-SECTION
# ==============================================================================

stopifnot(file.exists(cropland_final_file))
xw <- fread(crop_crosswalk_file)

avail <- arrow::open_dataset(cropland_final_file, format = "parquet")$schema$names
xw <- xw[cropshare_col %in% avail & yield_col %in% avail]
log_message(sprintf("Crops with both layers present: %d", nrow(xw)))
stopifnot(nrow(xw) > 0)

need <- intersect(c("grid_id", "cropland_frac", "pasture_frac",
                    "cropshare_total", xw$cropshare_col, xw$yield_col), avail)

cs <- setDT(arrow::read_parquet(cropland_final_file,
                                col_select = arrow::all_of(need)))
log_message(sprintf("  cropland cross-section: %s cells, %d columns",
                    format(nrow(cs), big.mark = ","), ncol(cs)))

dt <- merge(dt, cs, by = "grid_id", all.x = TRUE)
rm(cs); gc()

n_nocrop <- dt[is.na(cropshare_total), .N]
log_message(sprintf("  cells with no cropland row: %s (%.2f%%)",
                    format(n_nocrop, big.mark = ","), 100 * n_nocrop / nrow(dt)))
if (n_nocrop > 0.05 * nrow(dt)) {
  log_message("  WARNING: >5% of cells unmatched. Did Stage 1b run AFTER 0c finished?")
}

# ==============================================================================
# 3. PRICES AND TRANSPORT COST
# ==============================================================================

# GADM assigns Z0x codes to disputed areas, and they carry no ISO3 economic
# data. Z07 is labelled "India" and holds 2,985 forested cells; without a remap
# it would take the global-median price AND the global-median tau despite India
# having both observed. XPI / XSP / XCL are uninhabited islands (43 cells total,
# no agriculture) and keep the fallback - not worth a judgement call.
# country_econ is used ONLY for the price and tau joins; the fixed effect stays
# on country_iso3 so the disputed area is not silently folded into India.
GADM_ECON_REMAP <- c(Z07 = "IND")
dt[, country_econ := country_iso3]
dt[country_iso3 %in% names(GADM_ECON_REMAP),
   country_econ := GADM_ECON_REMAP[country_iso3]]
n_remap <- dt[country_econ != country_iso3, .N]
if (n_remap > 0) {
  log_message(sprintf("Remapped %s cells from GADM disputed codes for economic joins: %s",
                      format(n_remap, big.mark = ","),
                      paste(names(GADM_ECON_REMAP), "->", GADM_ECON_REMAP,
                            collapse = ", ")))
}

stopifnot(file.exists(crop_price_pre_file), file.exists(crop_price_glb_file),
          file.exists(trucking_file))
pre <- setDT(read_parquet(crop_price_pre_file))
glb <- setDT(read_parquet(crop_price_glb_file))
tau <- setDT(read_parquet(trucking_file))

dt <- merge(dt, tau[, .(country_econ = country_iso3, tau_usd_ton_km, tau_source)],
            by = "country_econ", all.x = TRUE)
# A missing tau makes EVERY crop in the cell unpriceable, so the cell is dropped
# by the n_crops_priced guard - silently removing whole countries. 0e now builds
# tau over all ISO3 codes so this should be zero; the fallback keeps a stale
# lookup table from deleting a country's worth of observations without saying so.
n_notau <- dt[is.na(tau_usd_ton_km), .N]
log_message(sprintf("Cells without a trucking price: %s (%.2f%%)",
                    format(n_notau, big.mark = ","), 100 * n_notau / nrow(dt)))
if (n_notau > 0) {
  miss_ctry <- sort(unique(dt$country_econ[is.na(dt$tau_usd_ton_km)]))
  log_message(sprintf("  WARNING: no tau for %d countries: %s",
                      length(miss_ctry), paste(miss_ctry, collapse = ", ")))
  log_message("  Falling back to the global median tau. Re-run 0e to fix properly.")
  dt[is.na(tau_usd_ton_km), `:=`(tau_usd_ton_km = median(tau$tau_usd_ton_km, na.rm = TRUE),
                                 tau_source = "global_median_fallback")]
}

# tau_r: USD per TONNE hauled out of this cell.
dt[, tau_r := tau_usd_ton_km * travel_time_cities / 60]

# ==============================================================================
# 4. BUILD THE PRESSURE MEASURES
# ==============================================================================
# Accumulated crop by crop rather than reshaped to long: 2M cells x 44 crops is
# 88M rows, and every operation here is a rowwise product that does not need it.

log_message("Building pressure measures...")

price_w <- dcast(pre, country_iso3 ~ crop, value.var = "price_usd_t")

# Country -> price-row index, resolved once. The earlier draft merged price_w
# onto dt inside the crop loop, which reallocated a 2M-row, 100-column table 44
# times; match() does the same join as a single integer index.
ci        <- match(dt$country_econ, price_w$country_iso3)
glb_price <- setNames(as.list(glb$price_usd_t), glb$crop)

n_nocountry <- sum(is.na(ci))
log_message(sprintf("Cells whose country has NO FAO price row: %s (%.1f%%) - %s",
                    format(n_nocountry, big.mark = ","),
                    100 * n_nocountry / nrow(dt),
                    paste(sort(unique(dt$country_econ[is.na(ci)])), collapse = ", ")))
log_message("  these fall back to the global median crop price")
tau_r <- dt$tau_r
n     <- nrow(dt)

unpriced_crops <- character(0)    # crops with no price in any country
press_cell <- numeric(n)          # sum_c s_cr * (p_c - tau_r) * y_cr
press_max  <- rep(NA_real_, n)    # max_c (p_c - tau_r) * y_cr
n_priced   <- integer(n)

for (i in seq_len(nrow(xw))) {
  crop <- xw$earthstat[i]
  if (!crop %in% names(price_w)) next

  p <- price_w[[crop]][ci]            # USD/t, country-specific

  # Coalesce onto the global median wherever the country reports no price at
  # all. COD, MMR, PNG and GAB are absent from FAOSTAT producer prices entirely
  # - four of the largest tropical-forest countries - and without this they
  # would enter with pressure = 0, which reads as "no agriculture" when it means
  # "no FAO submission". That missingness tracks state capacity, which tracks
  # protection, so the silent zero would bias beta directly.
  #
  # A crop can also have NO price anywhere: plantain is in EarthStat and in the
  # crosswalk but FAOSTAT reports no producer price for it in any country, so
  # glb_price[["plantain"]] is NULL. Assigning NULL into a subset is an error
  # ("replacement has length zero"), which is what killed job 1447204. Such a
  # crop contributes nothing to any cell, so skip it and record it.
  gp <- glb_price[[crop]]
  if (!is.null(gp) && is.finite(gp)) {
    p[is.na(p)] <- gp
  }
  if (all(is.na(p))) {
    unpriced_crops <- c(unpriced_crops, crop)
    next
  }
  y <- dt[[xw$yield_col[i]]]          # t/ha on the crop's harvested area
  s <- dt[[xw$cropshare_col[i]]]      # harvested area / cell area

  pi_c <- (p - tau_r) * y             # USD per ha OF THIS CROP

  ok <- !is.na(pi_c) & !is.na(s)
  press_cell[ok] <- press_cell[ok] + s[ok] * pi_c[ok]
  n_priced[ok]   <- n_priced[ok] + 1L

  # na.rm keeps the running max defined from the first non-NA crop onward.
  press_max <- pmax(press_max, pi_c, na.rm = TRUE)
}

if (length(unpriced_crops) > 0) {
  log_message(sprintf("Crops with NO price in any country, excluded from pressure (%d): %s",
                      length(unpriced_crops), paste(unpriced_crops, collapse = ", ")))
}

dt[, `:=`(pressure_cell  = press_cell,
          pressure_max   = press_max,
          n_crops_priced = n_priced)]
rm(press_cell, press_max, n_priced, ci, tau_r); gc()

# Normalized variant: USD per hectare of CROPLAND rather than of cell. Undefined
# where nothing is harvested, which is most of the undisturbed sample - keeping
# both is the point. This is pressure_cell rescaled, not a separate sum.
dt[, pressure_crop := fifelse(is.na(cropshare_total) | cropshare_total <= 0,
                              NA_real_, pressure_cell / cropshare_total)]

# Report in $100/ha so logit coefficients are readable.
for (v in c("pressure_cell", "pressure_crop", "pressure_max")) {
  dt[, (paste0(v, "_100")) := get(v) / 100]
}

# ==============================================================================
# 5. ESTIMATION SAMPLE
# ==============================================================================

dt[, forested_base := !is.na(undist_base) & undist_base > FOREST_THRESHOLD]

# n_crops_priced == 0 means NOTHING was priced for this cell, so pressure_cell
# is a structural zero rather than a measured one. The global fallback above
# should make this empty; the filter stays as a guard so such cells can never
# enter the regression disguised as genuine zeros.
n_unpriced <- dt[forested_base == TRUE & n_crops_priced == 0, .N]
if (n_unpriced > 0) {
  log_message(sprintf("Dropping %s candidate cells with ZERO priced crops",
                      format(n_unpriced, big.mark = ",")))
}

est <- dt[forested_base == TRUE & n_crops_priced > 0 &
          !is.na(pressure_cell) & is.finite(pressure_cell) &
          !is.na(ever_protected) & !is.na(area_km2) &
          !is.na(centroid_lat) & !is.na(centroid_lon)]

log_message("----------------------------------------")
log_message(sprintf("Estimation sample: %s cells (%.1f%% of %s)",
                    format(nrow(est), big.mark = ","),
                    100 * nrow(est) / nrow(dt),
                    format(nrow(dt), big.mark = ",")))
log_message(sprintf("  ever protected: %s (%.1f%%)",
                    format(sum(est$ever_protected), big.mark = ","),
                    100 * mean(est$ever_protected)))
log_message(sprintf("  countries: %d", uniqueN(est$country_iso3)))

# The degeneracy this specification has to own up to.
z <- est[, mean(pressure_cell == 0)]
log_message(sprintf("  pressure_cell EXACTLY ZERO: %.1f%% of the sample", 100 * z))
log_message(sprintf("  pressure_cell < 0 (transport exceeds price): %.1f%%",
                    100 * est[, mean(pressure_cell < 0)]))
if (z > 0.9) {
  log_message("  WARNING: >90% of the sample is at zero, so beta is identified")
  log_message("  off a thin non-zero tail. Read pressure_crop, and GAEZ, instead.")
}

for (v in c("pressure_cell", "pressure_crop", "pressure_max")) {
  x <- est[[v]]; x <- x[is.finite(x)]
  if (!length(x)) next
  log_message(sprintf("  %-14s n=%s  mean %9.2f  sd %9.2f  p50 %8.2f  p99 %9.2f",
                      v, format(length(x), big.mark = ","),
                      mean(x), sd(x), median(x), quantile(x, .99)))
}

# ==============================================================================
# 6. SPEC LADDER
# ==============================================================================
# travel_time_cities and dist_to_city_km are deliberately EXCLUDED from the
# controls: tau_r is built from travel time, so conditioning on it absorbs the
# transport leg of the very object being measured.

# elevation enters as asinh, not log(x+1): the DEM carries below-sea-level cells
# and log(elevation_m + 1) returns NaN for anything under -1 m, which fixest
# drops as an NA. That silently removed ~660 observations in job 1476554 with
# only a "NaNs produced" warning. asinh is defined on the whole line and behaves
# like log for large positive values, so the transform is unchanged where it
# matters.
cv <- function() vcov_conley(lat = ~centroid_lat, lon = ~centroid_lon,
                             cutoff = CONLEY_KM)

log_message("Estimating...")

m1 <- feglm(ever_protected ~ pressure_cell_100, est,
            family = binomial(), weights = ~area_km2, vcov = cv())

m2 <- feglm(ever_protected ~ pressure_cell_100 | country_iso3, est,
            family = binomial(), weights = ~area_km2, vcov = cv())

m3 <- feglm(ever_protected ~ pressure_cell_100 + slope_degrees +
              terrain_ruggedness + asinh(elevation_m) +
              log(pop_density_1990 + 1) | country_iso3, est,
            family = binomial(), weights = ~area_km2, vcov = cv())

m4 <- feglm(ever_protected ~ pressure_cell_100 + slope_degrees +
              terrain_ruggedness + asinh(elevation_m) +
              log(pop_density_1990 + 1) + log(agc_2010 + 1) +
              log(biodiv_thr + 1) | country_iso3, est,
            family = binomial(), weights = ~area_km2, vcov = cv())

models <- list("Pooled" = m1, "+ country FE" = m2,
               "+ cost side" = m3, "+ benefit side" = m4)

etable(models, digits = 4)

# Robustness: the cropland-normalized measure, where it exists.
est_c <- est[is.finite(pressure_crop)]
log_message(sprintf("pressure_crop defined on %s cells (%.1f%% of sample)",
                    format(nrow(est_c), big.mark = ","),
                    100 * nrow(est_c) / max(nrow(est), 1)))
if (nrow(est_c) > 1000) {
  m5 <- feglm(ever_protected ~ pressure_crop_100 + slope_degrees +
                terrain_ruggedness + asinh(elevation_m) +
                log(pop_density_1990 + 1) | country_iso3, est_c,
              family = binomial(), weights = ~area_km2, vcov = cv())
  etable(list("cropland-normalized" = m5), digits = 4)
}

# ==============================================================================
# 7. OUTPUT
# ==============================================================================

tex_path <- file.path(output_dir, "protection_logit.tex")
etable(models, file = tex_path, replace = TRUE, digits = 4,
       title = paste("Protection selection on \\$-denominated agricultural",
                     "pressure (logit, area-weighted, Conley 50km SEs)"),
       label = "tab:protection_logit")
log_message(sprintf("Wrote %s", tex_path))

sum_path <- file.path(output_dir, "pressure_summary.txt")
sink(sum_path)
cat("PRESSURE MEASURE DIAGNOSTICS\n============================\n")
cat(sprintf("Generated: %s\n\n", Sys.time()))
cat(sprintf("Estimation sample: %s cells, %d countries, %.1f%% ever protected\n",
            format(nrow(est), big.mark = ","), uniqueN(est$country_iso3),
            100 * mean(est$ever_protected)))
cat(sprintf("Crops priced: %d\n", nrow(xw)))
cat(sprintf("pressure_cell exactly zero: %.1f%%\n", 100 * z))
cat(sprintf("pressure_cell negative:     %.1f%%\n\n",
            100 * est[, mean(pressure_cell < 0)]))
cat("Yields are EarthStat ACTUAL (zero where the crop is not grown), so this\n")
cat("measures demonstrated agricultural value, not agronomic potential.\n\n")
print(est[, .(mean_pressure_cell = mean(pressure_cell),
              mean_pressure_max  = mean(pressure_max, na.rm = TRUE),
              mean_tau_r         = mean(tau_r, na.rm = TRUE),
              n = .N),
          by = .(ever_protected)])
sink()
log_message(sprintf("Wrote %s", sum_path))

log_message(sprintf("Completed in %.1f minutes",
                    as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
