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
#
# SAMPLE: forested-at-baseline cells NOT already protected by 2000. Designations
# on or before 2000 predate the 1991-2000 price window and are dropped, so the
# regressor always precedes the decision it is meant to explain.
# Outputs: output/analysis/
#            protection_logit.tex        coefficients + dollar equivalents
#            binscatter_pressure.png     binsreg, dots only
#
# Usage: Rscript code/analysis/6_protection_logit.R
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(duckdb)
  library(DBI)
  library(arrow)
  library(fixest)
  library(ggplot2)
  library(scales)
  library(binsreg)
  library(sf)
  library(rnaturalearth)
})

here::i_am('code/analysis/6_protection_logit.R')
source("code/build/BUILD_workspace.R")

output_dir <- file.path(project_root, "output", "analysis", "protection_logit")
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
    MAX(CASE WHEN year = %d THEN Undisturbed_TMF END)    AS undist_base,
    MAX(CASE WHEN year = 1990 THEN forest_cover END)      AS fc_1990,
    MAX(CASE WHEN year = 2000 THEN forest_cover END)      AS fc_2000
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

# GAEZ potential-yield columns for the same crops. gaez is pipe-separated where
# GAEZ ships no aggregate layer (rice = ricd|ricw); 1b has already materialised
# the composite as gaez_<earthstat-name>200a/b, so the column to ask for is
# always keyed on the EarthStat crop name, never on the raw GAEZ codes.
xw[, gaez_a := fifelse(nzchar(gaez), paste0("gaez_", earthstat, "200a"), NA_character_)]
xw[, gaez_b := fifelse(nzchar(gaez), paste0("gaez_", earthstat, "200b"), NA_character_)]
# ... except for single-code crops, where 1b builds no composite and the column
# carries the GAEZ code itself.
xw[nzchar(gaez) & !grepl("|", gaez, fixed = TRUE),
   `:=`(gaez_a = paste0("gaez_", gaez, "200a"),
        gaez_b = paste0("gaez_", gaez, "200b"))]

HAS_GAEZ <- any(grepl("^gaez_", avail))
if (HAS_GAEZ) {
  xw[, `:=`(gaez_a = fifelse(gaez_a %in% avail, gaez_a, NA_character_),
            gaez_b = fifelse(gaez_b %in% avail, gaez_b, NA_character_))]
  log_message(sprintf("GAEZ layers usable: %d crops (200a), %d crops (200b)",
                      sum(!is.na(xw$gaez_a)), sum(!is.na(xw$gaez_b))))
} else {
  log_message("No gaez_* columns in the cross-section; IV section will be skipped.")
  xw[, `:=`(gaez_a = NA_character_, gaez_b = NA_character_)]
}

need <- intersect(c("grid_id", "cropland_frac", "pasture_frac",
                    "cropshare_total", xw$cropshare_col, xw$yield_col,
                    na.omit(xw$gaez_a), na.omit(xw$gaez_b)), avail)

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

# ------------------------------------------------------------------------------
# POTENTIAL profit, the instrument. Same prices, same transport wedge, GAEZ
# potential yield in place of EarthStat realised yield.
# ------------------------------------------------------------------------------
#   pi_pot_cr = (p_c - tau_r) * y^GAEZ_cr / 1000       GAEZ is kg/ha, prices $/t
#
# TWO weightings, because the choice is not innocuous:
#
#   pressure_pot_cell  sum_c cropshare_cr * pi_pot_cr
#       The literal "only the yield leg is instrumented" object. But it inherits
#       the EarthStat weights, so it is ZERO wherever nothing is grown - the same
#       42.9% of the sample where the endogenous regressor is zero. It therefore
#       has no first-stage power in exactly the forest interior the design is
#       about. Built for comparison, not as the headline.
#
#   pressure_pot_max   max_c pi_pot_cr
#       Best-use potential profit. No planted-area weights at all, so it is
#       defined on every cell GAEZ covers, cleared or not. This is the instrument
#       with power where it is needed.
#
# 200a (agro-climatic, climate only) is the default: 200b folds in soil and
# terrain, which are on the RHS as controls, so instrumenting with it leans on
# variation the second stage is already conditioning away.

if (HAS_GAEZ) {
  gz <- xw[!is.na(gaez_a)]
  log_message(sprintf("Building potential-profit instruments from %d GAEZ crops...",
                      nrow(gz)))

  ci2       <- match(dt$country_econ, price_w$country_iso3)
  tau_v     <- dt$tau_r
  pot_cell  <- numeric(nrow(dt))
  pot_max   <- rep(NA_real_, nrow(dt))
  pot_cellb <- numeric(nrow(dt))
  n_pot     <- integer(nrow(dt))

  for (i in seq_len(nrow(gz))) {
    crop <- gz$earthstat[i]
    if (!crop %in% names(price_w)) next
    p <- price_w[[crop]][ci2]
    gp <- glb_price[[crop]]
    if (!is.null(gp) && is.finite(gp)) p[is.na(p)] <- gp
    if (all(is.na(p))) next

    s  <- dt[[gz$cropshare_col[i]]]
    ya <- dt[[gz$gaez_a[i]]] / 1000            # kg/ha -> t/ha
    pia <- (p - tau_v) * ya

    ok <- !is.na(pia) & !is.na(s)
    pot_cell[ok] <- pot_cell[ok] + s[ok] * pia[ok]
    n_pot[ok]    <- n_pot[ok] + 1L
    pot_max      <- pmax(pot_max, pia, na.rm = TRUE)

    if (!is.na(gz$gaez_b[i])) {
      yb  <- dt[[gz$gaez_b[i]]] / 1000
      pib <- (p - tau_v) * yb
      okb <- !is.na(pib) & !is.na(s)
      pot_cellb[okb] <- pot_cellb[okb] + s[okb] * pib[okb]
    }
  }

  # NOTE: these potential-profit measures are currently UNUSED. They were built
  # to instrument realised profitability; that approach has been removed. They
  # are kept because they are cheap to construct and are the natural candidate
  # for entering the ladder DIRECTLY as an exogenous measure of land value -
  # unlike the EarthStat measure they are defined on every cell and do not
  # collapse when forest state is conditioned on. Delete this block if that is
  # not wanted.
  dt[, `:=`(pressure_pot_cell  = pot_cell,
            pressure_pot_max   = pot_max,
            pressure_pot_cellb = pot_cellb,
            n_crops_pot        = n_pot)]
  rm(pot_cell, pot_max, pot_cellb, n_pot, ci2, tau_v); gc()

  for (v in c("pressure_pot_cell", "pressure_pot_max", "pressure_pot_cellb")) {
    dt[, (paste0(v, "_100")) := get(v) / 100]
  }
}
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

# ------------------------------------------------------------------------------
# Pre-period forest state, both in PERCENTAGE POINTS so the dollar equivalents
# read per pp rather than per unit-share.
#   forest_1990_pp   forest cover at the start of the window
#   defor_9000_pp    forest cover LOST over 1990-2000; positive = loss
# Both are measured entirely before DESIG_CUTOFF, so they precede every
# designation in the estimation sample and are not post-treatment. They are
# still endogenous in the weaker sense that prior clearing and future protection
# share unobserved drivers - they absorb variation, they do not identify it.
# ------------------------------------------------------------------------------
dt[, forest_1990_pp := fc_1990 * 100]
dt[, defor_9000_pp  := (fc_1990 - fc_2000) * 100]

# n_crops_priced == 0 means NOTHING was priced for this cell, so pressure_cell
# is a structural zero rather than a measured one. The global fallback above
# should make this empty; the filter stays as a guard so such cells can never
# enter the regression disguised as genuine zeros.
n_unpriced <- dt[forested_base == TRUE & n_crops_priced == 0, .N]
if (n_unpriced > 0) {
  log_message(sprintf("Dropping %s candidate cells with ZERO priced crops",
                      format(n_unpriced, big.mark = ",")))
}

# ------------------------------------------------------------------------------
# RISK SET: only designations that START after the price window closes.
# ------------------------------------------------------------------------------
# Prices are the 1991-2000 mean (PRICE_PRE_YEARS). A cell already protected by
# 2000 was not at risk over the window those prices describe - its designation
# was decided under earlier prices, and it cannot respond to the regressor. It
# is neither a valid treated unit nor a valid control, so it is DROPPED rather
# than recoded to zero, which would put already-protected land in the control
# group and bias beta toward zero.
#
# Cells protected with an UNKNOWN desig_year are also dropped: they cannot be
# placed on either side of the cutoff, and assuming they are post-2000 would
# manufacture treatment.
#
# Surviving definition:
#   ever_protected = 1  designated after DESIG_CUTOFF
#   ever_protected = 0  never protected
DESIG_CUTOFF <- max(PRICE_PRE_YEARS)

dt[, prot_pre  := ever_protected == 1 & !is.na(desig_year) & desig_year <= DESIG_CUTOFF]
dt[, prot_post := ever_protected == 1 & !is.na(desig_year) & desig_year >  DESIG_CUTOFF]
dt[, prot_unk  := ever_protected == 1 & is.na(desig_year)]

cand <- dt[forested_base == TRUE & n_crops_priced > 0 &
           !is.na(pressure_cell) & is.finite(pressure_cell) &
           !is.na(ever_protected) & !is.na(area_km2) &
           !is.na(centroid_lat) & !is.na(centroid_lon) &
           !is.na(forest_1990_pp) & !is.na(defor_9000_pp)]

log_message("----------------------------------------")
log_message(sprintf("Candidate cells before the %d cutoff: %s",
                    DESIG_CUTOFF, format(nrow(cand), big.mark = ",")))
log_message(sprintf("  protected on/before %d (dropped, not at risk): %s",
                    DESIG_CUTOFF, format(cand[prot_pre == TRUE, .N], big.mark = ",")))
log_message(sprintf("  protected, desig_year unknown (dropped):       %s",
                    format(cand[prot_unk == TRUE, .N], big.mark = ",")))
log_message(sprintf("  protected after %d (treated):                %s",
                    DESIG_CUTOFF, format(cand[prot_post == TRUE, .N], big.mark = ",")))
log_message(sprintf("  never protected (control):                    %s",
                    format(cand[ever_protected == 0, .N], big.mark = ",")))

est <- cand[prot_pre == FALSE & prot_unk == FALSE]
est[, ever_protected := as.integer(prot_post)]

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
#
# CONTROLS ARE LINEAR, deliberately. The point of this table is the marginal
# rate of substitution against the priced regressor: with x_k entering linearly,
# -gamma_k / beta_p is read straight off as "dollars per hectare of forgone
# agricultural profit that the planner behaves as if one extra unit of x_k is
# worth". Under log(x_k) that ratio would be dollars per LOG POINT, which has no
# natural unit and cannot be quoted per tCO2 or per metre.
#
# elevation enters raw for the same reason. The earlier asinh() guarded against
# NaN from below-sea-level cells under log(x+1); linear needs no guard.

# ------------------------------------------------------------------------------
# BIOME, for the country x biome fixed-effect variant.
#
# Read from the lookup 1_stacked_did.R already built and cached; nothing new is
# computed here. If the cache is absent that variant is skipped rather than
# built, since constructing it is a separate decision.
#
# ADM-2 is NOT available: the grid was cut on GADM ADM_0 only, so the panel
# carries no sub-national identifier. That variant needs a centroid join to
# GADM ADM_2 (47,217 polygons) which has not been made.
# ------------------------------------------------------------------------------

biome_cache <- file.path(project_root, "Data/build/final/grid_biome.parquet")
if (file.exists(biome_cache)) {
  est <- merge(est, as.data.table(arrow::read_parquet(biome_cache)),
               by = "grid_id", all.x = TRUE)
  log_message(sprintf("Biome attached from cache: %d biomes, %s cells unassigned",
                      uniqueN(est$biome, na.rm = TRUE),
                      format(sum(is.na(est$biome)), big.mark = ",")))
} else {
  log_message("Biome cache absent; the country x biome variant will be skipped.")
}

cv <- function() vcov_conley(lat = ~centroid_lat, lon = ~centroid_lon,
                             cutoff = CONLEY_KM)

log_message("Estimating...")

# Columns (4) and (5) are ALTERNATIVES, not a nesting. Forest state and carbon
# are close to collinear - above-ground carbon is essentially forest extent
# priced in tCO2 - so entering both put two measures of the same thing on the
# right-hand side, and forest_1990_pp went from 0.0277*** to 0.0036 n.s. between
# them. Splitting the blocks lets each be read on its own:
#   (4) forest block   forest_1990_pp + defor_9000_pp
#   (5) carbon block   agc_2010 + biodiv_thr
# Both sit on top of (3), so they are directly comparable to each other and to
# the geography-only column.
BLOCK_GEO    <- c("slope_degrees", "elevation_m", "pop_density_1990")
BLOCK_FOREST <- c("forest_1990_pp", "defor_9000_pp")
BLOCK_CARBON <- c("agc_2010", "biodiv_thr")

mk_f <- function(rhs, fe) {
  rhs_txt <- paste(c("pressure_cell_100", rhs), collapse = " + ")
  as.formula(if (is.null(fe)) paste("ever_protected ~", rhs_txt)
             else paste("ever_protected ~", rhs_txt, "|", fe))
}

fit <- function(f) feglm(f, est, family = binomial(), weights = ~area_km2,
                         vcov = cv())

# Fixed-effect variants. Each writes its own table; nothing is chosen here.
#   country        country_iso3          the original
#   (ADM-2 omitted: no sub-national identifier exists in the panel)
#   country_biome  country_iso3 ^ biome  matches the stacked DiD stratum
FE_SPECS <- list(
  list(key = "countryFE",      fe = "country_iso3",           label = "Country"),
  list(key = "country_biomeFE", fe = "country_iso3^biome",    label = "Country $\\times$ biome")
)

# Defined at top level, not inside run_fe_spec(): the IV section at the end of
# the script also needs it, and a local copy left it invisible there ("object
# 'LABEL' not found" killed job 1621938 after every model had already been fit).
LABEL <- c(
  pressure_cell_100    = "Agricultural profitability (\\$100/ha)",
  fit_pressure_cell_100 = "Agricultural profitability (\\$100/ha)",
  pressure_pot_max_100 = "Potential profitability, best use (\\$100/ha)",
  slope_degrees        = "Slope (degrees)",
  elevation_m          = "Elevation (m)",
  pop_density_1990     = "Population density, 1990 (per km$^2$)",
  forest_1990_pp       = "Forest cover, 1990 (pp)",
  defor_9000_pp        = "Forest lost 1990--2000 (pp)",
  agc_2010             = "Above-ground carbon (tCO$_2$/ha)",
  biodiv_thr           = "Threatened species richness"
)

run_fe_spec <- function(FE_KEY, FE_FE, FE_LABEL) {
  log_message(sprintf("--- fixed effects: %s ---", FE_LABEL))
  m1 <- fit(mk_f(character(0), NULL))
  m2 <- fit(mk_f(character(0), FE_FE))
  m3 <- fit(mk_f(BLOCK_GEO, FE_FE))
  m4 <- fit(mk_f(c(BLOCK_GEO, BLOCK_FOREST), FE_FE))
  m5 <- fit(mk_f(c(BLOCK_GEO, BLOCK_CARBON), FE_FE))
  mods <- list(m1, m2, m3, m4, m5)

  log_message("Ladder, coefficient on agricultural profitability ($100/ha):")
  for (i in seq_along(mods)) {
    bb <- coef(mods[[i]])[["pressure_cell_100"]]
    ss <- sqrt(vcov(mods[[i]])["pressure_cell_100", "pressure_cell_100"])
    log_message(sprintf("  (%d) %8.4f  (SE %.4f, t %6.2f)", i, bb, ss, bb / ss))
  }

  # ==============================================================================
  # 7. DOLLAR-EQUIVALENT VALUATIONS (marginal rate of substitution)
  # ==============================================================================
  # In  Pr(protect) = Lambda(alpha + beta_p * pressure + sum_k gamma_k x_k),
  # holding the index fixed:
  #
  #     beta_p * d(pressure) + gamma_k * dx_k = 0
  #     =>  d(pressure)/dx_k = -gamma_k / beta_p
  #
  # pressure is scaled in $100/ha, so multiplying by 100 gives USD/ha per unit of
  # x_k: the agricultural profit the siting process behaved as if it was willing
  # to forgo for one more unit of x_k. Revealed preference, not social value, and
  # only as credible as beta_p.
  #
  # Delta method on r = -gamma/beta:
  #   dr/dgamma = -1/beta ,  dr/dbeta = gamma/beta^2
  # The ratio blows up as beta_p -> 0, so beta_p's t is reported in the table.

  dollar_equiv <- function(model, price_var = "pressure_cell_100", scale = 100) {
    b  <- coef(model); V <- vcov(model)
    bp <- b[[price_var]]; vp <- V[price_var, price_var]

    rbindlist(lapply(setdiff(names(b), price_var), function(k) {
      g <- b[[k]]; vg <- V[k, k]; cg <- V[k, price_var]
      dg <- -1 / bp * scale
      db <-  g / bp^2 * scale
      data.table(term = k, usd_per_unit = -g / bp * scale,
                 se = sqrt(max(dg^2 * vg + db^2 * vp + 2 * dg * db * cg, 0)))
    }))
  }

  de <- dollar_equiv(m5)

  # agc_2010 is above-ground biomass CARBON (tC/ha). tCO2 = tC * 44/12, so both
  # the coefficient and its dollar value are rescaled by 12/44 to read per tCO2 -
  # same convention as 5_protection_lpm.R's committed-carbon calculation.
  C_TO_CO2 <- 44 / 12
  de[term == "agc_2010", `:=`(usd_per_unit = usd_per_unit / C_TO_CO2,
                              se           = se / C_TO_CO2)]

  # ------------------------------------------------------------------------------
  # Logit elasticity of the protection probability w.r.t. agricultural profit.
  #   p = Lambda(eta),  d p / d x = beta * p (1-p)
  #   elasticity = (d p / d x) * (x / p) = beta * x * (1 - p)
  # Evaluated at the area-weighted sample means, matching the estimation weights.
  # ------------------------------------------------------------------------------
  xbar_100 <- est[, weighted.mean(pressure_cell_100, area_km2)]
  pbar     <- est[, weighted.mean(ever_protected,    area_km2)]
  beta_p   <- coef(m5)[["pressure_cell_100"]]
  se_bp    <- sqrt(vcov(m5)["pressure_cell_100", "pressure_cell_100"])
  elast    <- beta_p * xbar_100 * (1 - pbar)
  elast_se <- abs(se_bp * xbar_100 * (1 - pbar))

  log_message("----------------------------------------")
  log_message(sprintf("beta_p = %.4f (t = %.2f)", beta_p, beta_p / se_bp))
  log_message(sprintf("elasticity of P(protect) wrt ag profit = %.4f (SE %.4f)",
                      elast, elast_se))
  for (i in seq_len(nrow(de))) {
    log_message(sprintf("  %-18s %10.2f  (SE %8.2f)  USD/ha",
                        de$term[i], de$usd_per_unit[i], de$se[i]))
  }

  # ==============================================================================
  # 8. TABLE
  # ==============================================================================
  # ONE table. Column (1) is the fitted logit in raw units; column (2) is the same
  # row re-expressed as dollars per hectare of forgone agricultural profit. The
  # price variable is the numeraire, so it has no entry in (2).


  stars <- function(t) if (!is.finite(t)) "" else
                       if (abs(t) > 2.576) "^{***}" else
                       if (abs(t) > 1.960) "^{**}"  else
                       if (abs(t) > 1.645) "^{*}"   else ""

  # agc_2010 is reported per tCO2 in EVERY column, coefficient and dollar alike,
  # so the units never mix within a row. tCO2 = tC * 44/12.
  rescale <- function(k, v) if (k == "agc_2010") v / C_TO_CO2 else v

  cell <- function(m, k) {
    b <- coef(m)
    if (!k %in% names(b)) return(list("", ""))
    se <- sqrt(vcov(m)[k, k])
    list(sprintf("$%.4f%s$", rescale(k, b[[k]]), stars(b[[k]] / se)),
         sprintf("$(%.4f)$", rescale(k, se)))
  }

  NCOL <- length(mods) + 1L   # coefficient columns + the dollar column

  tex <- file.path(output_dir, sprintf("protection_logit_%s.tex", FE_KEY))
  f <- file(tex, "w")
  cat("\\begin{table}[htbp]\n\\centering\n",
      "\\caption{\\label{tab:protection_logit_", FE_KEY, "} Protected-area siting and agricultural profitability.",
      " Area-weighted logits of post-2000 designation on cells forested at baseline and not already",
      " protected by ", DESIG_CUTOFF, "; Conley (50 km) standard errors in parentheses.",
      " Column (", NCOL, ") re-expresses the column (", NCOL - 1L, ") coefficients as",
      " $-\\gamma_k/\\beta_p$: the agricultural profit per hectare the siting process behaves as if",
      " it will forgo for one more unit of the attribute, with delta-method standard errors.}\n",
      "\\begin{tabular}{l", paste(rep("c", NCOL), collapse = ""), "}\n\\toprule\n",
      sep = "", file = f)
  cat(" & ", paste(sprintf("(%d)", seq_len(NCOL)), collapse = " & "), " \\\\\n",
      sprintf(" & \\multicolumn{%d}{c}{Coefficient} & \\$/ha equivalent \\\\\n",
              length(mods)),
      # coefficient columns occupy tabular columns 2..(1+k); the dollar column is
      # tabular column 2+k, one further right than NCOL suggests because column 1
      # is the row label.
      sprintf("\\cmidrule(lr){2-%d}\\cmidrule(lr){%d-%d}\n",
              length(mods) + 1L, length(mods) + 2L, length(mods) + 2L),
      sep = "", file = f)

  for (k in names(LABEL)) {
    top <- character(0); bot <- character(0)
    for (m in mods) { cc <- cell(m, k); top <- c(top, cc[[1]]); bot <- c(bot, cc[[2]]) }
    if (all(top == "")) next
    dcol <- if (k == "pressure_cell_100") "---" else {
      r <- de[term == k]
      if (nrow(r) == 0) "" else sprintf("$%.2f$", r$usd_per_unit)
    }
    dse <- if (k == "pressure_cell_100" || nrow(de[term == k]) == 0) "" else
             sprintf("$(%.2f)$", de[term == k, se])
    cat(LABEL[[k]], " & ", paste(top, collapse = " & "), " & ", dcol, " \\\\\n",
        " & ", paste(bot, collapse = " & "), " & ", dse, " \\\\\n", sep = "", file = f)
  }

  blank <- paste(rep("", length(mods)), collapse = " & ")
  cat("\\midrule\n",
      FE_LABEL, " fixed effects & ", paste(c("No", rep("Yes", length(mods) - 1L)),
                                        collapse = " & "), " & \\\\\n",
      "Observations & ", paste(sapply(mods, function(m) format(m$nobs, big.mark = ",")),
                               collapse = " & "), " & \\\\\n",
      "Pseudo $R^2$ & ", paste(sapply(mods, function(m) sprintf("%.4f", fitstat(m, "pr2")$pr2)),
                               collapse = " & "), " & \\\\\n",
      "\\midrule\n",
      "Elasticity of $\\Pr(\\text{protect})$ w.r.t. profit & ",
      paste(rep("", length(mods) - 1L), collapse = " & "),
      sprintf(" & $%.4f$ & \\\\\n", elast),
      " & ", paste(rep("", length(mods) - 1L), collapse = " & "),
      sprintf(" & $(%.4f)$ & \\\\\n", elast_se),
      "\\bottomrule\n\\end{tabular}\n\\end{table}\n", sep = "", file = f)
  close(f)
  log_message(sprintf("Wrote %s", tex))

  invisible(list(m5 = m5, mods = mods))
}

results <- list()
for (s in FE_SPECS) {
  results[[s$key]] <- tryCatch(run_fe_spec(s$key, s$fe, s$label),
    error = function(e) { log_message(sprintf("  FAILED (%s): %s", s$key, conditionMessage(e))); NULL })
}

# Downstream figures use the country-FE fit, the original specification.
stopifnot(!is.null(results$countryFE))
m5   <- results$countryFE$m5
mods <- results$countryFE$mods

# ==============================================================================
# 9. BINSCATTER
# ==============================================================================
# binsreg, dots only. It handles the 42.9% mass at zero itself: the repeated
# knots collapse into a single dot and it warns which bins it dropped, which is
# the honest treatment - equal-count quantile bins would otherwise scatter one
# number across a third of the plot.

X_LABEL <- "Agricultural profitability, $/ton"

bsr <- binsreg(y = est$ever_protected, x = est$pressure_cell,
               weights = est$area_km2, nbins = 20)
dots <- as.data.table(bsr$data.plot[[1]]$data.dots)
log_message(sprintf("binsreg returned %d dots (20 requested; repeated knots collapse the zero mass)",
                    nrow(dots)))

p <- ggplot(dots, aes(x = x, y = fit)) +
  geom_point(colour = "steelblue", size = 2.4) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = X_LABEL, y = "Probability of protection") +
  theme_classic(base_size = 12) +
  theme(panel.grid = element_blank())

ggsave(file.path(output_dir, "binscatter_pressure.png"), p,
       width = 6.5, height = 4.5, dpi = 300)
log_message(sprintf("Wrote %s",
                    file.path(output_dir, "binscatter_pressure.png")))

# ==============================================================================
# 10. HETEROGENEITY IN THE PROFITABILITY SLOPE
# ==============================================================================
# Two cuts of the column (5) specification, both letting beta_p vary:
#   by COUNTRY  - one slope per country, country FE still absorbing levels
#   by COHORT   - one slope per designation-year block, estimated stacked
#
# Everything else in the specification is held at its column (5) form so the
# heterogeneous slopes are comparable to the pooled number.

MIN_TREATED_COHORT  <- 500    # cohorts are stacked; below this the stack is pointless

RHS_CTRL <- c("slope_degrees", "elevation_m", "pop_density_1990",
              "forest_1990_pp", "defor_9000_pp", "agc_2010", "biodiv_thr")

# ------------------------------------------------------------------------------
# 10a. COUNTRY-LEVEL SLOPES
# ------------------------------------------------------------------------------
# i(country_iso3, pressure_cell_100, ref = NULL) gives one slope per country;
# the country FE still absorbs the level, so these are within-country slopes and
# the pooled beta_p is their (precision-weighted) centre of mass.
#
# NO minimum-sample screen. Countries are excluded only where a slope does not
# EXIST, never because it would be imprecise:
#   n_treated == 0 or == n   the outcome is constant within the country, so the
#                            country FE absorbs it entirely and fixest drops it
#   sd(pressure) == 0        no regressor variation, so no slope to estimate
# Everything else is kept however small. Precision is a separate problem from
# identification, and a small country's noisy estimate is information to be
# shrunk later, not information to discard now - 5_protection_lpm.R already
# carries a normal-normal shrinkage step that does exactly this for country WTP.
# Earlier thresholds (n >= 2000, treated >= 100) admitted 18 of 97 countries.

ctry_n <- est[, .(n = .N, n_treated = sum(ever_protected),
                  sd_press = sd(pressure_cell_100)), by = country_iso3]

ctry_n[, identified := n_treated > 0 & n_treated < n &
                       !is.na(sd_press) & sd_press > 0]
keep_ctry <- ctry_n[identified == TRUE, country_iso3]

log_message(sprintf("Country slopes: %d of %d countries identified (no size screen)",
                    length(keep_ctry), nrow(ctry_n)))
log_message(sprintf("  excluded: %d never protected, %d always protected, %d no pressure variation",
                    ctry_n[n_treated == 0, .N],
                    ctry_n[n_treated > 0 & n_treated >= n, .N],
                    ctry_n[n_treated > 0 & n_treated < n &
                           (is.na(sd_press) | sd_press <= 0), .N]))
log_message(sprintf("  smallest retained country: %d cells, %d treated",
                    ctry_n[identified == TRUE, min(n)],
                    ctry_n[identified == TRUE, min(n_treated)]))

est_c <- est[country_iso3 %chin% keep_ctry]

# i(f, x) with a CONTINUOUS x returns a slope for every level - no reference
# category is needed or wanted, because these are slopes rather than dummies and
# the country FE already absorbs the levels. Do NOT add ref = NULL: fixest
# 0.13.2 throws "node stack overflow" on it, which reads like a data problem but
# is purely the argument.
f_ctry <- as.formula(paste0(
  "ever_protected ~ i(country_iso3, pressure_cell_100) + ",
  paste(RHS_CTRL, collapse = " + "), " | country_iso3"))

# glm.iter raised from the default 25: with ~60 country-specific slopes the IRLS
# hit the iteration cap in job 1621938 ("Absence of convergence") and the VCOV
# came back non-positive-semi-definite and silently "fixed", which makes every
# reported SE untrustworthy. Convergence is checked explicitly below rather than
# left to a warning that scrolls past.
m_ctry <- feglm(f_ctry, est_c, family = binomial(), weights = ~area_km2,
                vcov = cv(), glm.iter = 100)

if (!isTRUE(m_ctry$convStatus)) {
  log_message("  WARNING: country-slope logit did NOT converge - slopes and SEs are unreliable")
} else {
  log_message(sprintf("  converged in %s IRLS iterations",
                      if (is.null(m_ctry$iterations)) "?" else m_ctry$iterations))
}

cb <- coef(m_ctry); cs <- sqrt(diag(vcov(m_ctry)))
sel <- grep("^country_iso3::", names(cb))
# Names come back as "country_iso3::BRA:pressure_cell_100". Stripped with two
# plain substitutions rather than a backreference - shorter, and it cannot be
# silently mangled by escaping.
ctry_lab <- sub("^country_iso3::", "", names(cb)[sel])
ctry_lab <- sub(":pressure_cell_100$", "", ctry_lab)

beta_ctry <- data.table(country_iso3 = ctry_lab,
                        beta = as.numeric(cb[sel]), se = as.numeric(cs[sel]))
stopifnot(!any(grepl("[:]", beta_ctry$country_iso3)))
beta_ctry <- merge(beta_ctry, ctry_n, by = "country_iso3")
beta_ctry[, t := beta / se]

log_message(sprintf("  slopes estimated for %d countries; %d negative, %d significant at 5%%",
                    nrow(beta_ctry), beta_ctry[beta < 0, .N],
                    beta_ctry[abs(t) > 1.96, .N]))
log_message(sprintf("  precision: %d with |t| < 1 (essentially uninformative), median |t| = %.2f",
                    beta_ctry[abs(t) < 1, .N], median(abs(beta_ctry$t), na.rm = TRUE)))
log_message(sprintf("  treated-cell count across mapped countries: min %d, median %d",
                    min(beta_ctry$n_treated), as.integer(median(beta_ctry$n_treated))))

world <- ne_countries(scale = "medium", returnclass = "sf")

# ------------------------------------------------------------------------------
# Winsorize the country slopes before mapping.
# ------------------------------------------------------------------------------
# With no size screen, a country identified off a handful of treated cells can
# return an enormous slope. Unwinsorized, one such estimate sets the colour
# scale and every other country collapses toward the midpoint - the map goes
# white and reads as "no heterogeneity" when the opposite is true.
#
# The estimates THEMSELVES are winsorized, not merely the scale limits, so the
# legend numbers correspond to what is actually drawn. Symmetric quantiles of
# the signed distribution, so the zero midpoint stays meaningful and positive
# and negative slopes are treated alike. What got clipped is logged.
WINSOR_P <- 0.05

lo_w <- as.numeric(quantile(beta_ctry$beta, WINSOR_P,     na.rm = TRUE))
hi_w <- as.numeric(quantile(beta_ctry$beta, 1 - WINSOR_P, na.rm = TRUE))
lim  <- max(abs(c(lo_w, hi_w)))          # symmetric about zero

beta_ctry[, beta_w := pmax(pmin(beta, lim), -lim)]
n_clip <- beta_ctry[abs(beta) > lim, .N]
log_message(sprintf("  winsorized at +/-%.3f (%.0f%% tails): %d of %d countries clipped",
                    lim, 100 * WINSOR_P, n_clip, nrow(beta_ctry)))
if (n_clip > 0) {
  log_message(sprintf("    clipped: %s",
                      paste(sprintf("%s (%.2f)",
                                    beta_ctry[abs(beta) > lim, country_iso3],
                                    beta_ctry[abs(beta) > lim, beta]),
                            collapse = ", ")))
}
log_message(sprintf("  raw slope range: [%.2f, %.2f]; mapped range: [%.2f, %.2f]",
                    min(beta_ctry$beta), max(beta_ctry$beta), -lim, lim))

wmap <- merge(world, beta_ctry, by.x = "iso_a3", by.y = "country_iso3",
              all.x = TRUE)

p_map <- ggplot(wmap) +
  geom_sf(aes(fill = beta_w), colour = "gray70", linewidth = 0.1) +
  scale_fill_gradient2(low = "#08306b", mid = "white", high = "#67000d",
                       midpoint = 0, limits = c(-lim, lim), na.value = "gray92",
                       name = expression(beta[p])) +
  coord_sf(ylim = c(-35, 35), expand = FALSE) +
  theme_void() +
  theme(legend.position = "right")

# 360 degrees of longitude against 70 of latitude is a ~5:1 frame; anything
# taller letterboxes the map inside a half-empty canvas.
ggsave(file.path(output_dir, "beta_by_country_map.png"), p_map,
       width = 11, height = 3.2, dpi = 300, bg = "white")
log_message(sprintf("Wrote %s", file.path(output_dir, "beta_by_country_map.png")))

# ------------------------------------------------------------------------------
# 10b. COHORT-LEVEL SLOPES
# ------------------------------------------------------------------------------
# Cohort is defined only for TREATED cells - a never-protected cell has no
# designation year - so this cannot be a single interacted regression. Each
# cohort is estimated on its own stack: {cells designated in that cohort} plus
# the FULL never-protected control pool.
#
# The control pool is therefore REUSED across cohorts. That is standard for a
# stacked design and leaves each beta unbiased, but the cohort estimates are not
# independent of one another, so the plotted intervals should not be read as if
# a difference between two cohorts were a two-sample test.

est[, cohort := cut(desig_year,
                    breaks = c(DESIG_CUTOFF, 2005, 2010, 2015, Inf),
                    labels = c("2001-2005", "2006-2010", "2011-2015", "2016+"),
                    right = TRUE)]

controls_pool <- est[ever_protected == 0]
coh_levels <- levels(est$cohort)
coh_n <- est[ever_protected == 1, .N, by = cohort][order(cohort)]
log_message("Treated cells by cohort:")
for (i in seq_len(nrow(coh_n))) {
  log_message(sprintf("  %-10s %s", as.character(coh_n$cohort[i]),
                      format(coh_n$N[i], big.mark = ",")))
}

f_coh <- as.formula(paste0("ever_protected ~ pressure_cell_100 + ",
                           paste(RHS_CTRL, collapse = " + "), " | country_iso3"))

beta_coh <- rbindlist(lapply(coh_levels, function(cl) {
  n_tr <- est[ever_protected == 1 & cohort == cl, .N]
  if (is.na(n_tr) || n_tr < MIN_TREATED_COHORT) {
    log_message(sprintf("  skipping cohort %s: only %s treated cells",
                        cl, format(n_tr, big.mark = ",")))
    return(NULL)
  }
  stack <- rbind(est[ever_protected == 1 & cohort == cl], controls_pool)
  mm <- feglm(f_coh, stack, family = binomial(), weights = ~area_km2,
              vcov = cv())
  data.table(cohort = cl,
             beta = coef(mm)[["pressure_cell_100"]],
             se = sqrt(vcov(mm)["pressure_cell_100", "pressure_cell_100"]),
             n_treated = n_tr)
}))

if (nrow(beta_coh) > 0) {
  beta_coh[, `:=`(lo = beta - 1.96 * se, hi = beta + 1.96 * se)]
  beta_coh[, cohort := factor(cohort, levels = coh_levels)]

  log_message("Cohort slopes on agricultural profitability:")
  for (i in seq_len(nrow(beta_coh))) {
    log_message(sprintf("  %-10s %8.4f  [%7.4f, %7.4f]  n_treated %s",
                        as.character(beta_coh$cohort[i]), beta_coh$beta[i],
                        beta_coh$lo[i], beta_coh$hi[i],
                        format(beta_coh$n_treated[i], big.mark = ",")))
  }

  # Intervals are the point of a coefficient plot, so they stay - unlike the
  # binscatter, where they were decoration.
  p_coh <- ggplot(beta_coh, aes(x = cohort, y = beta)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "gray50") +
    geom_hline(yintercept = coef(m5)[["pressure_cell_100"]],
               linetype = "dashed", linewidth = 0.3, colour = "gray50") +
    geom_linerange(aes(ymin = lo, ymax = hi), colour = "steelblue",
                   linewidth = 0.6) +
    geom_point(colour = "steelblue", size = 2.4) +
    labs(x = "Designation cohort",
         y = "Coefficient on agricultural profitability") +
    theme_classic(base_size = 12) +
    theme(panel.grid = element_blank())

  ggsave(file.path(output_dir, "beta_by_cohort.png"), p_coh,
         width = 6.5, height = 4.5, dpi = 300, bg = "white")
  log_message(sprintf("Wrote %s", file.path(output_dir, "beta_by_cohort.png")))
}

log_message(sprintf("Completed in %.1f minutes",
                    as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
