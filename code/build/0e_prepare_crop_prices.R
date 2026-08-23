# ==============================================================================
# STAGE 0e: PREPARE CROP PRICES AND TRANSPORT COSTS
# ==============================================================================
# Turns two shared non-spatial files into country-keyed tables the analysis can
# join onto the cropland cross-section. No rasters, no grid - runs in seconds on
# a login node.
#
# Input:  <GF>/non-spatial/prices/Prices_E_All_Data_NOFLAG.csv     (FAOSTAT)
#         <GF>/non-spatial/trade costs/SED Trucking Cost/*.csv     (World Bank)
#         Data/lookup/crop_crosswalk_earthstat_fao.csv             (this repo)
# Output: Data/lookup/fao_crop_prices.parquet       country x crop x year
#         Data/lookup/crop_price_preperiod.parquet  country x crop, pre-period
#         Data/lookup/trucking_cost.parquet         country
#
# WHY A PRE-PERIOD MEAN
# ---------------------
# The static protection logit needs a price that is PREDETERMINED relative to
# the designation decision. Most WDPA designations in this panel post-date 2000,
# so prices are averaged over PRICE_PRE_YEARS (1991-2000) and held fixed. Using
# a contemporaneous price would let the outcome move the regressor.
#
# UNITS AND DEFLATION
# -------------------
# FAOSTAT producer prices are NOMINAL USD/tonne. They are NOT deflated here.
# Averaging over a fixed decade makes the result "average 1990s USD/tonne",
# which is internally consistent across crops and countries - all of them share
# the same window - and is the right object for a single cross-section. Any
# time-varying specification (a designation hazard) MUST deflate first; the
# year-level table is written out precisely so that stays possible.
#
# Usage: Rscript code/build/0e_prepare_crop_prices.R
# ==============================================================================

here::i_am('code/build/0e_prepare_crop_prices.R')
source("code/build/BUILD_workspace.R")

suppressPackageStartupMessages(library(countrycode))

start_time <- Sys.time()
log_job_start("0e_prepare_crop_prices.R", task_id = 1)

PRICE_PRE_YEARS <- 1991:2000

stopifnot(file.exists(fao_prices_path),
          file.exists(trucking_cost_path),
          file.exists(crop_crosswalk_file))

# ==============================================================================
# CROP CROSSWALK
# ==============================================================================

xw <- fread(crop_crosswalk_file)
log_message(sprintf("Crosswalk: %d crops (%d flagged key)",
                    nrow(xw), sum(xw$key_crop == 1)))

# ==============================================================================
# FAO PRODUCER PRICES
# ==============================================================================

log_message("Reading FAOSTAT producer prices...")

price <- fread(fao_prices_path, encoding = "Latin-1")

price <- price[Element == "Producer Price (USD/tonne)"]
log_message(sprintf("  %s USD/tonne rows, %d areas, %d items",
                    format(nrow(price), big.mark = ","),
                    uniqueN(price$Area), uniqueN(price$Item)))

year_cols <- grep("^Y[0-9]{4}$", names(price), value = TRUE)
price <- melt(price, id.vars = c("Area", "Item"), measure.vars = year_cols,
              variable.name = "year", value.name = "price_usd_t",
              variable.factor = FALSE)
price[, year := as.integer(sub("^Y", "", year))]
price <- price[!is.na(price_usd_t) & price_usd_t > 0]

# Join on the LEGACY item name - this file predates the 2021 CPC rename. All 44
# crosswalk entries resolve; a drop to zero here means the price file was
# swapped for a current-vintage one and fao_item_current is the join key.
matched <- xw[fao_item_legacy %in% unique(price$Item)]
log_message(sprintf("  Crosswalk items found in price file: %d/%d",
                    nrow(matched), nrow(xw)))
if (nrow(matched) == 0) {
  stop("No crosswalk items matched. Check whether the price file vintage uses ",
       "fao_item_current instead of fao_item_legacy.")
}

price <- merge(price, xw[, .(fao_item_legacy, crop = earthstat)],
               by.x = "Item", by.y = "fao_item_legacy")

# FAOSTAT "Area" includes aggregates (World, Africa, ...); countrycode returns
# NA for those, which is exactly the filter we want.
price[, country_iso3 := suppressWarnings(
  countrycode(Area, "country.name.en", "iso3c"))]
n_unmapped <- uniqueN(price[is.na(country_iso3)]$Area)
if (n_unmapped > 0) {
  log_message(sprintf("  Dropping %d unmapped areas (FAO aggregates): %s",
                      n_unmapped,
                      paste(head(unique(price[is.na(country_iso3)]$Area), 8),
                            collapse = ", ")))
}
price <- price[!is.na(country_iso3)]

prices_year <- price[, .(price_usd_t = mean(price_usd_t, na.rm = TRUE)),
                     by = .(country_iso3, crop, year)]
setorder(prices_year, country_iso3, crop, year)

log_message(sprintf("Year-level prices: %s rows, %d countries, %d crops, %d-%d",
                    format(nrow(prices_year), big.mark = ","),
                    uniqueN(prices_year$country_iso3),
                    uniqueN(prices_year$crop),
                    min(prices_year$year), max(prices_year$year)))

# ==============================================================================
# PRE-PERIOD PRICE, WITH A GLOBAL FALLBACK
# ==============================================================================
# Country x crop coverage is thin: a tropical-forest country typically reports
# prices for a handful of its crops. Without a fallback the profit index would
# silently drop those crops from that country's sum, making the index a function
# of FAO REPORTING COVERAGE as much as of agronomy - and reporting correlates
# with state capacity, which correlates with protection. The fallback is the
# crop's global median over the same window, and price_source records which cells
# got it so the sensitivity can be tested.

pre <- prices_year[year %in% PRICE_PRE_YEARS,
                   .(price_usd_t = mean(price_usd_t, na.rm = TRUE),
                     n_years = .N),
                   by = .(country_iso3, crop)]

global_pre <- prices_year[year %in% PRICE_PRE_YEARS,
                          .(price_global = median(price_usd_t, na.rm = TRUE)),
                          by = crop]
global_out <- copy(global_pre)

countries <- sort(unique(prices_year$country_iso3))
grid <- CJ(country_iso3 = countries, crop = sort(unique(xw$earthstat)))

pre <- merge(grid, pre, by = c("country_iso3", "crop"), all.x = TRUE)
pre <- merge(pre, global_pre, by = "crop", all.x = TRUE)

pre[, price_source := fifelse(!is.na(price_usd_t), "observed",
                       fifelse(!is.na(price_global), "global_median", "missing"))]
pre[is.na(price_usd_t), price_usd_t := price_global]
pre[is.na(n_years), n_years := 0L]
pre[, price_global := NULL]

# Written as its own table because the country x crop grid above only covers
# countries FAOSTAT actually reports. Four major tropical-forest countries -
# COD, MMR, PNG, GAB - report NO producer prices at all and so appear nowhere in
# `pre`. A consumer joining only on country would hand those cells a silent
# zero, which is not "no agriculture" but "no FAO submission", and FAO reporting
# tracks state capacity, which tracks protection. Consumers must coalesce onto
# this table for any country missing from `pre`.
setnames(global_out, "price_global", "price_usd_t")
log_message(sprintf("Global fallback prices: %d crops", nrow(global_out)))

log_message(sprintf("Pre-period prices (%d-%d): %s country-crop cells",
                    min(PRICE_PRE_YEARS), max(PRICE_PRE_YEARS),
                    format(nrow(pre), big.mark = ",")))
for (src in c("observed", "global_median", "missing")) {
  n <- pre[price_source == src, .N]
  log_message(sprintf("    %-14s %6s (%.1f%%)", src, format(n, big.mark = ","),
                      100 * n / nrow(pre)))
}

# ==============================================================================
# TRUCKING COST
# ==============================================================================

log_message("Reading World Bank trucking costs...")

truck <- fread(trucking_cost_path)
setnames(truck, "Average/Median trucking price (USD per ton-km)", "tau_usd_ton_km")

truck[, country_iso3 := suppressWarnings(
  countrycode(Country, "country.name.en", "iso3c"))]
truck <- truck[!is.na(country_iso3) & !is.na(tau_usd_ton_km),
               .(tau_usd_ton_km = median(tau_usd_ton_km, na.rm = TRUE)),
               by = country_iso3]

log_message(sprintf("  %d countries with an observed trucking price",
                    nrow(truck)))

# Imputation. The companion (4_build.R:1141-1151, 1522-1528) does neighbour-mean
# then continent-median; there is no contiguity table in this repo, so this does
# continent-median then global-median and records which. Coarser, and stated as
# such rather than presented as the same procedure.
tau <- merge(CJ(country_iso3 = countries), truck, by = "country_iso3", all.x = TRUE)
tau[, continent := suppressWarnings(countrycode(country_iso3, "iso3c", "continent"))]

tau[, tau_continent := median(tau_usd_ton_km, na.rm = TRUE), by = continent]
tau_global <- median(tau$tau_usd_ton_km, na.rm = TRUE)

tau[, tau_source := fifelse(!is.na(tau_usd_ton_km), "observed",
                     fifelse(!is.na(tau_continent), "continent_median", "global_median"))]
tau[is.na(tau_usd_ton_km), tau_usd_ton_km := tau_continent]
tau[is.na(tau_usd_ton_km), tau_usd_ton_km := tau_global]
tau[, tau_continent := NULL]

for (src in c("observed", "continent_median", "global_median")) {
  n <- tau[tau_source == src, .N]
  log_message(sprintf("    %-18s %4d countries", src, n))
}
log_message(sprintf("  tau range: %.4f - %.4f USD/ton-km",
                    min(tau$tau_usd_ton_km), max(tau$tau_usd_ton_km)))

# ==============================================================================
# WRITE
# ==============================================================================

write_atomic(prices_year, crop_prices_file)
write_atomic(pre,         crop_price_pre_file)
write_atomic(global_out,  crop_price_glb_file)
write_atomic(tau,         trucking_file)

log_message(sprintf("Wrote:\n  %s\n  %s\n  %s\n  %s",
                    crop_prices_file, crop_price_pre_file,
                    crop_price_glb_file, trucking_file))

log_job_end(start_time)
