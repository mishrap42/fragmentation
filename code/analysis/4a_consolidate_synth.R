# ==============================================================================
# 4a_consolidate_synth.R
# Consolidate per-country synthetic control outputs from 4_synthetic_control.R
# into the project-level datasets for downstream processing.
#
#   synth_project_effects.csv - project ID, country, treatment effect, project
#                               area (+ inference and fit diagnostics)
#   synth_dynamic_effects.csv - project ID, country, year, event_time, and the
#                               per-period effect tau_t (= real_y - synth_y)
#   synth_stacked_panel.csv   - stacked treated+control panel (all cohorts) for
#                               a feols event study
#
# Usage: Rscript 4a_consolidate_synth.R
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

if (Sys.getenv("SLURM_SUBMIT_DIR") != "") {
  project_root <- Sys.getenv("SLURM_SUBMIT_DIR")
} else {
  project_root <- here::here()
}

synth_dir   <- file.path(project_root, "Data/build/intermediate/synth")
effects_dir <- file.path(synth_dir, "effects")
dynamic_dir <- file.path(synth_dir, "dynamic")
stacked_dir <- file.path(synth_dir, "stacked")

eff_files <- list.files(effects_dir, pattern = "^effects_.*\\.csv$",
                        full.names = TRUE)
dyn_files <- list.files(dynamic_dir, pattern = "^dynamic_.*\\.csv$",
                        full.names = TRUE)
stk_files <- list.files(stacked_dir, pattern = "^stacked_.*\\.csv$",
                        full.names = TRUE)

cat("Found", length(eff_files), "effects files in", effects_dir, "\n")
cat("Found", length(dyn_files), "dynamic files in", dynamic_dir, "\n")
cat("Found", length(stk_files), "stacked files in", stacked_dir, "\n")

if (length(eff_files) == 0) {
  stop("No per-country effects files found. Run 4_synthetic_control.R first.")
}

# --- Dataset 1: project-level treatment effects -------------------------------
effects <- rbindlist(lapply(eff_files, fread), use.names = TRUE, fill = TRUE)
setorder(effects, country_iso3, wdpa_pid)

eff_out <- file.path(synth_dir, "synth_project_effects.csv")
fwrite(effects, eff_out)
cat("Wrote", nrow(effects), "project effects ->", eff_out, "\n")

# --- Dataset 2: dynamic (per-period) effects ----------------------------------
dynamic <- if (length(dyn_files))
  rbindlist(lapply(dyn_files, fread), use.names = TRUE, fill = TRUE) else
  data.table()
if (nrow(dynamic)) setorder(dynamic, country_iso3, wdpa_pid, year)

dyn_out <- file.path(synth_dir, "synth_dynamic_effects.csv")
fwrite(dynamic, dyn_out)
cat("Wrote", nrow(dynamic), "dynamic-effect rows ->", dyn_out, "\n")

# --- Dataset 3: stacked treated+control panel ---------------------------------
stacked <- if (length(stk_files))
  rbindlist(lapply(stk_files, fread), use.names = TRUE, fill = TRUE) else
  data.table()
if (nrow(stacked)) setorder(stacked, country_iso3, wdpa_pid, unit_id, year)

stk_out <- file.path(synth_dir, "synth_stacked_panel.csv")
fwrite(stacked, stk_out)
cat("Wrote", nrow(stacked), "stacked-panel rows ->", stk_out, "\n")

# --- Quick summary ------------------------------------------------------------
cat("\n=== Synthetic control summary ===\n")
cat("Projects:", nrow(effects), "across", uniqueN(effects$country_iso3),
    "countries\n")
if (nrow(effects) > 0) {
  cat(sprintf("Mean ATT: %.4f | Median ATT: %.4f\n",
              mean(effects$treatment_effect, na.rm = TRUE),
              median(effects$treatment_effect, na.rm = TRUE)))
  cat("Projects with Fisher p < 0.10:",
      sum(effects$fishers_exact_pvalue < 0.10, na.rm = TRUE), "\n")
}
