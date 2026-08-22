# ==============================================================================
# install_packages.R — top up the shared group R library for this project.
#
# The tree at /resnick/groups/MishraLab/Rlib/4.5 is built by the Global Forest
# Repo and already carries sf, terra, arrow, exactextractr, data.table, fixest,
# Hmisc, here and ncdf4. These are the packages the Fragmentation pipeline adds.
#
# Run ONCE on a login node, from the project root:
#
#     source code/bash/config.sh && frag_load_modules && \
#       Rscript code/build/install_packages.R
#
# PREREQUISITE: clear the stale 00LOCK-* directories left by the interrupted
# 2026-07-30/31 installs, or install.packages() will refuse to proceed:
#
#     rm -rf /resnick/groups/MishraLab/Rlib/4.5/00LOCK-*
# ==============================================================================

lib <- Sys.getenv("R_LIBS_USER", unset = "/resnick/groups/MishraLab/Rlib/4.5")
if (!dir.exists(lib)) dir.create(lib, recursive = TRUE)
cat("Installing into:", lib, "\n")

# duckdb is the critical one: every analysis script reads the ~67M-row panel
# through it, because arrow trips the thrift size limit on that parquet.
pkgs <- c(
  "duckdb",      # panel reads in code/analysis/*.R
  "tidysynth",   # 4_synthetic_control.R
  "MatchIt",     # 3_matching.R, 3a_matching.R
  "cobalt",      # balance diagnostics alongside MatchIt
  "RANN",        # nearest-neighbour lookups in the matching + frontier code
  "janitor",
  "foreign",
  "ggnewscale"
)

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) == 0L) {
  cat("Nothing to do — all packages already present.\n")
} else {
  cat("Missing:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, lib = lib, repos = "https://cloud.r-project.org")
}

# Verify, and fail loudly so a broken build does not surface later inside a
# 24-hour array job.
still <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(still) > 0L) {
  stop("Failed to install: ", paste(still, collapse = ", "))
}
cat("All packages available.\n")
print(sessionInfo())
