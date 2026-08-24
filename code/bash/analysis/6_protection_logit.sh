#!/bin/bash -l
#SBATCH --time=06:00:00
#SBATCH --job-name=PROT_LOGIT
#SBATCH --account=mishralab
#SBATCH --partition=expansion
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=180G
#SBATCH --output=LOGS/%x.%A.out
#SBATCH --error=LOGS/%x.%A.err

# ==============================================================================
# 6_protection_logit.sh - PA selection on $-denominated agricultural pressure
# Single job (no array). Joins the cell-year panel (collapsed via DuckDB), the
# EarthStat cropland cross-section, FAO producer prices and country trucking
# costs; builds
#     pressure_r = sum_c cropshare_cr * (p_c - tau_r) * y_cr   [USD/ha of cell]
# and fits a spec ladder of area-weighted logits with Conley 50km SEs. Writes:
#   - protection_logit.tex   (spec ladder)
#   - pressure_summary.txt   (measure diagnostics)
#
# Memory: 180G rather than the LPM's 64G. This job holds the ~2M-cell collapse
# AND the 44-crop slice of the cropland cross-section (88 float columns) at the
# same time, which the LPM does not. Raised from 128G, and walltime 3h -> 6h,
# when the heterogeneity section was added: it fits ten models with Conley SEs
# rather than five, and the country-interacted model carries ~100 extra
# coefficients.
#
# Prerequisites, in order:
#   Stage 0c  -> Data/build/yields/*.parquet
#   Stage 1b  -> Data/build/final/TMF_5km_cropland.parquet
#   Stage 0e  -> Data/lookup/{crop_price_preperiod,crop_price_global,trucking_cost}.parquet
#   Stage 6   -> Data/build/final/TMF_5km_panel.parquet
# ==============================================================================

# sbatch runs a COPY of this script from /var/spool/slurmd/<job>/, so
# ${BASH_SOURCE[0]} does NOT point at code/bash/ under SLURM. Resolve from
# SLURM_SUBMIT_DIR (submit from the project root), and fall back to the script
# location for plain `bash <script>` runs. The rerun_* helpers submit a sed'd
# copy from $(mktemp), which is exactly why the fallback cannot be trusted here.
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  HERE="$SLURM_SUBMIT_DIR/code/bash"
else
  HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fi
[ -f "$HERE/config.sh" ] || { echo "config.sh not found at $HERE" >&2; exit 1; }
source "$HERE/config.sh"
frag_load_modules || exit 1
frag_ensure_logs

echo "========================================"
echo "SLURM Job Information"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Submit Directory: $SLURM_SUBMIT_DIR"
echo "Node: $SLURMD_NODENAME"
echo "Start Time: $(date)"
echo "========================================"

if [[ -n "$SLURM_SUBMIT_DIR" ]]; then
    PROJECT_ROOT="$SLURM_SUBMIT_DIR"
else
    PROJECT_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
fi

echo "Project root: $PROJECT_ROOT"
cd "$PROJECT_ROOT" || exit 1

R_SCRIPT="code/analysis/6_protection_logit.R"

if [[ -f "$R_SCRIPT" ]]; then
    echo "Running R script: $R_SCRIPT"
    Rscript "$R_SCRIPT"
    EXIT_CODE=$?
    echo "R script exited with code: $EXIT_CODE"
else
    echo "Error: R script not found at $R_SCRIPT"
    exit 1
fi

echo "========================================"
echo "End Time: $(date)"
echo "========================================"

exit $EXIT_CODE
