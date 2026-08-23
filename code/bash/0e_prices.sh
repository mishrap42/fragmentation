#!/bin/bash -l
#SBATCH --time=00:30:00
#SBATCH --job-name=TMF_PRICES
#SBATCH --account=mishralab
#SBATCH --partition=expansion
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --output=LOGS/%x.%A.out
#SBATCH --error=LOGS/%x.%A.err

# ==============================================================================
# STAGE 0e: Crop Prices and Transport Costs
# Reads FAOSTAT producer prices + World Bank trucking costs from the shared
# Global Forest Repo tree and writes three small lookup tables to Data/lookup/:
#   fao_crop_prices.parquet       country x crop x year
#   crop_price_preperiod.parquet  country x crop, 1991-2000 mean
#   crop_price_global.parquet     crop, global median (fallback)
#   trucking_cost.parquet         country, USD per ton-km
#
# No rasters, no grid, no dependency on any other stage - it reads two CSVs and
# a crosswalk. Sized at 30 min / 16G purely for headroom; it runs in seconds and
# is fine on a login node. It is a batch script only so the run is logged like
# everything else.
#
# PREREQUISITE: none
# ==============================================================================

# sbatch runs a COPY of this script from /var/spool/slurmd/<job>/, so
# ${BASH_SOURCE[0]} does NOT point at code/bash/ under SLURM. Resolve from
# SLURM_SUBMIT_DIR (submit from the project root), and fall back to the script
# location for plain `bash <script>` runs. The rerun_* helpers submit a sed'd
# copy from $(mktemp), which is exactly why the fallback cannot be trusted here.
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  HERE="$SLURM_SUBMIT_DIR/code/bash"
else
  HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
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
    PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
fi

echo "Project root: $PROJECT_ROOT"
cd "$PROJECT_ROOT" || exit 1

R_SCRIPT="code/build/0e_prepare_crop_prices.R"

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
