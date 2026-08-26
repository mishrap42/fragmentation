#!/bin/bash -l
#SBATCH --time=02:00:00
#SBATCH --job-name=TMF_GAEZ
#SBATCH --account=mishralab
#SBATCH --partition=expansion
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --array=1-344
#SBATCH --mem=8G
#SBATCH --output=LOGS/%x.%A_%a.out
#SBATCH --error=LOGS/%x.%A_%a.err

# ==============================================================================
# STAGE 0d: GAEZ v4 Attainable-Yield Extraction
# Per sub-tile, extracts 150 rasters onto the grid:
#   75 x <crop>200a_yld  -> gaez_<crop>200a   agro-climatic,   kg/ha DM
#   75 x <crop>200b_yld  -> gaez_<crop>200b   agro-ecological, kg/ha DM
#
# 150 rasters against 0c's 346, on the same grid and the same 5-arcmin WGS84
# footprint, so 2h is generous against 0c's measured cost.
#
# Array size: 344 sub-tiles (same as Stage 0). 40 of them have no grid geometry
# and will exit non-zero - those sub-tiles are empty ocean, not failures, so do
# NOT chain anything off this with --dependency=afterok.
#
# PREREQUISITE: Run Stage 0 first to create grid files
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
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
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

R_SCRIPT="code/build/0d_extract_gaez.R"

if [[ -f "$R_SCRIPT" ]]; then
    echo "Running R script: $R_SCRIPT with task ID: $SLURM_ARRAY_TASK_ID"
    Rscript "$R_SCRIPT" "$SLURM_ARRAY_TASK_ID"
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
