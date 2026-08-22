#!/bin/bash -l
#SBATCH --time=04:00:00
#SBATCH --job-name=TMF_COVARIATES
#SBATCH --account=mishralab
#SBATCH --partition=expansion
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --array=1-86
#SBATCH --mem=120G
#SBATCH --output=LOGS/%x.%A_%a.out
#SBATCH --error=LOGS/%x.%A_%a.err

# ==============================================================================
# ORDERING: THIS STAGE RUNS AFTER STAGE 4, NOT AFTER STAGE 0.
#
# 5_extract_covariates.R computes signed distance to the forest frontier and
# stopifnot()s on the stage 4 output:
#     Error: file.exists(frontier_file) is not TRUE
# It gets all the way through cities/population/biomass extraction first - about
# an hour per tile - before dying at that last step, so submitting it early
# looks like it is working and then wastes the whole run. Observed 2026-08-19:
# 85 tasks burned this way. Correct position: 0 -> 1 -> 2a -> 3 -> 4 -> 5 -> 6.
# ==============================================================================

# ==============================================================================
# STAGE 5: Static Covariate Extraction
# Extracts static covariates (population, biomass, elevation, etc.) for each tile
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

# Print job information
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

# Determine project root
if [[ -n "$SLURM_SUBMIT_DIR" ]]; then
    PROJECT_ROOT="$SLURM_SUBMIT_DIR"
else
    PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
fi

echo "Project root: $PROJECT_ROOT"

# Change to project directory
cd "$PROJECT_ROOT" || exit 1

# Set R script path
R_SCRIPT="code/build/5_extract_covariates.R"

if [[ -f "$R_SCRIPT" ]]; then
    echo "Running R script: $R_SCRIPT with task ID: $SLURM_ARRAY_TASK_ID"
    Rscript "$R_SCRIPT" "$SLURM_ARRAY_TASK_ID"
    EXIT_CODE=$?
    echo "R script exited with code: $EXIT_CODE"
else
    echo "Error: R script not found at $R_SCRIPT"
    echo "Current working directory: $(pwd)"
    ls -la code/build/
    exit 1
fi

echo "========================================"
echo "End Time: $(date)"
echo "========================================"

exit $EXIT_CODE
