#!/bin/bash -l
#SBATCH --time=03:00:00
#SBATCH --job-name=VALIDATE_YG
#SBATCH --account=mishralab
#SBATCH --partition=expansion
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --output=LOGS/%x.%j.out
#SBATCH --error=LOGS/%x.%j.err

# ==============================================================================
# VALIDATE YIELDS (Stage 0c) + GFED (Stage 1a)
#
# Value-level validation, not just file counts: diagnose_pipeline.R's
# file.exists() check passed while every yield file was a 964-byte empty.
# Runs independently of stages 2-6; needs only 0c, 1a and the stage 0 grids.
#
# --mem=64G: yields load fully (~0.12 GB on disk); GFED is ~100M rows across
# 304 files and is aggregated per file rather than rbound, so peak memory is
# one sub-tile.
# ==============================================================================

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

R_SCRIPT="code/build/validate_yields_gfed.R"

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
