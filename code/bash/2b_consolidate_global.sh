#!/bin/bash -l
#SBATCH --time=04:00:00
#SBATCH --job-name=TMF_CONS_GLOBAL
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --output=LOGS/%x.%A.out
#SBATCH --error=LOGS/%x.%A.err

# ==============================================================================
# STAGE 2B: Global Consolidation
# Merges all tile-level consolidated files into single global dataset
# ==============================================================================

# Create LOGS directory if it doesn't exist
mkdir -p LOGS

# Print job information
echo "========================================"
echo "SLURM Job Information"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Submit Directory: $SLURM_SUBMIT_DIR"
echo "Node: $SLURMD_NODENAME"
echo "Start Time: $(date)"
echo "========================================"

# Determine project root
if [[ -n "$SLURM_SUBMIT_DIR" ]]; then
    PROJECT_ROOT="$SLURM_SUBMIT_DIR"
else
    PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
fi

echo "Project root: $PROJECT_ROOT"
cd "$PROJECT_ROOT" || exit 1

R_SCRIPT="Code/2b_consolidate_global.R"

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
