#!/bin/bash -l
#SBATCH --time=04:00:00
#SBATCH --job-name=TMF_WDPA
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --array=1-87
#SBATCH --mem=64G
#SBATCH --output=LOGS/%x.%A_%a.out
#SBATCH --error=LOGS/%x.%A_%a.err

# ==============================================================================
# STAGE 5: WDPA Protected Area Extraction
# Tasks 1-86: Extract WDPA for each tile
# Task 87: Consolidate all tile outputs
# ==============================================================================

mkdir -p LOGS

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
    PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
fi

echo "Project root: $PROJECT_ROOT"
cd "$PROJECT_ROOT" || exit 1

R_SCRIPT="Code/5_extract_WDPA.R"

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
