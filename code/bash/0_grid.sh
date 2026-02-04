#!/bin/bash -l
#SBATCH --time=01:00:00
#SBATCH --job-name=TMF_GRID
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --array=1-86
#SBATCH --mem=32G
#SBATCH --output=LOGS/%x.%A_%a.out
#SBATCH --error=LOGS/%x.%A_%a.err

# ==============================================================================
# STAGE 0: Grid Generation
# Creates 1km x 1km equal-area grid cells for each TMF tile
# ==============================================================================

# Create LOGS directory if it doesn't exist
mkdir -p LOGS

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
    # Fallback for local testing
    PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
fi

echo "Project root: $PROJECT_ROOT"

# Change to project directory
cd "$PROJECT_ROOT" || exit 1

# Set R script path
R_SCRIPT="Code/0_create_grid.R"

if [[ -f "$R_SCRIPT" ]]; then
    echo "Running R script: $R_SCRIPT with task ID: $SLURM_ARRAY_TASK_ID"
    Rscript "$R_SCRIPT" "$SLURM_ARRAY_TASK_ID"
    EXIT_CODE=$?
    echo "R script exited with code: $EXIT_CODE"
else
    echo "Error: R script not found at $R_SCRIPT"
    echo "Current working directory: $(pwd)"
    ls -la Code/
    exit 1
fi

echo "========================================"
echo "End Time: $(date)"
echo "========================================"

exit $EXIT_CODE
