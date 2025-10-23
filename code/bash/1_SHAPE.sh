#!/bin/bash -l
#SBATCH --time=02:00:00
#SBATCH --job-name=SHAPE
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --array=1
#SBATCH --output=LOGS/%x.%A_%a.out
#SBATCH --error=LOGS/%x.%A_%a.err

# Create LOGS directory if it doesn't exist
mkdir -p LOGS

# Print job information
echo "========================================"
echo "SLURM Job Information"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Submit Directory: $SLURM_SUBMIT_DIR"
echo "Start Time: $(date)"
echo "========================================"

# SLURM_SUBMIT_DIR is the directory where sbatch was called from
if [[ -n "$SLURM_SUBMIT_DIR" ]]; then
    echo "SLURM submit directory: $SLURM_SUBMIT_DIR"
    # Assume script is submitted from project root or adjust path accordingly
    R_SCRIPT="${SLURM_SUBMIT_DIR}/code/build/1_calculate_shapes.R"
else
    # Fallback if not running under SLURM
    R_SCRIPT="/Volumes/f007z52/Fragmentation/code/build/1_calculate_shapes.R"
fi

echo "Array job task ID: $SLURM_ARRAY_TASK_ID"
echo "Looking for R script at: $R_SCRIPT"

if [[ -f "$R_SCRIPT" ]]; then
    echo "Running R script: $R_SCRIPT with task ID: $SLURM_ARRAY_TASK_ID"
    Rscript "$R_SCRIPT" "$SLURM_ARRAY_TASK_ID"
else
    echo "Error: R script not found at $R_SCRIPT"
    echo "Current working directory: $(pwd)"
    echo "Available environment variables:"
    env | grep -E "(SLURM|HOME)" | sort
    exit 1
fi