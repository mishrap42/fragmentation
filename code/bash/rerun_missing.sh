#!/bin/bash -l
# ==============================================================================
# Rerun Missing Pipeline Jobs
# ==============================================================================
# Usage: bash code/bash/rerun_missing.sh <stage> [--delete-only]
#   stage: 0, 0c, 1, 2a, 3, 4, 5, or 6
#   --delete-only: Only delete files, don't submit jobs (for manual rerun)
#
# This script DELETES existing output files for the listed jobs before
# resubmitting, ensuring a clean rerun. Use this when you need to refresh
# data that was created with bad parameters or geometry errors.
#
# Pipeline Stages:
#   0  - Grid Creation (N_SUB_TILES jobs)
#   0c - Yield Extraction (N_SUB_TILES jobs, optional)
#   1  - TMF Extraction (N_EXTRACTION_JOBS)
#   2a - Tile Consolidation (N_TMF_TILES jobs)
#   3  - Interior Classification (N_TMF_TILES jobs)
#   4  - Frontier Calculation (N_TMF_TILES jobs)
#   5  - Covariate Extraction (N_TMF_TILES jobs)
#   6  - Final Assembly (1 job)
#
# First run: Rscript code/build/diagnose_pipeline.R
# Or create LOGS/rerun/stage*_missing.txt manually with job IDs to rerun
# ==============================================================================

set -e

# Resource flags and the LOGS directory come from config.sh. A parent's
# account/partition/qos are NOT inherited by the jobs it submits, so every
# sbatch below must forward $SBATCH_COMMON explicitly.
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
[ -f "$HERE/config.sh" ] || { echo "config.sh not found at $HERE" >&2; exit 1; }
source "$HERE/config.sh"

# Determine project root and cd to it
if [[ -n "$SLURM_SUBMIT_DIR" ]]; then
    PROJECT_ROOT="$SLURM_SUBMIT_DIR"
else
    PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
fi
cd "$PROJECT_ROOT" || exit 1
echo "Project root: $PROJECT_ROOT"

# Configuration - must match BUILD_workspace.R
BUILD_DATA_PATH="Data/build"
GRID_PATH="${BUILD_DATA_PATH}/grids"
TMF_PATH="${BUILD_DATA_PATH}/tmf"
CONSOLIDATED_PATH="${BUILD_DATA_PATH}/tmf_consolidated"
INTERIOR_PATH="${BUILD_DATA_PATH}/classification"
FRONTIER_PATH="${BUILD_DATA_PATH}/classification"
COVARIATES_PATH="${BUILD_DATA_PATH}/covariates"

# TMF years configuration (must match BUILD_workspace.R)
TMF_YEAR_START=1990
TMF_YEAR_END=2023
N_TMF_YEARS=$((TMF_YEAR_END - TMF_YEAR_START + 1))  # 34 years

# Function to convert Stage 1 job_id to (tile_id, year)
job_to_year_tile() {
    local job_id=$1
    local year_idx=$(( (job_id - 1) % N_TMF_YEARS ))
    local tile_idx=$(( (job_id - 1) / N_TMF_YEARS + 1 ))
    local year=$((TMF_YEAR_START + year_idx))
    echo "$tile_idx $year"
}

if [[ -z "$1" ]]; then
    echo "Usage: bash code/bash/rerun_missing.sh <stage> [--delete-only]"
    echo "  stage: 0, 0c, 1, 2a, 3, 4, 5, or 6"
    echo "  --delete-only: Only delete files, don't submit jobs"
    echo ""
    echo "First run diagnosis to find missing jobs:"
    echo "  Rscript code/build/diagnose_pipeline.R"
    echo ""
    echo "Or create LOGS/rerun/stage*_missing.txt manually with job IDs to refresh"
    exit 1
fi

STAGE="$1"
DELETE_ONLY=false
if [[ "$2" == "--delete-only" ]]; then
    DELETE_ONLY=true
fi

# Find job file
JOB_FILE="LOGS/rerun/stage${STAGE}_missing.txt"

if [[ ! -f "$JOB_FILE" ]]; then
    echo "No job file found: $JOB_FILE"
    echo ""
    echo "Run diagnosis first:"
    echo "  sbatch code/bash/diagnose_pipeline.sh"
    exit 1
fi

echo "Found: $JOB_FILE"
N_JOBS=$(wc -l < "$JOB_FILE" | tr -d ' ')

if [[ "$N_JOBS" -eq 0 ]]; then
    echo "No jobs to rerun for stage $STAGE"
    exit 0
fi

echo ""
echo "========================================"
echo "Rerunning Stage $STAGE: $N_JOBS jobs"
echo "========================================"

# Read job IDs into array
mapfile -t JOB_IDS < "$JOB_FILE"

# Convert to comma-separated list for SLURM array
JOB_LIST=$(IFS=,; echo "${JOB_IDS[*]}")

# ==============================================================================
# DELETE EXISTING OUTPUT FILES
# ==============================================================================
# This ensures jobs will run fresh, even if output files already exist
# ==============================================================================

echo ""
echo "Deleting existing output files for ${#JOB_IDS[@]} jobs..."

deleted_count=0

for job_id in "${JOB_IDS[@]}"; do
    case "$STAGE" in
        0)
            # Stage 0: Grid files (parquet + gpkg)
            grid_pq=$(printf "%s/grid_sub_%04d.parquet" "$GRID_PATH" "$job_id")
            grid_gpkg=$(printf "%s/grid_sub_%04d.gpkg" "$GRID_PATH" "$job_id")
            for f in "$grid_pq" "$grid_gpkg"; do
                if [[ -f "$f" ]]; then
                    rm -f "$f"
                    deleted_count=$((deleted_count + 1))
                    echo "  Deleted: $f"
                fi
            done
            ;;
        0c)
            # Stage 0c: Yield files
            yield_file=$(printf "%s/yields_sub_%04d.parquet" "$GRID_PATH" "$job_id")
            if [[ -f "$yield_file" ]]; then
                rm -f "$yield_file"
                deleted_count=$((deleted_count + 1))
                echo "  Deleted: $yield_file"
            fi
            ;;
        1)
            # Stage 1: TMF extraction files (need to convert job_id to tile_id + year)
            read tile_id year <<< $(job_to_year_tile "$job_id")
            tmf_file=$(printf "%s/tmf_%03d_%d.csv.gz" "$TMF_PATH" "$tile_id" "$year")
            if [[ -f "$tmf_file" ]]; then
                rm -f "$tmf_file"
                deleted_count=$((deleted_count + 1))
                echo "  Deleted: $tmf_file"
            fi
            ;;
        2a)
            # Stage 2a: Consolidated tile files
            cons_file=$(printf "%s/tile_%03d.parquet" "$CONSOLIDATED_PATH" "$job_id")
            if [[ -f "$cons_file" ]]; then
                rm -f "$cons_file"
                deleted_count=$((deleted_count + 1))
                echo "  Deleted: $cons_file"
            fi
            ;;
        3)
            # Stage 3: Interior classification files
            int_file=$(printf "%s/interior_tile_%03d.parquet" "$INTERIOR_PATH" "$job_id")
            if [[ -f "$int_file" ]]; then
                rm -f "$int_file"
                deleted_count=$((deleted_count + 1))
                echo "  Deleted: $int_file"
            fi
            ;;
        4)
            # Stage 4: Frontier files
            front_file=$(printf "%s/frontier_tile_%03d.parquet" "$FRONTIER_PATH" "$job_id")
            if [[ -f "$front_file" ]]; then
                rm -f "$front_file"
                deleted_count=$((deleted_count + 1))
                echo "  Deleted: $front_file"
            fi
            ;;
        5)
            # Stage 5: Covariate files
            cov_file=$(printf "%s/covariates_tile_%03d.parquet" "$COVARIATES_PATH" "$job_id")
            if [[ -f "$cov_file" ]]; then
                rm -f "$cov_file"
                deleted_count=$((deleted_count + 1))
                echo "  Deleted: $cov_file"
            fi
            ;;
        6)
            # Stage 6: Final assembly (single file)
            final_file="Data/output/TMF_5km_panel.parquet"
            if [[ -f "$final_file" ]]; then
                rm -f "$final_file"
                deleted_count=$((deleted_count + 1))
                echo "  Deleted: $final_file"
            fi
            ;;
    esac
done

echo ""
echo "Deleted $deleted_count files"

if [[ "$DELETE_ONLY" == true ]]; then
    echo ""
    echo "========================================"
    echo "--delete-only: Skipping job submission"
    echo "========================================"
    exit 0
fi

# Determine which script to run
case "$STAGE" in
    0)
        SCRIPT="code/bash/0_grid.sh"
        ;;
    0c)
        SCRIPT="code/bash/0c_yields.sh"
        ;;
    1)
        SCRIPT="code/bash/1_extract_TMF.sh"
        ;;
    2a)
        SCRIPT="code/bash/2a_consolidate_tile.sh"
        ;;
    3)
        SCRIPT="code/bash/3_interior.sh"
        ;;
    4)
        SCRIPT="code/bash/4_frontier.sh"
        ;;
    5)
        SCRIPT="code/bash/5_covariates.sh"
        ;;
    6)
        SCRIPT="code/bash/6_assemble.sh"
        ;;
    *)
        echo "Unknown stage: $STAGE"
        echo "Valid stages: 0, 0c, 1, 2a, 3, 4, 5, 6"
        exit 1
        ;;
esac

if [[ ! -f "$SCRIPT" ]]; then
    echo "Script not found: $SCRIPT"
    exit 1
fi

echo "Script: $SCRIPT"
echo "Job IDs: $JOB_LIST"
echo ""

# Submit with specific array indices
# We need to modify the array specification
# Create a temporary script with the correct array
TEMP_SCRIPT=$(mktemp)

# Copy the script but replace the array line
sed "s/#SBATCH --array=.*/#SBATCH --array=$JOB_LIST/" "$SCRIPT" > "$TEMP_SCRIPT"

echo "Submitting jobs..."
JOB_ID=$(sbatch --parsable $SBATCH_COMMON "$TEMP_SCRIPT")

rm "$TEMP_SCRIPT"

echo ""
echo "========================================"
echo "Submitted job array: $JOB_ID"
echo "Jobs: $JOB_LIST"
echo "========================================"
echo "Monitor with: squeue -u \$USER"
echo "Logs in: LOGS/"
