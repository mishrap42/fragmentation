#!/bin/bash
# ==============================================================================
# TMF 5km Pipeline - Master Run Script
# ==============================================================================
# This script submits all pipeline stages in order with proper dependencies.
#
# PREREQUISITE: Run 0a_wdpa.sh first to create wdpa_cleaned.gpkg
#   sbatch code/bash/0a_wdpa.sh
#   (Wait for completion before running this script)
#
# Usage:
#   ./code/bash/run_pipeline.sh           # Run full pipeline (stages 0-5)
#   ./code/bash/run_pipeline.sh 3         # Start from stage 3
#   ./code/bash/run_pipeline.sh 1 3       # Run stages 1-3 only
#
# Pipeline stages (Stage 5 removed - WDPA now in grid):
#   0  - Grid Generation (creates 5km cells cut on country + WDPA boundaries)
#   0c - Yield Extraction (runs automatically after Stage 0)
#   1  - TMF Extraction
#   2  - Tile Consolidation
#   3  - Interior Classification
#   4  - Frontier Calculation
#   5  - Final Assembly (was Stage 6)
# ==============================================================================

set -e

# Resource flags and the LOGS directory come from config.sh. A parent's
# account/partition/qos are NOT inherited by the jobs it submits, so every
# sbatch below must forward $SBATCH_COMMON explicitly.
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
[ -f "$HERE/config.sh" ] || { echo "config.sh not found at $HERE" >&2; exit 1; }
source "$HERE/config.sh"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Parse arguments
START_STAGE=${1:-0}
END_STAGE=${2:-5}

echo "========================================"
echo "TMF 5km Data Pipeline"
echo "========================================"
echo "Starting from stage: $START_STAGE"
echo "Ending at stage: $END_STAGE"
echo "========================================"

# Change to project root
cd "$(dirname "$0")/.."
PROJECT_ROOT=$(pwd)
echo "Project root: $PROJECT_ROOT"

# Check that all bash scripts exist
for stage in $(seq $START_STAGE $END_STAGE); do
    case $stage in
        0) script="code/bash/0_grid.sh" ;;
        1) script="code/bash/1_extract_TMF.sh" ;;
        2) script="code/bash/2a_consolidate_tile.sh" ;;
        3) script="code/bash/3_interior.sh" ;;
        4) script="code/bash/4_frontier.sh" ;;
        5) script="code/bash/6_assemble.sh" ;;
    esac

    if [[ ! -f "$script" ]]; then
        echo -e "${RED}Error: Script not found: $script${NC}"
        exit 1
    fi
done

# Check WDPA prerequisite
if [[ $START_STAGE -le 0 ]]; then
    WDPA_FILE="Data/build/wdpa_cleaned.gpkg"
    if [[ ! -f "$WDPA_FILE" ]]; then
        echo -e "${YELLOW}WARNING: WDPA file not found: $WDPA_FILE${NC}"
        echo -e "${YELLOW}Run 'sbatch code/bash/0a_wdpa.sh' first and wait for completion.${NC}"
        read -p "Continue anyway? (y/N) " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            exit 1
        fi
    fi
fi

# LOGS must exist BEFORE the first sbatch: SLURM opens the --output file
# before the job script runs.
frag_ensure_logs

# Function to submit job and get job ID
submit_job() {
    local script=$1
    local dependency=$2
    local description=$3

    echo -e "${YELLOW}Submitting: $description${NC}"

    if [[ -n "$dependency" ]]; then
        job_id=$(sbatch --parsable $SBATCH_COMMON --dependency=afterok:$dependency "$script")
    else
        job_id=$(sbatch --parsable $SBATCH_COMMON "$script")
    fi

    echo -e "${GREEN}  Submitted job: $job_id${NC}"
    echo "$job_id"
}

# Track job IDs for dependencies
PREV_JOB=""

# ==============================================================================
# Stage 0: Grid Generation (5km cells cut on country + WDPA boundaries)
# ==============================================================================
if [[ $START_STAGE -le 0 && $END_STAGE -ge 0 ]]; then
    JOB_0=$(submit_job "code/bash/0_grid.sh" "$PREV_JOB" "Stage 0: Grid Generation (344 sub-tile jobs)")
    PREV_JOB=$JOB_0

    # Stage 0c: EarthStat Extraction (runs after grid generation)
    if [[ -f "code/bash/0c_yields.sh" ]]; then
        JOB_0C=$(submit_job "code/bash/0c_yields.sh" "$JOB_0" "Stage 0c: EarthStat Extraction (344 sub-tile jobs)")
        # Note: Stage 1 depends on Stage 0, not 0c (yields are optional metadata)

        # Stage 1b: consolidate 0c into the cropland cross-section. Chained here
        # rather than near Stage 6 because it belongs to the extraction side of
        # the build: 0c is its only input, and nothing downstream consumes it.
        # It runs in parallel with the whole TMF chain instead of behind it.
        if [[ -f "code/bash/1b_assemble_cropland.sh" ]]; then
            submit_job "code/bash/1b_assemble_cropland.sh" "$JOB_0C" \
                "Stage 1b: Cropland Cross-Section (1 job)" > /dev/null
        fi
    fi
fi

# ==============================================================================
# Stage 1: TMF Extraction
# ==============================================================================
if [[ $START_STAGE -le 1 && $END_STAGE -ge 1 ]]; then
    JOB_1=$(submit_job "code/bash/1_extract_TMF.sh" "$PREV_JOB" "Stage 1: TMF Extraction (2924 jobs)")
    PREV_JOB=$JOB_1
fi

# ==============================================================================
# Stage 2a: Tile Consolidation
# ==============================================================================
if [[ $START_STAGE -le 2 && $END_STAGE -ge 2 ]]; then
    JOB_2A=$(submit_job "code/bash/2a_consolidate_tile.sh" "$PREV_JOB" "Stage 2a: Tile Consolidation (86 jobs)")
    PREV_JOB=$JOB_2A

    # NOTE: Stage 2b (Global Consolidation) is SKIPPED - no longer required.
    # Stage 6 now reads directly from tile-level parquet files using Arrow,
    # which avoids R's 2.1 billion row limit (~3B rows globally).
fi

# ==============================================================================
# Stage 3: Interior Classification
# ==============================================================================
if [[ $START_STAGE -le 3 && $END_STAGE -ge 3 ]]; then
    JOB_3=$(submit_job "code/bash/3_interior.sh" "$PREV_JOB" "Stage 3: Interior Classification (86 jobs)")
    PREV_JOB=$JOB_3
fi

# ==============================================================================
# Stage 4: Frontier Calculation
# ==============================================================================
if [[ $START_STAGE -le 4 && $END_STAGE -ge 4 ]]; then
    JOB_4=$(submit_job "code/bash/4_frontier.sh" "$PREV_JOB" "Stage 4: Frontier Calculation (86 jobs)")
    PREV_JOB=$JOB_4
fi

# ==============================================================================
# Stage 5: Final Assembly (formerly Stage 6 - WDPA stage removed)
# ==============================================================================
if [[ $START_STAGE -le 5 && $END_STAGE -ge 5 ]]; then
    JOB_5=$(submit_job "code/bash/6_assemble.sh" "$PREV_JOB" "Stage 5: Final Assembly (1 job)")
fi

echo ""
echo "========================================"
echo -e "${GREEN}All jobs submitted!${NC}"
echo "========================================"
echo "Monitor progress with: squeue -u \$USER"
echo "Check logs in: LOGS/"
echo "========================================"
