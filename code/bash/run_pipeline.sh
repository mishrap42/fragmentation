#!/bin/bash
# ==============================================================================
# TMF 1km Pipeline - Master Run Script
# ==============================================================================
# This script submits all pipeline stages in order with proper dependencies.
#
# Usage:
#   ./bash/run_pipeline.sh           # Run full pipeline
#   ./bash/run_pipeline.sh 3         # Start from stage 3
#   ./bash/run_pipeline.sh 1 3       # Run stages 1-3 only
# ==============================================================================

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Parse arguments
START_STAGE=${1:-0}
END_STAGE=${2:-6}

echo "========================================"
echo "TMF 1km Data Pipeline"
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
        0) script="bash/0_grid.sh" ;;
        1) script="bash/1_extract_TMF.sh" ;;
        2) script="bash/2a_consolidate_tile.sh" ;;
        3) script="bash/3_interior.sh" ;;
        4) script="bash/4_frontier.sh" ;;
        5) script="bash/5_wdpa.sh" ;;
        6) script="bash/6_assemble.sh" ;;
    esac

    if [[ ! -f "$script" ]]; then
        echo -e "${RED}Error: Script not found: $script${NC}"
        exit 1
    fi
done

# Create LOGS directory
mkdir -p LOGS

# Function to submit job and get job ID
submit_job() {
    local script=$1
    local dependency=$2
    local description=$3

    echo -e "${YELLOW}Submitting: $description${NC}"

    if [[ -n "$dependency" ]]; then
        job_id=$(sbatch --dependency=afterok:$dependency "$script" | awk '{print $4}')
    else
        job_id=$(sbatch "$script" | awk '{print $4}')
    fi

    echo -e "${GREEN}  Submitted job: $job_id${NC}"
    echo "$job_id"
}

# Track job IDs for dependencies
PREV_JOB=""

# ==============================================================================
# Stage 0: Grid Generation
# ==============================================================================
if [[ $START_STAGE -le 0 && $END_STAGE -ge 0 ]]; then
    JOB_0=$(submit_job "bash/0_grid.sh" "$PREV_JOB" "Stage 0: Grid Generation (86 jobs)")
    PREV_JOB=$JOB_0
fi

# ==============================================================================
# Stage 1: TMF Extraction
# ==============================================================================
if [[ $START_STAGE -le 1 && $END_STAGE -ge 1 ]]; then
    JOB_1=$(submit_job "bash/1_extract_TMF.sh" "$PREV_JOB" "Stage 1: TMF Extraction (3010 jobs)")
    PREV_JOB=$JOB_1
fi

# ==============================================================================
# Stage 2a: Tile Consolidation
# ==============================================================================
if [[ $START_STAGE -le 2 && $END_STAGE -ge 2 ]]; then
    JOB_2A=$(submit_job "bash/2a_consolidate_tile.sh" "$PREV_JOB" "Stage 2a: Tile Consolidation (86 jobs)")

    # Stage 2b depends on 2a
    JOB_2B=$(submit_job "bash/2b_consolidate_global.sh" "$JOB_2A" "Stage 2b: Global Consolidation (1 job)")
    PREV_JOB=$JOB_2B
fi

# ==============================================================================
# Stage 3: Interior Classification
# ==============================================================================
if [[ $START_STAGE -le 3 && $END_STAGE -ge 3 ]]; then
    JOB_3=$(submit_job "bash/3_interior.sh" "$PREV_JOB" "Stage 3: Interior Classification (86 jobs)")
    PREV_JOB=$JOB_3
fi

# ==============================================================================
# Stage 4: Frontier Calculation
# ==============================================================================
if [[ $START_STAGE -le 4 && $END_STAGE -ge 4 ]]; then
    JOB_4=$(submit_job "bash/4_frontier.sh" "$PREV_JOB" "Stage 4: Frontier Calculation (86 jobs)")
    PREV_JOB=$JOB_4
fi

# ==============================================================================
# Stage 5: WDPA Extraction
# ==============================================================================
if [[ $START_STAGE -le 5 && $END_STAGE -ge 5 ]]; then
    JOB_5=$(submit_job "bash/5_wdpa.sh" "$PREV_JOB" "Stage 5: WDPA Extraction (87 jobs)")
    PREV_JOB=$JOB_5
fi

# ==============================================================================
# Stage 6: Final Assembly
# ==============================================================================
if [[ $START_STAGE -le 6 && $END_STAGE -ge 6 ]]; then
    JOB_6=$(submit_job "bash/6_assemble.sh" "$PREV_JOB" "Stage 6: Final Assembly (1 job)")
fi

echo ""
echo "========================================"
echo -e "${GREEN}All jobs submitted!${NC}"
echo "========================================"
echo "Monitor progress with: squeue -u \$USER"
echo "Check logs in: LOGS/"
echo "========================================"
