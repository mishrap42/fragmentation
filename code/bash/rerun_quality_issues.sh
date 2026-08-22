#!/bin/bash -l
# ==============================================================================
# Rerun tiles with data quality issues
# ==============================================================================
# Usage: bash code/bash/rerun_quality_issues.sh [stage]
#   stage: 1, 2a, or "all" (default: all)
#
# First run: bash code/bash/diagnose_data_quality.sh
# ==============================================================================

set -e

# Resource flags and the LOGS directory come from config.sh. A parent's
# account/partition/qos are NOT inherited by the jobs it submits, so every
# sbatch below must forward $SBATCH_COMMON explicitly.
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
[ -f "$HERE/config.sh" ] || { echo "config.sh not found at $HERE" >&2; exit 1; }
source "$HERE/config.sh"

STAGE="${1:-all}"

STAGE1_FILE="LOGS/rerun/stage1_quality_issues.txt"
STAGE2A_FILE="LOGS/rerun/stage2a_quality_issues.txt"

# Check files exist
if [[ ! -f "$STAGE1_FILE" ]] && [[ ! -f "$STAGE2A_FILE" ]]; then
    echo "No quality issues files found."
    echo "Run diagnosis first: bash code/bash/diagnose_data_quality.sh"
    exit 1
fi

submit_stage1() {
    if [[ ! -f "$STAGE1_FILE" ]]; then
        echo "No Stage 1 issues file: $STAGE1_FILE"
        return
    fi

    N_JOBS=$(wc -l < "$STAGE1_FILE" | tr -d ' ')
    if [[ "$N_JOBS" -eq 0 ]]; then
        echo "No Stage 1 jobs to rerun"
        return
    fi

    echo "========================================"
    echo "Stage 1: Rerunning $N_JOBS TMF extraction jobs"
    echo "========================================"

    # Read job IDs
    mapfile -t JOB_IDS < "$STAGE1_FILE"
    JOB_LIST=$(IFS=,; echo "${JOB_IDS[*]}")

    # Create temp script with updated array
    TEMP_SCRIPT=$(mktemp)
    sed "s/#SBATCH --array=.*/#SBATCH --array=$JOB_LIST/" code/bash/1_extract_TMF.sh > "$TEMP_SCRIPT"

    STAGE1_JOB=$(sbatch --parsable $SBATCH_COMMON "$TEMP_SCRIPT")
    rm "$TEMP_SCRIPT"

    echo "Submitted Stage 1 job array: $STAGE1_JOB"
    echo "$STAGE1_JOB" > LOGS/rerun/stage1_quality_jobid.txt
}

submit_stage2a() {
    if [[ ! -f "$STAGE2A_FILE" ]]; then
        echo "No Stage 2a issues file: $STAGE2A_FILE"
        return
    fi

    N_JOBS=$(wc -l < "$STAGE2A_FILE" | tr -d ' ')
    if [[ "$N_JOBS" -eq 0 ]]; then
        echo "No Stage 2a jobs to rerun"
        return
    fi

    echo "========================================"
    echo "Stage 2a: Rerunning $N_JOBS tile consolidation jobs"
    echo "========================================"

    # Read tile IDs
    mapfile -t TILE_IDS < "$STAGE2A_FILE"
    TILE_LIST=$(IFS=,; echo "${TILE_IDS[*]}")

    # Create temp script with updated array
    TEMP_SCRIPT=$(mktemp)
    sed "s/#SBATCH --array=.*/#SBATCH --array=$TILE_LIST/" code/bash/2a_consolidate_tile.sh > "$TEMP_SCRIPT"

    # Add dependency if Stage 1 was submitted
    if [[ -f "LOGS/rerun/stage1_quality_jobid.txt" ]]; then
        STAGE1_JOB=$(cat LOGS/rerun/stage1_quality_jobid.txt)
        STAGE2A_JOB=$(sbatch --parsable $SBATCH_COMMON --dependency=afterok:$STAGE1_JOB "$TEMP_SCRIPT")
        echo "Submitted Stage 2a (depends on $STAGE1_JOB): $STAGE2A_JOB"
    else
        STAGE2A_JOB=$(sbatch --parsable $SBATCH_COMMON "$TEMP_SCRIPT")
        echo "Submitted Stage 2a: $STAGE2A_JOB"
    fi

    rm "$TEMP_SCRIPT"
}

case "$STAGE" in
    1)
        submit_stage1
        ;;
    2a)
        submit_stage2a
        ;;
    all)
        submit_stage1
        submit_stage2a
        ;;
    *)
        echo "Unknown stage: $STAGE"
        echo "Usage: bash code/bash/rerun_quality_issues.sh [1|2a|all]"
        exit 1
        ;;
esac

echo ""
echo "========================================"
echo "Monitor with: squeue -u \$USER"
echo "========================================"
