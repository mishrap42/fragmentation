#!/bin/bash -l
#SBATCH --time=96:00:00
#SBATCH --job-name=SYNTH_PA
#SBATCH --account=mishralab
#SBATCH --partition=expansion
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --array=1
#SBATCH --mem=64G
#SBATCH --output=LOGS/%x.%A_%a.out
#SBATCH --error=LOGS/%x.%A_%a.err

# ==============================================================================
# 4_synthetic_control.sh - Per-project synthetic control (tidysynth)
# Each array task: one country. Loops over the country's protected areas,
# fits a synthetic control for each PA taken as a whole against the country's
# never-protected cells, and records project-level treatment effects + donor
# weights with placebo-based inference.
#
# Array size: 128 (overestimate; R script handles out-of-range gracefully).
# NOTE: with no donor cap (default), placebo inference runs one SC optimization
#       per donor cell. Large countries can be slow / memory-heavy. To cap the
#       donor pool, export SYNTH_MAX_DONORS before sbatch (e.g. =500).
# ==============================================================================

# sbatch runs a COPY of this script from /var/spool/slurmd/<job>/, so
# ${BASH_SOURCE[0]} does NOT point at code/bash/ under SLURM. Resolve from
# SLURM_SUBMIT_DIR (submit from the project root), and fall back to the script
# location for plain `bash <script>` runs. The rerun_* helpers submit a sed'd
# copy from $(mktemp), which is exactly why the fallback cannot be trusted here.
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  HERE="$SLURM_SUBMIT_DIR/code/bash"
else
  HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
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
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Submit Directory: $SLURM_SUBMIT_DIR"
echo "Node: $SLURMD_NODENAME"
echo "Start Time: $(date)"
echo "SYNTH_OUTCOME: ${SYNTH_OUTCOME:-forest_cover (default)}"
echo "SYNTH_MAX_DONORS: ${SYNTH_MAX_DONORS:-Inf (no cap, default)}"
echo "========================================"

if [[ -n "$SLURM_SUBMIT_DIR" ]]; then
    PROJECT_ROOT="$SLURM_SUBMIT_DIR"
else
    PROJECT_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
fi

echo "Project root: $PROJECT_ROOT"
cd "$PROJECT_ROOT" || exit 1

R_SCRIPT="code/analysis/4_synthetic_control.R"

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
