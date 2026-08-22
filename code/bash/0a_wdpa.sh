#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH --job-name=TMF_WDPA_CLEAN
#SBATCH --account=mishralab
#SBATCH --partition=expansion
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --output=LOGS/%x.%j.out
#SBATCH --error=LOGS/%x.%j.err
#SBATCH --requeue
#SBATCH --open-mode=append

# ==============================================================================
# RESOURCE / RESILIENCE NOTES
#
#   --time=48:00:00  Job 1067292 died at the old 4h limit still inside the
#                    conversion, so 4h is known-insufficient and the true cost
#                    is not yet known. expansion allows 14 days; there is no
#                    reason to run this tight. Phase A tuning should cut it
#                    substantially, but the limit stays generous until we have
#                    one completed run to size against.
#
#   --mem=120G       The validation phase (st_read + st_is_valid +
#                    st_make_valid on ~255k polygons) has never actually run -
#                    the job has never got past conversion. It is the peak and
#                    it is untested, so 64G -> 120G. If this pends forever,
#                    check the partition ceiling: sinfo -p expansion -o "%n %m"
#
#   --requeue        Node failure or preemption re-runs the job automatically.
#                    Cheap now that the R script checkpoints: a requeue after
#                    Phase A resumes at Phase B rather than restarting.
#
#   --open-mode=append   Without this a requeue truncates the log, destroying
#                    the record of why the first attempt died.
# ==============================================================================

# ==============================================================================
# STAGE 0a: WDPA Preprocessing
# Run ONCE before the main pipeline to create cleaned WDPA file
# This is a single job (no array) that takes 1-2 hours
# ==============================================================================

# sbatch runs a COPY of this script from /var/spool/slurmd/<job>/, so
# ${BASH_SOURCE[0]} does NOT point at code/bash/ under SLURM. Resolve from
# SLURM_SUBMIT_DIR (submit from the project root), and fall back to the script
# location for plain `bash <script>` runs. The rerun_* helpers submit a sed'd
# copy from $(mktemp), which is exactly why the fallback cannot be trusted here.
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

R_SCRIPT="code/build/0a_clean_wdpa.R"

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
