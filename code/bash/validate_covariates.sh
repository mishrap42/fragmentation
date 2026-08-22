#!/bin/bash -l
#SBATCH --time=01:00:00
#SBATCH --job-name=VALIDATE_COV
#SBATCH --account=mishralab
#SBATCH --partition=expansion
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --output=LOGS/%x.%j.out
#SBATCH --error=LOGS/%x.%j.err

# ==============================================================================
# Validate Stage 5 Covariate Extraction
# Creates visualizations and summary statistics for all covariate columns
# ==============================================================================
#
# Usage:
#   sbatch code/bash/validate_covariates.sh           # Validate all tiles
#   sbatch code/bash/validate_covariates.sh 1         # Validate tile 1 only
#
# Output: LOGS/validation/covariates/
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
echo "COVARIATE VALIDATION"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Start Time: $(date)"
echo "========================================"

# Determine project root
if [[ -n "$SLURM_SUBMIT_DIR" ]]; then
    PROJECT_ROOT="$SLURM_SUBMIT_DIR"
else
    PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
fi

cd "$PROJECT_ROOT" || exit 1
echo "Project root: $PROJECT_ROOT"

# Get optional tile argument
TILE_ARG="${1:-}"

R_SCRIPT="code/build/validate_covariates.R"

if [[ -f "$R_SCRIPT" ]]; then
    if [[ -n "$TILE_ARG" ]]; then
        echo "Validating tile: $TILE_ARG"
        Rscript "$R_SCRIPT" "$TILE_ARG"
    else
        echo "Validating all tiles"
        Rscript "$R_SCRIPT"
    fi
    EXIT_CODE=$?
    echo "R script exited with code: $EXIT_CODE"
else
    echo "Error: R script not found at $R_SCRIPT"
    exit 1
fi

echo "========================================"
echo "End Time: $(date)"
echo "========================================"

# List output files
echo ""
echo "Output files:"
ls -la LOGS/validation/covariates/ 2>/dev/null || echo "No output files found"

exit $EXIT_CODE
