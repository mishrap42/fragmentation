# config.sh — single source of truth for the Fragmentation (TMF 5km) SLURM
# jobs on the Caltech Resnick HPC cluster. Sourced by every script in
# code/bash/ and code/bash/analysis/.
#
# Modelled on the Global Forest Repo's code/bashfiles/config.sh, and it shares
# that project's R library tree (see R_LIBS_GROUP below) — the two pipelines
# use the same system-R + system-geo toolchain, so one build serves both.
#
# Usage, at the top of a job script (after the #SBATCH block):
#
#     if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
#       HERE="$SLURM_SUBMIT_DIR/code/bash"
#     else
#       HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
#     fi
#     source "$HERE/config.sh"
#     frag_load_modules
#
# Submitter scripts additionally pass $SBATCH_COMMON to every child sbatch.

# ---------------------------------------------------------------------------
# Cluster
# ---------------------------------------------------------------------------
REPO=/resnick/groups/MishraLab/Fragmentation   # already space-free; no symlink
ACCOUNT=mishralab
PARTITION=expansion                            # general CPU, 14-day wall limit
QOS=normal                                     # debug|normal|long

# Resource flags every job needs. Resnick REQUIRES these; the previous
# (Dartmouth) versions of these scripts carried none, which is why they
# cannot be submitted here unmodified. Submitter scripts must forward this to
# the children they sbatch — a parent's account/partition is NOT inherited.
SBATCH_COMMON="--account=$ACCOUNT --partition=$PARTITION --qos=$QOS"

# ---------------------------------------------------------------------------
# Toolchain — the RHEL9 SYSTEM stack, with udunits as the only module.
# See code/bashfiles/config.sh in the Global Forest Repo for the full rationale
# (system R 4.5.1 + system gdal/geos/proj/netcdf headers; the spack gdal and
# netcdf-c modules are phantoms). The one gap is /usr/include/udunits2.h,
# which is absent, so `units` — and therefore `sf` — needs the module below.
# ---------------------------------------------------------------------------
UDUNITS_MODULE=udunits/2.2.28-gcc-11.3.1-opodm73

# Shared R library, group space. Binaries here are linked against the SYSTEM
# gdal/geos/proj/netcdf and the module udunits, so this tree is NOT portable
# from another cluster. Built by the Global Forest Repo's install_packages.R;
# this project's extra packages come from code/build/install_packages.R.
R_LIBS_GROUP=/resnick/groups/MishraLab/Rlib/4.5

frag_load_modules() {
  module purge
  module load "$UDUNITS_MODULE"

  # The udunits modulefile prepends PATH and CMAKE_PREFIX_PATH ONLY — it sets
  # no CPATH, LD_LIBRARY_PATH or PKG_CONFIG_PATH. Derive the spack prefix from
  # the PATH entry the module just added, then hand `units` the locations
  # directly (its configure reads UDUNITS2_INCLUDE/UDUNITS2_LIBS) and put the
  # library on the runtime path for library(units) later.
  local udunits_bin udunits_prefix
  udunits_bin=$(echo "$PATH" | tr ':' '\n' | grep -m1 'udunits')
  udunits_prefix=${udunits_bin%/bin}

  if [ -z "$udunits_prefix" ] || [ ! -f "$udunits_prefix/include/udunits2.h" ]; then
    echo "ERROR: udunits2.h not found under '$udunits_prefix'." >&2
    echo "       Module $UDUNITS_MODULE loaded but its layout changed;" >&2
    echo '       the units package (and therefore sf) will not build.' >&2
    return 1
  fi

  export UDUNITS2_INCLUDE="$udunits_prefix/include"
  export UDUNITS2_LIBS="$udunits_prefix/lib"
  export LD_LIBRARY_PATH="$udunits_prefix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
  export PKG_CONFIG_PATH="$udunits_prefix/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"

  # Everything else is the RHEL9 system stack: /usr/bin/R 4.5.1 and the
  # gdal/geos/proj/netcdf headers in /usr/include. Nothing to load.
  export R_LIBS_USER="$R_LIBS_GROUP"
  mkdir -p "$R_LIBS_USER"

  echo "== modules =="
  module list 2>&1
  echo "== R            : $(command -v R)  $(R --version 2>/dev/null | head -1)"
  echo "== GDAL         : $(gdalinfo --version 2>/dev/null)"
  echo "== UDUNITS2     : $udunits_prefix"
  echo "== R_LIBS_USER  : $R_LIBS_USER"
}

# Create the log directory BEFORE submitting. SLURM opens the --output file
# before the job script runs, so a `mkdir -p LOGS` inside the script is too
# late: if LOGS/ is missing the job dies with no output and no clear cause.
frag_ensure_logs() {
  mkdir -p "${1:-LOGS}"
}
