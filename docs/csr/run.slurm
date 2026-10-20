#!/bin/bash
#SBATCH --job-name=csr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --account=a_imb_ccbcd
#SBATCH --mail-type=END,FAIL
# No --mail-user here on purpose: this repo is public. See README.md.
#SBATCH --output=protdft/%x-%j.out
#SBATCH --error=protdft/%x-%j.err
#
# One job script for every mcube stage. Resources come from the sbatch command
# line, which overrides the defaults above, so there is nothing to edit per run.
#
#   sbatch -J <name> [resources] run.slurm <profile> <stage> <section> [start.pdb]
#
#   profile  xtb | dft   threading only, see below
#   stage    qmmm | scan | neb | ts
#   section  the config.yaml section, e.g. scan_znco
#   start    optional --pdb; omit to let mcube take it from the previous stage
#
# Measured efficiency on this system (seff), which is why profile matters:
#   DFT single points, 16 CPUs : 96% CPU, 56 GB   -> scales, give it cores
#   xTB/MM MD,          8 CPUs : 30% CPU, 0.5 GB  -> does not scale past ~4
# So run DFT wide, and make xTB work go faster by running points in parallel
# rather than by adding threads.
#
#   sbatch -J coscan --cpus-per-task=8  --mem=16G run.slurm xtb scan scan_znco start.pdb
#   sbatch -J spdft  --cpus-per-task=64 --mem=192G --time=12:00:00 \
#          run.slurm dft qmmm qmmm_sp_dft geom.pdb
set -euo pipefail
ulimit -s unlimited
export OMP_STACKSIZE=1G
export PYTHONUNBUFFERED=1

PROFILE=${1:?profile: xtb or dft}
STAGE=${2:?stage: qmmm, scan, neb or ts}
SECTION=${3:?config.yaml section}
START=${4:-}

case "$PROFILE" in
  xtb)
    # xTB runs as a separate process, so binding its threads is safe and is
    # worth 2-9x on SMT nodes. 8 allocated CPUs is 4 physical cores.
    export OMP_NUM_THREADS=$(( ${SLURM_CPUS_PER_TASK:-8} / 2 ))
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    ;;
  dft)
    # PySCF runs in-process and sets its own thread count from dft_cores or
    # sp_cores. A conflicting OMP_NUM_THREADS makes it warn and can clamp it,
    # and process-wide binding would pin its BLAS threads to one core.
    export OMP_NUM_THREADS=1
    ;;
  *) echo "unknown profile: $PROFILE (xtb or dft)" >&2; exit 2 ;;
esac

module load miniforge/26.1.0-0
source "$ROOTMINIFORGE/etc/profile.d/conda.sh"
conda activate ASH

ROOT=$PWD
OUT=protdft/$SECTION
mkdir -p "$OUT"

if [ -n "$START" ]; then
  START_ABS=$(readlink -f "$START")
  test -f "$START_ABS" || { echo "missing start structure: $START" >&2; exit 1; }
  ARGS=(--pdb "$START_ABS")
else
  ARGS=()
fi

cd "$OUT"
date
mcube "$STAGE" "$ROOT/config.yaml" --section "$SECTION" "${ARGS[@]}"
date
echo "[DONE] $SECTION"

# Always report the site afterwards: cheap, and it is the first thing anyone
# asks of a finished structure.
for PDB in QM_MM_optimized.pdb QM_MM_lastframe.pdb; do
  [ -f "$PDB" ] && python "$ROOT/../mcube/docs/csr/site.py" sphere "$PDB" \
      --site="$ROOT/../mcube/docs/csr/site.yaml" || true
done
