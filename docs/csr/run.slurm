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
#   profile  xtb | xtb-parallel | dft   thread placement only, see below
#   stage    preprocess | md | qmmm | sp | scan | neb | ts
#   section  the config.yaml section, e.g. scan_znco; for preprocess and md,
#            which have no sections, it just names the output directory
#   start    optional --pdb; omit to let mcube take it from the previous stage
#
# Measured efficiency on this system (seff), which is why profile matters:
#   DFT single points (sp),  16 : 96% CPU, 56 GB   -> scales, give it cores
#   xTB/MM MD,          8 CPUs : 30% CPU, 0.5 GB  -> does not scale past ~4
# So run DFT wide, and make xTB work go faster by running points in parallel
# rather than by adding threads.
#
#   sbatch -J coscan --cpus-per-task=8  --mem=16G run.slurm xtb scan scan_znco start.pdb
#   sbatch -J spdft  --cpus-per-task=64 --mem=192G --time=12:00:00 \
#          run.slurm dft sp sp_sp_dft geom.pdb
set -euo pipefail
ulimit -s unlimited
export OMP_STACKSIZE=1G
export PYTHONUNBUFFERED=1

PROFILE=${1:?profile: xtb, xtb-parallel or dft}
STAGE=${2:?stage: preprocess, md, qmmm, sp, scan, neb or ts}
SECTION=${3:?config.yaml section, or the stage name for preprocess/md}
START=${4:-}

# preprocess and md have no --section: they write to the one place in the
# config. For them the third argument only names the output directory.
case "$STAGE" in
  preprocess|md) SECTION_ARGS=();            OUT=${STAGE}_results ;;
  *)             SECTION_ARGS=(--section "$SECTION"); OUT=protdft/$SECTION ;;
esac

# Thread placement lives here, not in mcube: it has to be in the environment
# before ASH, OpenMM or xTB start, and the right policy differs per stage.
# mcube only reads a thread *count* (a stage's `cores`, else $OMP_NUM_THREADS,
# else $SLURM_CPUS_PER_TASK).
case "$PROFILE" in
  xtb)
    # One xTB process at a time: bind one thread per physical core. Worth 2-9x
    # on SMT nodes, where N allocated CPUs are N/2 physical cores.
    # Stages: md, qmmm, scan, and ts without freq_workers.
    export OMP_NUM_THREADS=$(( ${SLURM_CPUS_PER_TASK:-8} / 2 ))
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    ;;
  xtb-parallel)
    # Several xTB processes at once (neb with parallel images, ts with
    # freq_workers > 1). Binding here is actively harmful: every worker would
    # pin its threads to the same first cores. Leave placement to the OS.
    export OMP_NUM_THREADS=$(( ${SLURM_CPUS_PER_TASK:-8} / 2 ))
    unset OMP_PROC_BIND OMP_PLACES
    ;;
  dft)
    # PySCF runs in-process and sets its own thread count from dft_cores or
    # sp_cores. A conflicting OMP_NUM_THREADS makes it warn and can clamp it,
    # and process-wide binding pins its BLAS threads to one core (measured:
    # ~18% CPU per thread). Stages: sp, and anything with qm_engine: dft.
    export OMP_NUM_THREADS=1
    unset OMP_PROC_BIND OMP_PLACES
    ;;
  *) echo "unknown profile: $PROFILE (xtb, xtb-parallel or dft)" >&2; exit 2 ;;
esac

module load miniforge/26.1.0-0
source "$ROOTMINIFORGE/etc/profile.d/conda.sh"
conda activate ASH

ROOT=$PWD
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
mcube "$STAGE" "$ROOT/config.yaml" "${SECTION_ARGS[@]}" "${ARGS[@]}"
date
echo "[DONE] $SECTION"

# Always report the site afterwards: cheap, and it is the first thing anyone
# asks of a finished structure.
for PDB in QM_MM_optimized.pdb QM_MM_lastframe.pdb; do
  [ -f "$PDB" ] && python "$ROOT/../mcube/docs/csr/site.py" sphere "$PDB" \
      --site="$ROOT/../mcube/docs/csr/site.yaml" || true
done
