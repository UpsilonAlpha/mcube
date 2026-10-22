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
# One job script for every mcube stage.
#
#   sbatch -J <name> [resources] run.slurm <stage> <section> [start.pdb]
#
#   stage    preprocess | md | qmmm | sp | scan | neb | ts
#   section  the config.yaml section, e.g. scan_znco or md_nac. preprocess has
#            none, so there it just names the output directory
#   start    optional --pdb; omit to let mcube take it from the previous stage
#
# Examples:
#
#   sbatch -J prep   --cpus-per-task=2  --mem=4G                 run.slurm preprocess preprocess
#   sbatch -J md     --cpus-per-task=16 --mem=16G                run.slurm md md
#   sbatch -J mdnac  --cpus-per-task=16 --mem=16G                run.slurm md md_nac
#   sbatch -J react  --cpus-per-task=8  --mem=16G                run.slurm qmmm qmmm_reactant
#   sbatch -J spR    --cpus-per-task=32 --mem=96G --time=12:00:00 run.slurm sp sp_reactant \
#          protdft/qmmm_reactant/QM_MM_optimized.pdb
#   sbatch -J coscan --cpus-per-task=8  --mem=16G                run.slurm scan scan_znco start.pdb
#   sbatch -J neb1   --cpus-per-task=16 --mem=32G                run.slurm neb neb_step1
#   sbatch -J ts1    --cpus-per-task=8  --mem=16G                run.slurm ts ts_saddle1 saddlepoint.pdb
#
# Thread placement is not an argument: `mcube env` derives it from the stage and
# the config (see thread_policy in cli.py) and this script applies it, so a stage
# cannot be launched under the wrong policy. It also warns when the allocation is
# bigger than the stage can use.
set -euo pipefail
ulimit -s unlimited
export OMP_STACKSIZE=1G
export PYTHONUNBUFFERED=1

STAGE=${1:?stage: preprocess, md, qmmm, sp, scan, neb or ts}
SECTION=${2:?config.yaml section, or the stage name for preprocess/md}
START=${3:-}

# preprocess has no --section: it writes to the one place in the config, and the
# second argument only names its output directory. md takes one, and keeps its
# own `<section>_results` convention so `md md` still writes md_results.
case "$STAGE" in
  preprocess) SECTION_ARGS=();                     OUT=preprocess_results ;;
  md)         SECTION_ARGS=(--section "$SECTION"); OUT=${SECTION}_results ;;
  *)          SECTION_ARGS=(--section "$SECTION"); OUT=protdft/$SECTION ;;
esac

module load miniforge/26.1.0-0
source "$ROOTMINIFORGE/etc/profile.d/conda.sh"
conda activate ASH

ROOT=$PWD

# Thread placement, worked out from the stage and the config and applied here,
# before anything starts. Nothing to choose and nothing to remember.
eval "$(mcube env "$ROOT/config.yaml" "$STAGE" "${SECTION_ARGS[@]}")"

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
