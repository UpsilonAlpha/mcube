#!/bin/bash
#SBATCH --account=a_imb_ccbcd
#SBATCH --nodes=1 --ntasks=1
#SBATCH --output=protdft/%x-%j.out --error=protdft/%x-%j.err
#
# A block of umbrella windows, one after another, from the run directory.
#
#   sbatch -J u_asp_0 --cpus-per-task=8 --mem=64G -t 3-00:00:00 \
#          --export=ALL,SECTION=umbrella_asp,FIRST=0,LAST=5 \
#          $S/mcube/docs/csr/umbrella_windows.slurm
#
# This is `docs/channel/system/umbrella_windows.slurm` with the resource
# directives REMOVED, so the caller supplies them on the sbatch line exactly as
# `run.slurm` does. The channel version hardcodes `--partition=gpu_cuda
# --gres=gpu:1`, which is right for its MM windows along a channel's z and wrong
# for these: a QM/MM window spends its time in xTB, and a GPU does not help it
# (user's standing rule). The two scripts should eventually be one; until then
# this is the only difference, and it is deliberate.
#
# Resources that match the measured rate (D157): 2.57 ps/h on 8 cores, so a
# 105 ps window (45 ps equilibration + 60 ps data) takes **40.9 h**. One window
# per job fits a 3-day allocation with room; two do not. Give each block at most
# one window unless the allocation is longer, and run blocks in parallel -
# 26 windows is about two days of wall clock, not two days times 26.
#
# `umbrella` takes a SUBCOMMAND (`pull` / `window --window N` / `wham`), which
# `run.slurm` cannot pass, which is why umbrella needs its own launcher at all.
# A window that already has results is skipped, so a block can be resubmitted
# after a timeout and will pick up where it stopped.
set -euo pipefail
ulimit -s unlimited
export OMP_STACKSIZE=1G

S=/scratch/user/uquadhik
module load miniforge/26.1.0-0
# shellcheck disable=SC1091
source "$ROOTMINIFORGE/etc/profile.d/conda.sh"
conda activate ASH
export PYTHONPATH=$S/mcube/src PYTHONUNBUFFERED=1
PY=$CONDA_PREFIX/bin/python

SECTION=${SECTION:-umbrella}
: "${FIRST:?give FIRST and LAST via --export}"
: "${LAST:?give FIRST and LAST via --export}"

# Thread placement from the stage and the config, applied before anything runs,
# the same way run.slurm does it.
eval "$($PY -m mcube env config.yaml umbrella)"
mkdir -p protdft

FAIL=0
for N in $(seq "$FIRST" "$LAST"); do
    W=$(printf '%03d' "$N")
    if $PY - "$SECTION" "$W" <<'PY'
import sys, yaml
section, key = sys.argv[1], sys.argv[2]
d = yaml.safe_load(open("config.yaml")) or {}
w = ((d.get(section) or {}).get("results") or {}).get("windows") or {}
sys.exit(0 if (w.get(key) or {}).get("done") else 1)
PY
    then
        echo "== $SECTION window $N: done already, skipping"
        continue
    fi
    echo "== $SECTION window $N  $(date)"
    $PY -m mcube umbrella config.yaml window --window "$N" --section "$SECTION" \
        || FAIL=$((FAIL + 1))
done
echo "failed windows: $FAIL"
exit "$FAIL"
