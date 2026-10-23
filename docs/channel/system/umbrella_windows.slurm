#!/bin/bash
#SBATCH --partition=gpu_cuda --qos=gpu --gres=gpu:1 --account=a_imb_ccbcd
#SBATCH --cpus-per-task=8 --mem=32G --time=24:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# A block of MM umbrella windows, one after another on one GPU (the per-user GPU QOS allows
# 4 jobs, so windows are batched rather than one job each). From the ion's run directory:
#
#   sbatch -J w_LA_0 --export=ALL,FIRST=0,LAST=22 $S/mcube/docs/channel/system/umbrella_windows.slurm
#
# A window that already has results is skipped, so a block can be resubmitted after a timeout.

S=/scratch/user/uquadhik
module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH
export PYTHONPATH=$S/mcube/src PYTHONUNBUFFERED=1
PY=$CONDA_PREFIX/bin/python            # ASH's interpreter by path (see bulk/bulk.slurm)
eval "$($PY -m mcube env config.yaml umbrella)"
FAIL=0
for N in $(seq $FIRST $LAST); do
    W=$(printf '%03d' $N)
    if $PY -c "import yaml,sys; d=yaml.safe_load(open('config.yaml')) or {}; w=((d.get('umbrella') or {}).get('results') or {}).get('windows') or {}; sys.exit(0 if (w.get('$W') or {}).get('done') else 1)"; then
        echo "== window $N done already"; continue
    fi
    echo "== window $N"
    $PY -m mcube umbrella config.yaml window --window $N || FAIL=$((FAIL + 1))
done
echo "failed windows: $FAIL"
exit $FAIL
