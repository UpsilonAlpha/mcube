#!/bin/bash
#SBATCH --partition=general --account=a_imb_ccbcd --cpus-per-task=32 --mem=96G --time=08:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# g-xTB vs GFN2-xTB vs r2SCAN/def2-SVP on La3+/Lu3+ aqua clusters cut from the failing
# xTB/MM trajectories (cluster_test.py). Smoke test, xTB only, one frame:
#
#   sbatch -c 8 --mem=16G -t 0:30:00 -J gx_smoke --export=ALL,OPTS="--frames 1 --skip-dft --skip-opt" cluster.slurm
#   sbatch -J gx_test cluster.slurm

S=/scratch/user/uquadhik
module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH
export PYTHONUNBUFFERED=1 PYSCF_MAX_MEMORY=80000
PY=$CONDA_PREFIX/bin/python            # ASH's interpreter by path (see bulk/bulk.slurm)
OUT=${OUT:-$S/channel/gxtb/${SLURM_JOB_NAME}}
$PY $S/mcube/docs/channel/gxtb/cluster_test.py \
    --run LA=$S/channel/bulk/LA_1264 --run LU=$S/channel/bulk/LU_1264 \
    --gxtb $S/opt/gxtb/xtb-6.7.1/bin/xtb --gfn2 $CONDA_PREFIX/bin/xtb \
    --basis-file $S/channel/dft/def2-svp.nw --cores $SLURM_CPUS_PER_TASK --out $OUT $OPTS
RC=$?
echo "RC=$RC"
exit $RC
