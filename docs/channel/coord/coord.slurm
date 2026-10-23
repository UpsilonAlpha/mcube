#!/bin/bash
#SBATCH --partition=general --account=a_imb_ccbcd --cpus-per-task=32 --mem=96G --time=24:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# 12-6-4 vs r2SCAN for water -> carboxylate exchange around Ca2+, La3+, Lu3+
# (carboxylate_exchange.py). Extra options in OPTS, e.g. OPTS="--set iod".
#
#   sbatch -J coord_check docs/channel/coord/coord.slurm

S=/scratch/user/uquadhik
module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH
export PYTHONPATH=$S/mcube/src PYTHONUNBUFFERED=1 PYSCF_MAX_MEMORY=80000 OMP_NUM_THREADS=1
PY=$CONDA_PREFIX/bin/python            # ASH's interpreter by path (see bulk/bulk.slurm)
$PY $S/mcube/docs/channel/coord/carboxylate_exchange.py --basis-file $S/channel/dft/def2-svp.nw \
    --cores $SLURM_CPUS_PER_TASK --out $S/channel/coord/${SLURM_JOB_NAME} $OPTS
RC=$?
echo "RC=$RC"
exit $RC
