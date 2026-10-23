#!/bin/bash
#SBATCH --partition=general --account=a_imb_ccbcd --cpus-per-task=8 --mem=16G --time=24:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# The Lu3+ 12-6-4 box, Amber-built, run without mcube in OpenMM or in AmberTools' sander
# (amber_crosscheck.py). sander here is serial, so its job uses 1 core:
#
#   sbatch -J xc_LU_openmm amber_crosscheck.slurm LU openmm
#   sbatch -J xc_LU_sander -c 1 --mem=4G amber_crosscheck.slurm LU sander

ION=$1
ENGINE=$2
S=/scratch/user/uquadhik
module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH
export PYTHONPATH=$S/mcube/src PYTHONUNBUFFERED=1 OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
PY=$CONDA_PREFIX/bin/python            # ASH's interpreter by path (see bulk.slurm)
export PATH=$CONDA_PREFIX/bin:$PATH    # tleap, parmed, sander from the same env
$PY $S/mcube/docs/channel/bulk/amber_crosscheck.py --pdb $S/channel/bulk/${ION}_1264_lsm/system.pdb \
    --ion $ION --engine $ENGINE --out $S/channel/bulk/xc_${ION}_${ENGINE} $OPTS
RC=$?
echo "RC=$RC"
exit $RC
