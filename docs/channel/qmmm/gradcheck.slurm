#!/bin/bash
#SBATCH --partition=general --account=a_imb_ccbcd --cpus-per-task=32 --mem=96G --time=3:00:00
# 32 CPUs: the job configs set dft_cores 32 (exchange.py frames --cores)
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# Finite-difference check of one QM/MM job's gradient (gradcheck.py):
#   sbatch -J gradcheck --export=ALL,JOB=$S/channel/qmmm/water/LA_f0024/LA docs/channel/qmmm/gradcheck.slurm
# ARGS passes extra options, e.g. ARGS="--h 0.002" (no commas).

S=/scratch/user/uquadhik
module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH
export PYTHONPATH=$S/mcube/src PYTHONUNBUFFERED=1 OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
$CONDA_PREFIX/bin/python $S/mcube/docs/channel/qmmm/gradcheck.py --job "$JOB" $ARGS
echo "RC=$?"
