#!/bin/bash
#SBATCH --partition=general --account=a_imb_ccbcd --cpus-per-task=32 --mem=96G --time=03:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# Time DFT/MM gradients (r2SCAN/def2-SVP, density fitting) for an Ln3+ ion and its first shell
# in a bulk box, against one xTB/MM gradient. Run at two core counts to see the scaling:
#
#   sbatch -c 16 --mem=64G -J dft16 docs/channel/dft/bench.slurm LA 1264
#   sbatch -c 32 --mem=96G -J dft32 docs/channel/dft/bench.slurm LA 1264
#
# Extra bench_dft.py options go in BENCH_OPTS (e.g. "--basis-file lu.nw --ecp-file lu.nw").

ION=$1
SET=$2
S=/scratch/user/uquadhik
module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH
export PYTHONPATH=$S/mcube/src PYTHONUNBUFFERED=1
PY=$CONDA_PREFIX/bin/python            # ASH's interpreter by path (see bulk.slurm)
python() { $PY "$@"; }
# PySCF sets its own thread count from --cores; xTB is bound like the qmmm stage
export OMP_NUM_THREADS=1
python $S/mcube/docs/channel/dft/bench_dft.py --pdb $S/channel/bulk/${ION}_${SET}/system.pdb \
    --ion $ION --set $SET --cores $SLURM_CPUS_PER_TASK \
    --out $S/channel/dft/bench_${ION}_${SLURM_CPUS_PER_TASK} $BENCH_OPTS
RC=$?
echo "RC=$RC"
exit $RC
