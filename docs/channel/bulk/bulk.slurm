#!/bin/bash
#SBATCH --partition=general --account=a_imb_ccbcd --cpus-per-task=8 --mem=16G --time=06:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# Bulk-water check of one ion model, in three steps, all CPU:
#   1. MM: bulk_ion.py run (NPT TIP3P box, mcube's metal parameters, shell analysis)
#   2. xTB/MM MD of the same box through mcube's own stages: preprocess (skip_solvate)
#      detects the shell and makes ion + shell the QM region; qmmm runs QM/MM MD only
#   3. the same shell analysis on the QM/MM trajectory
#
#   sbatch -J bulk_LA_iod docs/channel/bulk/bulk.slurm LA iod
#
# Writes $RUN/report.json (MM) and $RUN/qmmm_report/report.json (xTB/MM),
# RUN=/scratch/user/uquadhik/channel/bulk/<ION>_<SET>. NS and QMPS set the lengths; RUN can be
# overridden, e.g. for a smoke test: RUN=$S/channel/bulk/smoke NS=0.01 QMPS=0.05 bash bulk.slurm LA iod

ION=$1
SET=$2
NS=${NS:-4}
QMPS=${QMPS:-10}
EVERY=20                      # QM/MM steps between saved frames (1 fs steps)
DT=0.02                       # ps between saved QM/MM frames = EVERY * 1 fs
S=/scratch/user/uquadhik
B=$S/mcube/docs/channel/bulk/bulk_ion.py
RUN=${RUN:-$S/channel/bulk/${ION}_${SET}}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}   # interactive smoke tests: 2 threads

module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH
export PYTHONPATH=$S/mcube/src PYTHONUNBUFFERED=1

echo "== MM $ION $SET, $NS ns"
python $B run --ion $ION --set $SET --ns $NS --out $RUN || exit $?

echo "== xTB/MM, $QMPS ps"
python $B qmmm-config --run $RUN --ion $ION --ps $QMPS --every $EVERY --cores 4 || exit $?
cd $RUN/qmmm || exit 1
eval "$(python -m mcube env config.yaml preprocess)"
python -m mcube preprocess config.yaml || exit $?
eval "$(python -m mcube env config.yaml qmmm)"
python -m mcube qmmm config.yaml
RC=$?

echo "== shell analysis of the QM/MM trajectory"
python $B analyse --top $RUN/system.pdb --traj QM_MM.dcd --ion $ION --dt $DT --skip 2 \
    --out $RUN/qmmm_report
echo "RC=$RC"
exit $RC
