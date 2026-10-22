#!/bin/bash
#SBATCH --partition=general --account=a_imb_ccbcd --cpus-per-task=8 --mem=16G --time=12:00:00
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
# Writes $RUN/report.json (MM) and $RUN/$QMDIR/../qmmm_report.../report.json (xTB/MM),
# RUN=/scratch/user/uquadhik/channel/bulk/<ION>_<SET>. Settings by --export:
#   NS, EQUIL, QMPS   MM length (ns), MM equilibration (ps), QM/MM length (ps)
#   SKIP_MM=1         reuse $RUN's finished MM box and report; run only the QM/MM leg
#   QMDIR             subdirectory for the QM/MM leg (default qmmm); its report goes to
#                     $RUN/${QMDIR}_report
#   XTB_ETEMP         xTB electronic temperature (K); SHELL2 a QM cutoff (A) that takes in
#                     the second shell; PC_HARDNESS the xTB hardness of the MM point charges
# MM on these 4.1k-atom boxes ran ~15-20 ns/day on 8 CPUs, so 4 ns needs ~5-6 h; the limit is 12 h.
# Smoke test:
#   sbatch --export=ALL,RUN=$S/channel/bulk/smoke,NS=0.01,EQUIL=10,QMPS=0.2 -t 1:00:00 bulk.slurm LA iod

ION=$1
SET=$2
NS=${NS:-4}
QMPS=${QMPS:-10}
EQUIL=${EQUIL:-200}           # ps of MM NPT before sampling
QMDIR=${QMDIR:-qmmm}
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

if [ "${SKIP_MM:-0}" = 1 ]; then
    echo "== MM skipped: reusing $RUN"
    [ -f $RUN/report.json ] && [ -f $RUN/system.pdb ] || { echo "no MM box in $RUN"; exit 1; }
else
    echo "== MM $ION $SET, $NS ns"
    python $B run --ion $ION --set $SET --ns $NS --equil $EQUIL --out $RUN || exit $?
fi

echo "== xTB/MM, $QMPS ps, in $QMDIR (etemp ${XTB_ETEMP:-default}, shell2 ${SHELL2:-no})"
OPTS=""
[ -n "$XTB_ETEMP" ] && OPTS="$OPTS --etemp $XTB_ETEMP"
[ -n "$PC_HARDNESS" ] && OPTS="$OPTS --pc-hardness $PC_HARDNESS"
[ -n "$SHELL2" ] && OPTS="$OPTS --shell2 $SHELL2"
python $B qmmm-config --run $RUN --ion $ION --ps $QMPS --every $EVERY --cores 4 --dir $QMDIR \
    $OPTS || exit $?
cd $RUN/$QMDIR || exit 1
eval "$(python -m mcube env config.yaml preprocess)"
python -m mcube preprocess config.yaml || exit $?
eval "$(python -m mcube env config.yaml qmmm)"
python -m mcube qmmm config.yaml
RC=$?

echo "== shell analysis of the QM/MM trajectory"
python $B analyse --top $RUN/system.pdb --traj QM_MM.dcd --ion $ION --dt $DT --skip 2 \
    --out $RUN/${QMDIR}_report
echo "RC=$RC"
exit $RC
