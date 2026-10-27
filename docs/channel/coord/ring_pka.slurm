#!/bin/bash
#SBATCH --partition=general --account=a_imb_ccbcd --cpus-per-task=4 --mem=32G --time=4:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# PROPKA and the carboxylate environment of the ring Glu (ring_pka.py) for the apo, Ca and the four
# La/Lu residence trajectories. Writes $S/channel/pka/<label>/report.json.
#   sbatch -J ring_pka docs/channel/coord/ring_pka.slurm
# Needs the venv $S/envs/pka (ASH python with system site packages, plus `pip install propka`).

S=/scratch/user/uquadhik
SYS=$S/channel/system
PY=$S/envs/pka/bin/python
export PATH=$S/envs/pka/bin:$PATH
cd $S/channel || exit 1
R=$S/mcube/docs/channel/coord/ring_pka.py

$PY $R --top $SYS/apo/NVTtrajectory_firstframe.pdb --traj $SYS/apo/NVTtrajectory.dcd --stride 2 --label apo --out pka
$PY $R --top $SYS/CA/NVTtrajectory_firstframe.pdb --traj $SYS/CA/NVTtrajectory.dcd --stride 2 --label ca_nvt --out pka
for d in LA/umbrella_res_035:la_def LA/umbrella_res_chbe_035:la_chbe LU/umbrella_res_035:lu_def LU/umbrella_res_chbe_035:lu_chbe; do
    $PY $R --top $SYS/${d%%:*}/start.pdb --traj $SYS/${d%%:*}/window.dcd --label ${d##*:} --out pka
done
echo "RC=$?"
