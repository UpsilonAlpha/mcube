#!/bin/bash
#SBATCH --job-name=csr_md3pose
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --account=a_imb_ccbcd
#SBATCH --mail-type=END,FAIL
# No --mail-user here on purpose: this repo is public. See docs/csr/README.md.
#SBATCH --output=protdft/md3pose-%j.out
#SBATCH --error=protdft/md3pose-%j.err

# D59 pose-robustness test: 25 ps of unconstrained QM/MM MD starting from the
# far end of the carbonyl scan (Zn-O2x 4.00 A), to find out whether the
# substrate pose is a real preference or an inheritance from the rigid MM
# constraint that pinned the carbonyl at 1.94 A through equilibration.
#
# Nothing is constrained here: qmmm never reads preprocess's metal constraints,
# so where the carbonyl ends up is the answer.
set -euo pipefail
ulimit -s unlimited
export OMP_STACKSIZE=1G
export PYTHONUNBUFFERED=1

# 8 CPUs on an SMT node is 4 physical cores; bind xTB's 4 threads to them.
export OMP_NUM_THREADS=4
export OMP_PROC_BIND=close
export OMP_PLACES=cores

module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH

ROOT=$PWD
START=$ROOT/protdft/coscan/scan_08.pdb
test -f "$START" || { echo "missing: $START (did the carbonyl scan finish?)" >&2; exit 1; }
mkdir -p protdft/md3_pose
cd protdft/md3_pose
mcube qmmm ../../config.yaml --pdb "$START" --section qmmm_md3_pose
echo "[DONE] md3pose"

# Where did the site end up? Compare against the start.
python "$ROOT/../mcube/docs/csr/zn_sphere.py" "$START" QM_MM_lastframe.pdb || true
