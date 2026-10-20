#!/bin/bash
#SBATCH --job-name=csr_coscan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --account=a_imb_ccbcd
#SBATCH --mail-type=END,FAIL
# No --mail-user here on purpose: this repo is public. See docs/csr/README.md.
#SBATCH --output=protdft/coscan-%j.out
#SBATCH --error=protdft/coscan-%j.err

# D58: relaxed xTB/MM scan pulling the substrate carbonyl off the zinc
# (Zn-O2x 2.15 -> 4.00 A), to see whether the attack geometry improves as the
# site opens. Same start structure as the Zn-Asp92 scan so the two are
# comparable. 4 bound threads is the benchmarked optimum for xTB.
set -euo pipefail
ulimit -s unlimited
export OMP_STACKSIZE=1G
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=4
export OMP_PROC_BIND=close
export OMP_PLACES=cores

module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH

ROOT=$PWD
START=$ROOT/protdft/zncheck/relax/QM_MM_optimized.pdb
test -f "$START" || { echo "missing: $START" >&2; exit 1; }
mkdir -p protdft/coscan
cd protdft/coscan
mcube scan ../../config.yaml --section scan_znco --pdb "$START"
echo "[DONE] coscan"

# Attack angles are not available through `monitor`, so derive them from the
# structures. zn_sphere.py prints the Zn sphere and both O-C4x-O2x angles.
python "$ROOT/../mcube/docs/csr/zn_sphere.py" scan_*.pdb || true
