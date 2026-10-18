#!/bin/bash
#SBATCH --job-name=csr_md2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --account=a_imb_ccbcd
#SBATCH --mail-type=END,FAIL
# No --mail-user here on purpose: this repo is public, so an address committed
# here would be scraped. See docs/csr/README.md.
#SBATCH --output=protdft/md2-%j.out
#SBATCH --error=protdft/md2-%j.err

# D50: 25 ps QM/MM MD from the corrected resting-state reactant (D48), to ask
# whether a free Asp92 recruits a water into the donating position - the thing
# D49 found missing everywhere in the original trajectory.
set -euo pipefail
ulimit -s unlimited
export OMP_STACKSIZE=1G
export PYTHONUNBUFFERED=1

# 8 CPUs on an SMT node is 4 physical cores; bind xTB's 4 threads to them.
# mcube binds xTB only, through a wrapper, so PySCF and BLAS stay unpinned
# (see docs/running-on-bunya.md).
export OMP_NUM_THREADS=4
export OMP_PROC_BIND=close
export OMP_PLACES=cores

module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH

START=$PWD/protdft/cleanR/ohH/QM_MM_optimized.pdb
test -f "$START" || { echo "missing start structure: $START" >&2; exit 1; }
mkdir -p protdft/md2_ohH
cd protdft/md2_ohH
mcube qmmm ../../config.yaml --pdb "$START" --section qmmm_md2_ohH
echo "[DONE] md2; now run joint_nac.py on QM_MM.dcd with QM_MM_lastframe.pdb"
