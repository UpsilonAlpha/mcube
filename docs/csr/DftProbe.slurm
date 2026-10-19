#!/bin/bash
#SBATCH --job-name=csr_dftprobe
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --account=a_imb_ccbcd
#SBATCH --mail-type=END,FAIL
# No --mail-user here on purpose: this repo is public. See docs/csr/README.md.
#SBATCH --output=protdft/dftprobe-%j.out
#SBATCH --error=protdft/dftprobe-%j.err

# D55: how long does one r2SCAN-D4/def2-TZVP QM/MM gradient take on the real
# 119-atom QM region? Three optimisation steps, then stop. Read the per-step
# wall time before committing to the two full optimisations.
#
# No thread binding here: PySCF runs in-process, so binding would pin its BLAS
# threads to one core (docs/running-on-bunya.md).
set -euo pipefail
ulimit -s unlimited
export OMP_STACKSIZE=1G
export PYTHONUNBUFFERED=1
module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH

START=$PWD/md_results/NVTtrajectory_lastframe.pdb
test -f "$START" || { echo "missing: $START" >&2; exit 1; }
mkdir -p protdft/dftopt/probe
cd protdft/dftopt/probe
date
mcube qmmm ../../../config.yaml --pdb "$START" --section dftopt_probe
date
echo "[DONE] probe - per-step times are in optimization_energies.log"
