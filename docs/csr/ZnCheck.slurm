#!/bin/bash
#SBATCH --job-name=csr_zncheck
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=2-00:00:00
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --account=a_imb_ccbcd
#SBATCH --mail-type=END,FAIL
# No --mail-user here on purpose: this repo is public. See docs/csr/README.md.
#SBATCH --output=protdft/zncheck-%j.out
#SBATCH --error=protdft/zncheck-%j.err

# D53: is the six-coordinate zinc real, or the MM over-coordination artefact that
# Yang et al. (JCIM 2021, PMID 34748329) document? Relax, walk Zn-Asp92 OD2 from
# 2.05 to 3.05 A at xTB/MM, then recompute four of those geometries at DFT/MM.
# If the two methods put the minimum in different places, the mechanism in D52
# rests on an artefact.
set -euo pipefail
ulimit -s unlimited
export OMP_STACKSIZE=1G
export PYTHONUNBUFFERED=1
module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH

ROOT=$PWD
START=$ROOT/protdft/md2_ohH/QM_MM_lastframe.pdb
test -f "$START" || { echo "missing: $START" >&2; exit 1; }

# 1. Relax the MD end point at xTB/MM so the scan starts from a minimum,
#    not a thermal snapshot. 4 bound threads is the benchmarked optimum for xTB.
mkdir -p protdft/zncheck/relax
( export OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=cores
  cd protdft/zncheck/relax
  mcube qmmm ../../../config.yaml --pdb "$START" --section qmmm_znrelax )
echo "[DONE] relax"

# 2. Relaxed scan along Zn-Asp92 OD2.
mkdir -p protdft/zncheck/scan
( export OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=cores
  cd protdft/zncheck/scan
  mcube scan ../../../config.yaml --section scan_znasp \
      --pdb ../relax/QM_MM_optimized.pdb )
echo "[DONE] scan"

# 3. DFT/MM single points at four points along it: 2.05, 2.43, 2.80, 3.05 A.
#    All 16 cores here; PySCF runs in-process and must not be thread-bound.
for P in 00 03 06 08; do
  SRC=$ROOT/protdft/zncheck/scan/scan_$P.pdb
  if [ ! -f "$SRC" ]; then echo "no scan_$P.pdb, skipping"; continue; fi
  mkdir -p protdft/zncheck/sp_$P
  ( cd protdft/zncheck/sp_$P
    mcube qmmm ../../../config.yaml --pdb "$SRC" --section qmmm_sp_dft_$P ) \
      || echo "sp_$P FAILED"
  echo "[DONE] sp_$P"
done
echo "[ALL DONE] compare rel_kcal_mol in scan/scan_energies.csv (xTB) against the"
echo "sp_energy_disp values in the qmmm_sp_dft results for each point (DFT)."
