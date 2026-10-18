#!/bin/bash
#SBATCH --job-name=csr_gpubench
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --partition=gpu_cuda
#SBATCH --qos=gpu
#SBATCH --account=a_imb_ccbcd
#SBATCH --output=protdft/gpubench-%j.out
#SBATCH --error=protdft/gpubench-%j.err

# D45: DFT/MM gradient timing on a GPU. A few optimisation steps per section; the
# per-step wall time in each optimization_energies.log is the gradient cost.
set -euo pipefail
ulimit -s unlimited
export OMP_STACKSIZE=1G
export PYTHONUNBUFFERED=1
module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH

nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
START=$PWD/md_results/NVTtrajectory_lastframe.pdb

for SEC in qmmm_gpu_wb97x qmmm_gpu_r2scan qmmm_gpu_r2scan_pc; do
  mkdir -p protdft/gpubench/$SEC
  echo "=== $SEC start $(date +%T)"
  # keep going if one level fails, so the others still give numbers
  (cd protdft/gpubench/$SEC && mcube qmmm ../../../config.yaml --pdb "$START" --section $SEC) \
      || echo "=== $SEC FAILED"
  echo "=== $SEC end $(date +%T)"
done
