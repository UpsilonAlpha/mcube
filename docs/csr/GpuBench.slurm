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

# cupy's cuda12x wheels need the CUDA 12 runtime, and Bunya has no CUDA 12 module
# (only cuda11.7), so the libraries come from the nvidia-*-cu12 pip wheels. Put their
# lib directories on the loader path, or cupy fails with
# "ImportError: libcusolver.so.11: cannot open shared object file" (job 30054626).
NVLIB=$(python -c "import site,glob,os; p=site.getsitepackages()[0]; print(':'.join(sorted(glob.glob(os.path.join(p,'nvidia','*','lib')))))")
export LD_LIBRARY_PATH="$NVLIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
python -c "import cupy; cupy.zeros(3).sum(); print('cupy OK', cupy.__version__, cupy.cuda.runtime.getDeviceCount(), 'device(s)')"

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
