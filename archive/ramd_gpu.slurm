#!/bin/bash
# preprocess -> md -> ramd on one NVIDIA GPU (Bunya gpu_cuda). Run from the system's
# directory, next to config.yaml:   sbatch /path/to/mcube/scripts/ramd_gpu.slurm
# STAGES="ramd" sbatch ... runs only some stages; --array=0-3 with RAMD_PER_REPLICA=1
# runs one replica per array task.
#SBATCH --job-name=mcube_ramd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=gpu_cuda
#SBATCH --qos=gpu
#SBATCH --gres=gpu:1
#SBATCH --account=a_imb_ccbcd
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PYTHONUNBUFFERED=1
module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH

STAGES=${STAGES:-"preprocess md ramd"}
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv

# Fail early if OpenMM cannot use the GPU (e.g. driver older than the env's CUDA)
python - <<'EOF'
import openmm
p = openmm.Platform.getPlatformByName("CUDA")
s = openmm.System(); s.addParticle(1.0)
ctx = openmm.Context(s, openmm.VerletIntegrator(0.001), p)
print("OpenMM", openmm.__version__, "CUDA device:", p.getPropertyValue(ctx, "DeviceName"))
EOF

for stage in $STAGES; do
    echo "[STEP] mcube $stage"
    case $stage in
        preprocess) mkdir -p preprocess_results && (cd preprocess_results && mcube preprocess ../config.yaml) ;;
        md)         mkdir -p md && (cd md && mcube md ../config.yaml) ;;
        ramd)
            mkdir -p ramd
            if [[ -n "${RAMD_PER_REPLICA:-}" && -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
                (cd ramd && mcube ramd ../config.yaml --replicas "$SLURM_ARRAY_TASK_ID")
            else
                (cd ramd && mcube ramd ../config.yaml)
            fi ;;
        *) echo "unknown stage $stage"; exit 1 ;;
    esac
done
echo "[DONE]"
