#!/bin/bash
#SBATCH --partition=gpu_cuda --qos=gpu --gres=gpu:1 --account=a_imb_ccbcd
#SBATCH --cpus-per-task=8 --mem=32G --time=12:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# Build and equilibrate the apo channel in POPE: preprocess (CPU part) then md (GPU).
# Run from the run directory holding config.yaml and channel_cropped_oriented.pdb:
#
#   sbatch -J apo $S/mcube/docs/channel/system/apo.slurm

S=/scratch/user/uquadhik
module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH
export PYTHONPATH=$S/mcube/src PYTHONUNBUFFERED=1
PY=$CONDA_PREFIX/bin/python            # ASH's interpreter by path (see bulk/bulk.slurm)
$PY -c "import openmm; print([openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])"
eval "$($PY -m mcube env config.yaml preprocess)"
$PY -m mcube preprocess config.yaml || exit $?
eval "$($PY -m mcube env config.yaml md)"
$PY -m mcube md config.yaml
RC=$?
echo "RC=$RC"
exit $RC
