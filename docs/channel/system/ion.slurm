#!/bin/bash
#SBATCH --partition=gpu_cuda --qos=gpu --gres=gpu:1 --account=a_imb_ccbcd
#SBATCH --cpus-per-task=8 --mem=32G --time=12:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# One ion's system: preprocess (skip_solvate), short md with the shell restrained, then
# the umbrella pull. Run from the ion's run directory (config.yaml, channel_<ION>.pdb):
#
#   sbatch -J ion_LA $S/mcube/docs/channel/system/ion.slurm
#   sbatch -J pull_LA --export=ALL,STAGES=pull $S/mcube/docs/channel/system/ion.slurm
#
# Then the windows: umbrella_windows.slurm.

S=/scratch/user/uquadhik
module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH
export PYTHONPATH=$S/mcube/src PYTHONUNBUFFERED=1
PY=$CONDA_PREFIX/bin/python            # ASH's interpreter by path (see bulk/bulk.slurm)
# STAGES picks what runs (default all three), e.g. STAGES=pull after md has finished
STAGES=${STAGES:-"preprocess md pull"}
for STAGE in $STAGES; do
    case $STAGE in
        preprocess|md) eval "$($PY -m mcube env config.yaml $STAGE)"
                       $PY -m mcube $STAGE config.yaml || exit $? ;;
        pull)          eval "$($PY -m mcube env config.yaml umbrella)"
                       $PY -m mcube umbrella config.yaml pull || exit $? ;;
    esac
done
RC=0
echo "RC=$RC"
exit $RC
