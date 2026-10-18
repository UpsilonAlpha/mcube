#!/bin/bash
#SBATCH --job-name=csr_cleanR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=2-00:00:00
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --account=a_imb_ccbcd
#SBATCH --mail-type=END,FAIL
# No --mail-user here on purpose: this repo is public, so an address committed
# here would be scraped. SLURM mails the submitting user's configured address by
# default. If your cluster has no default, either submit with an explicit
# address, which overrides the directives:
#     sbatch --mail-user=you@example.com THIS_SCRIPT.slurm
# or add a --mail-user line to your own working copy, which is not in git.
#SBATCH --output=protdft/cleanR-%j.out
#SBATCH --error=protdft/cleanR-%j.err

# D42: is the Lys226 -> OH- proton transfer still favourable in a clean
# four-coordinate sphere? Both states run sequentially from the same NVT end
# point and share lysH's active region, so only the proton position differs.
set -euo pipefail
ulimit -s unlimited
export OMP_STACKSIZE=1G
export PYTHONUNBUFFERED=1
module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH

START=$PWD/md_results/NVTtrajectory_lastframe.pdb
mkdir -p protdft/cleanR/lysH protdft/cleanR/ohH

(cd protdft/cleanR/lysH && mcube qmmm ../../../config.yaml --pdb "$START" --section qmmm_cleanR_lysH)
echo "[DONE] lysH"
(cd protdft/cleanR/ohH && mcube qmmm ../../../config.yaml --pdb "$START" --section qmmm_cleanR_ohH \
    --actatoms ../lysH/active_atoms)
echo "[DONE] ohH"
