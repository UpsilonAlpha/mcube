#!/bin/bash
#SBATCH --partition=general --account=a_imb_ccbcd --cpus-per-task=32 --mem=96G --time=24:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# Larger-basis + counterpoise single points (cp_tzvp.py) on the geometries a
# carboxylate_exchange.py run saved in <run>/geom/. Chain it after the geometry run:
#
#   sbatch -J coord_cp --dependency=afterok:<geometry job> docs/channel/coord/coord_cp.slurm

S=/scratch/user/uquadhik
module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH
export PYTHONPATH=$S/mcube/src PYTHONUNBUFFERED=1 PYSCF_MAX_MEMORY=80000 OMP_NUM_THREADS=1
PY=$CONDA_PREFIX/bin/python            # ASH's interpreter by path (see bulk/bulk.slurm)
RUN=${RUN:-$S/channel/coord/coord_check}
BASIS=${BASIS:-$S/channel/dft/def2-tzvp_lcecp1.nw}
$PY $S/mcube/docs/channel/coord/cp_tzvp.py --run $RUN --basis-file $BASIS \
    --cores $SLURM_CPUS_PER_TASK --out $S/channel/coord/${SLURM_JOB_NAME}
RC=$?
echo "RC=$RC"
exit $RC
