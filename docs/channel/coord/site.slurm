#!/bin/bash
#SBATCH --partition=general --account=a_imb_ccbcd --cpus-per-task=32 --mem=96G --time=24:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# Ion exchange free energies at the constriction site (site_exchange.py), one site cluster per
# job, or the acetate calibration. Sites are built beforehand with `site_exchange.py build`.
#
#   sbatch -J site_la_def_ring --export=ALL,ACTION=run,SITE=$S/channel/site/la_def_ring docs/channel/coord/site.slurm
#   sbatch -J site_calib --export=ALL,ACTION=calibrate docs/channel/coord/site.slurm
#
# Resubmitting continues: finished optimisations and single points are read from report.json.
# SMD=1 runs on the SMD-enabled PySCF venv ($S/envs/smd, build_pyscf_smd.slurm); pass
# --solvation smd in OPTS, and for calibrate CALIB_OUT=<fresh dir>.

S=/scratch/user/uquadhik
module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH
export PYTHONPATH=$S/mcube/src PYTHONUNBUFFERED=1 PYSCF_MAX_MEMORY=80000 OMP_NUM_THREADS=1
PY=$CONDA_PREFIX/bin/python            # ASH's interpreter by path (see bulk/bulk.slurm)
if [ "${SMD:-0}" = 1 ]; then
    module load foss/2023a             # the runtime libraries the SMD build linked against
    PY=$S/envs/smd/bin/python
fi
X=$S/mcube/docs/channel/coord/site_exchange.py
D=$S/channel/dft
COMMON="--aqua $S/channel/coord/coord_check/geom --basis-tzvp $D/def2-tzvp_lcecp1.nw --cores $SLURM_CPUS_PER_TASK"
if [ "$ACTION" = calibrate ]; then
    $PY $X calibrate $COMMON --out ${CALIB_OUT:-$S/channel/site/calibrate} $OPTS
else
    $PY $X run $COMMON --site $SITE --basis-svp $D/def2-svp_lcecp1.nw $OPTS
fi
RC=$?
echo "RC=$RC"
exit $RC
