#!/bin/bash
#SBATCH --partition=general --account=a_imb_ccbcd --cpus-per-task=32 --mem=96G --time=24:00:00
#SBATCH --output=%x-%A_%a.out --error=%x-%A_%a.err
#
# One QM/MM vertical-swap job (docs/channel/qmmm/exchange.py): one frame with one ion, through
# mcube's own stages on the job directory's config.yaml:
#   preprocess (skip_solvate) -> qmmm (no QM/MM MD; DFT/MM optimisation of the 8 A region) -> sp
#
# As an array over a job list that `exchange.py frames` wrote (submit.sh does this):
#   sbatch -J qx_site --array=1-30%10 --export=ALL,LIST=$S/channel/qmmm/jobs_site.txt exchange.slurm
# or one job directory:
#   sbatch -J qx_one --export=ALL,JOB=$S/channel/qmmm/site/LA_f0399/LU exchange.slurm
#
# Resubmitting is safe: a stage whose results are already in config.yaml is skipped (a qmmm
# that was killed mid-optimisation starts that optimisation again). Each job's own log is
# <job>/job.log. --cpus-per-task must match `cores`/`dft_cores`/`sp_cores` in the configs
# (exchange.py frames --cores, default 32); never put commas inside an --export value.

S=/scratch/user/uquadhik
if [ -z "$JOB" ]; then
    [ -n "$LIST" ] && [ -n "$SLURM_ARRAY_TASK_ID" ] || { echo "set JOB, or LIST with --array"; exit 1; }
    JOB=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
fi
[ -f "$JOB/config.yaml" ] || { echo "no config.yaml in '$JOB'"; exit 1; }

module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH
export PYTHONPATH=$S/mcube/src PYTHONUNBUFFERED=1
PY=$CONDA_PREFIX/bin/python            # ASH's interpreter by path (see bulk/bulk.slurm)
$PY -c "import openmm, ash" || { echo "ASH env not usable: $PY"; exit 1; }
python() { $PY "$@"; }

cd "$JOB" || exit 1
echo "== $(date) job $SLURM_JOB_ID ($SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID) on $(hostname): $JOB" | tee -a job.log
exec >> job.log 2>&1

# done <section> <key>: is <section>.results.<key> set in config.yaml?
done_() {
    python -c "import yaml,sys; d=yaml.safe_load(open('config.yaml')) or {}; \
sys.exit(0 if ((d.get('$1') or {}).get('results') or {}).get('$2') is not None else 1)"
}

if done_ sp sp_energy; then echo "== finished already"; exit 0; fi
if ! done_ preprocess processed_pdb; then
    echo "== preprocess $(date)"
    eval "$(python -m mcube env config.yaml preprocess)"
    python -m mcube preprocess config.yaml || exit $?
fi
if ! done_ qmmm optimized_pdb; then
    echo "== qmmm (optimisation only) $(date)"
    eval "$(python -m mcube env config.yaml qmmm)"
    python -m mcube qmmm config.yaml || exit $?
fi
echo "== sp $(date)"
eval "$(python -m mcube env config.yaml sp)"
python -m mcube sp config.yaml
RC=$?
echo "== $(date) RC=$RC"
exit $RC
