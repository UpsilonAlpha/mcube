#!/bin/bash
# Submit `relax10` - a 10 ns continuation of the repaired structure (D225).
#
# NOT 100 ns. `convergence.py` showed eq6's 10 ns is converged (RMSD plateaus by
# ~3 ns, RMSF within 3% of final by 5 ns, T 300.31 +/- 1.45) while asppull3's
# 2 ns is not (RMSD still climbing, box volume drift 8x eq6's). So this matches
# the length eq6 demonstrated is sufficient for this system, and no more.
#
# `gentle: false` so it continues rather than re-annealing from 1 K as every
# pull stage did; NPT left ON because the box is still relaxing after the pull.
#
# In the repo rather than typed, because the code-server terminal truncates long
# single-line commands silently (bunya-working-rules).
set -euo pipefail

R=/scratch/user/uquadhik/mcube/docs/csr
W=/scratch/user/uquadhik/CSR/CSR_mich

module load miniforge/26.1.0-0
# shellcheck disable=SC1091
source "$ROOTMINIFORGE/etc/profile.d/conda.sh"
conda activate ASH

cd "$W"
echo "=== cwd: $PWD"
echo "=== date: $(date)"
echo "=== quota BEFORE (rule 8: an over-quota write fails silently)"
rquota

echo "=== sections.py relax10"
python "$R/sections.py" "$R/site.yaml" relax10

echo "=== what landed in config.yaml"
python - <<'PY'
import yaml
d = yaml.safe_load(open("config.yaml"))
s = d.get("md_relax10")
if not s:
    raise SystemExit("md_relax10 is not in config.yaml - sections.py refused it")
i = s["inputs"]
for k in ("pdb_file", "nvt_time", "nvt_timestep", "gentle", "skip_npt",
          "coordination_redetect", "traj_frequency", "platform"):
    print(f"    {k}: {i.get(k)}")
for r in i.get("restraints", []):
    print("    restraint", r)
PY

# 10 ns at 2 fs on one H100; eq6 did exactly this in 43.7 min, so 3 h is ample.
echo "=== sbatch"
sbatch -J relax10 -p gpu_cuda --qos=gpu --gres=gpu:h100:1 \
  --cpus-per-task=4 --mem=32G -t 03:00:00 \
  "$R/run.slurm" md md_relax10

echo "=== queue"
squeue -u "$USER" -o "%.10i %.9P %.10j %.2t %.10M %.20R" || true
