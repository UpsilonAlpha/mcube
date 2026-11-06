#!/bin/bash
# Submit `asppull` - stage 1 of the Asp92 second-shell repair (D218/D220).
#
# Pulls Zn-CG from eq6's 3.72 A towards 6DQH's 5.15 A, halfway first, with eq6's
# own restraint list minus the Asp92-Zn pair and coordination_redetect OFF.
# `md`, not `scan`: a relaxed scan cannot converge here (D134/D220).
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

echo "=== sections.py asppull"
python "$R/sections.py" "$R/site.yaml" asppull

echo "=== what landed in config.yaml"
python - <<'PY'
import yaml
d = yaml.safe_load(open("config.yaml"))
s = d.get("md_asppull")
if not s:
    raise SystemExit("md_asppull is not in config.yaml - sections.py refused it")
for r in s["inputs"]["restraints"]:
    print("    restraint", r)
for k in ("pdb_file", "nvt_time", "nvt_timestep", "traj_frequency",
          "coordination_redetect", "platform"):
    print(f"    {k}: {s['inputs'].get(k)}")
PY

# 1 ns at 2 fs. eq6 did 10 ns in 43.7 min on one H100, so this is ~5 min; 1 h is
# ample. Same resources as every md job in this project (D136).
echo "=== sbatch"
sbatch -J asppull -p gpu_cuda --qos=gpu --gres=gpu:h100:1 \
  --cpus-per-task=4 --mem=32G -t 01:00:00 \
  "$R/run.slurm" md md_asppull

echo "=== queue"
squeue -u "$USER" -o "%.10i %.9P %.10j %.2t %.10M %.20R" || true
