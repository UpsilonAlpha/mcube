#!/bin/bash
# Submit `eq6` - eq5 plus the fifth-donor and salt-bridge restraints (D210, D212).
#
# In the repo rather than typed, because the code-server terminal truncates long
# single-line commands silently and these paths are well over 100 characters
# (bunya-working-rules). Writes everything it does to stdout; the caller
# redirects that to a file under /scratch and reads it back over
# vscode-remote-resource, because the terminal is a canvas and cannot be read.
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

# sections.py refuses to write unless every structure exists and every selector
# resolves in it. That is the gate on the three restraint distances corrected in
# D212 and on `A:226:NZ`, which no section in this project has ever named.
echo "=== sections.py eq6"
python "$R/sections.py" "$R/site.yaml" eq6

echo "=== what landed in config.yaml"
python - <<'PY'
import yaml
d = yaml.safe_load(open("config.yaml"))
s = d.get("md_eq6")
if not s:
    raise SystemExit("md_eq6 is not in config.yaml - sections.py refused it")
for r in s["inputs"]["restraints"]:
    print("    restraint", r)
for k in ("pdb_file", "nvt_time", "nvt_timestep", "coordination_redetect", "platform"):
    print(f"    {k}: {s['inputs'].get(k)}")
PY

# md on a GPU: run.slurm's own defaults are the CPU general partition, so the
# resources are given here, as every md job in this project has been (D136).
# 10 ns at a 2 fs timestep; eq5's own length. 6 h is ample and well inside the
# QOS limit.
echo "=== sbatch"
sbatch -J eq6 -p gpu_cuda --qos=gpu --gres=gpu:h100:1 \
  --cpus-per-task=4 --mem=32G -t 06:00:00 \
  "$R/run.slurm" md md_eq6

echo "=== queue"
squeue -u "$USER" -o "%.10i %.9P %.10j %.2t %.10M %.20R" || true
