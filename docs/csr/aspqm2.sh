#!/bin/bash
# Submit `aspqm2` - one 50 ps xTB QM/MM chunk from the repaired structure (D224).
#
# `time` is the TOTAL ps wanted and `resume: true` continues rather than
# overwrites, so the 300 ps the pipeline calls for is reached by raising
# `aspqm2.inputs.time` in site.yaml - 50, 100, 150, 200, 250, 300 - and running
# this script again. ~19.5 h per chunk at the measured 2.57 ps/h.
#
# CPU, not GPU: xTB is the bottleneck in QM/MM and a GPU does not help.
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

echo "=== sections.py aspqm2"
python "$R/sections.py" "$R/site.yaml" aspqm2

echo "=== what landed in config.yaml"
python - <<'PY'
import yaml
d = yaml.safe_load(open("config.yaml"))
s = d.get("qmmm_aspqm2")
if not s:
    raise SystemExit("qmmm_aspqm2 is not in config.yaml - sections.py refused it")
i = s["inputs"]
for k in ("pdb_file", "charge", "time", "resume",
          "skip_equilibrate", "skip_optimize", "rigid_water"):
    print(f"    {k}: {i.get(k)}")
print("    qm_residues:", i.get("qm_residues"))
# Nothing should be restrained in this leg. Say so out loud rather than assume.
for k in ("restraints", "opt_constraints", "constraints"):
    if i.get(k):
        print(f"    [WARNING] {k} is set: {i[k]}")
# The whole point of this section is that it starts from the END of the 10 ns
# NVT, not from the unrelaxed 2 ns frame that job 32384737 was cancelled for
# using. Refuse rather than print, and refuse if that frame is not on disk yet
# (relax10 still running, or failed) - D228.
import os
want = "md_relax10_results/NVTtrajectory_lastframe.pdb"
if i.get("pdb_file") != want:
    raise SystemExit(f"pdb_file is {i.get('pdb_file')!r}, expected {want!r}")
if not os.path.isfile(want):
    raise SystemExit(f"{want} does not exist - relax10 has not finished")
if i.get("charge") != 2:
    raise SystemExit(f"charge is {i.get('charge')}, expected 2")
print("    [OK] starts from relax10's last frame, charge 2, nothing restrained")
PY

echo "=== sbatch"
sbatch -J aspqm2 -p general \
  --cpus-per-task=8 --mem=16G -t 24:00:00 \
  "$R/run.slurm" qmmm qmmm_aspqm2

echo "=== queue"
squeue -u "$USER" -o "%.10i %.9P %.10j %.2t %.10M %.20R" || true
