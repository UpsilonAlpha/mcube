#!/bin/bash
# Submit `aspdft` - the def2-SVP QM/MM geometry optimisation (D233).
#
# The instrument of record for Asp92's position. GFN2 chelates the carboxylate
# onto Zn1 from a crystal-consistent start (D232), so only DFT can say where it
# sits. No QM/MM equilibration: the user's instruction, and the only affordable
# choice at this level.
#
# It starts from `relax10`'s last frame, NOT from the GFN2 chunk's. Starting
# from a geometry that already has Asp92 on the metal would hand DFT the answer
# and D206 rule 3 would exclude the result. The gate below refuses anything else.
#
# `opt_maxiter: 25` so the job reaches its own maxiter inside the wall even if a
# step costs 50 min: ASH raises there, mcube keeps the last geometry and writes
# `QM_MM_optimized.pdb` with `converged: false`, whereas a walltime kill writes
# no PDB. The next chunk starts from that PDB with a maxiter sized from the rate
# this one measures.
#
# 32 cores, not 8: DFT scales with cores where xTB does not.
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

echo "=== sections.py aspdft"
python "$R/sections.py" "$R/site.yaml" aspdft

echo "=== what landed in config.yaml"
python - <<'PY'
import os
import yaml

d = yaml.safe_load(open("config.yaml"))
s = d.get("qmmm_aspdft")
if not s:
    raise SystemExit("qmmm_aspdft is not in config.yaml - sections.py refused it")
i = s["inputs"]
for k in ("pdb_file", "charge", "qm_engine", "dft_functional", "dft_basis",
          "dft_cores", "skip_equilibrate", "skip_optimize", "opt_maxiter",
          "rigid_water"):
    print(f"    {k}: {i.get(k)}")
print("    qm_residues:", i.get("qm_residues"))
print("    actradius:", i.get("actradius", (d.get("qmmm") or {}).get("inputs", {}).get("actradius")))

# D134: a $set block appears the moment a qmmm run is given opt_constraints, and
# convergence then cannot be reached. Nothing may be constrained here - that is
# also the scientific point, since a restrained carboxylate measures the
# restraint. Refuse rather than warn.
for k in ("restraints", "opt_constraints", "constraints"):
    if i.get(k):
        raise SystemExit(f"{k} is set ({i[k]}) - this leg must constrain nothing (D134)")
if i.get("rigid_water") is not False:
    raise SystemExit("rigid_water must be false: a $set water distance gates convergence (D137)")

# The start has to be the crystal-consistent MM frame, not a GFN2 geometry that
# already has Asp92 on the metal (D232, D206 rule 3).
want = "md_relax10_results/NVTtrajectory_lastframe.pdb"
if i.get("pdb_file") != want:
    raise SystemExit(f"pdb_file is {i.get('pdb_file')!r}, expected {want!r}")
if not os.path.isfile(want):
    raise SystemExit(f"{want} does not exist")
if i.get("charge") != 2:
    raise SystemExit(f"charge is {i.get('charge')}, expected 2")
if str(i.get("qm_engine", "xtb")).lower() not in ("dft", "pyscf"):
    raise SystemExit(f"qm_engine is {i.get('qm_engine')!r} - this leg is the DFT one")
if str(i.get("dft_basis", "")).lower() != "def2-svp":
    raise SystemExit(f"dft_basis is {i.get('dft_basis')!r}, expected def2-svp")
if i.get("skip_equilibrate") is not True or i.get("skip_optimize") is not False:
    raise SystemExit("this leg optimises and does not equilibrate: "
                     f"skip_equilibrate={i.get('skip_equilibrate')}, "
                     f"skip_optimize={i.get('skip_optimize')}")
print("    [OK] def2-SVP optimisation from relax10's last frame, charge 2, "
      "nothing constrained, no QM/MM MD")
PY

echo "=== sbatch"
sbatch -J aspdft -p general \
  --cpus-per-task=32 --mem=120G -t 24:00:00 \
  "$R/run.slurm" qmmm qmmm_aspdft

echo "=== queue"
squeue -u "$USER" -o "%.10i %.9P %.10j %.2t %.10M %.20R" || true
