"""Build CSR_mich, the Michaelis rebuild's run directory (D196).

    python build-mich.py [--force]

Run from anywhere; every path below is absolute. It creates the directory,
copies the four files a run needs, and writes `config.yaml`. It does NOT
submit anything and it does NOT run `sections.py`: `eq5` names
`preprocess_results/finalsystem.pdb` so that its five selectors are checked
against the real solvated system, and that file does not exist until
`preprocess` has run. The order is therefore

    1. this script
    2. sbatch -J mich_pre --cpus-per-task=2 --mem=8G run.slurm preprocess preprocess
    3. read the coordination report - Asp92 OD2 near 2.80, the carbonyl near
       2.60, and `removed_residues` naming D:1:HOH
    4. python sections.py site.yaml eq5      (checks the five selectors)
    5. sbatch -J mich_eq5 -p gpu_cuda --qos=gpu --gres=gpu:1 \
              --cpus-per-task=16 --mem=32G run.slurm md md_eq5

`run.slurm` requires the section as its second argument even for
preprocess, where it only names the output directory, and it writes its
own output to `protdft/`, so that directory has to exist before the
first sbatch. Both are done here.

Why a script rather than typed commands: this directory's provenance is part
of the result. The config it writes says which structure it came from, what
was removed from it, and what was changed relative to CSR_crystal_v2, and a
file in the repository can be read afterwards where a shell history cannot.

What differs from CSR_crystal_v2, and why:

  remove_residues: [D:1:HOH]   the water the substrate displaces. The input
                               has exactly B:1:ZN, C:1:mer (51 atoms),
                               C:2:OH (2) and D:1:HOH (3), and the QM/MM
                               frame measured that water at 1.98 A on the
                               metal - so removing it leaves the five of D196
                               rather than six.
  coordination_force: 500      D191 measured a 2.80 A restraint sitting at a
                               mean of 2.44; offset goes as F/k.
  ligand: null                 the ligand is NOT re-parameterised. The XML
  ligand_xmlfile: openff_MER   CSR_crystal_v2 already produced is copied in
                               and used as given, which avoids both the cost
                               and OpenFF's atom renaming (D169) - the input
                               structure's ligand atom names are the ones that
                               XML was generated for.
  extra_qm_residues            Arg90 and Lys226 dropped: D190 put the whole
                               RQK set in MM, and leaving them here would make
                               `preprocess.results.qm_atoms` disagree with
                               every section that runs.
"""
import argparse
import shutil
import sys
from pathlib import Path

import yaml

SCRATCH = Path("/scratch/user/uquadhik")
# The finished run directories were moved into `CSR/` on 2026-11-03 to clear the
# top of scratch. `CSR_mich` is the one that stayed out, because its MD was
# running at the time and Slurm had recorded the old WorkDir; move it in once
# that job is done and then RUNS is the only line that has to change.
RUNS = SCRATCH / "CSR"
SRC = RUNS / "CSR_crystal_v2"
DST = SCRATCH / "CSR_mich"
STRUCTURE = "csr_crystal_restrained_end.pdb"
LIGAND_XML = "openff_MER.xml"

# Chain, residue, name of every non-protein residue the input must contain for
# this build to be the one described above. Checked, not assumed: the whole
# point of the rebuild is which oxygen is on the metal.
EXPECTED = {"B:1:ZN": 1, "C:1:mer": 51, "C:2:OH": 2, "D:1:HOH": 3}


def non_protein(pdb: Path) -> dict:
    counts: dict = {}
    for line in pdb.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")) or line[21] == "A":
            continue
        key = "%s:%s:%s" % (line[21], line[22:26].strip(), line[17:20].strip())
        counts[key] = counts.get(key, 0) + 1
    return counts


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing CSR_mich/config.yaml")
    args = ap.parse_args()

    if DST.exists() and not args.force:
        print(f"{DST} already exists. Re-run with --force to rewrite its "
              "config.yaml, or move it aside first.")
        return 1

    structure = SRC / "preprocess_results" / STRUCTURE
    ligand_xml = SRC / "preprocess_results" / LIGAND_XML
    hydroxide = SRC / "hydroxide.xml"
    slurm = RUNS / "CSR_wt" / "run.slurm"
    for path in (structure, ligand_xml, hydroxide, slurm, SRC / "config.yaml"):
        if not path.is_file():
            print(f"missing: {path}")
            return 1

    found = non_protein(structure)
    if found != EXPECTED:
        print("The input structure is not the one this build is described "
              f"against.\n  wanted {EXPECTED}\n  found  {found}")
        return 1
    print(f"[OK] {STRUCTURE}: " +
          ", ".join(f"{k} x{v}" for k, v in sorted(found.items())))

    src_config = yaml.safe_load(SRC.joinpath("config.yaml").read_text())
    pre = dict(src_config["preprocess"]["inputs"])
    md = dict(src_config["md"]["inputs"])

    DST.mkdir(parents=True, exist_ok=True)
    # run.slurm's --output is protdft/%x-%j.out, and sbatch fails
    # outright if that directory is missing.
    (DST / "protdft").mkdir(exist_ok=True)
    (DST / "preprocess_results").mkdir(exist_ok=True)
    for path in (structure, ligand_xml, hydroxide, slurm):
        shutil.copy2(path, DST / path.name)
        print(f"[COPY] {path.name}")

    pre["protein_pdb"] = (DST / STRUCTURE).as_posix()
    pre["ligand"] = None
    pre["ligand_xmlfile"] = (DST / LIGAND_XML).as_posix()
    pre["remove_residues"] = ["D:1:HOH"]
    pre["coordination_force"] = 500.0
    pre["extra_qm_residues"] = {"C": [1]}
    pre.pop("results", None)

    md["pdb_file"] = None
    md["nvt_timestep"] = 0.002

    config = {"preprocess": {"inputs": pre}, "md": {"inputs": md}}
    out = DST / "config.yaml"
    out.write_text(yaml.safe_dump(config, sort_keys=False))
    print(f"[WRITE] {out}")
    print(f"\nNext, from {DST}:")
    print("  sbatch -J mich_pre --cpus-per-task=2 --mem=8G "
          "run.slurm preprocess preprocess")
    return 0


if __name__ == "__main__":
    sys.exit(main())
