"""Build the hydrolysed-product run directory from the reactant's.

    python make_product_run.py <reactant_dir> <product_dir> <ligand.pdb>

The product complex is the reactant system with meropenem replaced by
hydrolysed meropenem at the pose transferred from SMB-1's 5AXO (D156, D159).
Everything else - the protein, the zinc and its placement, the force fields,
the histidine variants, the 2.5 A detection cutoff - is the reactant's, because
the two states are only comparable if nothing else differs.

What this changes, and why each one:

  protein_pdb     the reactant's start structure with the substrate removed and
                  **the hydroxide turned back into a water**. The nucleophile is
                  part of the product now, so the OH- residue cannot stay; but
                  the position is a crystal water in 6NC5 (HOH493, 2.14 A from
                  the metal) and 5AXO keeps a water on Zn1 at 1.82 A, so
                  deleting it outright would leave a three-coordinate zinc that
                  neither crystal supports.
  ligand          the grafted product, our atom names on 5AXO's coordinates.
  ligand_smiles   the product as the zwitterion: protonated pyrrolidine, both
                  carboxylates ionised, and the proton the hydrolysis delivers
                  on the leaving nitrogen. Net -1, which is D13's figure and
                  what `hyd_zw.pdb` actually carries.
  charge          0 -> -1, the ligand's.

Results are dropped from every section: a results block copied into a new
directory points at the other run's files, and a stage that reads one would
start from the wrong structure and converge perfectly well.
"""
from __future__ import annotations

import sys
from pathlib import Path

import yaml

# The hydrolysed product, zwitterionic: C[C@@H](O)... with [NH2+] and both
# carboxylates as [O-]. Derived from CSR_hydrolysed's neutral-amine SMILES by
# the two changes D13 specifies, and checked against hyd_zw.pdb's own
# hydrogens - N1x-H03, N2x-H01/H10x, all four carboxylate oxygens bare.
PRODUCT_SMILES = ("C[C@@H](O)[C@@H](C(=O)[O-])[C@@H]1NC(C(=O)[O-])="
                  "C(S[C@@H]2C[NH2+][C@H](C(=O)N(C)C)C2)[C@@H]1C")

HYDROXIDE_RESNAMES = ("OH", "HOH")


def product_protein(src: Path, dst: Path, ligand_resname: str = "mer",
                    hydroxide_resname: str = "OH") -> dict:
    """Write `src` without the substrate and with the hydroxide as a water.

    Only the hydroxide's oxygen is kept; `preprocess` adds the second hydrogen,
    which is the same route every other water in the system takes.
    """
    kept, dropped, converted = [], 0, 0
    for line in src.read_text().splitlines(keepends=True):
        if not line.startswith(("ATOM", "HETATM")):
            if not line.startswith("END"):
                kept.append(line)
            continue
        resname = line[17:20].strip()
        name = line[12:16].strip()
        if resname == ligand_resname:
            dropped += 1
            continue
        if resname == hydroxide_resname:
            if name.startswith("H"):          # the hydroxide's own proton
                dropped += 1
                continue
            kept.append(line[:17] + "HOH" + line[20:])
            converted += 1
            continue
        kept.append(line)
    kept.append("END\n")
    dst.write_text("".join(kept))
    return {"dropped": dropped, "converted": converted, "kept": len(kept) - 1}


def product_config(reactant_cfg: dict, protein: str, ligand: str,
                   resname: str) -> dict:
    """The reactant's config with the product's inputs and no results."""
    out = {"preprocess": {"inputs": dict(
        (reactant_cfg.get("preprocess") or {}).get("inputs") or {})}}
    inputs = out["preprocess"]["inputs"]
    inputs.update({"protein_pdb": protein, "ligand": ligand,
                   "ligand_resname": resname, "ligand_smiles": PRODUCT_SMILES,
                   "charge": -1, "ligand_xmlfile": None})
    for section in ("md", "qmmm"):
        body = reactant_cfg.get(section)
        if body and body.get("inputs"):
            out[section] = {"inputs": dict(body["inputs"])}
    return out


def main(reactant_dir, product_dir, ligand_pdb):
    src, dst = Path(reactant_dir), Path(product_dir)
    cfg = yaml.safe_load((src / "config.yaml").read_text())
    start = (cfg.get("preprocess", {}).get("inputs", {})
             .get("protein_pdb", "")).lstrip("./")
    start_path = src / start
    if not start_path.is_file():
        start_path = src / "preprocess_results" / Path(start).name
    if not start_path.is_file():
        raise SystemExit("cannot find the reactant's protein_pdb (%s)" % start)

    dst.mkdir(parents=True, exist_ok=True)
    protein = dst / "product_protein.pdb"
    counts = product_protein(start_path, protein)
    print("protein from %s: %d lines kept, %d atoms dropped, hydroxide -> "
          "water on %d atom(s)" % (start_path, counts["kept"],
                                   counts["dropped"], counts["converted"]))
    if counts["converted"] != 1:
        raise SystemExit(
            "expected exactly one hydroxide oxygen to become a water, found "
            "%d. The start structure is not what this script assumes; nothing "
            "else was written." % counts["converted"])

    ligand = Path(ligand_pdb)
    (dst / ligand.name).write_text(ligand.read_text())
    out = product_config(cfg, "./%s" % protein.name, "./%s" % ligand.name,
                         "hyd")
    (dst / "config.yaml").write_text(
        yaml.safe_dump(out, sort_keys=False, default_flow_style=False))
    print("wrote %s, %s and %s" % (protein, dst / ligand.name,
                                   dst / "config.yaml"))
    print("ligand charge -1, smiles %s" % PRODUCT_SMILES)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
