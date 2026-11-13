"""Build the hydrolysed-product run directory from the reactant's.

    python make_product_run.py <reactant_dir> <product_dir> <ligand.pdb>

The product complex is the reactant system with meropenem replaced by
hydrolysed meropenem at the pose transferred from SMB-1's 5AXO (¶156, ¶159).
Everything else - the protein, the zinc and its placement, the force fields,
the histidine variants, the detection cutoff - is the reactant's, because
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
                  on the leaving nitrogen. Net -1, which is ¶13's figure and
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
# the two changes ¶13 specifies, and checked against hyd_zw.pdb's own
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
    # Where the existing waters live, so the converted hydroxide can join them.
    # It cannot stay in the ligand's chain: removing the substrate leaves that
    # chain starting at residue 2, and mcube refuses a file not numbered 1..N
    # per chain because ASH renumbers them and every selector would then point
    # somewhere else. It is a water now, so the water chain is where it goes.
    waters = [l for l in src.read_text().splitlines(keepends=True)
              if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == "HOH"]
    water_chain = waters[0][21] if waters else None
    next_water = max((int(l[22:26]) for l in waters), default=0) + 1

    kept, dropped, converted, carried = [], 0, 0, None
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
            # " O01" is the hydroxide's name, and `amber14/tip3p.xml` types a
            # water by its atom names: a HOH whose oxygen is called O01 has no
            # template and the system cannot be built (¶163).
            moved = line[:12] + " O  " + line[16:17] + "HOH" + line[20:]
            if water_chain:
                moved = (moved[:21] + water_chain + "%4d" % next_water
                         + moved[26:])
            carried = moved
            converted += 1
            continue
        kept.append(line)
    # After the other waters, not where the hydroxide sat: the numbering check
    # walks the file in order, so a chain whose residues appear 2 then 1 reads
    # as "starts at 2" and is refused (¶162).
    if carried is not None:
        last = max((i for i, l in enumerate(kept)
                    if l.startswith(("ATOM", "HETATM"))), default=len(kept) - 1)
        kept.insert(last + 1, carried)
    kept.append("END\n")
    dst.write_text("".join(kept))
    return {"dropped": dropped, "converted": converted, "kept": len(kept) - 1}


def product_config(reactant_cfg: dict, protein: str, ligand: str,
                   resname: str) -> dict:
    """The reactant's config with the product's inputs and no results."""
    out = {"preprocess": {"inputs": dict(
        (reactant_cfg.get("preprocess") or {}).get("inputs") or {})}}
    inputs = out["preprocess"]["inputs"]
    # `merge_ligand: false` is right for the reactant, whose substrate is
    # already inside `csr_crystal_start.pdb` - there the ligand file only
    # supplies parameters. Here the substrate has been removed on purpose and
    # the product exists only in its own file, so inheriting `false` builds a
    # system with no ligand at all: preprocess reported [DONE], the zinc shell
    # was correct, and nothing in it was the product (¶165).
    inputs.update({"protein_pdb": protein, "ligand": ligand,
                   "ligand_resname": resname, "ligand_smiles": PRODUCT_SMILES,
                   "charge": -1, "ligand_xmlfile": None,
                   "merge_ligand": True})
    retarget_ligand_chain(inputs, LIGAND_CHAIN)
    for section in ("md", "qmmm"):
        body = reactant_cfg.get(section)
        if body and body.get("inputs"):
            out[section] = {"inputs": retarget_ligand_xml(
                dict(body["inputs"]), resname)}
    return out


FF_KEYS = ("forcefields", "default_forcefields")

# Where ASH's merge puts the ligand. `product_protein.pdb` holds the protein in
# A, the zinc in B and the waters in W; the merge renames W to C and appends the
# ligand as a chain of its own, which for this system is D. The reactant's
# selectors say C:1 because there the substrate was inside the protein file, so
# inherited unchanged they name a **water** - and nothing complains, because a
# water at C:1 exists. That cost two preprocess runs: the ligand stayed out of
# the QM region and the 1.90 A water lost its restraint, since the exclusion
# meant for the ligand landed on it (¶165).
LIGAND_CHAIN = "D"

# Asp92 is in the reactant's QM region, but only because that run was built
# with the old 2.8 A default cutoff and its OD2 sits at 2.80: at the 2.5 A this
# config now sets, detection misses it. It is the general base, so it belongs in
# both QM regions by intent rather than by a cutoff that happens to reach.
GENERAL_BASE = 92

SELECTOR_KEYS = ("extra_qm_residues", "constraint_exclude_residues",
                 "exclude_residues")


def retarget_ligand_chain(inputs: dict, chain: str) -> dict:
    """Move the reactant's ligand selectors onto the merged ligand's chain."""
    old = (inputs.get("ligand_chain") or "C")
    for key in SELECTOR_KEYS:
        block = inputs.get(key)
        if not isinstance(block, dict) or old not in block:
            continue
        block = dict(block)
        block[chain] = block.pop(old)
        inputs[key] = block
    extra = inputs.get("extra_qm_residues")
    if isinstance(extra, dict) and "A" in extra:
        if GENERAL_BASE not in extra["A"]:
            inputs["extra_qm_residues"] = dict(
                extra, A=sorted(list(extra["A"]) + [GENERAL_BASE]))
    return inputs


def retarget_ligand_xml(inputs: dict, resname: str) -> dict:
    """Point inherited force-field lists at this ligand's openff XML.

    The reactant's `qmmm.inputs.forcefields` names `openff_MER.xml`, which is
    the substrate's parameters: inherited unchanged it would run the product
    with the reactant's ligand typed in, and OpenMM would either fail to find a
    template or quietly use the wrong one. The name follows the resname, so it
    is derivable rather than a thing to remember.
    """
    want = "openff_%s.xml" % resname
    for key in FF_KEYS:
        entries = inputs.get(key)
        if not entries:
            continue
        inputs[key] = [
            ("%s/%s" % (e.rsplit("/", 1)[0], want) if "/" in e else want)
            if Path(e).name.startswith("openff_") and e.endswith(".xml") else e
            for e in entries]
    return inputs


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

    # Into preprocess_results/, not the run directory: `run.slurm` does
    # `cd "$OUT"` before calling mcube, so preprocess reads its `protein_pdb`
    # and `ligand` relative to **that** directory, not to the config's. It is
    # why the reactant's own `csr_crystal_start.pdb` lives there too, and
    # writing them one level up cost job 31299717 an opaque RDKit "Bad input
    # file" (¶161).
    here = dst / "preprocess_results"
    here.mkdir(parents=True, exist_ok=True)
    protein = here / "product_protein.pdb"
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
    (here / ligand.name).write_text(ligand.read_text())
    out = product_config(cfg, "./%s" % protein.name, "./%s" % ligand.name,
                         "hyd")
    (dst / "config.yaml").write_text(
        yaml.safe_dump(out, sort_keys=False, default_flow_style=False))
    print("wrote %s, %s and %s" % (protein, here / ligand.name,
                                   dst / "config.yaml"))
    print("ligand charge -1, smiles %s" % PRODUCT_SMILES)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
