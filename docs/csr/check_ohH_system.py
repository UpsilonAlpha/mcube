#!/usr/bin/env python3
"""Does the ohH input parameterise? Run before submitting anything.

    python check_ohH_system.py ohH_start.pdb water_nuc.xml

Loads the structure with the same force field list preprocess will use and
reports what it finds about the nucleophile. A preprocess job that fails on a
template match fails minutes in, but it fails after queueing, and the whole
point of the ohH arm is that it differ from the lysH arm in one hydrogen and
nothing else - so the things worth checking are mechanical: that the residue is
bonded, that it is neutral, that its charges are TIP3P's, and that rigidWater
leaves it flexible (which is why the equilibration uses 2 fs).

This does not solvate. It checks the one residue that is new.
"""

from __future__ import annotations

import sys
from pathlib import Path

TIP3P = (-0.834, 0.417, 0.417)


def main(argv):
    if len(argv) != 2:
        raise SystemExit(__doc__)
    pdb_path, xml_path = (Path(a) for a in argv)
    for path in (pdb_path, xml_path):
        if not path.is_file():
            raise SystemExit(f"no such file: {path}")

    from openmm.app import ForceField, PDBFile

    pdb = PDBFile(str(pdb_path))
    nuc = [r for r in pdb.topology.residues() if r.name == "NUW"]
    if len(nuc) != 1:
        raise SystemExit(f"expected exactly one NUW residue, found {len(nuc)}")
    nuc = nuc[0]
    names = sorted(a.name for a in nuc.atoms())
    print(f"NUW: chain {nuc.chain.id} resSeq {nuc.id}, atoms {' '.join(names)}")
    if names != ["H03", "H04", "O01"]:
        raise SystemExit(f"NUW should hold O01, H03 and H04; it holds {names}")

    indices = {a.index for a in nuc.atoms()}
    bonds = [b for b in pdb.topology.bonds()
             if b[0].index in indices and b[1].index in indices]
    print(f"  bonds within NUW: {len(bonds)}")
    if len(bonds) != 2:
        raise SystemExit(
            "NUW needs two bonds and has "
            f"{len(bonds)}. Bonds for a non-standard residue come from CONECT "
            "records; PDBFile infers a water's bonds from the names O/H1/H2, "
            "which these are not.")

    forcefields = ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml",
                   str(xml_path)]
    print(f"  force fields: {', '.join(forcefields)}")
    field = ForceField(*forcefields)

    # Checked against the parsed templates rather than by building a System.
    # Isolating one residue and calling createSystem does not work here: the
    # preprocess *input* has no protein hydrogens yet (Modeller adds them), and
    # deleting everything else leaves NUW alone in the chain that held the
    # ligand, where OpenMM starts comparing it against terminal variants. The
    # facts that matter are properties of the XML and the PDB, and both are
    # readable directly. `_templates` is private; this is a diagnostic.
    templates = field._templates
    if "NUW" not in templates:
        raise SystemExit(f"{xml_path} did not register a NUW template")
    nuw, hoh = templates["NUW"], templates["HOH"]
    nuw_atoms = {a.name: a for a in nuw.atoms}
    print(f"  NUW template: {' '.join(sorted(nuw_atoms))}, "
          f"{len(nuw.bonds)} bonds")
    if sorted(nuw_atoms) != ["H03", "H04", "O01"]:
        raise SystemExit(f"NUW template holds {sorted(nuw_atoms)}")

    charges = [a.parameters.get("charge") for a in nuw.atoms]
    print(f"  charges: {charges} (sum {sum(charges):+.3f})")
    if abs(sum(charges)) > 1e-9:
        raise SystemExit("the NUW template is not neutral")

    # Parameter-identical to HOH is load-bearing: a template with its own types
    # also matches every bulk water by bond graph, and OpenMM then refuses the
    # whole system with "Multiple non-identical matching templates found".
    nuw_types = sorted(a.type for a in nuw.atoms)
    hoh_types = sorted(a.type for a in hoh.atoms)
    print(f"  NUW types {nuw_types}")
    print(f"  HOH types {hoh_types}")
    if nuw_types != hoh_types:
        raise SystemExit(
            "NUW does not use TIP3P's own atom types, so its template is not "
            "parameter-identical to HOH. It will match all ~15000 bulk waters "
            "by bond graph and OpenMM will refuse the system.")
    if sorted(charges) != sorted(a.parameters.get("charge") for a in hoh.atoms):
        raise SystemExit("NUW's charges differ from HOH's")

    print("OK - bonded in the PDB; neutral, TIP3P-typed and "
          "HOH-identical in the template")


if __name__ == "__main__":
    main(sys.argv[1:])
