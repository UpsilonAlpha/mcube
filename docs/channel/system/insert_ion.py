#!/usr/bin/env python3
"""
Put one ion into the bulk water of the equilibrated apo channel box (system build step 2 in
channel-us-plan.md): the water whose oxygen is nearest (x, y, z) is replaced by the ion, and
Na+ ions are removed (those furthest from the protein) to keep the box neutral. x, y default
to the pore axis, the protein's Cα centroid in xy.

The ion goes in as a chain of its own at the end, so mcube's `ion: "<chain>:1:<ION>"` selects
it; the chain letter OpenMM writes is printed and stored in the report.

  python insert_ion.py --pdb ../apo/md_lastframe.pdb --ion LA --z -30 --out channel_LA.pdb
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

ION = {"CA": ("Ca", 2), "LA": ("La", 3), "LU": ("Lu", 3)}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pdb", required=True)
    p.add_argument("--ion", required=True, choices=sorted(ION))
    p.add_argument("--z", type=float, required=True, help="Å, in the PDB's frame")
    p.add_argument("--x", type=float)
    p.add_argument("--y", type=float)
    p.add_argument("--protein-chain", default="A")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    from openmm import unit
    from openmm.app import Element, Modeller, PDBFile, Topology
    pdb = PDBFile(args.pdb)
    top = pdb.topology
    xyz = np.array(pdb.getPositions().value_in_unit(unit.angstrom))
    prot = [a.index for a in top.atoms() if a.residue.chain.id == args.protein_chain and a.name == "CA"]
    cx, cy = xyz[prot, 0].mean(), xyz[prot, 1].mean()
    target = np.array([args.x if args.x is not None else cx, args.y if args.y is not None else cy, args.z])

    waters = [r for r in top.residues() if r.name in ("HOH", "WAT")]
    o = {r: next(a.index for a in r.atoms() if a.element.symbol == "O") for r in waters}
    water = min(waters, key=lambda r: np.linalg.norm(xyz[o[r]] - target))
    site = xyz[o[water]].copy()

    sym, charge = ION[args.ion]
    sodium = [r for r in top.residues() if r.name in ("NA", "Na+")]
    if len(sodium) < charge:
        raise SystemExit(f"only {len(sodium)} Na+ in the box; cannot neutralise a {charge}+ ion")
    heavy = np.array([a.index for a in top.atoms() if a.residue.chain.id == args.protein_chain])
    far = sorted(sodium, key=lambda r: -np.min(np.linalg.norm(
        xyz[heavy] - xyz[next(iter(r.atoms())).index], axis=1)))[:charge]

    mod = Modeller(top, pdb.positions)
    mod.delete([water] + far)
    ion_top = Topology()
    chain = ion_top.addChain()
    res = ion_top.addResidue(args.ion, chain)
    ion_top.addAtom(args.ion, Element.getBySymbol(sym), res)
    mod.add(ion_top, [site / 10.0] * unit.nanometer)
    with open(args.out, "w") as fh:
        PDBFile.writeFile(mod.topology, mod.positions, fh)

    chains = list(mod.topology.chains())
    letter = None
    for line in open(args.out):
        if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == args.ion:
            letter = line[21]
    report = {"ion": args.ion, "replaced_water_O_A": [round(float(v), 2) for v in site],
              "pore_axis_xy_A": [round(float(cx), 2), round(float(cy), 2)],
              "removed_na": len(far), "chains": len(chains), "ion_chain": letter,
              "umbrella_ion": f"{letter}:1:{args.ion}",
              "box_A": [round(float(v), 2) for v in pdb.topology.getUnitCellDimensions().value_in_unit(unit.angstrom)]}
    Path(args.out).with_suffix(".json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
