#!/usr/bin/env python3
"""
Put one ion into the bulk water of the equilibrated apo channel box (system build step 2 in
channel-us-plan.md): the water whose oxygen is nearest (x, y, z) is replaced by the ion, and
Na+ ions are removed (those furthest from the protein) to keep the box neutral. x, y default
to the pore axis, the protein's Cα centroid in xy.

The ion goes in as a chain of its own at the end, so mcube's `ion: "<chain>:1:<ION>"` selects
it; the chain letter OpenMM writes is printed and stored in the report.

Bonds come from --topology (default: finalsystem.pdb next to --pdb, the file the membrane
builder wrote). md's imaged last frame has no CONECT records, and without them POPE has no
bonds when the file is read back, so md fails with "No template found for residue (POP)".
Positions and box come from --pdb; the two files must list the same atoms in the same order.
The imaged last frame also splits molecules across the periodic boundary (1300 lipid bonds
were 3-169 Å long), so every molecule is made whole first: walking its bonds, each atom is put
at the periodic image nearest the atom it is bonded to. The report gives the longest bond
before and after.

  python insert_ion.py --pdb ../apo/md_lastframe.pdb --ion LA --z -30 --out channel_LA.pdb
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

ION = {"CA": ("Ca", 2), "LA": ("La", 3), "LU": ("Lu", 3)}


def make_whole(xyz, top, box):
    """Undo periodic wrapping inside every molecule (orthorhombic box), by bond-graph walk."""
    from collections import defaultdict, deque
    nbrs = defaultdict(list)
    for a, b in top.bonds():
        nbrs[a.index].append(b.index)
        nbrs[b.index].append(a.index)
    seen = np.zeros(len(xyz), bool)
    out = xyz.copy()
    for start in range(len(xyz)):
        if seen[start]:
            continue
        seen[start] = True
        queue = deque([start])
        while queue:
            i = queue.popleft()
            for j in nbrs[i]:
                if not seen[j]:
                    d = out[j] - out[i]
                    out[j] = out[i] + d - box * np.round(d / box)
                    seen[j] = True
                    queue.append(j)
    return out


def longest_bond(xyz, top):
    return max(float(np.linalg.norm(xyz[a.index] - xyz[b.index])) for a, b in top.bonds())


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pdb", required=True)
    p.add_argument("--ion", required=True, choices=sorted(ION))
    p.add_argument("--z", type=float, required=True, help="Å, in the PDB's frame")
    p.add_argument("--x", type=float)
    p.add_argument("--y", type=float)
    p.add_argument("--protein-chain", default="A")
    p.add_argument("--topology", help="PDB with the bonds (default: finalsystem.pdb beside --pdb)")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    from openmm import unit
    from openmm.app import Element, Modeller, PDBFile, Topology
    pdb = PDBFile(args.pdb)
    topo_file = args.topology or str(Path(args.pdb).resolve().parent / "finalsystem.pdb")
    bonded = PDBFile(topo_file)
    top = bonded.topology
    names = [(a.name, a.residue.name, a.residue.chain.id) for a in pdb.topology.atoms()]
    if names != [(a.name, a.residue.name, a.residue.chain.id) for a in top.atoms()]:
        raise SystemExit(f"{topo_file} and {args.pdb} do not list the same atoms in the same order")
    top.setPeriodicBoxVectors(pdb.topology.getPeriodicBoxVectors())
    box = np.array(pdb.topology.getUnitCellDimensions().value_in_unit(unit.angstrom))
    xyz = np.array(pdb.getPositions().value_in_unit(unit.angstrom))
    bond_before = longest_bond(xyz, top)
    xyz = make_whole(xyz, top, box)
    bond_after = longest_bond(xyz, top)
    if bond_after > 3.0:
        raise SystemExit(f"a bond is still {bond_after:.1f} A after making molecules whole")
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

    mod = Modeller(top, (xyz / 10.0) * unit.nanometer)
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
              "topology": topo_file, "bonds": mod.topology.getNumBonds(),
              "longest_bond_A": {"input": round(bond_before, 2), "written": round(bond_after, 2)},
              "umbrella_ion": f"{letter}:1:{args.ion}",
              "box_A": [round(float(v), 2) for v in pdb.topology.getUnitCellDimensions().value_in_unit(unit.angstrom)]}
    Path(args.out).with_suffix(".json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
