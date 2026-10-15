#!/usr/bin/env python3
"""
Build a solvated NEB product that shares the reactant's solvent box.

NEB needs reactant and product with identical atoms in identical order, including
water and ions. This takes the solvated reactant (e.g. preprocess finalsystem.pdb)
and replaces the solute coordinates with those of an unsolvated product structure
whose solute atoms are in the same order (same names and elements).

Solute = every atom that is not water or a solvent ion. The product is aligned
onto the reactant using the solute atoms outside `--reactive-resnames`, in case
solvation moved the solute. Waters and ions keep the reactant's coordinates,
so check the reported clashes and QM/MM-optimise both endpoints afterwards.

Usage:
    python build_product_box.py finalsystem.pdb hydrolysed_lastframe.pdb product_solvated.pdb \
        --reactive-resnames MER OH HYD mer
"""

import argparse

import numpy as np
from openmm import unit
from openmm.app import PDBFile

SOLVENT = {"HOH", "WAT", "NA", "CL", "K", "Na+", "Cl-"}


def kabsch(mobile, target):
    """Rotation and translation that best superpose mobile onto target."""
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    h = (mobile - mc).T @ (target - tc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, tc - mc @ rot.T


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("reactant_solvated", help="Solvated reactant PDB (defines topology, solvent and box)")
    p.add_argument("product_solute", help="Product PDB with the same solute atoms in the same order")
    p.add_argument("output", help="Output solvated product PDB")
    p.add_argument("--reactive-resnames", nargs="+", default=["MER", "HYD", "OH", "mer", "hyd"],
                   help="Residues left out of the alignment (they move in the reaction)")
    p.add_argument("--clash", type=float, default=1.5, help="Water/ion clash distance (A)")
    args = p.parse_args()

    react = PDBFile(args.reactant_solvated)
    prod = PDBFile(args.product_solute)

    r_atoms = [a for a in react.topology.atoms() if a.residue.name not in SOLVENT]
    p_atoms = [a for a in prod.topology.atoms() if a.residue.name not in SOLVENT]

    if len(r_atoms) != len(p_atoms):
        raise SystemExit(f"Solute atom counts differ: reactant {len(r_atoms)}, product {len(p_atoms)}")
    for ra, pa in zip(r_atoms, p_atoms):
        if ra.name != pa.name or ra.element != pa.element:
            raise SystemExit(f"Atom order differs at reactant atom {ra.index} ({ra.residue.name} {ra.name}) "
                             f"vs product atom {pa.index} ({pa.residue.name} {pa.name})")

    r_xyz = react.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    p_xyz = prod.getPositions(asNumpy=True).value_in_unit(unit.angstrom)

    r_idx = np.array([a.index for a in r_atoms])
    p_idx = np.array([a.index for a in p_atoms])
    fixed = np.array([a.residue.name not in args.reactive_resnames for a in r_atoms])

    rot, trans = kabsch(p_xyz[p_idx][fixed], r_xyz[r_idx][fixed])
    p_solute = p_xyz[p_idx] @ rot.T + trans
    rmsd_fixed = np.sqrt(((p_solute[fixed] - r_xyz[r_idx][fixed]) ** 2).sum(axis=1).mean())
    moved = np.linalg.norm(p_solute - r_xyz[r_idx], axis=1)
    print(f"Aligned on {fixed.sum()} non-reactive solute atoms: RMSD {rmsd_fixed:.3f} A")
    print(f"Largest solute displacement reactant->product: {moved.max():.2f} A "
          f"(atom {r_idx[moved.argmax()]}, {r_atoms[moved.argmax()].residue.name} {r_atoms[moved.argmax()].name})")

    new_xyz = r_xyz.copy()
    new_xyz[r_idx] = p_solute

    # Solvent atoms now closer than --clash to any product solute atom
    solvent_idx = np.array([a.index for a in react.topology.atoms() if a.residue.name in SOLVENT])
    clashes = []
    if len(solvent_idx):
        moved_idx = r_idx[moved > 0.1]
        for i in moved_idx:
            d = np.linalg.norm(new_xyz[solvent_idx] - new_xyz[i], axis=1)
            for j in solvent_idx[d < args.clash]:
                clashes.append((i, j, np.linalg.norm(new_xyz[j] - new_xyz[i])))
    atoms = list(react.topology.atoms())
    print(f"Solvent clashes (< {args.clash} A) with moved solute atoms: {len(clashes)}")
    for i, j, d in sorted(clashes, key=lambda c: c[2])[:20]:
        print(f"  solute {i} ({atoms[i].residue.name} {atoms[i].name}) - "
              f"{atoms[j].residue.name}{atoms[j].residue.id} {atoms[j].name} ({j}): {d:.2f} A")

    with open(args.output, "w") as f:
        PDBFile.writeFile(react.topology, new_xyz * unit.angstrom, f, keepIds=True)
    print(f"Wrote {args.output} ({len(atoms)} atoms, same topology and box as {args.reactant_solvated})")


if __name__ == "__main__":
    main()
