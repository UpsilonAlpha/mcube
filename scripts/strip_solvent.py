#!/usr/bin/env python3
"""
Remove bulk water and monatomic ions from a (solvated) PDB, keeping chosen
waters, so an MD end point can be fed back into `mcube preprocess` as
`protein_pdb` (it re-solvates and re-adds ions). Atom names, residue ids,
chains, hydrogens and ligand CONECT records are kept (written with OpenMM).

Usage:
    python strip_solvent.py in.pdb out.pdb [--keep D] [--keep C:2]

--keep takes a chain ("D") or a residue ("chain:resid"); those waters/ions stay.
"""

import argparse

from openmm.app import PDBFile, Modeller

WATER_IONS = {"HOH", "WAT", "SOL", "TIP3", "NA", "CL", "K", "Na+", "Cl-", "K+", "NA+", "CL-"}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pdb_in")
    p.add_argument("pdb_out")
    p.add_argument("--keep", action="append", default=[], help="Chain or chain:resid to keep")
    args = p.parse_args()

    chains = {k for k in args.keep if ":" not in k}
    residues = {tuple(k.split(":")) for k in args.keep if ":" in k}
    pdb = PDBFile(args.pdb_in)
    modeller = Modeller(pdb.topology, pdb.positions)
    drop = [r for r in modeller.topology.residues()
            if r.name in WATER_IONS and r.chain.id not in chains and (r.chain.id, r.id.strip()) not in residues]
    kept = [f"{r.name} {r.chain.id}:{r.id.strip()}" for r in modeller.topology.residues()
            if r.name in WATER_IONS and r not in drop]
    modeller.delete(drop)
    with open(args.pdb_out, "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f, keepIds=True)
    print(f"Removed {len(drop)} water/ion residues; kept {kept}; "
          f"{modeller.topology.getNumAtoms()} atoms -> {args.pdb_out}")


if __name__ == "__main__":
    main()
