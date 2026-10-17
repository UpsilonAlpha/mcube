#!/usr/bin/env python3
"""
Post-SCF DFT-D3 dispersion energy of the QM region of a QM/MM structure, for
single points that were run without dispersion (the correction is additive).

The QM region is rebuilt as mcube/ASH build it: the QM atoms (from the config
section's results, else preprocess) plus a hydrogen link atom on every QM-MM
bond, placed along the bond at ASH's standard distance (C-H 1.09, N-H 0.99,
O-H 0.98 Å). Needs `dftd3-python` (simple-dftd3).

Usage:
    python d3_correction.py structure.pdb config.yaml [--section qmmm_x] \
        [--method wb97x] [--damping d3bj]
"""

import argparse

import numpy as np
import yaml
from openmm import unit
from openmm.app import PDBFile

LINK_DISTANCE = {"C": 1.09, "N": 0.99, "O": 0.98}
ANGSTROM_TO_BOHR = 1.0 / 0.52917721092


def qm_region(pdb_file, qm_atoms):
    pdb = PDBFile(pdb_file)
    xyz = np.array(pdb.getPositions(asNumpy=True).value_in_unit(unit.angstrom))
    atoms = list(pdb.topology.atoms())
    qm = set(qm_atoms)
    symbols = [atoms[i].element.symbol for i in sorted(qm)]
    coords = [xyz[i] for i in sorted(qm)]
    links = []
    for a, b in pdb.topology.bonds():
        i, j = a.index, b.index
        if (i in qm) == (j in qm):
            continue
        q, m = (i, j) if i in qm else (j, i)
        d = LINK_DISTANCE.get(atoms[q].element.symbol, 1.09)
        v = xyz[m] - xyz[q]
        coords.append(xyz[q] + d * v / np.linalg.norm(v))
        symbols.append("H")
        links.append((q, m))
    return symbols, np.array(coords), links


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pdb")
    p.add_argument("config")
    p.add_argument("--section", help="Config section whose results hold qm_atoms (default: preprocess)")
    p.add_argument("--method", default="wb97x", help="Functional name in simple-dftd3")
    p.add_argument("--damping", default="d3bj", choices=["d3bj", "d3zero"])
    args = p.parse_args()

    cfg = yaml.safe_load(open(args.config))
    res = (cfg.get(args.section) or {}).get("results") or {} if args.section else {}
    qm_atoms = res.get("qm_atoms") or cfg["preprocess"]["results"]["qm_atoms"]
    symbols, coords, links = qm_region(args.pdb, qm_atoms)

    from mcube.qmmm import d3_energy
    e = d3_energy(symbols, coords * ANGSTROM_TO_BOHR, args.method, args.damping)
    print(f"{len(qm_atoms)} QM atoms + {len(links)} link atoms {links}")
    print(f"{args.damping} dispersion ({args.method}): {e:.8f} Eh ({e * 627.509474:.3f} kcal/mol)")


if __name__ == "__main__":
    main()
