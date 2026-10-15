#!/usr/bin/env python3
"""
Make every residue's atoms contiguous in a PDB, and apply the same atom
permutation to other PDBs with the same atoms in the same order.

PDB readers (OpenMM, pdbfixer, RDKit) start a new residue whenever the residue
number changes between consecutive lines, so a residue whose atoms are
interleaved with another residue's (e.g. hydroxide atoms placed inside a ligand
to match a reaction product's atom order) is read as several fragments.

The permutation is taken from the first (reference) file: residues in order of
first appearance, atoms within a residue in their original order. Every file
gets that same permutation, so atom i still corresponds to atom i across files
(names and elements are checked). Serial numbers are renumbered and CONECT
records remapped. Output: <name>_contig.pdb next to each input.

Usage:
    python make_residues_contiguous.py reactant.pdb [product.pdb ...]
"""

import argparse
from pathlib import Path


def read(path):
    lines = Path(path).read_text().splitlines()
    atoms = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
    conect = [l for l in lines if l.startswith("CONECT")]
    header = [l for l in lines if l.startswith(("CRYST1", "REMARK", "HEADER", "TITLE"))]
    return header, atoms, conect


def residue_key(line):
    # chain, resSeq, insertion code, resName
    return line[21], line[22:26], line[26], line[17:20]


def ident(line):
    # atom name and element
    return line[12:16].strip(), line[76:78].strip()


def permutation(atoms):
    order, groups = [], {}
    for i, line in enumerate(atoms):
        key = residue_key(line)
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(i)
    return [i for key in order for i in groups[key]]


def write(path, header, atoms, conect, perm):
    old_serial = [int(l[6:11]) for l in atoms]
    new_serial = {old_serial[old]: new + 1 for new, old in enumerate(perm)}

    out = list(header)
    prev_chain = None
    for new, old in enumerate(perm):
        line = atoms[old]
        if prev_chain is not None and line[21] != prev_chain:
            out.append("TER")
        prev_chain = line[21]
        out.append(f"{line[:6]}{new + 1:5d}{line[11:]}")
    out.append("TER")

    for line in conect:
        fields = [int(line[i:i + 5]) for i in range(6, len(line.rstrip()), 5) if line[i:i + 5].strip()]
        if fields and all(f in new_serial for f in fields):
            out.append("CONECT" + "".join(f"{new_serial[f]:5d}" for f in fields))
    out.append("END")
    Path(path).write_text("\n".join(out) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("reference", help="PDB that defines the permutation (e.g. the reactant)")
    p.add_argument("others", nargs="*", help="PDBs with the same atoms in the same order (e.g. the product)")
    args = p.parse_args()

    ref_header, ref_atoms, ref_conect = read(args.reference)
    perm = permutation(ref_atoms)
    moved = sum(1 for new, old in enumerate(perm) if new != old)
    print(f"{args.reference}: {len(ref_atoms)} atoms, {moved} change position")

    for path in [args.reference] + args.others:
        header, atoms, conect = read(path)
        if len(atoms) != len(ref_atoms):
            raise SystemExit(f"{path} has {len(atoms)} atoms, reference has {len(ref_atoms)}")
        for i, (a, b) in enumerate(zip(ref_atoms, atoms)):
            if ident(a) != ident(b):
                raise SystemExit(f"{path} atom {i + 1} is {ident(b)}, reference has {ident(a)}")
        out = Path(path).with_name(Path(path).stem + "_contig.pdb")
        write(out, header, atoms, conect, perm)
        print(f"Wrote {out}")


if __name__ == "__main__":
    main()
