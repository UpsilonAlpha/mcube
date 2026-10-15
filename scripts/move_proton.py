#!/usr/bin/env python3
"""
Move a proton to a different heavy atom within one residue of a PDB, keeping
the atom order (and so the atom correspondence between files) unchanged.

Used to change a ligand's protonation state, e.g. meropenem COOH/amine ->
zwitterion COO-/NH2+ by moving H01 from O3x to N2x. The new position points
away from the new parent's other bonded neighbours (tetrahedral-ish for sp3).
CONECT records are updated. Optionally writes the residue on its own
(renumbered, internal CONECT only) for re-parameterisation.

Usage:
    python move_proton.py in.pdb out.pdb --resname mer --proton H01 --to N2x \
        [--extract ligand.pdb]
"""

import argparse

import numpy as np

BOND = {"N": 1.02, "O": 0.97, "S": 1.34, "C": 1.09}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pdb")
    p.add_argument("out")
    p.add_argument("--resname", required=True)
    p.add_argument("--proton", required=True, help="Name of the H atom to move")
    p.add_argument("--to", required=True, help="Name of the new parent atom (same residue)")
    p.add_argument("--extract", help="Also write the residue alone to this PDB")
    args = p.parse_args()

    lines = open(args.pdb).read().splitlines()
    atom_idx = [i for i, l in enumerate(lines) if l.startswith(("ATOM", "HETATM"))]
    serial = {int(lines[i][6:11]): i for i in atom_idx}

    def find(name):
        hits = [i for i in atom_idx if lines[i][17:20].strip() == args.resname and lines[i][12:16].strip() == name]
        if len(hits) != 1:
            raise SystemExit(f"Expected one {args.resname} {name}, found {len(hits)}")
        return hits[0]

    h_line, n_line = find(args.proton), find(args.to)
    h_ser, n_ser = int(lines[h_line][6:11]), int(lines[n_line][6:11])

    bonds = {}
    for l in lines:
        if l.startswith("CONECT"):
            f = [int(l[i:i + 5]) for i in range(6, len(l.rstrip()), 5) if l[i:i + 5].strip()]
            bonds.setdefault(f[0], set()).update(f[1:])
            for b in f[1:]:
                bonds.setdefault(b, set()).add(f[0])

    old_parent = bonds.get(h_ser, set())
    if len(old_parent) != 1:
        raise SystemExit(f"{args.proton} should have one bonded partner in CONECT, has {sorted(old_parent)}")
    old_ser = next(iter(old_parent))

    xyz = lambda i: np.array([float(lines[i][30:38]), float(lines[i][38:46]), float(lines[i][46:54])])
    n_xyz = xyz(n_line)
    neigh = [s for s in bonds.get(n_ser, set()) if s != h_ser]
    if not neigh:
        raise SystemExit(f"{args.to} has no bonded neighbours in CONECT")
    v = sum((xyz(serial[s]) - n_xyz) / np.linalg.norm(xyz(serial[s]) - n_xyz) for s in neigh)
    direction = -v / np.linalg.norm(v)
    elem = lines[n_line][76:78].strip() or args.to[0]
    new_h = n_xyz + BOND.get(elem, 1.0) * direction
    old_h = xyz(h_line)
    lines[h_line] = lines[h_line][:30] + "{:8.3f}{:8.3f}{:8.3f}".format(*new_h) + lines[h_line][54:]

    bonds[h_ser] = {n_ser}
    bonds[old_ser].discard(h_ser)
    bonds.setdefault(n_ser, set()).add(h_ser)

    # Closest non-bonded contacts of the moved proton
    contacts = []
    for i in atom_idx:
        s = int(lines[i][6:11])
        if s in (h_ser, n_ser):
            continue
        d = np.linalg.norm(xyz(i) - new_h)
        if d < 2.0:
            contacts.append((d, lines[i][12:27].strip()))
    print(f"Moved {args.proton}: {lines[serial[old_ser]][12:16].strip()} -> {args.to}, "
          f"displacement {np.linalg.norm(new_h - old_h):.2f} A, new {args.to}-H {BOND.get(elem, 1.0):.2f} A")
    print("Contacts < 2.0 A from new H:", ", ".join(f"{n} {d:.2f}" for d, n in sorted(contacts)) or "none")

    out = [l for l in lines if not l.startswith(("CONECT", "END"))]
    for s in sorted(bonds):
        if bonds[s] and s in serial:
            nb = sorted(bonds[s])
            for k in range(0, len(nb), 4):
                out.append("CONECT" + f"{s:5d}" + "".join(f"{b:5d}" for b in nb[k:k + 4]))
    out.append("END")
    open(args.out, "w").write("\n".join(out) + "\n")
    print(f"Wrote {args.out}")

    if args.extract:
        res = [i for i in atom_idx if lines[i][17:20].strip() == args.resname]
        renum = {int(lines[i][6:11]): k + 1 for k, i in enumerate(res)}
        lig = [f"{lines[i][:6]}{renum[int(lines[i][6:11])]:5d}{lines[i][11:]}" for i in res]
        for s in sorted(renum):
            nb = sorted(b for b in bonds.get(s, ()) if b in renum)
            if nb:
                lig.append("CONECT" + f"{renum[s]:5d}" + "".join(f"{renum[b]:5d}" for b in nb))
        lig.append("END")
        open(args.extract, "w").write("\n".join(lig) + "\n")
        print(f"Wrote {args.extract} ({len(res)} atoms)")


if __name__ == "__main__":
    main()
