#!/usr/bin/env python3
"""WHERE, radially, did a relaxation happen? Two PDBs, same atom order.

    python shell_shift.py REF.pdb MOVED.pdb [--anchor B:1:ZN] [--tol 0.1]
                          [--show 20] [--shells 3,4.5,6,8]

¶341 measured `saddle1` draining 18.41 kcal/mol out of the MM term at a FIXED
reaction coordinate, and ¶332 asks whether the published barrier should be
scored in a common FROZEN environment. Those two questions meet at one number
this project has never measured: **how far from the metal that relaxation
lives.**

It decides ¶332 and the two answers are opposite:

  far field   Ryde (JCTC 2017, 13, 5745) prices a far-from-site hydrogen bond at
              ~20 kJ/mol and calls the resulting energy change "most likely
              irrelevant for the reaction". Calixto (JCTC 2017, 13, 5486)
              prices freezing everything outside a 6.00 A free region at 0.8
              kcal/mol average, 1.6 maximum. If the relaxation is out there,
              freezing removes a large error for a small one and ¶332 settles
              as "freeze".
  first shell ¶333 measured the NEB's movers and every one was first or second
              shell - the ligand, the hydroxide, Asp92, Ser186, His88, His164,
              Asn225, the zinc itself, and six ordered waters. Energy carried
              THERE is the reaction's own environmental response, and freezing
              it away discards chemistry rather than noise. ¶332 then cannot
              settle as "freeze" at any radius.

So this reports displacement per radial shell rather than per image. The shell
edges default to 3 / 4.5 / 6 / 8 A because each one is a threshold someone has
already argued for: 3 A is the user's standing restraint rule ("the only thing
that needs to be restrained are coordinating atoms in the metal's direct
environment - roughly 3 angstroms"), 4.5 A is the first solvation shell past the
ligands, 6.00 A is Calixto's measured free-region threshold, and 8 A is this
project's `actradius`, outside which the atoms are frozen by construction and
must show zero.

`mep_shift.py` answers "which residues moved" for a multi-frame NEB path. This
answers "at what radius" for a single before/after pair, which is the form every
optimisation handoff in this project takes.
"""
from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict


def read_pdb(path):
    """[(chain, resseq, resname, atomname, x, y, z), ...] in file order."""
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                out.append((
                    line[21], line[22:26].strip(), line[17:20].strip(),
                    line[12:16].strip(),
                    float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return out


def category(resname):
    if resname in ("HOH", "WAT", "TIP3"):
        return "water"
    if resname in ("ZN", "MER", "OH", "mer"):
        return "ligand/metal"
    return "protein"


def find_anchor(atoms, spec):
    """`B:1:ZN` -> that atom's coordinates, or the residue's first atom."""
    want = spec.split(":")
    for a in atoms:
        if [a[0], a[1], a[2]][:len(want)] == want[:3] or \
           (len(want) == 3 and (a[0], a[1], a[3]) == tuple(want)):
            return a[4:7]
    raise SystemExit("anchor %s not found" % spec)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("ref", help="the structure the other one started from")
    p.add_argument("moved", help="the structure it relaxed to, SAME atom order")
    p.add_argument("--anchor", default="B:1:ZN",
                   help="chain:resseq:name or chain:resseq:resname (default B:1:ZN)")
    p.add_argument("--tol", type=float, default=0.1,
                   help="displacement in Angstrom to count as moved (default 0.1)")
    p.add_argument("--show", type=int, default=20, help="residues to list")
    p.add_argument("--shells", default="3,4.5,6,8",
                   help="comma-separated shell edges in Angstrom (default 3,4.5,6,8)")
    args = p.parse_args(argv)

    a = read_pdb(args.ref)
    b = read_pdb(args.moved)
    if len(a) != len(b):
        raise SystemExit("atom count mismatch: %d vs %d" % (len(a), len(b)))

    anchor = find_anchor(a, args.anchor)
    edges = [0.0] + [float(x) for x in args.shells.split(",")] + [float("inf")]
    rows = []
    for x, y in zip(a, b):
        rows.append((math.dist(x[4:7], y[4:7]), math.dist(x[4:7], anchor), x))

    print("%d atoms, anchor %s, tol %.3g A" % (len(a), args.anchor, args.tol))
    print("%s\n  vs\n%s\n" % (args.ref, args.moved))
    print("%-12s %8s %8s %8s %9s %10s" %
          ("shell/A", "atoms", "moved", "max/A", "mean/A", "sum d^2"))
    for lo, hi in zip(edges, edges[1:]):
        g = [r for r in rows if lo <= r[1] < hi]
        if not g:
            continue
        mv = [r for r in g if r[0] > args.tol]
        label = "%.3g-%.3g" % (lo, hi) if hi != float("inf") else ">%.3g" % lo
        print("%-12s %8d %8d %8.3f %9.4f %10.2f" %
              (label, len(g), len(mv), max(r[0] for r in g),
               sum(r[0] for r in g) / len(g), sum(r[0] ** 2 for r in g)))

    mv = [r for r in rows if r[0] > args.tol]
    print("\n%d of %d atoms moved past %.3g A, largest %.3f A"
          % (len(mv), len(rows), args.tol, max([r[0] for r in rows], default=0.0)))
    per_cat = defaultdict(float)
    for d, dz, x in mv:
        per_cat[category(x[2])] = max(per_cat[category(x[2])], d)
    if per_cat:
        print("largest by class:  %s"
              % ", ".join("%s %.2f" % kv for kv in sorted(per_cat.items())))

    per_res = {}
    for d, dz, x in mv:
        k = (x[0], x[1], x[2])
        if d > per_res.get(k, (0.0, 0.0))[0]:
            per_res[k] = (d, dz)
    print("\n%-16s %8s %9s  %s" % ("residue", "max/A", "anchor/A", "class"))
    for k, v in sorted(per_res.items(), key=lambda kv: -kv[1][0])[:args.show]:
        print("%-16s %8.3f %9.2f  %s"
              % ("%s:%s:%s" % k, v[0], v[1], category(k[2])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
