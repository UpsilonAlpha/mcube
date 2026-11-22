#!/usr/bin/env python
"""Do two NEB endpoints differ ONLY inside the active region?

ASH's NEB with `ActiveRegion=True` does not freeze the inactive atoms at one
geometry. `KnarrCalculator.Compute` rebuilds each image's full coordinates from
whichever endpoint the image index is nearest - the reactant's for the first
half, the product's for the second - so if the two endpoints' inactive atoms sit
anywhere different, the frozen environment SWITCHES at the midpoint of the path.
That is a step change in the energy, placed where the climbing image goes, and
nothing in the output says so. (The comment in ASH's own `__init__` says
"Inactive part of reactant will be used for all images"; the code does not do
that.)

So the condition a NEB needs is stronger than "the same active region in both
runs": the two endpoints must carry the SAME COORDINATES for every atom outside
it. Then the midpoint switch is a no-op and the discontinuity does not exist.

A product written by `scan` from the reactant satisfies this by construction, as
long as both runs resolved the same active region - the scan never moved an atom
outside it. This script measures that instead of assuming it, because the radius
is re-derived per structure (`resolve_actatoms`) and a structure that drifts can
bring a different residue inside it.

    python frozen_check.py REACTANT.pdb PRODUCT.pdb [--actatoms FILE] [--tol 1e-3]

`--actatoms` is the `active_atoms` file either run wrote (`results.actatoms_file`);
indices are 0-based, as ASH writes and reads them. Without it, the script still
reports every atom that moved, which says what the active region would have to
contain for the path to be clean.
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import OrderedDict


def read_pdb(path):
    """(key, xyz) per ATOM/HETATM record, in file order - the order ASH indexes."""
    keys, xyz = [], []
    with open(path) as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                keys.append((line[21], line[22:27].strip(), line[17:20].strip(),
                             line[12:16].strip()))
                xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return keys, xyz


def read_actatoms(path):
    with open(path) as handle:
        return set(int(tok) for tok in handle.read().split())


def label(key):
    chain, resseq, resname, atom = key
    return f"{chain}:{resseq}:{resname}:{atom}"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("reactant")
    parser.add_argument("product")
    parser.add_argument("--actatoms", help="the `active_atoms` file a run wrote (0-based)")
    parser.add_argument("--tol", type=float, default=1e-3,
                        help="displacement counted as motion, A (default 1e-3: "
                             "PDB coordinates carry 3 decimals, so anything "
                             "smaller is rounding)")
    parser.add_argument("--show", type=int, default=12, help="moved atoms to list")
    args = parser.parse_args(argv)

    left_keys, left_xyz = read_pdb(args.reactant)
    right_keys, right_xyz = read_pdb(args.product)

    print(f"reactant {args.reactant}: {len(left_keys)} atoms")
    print(f"product  {args.product}: {len(right_keys)} atoms")
    if len(left_keys) != len(right_keys):
        print(f"[FAIL] different atom counts: NEB needs one system in two "
              f"geometries, so this pair cannot be a path at all")
        return 2
    bad = [i for i, (a, b) in enumerate(zip(left_keys, right_keys)) if a != b]
    if bad:
        print(f"[FAIL] same count, different atom order ({len(bad)} positions). "
              f"First: index {bad[0]}, reactant {label(left_keys[bad[0]])}, "
              f"product {label(right_keys[bad[0]])}")
        return 2
    print(f"[OK]   same atoms in the same order")

    disp = [math.dist(a, b) for a, b in zip(left_xyz, right_xyz)]
    moved = [i for i, d in enumerate(disp) if d > args.tol]
    print(f"[INFO] {len(moved)} of {len(disp)} atoms moved by more than "
          f"{args.tol} A; largest displacement {max(disp):.3f} A")

    if args.actatoms is None:
        print(chr(10) + "No --actatoms given, so the verdict is not computable. "
              "What moved, by residue:")
        residues = OrderedDict()
        for i in moved:
            chain, resseq, resname, _ = left_keys[i]
            residues.setdefault(f"{chain}:{resseq}:{resname}", []).append(disp[i])
        for name, values in list(residues.items())[:40]:
            print(f"    {name:18s} {len(values):3d} atoms, max {max(values):.3f} A")
        if len(residues) > 40:
            print(f"    ... and {len(residues) - 40} more residues")
        return 1

    active = read_actatoms(args.actatoms)
    print(f"[INFO] active region {args.actatoms}: {len(active)} atoms")
    outside = sorted(i for i in moved if i not in active)
    inside_max = max([disp[i] for i in range(len(disp)) if i in active], default=0.0)
    print(f"[INFO] largest displacement inside the active region: {inside_max:.3f} A")

    if not outside:
        print(chr(10) + "[PASS] every atom that moved is in the active region, so the "
              "two endpoints carry identical coordinates for every frozen atom. "
              "ASH's midpoint switch from the reactant's inactive geometry to the "
              "product's is a no-op, and the path has no discontinuity from it.")
        return 0

    print(chr(10) + f"[FAIL] {len(outside)} atoms OUTSIDE the active region differ "
          f"between the endpoints, by up to "
          f"{max(disp[i] for i in outside):.3f} A.")
    print("These are frozen during the NEB, and ASH takes them from the reactant "
          "for the first half of the path and from the product for the second, so "
          "the frozen environment changes under the climbing image. The barrier "
          "would include that step. Either widen the active region to contain "
          "them, or rebuild the product from the reactant so nothing outside the "
          "region ever moves.")
    for i in outside[:args.show]:
        print(f"    index {i:6d}  {label(left_keys[i]):24s} {disp[i]:.3f} A")
    if len(outside) > args.show:
        print(f"    ... and {len(outside) - args.show} more")
    return 3


if __name__ == "__main__":
    sys.exit(main())
