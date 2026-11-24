#!/usr/bin/env python3
"""What moved, by residue, between a NEB image and its own frozen-source endpoint.

    python mep_shift.py knarr_MEP_FULL.xyz REF.pdb [--actatoms FILE] [--qm FILE]
                        [--tol 0.1] [--show 12]

D331 measured the MM term of `step1`'s path relaxing 18.6-21.2 kcal/mol below
both endpoints, with both endpoints verified to 0.06 kcal/mol of their own
converged optimisations (D330). That says the images find arrangements of the
665-atom active region that neither endpoint optimisation reached. **This says
WHICH atoms.**

**The comparison has to be chosen carefully or it measures the wrong thing.**
ASH's ActiveRegion NEB rebuilds every image's FULL coordinates from the nearer
endpoint's frozen environment and switches at the midpoint (D309), so for 14
images, images 1-6 carry the REACTANT's frozen atoms and images 7-12 carry the
PRODUCT's. Differencing image 12 against image 0 therefore mixes two things:
the active region's relaxation, which is what we want, and the 653-atom
endpoint-to-endpoint difference in the frozen part, which is not. So each image
is compared against its OWN frozen source - 1-6 against image 0, 7-12 against
image 13 - and then every atom that differs is an active-region atom by
construction.

Prints, per image: how many atoms moved past `--tol`, the largest, and the top
residues by maximum displacement, split into water / protein / ligand+metal and
marked QM where a QM atom list is given. The point is to tell a water network
reorganising from a side chain or a loop moving, because those have different
consequences and different literature.
"""
from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict


def read_pdb_labels(path):
    """(chain, resseq, resname, atomname) per ATOM/HETATM, in file order."""
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                out.append((line[21], line[22:26].strip(),
                            line[17:20].strip(), line[12:16].strip()))
    return out


def read_xyz_frames(path):
    """[(comment, [(x,y,z), ...]), ...] for a multi-frame xyz."""
    frames = []
    with open(path) as fh:
        while True:
            head = fh.readline()
            if not head.strip():
                break
            n = int(head.split()[0])
            comment = fh.readline().rstrip()
            coords = []
            for _ in range(n):
                parts = fh.readline().split()
                coords.append((float(parts[1]), float(parts[2]), float(parts[3])))
            frames.append((comment, coords))
    return frames


def read_intlist(path):
    with open(path) as fh:
        return {int(x) for x in fh.read().split()}


def category(resname):
    if resname in ("HOH", "WAT", "TIP3"):
        return "water"
    if resname in ("ZN", "MER", "OH", "mer"):
        return "ligand/metal"
    return "protein"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("mep")
    p.add_argument("pdb", help="a PDB with the SAME atom order (either endpoint)")
    p.add_argument("--actatoms", help="the run's `active_atoms` file (0-based)")
    p.add_argument("--qm", help="file of QM atom indices, to mark them")
    p.add_argument("--tol", type=float, default=0.1,
                   help="displacement in Angstrom to count as moved (default 0.1)")
    p.add_argument("--show", type=int, default=12, help="residues to list per image")
    args = p.parse_args(argv)

    labels = read_pdb_labels(args.pdb)
    frames = read_xyz_frames(args.mep)
    n_img = len(frames)
    if not frames:
        raise SystemExit("no frames in %s" % args.mep)
    if len(frames[0][1]) != len(labels):
        raise SystemExit("atom count mismatch: %s has %d, %s has %d"
                         % (args.mep, len(frames[0][1]), args.pdb, len(labels)))

    active = read_intlist(args.actatoms) if args.actatoms else None
    qm = read_intlist(args.qm) if args.qm else set()

    last = n_img - 1
    print("%d images, %d atoms. Each interior image is compared against its own"
          % (n_img, len(labels)))
    print("frozen source: images 1-%d against image 0, images %d-%d against image %d."
          % (last // 2, last // 2 + 1, last - 1, last))
    if active is not None:
        print("Active region: %d atoms." % len(active))
    print()

    for img in range(1, last):
        ref = 0 if img <= last // 2 else last
        a, b = frames[ref][1], frames[img][1]
        moved = []
        for i, (p1, p2) in enumerate(zip(a, b)):
            d = math.dist(p1, p2)
            if d > args.tol:
                moved.append((d, i))
        if not moved:
            print("image %-2d vs %-2d   nothing moved past %.3g A" % (img, ref, args.tol))
            continue
        outside = [i for _, i in moved if active is not None and i not in active]
        per_res = defaultdict(float)
        per_cat = defaultdict(float)
        for d, i in moved:
            ch, seq, rn, _ = labels[i]
            per_res[(ch, seq, rn)] = max(per_res[(ch, seq, rn)], d)
            per_cat[category(rn)] = max(per_cat[category(rn)], d)
        top = sorted(per_res.items(), key=lambda kv: -kv[1])[:args.show]
        print("image %-2d vs %-2d   %d atoms moved > %.3g A, largest %.3f, "
              "%d residues%s"
              % (img, ref, len(moved), args.tol, max(d for d, _ in moved),
                 len(per_res),
                 ", %d OUTSIDE the active region" % len(outside) if outside else ""))
        cats = ", ".join("%s %.2f" % (k, v) for k, v in sorted(per_cat.items()))
        print("      largest by class:  %s" % cats)
        for (ch, seq, rn), d in top:
            tag = "  [QM]" if any(labels[i][:3] == (ch, seq, rn) and i in qm
                                  for _, i in moved) else ""
            print("      %-14s %.3f%s" % ("%s:%s:%s" % (ch, seq, rn), d, tag))
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
