#!/usr/bin/env python3
"""
Polar contacts between chosen residues and a target group over a trajectory,
with occupancy: e.g. which Lys/Arg/Asp atoms touch the substrate, the
hydroxide or the metal-bound water, and how often.

For every pair (residue atom, target atom) with both atoms polar (N, O, S by
default) that comes within --cutoff Å in any frame, it reports the fraction of
frames within the cutoff and the min/mean distance. --frames restricts the
analysis to listed frames (e.g. the near-attack frames), so subsets can be
compared. --near lists every residue of a kind within a distance of the target
(e.g. all arginines within 10 Å) before the contact table.

Usage:
    python contact_report.py topology.pdb traj.dcd --target "resname mer or resname OH" \
        --sel "resname ARG LYS ASP GLU GLN HIS" [--cutoff 3.5] [--near 10] \
        [--frames frames.txt] [--stride 1] [--out contacts.csv]
"""

import argparse

import mdtraj
import numpy as np

POLAR = {"N", "O", "S"}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("topology")
    p.add_argument("trajectory")
    p.add_argument("--target", required=True, help="mdtraj selection of the target group")
    p.add_argument("--sel", required=True, help="mdtraj selection of residues to test")
    p.add_argument("--cutoff", type=float, default=3.5, help="Contact distance (Å)")
    p.add_argument("--near", type=float, default=0.0, help="Also list --sel residues within this distance (Å)")
    p.add_argument("--all-elements", action="store_true", help="Include C/H atoms, not only N/O/S")
    p.add_argument("--frames", help="File of 0-based frame indices (one per line or whitespace separated)")
    p.add_argument("--stride", type=int, default=1)
    p.add_argument("--min-occupancy", type=float, default=0.05)
    p.add_argument("--out", default="contacts.csv")
    args = p.parse_args()

    traj = mdtraj.load(args.trajectory, top=args.topology, stride=args.stride)
    if args.frames:
        idx = [int(x) // args.stride for x in open(args.frames).read().split()]
        traj = traj[idx]
    top = traj.topology
    keep = (lambda a: True) if args.all_elements else (lambda a: a.element is not None and a.element.symbol in POLAR)
    target = [i for i in top.select(args.target) if keep(top.atom(i))]
    tset = set(top.select(args.target))
    cand = [i for i in top.select(args.sel) if keep(top.atom(i)) and i not in tset]
    if not target or not cand:
        raise SystemExit("empty target or selection")
    lab = lambda i: f"{top.atom(i).residue.name}{top.atom(i).residue.resSeq}:{top.atom(i).name}"  # noqa: E731
    print(f"{traj.n_frames} frames; target {len(target)} atoms ({args.target}); {len(cand)} candidate atoms")

    if args.near:
        tcom = mdtraj.compute_center_of_mass(traj.atom_slice(sorted(tset)))
        res = {}
        for a in cand:
            r = top.atom(a).residue
            d = np.linalg.norm(traj.xyz[:, a] - tcom, axis=1) * 10.0
            res.setdefault(r, []).append(d.mean())
        near = sorted((min(v), f"{r.name}{r.resSeq}") for r, v in res.items() if min(v) <= args.near)
        print(f"Residues with a polar atom within {args.near} Å of the target centre (mean distance): "
              + ", ".join(f"{n} {d:.1f}" for d, n in near))

    pairs = [(a, t) for a in cand for t in target]
    d = mdtraj.compute_distances(traj, pairs) * 10.0
    within = d < args.cutoff
    occ = within.mean(0)
    order = np.argsort(-occ)
    rows = [(pairs[k], occ[k], d[:, k].min(), d[:, k].mean()) for k in order if occ[k] >= args.min_occupancy]
    with open(args.out, "w") as f:
        f.write("residue_atom,target_atom,occupancy,min_A,mean_A\n")
        for (a, t), o, mn, me in rows:
            f.write(f"{lab(a)},{lab(t)},{o:.3f},{mn:.2f},{me:.2f}\n")
    print(f"Contacts within {args.cutoff} Å (occupancy >= {args.min_occupancy:.0%}):")
    for (a, t), o, mn, me in rows[:40]:
        print(f"  {lab(a):>16} .. {lab(t):<14} {o:6.1%}   min {mn:.2f}  mean {me:.2f} Å")
    if not rows:
        print("  none")


if __name__ == "__main__":
    main()
