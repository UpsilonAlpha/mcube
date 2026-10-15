#!/usr/bin/env python3
"""
Metal-site and stability report for an MD trajectory (any system).

- Metal coordination shell: every heavy atom within --cutoff of the metal in any
  frame, with the fraction of frames it is bound and its min/mean/max distance.
- Extra distances (--pair A B) and angles (--angle A B C) per frame.
- Protein CA RMSD (to frame 0) and the RMSD of --rmsd-sel groups after fitting
  on CA, to spot drift of a ligand.

Atoms are given as an mdtraj selection ("resname MER and name O2x"), a
"chain:resid:atom" selector ("B:1:ZN", chain letters as in the PDB) or a
0-based index. Writes site_frames.csv and site_report.txt into --outdir.

Usage:
    python site_report.py topology.pdb traj.dcd --outdir site \
        [--metal "element Zn"] [--cutoff 2.8] [--skip 0] \
        [--pair "C:2:O01" "C:1:C4x"] [--angle "C:2:O01" "C:1:C4x" "C:1:O2x"] \
        [--rmsd-sel "resname MER and not element H"]
"""

import argparse
import re
from pathlib import Path

import mdtraj
import numpy as np

SELECTOR = re.compile(r"^\s*([^:\s]+)\s*:\s*(-?\d+)\s*:\s*([^:\s]+)\s*$")


def chain_letter(chain):
    return getattr(chain, "chain_id", None) or str(chain.index)


def atom_label(top, i):
    a = top.atom(i)
    r = a.residue
    return f"{chain_letter(r.chain)}:{r.resSeq}:{a.name} ({r.name}, #{i})"


def resolve(top, spec):
    spec = str(spec)
    if spec.strip().isdigit():
        return int(spec)
    m = SELECTOR.match(spec)
    if m:
        ch, rid, name = m.group(1), int(m.group(2)), m.group(3)
        hits = [a.index for a in top.atoms
                if a.name == name and a.residue.resSeq == rid and chain_letter(a.residue.chain) == ch]
    else:
        hits = list(top.select(spec))
    if len(hits) != 1:
        raise SystemExit(f"'{spec}' matched {len(hits)} atoms, expected 1")
    return int(hits[0])


def stats(x):
    return f"{np.min(x):.2f} / {np.mean(x):.2f} / {np.max(x):.2f}"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("topology")
    p.add_argument("trajectory")
    p.add_argument("--outdir", default="site")
    p.add_argument("--metal", default="element Zn", help="Metal atom (one atom)")
    p.add_argument("--cutoff", type=float, default=2.8, help="Coordination cutoff (A)")
    p.add_argument("--pair", nargs=2, action="append", default=[], metavar=("A", "B"))
    p.add_argument("--angle", nargs=3, action="append", default=[], metavar=("A", "B", "C"))
    p.add_argument("--rmsd-sel", action="append", default=[], help="Group RMSD after CA fit")
    p.add_argument("--skip", type=int, default=0, help="Ignore the first N frames")
    args = p.parse_args()

    traj = mdtraj.load(args.trajectory, top=args.topology)[args.skip:]
    top = traj.topology
    n = traj.n_frames
    metal = resolve(top, args.metal)
    lines = [f"{n} frames from {args.trajectory} (skip {args.skip}); metal {atom_label(top, metal)}",
             f"Coordination shell (heavy atoms within {args.cutoff} A in any frame): bound fraction, "
             "min / mean / max distance (A)"]

    heavy = [a.index for a in top.atoms if a.element is not None and a.element.symbol != "H" and a.index != metal]
    near = mdtraj.compute_neighbors(traj, args.cutoff / 10.0, [metal], haystack_indices=heavy)
    candidates = sorted(set(int(i) for frame in near for i in frame))
    cols = {}
    if candidates:
        d = mdtraj.compute_distances(traj, [[metal, c] for c in candidates]) * 10.0
        for k, c in enumerate(candidates):
            frac = float(np.mean(d[:, k] <= args.cutoff))
            lines.append(f"  {atom_label(top, c):45s} {frac:6.1%}   {stats(d[:, k])}")
            cols[f"M-{c}"] = d[:, k]
        cn = np.sum(d <= args.cutoff, axis=1)
        lines.append(f"Coordination number: min {cn.min()}, mean {cn.mean():.2f}, max {cn.max()}")
        cols["CN"] = cn

    for a, b in args.pair:
        i, j = resolve(top, a), resolve(top, b)
        x = mdtraj.compute_distances(traj, [[i, j]])[:, 0] * 10.0
        cols[f"d {a}-{b}"] = x
        lines.append(f"Distance {atom_label(top, i)} - {atom_label(top, j)}: {stats(x)} A; "
                     f"first {x[0]:.2f}, last {x[-1]:.2f}")
    for a, b, c in args.angle:
        i, j, k = resolve(top, a), resolve(top, b), resolve(top, c)
        x = np.degrees(mdtraj.compute_angles(traj, [[i, j, k]])[:, 0])
        cols[f"ang {a}-{b}-{c}"] = x
        lines.append(f"Angle {a}-{b}-{c}: {stats(x)} deg; first {x[0]:.1f}, last {x[-1]:.1f}")

    ca = top.select("name CA and protein")
    if len(ca):
        traj.superpose(traj, 0, atom_indices=ca)
        r = np.sqrt(np.mean(np.sum((traj.xyz[:, ca] - traj.xyz[0, ca]) ** 2, axis=2), axis=1)) * 10.0
        cols["CA RMSD"] = r
        lines.append(f"Protein CA RMSD to first frame: {stats(r)} A; last {r[-1]:.2f}")
        for sel in args.rmsd_sel:
            g = top.select(sel)
            if not len(g):
                raise SystemExit(f"--rmsd-sel '{sel}' matched no atoms")
            rg = np.sqrt(np.mean(np.sum((traj.xyz[:, g] - traj.xyz[0, g]) ** 2, axis=2), axis=1)) * 10.0
            cols[f"RMSD {sel}"] = rg
            lines.append(f"RMSD of '{sel}' ({len(g)} atoms) after CA fit: {stats(rg)} A; last {rg[-1]:.2f}")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    names = list(cols)
    with open(out / "site_frames.csv", "w") as f:
        f.write("frame," + ",".join(names) + "\n")
        for t in range(n):
            f.write(f"{t + args.skip}," + ",".join(f"{cols[k][t]:.3f}" for k in names) + "\n")
    text = "\n".join(lines)
    (out / "site_report.txt").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
