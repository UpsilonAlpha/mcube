#!/usr/bin/env python3
"""
Crop an AF2 model of the channel and orient it for mcube's membrane builder, which never
moves the protein: the barrel axis on z and the hydrophobic belt centred on z = 0.

  - crop: drop residues by the channel's own numbering (default 110-130, the N-terminal
    segment that sits in the lumen) and renumber the rest to that numbering
  - axis: the β-strand directions, sign-aligned and averaged (as in pore.py)
  - sign: both termini of an outer-membrane β-barrel are periplasmic, so the first strand's
    start goes to z < 0 (periplasm below, extracellular side above)
  - centre (--centre charge-gap, the default): z = 0 is the middle of the longest stretch of
    the barrel wall with no charged lipid-facing residue (Lys/Arg/Asp/Glu with Cβ further from
    the axis than Cα and Cα at least 9 Å out). Charged side chains snorkel to the interfaces,
    so that stretch is the hydrophobic belt. (--centre kd: the --slab-thick slab with the
    highest mean Kyte-Doolittle score of lipid-facing strand residues. On the channel model it
    put the belt 6 Å too low, pulled down by polar residues near the periplasmic turns.)
    The lipid-facing Tyr/Trp/Phe are reported: the Tyr/Trp girdles should sit near the edges.

  python orient.py --pdb ranked_0.pdb --crop 110-130 --out channel_cropped_oriented.pdb
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

KD = {"ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5, "GLN": -3.5, "GLU": -3.5,
      "GLY": -0.4, "HIS": -3.2, "ILE": 4.5, "LEU": 3.8, "LYS": -3.9, "MET": 1.9, "PHE": 2.8,
      "PRO": -1.6, "SER": -0.8, "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2}


def strands(dssp, min_len=5):
    out, start = [], None
    for i, s in enumerate(list(dssp) + ["-"]):
        if s == "E" and start is None:
            start = i
        elif s != "E" and start is not None:
            if i - start >= min_len:
                out.append((start, i - 1))
            start = None
    return out


def rotation_to_z(axis):
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(axis, z)
    c = float(axis @ z)
    if np.linalg.norm(v) < 1e-8:
        return np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1 / (1 + c))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pdb", required=True)
    p.add_argument("--first-residue", type=int, default=110, help="channel number of residue 1")
    p.add_argument("--crop", default="110-130", help="channel residues to remove, e.g. 110-130 ('' for none)")
    p.add_argument("--slab-thick", type=float, default=24.0, help="hydrophobic slab for --centre kd (Å)")
    p.add_argument("--centre", choices=["charge-gap", "kd"], default="charge-gap")
    p.add_argument("--number-from", type=int,
                   help="renumber the written residues from this (mcube wants 1..N per chain); "
                        "the report keeps the channel numbering and gives the offset")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    import mdtraj
    traj = mdtraj.load(args.pdb)
    top = traj.topology
    # channel numbering first, so crop and output share it
    for res in top.residues:
        res.resSeq = res.resSeq + args.first_residue - 1
    if args.crop:
        lo, hi = (int(v) for v in args.crop.split("-"))
        keep = [a.index for a in top.atoms if not lo <= a.residue.resSeq <= hi]
        traj = traj.atom_slice(keep)
        top = traj.topology

    xyz = traj.xyz[0] * 10.0
    dssp = mdtraj.compute_dssp(traj, simplified=True)[0]
    runs = strands(dssp)
    ca = {r.index: next(a.index for a in r.atoms if a.name == "CA") for r in top.residues}
    cb = {r.index: a.index for r in top.residues for a in r.atoms if a.name == "CB"}
    ref = xyz[ca[runs[0][1]]] - xyz[ca[runs[0][0]]]
    dirs = [(lambda d: d if d @ ref > 0 else -d)(xyz[ca[b]] - xyz[ca[a]]) for a, b in runs]
    axis = np.mean(dirs, axis=0)
    axis /= np.linalg.norm(axis)
    centre = xyz[[ca[i] for a, b in runs for i in range(a, b + 1)]].mean(0)
    # the first strand runs from the periplasm outward: make its direction +z
    if ref @ axis < 0:
        axis = -axis
    rot = rotation_to_z(axis)
    new = (xyz - centre) @ rot.T

    # lipid-facing strand residues and their z
    in_strand = {i for a, b in runs for i in range(a, b + 1)}
    radial = lambda p: np.linalg.norm(p[:2])
    facing = []
    for r in top.residues:
        if r.index in in_strand and r.index in cb and radial(new[cb[r.index]]) > radial(new[ca[r.index]]):
            facing.append((float(new[ca[r.index]][2]), r.name, r.resSeq))
    zs = np.array([f[0] for f in facing])
    kd = np.array([KD.get(f[1], 0.0) for f in facing])
    best_kd, best_score = 0.0, -1e9
    for zc in np.arange(zs.min() + args.slab_thick / 2, zs.max() - args.slab_thick / 2 + 0.01, 0.5):
        inside = np.abs(zs - zc) <= args.slab_thick / 2
        if inside.sum() < 10:
            continue
        score = kd[inside].mean()
        if score > best_score:
            best_kd, best_score = float(zc), float(score)
    # charged lipid-facing residues anywhere on the wall (not only in strands)
    charged = sorted(float(new[ca[r.index]][2]) for r in top.residues
                     if r.name in ("LYS", "ARG", "ASP", "GLU") and r.index in cb
                     and radial(new[ca[r.index]]) >= 9.0
                     and radial(new[cb[r.index]]) > radial(new[ca[r.index]]))
    gaps = [(b - a, a, b) for a, b in zip(charged, charged[1:])]
    gap_len, gap_lo, gap_hi = max(gaps)
    best_gap = (gap_lo + gap_hi) / 2
    best = best_gap if args.centre == "charge-gap" else best_kd
    new[:, 2] -= best

    arom = sorted((round(z - best, 1), f"{n}{s}") for z, n, s in facing if n in ("TYR", "TRP", "PHE"))
    traj.xyz[0] = new / 10.0
    traj.unitcell_vectors = None
    first = min(r.resSeq for r in top.residues)
    if args.number_from is not None:
        for k, r in enumerate(top.residues):
            r.resSeq = args.number_from + k
    traj.save_pdb(args.out)
    report = {"pdb": str(Path(args.pdb).resolve()), "out": str(Path(args.out).resolve()),
              "cropped": args.crop, "residues": top.n_residues,
              "first_residue": first,
              "written_numbering_offset": (first - args.number_from) if args.number_from is not None else 0,
              "strands": len(runs), "centre_method": args.centre, "belt_shift_A": round(best, 1),
              "charge_free_span_A": [round(gap_lo - best, 1), round(gap_hi - best, 1)],
              "charge_free_length_A": round(gap_len, 1),
              "kd_slab_centre_A": round(best_kd - best, 1), "belt_mean_KD": round(best_score, 2),
              "lipid_facing_aromatics_z": arom,
              "z_extent_A": [round(float(new[:, 2].min()), 1), round(float(new[:, 2].max()), 1)]}
    Path(args.out).with_suffix(".json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
