#!/usr/bin/env python3
"""
Does an AF2 model of the channel have a continuous pore, and where is it narrowest?

  - β-strands from DSSP (runs of >= 5 'E'), their count, and whether the first and last
    strand pair up (a closed barrel).
  - The barrel axis: the strands' directions, flipped to agree with the first strand and
    averaged; the centre is the mean strand Cα.
  - A HOLE-like radius profile along that axis: in each 1 Å slice the largest sphere that
    fits between the heavy atoms (van der Waals radii), its centre searched on an in-plane
    grid. The narrowest slice is the constriction; the residues lining it are listed.
  - pLDDT from the B-factor column, and pTM from the AF2 result pickle if given.

  python pore.py --pdb ranked_0.pdb [--result result_model_3_ptm_pred_0.pkl] --out pore.json
"""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}


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


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pdb", required=True)
    p.add_argument("--result", help="AF2 result_*.pkl, for pTM / ipTM")
    p.add_argument("--first-residue", type=int, default=110, help="sequence number of residue 1")
    p.add_argument("--half-length", type=float, default=30.0, help="profile from -L to +L (Å)")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    import mdtraj
    traj = mdtraj.load(args.pdb)
    top = traj.topology
    xyz = traj.xyz[0] * 10.0
    ca = np.array([a.index for a in top.atoms if a.name == "CA"])
    plddt = np.array([float(l[60:66]) for l in Path(args.pdb).read_text().splitlines()
                      if l.startswith("ATOM") and l[12:16].strip() == "CA"])

    dssp = mdtraj.compute_dssp(traj, simplified=True)[0]
    runs = strands(dssp)
    ref = xyz[ca[runs[0][1]]] - xyz[ca[runs[0][0]]]
    dirs = []
    for a, b in runs:
        d = xyz[ca[b]] - xyz[ca[a]]
        dirs.append(d if d @ ref > 0 else -d)
    axis = np.mean(dirs, axis=0)
    axis /= np.linalg.norm(axis)
    strand_ca = np.concatenate([ca[a:b + 1] for a, b in runs])
    centre = xyz[strand_ca].mean(axis=0)

    # closed barrel: does any residue of the last strand sit within H-bonding Cα distance
    # (< 5.5 Å) of one in the first?
    first = xyz[ca[runs[0][0]:runs[0][1] + 1]]
    last = xyz[ca[runs[-1][0]:runs[-1][1] + 1]]
    closure = float(np.min(np.linalg.norm(first[:, None] - last[None], axis=2)))

    heavy = np.array([a.index for a in top.atoms if a.element.symbol != "H"])
    radii = np.array([VDW.get(top.atom(i).element.symbol, 1.8) for i in heavy])
    pos = xyz[heavy] - centre
    z_atoms = pos @ axis
    u = np.cross(axis, [1.0, 0, 0] if abs(axis[0]) < 0.9 else [0, 1.0, 0])
    u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    # barrel radius: strand Cα distance from the axis; the probe stays inside 80 % of it,
    # or in a protein with no membrane around it the widest "pore" is the open space outside
    sc = xyz[strand_ca] - centre
    barrel_r = float(np.linalg.norm(sc - np.outer(sc @ axis, axis), axis=1).mean())
    grid = np.arange(-barrel_r, barrel_r + 0.01, 0.5)
    gu, gv = np.meshgrid(grid, grid)
    plane = np.stack([gu.ravel(), gv.ravel()], axis=1)
    plane = plane[np.linalg.norm(plane, axis=1) < 0.8 * barrel_r]

    profile, prev = [], np.zeros(2)
    for z in np.arange(-args.half_length, args.half_length + 0.01, 1.0):
        near = np.abs(z_atoms - z) < 8.0
        if not near.any():
            continue
        # search near the previous slice's centre so the path stays connected
        cand = plane[np.linalg.norm(plane - prev, axis=1) < 6.0]
        pts = z * axis + cand[:, :1] * u + cand[:, 1:] * v
        d = np.linalg.norm(pts[:, None, :] - pos[near][None], axis=2) - radii[near][None]
        free = d.min(axis=1)
        k = int(np.argmax(free))
        prev = cand[k]
        profile.append({"z": round(float(z), 1), "radius": round(float(free[k]), 2),
                        "centre_offset": round(float(np.linalg.norm(cand[k])), 2),
                        "_pt": pts[k]})

    inside = [s for s in profile if abs(s["z"]) <= 12.0]       # the membrane-spanning part
    narrow = min(inside, key=lambda s: s["radius"])
    d = np.linalg.norm(pos - narrow["_pt"], axis=1) - radii
    lining = sorted({(top.atom(heavy[i]).residue.resSeq + args.first_residue - 1,
                      top.atom(heavy[i]).residue.name) for i in np.where(d < 2.0)[0]})
    for s in profile:
        s.pop("_pt")

    report = {"pdb": str(Path(args.pdb).resolve()), "residues": int(len(ca)),
              "mean_plddt": round(float(plddt.mean()), 1),
              "plddt_ge90": int((plddt >= 90).sum()), "plddt_lt70": int((plddt < 70).sum()),
              "dssp": "".join(dssp),
              "strands": len(runs),
              "strand_ranges": [[a + args.first_residue, b + args.first_residue] for a, b in runs],
              "first_last_strand_min_ca_A": round(closure, 2),
              "barrel_closed": closure < 5.5,
              "barrel_radius_ca_A": round(barrel_r, 1),
              "constriction": {"z": narrow["z"], "radius_A": narrow["radius"],
                               "lining": [f"{n}{r}" for r, n in lining]},
              "min_radius_membrane_A": narrow["radius"],
              "profile": profile}
    if args.result:
        with open(args.result, "rb") as fh:
            res = pickle.load(fh)
        for key in ("ptm", "iptm", "ranking_confidence"):
            if key in res:
                report[key] = round(float(np.asarray(res[key])), 3)
    Path(args.out).write_text(json.dumps(report, indent=1))
    print(json.dumps({k: v for k, v in report.items() if k != "profile"}, indent=1))


if __name__ == "__main__":
    main()
