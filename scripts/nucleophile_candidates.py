#!/usr/bin/env python3
"""
Candidate nucleophiles for attack on an electrophilic carbon over a trajectory
(e.g. waters near a beta-lactam carbonyl carbon), and whether a base is ready
to take a proton from each.

For every candidate oxygen (all water O atoms, plus any --extra atoms such as a
metal-bound hydroxide) that comes within --cutoff Å of the carbon in any frame:
- distance d(O..C) and attack angle O..C=O (Bürgi-Dunitz ~95-120 deg);
- whether it is bound to the metal (d(O..metal) <= --metal-cutoff);
- H-bonds to each --base acceptor: heavy-atom distance <= --hb-dist and one of
  the candidate's hydrogens within --hb-hdist of the acceptor (H pointing at it);
- "attack frames": d <= --dmax and amin <= angle <= amax.
Summary per candidate (occupancy, minimum d, best angle, attack fraction,
metal-bound fraction, H-bond fractions, and attack frames *with* a base H-bond),
plus per-frame best candidate. Atoms: mdtraj selections, chain:resid:atom or
indices.

Usage:
    python nucleophile_candidates.py top.pdb traj.dcd --carbon C:1:C4x --oxygen C:1:O2x \
        --metal B:1:ZN --base A:92:OD1 --base A:92:OD2 --base A:226:NZ \
        --extra C:2:O01 [--cutoff 5] [--dmax 3.5] [--amin 95 --amax 120] [--outdir nuc]
"""

import argparse
import re
from pathlib import Path

import mdtraj
import numpy as np

SELECTOR = re.compile(r"^\s*([^:\s]+)\s*:\s*(-?\d+)\s*:\s*([^:\s]+)\s*$")


def chain_letter(chain):
    return getattr(chain, "chain_id", None) or str(chain.index)


def resolve(top, spec):
    spec = str(spec)
    if spec.strip().isdigit():
        return int(spec)
    m = SELECTOR.match(spec)
    if m:
        hits = [a.index for a in top.atoms if a.name == m.group(3) and a.residue.resSeq == int(m.group(2))
                and chain_letter(a.residue.chain) == m.group(1)]
    else:
        hits = list(top.select(spec))
    if len(hits) != 1:
        raise SystemExit(f"'{spec}' matched {len(hits)} atoms, expected 1")
    return int(hits[0])


def label(top, i):
    a = top.atom(i)
    return f"{a.residue.name}{chain_letter(a.residue.chain)}{a.residue.resSeq}:{a.name}({i})"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("topology")
    p.add_argument("trajectory")
    p.add_argument("--carbon", required=True)
    p.add_argument("--oxygen", required=True, help="Carbonyl O of the electrophile (for the attack angle)")
    p.add_argument("--metal", help="Metal atom (for metal-bound flags)")
    p.add_argument("--base", action="append", default=[], help="Base acceptor atom (repeatable)")
    p.add_argument("--extra", action="append", default=[], help="Extra candidate O atoms (e.g. a hydroxide)")
    p.add_argument("--water", default="water and element O", help="mdtraj selection of water oxygens")
    p.add_argument("--cutoff", type=float, default=5.0, help="Candidate search radius around the carbon (Å)")
    p.add_argument("--dmax", type=float, default=3.5)
    p.add_argument("--amin", type=float, default=95.0)
    p.add_argument("--amax", type=float, default=120.0)
    p.add_argument("--metal-cutoff", type=float, default=2.6)
    p.add_argument("--hb-dist", type=float, default=3.2)
    p.add_argument("--hb-hdist", type=float, default=2.5)
    p.add_argument("--stride", type=int, default=1)
    p.add_argument("--skip", type=int, default=0)
    p.add_argument("--outdir", default="nucleophiles")
    args = p.parse_args()

    traj = mdtraj.load(args.trajectory, top=args.topology, stride=args.stride)[args.skip:]
    top = traj.topology
    n = traj.n_frames
    c, o = resolve(top, args.carbon), resolve(top, args.oxygen)
    metal = resolve(top, args.metal) if args.metal else None
    bases = [resolve(top, b) for b in args.base]
    extra = [resolve(top, e) for e in args.extra]

    pool = sorted(set(int(i) for i in top.select(args.water)) | set(extra))
    near = mdtraj.compute_neighbors(traj, args.cutoff / 10.0, [c], haystack_indices=pool)
    cands = sorted(set(int(i) for f in near for i in f) | set(extra))
    print(f"{n} frames; {len(cands)} candidate O atoms within {args.cutoff} Å of {label(top, c)} in any frame")

    def hydrogens(oi):
        return [b.index for a, b in top.bonds if a.index == oi and b.element.symbol == "H"] + \
               [a.index for a, b in top.bonds if b.index == oi and a.element.symbol == "H"]

    d = mdtraj.compute_distances(traj, [[oi, c] for oi in cands]) * 10.0
    ang = np.degrees(mdtraj.compute_angles(traj, [[oi, c, o] for oi in cands]))
    attack = (d <= args.dmax) & (ang >= args.amin) & (ang <= args.amax)
    mbound = (mdtraj.compute_distances(traj, [[oi, metal] for oi in cands]) * 10.0 <= args.metal_cutoff) \
        if metal is not None else np.zeros_like(d, bool)

    hb = {}
    for b in bases:
        heavy = mdtraj.compute_distances(traj, [[oi, b] for oi in cands]) * 10.0 <= args.hb_dist
        hpoint = np.zeros_like(heavy)
        for k, oi in enumerate(cands):
            hs = hydrogens(oi)
            if hs:
                dh = mdtraj.compute_distances(traj, [[h, b] for h in hs]) * 10.0
                hpoint[:, k] = (dh <= args.hb_hdist).any(axis=1)
        hb[b] = heavy & hpoint
    any_hb = np.any(np.stack(list(hb.values())), axis=0) if hb else np.zeros_like(d, bool)

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for k, oi in enumerate(cands):
        within = d[:, k] <= args.cutoff
        rows.append({
            "candidate": label(top, oi),
            "occupancy": within.mean(),
            "d_min": d[:, k].min(),
            "d_mean_when_near": d[within, k].mean() if within.any() else np.nan,
            "best_angle_at_d<=4": ang[d[:, k] <= 4.0, k].max() if (d[:, k] <= 4.0).any() else np.nan,
            "attack": attack[:, k].mean(),
            "metal_bound": mbound[:, k].mean(),
            **{f"hb_{label(top, b)}": hb[b][:, k].mean() for b in bases},
            "attack_with_base_hb": (attack[:, k] & any_hb[:, k]).mean(),
            "attack_frames": " ".join(str(f * args.stride + args.skip) for f in np.where(attack[:, k])[0][:40]),
        })
    rows.sort(key=lambda r: (-r["attack_with_base_hb"], -r["attack"], r["d_min"]))
    keys = list(rows[0].keys())
    with open(out / "candidates.csv", "w") as f:
        f.write(",".join(keys) + "\n")
        for r in rows:
            f.write(",".join(f"{r[k]:.3f}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")

    # Per frame: the candidate closest to an ideal attack geometry
    score = np.abs(ang - 107.0) / 10.0 + np.maximum(d - 2.9, 0)
    best = score.argmin(axis=1)
    with open(out / "best_per_frame.csv", "w") as f:
        f.write("frame,candidate,d,angle,metal_bound," + ",".join(f"hb_{label(top, b)}" for b in bases) + "\n")
        for t in range(n):
            k = best[t]
            f.write(f"{t * args.stride + args.skip},{label(top, cands[k])},{d[t, k]:.2f},{ang[t, k]:.1f},"
                    f"{int(mbound[t, k])}," + ",".join(str(int(hb[b][t, k])) for b in bases) + "\n")

    show = [r for r in rows if r["occupancy"] >= 0.05 or r["attack"] > 0][:20]
    print(f"Candidates (occupancy >= 5% or any attack frame); attack = d <= {args.dmax} Å and "
          f"{args.amin:g}-{args.amax:g} deg; H-bond = O..acceptor <= {args.hb_dist} Å with an H <= {args.hb_hdist} Å")
    for r in show:
        hbs = " ".join(f"{k[3:]}={v:.0%}" for k, v in r.items() if k.startswith("hb_") and v > 0)
        print(f"  {r['candidate']:<24} near {r['occupancy']:5.0%}  dmin {r['d_min']:.2f}  "
              f"best angle(d<=4) {r['best_angle_at_d<=4']:.0f}  attack {r['attack']:5.1%}  "
              f"attack+baseHB {r['attack_with_base_hb']:5.1%}  metal {r['metal_bound']:4.0%}  {hbs}")
    print(f"Written {out}/candidates.csv and {out}/best_per_frame.csv")


if __name__ == "__main__":
    main()
