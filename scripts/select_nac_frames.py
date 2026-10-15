#!/usr/bin/env python3
"""
Select near-attack-conformation (NAC) frames from an MD trajectory.

For each frame, computes the nucleophile-electrophile distance and the
nucleophile-C=O angle (Burgi-Dunitz). It writes all values to a CSV and saves
the frames that meet the criteria as full-system PDBs (same atom order as the
topology), spread out in time.

Defaults match the CSR meropenem system (hydroxide O01 attacking C4x=O2x).

Usage:
    python select_nac_frames.py topology.pdb traj.dcd --outdir nac_frames \
        [--nuc "resname OH and name O01"] [--carbon "resname mer and name C4x"] \
        [--oxygen "resname mer and name O2x"] [--dmax 3.2] [--amin 95 --amax 120] \
        [--n 10] [--min-separation 50] [--skip 0]
"""

import argparse
from pathlib import Path

import mdtraj
import numpy as np
from openmm import unit
from openmm.app import PDBFile


def one(top, sel):
    idx = top.select(sel)
    if len(idx) != 1:
        raise SystemExit(f"Selection '{sel}' matched {len(idx)} atoms, expected 1")
    return int(idx[0])


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("topology")
    p.add_argument("trajectory")
    p.add_argument("--outdir", default="nac_frames")
    p.add_argument("--nuc", default="resname OH and name O01")
    p.add_argument("--carbon", default="resname mer MER and name C4x")
    p.add_argument("--oxygen", default="resname mer MER and name O2x")
    p.add_argument("--dmax", type=float, default=3.2, help="Max nucleophile-C distance (A)")
    p.add_argument("--amin", type=float, default=95.0, help="Min Nu-C=O angle (deg)")
    p.add_argument("--amax", type=float, default=120.0, help="Max Nu-C=O angle (deg)")
    p.add_argument("--n", type=int, default=10, help="Max number of frames to write")
    p.add_argument("--min-separation", type=int, default=50, help="Min frames between selected frames")
    p.add_argument("--skip", type=int, default=0, help="Ignore the first N frames (equilibration)")
    args = p.parse_args()

    traj = mdtraj.load(args.trajectory, top=args.topology)
    top = traj.topology
    nu, c, o = one(top, args.nuc), one(top, args.carbon), one(top, args.oxygen)

    dist = mdtraj.compute_distances(traj, [[nu, c]])[:, 0] * 10.0
    angle = np.degrees(mdtraj.compute_angles(traj, [[nu, c, o]])[:, 0])

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "nac_geometry.csv", "w") as f:
        f.write("frame,time_ps,dist_A,angle_deg\n")
        for i, (d, a) in enumerate(zip(dist, angle)):
            f.write(f"{i},{traj.time[i]:.1f},{d:.3f},{a:.1f}\n")

    ok = np.where((dist <= args.dmax) & (angle >= args.amin) & (angle <= args.amax)
                  & (np.arange(len(dist)) >= args.skip))[0]
    print(f"{len(traj)} frames; distance {dist.min():.2f}-{dist.max():.2f} A (mean {dist.mean():.2f}); "
          f"angle {angle.min():.0f}-{angle.max():.0f} deg (mean {angle.mean():.0f})")
    print(f"Frames meeting d <= {args.dmax} A and {args.amin} <= angle <= {args.amax} deg: {len(ok)} "
          f"({100 * len(ok) / len(dist):.1f}%)")

    # Spread the chosen frames out in time: closest to ideal first, then enforce separation
    ideal = np.abs(angle[ok] - 107.0) / 10.0 + np.abs(dist[ok] - 2.9)
    chosen = []
    for i in (ok[np.argsort(ideal)] if args.n > 0 else []):  # --n 0: only write the CSV
        if all(abs(i - j) >= args.min_separation for j in chosen):
            chosen.append(int(i))
        if len(chosen) == args.n:
            break
    # Written with OpenMM and the original topology, so atom names, residues and
    # ligand CONECT records match the input exactly (it must be an OpenMM-readable PDB)
    omm_top = PDBFile(args.topology).topology
    for i in sorted(chosen):
        name = out / f"nac_frame{i:05d}.pdb"
        with open(name, "w") as f:
            PDBFile.writeFile(omm_top, traj.xyz[i] * 10.0 * unit.angstrom, f, keepIds=True)
        print(f"  frame {i} (t = {traj.time[i]:.0f} ps): d = {dist[i]:.2f} A, angle = {angle[i]:.0f} deg -> {name}")
    if args.n > 0 and len(chosen) < 5:
        print("[WARNING] Fewer than 5 NAC frames; consider longer sampling or a restraint (see decision log D16).")


if __name__ == "__main__":
    main()
