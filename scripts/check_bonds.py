#!/usr/bin/env python3
"""
Check that bond geometry is preserved along a trajectory, for a selection of
atoms (by default the QM region saved by `mcube preprocess`).

Per frame it flags:
- stretched/compressed bonds: |d - d_ref| > --tol (topology bonds within the
  selection; d_ref = length in the first frame, or in --ref);
- new contacts: two non-bonded heavy atoms of the selection closer than
  --contact Å (separated by more than two bonds), i.e. a bond forming;
- proton transfers: an H closer to another heavy atom than to the atom it is
  bonded to.
Metal–ligand pairs are not covalent bonds in the topology, so they are not
checked here (use site_report.py for the coordination shell).

Writes bond_frames.csv (frame, n_stretched, n_contacts, n_htransfer, intact)
and prints a summary with the worst bonds.

Usage:
    python check_bonds.py topology.pdb traj.dcd [--config config.yaml | --sel "resname MER"]
        [--tol 0.25] [--contact 1.75] [--stride 1] [--outdir bonds]
"""

import argparse
from pathlib import Path

import mdtraj
import numpy as np
import yaml


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("topology")
    p.add_argument("trajectory")
    p.add_argument("--config", help="Use preprocess.results.qm_atoms from this mcube config")
    p.add_argument("--sel", help="mdtraj selection instead of the QM atoms")
    p.add_argument("--ref", help="Reference structure for bond lengths (default: first frame)")
    p.add_argument("--tol", type=float, default=0.25, help="Allowed bond-length change (Å)")
    p.add_argument("--contact", type=float, default=1.75, help="Non-bonded heavy-atom contact flag (Å)")
    p.add_argument("--stride", type=int, default=1)
    p.add_argument("--outdir", default="bonds")
    args = p.parse_args()

    traj = mdtraj.load(args.trajectory, top=args.topology, stride=args.stride)
    top = traj.topology
    if args.sel:
        atoms = set(int(i) for i in top.select(args.sel))
    elif args.config:
        atoms = set((yaml.safe_load(open(args.config))["preprocess"]["results"])["qm_atoms"])
    else:
        raise SystemExit("give --config (QM atoms) or --sel")

    bonds = sorted({tuple(sorted((a.index, b.index))) for a, b in top.bonds
                    if a.index in atoms and b.index in atoms})
    if not bonds:
        raise SystemExit("no topology bonds inside the selection (does the topology PDB have CONECT records?)")
    label = lambda i: f"{top.atom(i).residue}:{top.atom(i).name}({i})"  # noqa: E731
    d = mdtraj.compute_distances(traj, bonds) * 10.0
    ref = (mdtraj.compute_distances(mdtraj.load(args.ref), bonds)[0] * 10.0) if args.ref else d[0]
    dev = d - ref
    stretched = np.abs(dev) > args.tol

    # Bond graph within the selection, for 1-2/1-3 exclusions and H partners
    nbr = {i: set() for i in atoms}
    for a, b in bonds:
        nbr[a].add(b)
        nbr[b].add(a)
    heavy = sorted(i for i in atoms if top.atom(i).element.symbol != "H")
    close = set()
    for a in heavy:
        close |= {(a, b) for b in nbr[a] | {c for n in nbr[a] for c in nbr[n]}}
    pairs = [(a, b) for k, a in enumerate(heavy) for b in heavy[k + 1:] if (a, b) not in close]
    cd = mdtraj.compute_distances(traj, pairs) * 10.0 if pairs else np.zeros((traj.n_frames, 0))
    contacts = cd < args.contact

    hs = [i for i in atoms if top.atom(i).element.symbol == "H" and len(nbr[i]) == 1]
    htrans = np.zeros((traj.n_frames, len(hs)), bool)
    h_detail = []
    for k, h in enumerate(hs):
        x = next(iter(nbr[h]))
        others = [y for y in heavy if y != x]
        dx = mdtraj.compute_distances(traj, [[h, x]])[:, 0]
        dy = mdtraj.compute_distances(traj, [[h, y] for y in others])
        j = dy.argmin(axis=1)
        htrans[:, k] = dy[np.arange(traj.n_frames), j] < dx
        if htrans[:, k].any():
            f = int(np.argmax(htrans[:, k]))
            h_detail.append(f"{label(h)} leaves {label(x)} for {label(others[j[f]])} "
                            f"in {htrans[:, k].mean():.0%} of frames (first frame {f * args.stride})")

    n_s, n_c, n_h = stretched.sum(1), contacts.sum(1), htrans.sum(1)
    intact = (n_s == 0) & (n_c == 0) & (n_h == 0)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "bond_frames.csv", "w") as f:
        f.write("frame,n_stretched,n_contacts,n_htransfer,intact\n")
        for t in range(traj.n_frames):
            f.write(f"{t * args.stride},{n_s[t]},{n_c[t]},{n_h[t]},{int(intact[t])}\n")

    print(f"{traj.n_frames} frames (stride {args.stride}); {len(atoms)} atoms, {len(bonds)} bonds, "
          f"{len(pairs)} non-bonded heavy pairs, {len(hs)} H")
    print(f"Frames with all bonds intact: {intact.sum()} / {traj.n_frames}")
    worst = np.argsort(-np.abs(dev).max(0))[:8]
    print(f"Largest bond-length changes (tolerance {args.tol} Å):")
    for k in worst:
        a, b = bonds[k]
        print(f"  {label(a)}-{label(b)}: ref {ref[k]:.3f}, range {d[:, k].min():.3f}-{d[:, k].max():.3f} Å, "
              f"out of tolerance in {stretched[:, k].mean():.0%} of frames")
    if contacts.any():
        for k in np.where(contacts.any(0))[0][:10]:
            a, b = pairs[k]
            print(f"  NEW CONTACT {label(a)}..{label(b)}: min {cd[:, k].min():.2f} Å in {contacts[:, k].mean():.0%} of frames")
    else:
        print(f"No non-bonded heavy-atom contacts below {args.contact} Å")
    print("\n".join(h_detail) if h_detail else "No proton transfers")


if __name__ == "__main__":
    main()
