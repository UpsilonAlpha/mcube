#!/usr/bin/env python3
"""Named active-site contacts over a QM/MM trajectory.

    python traj_contacts.py QM_MM.dcd QM_MM_firstframe.pdb
    python traj_contacts.py --list QM_MM_firstframe.pdb

D107 measured four of these by hand for the `nacqm` run; this makes the set
reproducible and adds the arginines. The hypothesis under test (D107, D110) is
that Arg90 and Lys226, as two monovalent cations at metal-ligand positions,
supply the +2 the absent Zn2 would have carried - anchoring the substrate
electrostatically while being unable to do Zn2's catalytic jobs. So the
quantities that matter are each cation's reach to the substrate's anionic
centres, over a trajectory rather than in one frame: a single optimised frame
cannot tell a spectator from a residue that is merely elsewhere at that instant,
which is the error D107 had to correct.

Distances are minimum-image (the box vectors are in the DCD) and in angstroms.
A group on either side means the closest pair is taken per frame, and the atoms
that achieved it are named, so "Arg229 reaches the substrate" does not hide
which oxygen it reached.
"""

from __future__ import annotations

import sys
from pathlib import Path

CUTOFF = 3.5  # A; salt bridge / H-bond heavy-atom contact


def atoms(top, names=None, resseq=None, resname=None, element=None):
    """Atom indices matching every constraint given. Raises if none match."""
    # Waters carry the same resSeq as protein residues here (a HOH 226 sits
    # beside Lys226) and this refuses only an *empty* selection, not an
    # over-matching one, so probes name the residue as well as the number.
    picked = []
    for atom in top.atoms:
        if resseq is not None and atom.residue.resSeq != resseq:
            continue
        if resname is not None and atom.residue.name != resname:
            continue
        if names is not None and atom.name not in names:
            continue
        if element is not None and atom.element.symbol != element:
            continue
        picked.append(atom.index)
    if not picked:
        raise ValueError(
            f"nothing matched names={names} resSeq={resseq} resname={resname} "
            f"element={element}. Run --list to see what the topology holds; a "
            "silently empty selection would be reported as a missing contact.")
    return picked


def probes(top):
    """(label, indices_a, indices_b) for every contact we want."""
    mer_o = dict(resname="mer", element="O")
    carboxylate = dict(resname="mer", names={"O3x", "O4x"})
    guanidinium = {"NE", "NH1", "NH2"}
    nucleophile = dict(resname="OH", names={"O01"})

    spec = [
        # Lys226 = BBL K263, at the vacant Zn2 ligand position (D107)
        ("Lys226 NZ - O01 nucleophile",
         dict(resseq=226, resname="LYS", names={"NZ"}), nucleophile),
        ("Lys226 NZ - mer O3x",
         dict(resseq=226, resname="LYS", names={"NZ"}), dict(resname="mer", names={"O3x"})),
        ("Lys226 NZ - mer O4x",
         dict(resseq=226, resname="LYS", names={"NZ"}), dict(resname="mer", names={"O4x"})),
        ("Lys226 NZ - mer N1x (leaving N)",
         dict(resseq=226, resname="LYS", names={"NZ"}), dict(resname="mer", names={"N1x"})),
        ("Lys226 NZ - any mer O",
         dict(resseq=226, resname="LYS", names={"NZ"}), mer_o),
        # Arg90 = BBL R118, at a Zn1 ligand position
        ("Arg90 guanidinium - carboxylate O3x/O4x",
         dict(resseq=90, resname="ARG", names=guanidinium), carboxylate),
        ("Arg90 guanidinium - any mer O",
         dict(resseq=90, resname="ARG", names=guanidinium), mer_o),
        ("Arg90 guanidinium - O01 nucleophile",
         dict(resseq=90, resname="ARG", names=guanidinium), nucleophile),
        # Arg229, in the QM region for the nacarg run only (D110)
        ("Arg229 guanidinium - carboxylate O3x/O4x",
         dict(resseq=229, resname="ARG", names=guanidinium), carboxylate),
        ("Arg229 guanidinium - any mer O",
         dict(resseq=229, resname="ARG", names=guanidinium), mer_o),
        ("Arg229 guanidinium - O01 nucleophile",
         dict(resseq=229, resname="ARG", names=guanidinium), nucleophile),
        # Asp92 = D120, the proposed general base; monodentate to the one zinc
        ("Asp92 OD1 - O01 nucleophile",
         dict(resseq=92, resname="ASP", names={"OD1"}), nucleophile),
        ("Asp92 OD1/OD2 - O01 nucleophile",
         dict(resseq=92, resname="ASP", names={"OD1", "OD2"}), nucleophile),
        # The attack coordinate and the metal, to sit beside the log's table
        ("O01 - mer C4x (attack)",
         nucleophile, dict(resname="mer", names={"C4x"})),
        ("Zn - mer O2x",
         dict(resname="ZN"), dict(resname="mer", names={"O2x"})),
        ("Zn - O01 nucleophile",
         dict(resname="ZN"), nucleophile),
    ]
    out = []
    for label, a, b in spec:
        out.append((label, atoms(top, **a), atoms(top, **b)))
    return out


def list_topology(pdb):
    import mdtraj as md

    top = md.load(pdb).topology
    print(f"{top.n_atoms} atoms, {top.n_residues} residues")
    print("\nNon-protein, non-water residues and their atoms:")
    skip = {"HOH", "WAT", "NA", "CL"}
    for res in top.residues:
        if res.is_protein or res.name in skip:
            continue
        names = " ".join(a.name for a in res.atoms)
        print(f"  {res.name:5s} resSeq {res.resSeq:<5d} {names}")
    print("\nActive-site residues by number:")
    for num in (90, 92, 93, 226, 229):
        hit = [r for r in top.residues if r.resSeq == num]
        for res in hit:
            print(f"  resSeq {num}: {res.name}, "
                  f"{' '.join(a.name for a in res.atoms)}")


def main(argv):
    if "--list" in argv:
        argv.remove("--list")
        if len(argv) != 1:
            raise SystemExit("usage: traj_contacts.py --list <topology.pdb>")
        return list_topology(argv[0])

    if len(argv) != 2:
        raise SystemExit(__doc__)
    traj_path, top_path = argv
    for path in (traj_path, top_path):
        if not Path(path).is_file():
            raise SystemExit(f"no such file: {path}")

    import mdtraj as md
    import numpy as np

    traj = md.load(traj_path, top=top_path)
    print(f"{traj.n_frames} frames, {traj.n_atoms} atoms, "
          f"periodic={traj.unitcell_lengths is not None}")
    if traj.n_frames < 2:
        print("[WARN] fewer than two frames: nothing to average over")

    rows = probes(traj.topology)
    # One compute_distances call for every pair in every probe, then split it
    # back up: mdtraj's minimum-image code is the slow part, not the bookkeeping.
    pairs, spans = [], []
    for _, a, b in rows:
        start = len(pairs)
        pairs.extend([(i, j) for i in a for j in b])
        spans.append((start, len(pairs)))
    dist = md.compute_distances(traj, np.array(pairs), periodic=True) * 10.0

    names = {atom.index: f"{atom.residue.name}{atom.residue.resSeq}:{atom.name}"
             for atom in traj.topology.atoms}
    width = max(len(label) for label, _, _ in rows)
    print(f"\n{'contact':{width}s}   min   max  mean   %<{CUTOFF}  closest pair")
    for (label, _, _), (start, stop) in zip(rows, spans):
        block = dist[:, start:stop]
        per_frame = block.min(axis=1)
        closest = pairs[start + int(block.min(axis=0).argmin())]
        frac = 100.0 * float((per_frame < CUTOFF).mean())
        print(f"{label:{width}s} {per_frame.min():5.2f} {per_frame.max():5.2f} "
              f"{per_frame.mean():5.2f} {frac:6.1f}  "
              f"{names[closest[0]]} - {names[closest[1]]}")


if __name__ == "__main__":
    main(sys.argv[1:])
