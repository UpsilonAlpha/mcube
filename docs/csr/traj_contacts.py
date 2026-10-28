#!/usr/bin/env python3
"""Named active-site contacts over a QM/MM trajectory.

    python traj_contacts.py QM_MM.dcd QM_MM_firstframe.pdb
    python traj_contacts.py --skip-missing <traj> <top>   # older ligand variants
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

# The nucleophile is residue OH in the hydroxide state and NUW once it has been
# protonated to water (docs/csr/make_water_nuc.py renames it, so that OpenMM
# matches it against water_nuc.xml rather than against every bulk HOH). The two
# arms of the protonation test therefore carry different residue names for the
# same oxygen, and the probe set has to find either.
NUCLEOPHILE_RESNAMES = ("OH", "NUW")


def nucleophile_resname(top):
    """Whichever of NUCLEOPHILE_RESNAMES this topology carries."""
    found = [name for name in NUCLEOPHILE_RESNAMES
             if any(res.name == name for res in top.residues)]
    if not found:
        raise ValueError(
            "no nucleophile residue: expected one of "
            f"{', '.join(NUCLEOPHILE_RESNAMES)}. Run --list to see what the "
            "topology holds.")
    if len(found) > 1:
        # Both at once means two nucleophiles, so every distance below would be
        # a minimum over them without saying so. That is a broken topology.
        raise ValueError(
            f"topology holds more than one nucleophile residue: {found}")
    return found[0]


def atoms(top, names=None, resseq=None, resname=None, element=None):
    """Atom indices matching every constraint given. Raises if none match."""
    # Waters carry the same resSeq as protein residues here (a HOH 226 sits
    # beside Lys226) and this refuses only an *empty* selection, not an
    # over-matching one, so probes name the residue as well as the number.
    picked, matched = [], set()
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
        matched.add(atom.name)
    if picked and names is not None:
        # A group selection survives when only some of its atoms exist, and then
        # measures a smaller group than its label claims - "carboxylate O3x/O4x"
        # against a ligand variant that has no O4x is a different quantity, not
        # a missing contact. Say so rather than letting the label stand for it.
        absent = sorted(names - matched)
        if absent:
            print(f"[NOTE] {resname or ''}{resseq or ''}: no "
                  f"{', '.join(absent)}; measuring the rest of {sorted(names)}")
    if not picked:
        raise ValueError(
            f"nothing matched names={names} resSeq={resseq} resname={resname} "
            f"element={element}. Run --list to see what the topology holds; a "
            "silently empty selection would be reported as a missing contact.")
    return picked


def probes(top, skip_missing=False):
    """(label, indices_a, indices_b) for every contact we want."""
    mer_o = dict(resname="mer", element="O")
    carboxylate = dict(resname="mer", names={"O3x", "O4x"})
    carbonyl = dict(resname="mer", names={"O2x"})
    guanidinium = {"NE", "NH1", "NH2"}
    nucleophile = dict(resname=nucleophile_resname(top), names={"O01"})

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
        # Gln93 = BBL Q121, a Zn2-site ligand position in canonical B3 and one of
        # the three B3-RQK substitutions. D121: a neutral amide standing in for the
        # metal as the oxyanion hole is documented across this subclass - Gln157 in
        # SMB-1, Asn176 in L1, Asn233 in CcrA - so the amide's reach to the
        # beta-lactam carbonyl is the measurement that tests the same role here.
        # Both amide heteroatoms are probed because the side-chain amide flips and
        # the crystal's assignment of OE1 against NE2 is not reliable either.
        ("Gln93 OE1 - mer O2x (carbonyl)",
         dict(resseq=93, resname="GLN", names={"OE1"}), carbonyl),
        ("Gln93 NE2 - mer O2x (carbonyl)",
         dict(resseq=93, resname="GLN", names={"NE2"}), carbonyl),
        ("Gln93 OE1/NE2 - mer O2x (carbonyl)",
         dict(resseq=93, resname="GLN", names={"OE1", "NE2"}), carbonyl),
        ("Gln93 NE2 - carboxylate O3x/O4x",
         dict(resseq=93, resname="GLN", names={"NE2"}), carboxylate),
        ("Gln93 NE2 - mer N1x (leaving N)",
         dict(resseq=93, resname="GLN", names={"NE2"}), dict(resname="mer", names={"N1x"})),
        ("Gln93 OE1/NE2 - O01 nucleophile",
         dict(resseq=93, resname="GLN", names={"OE1", "NE2"}), nucleophile),
        # D94 put Gln93 OE1 1.79 A from the zinc under MM, and D120 found Asp92 OD1
        # there instead. Whichever donor is in that position is in the carbonyl's,
        # so the amide's distance to the metal is measured beside its distance to
        # the substrate rather than inferred from it.
        ("Gln93 OE1 - Zn",
         dict(resseq=93, resname="GLN", names={"OE1"}), dict(resname="ZN")),
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
        # Meropenem can sit in the site either way round, and the two docked poses
        # this project has used do exactly that: the May pose presents the
        # beta-lactam carbonyl O2x to the metal at 2.55 A, the October one presents
        # the 6-alpha-hydroxyethyl O1x instead and leaves O2x at 5.43 A. Only the
        # first is a Michaelis complex (D121), so the orientation has to be a
        # measurement rather than something inferred from the pose's provenance.
        ("Zn - mer O1x (6a-OH)",
         dict(resname="ZN"), dict(resname="mer", names={"O1x"})),
        ("Zn - any mer O (which O faces the metal)",
         dict(resname="ZN"), mer_o),
        ("Zn - O01 nucleophile",
         dict(resname="ZN"), nucleophile),
    ]
    out, skipped = [], []
    for label, a, b in spec:
        try:
            out.append((label, atoms(top, **a), atoms(top, **b)))
        except ValueError:
            # Only with --skip-missing, and only an atom the system does not have.
            # The default stays strict because an empty selection silently reported
            # as a missing contact is the error this guard exists to prevent; but
            # the older systems carry earlier ligand variants (CSR_ligand's
            # meropenem has no O4x), and refusing to measure them at all would
            # mean comparing them by eye instead.
            if not skip_missing:
                raise
            skipped.append(label)
    if skipped:
        print("[SKIP] not present in this topology: " + "; ".join(skipped))
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
    skip_missing = "--skip-missing" in argv
    if skip_missing:
        argv.remove("--skip-missing")
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
          f"periodic={traj.unitcell_lengths is not None}, "
          f"nucleophile residue {nucleophile_resname(traj.topology)}")
    if traj.n_frames < 2:
        print("[WARN] fewer than two frames: nothing to average over")

    rows = probes(traj.topology, skip_missing=skip_missing)
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
