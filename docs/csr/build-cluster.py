# -*- coding: utf-8 -*-
"""Build the gas-phase cluster models for the GFN2-vs-DFT benchmark (P1/D187).

The question this exists to answer is narrow and was asked by the user on
2026-11-03: **does GFN2-xTB hold meropenem's 6-(1-hydroxyethyl) oxygen O1x on
the zinc when DFT would not?** `md_eq5` found O1x coordinated in 97.9% of
10,000 frames at 2.18 A mean (D197), closer than the restrained carbonyl and as
close as the anionic hydroxide, and nothing in the MM or QM/MM machinery can
say whether that is real:

- MM cannot, because a non-bonded Zn(2+) with no Pauli repulsion pulls any lone
  pair in, and MM MD is documented to overfavour octahedral zinc (PMID
  34748329).
- QM/MM MD cannot, because **no QM/MM run in this project has ever crossed the
  coordination boundary in either direction** (D197 corrections): held in, a
  donor stays in; left out, it goes further out. Ligand exchange at Zn(II) is
  slower than the tens of ps affordable, so the trajectory length is what would
  be measured, not the chemistry.

What is left is a thermodynamic comparison on fixed geometries, which is
Imhof's own protocol - the paper this project's method follows - and the reason
plan item 7 has been open since D9. Imhof found AM1/d made a 4 -> 5 coordination
change **7.55 kcal/mol** too favourable; the equivalent number for GFN2 is
unmeasured and no GFN2 benchmark for zinc coordination number in enzymes
appears to exist.

**The measurement.** Two geometries, O1x on the metal and O1x off it, each
scored at GFN2-xTB and at DFT. The quantity is the *difference of differences*:

    dd = [E_off - E_on](GFN2)  -  [E_off - E_on](DFT)

A negative `dd` means GFN2 stabilises the O1x-bound form relative to DFT, i.e.
it over-coordinates, and `dd` is then the size of the error in kcal/mol. Because
both methods see the SAME two geometries, strain and thermal noise enter both
brackets and cancel in the difference; what does not cancel is the method.

**Both geometries come from the same trajectory**, not from hand displacement:
`on` is the frame with the shortest Zn-O1x and `off` the frame with the longest.
O1x reaches 3.19 A inside the 2.5 A occupancy analysis and 3.74 A in the 3.5 A
contact analysis, so an unbound frame exists and need not be invented. Moving an
atom by hand would put a strained geometry into one bracket only, and a
single-point energy on it would report that strain as chemistry.

**What this does NOT do.** It does not optimise either cluster, so these are not
binding energies. `sp` is a single-point stage by construction (D188) and an
optimisation at both levels is a separate, larger job. The difference-of-
differences is still the method comparison, which is what is being asked; a
reader wanting a binding energy needs the optimisations as well, and that should
be said wherever this number is quoted.

**The truncation, read from CONECT records rather than guessed** (D187 is
explicit that this is the one thing not to infer from distances - inferring
bonding by distance already put a spurious O-O bond in a strained geometry
once). From `mer_crystal.pdb`'s 51 CONECT records the real topology is:

    beta-lactam ring (4):  C4x - C3x - C17x - N1x - C4x,  O2x on C4x
    pyrroline ring (5):    N1x - C5x - C7x - C15x - C17x, fused at N1x-C17x
    C3 carboxylate:        C6x(O3x, O4x) on C5x
    1-hydroxyethyl:        C3x - C2x(C1x methyl)(O1x hydroxyl)
    thio side chain:       C7x - S1x - C8x ... N3x(C12x, C13x)

so the fragment keeping BOTH oxygens of interest is the intact beta-lactam ring
plus the hydroxyethyl: C4x, O2x, N1x, C3x, C17x, C2x, O1x, C1x.

**The C3 carboxylate is deliberately left out, which departs from D187's scope**
(it expected a ~15-atom fragment of net -1, implying the carboxylate). Two
reasons. It sits 5.93 and 7.05 A from the zinc through the whole of `md_eq5`
(D197), so it is not part of this contact. And in a *gas-phase* cluster a bare
carboxylate 6 A from a dication adds an enormous charge-separation term to both
brackets; it would largely cancel in `dd`, but it also destabilises the SCF and
makes a 1-10 kcal/mol effect ride on a far bigger number. Leaving it out makes
the cluster neutral, which is the better-behaved calculation. The site part of
D187's scope is unchanged.

**Charge and multiplicity.** Zn(2+) + acetate(-1) + hydroxide(-1) + two neutral
4-methylimidazoles + a neutral fragment = **0**, multiplicity 1 for d10 Zn(II).
Every `sp` section states `charge` explicitly: a section left at an inherited
value converges its SCF on the wrong species while looking entirely normal
(D152).

Usage, from the run directory:

    python build-cluster.py scan            # find the two frames, write them
    python build-cluster.py cut             # cut and cap both clusters
    python build-cluster.py                 # both, in order
"""

import sys
from pathlib import Path

SCRATCH = Path("/scratch/user/uquadhik")
RUNS = SCRATCH / "CSR"
RUN = RUNS / "CSR_mich"
RESULTS = RUN / "md_eq5_results"
OUT = RUN / "cluster"

TOPOLOGY = RESULTS / "NVTtrajectory_firstframe.pdb"
TRAJECTORY = RESULTS / "NVTtrajectory_imaged.dcd"

ZN = ("B", 1, "ZN")
O1X = ("C", 1, "O1x")

# Which heavy atoms of each residue survive the cut. A name not in this list is
# removed, and any bond it had to a kept atom is capped with hydrogen.
KEEP = {
    ("B", 1): None,                                   # the zinc, whole
    ("C", 2): None,                                   # the hydroxide, whole
    ("A", 92): ["CB", "CG", "OD1", "OD2"],            # -> acetate
    ("A", 88): ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],   # -> 4-methylimidazole
    ("A", 164): ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    ("C", 1): ["C4x", "O2x", "N1x", "C3x", "C17x", "C2x", "O1x", "C1x"],
}

# Cap length by the element of the atom that KEEPS the bond.
CAP_LENGTH = {"C": 1.09, "N": 1.01, "O": 0.97, "S": 1.34}

CHARGE = 0
MULT = 1


def read_pdb(path):
    """`[(chain, resseq, resname, name, element, x, y, z, line)]`, ATOM/HETATM only."""
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        element = (line[76:78].strip() or line[12:16].strip()[0]).upper()
        rows.append((line[21], int(line[22:26]), line[17:20].strip(),
                     line[12:16].strip(), element,
                     float(line[30:38]), float(line[38:46]), float(line[46:54]),
                     line))
    return rows


def dist(a, b):
    return sum((a[i] - b[i]) ** 2 for i in range(3)) ** 0.5


def scan():
    """Write the shortest and longest Zn-O1x frames of the trajectory as PDBs.

    Both geometries must come from the same force field and the same run, or the
    bracket difference carries a build difference as well as a coordination one.
    """
    import MDAnalysis as mda

    OUT.mkdir(parents=True, exist_ok=True)
    u = mda.Universe(str(TOPOLOGY), str(TRAJECTORY))

    def one(sel):
        chain, resid, name = sel
        g = u.select_atoms(f"segid {chain} and resid {resid} and name {name}")
        if g.n_atoms != 1:
            g = u.select_atoms(f"resid {resid} and name {name}")
        if g.n_atoms != 1:
            raise SystemExit(f"{sel} selected {g.n_atoms} atoms, want exactly 1")
        return g[0]

    zn, o1x = one(ZN), one(O1X)
    best = {"min": (1e9, -1), "max": (-1.0, -1)}
    for ts in u.trajectory:
        d = float(((zn.position - o1x.position) ** 2).sum() ** 0.5)
        if d < best["min"][0]:
            best["min"] = (d, ts.frame)
        if d > best["max"][0]:
            best["max"] = (d, ts.frame)
    print(f"[INFO] {len(u.trajectory)} frames scanned")
    for key, label in (("min", "on"), ("max", "off")):
        d, frame = best[key]
        u.trajectory[frame]
        path = OUT / f"frame_o1x_{label}.pdb"
        u.atoms.write(str(path))
        print(f"[OK] Zn-O1x {d:.2f} A at frame {frame} -> {path.name}")
    return 0


def cut(frame_path, out_path):
    """Cut the cluster out of one frame and cap every severed bond."""
    rows = read_pdb(frame_path)

    # Hydrogens are assigned to a heavy atom by distance. H makes exactly one
    # bond and the next-nearest heavy atom is never within 1.3 A, so this is
    # unambiguous - unlike heavy-atom bonding, which is why that came from
    # CONECT instead.
    heavy = [r for r in rows if r[4] != "H"]
    hydrogens = [r for r in rows if r[4] == "H"]

    wanted, missing = [], []
    for (chain, resid), names in KEEP.items():
        here = [r for r in heavy if r[0] == chain and r[1] == resid]
        if not here:
            missing.append(f"  {chain}:{resid} - no such residue in {frame_path}")
            continue
        if names is None:
            wanted += here
            continue
        for name in names:
            hit = [r for r in here if r[3] == name]
            if not hit:
                have = "/".join(sorted(r[3] for r in here))
                missing.append(f"  {chain}:{resid}:{name} - not found; has {have}")
            wanted += hit
    if missing:
        raise SystemExit("build-cluster: the truncation names atoms that are not "
                         f"in {frame_path}:\n" + "\n".join(missing) +
                         "\nNothing was written. The names come from the CONECT "
                         "records of mer_crystal.pdb and the ff14SB residue "
                         "templates, so a miss means the wrong input frame.")

    kept_ids = {id(r) for r in wanted}
    out, caps = [], []

    for r in wanted:
        out.append(r)
        # every hydrogen attached to a kept heavy atom comes along
        for h in hydrogens:
            if dist(r[5:8], h[5:8]) < 1.3:
                out.append(h)

    # Severed heavy-heavy bonds: a kept atom with a near heavy neighbour in the
    # SAME residue that was dropped. Capping only within the residue is right
    # here because the cut set is whole residues apart from the truncations, and
    # the only bonds broken are the ones those truncations break.
    for r in wanted:
        if r[4] == "ZN":
            continue
        for other in heavy:
            if id(other) in kept_ids or other[0] != r[0] or other[1] != r[1]:
                continue
            d = dist(r[5:8], other[5:8])
            if d > 1.95:
                continue
            length = CAP_LENGTH.get(r[4], 1.09)
            vec = [(other[5 + i] - r[5 + i]) / d for i in range(3)]
            pos = [r[5 + i] + vec[i] * length for i in range(3)]
            caps.append((r, other, pos))

    lines, serial = [], 0
    def emit(chain, resid, resname, name, element, xyz):
        nonlocal serial
        serial += 1
        lines.append("HETATM%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00"
                     "          %2s\n" % (serial, name[:4], resname[:3], chain,
                                          resid, xyz[0], xyz[1], xyz[2],
                                          element[:2]))

    for r in out:
        emit(r[0], r[1], r[2], r[3], r[4], r[5:8])
    for i, (keeper, dropped, pos) in enumerate(caps, 1):
        emit(keeper[0], keeper[1], keeper[2], "HC%d" % i, "H", pos)
    lines.append("END\n")

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_path).write_text("".join(lines))

    counts = {}
    for r in out:
        counts[r[4]] = counts.get(r[4], 0) + 1
    counts["H"] = counts.get("H", 0) + len(caps)

    by = lambda sel: next((r for r in out if (r[0], r[1], r[3]) == sel), None)
    zn, o1x = by(ZN), by(O1X)
    o2x = by(("C", 1, "O2x"))
    print(f"[OK] {out_path.name}: {serial} atoms, charge {CHARGE}, mult {MULT}")
    print("[INFO]   formula " + " ".join(f"{k}{v}" for k, v in sorted(counts.items())))
    print(f"[INFO]   {len(caps)} bond(s) capped: " +
          ", ".join(f"{k[2]}{k[1]}:{k[3]}-{d[3]}" for k, d, _ in caps))
    if zn and o1x:
        print(f"[INFO]   Zn-O1x {dist(zn[5:8], o1x[5:8]):.2f} A", end="")
    if zn and o2x:
        print(f"   Zn-O2x {dist(zn[5:8], o2x[5:8]):.2f} A", end="")
    print()
    return serial


def main(argv):
    what = argv[1] if len(argv) > 1 else "all"
    if what in ("scan", "all"):
        scan()
    if what in ("cut", "all"):
        sizes = {}
        for label in ("on", "off"):
            src = OUT / f"frame_o1x_{label}.pdb"
            if not src.is_file():
                raise SystemExit(f"missing {src} - run `scan` first")
            sizes[label] = cut(src, OUT / f"cluster_o1x_{label}.pdb")
        if sizes["on"] != sizes["off"]:
            raise SystemExit(
                "build-cluster: the two clusters have different atom counts "
                f"({sizes}). The difference of differences is only a method "
                "comparison if both methods see the same two systems, so an "
                "asymmetric cut invalidates the benchmark. Most likely a cap "
                "was placed in one frame and not the other because a severed "
                "bond stretched past 1.95 A.")
        print(f"[DONE] both clusters are {sizes['on']} atoms")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
