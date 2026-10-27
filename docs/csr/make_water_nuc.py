#!/usr/bin/env python3
"""Build the ohH input structure: put the proton back on the nucleophile.

    python make_water_nuc.py csr_crystal_start.pdb csr_crystal_start_ohH.pdb

The lysH arm has a hydroxide (residue OH, atoms O01 + H03) with a cationic
Lys226. The ohH arm is the same site with that proton moved onto the
nucleophile, so this adds one hydrogen, H04, and renames the residue OH -> NUW
to match `water_nuc.xml`. Lys226 is turned neutral separately, by
`residue_variants: {A: {226: LYN}}` in preprocess - OpenMM's LYN variant drops
HZ3, so the two arms have the same atom count and the same total charge, and the
QM charge stays +3.

H04 goes where a metal-bound water's second hydrogen belongs: in the plane of
Zn-O01-H03, at the TIP3P geometry (0.9572 A, 104.52 degrees from H03), on
whichever side points further from the metal. That is a geometric choice, not a
mechanistic one - it does not aim the new proton at Asp92 or at anything else,
and the restrained MM equilibration that follows is what actually settles it.

Refuses rather than guesses: missing atoms, a degenerate Zn-O01-H03 geometry, or
a placement that lands on top of an existing atom all raise. The last matters
because a hydrogen dropped into a clash looks fine in a PDB and blows up a
minimisation later.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

OH_LEN = 0.9572          # A, TIP3P
HOH_ANGLE = 104.52       # degrees, TIP3P
CLASH = 1.5              # A; nothing but O01 may be closer than this to H04


def parse(line):
    """(name, resname, chain, resseq, x, y, z) from an ATOM/HETATM line."""
    return (line[12:16].strip(), line[17:20].strip(), line[21],
            line[22:26].strip(),
            float(line[30:38]), float(line[38:46]), float(line[46:54]))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    length = norm(a)
    if length < 1e-8:
        raise ValueError("cannot normalise a zero-length vector")
    return tuple(x / length for x in a)


def place_h04(zn, o01, h03):
    """The second hydrogen of a metal-bound water, in the Zn-O-H plane."""
    a = unit(sub(h03, o01))                 # O -> H03
    away = unit(sub(o01, zn))               # O -> away from the metal
    normal = cross(a, away)
    if norm(normal) < 1e-6:
        raise ValueError(
            "Zn, O01 and H03 are collinear, so the Zn-O01-H03 plane is "
            "undefined and there is no 'away from the metal' side to pick. "
            "Place H04 by hand for this structure.")
    normal = unit(normal)
    perp = unit(cross(normal, a))           # in plane, perpendicular to a
    theta = math.radians(HOH_ANGLE)
    options = []
    for sign in (1.0, -1.0):
        direction = tuple(math.cos(theta) * ai + sign * math.sin(theta) * pi
                          for ai, pi in zip(a, perp))
        options.append((dot(direction, away), direction))
    _, best = max(options)                  # the one pointing further from Zn
    return tuple(o + OH_LEN * d for o, d in zip(o01, best))


def parse_lys(spec):
    """"A:226" -> ("A", "226")."""
    parts = spec.split(":")
    if len(parts) != 2 or not all(parts):
        raise SystemExit(f"--lys wants CHAIN:RESID, not {spec!r}")
    return parts[0], parts[1]


def main(argv):
    lys = None
    if "--lys" in argv:
        at = argv.index("--lys")
        if at + 1 >= len(argv):
            raise SystemExit("--lys needs a CHAIN:RESID argument")
        lys = parse_lys(argv[at + 1])
        argv = argv[:at] + argv[at + 2:]
    if len(argv) != 2:
        raise SystemExit(__doc__)
    src, dest = Path(argv[0]), Path(argv[1])
    if not src.is_file():
        raise SystemExit(f"no such file: {src}")
    lines = src.read_text().splitlines(keepends=True)

    zn = o01 = h03 = None
    oh_index = oh_line = None
    nuc_residue = None
    lys_atoms = {}
    coords = []
    for index, line in enumerate(lines):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name, resname, chain, resseq, x, y, z = parse(line)
        coords.append(((chain, resseq, name), (x, y, z)))
        if resname == "ZN" or name == "ZN":
            zn = (x, y, z)
        if resname == "OH":
            nuc_residue = (chain, resseq)
            if name == "O01":
                o01, oh_index, oh_line = (x, y, z), index, line
            elif name == "H03":
                h03, oh_index = (x, y, z), index   # insert after H03
        if resname == "NUW":
            raise SystemExit(f"{src} already holds a NUW residue; this has "
                             "been run already.")
        if lys and (chain, resseq) == lys:
            lys_atoms[name] = index

    missing = [n for n, v in (("Zn", zn), ("OH O01", o01), ("OH H03", h03))
               if v is None]
    if missing:
        raise SystemExit(
            f"could not find {', '.join(missing)} in {src}. This script expects "
            "the lysH input: a residue named OH with atoms O01 and H03, and a "
            "zinc. Nothing was written.")

    # The lysine's proton is the one that moved, so it has to go. Its residue
    # keeps the name LYS: OpenMM matches the force-field template by bond
    # graph and finds LYN anyway, whereas renaming it would make
    # `build_qm_region` treat it as non-standard - the config's `amino_acids`
    # list overrides mcube's default and does not include LYN - and put its
    # whole backbone into the QM region, changing the region between the arms.
    drop_index = None
    if lys:
        if not lys_atoms:
            raise SystemExit(
                f"--lys {lys[0]}:{lys[1]} matched no atoms in {src}. Nothing "
                "was written.")
        if "HZ3" not in lys_atoms:
            have = " ".join(sorted(lys_atoms))
            raise SystemExit(
                f"{lys[0]}:{lys[1]} has no HZ3, so it is not a cationic "
                f"lysine to neutralise. It holds: {have}. Nothing was written.")
        if "HZ1" not in lys_atoms or "HZ2" not in lys_atoms:
            raise SystemExit(
                f"{lys[0]}:{lys[1]} is missing HZ1 or HZ2, so removing HZ3 "
                "would not leave a neutral lysine. Nothing was written.")
        drop_index = lys_atoms["HZ3"]

    h04 = place_h04(zn, o01, h03)

    # The clash check looks at OTHER residues only. H04's nearest neighbour is
    # always its own partner hydrogen, at 2*0.9572*sin(104.52/2) = 1.513 A,
    # which sits a hundredth of an angstrom above the cutoff - so counting the
    # nucleophile's own atoms would make every correct placement a near miss
    # and some of them a spurious refusal.
    near = sorted((norm(sub(h04, xyz)), key) for key, xyz in coords
                  if (key[0], key[1]) != nuc_residue)
    if not near:
        raise SystemExit(f"{src} holds nothing but the nucleophile.")
    distance, key = near[0]
    if distance < CLASH:
        raise SystemExit(
            f"H04 would sit {distance:.2f} A from {':'.join(key)}, inside the "
            f"{CLASH} A clash cutoff. Nothing was written - the placement "
            "needs doing by hand for this structure.")
    closest = near[0]

    # Build the new H04 line from the O01 line so the chain, residue number and
    # formatting are inherited rather than reconstructed.
    h04_line = (oh_line[:6] + oh_line[6:12] + f"{'H04':>4s}" + oh_line[16]
                + f"{'NUW':>3s}" + oh_line[20:30]
                + f"{h04[0]:8.3f}{h04[1]:8.3f}{h04[2]:8.3f}"
                + oh_line[54:])

    out = []
    dropped_serial = None
    for index, line in enumerate(lines):
        if index == drop_index:
            dropped_serial = int(line[6:11])
            continue
        if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == "OH":
            line = line[:17] + f"{'NUW':>3s}" + line[20:]
        out.append(line)
        if index == oh_index:
            out.append(h04_line)

    # Renumber serials so the inserted atom does not leave a duplicate, keeping
    # old -> new so the CONECT records can be rewritten. The nucleophile's bonds
    # come from CONECT (PDBFile infers bonds by atom name, and O01/H03/H04 are
    # not a standard water's names), so a stale CONECT would leave the new
    # hydrogen unbonded and the force field would reject the residue.
    remap = {}
    serial = 0
    h04_serial = None
    for index, line in enumerate(out):
        if line.startswith(("ATOM", "HETATM")):
            serial += 1
            old = line[6:11].strip()
            if line is h04_line or line[12:16].strip() == "H04":
                h04_serial = serial
            elif old:
                remap[int(old)] = serial
            out[index] = line[:6] + f"{serial:5d}" + line[11:]

    o01_serial = remap.get(int(oh_line[6:11]))
    conect_seen = False
    for index, line in enumerate(out):
        if not line.startswith("CONECT"):
            continue
        conect_seen = True
        fields = [line[pos:pos + 5].strip() for pos in range(6, 31, 5)]
        new = "CONECT"
        for field in fields:
            if not field:
                continue
            if dropped_serial is not None and int(field) == dropped_serial:
                continue    # the lysine proton we removed
            if int(field) not in remap:
                raise SystemExit(
                    f"CONECT refers to serial {field}, which is not an atom in "
                    f"{src}. Nothing was written.")
            new += f"{remap[int(field)]:5d}"
        # The O01 record gains the new hydrogen.
        if fields and fields[0] and remap[int(fields[0])] == o01_serial:
            new += f"{h04_serial:5d}"
        out[index] = new + "\n"

    if not conect_seen:
        raise SystemExit(
            f"{src} has no CONECT records, so the nucleophile's bonds cannot be "
            "expressed and the force field would reject residue NUW. Nothing "
            "was written.")

    # The reciprocal record goes before END: PDB readers stop there, so a
    # CONECT appended after it would be dropped without complaint and the new
    # hydrogen would arrive unbonded.
    reciprocal = f"CONECT{h04_serial:5d}{o01_serial:5d}\n"
    end = next((i for i, l in enumerate(out) if l.startswith("END")), len(out))
    out.insert(end, reciprocal)

    dest.write_text("".join(out))
    print(f"H04 placed at {h04[0]:.3f} {h04[1]:.3f} {h04[2]:.3f}")
    print(f"  O01-H04 {norm(sub(h04, o01)):.4f} A, "
          f"H03-O01-H04 {math.degrees(math.acos(max(-1.0, min(1.0, dot(unit(sub(h03, o01)), unit(sub(h04, o01))))))):.2f} deg")
    print(f"  Zn-H04 {norm(sub(h04, zn)):.3f} A "
          f"(Zn-H03 {norm(sub(h03, zn)):.3f} A)")
    print(f"  nearest other atom: {':'.join(closest[1])} at {closest[0]:.2f} A")
    if lys:
        print(f"  {lys[0]}:{lys[1]} HZ3 removed (serial {dropped_serial}); the "
              "residue keeps the name LYS and matches LYN by bond graph")
    print(f"  residue OH renamed NUW; wrote {dest} ({serial} atoms)")


if __name__ == "__main__":
    main(sys.argv[1:])
