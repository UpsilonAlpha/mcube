#!/usr/bin/env python3
"""Build the ohH input structure: put the proton back on the nucleophile.

    python make_water_nuc.py in.pdb out.pdb [--lys A:226]

The lysH arm has a hydroxide (residue OH, atoms O01 + H03) with a cationic
Lys226. The ohH arm is the same site with that proton moved onto the
nucleophile, so this adds one hydrogen, H04, and renames the residue OH -> NUW
to match `water_nuc.xml`.

`--lys A:226` completes the move by removing that lysine's HZ3. One hydrogen
leaves the lysine and one joins the nucleophile, so both arms keep the same atom
count, the same total charge, and the same QM charge of +3. The residue
deliberately keeps the name LYS: OpenMM matches the force-field template by bond
graph and finds LYN anyway, whereas renaming it would make `build_qm_region`
treat it as non-standard - the config's `amino_acids` list overrides mcube's
default and omits LYN - and put its whole backbone into the QM region, changing
the QM region between the two arms. Use `--lys` when starting from an
equilibrated frame, which is what the ohH run does (naclong's 100 ns MM
lastframe, the frame nacqm started from). A run that goes through preprocess can
set `residue_variants: {A: {226: LYN}}` instead, since Modeller adds hydrogens
there.

H04 goes where a metal-bound water's second hydrogen belongs: in the plane of
Zn-O01-H03, at the TIP3P geometry (0.9572 A, 104.52 degrees from H03), on
whichever side points further from the metal. That is a geometric choice, not a
mechanistic one - it does not aim the new proton at Asp92 or at anything else,
and the first picoseconds of QM/MM are what actually settle it.

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
CLASH = 1.9              # A. An H...H contact at 1.5 A is strongly repulsive -
                         # normal is 2.2-2.4 - and the first MD step would feel
                         # it. 1.9 is a floor, not a target; the cone scan
                         # reports the clearance it actually achieved.


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


def cone_positions(o01, h03, steps=72):
    """Every H04 position allowed by the TIP3P geometry, as (angle, xyz).

    With O01 and H03 fixed, the bond length and the H-O-H angle leave one free
    parameter: H04 lies on a cone about the O01->H03 axis. Picking a position in
    the Zn-O-H plane is only one point on that cone, and in a solvated frame it
    can land on top of the substrate - it did, 1.51 A from meropenem's H5x. So
    the cone is scanned instead and the clearest position taken.
    """
    axis = unit(sub(h03, o01))
    # Any vector not parallel to the axis, to start an orthonormal pair from.
    seed = (1.0, 0.0, 0.0)
    if abs(dot(seed, axis)) > 0.9:
        seed = (0.0, 1.0, 0.0)
    u = unit(cross(axis, seed))
    v = unit(cross(axis, u))
    theta = math.radians(HOH_ANGLE)
    out = []
    for step in range(steps):
        phi = 2.0 * math.pi * step / steps
        radial = tuple(math.cos(phi) * ui + math.sin(phi) * vi
                       for ui, vi in zip(u, v))
        direction = tuple(math.cos(theta) * ai + math.sin(theta) * ri
                          for ai, ri in zip(axis, radial))
        out.append((math.degrees(phi),
                    tuple(o + OH_LEN * d for o, d in zip(o01, direction))))
    return out


def sphere_directions(count=240):
    """Roughly even directions on a sphere (Fibonacci spiral)."""
    out = []
    golden = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(count):
        z = 1.0 - 2.0 * i / (count - 1.0)
        radius = math.sqrt(max(0.0, 1.0 - z * z))
        phi = golden * i
        out.append((math.cos(phi) * radius, math.sin(phi) * radius, z))
    return out


def reorient_water(zn, o01, h03, others, floor):
    """Both hydrogens placed afresh, with the oxygen left exactly where it is.

    Used when every position on the cone is blocked, which is what happens in
    naclong's frame: the restrained equilibration pressed the hydroxide against
    meropenem, and a second hydrogen does not fit at that orientation. H03's
    direction is a *hydroxide's* and has no claim on a water's geometry, while
    O01 is the nucleophile and the one coordinate the comparison rests on - so
    the oxygen stays put and the hydrogens are free. Among orientations that
    clear `floor`, the one keeping H03 closest to where it was wins, so the
    structure stays as near the input as the clash allows.
    """
    near = [(key, xyz) for key, xyz in others if norm(sub(xyz, o01)) < 6.0]
    original = unit(sub(h03, o01))
    away = unit(sub(o01, zn))
    theta = math.radians(HOH_ANGLE)
    best = None
    for d1 in sphere_directions():
        p1 = tuple(o + OH_LEN * d for o, d in zip(o01, d1))
        gap1 = min((norm(sub(p1, xyz)) for _, xyz in near), default=float("inf"))
        if gap1 < floor:
            continue
        seed = (1.0, 0.0, 0.0) if abs(dot((1.0, 0.0, 0.0), d1)) < 0.9 else (0.0, 1.0, 0.0)
        u = unit(cross(d1, seed))
        v = unit(cross(d1, u))
        for step in range(36):
            phi = 2.0 * math.pi * step / 36
            radial = tuple(math.cos(phi) * ui + math.sin(phi) * vi
                           for ui, vi in zip(u, v))
            d2 = tuple(math.cos(theta) * ai + math.sin(theta) * ri
                       for ai, ri in zip(d1, radial))
            p2 = tuple(o + OH_LEN * d for o, d in zip(o01, d2))
            gap2 = min((norm(sub(p2, xyz)) for _, xyz in near),
                       default=float("inf"))
            if gap2 < floor:
                continue
            clearance = min(gap1, gap2)
            # Stay as close to the input orientation as the clash allows, then
            # prefer hydrogens pointing away from the metal.
            score = (dot(d1, original), dot(d1, away) + dot(d2, away))
            if best is None or score > best[0]:
                best = (score, clearance, p1, p2)
    return best


def place_h04(zn, o01, h03, others=()):
    """The clearest allowed position for the water's second hydrogen.

    `others` is (key, xyz) for every atom that is not part of the nucleophile.
    The position chosen is the one whose closest approach to them is largest;
    among positions within 0.1 A of that, the one furthest from the metal wins,
    which keeps the chemical preference for both hydrogens pointing away from
    the cation as a tiebreak rather than as the rule.
    """
    if norm(sub(h03, o01)) < 1e-6:
        raise ValueError("O01 and H03 are at the same position")
    away = unit(sub(o01, zn))
    # Only atoms near the nucleophile can be the closest approach, and scanning
    # the cone against all 49000 would be needlessly slow.
    near = [(key, xyz) for key, xyz in others
            if norm(sub(xyz, o01)) < 6.0]
    scored = []
    for phi, xyz in cone_positions(o01, h03):
        if near:
            clearance, key = min((norm(sub(xyz, other)), k)
                                 for k, other in near)
        else:
            clearance, key = float("inf"), None
        scored.append((clearance, dot(unit(sub(xyz, o01)), away), phi, xyz, key))
    best = max(scored)[0]
    viable = [s for s in scored if s[0] >= best - 0.1]
    clearance, _, phi, xyz, key = max(viable, key=lambda s: s[1])
    return xyz, clearance, key, phi


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

    # The nucleophile's own atoms are excluded: H04's partner hydrogen is
    # always 2*0.9572*sin(104.52/2) = 1.513 A away by construction, so counting
    # it would make every correct placement look like a near miss.
    others = [(key, xyz) for key, xyz in coords
              if (key[0], key[1]) != nuc_residue]
    if not others:
        raise SystemExit(f"{src} holds nothing but the nucleophile.")
    h04, clearance, key, phi = place_h04(zn, o01, h03, others)
    h03_new = None
    if clearance < CLASH:
        print(f"[INFO] Every position on the cone clashes - the best is "
              f"{clearance:.2f} A from {':'.join(key)}, inside the {CLASH} A "
              "cutoff. The hydroxide is pressed against its surroundings, so "
              "reorienting the whole water with the oxygen fixed.")
        found = reorient_water(zn, o01, h03, others, CLASH)
        if found is None:
            raise SystemExit(
                f"no orientation of a water at O01 clears {CLASH} A. The "
                "nucleophile position itself has no room for two hydrogens in "
                "this frame, so it needs relaxing before the proton is added. "
                "Nothing was written.")
        _, clearance, h03_new, h04 = found
        turn = math.degrees(math.acos(max(-1.0, min(1.0, dot(
            unit(sub(h03_new, o01)), unit(sub(h03, o01)))))))
        print(f"[INFO] Water reoriented: H03 turned {turn:.0f} deg, clearance "
              f"now {clearance:.2f} A. O01 is untouched, which is the only "
              "coordinate the comparison rests on.")
        key = min((norm(sub(h04, xyz)), k) for k, xyz in others)[1]
    closest = (clearance, key)

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
            if h03_new is not None and line[12:16].strip() == "H03":
                line = (line[:30]
                        + f"{h03_new[0]:8.3f}{h03_new[1]:8.3f}{h03_new[2]:8.3f}"
                        + line[54:])
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
