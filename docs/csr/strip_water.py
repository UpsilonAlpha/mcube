"""Drop bulk water from a PDB, keeping the waters that coordinate a metal.

A solvated CSR frame is ~3.7 MB and 91% of its atoms are bulk water: 43,668 of
47,979 in `relax10`'s last frame. Everything this project measures - the metal,
its donors, Asp92, Lys226, the ligand - is in the other 9%, which is 341 kB. So
stripping the solvent is the difference between a file that can be moved off
Bunya through a read-only channel and one that cannot (¶251).

    python strip_water.py <in.pdb> <out.pdb> [--metal=B:1:ZN] [--shell=5.0]

Waters with any atom within `--shell` of the named metal are KEPT, because
whether a water holds the open coordination site is a result in this project
(¶191's resident waters, ¶214's four-coordinate sphere), not noise. Everything
else named HOH/WAT/SOL/TIP3 goes. CRYST1 is preserved so the box survives; so
are CONECT records for atoms that survive, since a grafted ligand without them
is unbonded to whatever reads it next (¶160).

Writes nothing if the output exists, and reports what it dropped.
"""
import math
import os
import sys

WATER = ("HOH", "WAT", "SOL", "TIP", "TIP3", "T3P", "SPC")


def _key(line):
    return (line[21:22], line[22:26].strip())


def _xyz(line):
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def _is_water(line):
    return line[17:20].strip().upper() in WATER


def _select(path, metal, shell):
    """Return (kept_lines, dropped_count, kept_water_count).

    Two passes: the first finds the metal and every water residue with an atom
    inside `shell`, the second emits. Two passes rather than one because a
    water's own atoms may appear before the metal does, and a single pass would
    have to buffer the whole file anyway.
    """
    want_chain, want_seq, want_name = (metal.split(":") + [None])[:3]
    centre, near = None, set()
    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if (line[21:22] == want_chain
                    and line[22:26].strip() == want_seq
                    and (want_name is None
                         or line[17:20].strip() == want_name
                         or line[12:16].strip() == want_name)):
                centre = _xyz(line)
                break
    if centre is None:
        raise SystemExit("no atom matches --metal=%s in %s" % (metal, path))

    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")) or not _is_water(line):
                continue
            x, y, z = _xyz(line)
            if math.dist((x, y, z), centre) <= shell:
                near.add(_key(line))

    kept, dropped, serials = [], 0, set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("CRYST1"):
                kept.append(line)
            elif line.startswith(("ATOM", "HETATM")):
                if _is_water(line) and _key(line) not in near:
                    dropped += 1
                else:
                    kept.append(line)
                    serials.add(line[6:11].strip())
            elif line.startswith("CONECT"):
                fields = [line[i:i + 5].strip() for i in range(6, len(line.rstrip()), 5)]
                if fields and all(f in serials for f in fields if f):
                    kept.append(line)
    kept.append("END\n")
    return kept, dropped, len(near)


def main(argv):
    paths = [a for a in argv if not a.startswith("--")]
    opts = dict(a.lstrip("-").split("=", 1) for a in argv if a.startswith("--"))
    if len(paths) != 2:
        raise SystemExit(__doc__)
    src, dst = paths
    if os.path.exists(dst):
        raise SystemExit("%s exists; refusing to overwrite" % dst)
    kept, dropped, near = _select(
        src, opts.get("metal", "B:1:ZN"), float(opts.get("shell", 5.0)))
    with open(dst, "w") as fh:
        fh.writelines(kept)
    print("%s -> %s" % (src, dst))
    print("  dropped %d bulk water atoms, kept %d water residues within %s A"
          % (dropped, near, opts.get("shell", 5.0)))
    print("  %.1f kB from %.1f MB"
          % (os.path.getsize(dst) / 1024, os.path.getsize(src) / 1048576))


if __name__ == "__main__":
    main(sys.argv[1:])
