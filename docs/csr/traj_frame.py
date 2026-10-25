#!/usr/bin/env python3
"""One frame of a multi-frame XYZ, written back onto a PDB template.

    python traj_frame.py <traj.xyz> <template.pdb> <out.pdb> [index]

`index` is a Python index into the frames, default -1 (the last). The template
supplies everything a PDB has and an XYZ does not - names, residues, chains,
occupancies - so it must be the same system in the same atom order, which is
checked rather than assumed: the atom count must match and every element must
agree with the template's atom name.

Written for geomeTRIC's `geometric_OPTtraj_Full.xyz`, which holds the whole
system at every optimisation step and reaches hundreds of megabytes, so the
file is streamed once and only the frame asked for is kept.
"""

from __future__ import annotations

import sys
from collections import deque


def frames(path):
    """Yield (comment, [(element, x, y, z)]) for each frame, one at a time."""
    with open(path) as fh:
        while True:
            head = fh.readline()
            if not head.strip():
                return
            n = int(head.split()[0])
            comment = fh.readline().rstrip("\n")
            rows = []
            for _ in range(n):
                parts = fh.readline().split()
                rows.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
            yield comment, rows


def pick(path, index: int):
    """The frame at `index`, counting frames without holding more than needed."""
    if index < 0:
        kept = deque(maxlen=-index)
        total = 0
        for frame in frames(path):
            kept.append(frame)
            total += 1
        if len(kept) < -index:
            raise IndexError(f"{path} has {total} frames; index {index} is out of range")
        return kept[0], total
    for i, frame in enumerate(frames(path)):
        if i == index:
            return frame, None
    raise IndexError(f"{path} has fewer than {index + 1} frames")


def energy_of(comment: str):
    words = comment.split()
    for key in ("Energy", "energy"):
        if key in words:
            try:
                return float(words[words.index(key) + 1])
            except (IndexError, ValueError):
                return None
    return None


def write_pdb(template: str, rows, out: str) -> int:
    """`template` with its coordinates replaced by `rows`. Returns atoms written."""
    atoms = [l for l in open(template) if l.startswith(("ATOM  ", "HETATM"))]
    if len(atoms) != len(rows):
        raise ValueError(f"{template} has {len(atoms)} atoms, the frame has {len(rows)}")
    for i, (line, row) in enumerate(zip(atoms, rows)):
        name = line[12:16].strip()
        el = (line[76:78].strip() or name[:1]).upper()
        if row[0].upper() not in (el, name[:1].upper()):
            raise ValueError(f"atom {i + 1}: frame says {row[0]}, template says {name} ({el}). "
                             "Different atom order - this is not the same system.")
    written = 0
    with open(out, "w") as fh:
        for line in open(template):
            if line.startswith(("ATOM  ", "HETATM")):
                _, x, y, z = rows[written]
                fh.write(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")
                written += 1
            else:
                fh.write(line)
    return written


def main(argv) -> int:
    if len(argv) < 4:
        print(__doc__)
        return 2
    traj, template, out = argv[1], argv[2], argv[3]
    index = int(argv[4]) if len(argv) > 4 else -1
    (comment, rows), total = pick(traj, index)
    e = energy_of(comment)
    print(f"[INFO] {traj}: frame {index}" + (f" of {total}" if total else "")
          + (f", energy {e}" if e is not None else ", no energy in the comment line"))
    n = write_pdb(template, rows, out)
    print(f"[INFO] wrote {out}: {n} atoms, coordinates from the frame, everything else "
          f"from {template}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
