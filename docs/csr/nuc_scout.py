"""Find the viable nucleophile and proton-acceptor pair in one optimised structure.

    python nuc_scout.py reactant.pdb

For every water oxygen near the electrophilic carbon, reports the attack geometry
and, for each of its hydrogens, every acceptor within reach. The point is to choose
scan coordinates from the structure that will actually be scanned, rather than from
a crystal structure or an MD frame: a second coordinate that starts at 4 A is not a
proton transfer, it drags a group across space.

A usable pair needs:
  * O...C within about 4.5 A and an O-C=O angle in the Burgi-Dunitz window (95-120),
  * and one of that water's own hydrogens already pointing at a base (< ~2.5 A).
If no single water satisfies both, say so: that is a result about the mechanism,
not a reason to force a grid.
"""
import math
import sys

C4 = ("C", "1", "C4x")
O2 = ("C", "1", "O2x")
N1 = ("C", "1", "N1x")
ZN = ("B", "1", "ZN")
BASES = {"Asp92 OD1": ("A", "92", "OD1"), "Asp92 OD2": ("A", "92", "OD2"),
         "mer N1x": N1, "mer O4x": ("C", "1", "O4x")}
MAX_O_C = 6.0
BD_LO, BD_HI = 95.0, 120.0


def load(path):
    at = []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            at.append((line[21], line[22:26].strip(), line[12:16].strip(),
                       line[17:20].strip(),
                       (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
    return at


def find(at, key):
    for i, a in enumerate(at):
        if (a[0], a[1], a[2]) == key:
            return i
    return None


def d(at, i, j):
    return math.dist(at[i][4], at[j][4])


def ang(at, i, j, k):
    v1 = [at[i][4][t] - at[j][4][t] for t in range(3)]
    v2 = [at[k][4][t] - at[j][4][t] for t in range(3)]
    c = sum(p * q for p, q in zip(v1, v2)) / (math.hypot(*v1) * math.hypot(*v2))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


for path in sys.argv[1:]:
    at = load(path)
    c4, o2, zn = find(at, C4), find(at, O2), find(at, ZN)
    print("==", path)
    cands = []
    for i, a in enumerate(at):
        if a[3] not in ("HOH", "WAT", "OH") or a[2][0] != "O":
            continue
        oc = d(at, i, c4)
        if oc > MAX_O_C:
            continue
        cands.append((oc, i))
    cands.sort()
    print("   %-18s %6s %7s %8s   %s" % ("water", "O-C4x", "angle", "Zn-O", "verdict"))
    viable = []
    for oc, i in cands:
        a = ang(at, i, c4, o2)
        zo = d(at, zn, i) if zn is not None else float("nan")
        tag = "%s:%s:%s" % (at[i][0], at[i][1], at[i][3])
        why = []
        if zo < 2.5:
            why.append("Zn-bound")
        if not (BD_LO <= a <= BD_HI):
            why.append("angle off")
        if oc > 4.5:
            why.append("far")
        verdict = "VIABLE" if not why else ", ".join(why)
        print("   %-18s %6.2f %6.0f %8.2f   %s" % (tag, oc, a, zo, verdict))
        if not why:
            viable.append(i)
    if not viable:
        print("   -> no water satisfies both distance and angle")
    # For every candidate (viable or not), where could its protons go?
    for oc, i in cands:
        hs = [j for j, b in enumerate(at) if b[2][0] == "H" and d(at, i, j) < 1.3]
        if not hs:
            continue
        print("   %s:%s:%s protons:" % (at[i][0], at[i][1], at[i][3]))
        for h in hs:
            acc = []
            for label, key in BASES.items():
                k = find(at, key)
                if k is not None:
                    acc.append((d(at, h, k), label))
            for k, b in enumerate(at):
                if b[2][0] in "ON" and k != i and abs(k - i) > 2 and d(at, h, k) < 2.6:
                    acc.append((d(at, h, k), "%s:%s:%s:%s" % (b[0], b[1], b[3], b[2])))
            acc.sort()
            shown = [x for x in acc if x[0] < 4.0][:4]
            print("      %-6s %s" % (at[h][2],
                  ", ".join("%s %.2f" % (lbl, dd) for dd, lbl in shown) or "nothing within 4 A"))
