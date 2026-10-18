"""Zinc coordination sphere and active-site protonation in one or more PDBs.

    python zn_sphere.py structure1.pdb [structure2.pdb ...]

Reports, for each structure: every N/O/S within 2.8 A of the Zn; the attack
geometry of each candidate nucleophile; and how many hydrogens sit on O01,
HOH485 and Lys226 NZ, which is what distinguishes hydroxide from water and
neutral lysine from ammonium.
"""
import math
import sys

# CSR_crystal_v2 selectors (chain:resid:atom)
ZN = ("B", "1", "ZN")
C4 = ("C", "1", "C4x")
O2 = ("C", "1", "O2x")
NUCS = {"O01 (Zn-bound, ex-HOH493)": ("C", "2", "O01"),
        "HOH485": ("D", "1", "O")}
PROTON_SITES = {"O01": ("C", "2", "O01"), "HOH485 O": ("D", "1", "O"),
                "Lys226 NZ": ("A", "226", "NZ")}
SPHERE_EXTRA = {"Asp92 OD2": ("A", "92", "OD2"), "mer O2x": ("C", "1", "O2x")}


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


def dist(at, i, j):
    return math.dist(at[i][4], at[j][4])


def angle(at, i, j, k):
    v1 = [at[i][4][t] - at[j][4][t] for t in range(3)]
    v2 = [at[k][4][t] - at[j][4][t] for t in range(3)]
    cos = sum(p * q for p, q in zip(v1, v2)) / (math.hypot(*v1) * math.hypot(*v2))
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))


for path in sys.argv[1:]:
    at = load(path)
    zn = find(at, ZN)
    print("==", path)
    if zn is None:
        print("   no Zn found")
        continue
    sphere = sorted((dist(at, zn, i), i) for i, a in enumerate(at)
                    if a[2][0] in "ONS" and i != zn and dist(at, zn, i) < 2.8)
    print("   Zn sphere (<2.8 A), %d ligands:" % len(sphere))
    for d, i in sphere:
        print("      %5.2f  %s:%s:%s:%s" % (d, at[i][0], at[i][1], at[i][3], at[i][2]))
    for label, key in SPHERE_EXTRA.items():
        i = find(at, key)
        if i is not None:
            print("   Zn-%-10s %5.2f" % (label, dist(at, zn, i)))
    c4, o2 = find(at, C4), find(at, O2)
    for label, key in NUCS.items():
        i = find(at, key)
        if i is not None and c4 is not None:
            print("   %-28s O-C4x %5.2f   angle O-C4x-O2x %3.0f"
                  % (label, dist(at, i, c4), angle(at, i, c4, o2)))
    for label, key in PROTON_SITES.items():
        i = find(at, key)
        if i is None:
            continue
        hs = [j for j, a in enumerate(at) if a[2][0] == "H" and dist(at, i, j) < 1.3]
        names = ", ".join("%s(%.2f)" % (at[j][2], dist(at, i, j)) for j in hs)
        print("   %-12s %d H   %s" % (label, len(hs), names))
