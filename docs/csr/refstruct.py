"""Measure a deposited structure, so a reference number can be re-derived.

    python refstruct.py shell    <pdb> [pdb ...]      what is on each metal
    python refstruct.py centre   <pdb> <resname> <resseq>
                                        a transition-state analogue's geometry

Every reference value quoted in `site.yaml`'s `michaelis_criteria` and in the
decision log came out of these two commands, and `reference-structures.md`
records the exact invocations. They are here rather than in a scratch directory
because the numbers end up in a paper: a criterion that cannot be re-derived
from a PDB ID is an assertion, not a measurement.

Deliberately independent of `site.py` and of `site.yaml`. These read deposited
structures with their own numbering and their own ligand codes, so there is no
site description to drive them, and nothing here should be able to inherit an
assumption from our own system's description.

Alternate conformations: only blank and A are read, so a partially occupied
side chain contributes one position rather than two.
"""
import math
import sys

AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER "
         "THR TRP TYR VAL MSE".split())


def rows(path):
    """(name, resname, chain, resseq, x, y, z, element) for each atom."""
    out = []
    for line in open(path, errors="replace"):
        if line.startswith(("ATOM", "HETATM")) and line[16] in (" ", "A"):
            out.append((line[12:16].strip(), line[17:20].strip(), line[21],
                        line[22:26].strip(),
                        float(line[30:38]), float(line[38:46]), float(line[46:54]),
                        line[76:78].strip().upper()))
    return out


def dist(a, b):
    return math.sqrt(sum((a[4 + i] - b[4 + i]) ** 2 for i in range(3)))


def angle(a, b, c):
    u = [a[4 + i] - b[4 + i] for i in range(3)]
    v = [c[4 + i] - b[4 + i] for i in range(3)]
    nu = math.sqrt(sum(x * x for x in u))
    nv = math.sqrt(sum(x * x for x in v))
    cos = sum(u[i] * v[i] for i in range(3)) / (nu * nv)
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))


def mode_shell(paths):
    """Donors within 3.2 A of every zinc, and any metal-metal separation.

    3.2 A rather than our own 2.6 A cutoff: the point of reading a deposited
    structure is to see what is there, including the borderline contact that a
    tighter cutoff would hide.
    """
    for path in paths:
        at = rows(path)
        zns = [a for a in at if a[7] == "ZN"]
        het = sorted({a[1] for a in at
                      if a[1] != "HOH" and a[7] != "ZN" and a[1] not in AA})
        print("\n==", path, "  %d Zn, ligands present: %s"
              % (len(zns), ",".join(het)[:120]))
        for z in zns:
            near = sorted((dist(z, a), a) for a in at
                          if a[7] in ("O", "N", "S") and dist(z, a) < 3.2)
            print("  Zn %s%s chain %s:" % (z[1], z[3], z[2]))
            for d, a in near:
                mark = "  <-- LIGAND" if a[1] in het else ""
                print("     %5.2f  %-4s %-3s %s%-4s%s"
                      % (d, a[0], a[1], a[2], a[3], mark))
        for i, z1 in enumerate(zns):
            for z2 in zns[i + 1:]:
                if dist(z1, z2) < 6:
                    print("  Zn-Zn %.2f (chain %s/%s)"
                          % (dist(z1, z2), z1[2], z2[2]))


def mode_centre(args):
    """A transition-state analogue's tetrahedral centre: its oxygens and angles.

    The centre is found by connectivity rather than by name - the P or C in the
    named residue that carries at least two oxygens within 1.75 A - so a
    phosphonate, a carboxylate and a tetrahedral carbon are all handled without
    being told which is which.
    """
    path, resn, seq = args[0], args[1], args[2]
    at = rows(path)
    zn = [a for a in at if a[7] == "ZN"]
    lig = [a for a in at if a[1] == resn and a[3] == seq]
    if not lig:
        raise SystemExit("no residue %s %s in %s" % (resn, seq, path))
    print("==", path, resn, seq)
    centres = [a for a in lig if a[7] in ("P", "C")
               and len([o for o in lig if o[7] == "O" and dist(a, o) < 1.75]) >= 2]
    for c in centres:
        oxy = sorted((dist(c, o), o) for o in lig if o[7] == "O" and dist(c, o) < 1.75)
        print("  centre %s: %d oxygens" % (c[0], len(oxy)))
        for dc, o in oxy:
            zd = min((dist(z, o) for z in zn), default=float("nan"))
            print("     %-4s  %s-O %.2f   Zn-O %.2f" % (o[0], c[0], dc, zd))
        for i in range(len(oxy)):
            for j in range(i + 1, len(oxy)):
                a, b = oxy[i][1], oxy[j][1]
                print("     angle %s-%s-%s = %.0f deg"
                      % (a[0], c[0], b[0], angle(a, c, b)))
        if len(oxy) >= 2:
            print("     O...O %.2f" % dist(oxy[0][1], oxy[1][1]))


if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1] not in ("shell", "centre"):
        raise SystemExit(__doc__)
    if sys.argv[1] == "shell":
        mode_shell(sys.argv[2:])
    else:
        if len(sys.argv) != 5:
            raise SystemExit("centre needs <pdb> <resname> <resseq>")
        mode_centre(sys.argv[2:])
