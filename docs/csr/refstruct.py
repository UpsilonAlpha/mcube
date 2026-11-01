"""Measure a deposited structure, so a reference number can be re-derived.

    python refstruct.py shell    <pdb> [pdb ...]      what is on each metal
    python refstruct.py centre   <pdb> <resname> <resseq>
                                        a transition-state analogue's geometry
    python refstruct.py transfer <template.pdb> <target.pdb> --ligand=X --zn=...
                                 --fit=a=b,... [--poses=p.pdbqt] [--out=f.pdb]
                                        put a deposited complex's ligand into
                                        our frame and compare our poses to it

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



def select(atoms, spec):
    """The one atom a `chain:resseq:name` selector names, or an error."""
    chain, seq, name = spec.split(":")
    hit = [a for a in atoms if a[2] == chain and a[3] == seq and a[0] == name]
    if len(hit) != 1:
        raise SystemExit("%r matches %d atoms, not 1" % (spec, len(hit)))
    return hit[0]


def kabsch(mobile, fixed):
    """(rotation, mobile centroid, fixed centroid, RMSD) fitting mobile onto fixed.

    numpy only for the SVD. Everything else here is deliberately plain, but
    writing a stable SVD by hand to keep that promise would be worse.
    """
    import numpy as np

    P = np.array(mobile, dtype=float)
    Q = np.array(fixed, dtype=float)
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    moved = (R @ (P - pc).T).T + qc
    rmsd = float(np.sqrt(((moved - Q) ** 2).sum(axis=1).mean()))
    return R, pc, qc, rmsd


def apply_fit(fit, xyz):
    import numpy as np

    R, pc, qc = fit[0], fit[1], fit[2]
    return (R @ (np.array(xyz, dtype=float) - pc)) + qc


def carboxylates(lig):
    """Carbons in `lig` carrying two oxygens within 1.75 A, with those oxygens.

    Same connectivity rule as `mode_centre`, and for the same reason: the atom
    names differ in every ligand code (LMP, DQM, MX0, MER) and a name table
    would be four chances to map the wrong atom.
    """
    out = []
    for a in lig:
        if a[7] != "C":
            continue
        oxy = [o for o in lig if o[7] == "O" and dist(a, o) < 1.75]
        if len(oxy) >= 2:
            out.append((a, oxy))
    return out


def amide_carbonyls(lig):
    """Carbons bonded to a nitrogen and to a terminal oxygen - the intact
    beta-lactam carbonyl, and nothing in a hydrolysed product, where that bond
    is the one the enzyme broke."""
    out = []
    for a in lig:
        if a[7] != "C":
            continue
        nitro = [n for n in lig if n[7] == "N" and dist(a, n) < 1.6]
        oxy = [o for o in lig if o[7] == "O" and dist(a, o) < 1.45
               and len([h for h in lig if h is not o and dist(o, h) < 1.75]) == 1]
        if nitro and oxy:
            out.append((a, oxy[0], nitro[0]))
    return out


def ligand_report(label, lig, zn, reference=None):
    """What this ligand presents to `zn`, by connectivity rather than by name."""
    print("  %s: %d atoms" % (label, len(lig)))
    near = sorted((dist(zn, a), a) for a in lig)[:4]
    for d, a in near:
        print("     Zn-%-4s %5.2f" % (a[0], d))
    for c, oxy in carboxylates(lig):
        ds = sorted(dist(zn, o) for o in oxy)
        print("     carboxylate at %-4s: Zn-O %s"
              % (c[0], ", ".join("%.2f" % x for x in ds)))
    # Which amide is the beta-lactam: in a carbapenem its carbonyl carbon is a
    # few angstroms from the C3 carboxylate on the same fused bicycle, while
    # meropenem's dimethylcarbamoyl amide is at the far end of the pyrrolidine
    # arm. A pdbqt pose has no usable atom names, so the ring has to be found
    # by distance or not at all.
    acids = [c for c, _ in carboxylates(lig)]
    amides = amide_carbonyls(lig)
    # 6 A gates it: in a hydrolysed product the only amide left is the distal
    # one (7.0 A from a carboxylate C in LMP), and marking that the beta-lactam
    # would label the one bond the enzyme has already broken.
    ring = min(((min((dist(c, a) for a in acids), default=float("inf")), c, o)
                for c, o, _ in amides), default=None)
    if ring and ring[0] > 6.0:
        ring = None
    # Abstain where the rule cannot be trusted: moxalactam carries a second
    # carboxylate on its side arm, and the nearest-carboxylate test then picks
    # the distal amide (6.52 A from the zinc) over the lactam the crystal puts
    # at 2.53 A. A wrong label is worse than none.
    if len(acids) != 1:
        ring = None
    for c, o, n in amides:
        near = min((dist(c, a) for a in acids), default=float("nan"))
        mark = ("  <-- beta-lactam" if ring and c is ring[1]
                else ("  (%d carboxylates: not marked)" % len(acids)
                      if len(acids) != 1 else ""))
        print("     amide C=O at %-4s: Zn-O %.2f (N %s, %.1f A from a "
              "carboxylate C)%s" % (c[0], dist(zn, o), n[0], near, mark))
    if reference:
        print("     centroid %.2f A from %s's"
              % (centroid_distance(lig, reference[1]), reference[0]))


def centroid(lig):
    n = float(len(lig))
    return [sum(a[4 + i] for a in lig) / n for i in range(3)]


def centroid_distance(a, b):
    ca, cb = centroid(a), centroid(b)
    return math.sqrt(sum((ca[i] - cb[i]) ** 2 for i in range(3)))


def models(path, resname=None):
    """[(label, atoms)] - one per MODEL for a multi-model file, else one entry.

    Reads .pdbqt as well as .pdb: the coordinate columns are the same and the
    extra trailing columns are never looked at.
    """
    out, current, label = [], [], None
    for line in open(path, errors="replace"):
        if line.startswith("MODEL"):
            current, label = [], "model %s" % line.split()[-1]
        elif line.startswith("ENDMDL"):
            out.append((label, current))
        elif line.startswith(("ATOM", "HETATM")) and line[16] in (" ", "A"):
            element = line[76:78].strip().upper()
            if not element:                      # pdbqt: its own type column
                element = line[12:16].strip()[0]
            current.append((line[12:16].strip(), line[17:20].strip(), line[21],
                            line[22:26].strip(),
                            float(line[30:38]), float(line[38:46]),
                            float(line[46:54]), element))
    if not out:
        out = [("", current)]
    if resname:
        out = [(lab, [a for a in at if a[1] == resname]) for lab, at in out]
    return [(lab, at) for lab, at in out if at]


def pdb_lines(lig, resname="LIG", chain="L", seq=1):
    for i, a in enumerate(lig, 1):
        yield ("HETATM%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
               % (i, a[0][:4], resname, chain, seq, a[4], a[5], a[6], a[7][:2]))


def mode_transfer(args):
    """Put a deposited complex's ligand into our frame, and compare poses to it.

        python refstruct.py transfer <template.pdb> <target.pdb> \
            --ligand=LMP --zn=B:1:ZN \
            --fit=A:72:NE2=A:88:NE2,A:150:NE2=A:164:NE2,... \
            [--poses=meropenem.pdbqt] [--out=transferred.pdb]

    The fit is a local one on the zinc site, not a fold alignment: the pairs
    are given explicitly because the two proteins have different numbering and
    guessing the correspondence is exactly the step that must not be guessed.
    The fit RMSD is printed, and a fit that does not close to a few tenths of
    an angstrom is a correspondence error rather than a result.
    """
    opts = {k.lstrip("-").split("=")[0]: k.split("=", 1)[1]
            for k in args if k.startswith("--")}
    paths = [a for a in args if not a.startswith("--")]
    template, target = paths[0], paths[1]
    tat, qat = rows(template), rows(target)

    pairs = [p.split("=") for p in opts["fit"].split(",")]
    mobile = [select(tat, t)[4:7] for t, _ in pairs]
    fixed = [select(qat, q)[4:7] for _, q in pairs]
    fit = kabsch(mobile, fixed)
    print("== %s -> %s" % (template, target))
    print("  fit on %d atoms, RMSD %.3f A" % (len(pairs), fit[3]))
    for (t, q), m, f in zip(pairs, mobile, fixed):
        moved = apply_fit(fit, m)
        print("     %-14s -> %-14s %.3f A"
              % (t, q, math.sqrt(sum((moved[i] - f[i]) ** 2 for i in range(3)))))

    zn = select(qat, opts["zn"])
    lig = [a for a in tat if a[1] == opts["ligand"]]
    if not lig:
        raise SystemExit("no %s in %s" % (opts["ligand"], template))
    moved = [a[:4] + tuple(apply_fit(fit, a[4:7])) + (a[7],) for a in lig]
    ligand_report("%s transferred" % opts["ligand"], moved, zn)

    if opts.get("out"):
        with open(opts["out"], "w") as fh:
            fh.writelines(pdb_lines(moved, opts["ligand"]))
            fh.write("END\n")
        print("  written to %s" % opts["out"])

    for path in (opts.get("poses") or "").split(",") if opts.get("poses") else []:
        print("  -- %s" % path)
        for label, pose in models(path, opts.get("poseres")):
            ligand_report(label or path, pose, zn,
                          ("the transferred ligand", moved))


if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1] not in ("shell", "centre", "transfer"):
        raise SystemExit(__doc__)
    if sys.argv[1] == "shell":
        mode_shell(sys.argv[2:])
    elif sys.argv[1] == "transfer":
        mode_transfer(sys.argv[2:])
    else:
        if len(sys.argv) != 5:
            raise SystemExit("centre needs <pdb> <resname> <resseq>")
        mode_centre(sys.argv[2:])
