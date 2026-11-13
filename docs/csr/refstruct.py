"""Measure a deposited structure, so a reference number can be re-derived.

    python refstruct.py shell    <pdb> [pdb ...]      what is on each metal
    python refstruct.py centre   <pdb> <resname> <resseq>
                                        a transition-state analogue's geometry
    python refstruct.py match <ours.pdb> <template.pdb> [--out=f.pdb]
                                        our atom names on a deposited ligand's
                                        coordinates, matched by bonding
    python refstruct.py pair <pdb> <sel>-<sel> [<sel>-<sel> ...]
                                        named distances, e.g. A:91:OD2-A:87:NE2
    python refstruct.py similarity <a.pdb>[:CHAIN] <b.pdb>[:CHAIN] ...
                                        identity and CA RMSD for every pair
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



def mode_pair(args):
    """Distances between named atoms of one deposited structure.

    `shell` answers what is on a metal, which needs a metal; this answers what
    is near what, which does not. It exists because the reference a wild-type
    model should be compared against may have no metal at all - CSR-1's own
    wild-type crystals (6DN4, 6DQ2) are metal-free, so the only way to ask
    where their Asp92 sits is against its own neighbours.

    An internal distance also needs no superposition, which is the point: two
    structures can be compared on it without a fit, so the answer cannot be an
    artefact of how they were aligned.
    """
    if len(args) < 2:
        raise SystemExit("pair needs <pdb> <chain:resseq:name>-<chain:resseq:name> ...")
    atoms = rows(args[0])
    print("%s (%d atoms read)" % (args[0], len(atoms)))
    for spec in args[1:]:
        if spec.count("-") != 1:
            raise SystemExit("%r is not <sel>-<sel>" % spec)
        left, right = spec.split("-")
        a, b = select(atoms, left), select(atoms, right)
        print("  %-22s %-22s %6.2f A   (%s %s / %s %s)"
              % (left, right, dist(a, b), a[1], a[3], b[1], b[3]))


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

    `--ligand` takes a residue name, or `CHAIN:RESSEQ[:RESNAME]` where the name
    is ambiguous - 6NC5 has six residues called ZN and only one of them is the
    site. `--zn` is optional: the target may be the metal-free structure, which
    is the case that makes this worth doing.

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

    # `--ligand` by residue NAME, or by `CHAIN:RESSEQ[:RESNAME]` when the name
    # alone is ambiguous. 6NC5 holds six zincs, all called ZN, and only ZN302
    # is the site we model - so a name is not enough to name a metal (¶184).
    want = opts["ligand"]
    if ":" in want:
        bits = want.split(":")
        chain, seq = bits[0], bits[1]
        lig = [a for a in tat if a[2] == chain and a[3] == seq
               and (len(bits) < 3 or a[1] == bits[2])]
    else:
        lig = [a for a in tat if a[1] == want]
    if not lig:
        raise SystemExit("no %s in %s" % (want, template))
    moved = [a[:4] + tuple(apply_fit(fit, a[4:7])) + (a[7],) for a in lig]

    # `--zn` is optional, because the structure worth transferring INTO may be
    # the one with no metal: CSR-1's wild-type crystals are metal-free, which
    # is the whole reason to move 6NC5's zinc into 6DN4's frame. Without it the
    # group is written out and measured with `pair`, which needs no metal.
    zn = select(qat, opts["zn"]) if opts.get("zn") else None
    if zn is not None:
        ligand_report("%s transferred" % want, moved, zn)
    else:
        print("  no --zn in the target: distances are not reported here. "
              "Write the group with --out, append it to the target, and "
              "measure with `pair` - which needs no metal.")

    if opts.get("out"):
        with open(opts["out"], "w") as fh:
            # The residue's own name, not `--ligand`: that may now be a
            # selector, and writing `A:302:ZN` into the four-character resname
            # column produces a file every other tool misparses (¶184).
            fh.writelines(pdb_lines(moved, moved[0][1]))
            fh.write("END\n")
        print("  written to %s" % opts["out"])

    if opts.get("poses") and zn is None:
        print("  --poses needs --zn: a pose is compared by how it sits on the "
              "metal, and the target has none.")
    for path in (opts.get("poses") or "").split(",") if opts.get("poses") and zn else []:
        print("  -- %s" % path)
        for label, pose in models(path, opts.get("poseres")):
            ligand_report(label or path, pose, zn,
                          ("the transferred ligand", moved))



# BLOSUM62, the 20 standard residues. Embedded rather than imported so this
# file keeps working in an environment without biopython, which is every
# environment this project has.
_B62_ORDER = "ARNDCQEGHILKMFPSTWYV"
_B62_ROWS = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
BLOSUM62 = {}
for _i, _row in enumerate(_B62_ROWS.strip().splitlines()):
    for _j, _v in enumerate(_row.split()):
        BLOSUM62[(_B62_ORDER[_i], _B62_ORDER[_j])] = int(_v)

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "MSE": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T",
    "TRP": "W", "TYR": "Y", "VAL": "V",
    # the protonation variants our own models carry
    "HID": "H", "HIE": "H", "HIP": "H", "ASH": "D", "GLH": "E", "LYN": "K",
    "CYX": "C",
}


def chain_residues(path, chain=None):
    """[(one-letter, CA xyz)] for one chain, in file order.

    Residues without a CA are dropped: they cannot be superposed and would
    make the sequence and the structure disagree about what was aligned.
    """
    seen, out = set(), []
    for a in rows(path):
        if a[0] != "CA" or a[1] not in THREE_TO_ONE:
            continue
        if chain and a[2] != chain:
            continue
        key = (a[2], a[3])
        if key in seen:
            continue
        seen.add(key)
        out.append((THREE_TO_ONE[a[1]], (a[4], a[5], a[6])))
    return out


def align(seq_a, seq_b, gap_open=-11, gap_extend=-1):
    """Global Needleman-Wunsch with affine gaps. [(i, j)] of aligned pairs.

    Affine rather than linear because a linear penalty tears loops apart to
    buy single matches, and the loops are exactly where two MBLs of 20-40%
    identity differ.
    """
    import numpy as np

    n, m = len(seq_a), len(seq_b)
    neg = -1e9
    M = np.full((n + 1, m + 1), neg)      # ends aligned
    X = np.full((n + 1, m + 1), neg)      # gap in b
    Y = np.full((n + 1, m + 1), neg)      # gap in a
    M[0, 0] = 0.0
    for i in range(1, n + 1):
        X[i, 0] = gap_open + (i - 1) * gap_extend
    for j in range(1, m + 1):
        Y[0, j] = gap_open + (j - 1) * gap_extend
    ptr = {}
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = BLOSUM62.get((seq_a[i - 1], seq_b[j - 1]), -1)
            best = max((M[i - 1, j - 1], "M"), (X[i - 1, j - 1], "X"),
                       (Y[i - 1, j - 1], "Y"))
            M[i, j] = best[0] + s
            ptr[("M", i, j)] = best[1]
            bx = max((M[i - 1, j] + gap_open, "M"), (X[i - 1, j] + gap_extend, "X"))
            X[i, j], ptr[("X", i, j)] = bx[0], bx[1]
            by = max((M[i, j - 1] + gap_open, "M"), (Y[i, j - 1] + gap_extend, "Y"))
            Y[i, j], ptr[("Y", i, j)] = by[0], by[1]
    state = max((M[n, m], "M"), (X[n, m], "X"), (Y[n, m], "Y"))[1]
    i, j, pairs = n, m, []
    while i > 0 and j > 0:
        nxt = ptr[(state, i, j)]
        if state == "M":
            pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif state == "X":
            i -= 1
        else:
            j -= 1
        state = nxt
    return list(reversed(pairs))


def superpose_pairs(a, b, pairs, cutoff=3.0, rounds=10):
    """Fit b's CAs onto a's over `pairs`, dropping outliers and refitting.

    Returns (all-pair RMSD, core size, core RMSD). The core is what a
    structural superposition actually reports: two enzymes of this family share
    a fold and differ in loops, and one flapping loop dominates an all-pair
    RMSD while saying nothing about whether the active sites agree.
    """
    import numpy as np

    keep = list(pairs)
    first = None
    for _ in range(rounds):
        P = [b[j][1] for _, j in keep]
        Q = [a[i][1] for i, _ in keep]
        R, pc, qc, rmsd = kabsch(P, Q)
        moved = np.array([apply_fit((R, pc, qc), p) for p in P])
        d = np.sqrt(((moved - np.array(Q)) ** 2).sum(axis=1))
        if first is None:
            first = float(np.sqrt((d ** 2).mean()))
        # Tighten gradually: one badly displaced loop drags the whole fit, so
        # a hard 3 A cut in the first round throws away residues that are only
        # misplaced because of it. Dropping at twice the current RMSD removes
        # the real outliers first and converges on the cutoff.
        limit = max(cutoff, 2.0 * rmsd)
        nxt = [pr for pr, dd in zip(keep, d) if dd <= limit]
        floor = min(20, len(pairs) // 2)
        if len(nxt) == len(keep) or len(nxt) < floor:
            return first, len(keep), rmsd
        keep = nxt
    return first, len(keep), rmsd


def mode_similarity(args):
    """How alike are two structures, in sequence and in fold?

        python refstruct.py similarity <a.pdb>[:CHAIN] <b.pdb>[:CHAIN] ...

    Every pair of the given structures is reported: identity over the aligned
    positions, then a CA superposition over those same positions with outliers
    dropped. The point is to answer "would these two hold a substrate the same
    way" with a number rather than with a subclass label.
    """
    loaded = []
    for spec in args:
        # rpartition, and only a single character counts as the chain: a
        # Windows path starts "C:/..." and partition() would read the drive
        # letter as the chain and the rest as the file.
        head, sep, tail = spec.rpartition(":")
        path, chain = (head, tail) if sep and len(tail) == 1 else (spec, "")
        res = chain_residues(path, chain or None)
        if not res:
            raise SystemExit("no residues with a CA in %s" % spec)
        loaded.append((spec, res))
        print("  %-28s %d residues" % (spec, len(res)))
    for x in range(len(loaded)):
        for y in range(x + 1, len(loaded)):
            (na, a), (nb, b) = loaded[x], loaded[y]
            pairs = align("".join(r[0] for r in a), "".join(r[0] for r in b))
            same = sum(1 for i, j in pairs if a[i][0] == b[j][0])
            ident = 100.0 * same / len(pairs) if pairs else 0.0
            full, core, core_rmsd = superpose_pairs(a, b, pairs)
            print("\n== %s  vs  %s" % (na, nb))
            print("   aligned %d positions, %d identical = %.1f%% identity"
                  % (len(pairs), same, ident))
            print("   CA RMSD over all aligned: %.2f A" % full)
            print("   core: %d of %d within 3.0 A, RMSD %.2f A"
                  % (core, len(pairs), core_rmsd))



# Covalent cutoffs, generous enough for a crystal ligand at 1.4-1.6 A
# resolution and tight enough not to bond across a ring.
_BOND_MAX = {("S", "C"): 1.95, ("S", "O"): 1.75, ("S", "N"): 1.85}
_BOND_DEFAULT = 1.75


def bond_max(e1, e2):
    return _BOND_MAX.get((e1, e2)) or _BOND_MAX.get((e2, e1)) or _BOND_DEFAULT



def conect_bonds(path):
    """{(i, j)} over ATOM/HETATM order, from the file's own CONECT records.

    None when the file has none. Our ligands carry them - they are the bonds
    the force field was built on - and reading them is the difference between
    using the topology and re-deriving it. A distance rule gets 27 atoms and
    28 bonds right and can still put one of them between the wrong pair.
    """
    index, order, conect = {}, 0, []
    for line in open(path, errors="replace"):
        if line.startswith(("ATOM", "HETATM")) and line[16] in (" ", "A"):
            serial = line[6:11].strip()
            if serial:
                index[int(serial)] = order
            order += 1
        elif line.startswith("CONECT"):
            conect.append(line.rstrip())
    if not conect:
        return None
    pairs = set()
    for line in conect:
        fields = [line[i:i + 5].strip() for i in range(6, min(len(line), 31), 5)]
        nums = [int(f) for f in fields if f]
        if len(nums) < 2:
            continue
        a = index.get(nums[0])
        for other in nums[1:]:
            b = index.get(other)
            if a is not None and b is not None and a != b:
                pairs.add((min(a, b), max(a, b)))
    return pairs


def heavy_graph(atoms, bonds=None):
    """(heavy atoms, adjacency sets). From `bonds` if given, else by distance.

    `bonds` is {(i, j)} over the full atom order, as `conect_bonds` returns it.
    Hydrogens are left out either way: a crystal ligand has none, and the match
    only has to be unambiguous over the heavy-atom skeleton.
    """
    keep = [i for i, a in enumerate(atoms) if a[7] != "H"]
    heavy = [atoms[i] for i in keep]
    position = {i: k for k, i in enumerate(keep)}
    adj = [set() for _ in heavy]
    if bonds is not None:
        for i, j in bonds:
            if i in position and j in position:
                adj[position[i]].add(position[j])
                adj[position[j]].add(position[i])
        return heavy, adj
    for i in range(len(heavy)):
        for j in range(i + 1, len(heavy)):
            if dist(heavy[i], heavy[j]) <= bond_max(heavy[i][7], heavy[j][7]):
                adj[i].add(j)
                adj[j].add(i)
    return heavy, adj


def isomorphisms(a_el, a_adj, b_el, b_adj, limit=64):
    """Every mapping a->b preserving element and bonding, up to `limit`.

    Backtracking with the most-connected atom first, which prunes hard on a
    fused bicycle. The count matters as much as the mapping: two carboxylate
    oxygens are interchangeable, and a caller that is told there are several
    mappings knows which parts of the molecule its choice was arbitrary for.
    """
    n = len(a_el)
    if n != len(b_el) or sorted(a_el) != sorted(b_el):
        return []
    order = sorted(range(n), key=lambda i: -len(a_adj[i]))
    found, mapping, used = [], {}, set()

    def step(k):
        if len(found) >= limit:
            return
        if k == n:
            found.append(dict(mapping))
            return
        i = order[k]
        for j in range(n):
            if j in used or b_el[j] != a_el[i] or len(b_adj[j]) != len(a_adj[i]):
                continue
            ok = True
            for nb, jj in mapping.items():
                if (nb in a_adj[i]) != (jj in b_adj[j]):
                    ok = False
                    break
            if ok:
                mapping[i] = j
                used.add(j)
                step(k + 1)
                del mapping[i]
                used.discard(j)

    step(0)
    return found



def hydrogen_parents(atoms, bonds=None):
    """{hydrogen index: heavy-atom index} over the full atom order.

    From `bonds` where the file has CONECT, else the nearest heavy atom within
    1.3 A, which is every X-H bond there is.
    """
    out = {}
    if bonds is not None:
        for i, j in bonds:
            for h, heavy in ((i, j), (j, i)):
                if atoms[h][7] == "H" and atoms[heavy][7] != "H":
                    out[h] = heavy
        return out
    heavy_idx = [k for k, a in enumerate(atoms) if a[7] != "H"]
    for k, a in enumerate(atoms):
        if a[7] != "H":
            continue
        near = sorted((dist(a, atoms[m]), m) for m in heavy_idx)
        if near and near[0][0] <= 1.3:
            out[k] = near[0][1]
    return out


def local_environment(index, adj, want=4):
    """`index` plus its neighbours and their neighbours, as a set.

    Enough atoms to define a rigid frame: a hydrogen's parent and one neighbour
    are two points, and two points leave a free rotation about their axis -
    which would spin the hydrogen anywhere on a cone.
    """
    env = {index} | set(adj[index])
    if len(env) < want:
        for nb in list(adj[index]):
            env |= set(adj[nb])
    return env


def transfer_hydrogens(atoms, parents, heavy_order, mapping, template_heavy, adj):
    """{hydrogen index: new xyz}, each placed by its parent's local frame.

    The heavy atoms move to the template's coordinates, which changes the
    conformation; a hydrogen left where it was would be a bond stretched across
    the molecule. Each one is carried over by fitting the small rigid group
    around its parent - the parent, its neighbours and theirs - onto the same
    atoms in their new positions, and applying that fit.
    """
    position = {atom_index: k for k, atom_index in enumerate(heavy_order)}
    moved = {}
    for h, parent in parents.items():
        if parent not in position:
            continue
        env = [i for i in local_environment(position[parent], adj)]
        if len(env) < 3:
            env = list(range(min(3, len(heavy_order))))
        src = [atoms[heavy_order[k]][4:7] for k in env]
        dst = [template_heavy[mapping[k]][4:7] for k in env]
        fit = kabsch(src, dst)
        moved[h] = tuple(apply_fit(fit, atoms[h][4:7]))
    return moved


def mode_match(args):
    """Put our ligand's atom names on a deposited ligand's coordinates.

        python refstruct.py match <ours.pdb> <template.pdb> [--out=f.pdb]

    Both files must hold one ligand and the same molecule. The heavy atoms are
    matched by graph isomorphism - element and bonding, never by name, because
    LMP and our `hyd` name the same 27 atoms completely differently - and the
    output is **our** atom names and order at the **template's** coordinates.

    Hydrogens come across too, each carried by a rigid fit of the small group
    around its parent, because mcube parameterises a ligand from its own file
    against `AddHs(MolFromSmiles(...))` and a structure without them cannot
    match that template. Leaving them at their old coordinates is not an
    option either: the heavy atoms have moved.

    The number of distinct mappings is printed. More than one means a symmetric
    group - the two oxygens of a carboxylate, the two methyls of an amide - and
    says which names were assigned arbitrarily.
    """
    opts = {k.lstrip("-").split("=")[0]: k.split("=", 1)[1]
            for k in args if k.startswith("--")}
    paths = [a for a in args if not a.startswith("--")]
    ours, template = rows(paths[0]), rows(paths[1])
    # Our ligand's CONECT is the topology its parameters were built on, so it
    # wins over any distance rule. The template is a crystal ligand with no
    # CONECT, where distance bonding is reliable - and if the two graphs
    # disagree, no isomorphism is found, which is the check rather than a
    # silent reinterpretation.
    our_bonds = conect_bonds(paths[0])
    a_heavy, a_adj = heavy_graph(ours, our_bonds)
    b_heavy, b_adj = heavy_graph(template, conect_bonds(paths[1]))
    print("  ours %d heavy atoms (bonds from %s), template %d"
          % (len(a_heavy), "CONECT" if our_bonds else "distance", len(b_heavy)))
    if len(a_heavy) != len(b_heavy):
        raise SystemExit(
            "different heavy-atom counts (%d and %d), so these are not the "
            "same molecule. Nothing was written." % (len(a_heavy), len(b_heavy)))

    maps = isomorphisms([a[7] for a in a_heavy], a_adj,
                        [b[7] for b in b_heavy], b_adj)
    if not maps:
        raise SystemExit(
            "no mapping preserves element and bonding, so the two ligands are "
            "not the same molecule - or one of them has a bond the distance "
            "cutoff did not find. Nothing was written.")
    print("  %d mapping%s found; using the first"
          % (len(maps), "" if len(maps) == 1 else "s"))
    best = maps[0]

    # Which names the symmetry left arbitrary, if any.
    varying = sorted({a_heavy[i][0] for m in maps[1:] for i in m
                      if m[i] != best[i]})
    if varying:
        print("  interchangeable under the symmetry: %s" % ", ".join(varying))

    heavy_order = [i for i, a in enumerate(ours) if a[7] != "H"]
    placed = {}
    for i, a in enumerate(a_heavy):
        t = b_heavy[best[i]]
        print("     %-5s <- %-5s" % (a[0], t[0]))
        placed[heavy_order[i]] = (t[4], t[5], t[6])

    # mcube parameterises a ligand from its own PDB and matches it against
    # `AddHs(MolFromSmiles(...))`, so the file needs its hydrogens: without
    # them the template carries 26 atoms the structure does not and the
    # substructure match cannot succeed, whatever the protonation. Five
    # preprocess jobs died on that (¶162).
    parents = hydrogen_parents(ours, our_bonds)
    placed.update(transfer_hydrogens(ours, parents, heavy_order, best,
                                     b_heavy, a_adj))
    missing = [i for i, a in enumerate(ours) if i not in placed]
    if missing:
        raise SystemExit(
            "%d atom(s) could not be placed, including %s - every hydrogen "
            "needs a heavy-atom parent. Nothing was written."
            % (len(missing), ours[missing[0]][0]))
    moved = [(a[0], a[1], a[2], a[3]) + placed[i] + (a[7],)
             for i, a in enumerate(ours)]
    print("  %d atoms placed (%d heavy, %d hydrogens carried on their "
          "parents' frames)" % (len(moved), len(a_heavy), len(parents)))

    if opts.get("out"):
        resname = moved[0][1]
        with open(opts["out"], "w") as fh:
            fh.writelines(pdb_lines(moved, resname,
                                    moved[0][2] or "L", int(moved[0][3] or 1)))
            # CONECT, because the bonds are the whole point of this file.
            # `preprocess` reads a ligand with `proximityBonding=False`, so
            # without them RDKit sees 27 unbonded atoms and refuses the file -
            # which is what job 31299591 died of, and the right behaviour: a
            # ligand whose bonds are guessed from distance is how this project
            # got a spurious O-O bond once already (¶21).
            for i, j in sorted(our_bonds or set()):
                fh.write("CONECT%5d%5d\n" % (i + 1, j + 1))
            if our_bonds is None:
                for i in range(len(a_adj)):
                    for j in sorted(a_adj[i]):
                        if j > i:
                            fh.write("CONECT%5d%5d\n"
                                     % (heavy_order[i] + 1, heavy_order[j] + 1))
            fh.write("END\n")
        print("  written to %s (all atoms and CONECT)" % opts["out"])


if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1] not in ("shell", "centre", "transfer", "match",
                                            "similarity", "pair"):
        raise SystemExit(__doc__)
    if sys.argv[1] == "shell":
        mode_shell(sys.argv[2:])
    elif sys.argv[1] == "pair":
        mode_pair(sys.argv[2:])
    elif sys.argv[1] == "transfer":
        mode_transfer(sys.argv[2:])
    elif sys.argv[1] == "similarity":
        mode_similarity(sys.argv[2:])
    elif sys.argv[1] == "match":
        mode_match(sys.argv[2:])
    else:
        if len(sys.argv) != 5:
            raise SystemExit("centre needs <pdb> <resname> <resseq>")
        mode_centre(sys.argv[2:])
