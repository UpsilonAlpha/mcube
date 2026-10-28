"""Active-site geometry analysis, driven by a site description.

    python site.py sphere       <pdb> [pdb ...]     metal sphere and protonation
    python site.py nucleophiles <pdb> [pdb ...]     attack geometry + proton acceptors
    python site.py shell        <traj> <topology>   metal coordination over a trajectory
    python site.py acceptors    <traj> <topology>   what CAN take the proton, counted
    python site.py joint        <traj> <topology>   attack AND donate, same molecule
    python site.py poses        <pdbqt> <receptor> [ligand.pdb]
                                                    docked poses vs the metal site
    python site.py graft        <pdbqt> <ligand.pdb> <receptor.pdb> <system.pdb>
                                                [1,3,4]  a pose into the system
    python site.py check        <pdb> [pdb ...] [--at=start|reactant]
                                                    accept or reject a structure

Every atom, every candidate base and every criterion comes from `site.yaml`
(override with --site), so this script carries no knowledge of a particular
enzyme and a new system needs a new YAML file rather than new code.

The modes answer different questions, in rough order of cost:

  sphere        what is coordinated, and where are the hydrogens?
  shell         over a trajectory, what is coordinated and how many of them, with
                every ligand's occupancy - the question `sphere` can only
                answer for the one frame it is given.
  nucleophiles  for one structure: which waters could attack, and for each of
                their protons, what is within hydrogen-bond reach?
  acceptors     over a trajectory, during the frames where a water is genuinely
                on the Burgi-Dunitz line, what is near its protons at all? Use
                this when you do not want to assume which base is involved.
  joint         over a trajectory, does any single water both reach the attack
                geometry and donate to a named base, directly or through one
                bridging water? A general-base mechanism needs both at once.
  poses         which docked poses are there, and which one are we actually on?
  graft         put another pose into the equilibrated system, so the poses can
                be compared in one system rather than in systems that differ.
  check         does this structure meet the criteria for a Michaelis complex
                at all? Cheapest of the lot, and the one to run first: it reads
                `michaelis_criteria` from the site file and exits non-zero if a
                criterion required at this stage is violated, so a start
                structure is accepted by a command rather than by eye.
"""
from __future__ import annotations

import collections
import math
import sys

import yaml


OPTS = {}   # --key=value from the command line, read by the modes


# --- structure handling ------------------------------------------------------

class Structure(list):
    """Atoms, carrying the unit cell so distances can cross the box edge."""
    cell = None


def load_pdb(path):
    """[(chain, resid, atom, resname, (x, y, z)), ...] in file order, with .cell."""
    out = Structure()
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            out.append((line[21], line[22:26].strip(), line[12:16].strip(),
                        line[17:20].strip(),
                        (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
        elif line.startswith("CRYST1"):
            lengths = [float(line[6:15]), float(line[15:24]), float(line[24:33])]
            angles = [float(line[33:40]), float(line[40:47]), float(line[47:54])]
            if all(abs(a - 90.0) < 1e-3 for a in angles):
                out.cell = lengths
            else:
                print("[WARNING] %s has a non-orthorhombic cell; distances are "
                      "measured without the minimum-image convention" % path)
    return out


# AutoDock types are not elements: aromatic carbon is "A", and the acceptor
# forms carry a trailing letter.
AUTODOCK_ELEMENT = {"A": "C", "NA": "N", "OA": "O", "SA": "S", "HD": "H"}


def load_pdbqt(path):
    """([pose heavy-atom coordinates, ...], [Vina score, ...], [element, ...]).

    Heavy atoms only. A PDBQT has no nonpolar hydrogens, so anything compared
    against it has to drop its own or the two are not the same molecule.
    """
    poses, scores, elements, current = [], [], [], []
    for line in open(path):
        if line.startswith("MODEL"):
            current = []
        elif line.startswith("REMARK VINA RESULT"):
            scores.append(float(line.split()[3]))
        elif line.startswith(("ATOM", "HETATM")):
            raw = line[77:79].strip() or line[12:16].strip()[:1]
            el = AUTODOCK_ELEMENT.get(raw, raw[0])
            if el != "H":
                current.append((float(line[30:38]), float(line[38:46]),
                                float(line[46:54])))
                if not poses:
                    elements.append(el)
        elif line.startswith("ENDMDL"):
            poses.append(current)
    if not poses and current:
        poses = [current]
    return poses, scores, elements


def index_of(atoms, selector):
    """'chain:resid:atom' -> index, or None."""
    chain, resid, name = selector.split(":")
    for i, a in enumerate(atoms):
        if a[0] == chain and a[1] == resid and a[2] == name:
            return i
    return None


def label(atoms, i):
    a = atoms[i]
    return "%s:%s:%s:%s" % (a[0], a[1], a[3], a[2])


def delta(atoms, i, j):
    """j - i, as the shortest vector across the periodic box if there is one.

    An MD frame is written with every atom wrapped into the cell, so a ligand
    sitting near a face comes out on the far side of the protein: measured
    naively, a 5.0 A contact reads as 76.7 A in an 80 A box. Only orthorhombic
    cells are handled, which is what the solvated boxes here are.
    """
    v = [atoms[j][4][t] - atoms[i][4][t] for t in range(3)]
    cell = getattr(atoms, "cell", None)
    if cell:
        v = [x - c * round(x / c) if c else x for x, c in zip(v, cell)]
    return v


def dist(atoms, i, j):
    return math.hypot(*delta(atoms, i, j))


def angle(atoms, i, j, k):
    v1, v2 = delta(atoms, j, i), delta(atoms, j, k)
    cos = sum(p * q for p, q in zip(v1, v2)) / (math.hypot(*v1) * math.hypot(*v2))
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))


def hydrogens_on(atoms, i, cutoff=1.3):
    return [j for j, a in enumerate(atoms) if a[2][0] == "H" and dist(atoms, i, j) < cutoff]


def waters_near(atoms, centre, cutoff):
    """Water oxygens within `cutoff` of atom index `centre`, nearest first."""
    out = [(dist(atoms, centre, i), i) for i, a in enumerate(atoms)
           if a[3] in ("HOH", "WAT", "OH") and a[2][0] == "O" and i != centre]
    return sorted(d for d in out if d[0] <= cutoff)


# --- modes -------------------------------------------------------------------

def mode_sphere(site, paths):
    bd = site["burgi_dunitz"]
    for path in paths:
        at = load_pdb(path)
        zn = index_of(at, site["metal"])
        print("==", path)
        if zn is None:
            print("   metal %s not found" % site["metal"])
            continue
        sphere = sorted((dist(at, zn, i), i) for i, a in enumerate(at)
                        if a[2][0] in "ONS" and i != zn and dist(at, zn, i) < 2.8)
        print("   %d ligands within 2.8 A:" % len(sphere))
        for d, i in sphere:
            print("      %5.2f  %s" % (d, label(at, i)))
        for name, sel in (site.get("monitor") or {}).items():
            i, j = index_of(at, sel[0]), index_of(at, sel[1])
            if i is not None and j is not None:
                print("   %-14s %5.2f" % (name, dist(at, i, j)))
        c4 = index_of(at, site["electrophile"])
        o2 = index_of(at, site["carbonyl_o"])
        for name, sel in (site.get("nucleophiles") or {}).items():
            i = index_of(at, sel)
            if i is not None and c4 is not None:
                print("   %-28s O-C %5.2f   angle %3.0f"
                      % (name, dist(at, i, c4), angle(at, i, c4, o2)))
        for name, sel in list((site.get("nucleophiles") or {}).items()) + \
                         [("Lys226 NZ", site["bases"].get("Lys226 NZ", ""))]:
            i = index_of(at, sel) if sel else None
            if i is None:
                continue
            hs = hydrogens_on(at, i)
            print("   %-28s %d H  %s" % (name, len(hs),
                  ", ".join("%s %.2f" % (at[j][2], dist(at, i, j)) for j in hs)))


def mode_nucleophiles(site, paths):
    bd = site["burgi_dunitz"]
    lo, hi = bd["angle"]
    for path in paths:
        at = load_pdb(path)
        zn = index_of(at, site["metal"])
        c4 = index_of(at, site["electrophile"])
        o2 = index_of(at, site["carbonyl_o"])
        print("==", path)
        print("   %-20s %6s %7s %8s   %s" % ("water", "O-C", "angle", "M-O", "verdict"))
        cands = waters_near(at, c4, 6.0)
        for oc, i in cands:
            a = angle(at, i, c4, o2)
            mo = dist(at, zn, i) if zn is not None else float("nan")
            why = []
            if mo < bd["metal_bound"]:
                why.append("metal-bound")
            if not (lo <= a <= hi):
                why.append("angle off")
            if oc > bd["max_distance"]:
                why.append("far")
            print("   %-20s %6.2f %6.0f %8.2f   %s"
                  % ("%s:%s:%s" % (at[i][0], at[i][1], at[i][3]), oc, a, mo,
                     "VIABLE" if not why else ", ".join(why)))
        if not cands:
            print("   no water within 6 A of the electrophile")
        for oc, i in cands:
            hs = hydrogens_on(at, i)
            if not hs:
                continue
            print("   %s:%s protons:" % (at[i][0], at[i][1]))
            for h in hs:
                near = []
                for name, sel in (site.get("bases") or {}).items():
                    k = index_of(at, sel)
                    if k is not None:
                        near.append((dist(at, h, k), name))
                for k, b in enumerate(at):
                    if b[2][0] in "ON" and k != i and abs(k - i) > 2 \
                            and dist(at, h, k) < bd["hbond"]:
                        near.append((dist(at, h, k), label(at, k)))
                near.sort()
                shown = [x for x in near if x[0] < 4.0][:4]
                print("      %-5s %s" % (at[h][2],
                      ", ".join("%s %.2f" % (n, d) for d, n in shown) or "nothing within 4 A"))


def _traj_setup(site, traj_file, top_file):
    import mdtraj as md
    top = md.load_topology(top_file)

    def one(selector):
        chain, resid, name = selector.split(":")
        for a in top.atoms:
            if a.name == name and str(a.residue.resSeq) == resid:
                return a.index
        raise SystemExit("not found in topology: " + selector)

    print("loading trajectory...")
    traj = md.load(traj_file, top=top_file)
    print("%d frames" % traj.n_frames)
    waters = [r for r in top.residues if r.name in ("HOH", "WAT")]
    w_o = {r.index: [a.index for a in r.atoms if a.element.symbol == "O"][0] for r in waters}
    w_h = {r.index: [a.index for a in r.atoms if a.element.symbol == "H"] for r in waters}
    return md, top, traj, one, waters, w_o, w_h


def _attack_masks(md, traj, site, one, w_o, w_h):
    """{residue index: boolean mask of frames on the Burgi-Dunitz line}."""
    import numpy as np
    bd = site["burgi_dunitz"]
    lo, hi = bd["angle"]
    zn, c4, o2 = one(site["metal"]), one(site["electrophile"]), one(site["carbonyl_o"])
    first = traj[0]
    pre = md.compute_distances(first, [[w_o[r], c4] for r in w_o])[0] * 10
    near = [r for r, d in zip(w_o, pre) if d < 12.0]
    masks = {}
    for r in near:
        d = md.compute_distances(traj, [[w_o[r], c4], [w_o[r], zn]]) * 10
        a = np.degrees(md.compute_angles(traj, [[w_o[r], c4, o2]])[:, 0])
        m = (d[:, 0] < bd["max_distance"]) & (d[:, 1] > bd["metal_bound"]) & (a > lo) & (a < hi)
        if m.any():
            masks[r] = (m, d[:, 0], a)
    return masks, near


def mode_acceptors(site, traj_file, top_file):
    import numpy as np
    md, top, traj, one, waters, w_o, w_h = _traj_setup(site, traj_file, top_file)
    masks, _ = _attack_masks(md, traj, site, one, w_o, w_h)
    if not masks:
        raise SystemExit("no water reaches the attack geometry in any frame")
    res = {r.index: r for r in waters}
    best = max(masks.items(), key=lambda kv: kv[1][0].sum())
    ri, (mask, dd, aa) = best
    r = res[ri]
    print("attacking water: %s%s, on the line in %d of %d frames\n"
          % (r.name, r.resSeq, int(mask.sum()), traj.n_frames))
    sub = traj[mask]
    heavy_all = [a.index for a in top.atoms if a.element.symbol in ("O", "N")]
    local = md.compute_distances(sub[0], [[w_o[ri], k] for k in heavy_all])[0] * 10
    heavy = [k for k, d in zip(heavy_all, local)
             if d < 8.0 and top.atom(k).residue.index != ri]
    print("%d O/N atoms within 8 A of the attacking oxygen\n" % len(heavy))
    counts, closest = collections.Counter(), {}
    for h in w_h[ri]:
        d = md.compute_distances(sub, [[h, k] for k in heavy]) * 10
        near = d < site["burgi_dunitz"]["hbond"]
        for col in np.flatnonzero(near.any(axis=0)):
            a = top.atom(heavy[col])
            name = "%s%s %s" % (a.residue.name, a.residue.resSeq, a.name)
            counts[name] += int(near[:, col].sum())
            closest[name] = min(closest.get(name, 9.9), float(d[:, col].min()))
    if not counts:
        print("Nothing comes within hydrogen-bond distance of either proton while attacking.")
        return
    print("%-24s %8s %10s" % ("acceptor", "frames", "closest"))
    for name, c in counts.most_common(12):
        print("%-24s %8d %9.2f A" % (name, c, closest[name]))
    print("\n(out of %d attacking frames)" % int(mask.sum()))


def mode_joint(site, traj_file, top_file):
    import numpy as np
    md, top, traj, one, waters, w_o, w_h = _traj_setup(site, traj_file, top_file)
    bd = site["burgi_dunitz"]
    masks, near = _attack_masks(md, traj, site, one, w_o, w_h)
    if not masks:
        raise SystemExit("no water reaches the attack geometry in any frame")
    base_idx = {n: one(s) for n, s in (site.get("bases") or {}).items()}
    res = {r.index: r for r in waters}
    n = traj.n_frames

    def donates(ri):
        m = np.zeros(n, dtype=bool)
        for h in w_h[ri]:
            d = md.compute_distances(traj, [[h, k] for k in base_idx.values()]) * 10
            m |= (d < bd["hbond"]).any(axis=1)
        return m

    print("\n%-16s %8s %9s %8s   %s"
          % ("water", "attack", "H->base", "BOTH", "best O-C / angle"))
    both_any = 0
    for ri, (mask, dd, aa) in sorted(masks.items(), key=lambda kv: -kv[1][0].sum()):
        hb = donates(ri)
        both = mask & hb
        both_any += int(both.sum())
        k = int(np.argmin(np.where(mask, dd, 1e9)))
        print("%-16s %8d %9d %8d   %.2f A / %.0f deg (frame %d)"
              % ("%s%s" % (res[ri].name, res[ri].resSeq), int(mask.sum()),
                 int(hb.sum()), int(both.sum()), dd[k], aa[k], k))
    if both_any:
        print("\nSome frames satisfy both: use one of those as the reactant.")
        return
    print("\nNo frame has both halves: the attacking water and the water that donates")
    print("to a base are never the same molecule.")
    print("\nTesting the relay: attacker -> bridging water -> base")
    bridges = {r: donates(r) for r in near if donates(r).any()}
    print("%d nearby waters donate to a base in at least one frame" % len(bridges))
    best = (0, None, None)
    for ri, (mask, _, _) in masks.items():
        for bi, bhb in bridges.items():
            if bi == ri:
                continue
            reach = np.zeros(n, dtype=bool)
            for h in w_h[ri]:
                reach |= (md.compute_distances(traj, [[h, w_o[bi]]])[:, 0] * 10) < bd["hbond"]
            j = mask & reach & bhb
            if j.sum() > best[0]:
                best = (int(j.sum()), ri, bi)
    if best[0] == 0:
        print("No frame has attacker, bridge and base aligned either, in any of the")
        print("%d frames. The route is not pre-organised; whichever is taken, the" % n)
        print("model must include the reorganisation that sets it up.")
    else:
        cnt, ri, bi = best
        print("%s%s -> %s%s -> base aligned in %d frames"
              % (res[ri].name, res[ri].resSeq, res[bi].name, res[bi].resSeq, cnt))


def mode_shell(site, traj_file, top_file):
    """The metal's coordination sphere over a trajectory, ligand by ligand.

    `sphere` answers this for one structure, which is the wrong unit: a single
    frame cannot tell a ligand that is bound from one that happens to be close
    at that instant, and the coordination *number* is the quantity semiempirical
    methods are known to get wrong (D114: AM1/d was out by 7.55 kcal/mol on a
    4 -> 5 coordinate change). So every donor is reported with the fraction of
    frames it spends inside the cutoff, and the per-frame coordination number is
    reported as a distribution rather than an average alone - a site flipping
    between 4 and 5 and a site sitting at 4.5 are different things.

    Donors just outside the cutoff are listed too. A ligand at 2.7 A with a 2.6
    cutoff is a borderline call about this site, not an absence, and hiding it
    would make the cutoff look like a measurement.
    """
    import numpy as np
    md, top, traj, one, waters, w_o, w_h = _traj_setup(site, traj_file, top_file)

    cfg = site.get("shell") or {}
    cutoff = float(cfg.get("cutoff", site["burgi_dunitz"]["metal_bound"]))
    donors = tuple(cfg.get("donor_elements", ("O", "N", "S")))
    margin = float(cfg.get("margin", 1.0))     # how far past the cutoff to list
    metals = cfg.get("metals") or [site["metal"]]

    for selector in metals:
        m = one(selector)
        print("\n== %s  (%s)" % (selector, traj_file))
        print("%d frames, cutoff %.2f A, donors %s, +%.1f A margin shown"
              % (traj.n_frames, cutoff, "/".join(donors), margin))

        # Prefilter on the first frame so the distance matrix stays small; the
        # window is generous because a ligand can exchange during the run.
        cand = [a.index for a in top.atoms
                if a.index != m and a.element.symbol in donors]
        pre = md.compute_distances(traj[0], [[m, i] for i in cand],
                                   periodic=True)[0] * 10.0
        cand = [i for i, d in zip(cand, pre) if d < cutoff + margin + 4.0]
        if not cand:
            print("   no donor atom within reach of the metal at all")
            continue

        dist = md.compute_distances(traj, [[m, i] for i in cand],
                                    periodic=True) * 10.0
        occ = (dist < cutoff).mean(axis=0) * 100.0
        shown = [k for k in range(len(cand))
                 if dist[:, k].min() < cutoff + margin]
        shown.sort(key=lambda k: (-occ[k], dist[:, k].mean()))

        print("   %-24s  occ%%    min  mean   max" % "ligand")
        for k in shown:
            tag = "" if occ[k] else "   (outside the cutoff in every frame)"
            print("   %-24s %5.1f  %5.2f %5.2f  %5.2f%s"
                  % (label_traj(top, cand[k]), occ[k], dist[:, k].min(),
                     dist[:, k].mean(), dist[:, k].max(), tag))

        cn = (dist < cutoff).sum(axis=1)
        counts = np.bincount(cn)
        spread = "   ".join("%d: %d frames (%.0f%%)"
                            % (n, c, 100.0 * c / traj.n_frames)
                            for n, c in enumerate(counts) if c)
        print("   coordination number: mean %.2f   %s" % (cn.mean(), spread))


def label_traj(top, index):
    a = top.atom(index)
    return "%s%d:%s" % (a.residue.name, a.residue.resSeq, a.name)


def receptor_metal(site):
    """The metal selector as the DOCKING RECEPTOR numbers it.

    The receptor is a separate file prepared before solvation, so its chains and
    residue numbers need not match the simulated system's - CSR-1's zinc is
    B:2:ZN in `CSR_nontruncated_1metal.pdb` and B:1:ZN in `finalsystem.pdb`.
    Falls back to `metal`, so a site whose two files agree needs no new key.
    """
    return (site.get("docking") or {}).get("receptor_metal") or site["metal"]

def mode_poses(site, args):
    """Where do docked poses sit relative to the metal, and which one did we use?

    Answers the question a single carried-forward pose cannot: are the
    alternatives in the same place, and is the structure we have actually the
    top-scored one? Vina's scoring is good to a couple of kcal/mol, so poses
    within that of each other are not distinguishable by score alone and the
    geometry has to decide.
    """
    pdbqt, receptor = args[0], args[1]
    reference = args[2] if len(args) > 2 else None

    def centroid(points):
        n = len(points)
        return tuple(sum(p[i] for p in points) / n for i in range(3))

    def overlap(a, b):
        """Mean nearest-neighbour distance, symmetrised.

        PDBQT and PDB order and name their atoms differently and PDBQT merges
        nonpolar hydrogens away, so there is no atom-to-atom correspondence to
        take an RMSD over. This asks the question that actually matters anyway:
        does this pose occupy the same space as that structure.
        """
        ab = sum(min(math.dist(p, q) for q in b) for p in a) / len(a)
        ba = sum(min(math.dist(q, p) for p in a) for q in b) / len(b)
        return (ab + ba) / 2

    metal = index_of(load_pdb(receptor), receptor_metal(site))
    if metal is None:
        raise SystemExit("metal %s not found in %s"
                         % (receptor_metal(site), receptor))
    zn = load_pdb(receptor)[metal][4]

    poses, scores, _ = load_pdbqt(pdbqt)

    ref_pts = ref_c = None
    if reference:
        atoms = load_pdb(reference)
        ref_pts = [a[4] for a in atoms if a[2][0] != "H"]
        ref_c = centroid(ref_pts)
        print("reference %s: %d heavy atoms of %d, centroid to metal %.2f A"
              % (reference, len(ref_pts), len(atoms), math.dist(ref_c, zn)))
    print("metal %s at %.2f %.2f %.2f in %s"
          % (receptor_metal(site), zn[0], zn[1], zn[2], receptor))
    print()

    print("%-6s %8s %7s %10s %11s %9s"
          % ("pose", "score", "heavy", "to metal", "centroid", "overlap"))
    for i, pts in enumerate(poses, 1):
        c = centroid(pts)
        score = "%8.3f" % scores[i - 1] if i <= len(scores) else "       -"
        cen = "%10.2f A" % math.dist(c, ref_c) if ref_c else "          -"
        ov = "%8.2f A" % overlap(pts, ref_pts) if ref_pts else "         -"
        print("%-6d %s %7d %9.2f A %s %s"
              % (i, score, len(pts), math.dist(c, zn), cen, ov))

    if ref_pts and poses:
        best = min(range(len(poses)), key=lambda i: overlap(poses[i], ref_pts))
        d = overlap(poses[best], ref_pts)
        print()
        print("closest pose to the reference: %d, overlap %.2f A" % (best + 1, d))
        if d < 0.3:
            print("  the reference is that pose, essentially unchanged")
        elif d < 1.5:
            print("  the reference is that pose after being moved or minimised")
        else:
            print("  the reference matches no pose; check where it came from")


def bond_graph(points, cutoff=1.9):
    """{i: {neighbours}} by distance. Crude, but a covalent bond is short."""
    graph = {i: set() for i in range(len(points))}
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            if math.dist(points[i], points[j]) < cutoff:
                graph[i].add(j)
                graph[j].add(i)
    return graph


def map_pose_to_reference(pose, pose_el, ref_pts, ref_el):
    """Which reference heavy atom is which pose heavy atom.

    Nearest-neighbour assignment is not safe even when the two structures are
    close: meropenem has atoms nearer to a neighbour's position than to their
    own, and the assignment collides. So match the molecules as graphs -- same
    elements, same bonds -- which is exact and does not care how far apart the
    two conformations are. RDKit is used only as the matcher; no bond orders
    are perceived and nothing is sanitised, so a docked beta-lactam cannot be
    misread on the way through.
    """
    if len(pose) != len(ref_pts):
        raise SystemExit("pose has %d heavy atoms, reference has %d"
                         % (len(pose), len(ref_pts)))
    try:
        from rdkit import Chem
    except ImportError:
        raise SystemExit("graft needs rdkit for the atom mapping")

    def build(points, els):
        mol = Chem.RWMol()
        for e in els:
            mol.AddAtom(Chem.Atom(e))
        for i, nbrs in bond_graph(points).items():
            for j in nbrs:
                if j > i:
                    mol.AddBond(i, j, Chem.BondType.SINGLE)
        out = mol.GetMol()
        Chem.FastFindRings(out)
        return out

    query, target = build(ref_pts, ref_el), build(pose, pose_el)
    match = target.GetSubstructMatch(query)
    if not match:
        raise SystemExit("the pose and the reference are not the same molecule "
                         "as graphs; refusing to guess atom identities")
    mapping = [None] * len(pose)
    for ref_i, pose_i in enumerate(match):
        mapping[pose_i] = ref_i
    if any(m is None for m in mapping):
        raise SystemExit("graph match did not cover every atom")
    print("atom mapping verified: %d atoms matched as a labelled graph"
          % len(mapping))
    return mapping


def kabsch(mobile, target):
    """Rigid transform taking `mobile` onto `target`; returns (apply, rmsd)."""
    import numpy as np
    a, b = np.asarray(mobile, float), np.asarray(target, float)
    ca, cb = a.mean(0), b.mean(0)
    u, _, vt = np.linalg.svd((a - ca).T @ (b - cb))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    rmsd = float(np.sqrt((((a - ca) @ rot.T - (b - cb)) ** 2).sum(1).mean()))
    return (lambda pts: [tuple(v) for v in
                         (np.asarray(pts, float) - ca) @ rot.T + cb], rmsd)


def frame_transform(docked_in, system, metal_selector, radius=12.0):
    """Put docking-frame coordinates into the simulated system's frame.

    A solvated system is built in its own box, so it is translated and may be
    rotated relative to the structure the ligand was docked into. Aligning on
    the residues around the metal rather than on the whole protein keeps the
    fit honest about the part the pose has to land in.
    """
    a, b = load_pdb(docked_in), load_pdb(system)
    mi = index_of(b, metal_selector)
    if mi is None:
        raise SystemExit("metal %s not found in %s" % (metal_selector, system))
    key = lambda x: (x[0], x[1], x[2])
    pa = {key(x): x[4] for x in a}
    pb = {key(x): x[4] for x in b}

    def fit(backbone_only):
        near = [key(x) for x in b
                if x[3] not in ("HOH", "WAT", "SOL")
                and (not backbone_only or x[2] in ("N", "CA", "C", "O"))
                and math.dist(x[4], b[mi][4]) < radius]
        shared = [k for k in near if k in pa]
        if len(shared) < 20:
            return None, None, len(shared)
        apply_, rmsd = kabsch([pa[k] for k in shared], [pb[k] for k in shared])
        return apply_, rmsd, len(shared)

    # Side chains rearrange during MD; backbone does not. Fitting on everything
    # lets that rearrangement masquerade as a frame difference, so fit on the
    # backbone and report the all-atom number only as a measure of how much the
    # site moved.
    apply_, rmsd, n = fit(True)
    _, all_rmsd, all_n = fit(False)
    if apply_ is None:
        raise SystemExit("only %d shared backbone atoms near the metal" % n)
    print("docking frame aligned onto the system on %d backbone atoms within "
          "%.0f A of the metal, RMSD %.2f A" % (n, radius, rmsd))
    print("  for comparison, all %d atoms there fit to %.2f A: that is how far "
          "the site itself moved" % (all_n, all_rmsd))
    if rmsd > 1.0:
        print("  WARNING: the backbone itself has moved; a docked pose may not "
              "be meaningful in this frame")
    return apply_


def place_hydrogens(ref_atoms, heavy_idx, new_heavy):
    """Carry each hydrogen across on the local frame of its parent heavy atom.

    Rebuilding hydrogens with a force field would mean sanitising a docked
    structure through a perception step that can get a beta-lactam wrong. The
    local geometry around a hydrogen does not change between poses -- only the
    torsions do -- so transferring it in the parent's frame is both exact for
    the bits that matter and free of chemistry perception.
    """
    def frame(origin, a, b):
        e1 = [a[k] - origin[k] for k in range(3)]
        n1 = math.sqrt(sum(v * v for v in e1)) or 1.0
        e1 = [v / n1 for v in e1]
        e2 = [b[k] - origin[k] for k in range(3)]
        dot = sum(e1[k] * e2[k] for k in range(3))
        e2 = [e2[k] - dot * e1[k] for k in range(3)]
        n2 = math.sqrt(sum(v * v for v in e2)) or 1.0
        e2 = [v / n2 for v in e2]
        e3 = [e1[1] * e2[2] - e1[2] * e2[1],
              e1[2] * e2[0] - e1[0] * e2[2],
              e1[0] * e2[1] - e1[1] * e2[0]]
        return e1, e2, e3

    old_heavy = [ref_atoms[i][4] for i in heavy_idx]
    out = {}
    for h, atom in enumerate(ref_atoms):
        if atom[2][0] != "H":
            continue
        pos = atom[4]
        order = sorted(range(len(heavy_idx)),
                       key=lambda k: math.dist(pos, old_heavy[k]))
        p, a, b = order[0], order[1], order[2]
        old_f = frame(old_heavy[p], old_heavy[a], old_heavy[b])
        new_f = frame(new_heavy[p], new_heavy[a], new_heavy[b])
        d = [pos[k] - old_heavy[p][k] for k in range(3)]
        local = [sum(d[k] * axis[k] for k in range(3)) for axis in old_f]
        out[h] = tuple(new_heavy[p][k] + sum(local[c] * new_f[c][k]
                                             for c in range(3))
                       for k in range(3))
    return out


def mode_graft(site, args):
    """Put a different docked pose into the system we already equilibrated.

    Re-solvating each pose from scratch would give every pose a different water
    count, and then the total energies are not comparable at all. Swapping the
    ligand into one equilibrated system keeps the protein and the waters fixed,
    so what differs between the runs is the thing under test.

        site.py graft <pdbqt> <ligand.pdb> <receptor.pdb> <system.pdb> [poses]

    `ligand.pdb` and `receptor.pdb` are the docking frame -- the ligand gives
    the atom names and hydrogens, the receptor gives the alignment onto the
    system. With no pose list it only measures, which is the cheap triage: a
    pose that puts nothing near the electrophile is not worth optimising.
    """
    pdbqt, ref_lig, receptor, system = args[0], args[1], args[2], args[3]
    wanted = ([int(p) for p in args[4].split(",")] if len(args) > 4 else [])

    poses, scores, pose_el = load_pdbqt(pdbqt)
    ref_atoms = load_pdb(ref_lig)
    heavy_idx = [i for i, a in enumerate(ref_atoms) if a[2][0] != "H"]
    ref_pts = [ref_atoms[i][4] for i in heavy_idx]
    ref_el = [ref_atoms[i][2][0] for i in heavy_idx]

    best = min(range(len(poses)),
               key=lambda i: sum(min(math.dist(p, q) for q in poses[i])
                                 for p in ref_pts))
    print("reference ligand matches pose %d; deriving atom names from it"
          % (best + 1))
    mapping = map_pose_to_reference(poses[best], pose_el, ref_pts, ref_el)

    # Everything from here on is in the system's frame, including the ligand
    # the names came from, so distances mean what they say.
    into_system = frame_transform(receptor, system, site["metal"])
    poses = [into_system(p) for p in poses]
    ref_atoms = [a[:4] + (c,) for a, c in
                 zip(ref_atoms, into_system([a[4] for a in ref_atoms]))]

    # The reference ligand sits in the system, so find it there by name.
    sys_atoms = load_pdb(system)
    lig_chain, lig_resid = ref_atoms[0][0], ref_atoms[0][1]
    lig_rows = [i for i, a in enumerate(sys_atoms)
                if a[0] == lig_chain and a[1] == lig_resid]
    if len(lig_rows) != len(ref_atoms):
        raise SystemExit("system has %d atoms for %s:%s, reference ligand has %d"
                         % (len(lig_rows), lig_chain, lig_resid, len(ref_atoms)))

    carb = index_of(ref_atoms, site["carbonyl_o"])
    if carb is None:
        raise SystemExit("carbonyl %s not found in %s"
                         % (site["carbonyl_o"], ref_lig))
    metal_i = index_of(sys_atoms, site["metal"])
    zn = sys_atoms[metal_i][4]
    bd = site["burgi_dunitz"]

    # Is the ligand in the system the same conformer as the pose, just moved?
    # Optimal superposition removes the placement and leaves the conformation,
    # which separates "this pose was repositioned" from "this is not that pose".
    sys_lig = [a[4] for a in sys_atoms
               if a[0] == lig_chain and a[1] == lig_resid and a[2][0] != "H"]

    print()
    print("%-6s %8s %10s %9s %12s %10s"
          % ("pose", "score", "Zn-carbonyl", "conformer", "waters <2.2A",
             "written"))
    for n, pose in enumerate(poses, 1):
        new_heavy = [None] * len(heavy_idx)
        for pi, ri in enumerate(mapping):
            new_heavy[ri] = pose[pi]

        pos_of = {heavy_idx[k]: new_heavy[k] for k in range(len(heavy_idx))}
        pos_of.update(place_hydrogens(ref_atoms, heavy_idx, new_heavy))

        znco = math.dist(pos_of[carb], zn)
        clashing = waters_within(sys_atoms, pos_of, lig_rows, 2.2)
        written = "-"
        if n in wanted:
            written = write_grafted(system, sys_atoms, lig_rows, ref_atoms,
                                    pos_of, clashing, n)
        score = "%8.3f" % scores[n - 1] if n <= len(scores) else "       -"
        conf = "        -"
        if len(sys_lig) == len(pose):
            ordered = [None] * len(pose)
            for pi, ri in enumerate(mapping):
                ordered[pi] = sys_lig[ri]
            _, crmsd = kabsch(ordered, pose)
            conf = "%7.2f A" % crmsd
        print("%-6d %s %9.2f A %9s %12d %10s"
              % (n, score, znco, conf, len(clashing), written))

    # The ligand already in the system is the control: if the docked pose it
    # came from sits somewhere else entirely, the pipeline moved it there.
    in_system = next((a[4] for a in sys_atoms
                      if a[0] == lig_chain and a[1] == lig_resid
                      and a[2] == site["carbonyl_o"].rsplit(":", 1)[1]), None)
    print()
    print("Zn-carbonyl under %.1f A is inner-sphere." % bd["metal_bound"])
    print("  the docked pose, aligned into this frame: %.2f A"
          % math.dist(ref_atoms[carb][4], zn))
    if in_system:
        print("  the ligand actually in the system:         %.2f A"
              % math.dist(in_system, zn))


def waters_within(sys_atoms, pos_of, lig_rows, cutoff):
    """Water residues with any atom within `cutoff` of a grafted ligand atom."""
    pts = list(pos_of.values())
    hit = set()
    for i, a in enumerate(sys_atoms):
        if i in set(lig_rows) or a[3] not in ("HOH", "WAT", "SOL"):
            continue
        if any(math.dist(a[4], p) < cutoff for p in pts):
            hit.add((a[0], a[1]))
    return hit


def write_grafted(system, sys_atoms, lig_rows, ref_atoms, pos_of, clashing, n):
    """The system PDB with the ligand replaced and clashing waters removed."""
    out = "%s_pose%d.pdb" % (system.rsplit(".", 1)[0], n)
    row_for = {row: k for k, row in enumerate(lig_rows)}
    with open(out, "w") as fh:
        i = -1
        for line in open(system):
            if not line.startswith(("ATOM", "HETATM")):
                fh.write(line)
                continue
            i += 1
            a = sys_atoms[i]
            if (a[0], a[1]) in clashing:
                continue
            if i in row_for:
                x, y, z = pos_of[row_for[i]]
                line = "%s%8.3f%8.3f%8.3f%s" % (line[:30], x, y, z, line[54:])
            fh.write(line)
    return out.rsplit("/", 1)[-1]


# --- entry point -------------------------------------------------------------

def coordination(site, at, metal_i):
    """Donors on the metal, using the same keys `shell` uses on a trajectory."""
    sh = site.get("shell") or {}
    cutoff = sh.get("cutoff", site["burgi_dunitz"]["metal_bound"])
    elements = tuple(sh.get("donor_elements") or ["O", "N", "S"])
    return sorted((dist(at, metal_i, i), i) for i, a in enumerate(at)
                  if i != metal_i and a[2][0] in elements
                  and dist(at, metal_i, i) <= cutoff)


def _criterion(site, at, c, stage):
    """(verdict, measured, detail) for one criterion against one structure.

    Verdicts: PASS, FAIL, WARN. A criterion is only ever FAIL when `stage` is in
    its `required_at`; otherwise a violation is a WARN, because a criterion that
    describes the reactant state is not a defect in a structure that has yet to
    be equilibrated into one.
    """
    required = stage in (c.get("required_at") or ["start", "reactant"])
    bad = "FAIL" if required else "WARN"
    lo, hi = c.get("range", (float("-inf"), float("inf")))

    if c.get("closest_to"):
        ref = index_of(at, c["closest_to"])
        if ref is None:
            return bad, "-", "%s not in the structure" % c["closest_to"]
        seen = [(dist(at, ref, i), s) for s in c["among"]
                for i in [index_of(at, s)] if i is not None]
        if not seen:
            return bad, "-", "none of %s present" % ", ".join(c["among"])
        seen.sort()
        d, nearest = seen[0]
        detail = ", ".join("%s %.2f" % (s, x) for x, s in seen)
        if nearest != c["must_be"]:
            return bad, "%s %.2f" % (nearest, d),                    "must be %s; %s" % (c["must_be"], detail)
        return "PASS", "%s %.2f" % (nearest, d), detail

    if c.get("coordination"):
        m = index_of(at, c["coordination"])
        if m is None:
            return bad, "-", "%s not in the structure" % c["coordination"]
        sphere = coordination(site, at, m)
        detail = ", ".join("%s %.2f" % (label(at, i), d) for d, i in sphere)
        n = len(sphere)
        return ("PASS" if lo <= n <= hi else bad), "%d" % n, detail

    idx = [index_of(at, s) for s in c["atoms"]]
    if any(i is None for i in idx):
        absent = [s for s, i in zip(c["atoms"], idx) if i is None]
        return bad, "-", "not in the structure: %s" % ", ".join(absent)
    if len(idx) == 3:
        v = angle(at, idx[0], idx[1], idx[2])
        shown = "%.0f" % v
    else:
        v = dist(at, idx[0], idx[1])
        shown = "%.2f" % v
    return ("PASS" if lo <= v <= hi else bad), shown, ""


def mode_check(site, args):
    """Accept or reject a candidate structure against `michaelis_criteria`.

    This mode exists because of D123: the pose everything ran on from October
    had the substrate bound the other way round, with a substituent oxygen on
    the metal and the carbonyl 5.43 A away, and the number sat in the decision
    log for two weeks without being read as a defect. Nothing here is new
    knowledge - every criterion was already measured and recorded. What was
    missing was a command that says yes or no, so a start structure is not
    accepted by eye. Exits non-zero if anything required at this stage fails.
    """
    stage = OPTS.get("at", "start")
    criteria = site.get("michaelis_criteria") or []
    if not criteria:
        raise SystemExit("no `michaelis_criteria` in the site file: nothing to check")
    worst = 0
    for path in args:
        at = load_pdb(path)
        print("==", path, "   (criteria required at: %s)" % stage)
        for c in criteria:
            verdict, measured, detail = _criterion(site, at, c, stage)
            want = c.get("range")
            window = ("%g-%g" % tuple(want)) if want else c.get("must_be", "")
            print("   %-4s %-26s %10s   want %-10s %s"
                  % (verdict, c["name"], measured, window, c.get("source", "")))
            if detail:
                print("          %s" % detail)
            if verdict == "FAIL":
                worst = 1
        print("   %s" % ("REJECTED" if worst else "accepted at stage '%s'" % stage))
    raise SystemExit(worst)


MODES = {"sphere": mode_sphere, "nucleophiles": mode_nucleophiles,
         "acceptors": mode_acceptors, "joint": mode_joint, "poses": mode_poses,
         "graft": mode_graft, "shell": mode_shell,
         "check": mode_check}

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    OPTS.update(o.lstrip("-").split("=", 1) for o in sys.argv[1:]
                if o.startswith("--") and "=" in o)
    site_file = next((a.split("=", 1)[1] for a in sys.argv[1:]
                      if a.startswith("--site=")), "site.yaml")
    if len(args) < 2 or args[0] not in MODES:
        raise SystemExit(__doc__)
    site = yaml.safe_load(open(site_file))["site"]
    mode, rest = args[0], args[1:]
    if mode in ("sphere", "nucleophiles", "poses", "graft", "check"):
        MODES[mode](site, rest)
    else:
        if len(rest) != 2:
            raise SystemExit("%s needs <trajectory> <topology>" % mode)
        MODES[mode](site, rest[0], rest[1])
