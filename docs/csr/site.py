"""Active-site geometry analysis, driven by a site description.

    python site.py sphere       <pdb> [pdb ...]     metal sphere and protonation
    python site.py nucleophiles <pdb> [pdb ...]     attack geometry + proton acceptors
    python site.py acceptors    <traj> <topology>   what CAN take the proton, counted
    python site.py joint        <traj> <topology>   attack AND donate, same molecule
    python site.py poses        <pdbqt> <receptor> [ligand.pdb]
                                                    docked poses vs the metal site
    python site.py graft        <pdbqt> <ligand.pdb> <system.pdb> [1,3,4]
                                                    a pose into the solvated system

Every atom, every candidate base and every criterion comes from `site.yaml`
(override with --site), so this script carries no knowledge of a particular
enzyme and a new system needs a new YAML file rather than new code.

The modes answer different questions, in rough order of cost:

  sphere        what is coordinated, and where are the hydrogens?
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
"""
from __future__ import annotations

import collections
import math
import sys

import yaml


# --- structure handling ------------------------------------------------------

def load_pdb(path):
    """[(chain, resid, atom, resname, (x, y, z)), ...] in file order."""
    out = []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            out.append((line[21], line[22:26].strip(), line[12:16].strip(),
                        line[17:20].strip(),
                        (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
    return out


def load_pdbqt(path):
    """([pose heavy-atom coordinates, ...], [Vina score, ...]).

    Heavy atoms only. A PDBQT has no nonpolar hydrogens, so anything compared
    against it has to drop its own or the two are not the same molecule.
    """
    poses, scores, current = [], [], []
    for line in open(path):
        if line.startswith("MODEL"):
            current = []
        elif line.startswith("REMARK VINA RESULT"):
            scores.append(float(line.split()[3]))
        elif line.startswith(("ATOM", "HETATM")):
            if (line[77:79].strip() or line[12:16].strip())[0] != "H":
                current.append((float(line[30:38]), float(line[38:46]),
                                float(line[46:54])))
        elif line.startswith("ENDMDL"):
            poses.append(current)
    if not poses and current:
        poses = [current]
    return poses, scores


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


def dist(atoms, i, j):
    return math.dist(atoms[i][4], atoms[j][4])


def angle(atoms, i, j, k):
    v1 = [atoms[i][4][t] - atoms[j][4][t] for t in range(3)]
    v2 = [atoms[k][4][t] - atoms[j][4][t] for t in range(3)]
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

    metal = index_of(load_pdb(receptor), site["metal"])
    if metal is None:
        raise SystemExit("metal %s not found in %s" % (site["metal"], receptor))
    zn = load_pdb(receptor)[metal][4]

    poses, scores = load_pdbqt(pdbqt)

    ref_pts = ref_c = None
    if reference:
        atoms = load_pdb(reference)
        ref_pts = [a[4] for a in atoms if a[2][0] != "H"]
        ref_c = centroid(ref_pts)
        print("reference %s: %d heavy atoms of %d, centroid to metal %.2f A"
              % (reference, len(ref_pts), len(atoms), math.dist(ref_c, zn)))
    print("metal %s at %.2f %.2f %.2f in %s"
          % (site["metal"], zn[0], zn[1], zn[2], receptor))
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


def map_pose_to_reference(pose, ref_pts):
    """Which reference heavy atom is which pose heavy atom.

    The pose this is derived from must be the one the reference came from, so
    the two differ by much less than a bond length and nearest-neighbour
    assignment is unambiguous. That is an assumption, so it is checked: the
    mapping has to be a bijection AND has to carry the pose's bond graph onto
    the reference's. If it does, the atom names are transferable to every other
    pose in the file, which all share the PDBQT's atom order.
    """
    if len(pose) != len(ref_pts):
        raise SystemExit("pose has %d heavy atoms, reference has %d"
                         % (len(pose), len(ref_pts)))
    mapping = [min(range(len(ref_pts)), key=lambda j: math.dist(p, ref_pts[j]))
               for p in pose]
    if len(set(mapping)) != len(mapping):
        raise SystemExit("mapping is not one-to-one; the reference is not this "
                         "pose, so atom names cannot be transferred")
    gp, gr = bond_graph(pose), bond_graph(ref_pts)
    for i, nbrs in gp.items():
        if {mapping[n] for n in nbrs} != gr[mapping[i]]:
            raise SystemExit("mapping does not preserve bonding at pose atom %d; "
                             "refusing to guess atom identities" % (i + 1))
    worst = max(math.dist(pose[i], ref_pts[mapping[i]]) for i in range(len(pose)))
    print("atom mapping verified: bijective, bond graph preserved, "
          "worst atom %.2f A" % worst)
    return mapping


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

        site.py graft <pdbqt> <reference_ligand.pdb> <system.pdb> [poses]

    With no pose list it only measures, which is the cheap triage: a pose that
    puts nothing near the electrophile is not worth optimising.
    """
    pdbqt, ref_lig, system = args[0], args[1], args[2]
    wanted = ([int(p) for p in args[3].split(",")] if len(args) > 3 else [])

    poses, scores = load_pdbqt(pdbqt)
    ref_atoms = load_pdb(ref_lig)
    heavy_idx = [i for i, a in enumerate(ref_atoms) if a[2][0] != "H"]
    ref_pts = [ref_atoms[i][4] for i in heavy_idx]

    best = min(range(len(poses)),
               key=lambda i: sum(min(math.dist(p, q) for q in poses[i])
                                 for p in ref_pts))
    print("reference ligand matches pose %d; deriving atom names from it"
          % (best + 1))
    mapping = map_pose_to_reference(poses[best], ref_pts)

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

    print()
    print("%-6s %8s %10s %12s %10s" % ("pose", "score", "Zn-carbonyl",
                                       "waters <%.1fA" % 2.2, "written"))
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
        print("%-6d %s %9.2f A %12d %10s" % (n, score, znco, len(clashing), written))

    print()
    print("Zn-carbonyl under %.1f A is inner-sphere; the reference pose sits at "
          "%.2f A" % (bd["metal_bound"], math.dist(ref_atoms[carb][4], zn)))


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

MODES = {"sphere": mode_sphere, "nucleophiles": mode_nucleophiles,
         "acceptors": mode_acceptors, "joint": mode_joint, "poses": mode_poses,
         "graft": mode_graft}

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    site_file = next((a.split("=", 1)[1] for a in sys.argv[1:]
                      if a.startswith("--site=")), "site.yaml")
    if len(args) < 2 or args[0] not in MODES:
        raise SystemExit(__doc__)
    site = yaml.safe_load(open(site_file))["site"]
    mode, rest = args[0], args[1:]
    if mode in ("sphere", "nucleophiles", "poses", "graft"):
        MODES[mode](site, rest)
    else:
        if len(rest) != 2:
            raise SystemExit("%s needs <trajectory> <topology>" % mode)
        MODES[mode](site, rest[0], rest[1])
