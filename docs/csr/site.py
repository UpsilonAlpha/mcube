"""Active-site geometry analysis, driven by a site description.

    python site.py sphere       <pdb> [pdb ...]     metal sphere and protonation
    python site.py nucleophiles <pdb> [pdb ...]     attack geometry + proton acceptors
    python site.py shell        <traj> <topology>   metal coordination over a trajectory
    python site.py contacts     <traj> <topology>   named pairs over a trajectory
                                [--profile=<substring>] [--bin=0.05]
    python site.py acceptors    <traj> <topology>   what CAN take the proton, counted
    python site.py joint        <traj> <topology>   attack AND donate, same molecule
    python site.py poses        <pdbqt> <receptor> [ligand.pdb]
                                                    docked poses vs the metal site
    python site.py graft        <pdbqt> <ligand.pdb> <receptor.pdb> <system.pdb>
                                                [1,3,4]  a pose into the system
                                       [--ref-frame=system]
    python site.py check        <pdb> [pdb ...] [--at=start|reactant]
                                       [--hypothesis=<name>]
                                                    accept or reject a structure
    python site.py protonate    <pdb> <chain:resid> --to=ASH [--out=<path>]
                                                    build the protonated variant

Every atom, every candidate base and every criterion comes from `site.yaml`
(override with --site), so this script carries no knowledge of a particular
enzyme and a new system needs a new YAML file rather than new code.

The modes answer different questions, in rough order of cost:

  sphere        what is coordinated, and where are the hydrogens?
  shell         over a trajectory, what is coordinated and how many of them, with
                every ligand's occupancy - the question `sphere` can only
                answer for the one frame it is given.
  contacts      over a trajectory, how often each pair named in `monitor` is
                within hydrogen-bonding reach, with the distance distribution.
                `shell` answers what is on the metal; this answers everything
                else - which residue holds the substrate's carboxylate, whether
                a glutamine reaches the carbonyl - and it answers it for the
                pairs the site description already names, so a contact is
                reported because the YAML names it and not because someone
                wrote a one-off script for this protein.
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
                one named set from `michaelis_hypotheses` and exits non-zero if
                a criterion required at this stage is violated, so a start
                structure is accepted by a command rather than by eye.
                `--hypothesis=` is required when more than one set is defined,
                because a criterion is only a criterion within a mechanism
                (D127). Its coordination criterion returns a COUNT, and four
                donors inside a cutoff can be the wrong four - where the
                identity of the sphere is the claim, read it off `shell`'s
                occupancies instead (D129, D130).
  protonate     not a measurement: the one mode that writes a structure. It
                renames a residue to its protonated Amber variant and builds
                the proton, so the alternative protonation state can be run at
                all. See `mode_protonate`.
"""
from __future__ import annotations

import collections
import math
import os
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
        """Atom index for "chain:resid:name", by residue number and atom name.

        The chain letter is only used to break a tie, because a trajectory
        topology does not have to carry the PDB's chain IDs - but a tie is
        exactly what used to go wrong silently: `D:1:O` matches the crystal
        water's oxygen AND residue 1's backbone O, and whichever came first in
        the topology was returned. An ambiguous selector now raises instead of
        guessing, and reports the candidates, because a wrong atom here is a
        measurement that looks right (D148).
        """
        chain, resid, name = selector.split(":")
        hits = [a for a in top.atoms
                if a.name == name and str(a.residue.resSeq) == resid]
        if not hits:
            raise SystemExit("not found in topology: " + selector)
        if len(hits) > 1 and chain:
            want = ord(chain[0].upper()) - ord("A")
            narrowed = [a for a in hits if a.residue.chain.index == want]
            if len(narrowed) == 1:
                return narrowed[0].index
        if len(hits) > 1:
            where = ", ".join("%s%d chain %d" % (a.residue.name,
                                                 a.residue.resSeq,
                                                 a.residue.chain.index)
                              for a in hits[:6])
            raise SystemExit("ambiguous selector %s: %d atoms match (%s)"
                             % (selector, len(hits), where))
        return hits[0].index

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

def mode_contacts(site, traj_file, top_file):
    """Every pair in `monitor`, over a trajectory: occupancy and distribution.

    `sphere` prints these distances for one frame, which cannot tell a contact
    that is held from one that happens to be there - the same objection
    `shell` exists for, applied to the pairs that are not on the metal. A
    residue "near the carboxylate" in one optimised structure and one at 100%
    occupancy over 5,000 frames are different claims, and only the second is
    worth writing down.

    Occupancy is the fraction of frames inside `contacts.cutoff` (3.5 A by
    default, the usual polar-contact convention for heavy-atom pairs). Pairs
    whose atoms are not in this topology are listed as missing rather than
    skipped silently, because a pair that is absent from the system is a fact
    about the model, not a typo to swallow.
    """
    import numpy as np
    md, top, traj, one, waters, w_o, w_h = _traj_setup(site, traj_file, top_file)

    cfg = site.get("contacts") or {}
    cutoff = float(cfg.get("cutoff", 3.5))
    pairs = site.get("monitor") or {}

    named, idx, missing = [], [], []
    for name, sel in pairs.items():
        try:
            i, j = one(sel[0]), one(sel[1])
        except SystemExit as e:
            missing.append((name, str(e)))
            continue
        named.append(name)
        idx.append([i, j])

    print()
    print("== %s  (%d frames, cutoff %.2f A)" % (traj_file, traj.n_frames, cutoff))
    if not idx:
        print("   no pair in `monitor` resolves in this topology")
        return
    dist = md.compute_distances(traj, idx, periodic=True) * 10.0
    occ = (dist < cutoff).mean(axis=0) * 100.0
    order = sorted(range(len(named)), key=lambda k: (-occ[k], dist[:, k].mean()))

    print("   %-26s  occ%%    min  mean   max" % "pair")
    for k in order:
        print("   %-26s %5.1f  %5.2f %5.2f  %5.2f"
              % (named[k], occ[k], dist[:, k].min(),
                 dist[:, k].mean(), dist[:, k].max()))
    for name, why in missing:
        print("   %-26s     -  %s" % (name, why))

    want = OPTS.get("profile")
    if want:
        _profile(np, named, dist, want, float(OPTS.get("bin", 0.05)),
                 float(cfg.get("temperature", 300.0)), traj.n_frames)


def _profile(np, named, dist, want, width, temperature, n_frames):
    """-RT ln P(r) for every pair whose name contains `want`.

    Unbiased sampling already carries the free energy along any distance it
    explored: the histogram of a 100 ns trajectory is a potential of mean force
    over the range it visited, for nothing but the cost of reading the frames.
    It is only as good as the sampling, so the count per bin is printed beside
    the energy - a bin with a handful of frames is noise wearing a number, and
    the top of a profile like this is always the thin end.

    This is not a substitute for umbrella sampling, which is what reaches the
    range the molecule does not visit on its own. It is the part of the answer
    that is already paid for.
    """
    rt = 0.0019872041 * temperature          # kcal/mol
    for k, name in enumerate(named):
        if want.lower() not in name.lower():
            continue
        d = dist[:, k]
        lo = np.floor(d.min() / width) * width
        hi = np.ceil(d.max() / width) * width
        edges = np.arange(lo, hi + width, width)
        counts, edges = np.histogram(d, bins=edges)
        centres = 0.5 * (edges[:-1] + edges[1:])
        keep = counts > 0
        # Jacobian: a distance shell grows as r^2, so the radial density is
        # divided by it before the log, or every profile tilts outward.
        p = counts[keep] / (centres[keep] ** 2)
        g = -rt * np.log(p / p.max())
        print()
        print("   -- %s: -RT ln P at %.0f K, %.2f A bins, %d frames"
              % (name, temperature, width, n_frames))
        print("      %6s %8s %10s" % ("r", "frames", "dG kcal"))
        for c, n, e in zip(centres[keep], counts[keep], g):
            print("      %6.2f %8d %10.2f" % (c, n, e))
        best = centres[keep][int(np.argmin(g))]
        print("      minimum at %.2f A; %d of %d frames in the range"
              % (best, counts.sum(), n_frames))


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

    `receptor.pdb` is the docking receptor and gives the alignment onto the
    system. `ligand.pdb` gives the atom names and the hydrogens, and has to be
    the same ligand the system holds, atom for atom. With no pose list it only
    measures, which is the cheap triage: a pose that puts nothing near the
    electrophile is not worth optimising.

    `--ref-frame=system` says the reference ligand is already in the system's
    frame rather than the docking frame. That is the case when the reference is
    extracted from the system itself, which is the only reference available
    when the docked run used an earlier variant of the ligand: CSR-1's docking
    predates meropenem's O4x, so the docking-frame copy is one atom short of
    the system's and cannot meet the atom-for-atom requirement from that side.
    Default `docking`, so existing calls are unchanged.
    """
    pdbqt, ref_lig, receptor, system = args[0], args[1], args[2], args[3]
    wanted = ([int(p) for p in args[4].split(",")] if len(args) > 4 else [])
    ref_in_system = OPTS.get("ref-frame", "docking") == "system"

    poses, scores, pose_el = load_pdbqt(pdbqt)
    ref_atoms = load_pdb(ref_lig)
    heavy_idx = [i for i, a in enumerate(ref_atoms) if a[2][0] != "H"]
    ref_pts = [ref_atoms[i][4] for i in heavy_idx]
    ref_el = [ref_atoms[i][2][0] for i in heavy_idx]

    # Everything from here on is in the system's frame, including the ligand
    # the names came from, so distances mean what they say. The pose-to-name
    # match has to run in whichever frame the reference is in, because it
    # compares geometry: a system-frame reference against docking-frame poses
    # compares two different coordinate systems and picks a pose at random.
    into_system = frame_transform(receptor, system, site["metal"])
    if ref_in_system:
        poses = [into_system(p) for p in poses]

    best = min(range(len(poses)),
               key=lambda i: sum(min(math.dist(p, q) for q in poses[i])
                                 for p in ref_pts))
    print("reference ligand matches pose %d; deriving atom names from it"
          % (best + 1))
    mapping = map_pose_to_reference(poses[best], pose_el, ref_pts, ref_el)

    if not ref_in_system:
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


# --- protonate ---------------------------------------------------------------
#
# The protonated Amber variants, as ff14SB actually defines them. Each entry was
# read off the force field rather than from a paper:
#
#   ForceField("amber14/protein.ff14SB.xml")._templates["ASH"]
#
#     ASH  N H CA HA CB HB2 HB3 CG OD1 OD2 HD2 C O      = ASP + HD2 on OD2
#     GLH  N H CA HA CB HB2 HB3 CG HG2 HG3 CD OE1 OE2 HE2 C O
#                                                       = GLU + HE2 on OE2
#     LYN  N H CA HA CB HB2 HB3 CG HG2 HG3 CD HD2 HD3 CE HE2 HE3 NZ HZ2 HZ3 C O
#                                                       = LYS - HZ1
#
# ASH and GLH are the same operation on different atom names, so both are built;
# LYN is not, and the reason is in OpenMM's data files rather than in the force
# field. Renaming a residue only works if something still knows how to bond it:
#
#   data/pdbNames.xml   <Residue name="ASP" alt1="ASH">
#                       <Residue name="GLU" alt1="GLH">   (no LYN entry)
#   data/residues.xml   ASP carries <Bond from="HD2" to="OD2"/>
#                       GLU carries <Bond from="HE2" to="OE2"/>
#                       ASH, GLH and LYN have no entry at all
#
# So an ASH in a file is renamed back to ASP on load and then bonded by ASP's
# definitions, which already include the HD2 the charged residue never has -
# the proton is bonded, by name, and the ASH *template* is matched afterwards by
# graph isomorphism. GLH goes the same way through GLU. A residue renamed LYN,
# though, is renamed by nothing and bonded by nothing: it has no residues.xml
# entry, so `createStandardBonds` gives it not one bond, and twenty-one loose
# atoms match no template. Neutralising a lysine therefore has to be done by
# leaving the residue named LYS and letting the graph match LYN, which is what
# `make_water_nuc.py --lys` does. It is a different operation, not this one, so
# `--to=LYN` refuses and says so rather than writing a structure that dies at
# `createSystem`.
#
# `parent` is the charged residue the input must currently be: refusing anything
# else is the point, because a residue that is already ASH, or is an ASN that
# someone mistyped the selector for, would otherwise be silently mangled.
# A proton placed from three heavy atoms knows nothing about the rest of the
# structure, so it can land inside a neighbour. 1.2 A is below any real contact:
# a strong O-H...O hydrogen bond is 1.6-2.0 A, and the shortest H...H in a
# protein structure is about 1.7 A, so anything under 1.2 A is an overlap and
# not a close contact. Reported either way - `subst_start.pdb` carried a 0.58 A
# overlap for weeks because nothing measured it (D152).
CLASH_A = 1.2


def proton_contacts(atoms, xyz, exclude, n=3):
    """The `n` nearest atoms to `xyz` outside the residue, nearest first.

    `exclude` is the set of line indices the residue occupies, so the proton is
    not reported as clashing with the oxygen it was just built on.
    """
    near = sorted((_norm(_sub(xyz, a[5])), i) for i, a in enumerate(atoms)
                  if a[0] not in exclude)
    return [(d, atoms[i]) for d, i in near[:n]]


PROTON_VARIANTS = {
    "ASH": {"parent": "ASP", "add": "HD2", "oxygen": "OD2", "carbon": "CG",
            "other": "OD1"},
    "GLH": {"parent": "GLU", "add": "HE2", "oxygen": "OE2", "carbon": "CD",
            "other": "OE1"},
    "LYN": {"parent": "LYS", "unsupported":
            "a neutral lysine cannot be built by renaming the residue. OpenMM "
            "has no LYN entry in residues.xml and no LYN alternative name in "
            "pdbNames.xml, so a residue renamed LYN keeps that name and comes "
            "out of createStandardBonds with no bonds at all - twenty-one "
            "loose atoms, matching no template. Remove HZ1 (ff14SB's LYN is "
            "LYS without HZ1) and leave the residue named LYS instead: OpenMM "
            "matches the template by bond graph and finds LYN anyway. "
            "`make_water_nuc.py --lys CHAIN:RESID` already does exactly that."},
}

OH_LENGTH = 0.98        # A. A carboxylic O-H; ff14SB's OH/HO pair sits here.
COH_ANGLE = 113.0       # degrees, C-O-H. The sp2 oxygen's angle, not 109.5.


def _sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _norm(a):
    return math.sqrt(_dot(a, a))


def _unit(a):
    length = _norm(a)
    if length < 1e-8:
        raise ValueError("cannot normalise a zero-length vector")
    return tuple(x / length for x in a)


def place_carboxyl_proton(carbon, other_o, oxygen,
                          length=OH_LENGTH, angle_deg=COH_ANGLE):
    """The proton on `oxygen`: in the carboxylate plane, anti to `other_o`.

    Everything comes from the three heavy atoms, so the result follows the
    residue wherever it sits in the box and whatever it has rotated to - no
    coordinate here is written down.

    `e1` runs from the hydroxyl oxygen to the carbon, `n` is the normal of the
    C/O/O plane, and `e2` is the in-plane direction perpendicular to `e1`, signed
    so that it points towards the *other* oxygen. Rotating `e1` by the C-O-H
    angle away from `e2` therefore lands the proton in the plane, at that angle,
    and on the far side of the C-O axis from the carbonyl oxygen: the anti
    conformer. Rotating towards `e2` would give syn, which is a real conformer
    but the less common one, and at a 1,3 separation across the carboxylate it
    puts the proton about 2.2 A from the carbonyl oxygen instead of 3.2 - a
    contact a minimisation has to undo before anything else happens.

    Raises on a degenerate (collinear) carboxylate, where the plane - and so the
    whole construction - is not defined.
    """
    e1 = _unit(_sub(carbon, oxygen))
    to_other = _sub(other_o, oxygen)
    n = _cross(e1, to_other)
    if _norm(n) < 1e-6:
        raise ValueError(
            "the three heavy atoms are collinear, so the carboxylate plane is "
            "undefined and there is no in-plane direction to build along")
    e2 = _unit(_cross(_unit(n), e1))
    if _dot(e2, to_other) < 0.0:
        e2 = tuple(-x for x in e2)
    theta = math.radians(angle_deg)
    direction = tuple(math.cos(theta) * a - math.sin(theta) * b
                      for a, b in zip(e1, e2))
    return tuple(o + length * d for o, d in zip(oxygen, direction))


def _pdb_atoms(lines):
    """[(index, name, resname, chain, resseq, (x, y, z)), ...] over ATOM lines."""
    out = []
    for index, line in enumerate(lines):
        if line.startswith(("ATOM", "HETATM")):
            out.append((index, line[12:16].strip(), line[17:20].strip(),
                        line[21], line[22:26].strip(),
                        (float(line[30:38]), float(line[38:46]),
                         float(line[46:54]))))
    return out


def _rename(line, variant):
    """The residue name field replaced, the rest of the line untouched."""
    return line[:17] + "%3s" % variant + line[20:]


def _hydrogen_line(model, name, xyz):
    """A new hydrogen's ATOM line, built from a line of the same residue.

    Inheriting chain, residue number, insertion code, occupancy and B-factor
    from a sibling rather than reconstructing them is what keeps the new atom in
    the same residue. The element columns (77-78) are then overwritten: the
    model is a heavy atom, and `PDBFile` believes that column over the atom
    name, so leaving it would declare the proton an oxygen - the mistake that
    once turned a water into an O-H-O that matched no template.
    """
    line = (model[:12] + "%-4s" % name + model[16]
            + model[17:30]
            + "%8.3f%8.3f%8.3f" % xyz
            + model[54:]).rstrip("\n").ljust(78)
    return line[:76] + " H" + line[78:] + "\n"


def mode_protonate(site, args):
    """Write the structure again with one residue in its protonated form.

        python site.py protonate in.pdb A:92 --to=ASH [--out=out.pdb]

    The structure is copied through line by line and only the named residue is
    touched: every other atom, the ordering, CRYST1, the waters, the metal and
    the CONECT records all arrive unchanged, because the question being asked is
    what one proton does and anything else that moved would be a confound.

    The proton's serial is one past the file's highest, so no existing line is
    renumbered and no existing CONECT has to be rewritten. `PDBFile` keys CONECT
    on serial numbers rather than on file order, so a serial out of sequence is
    read correctly; a duplicate one would not be. The proton is written
    immediately after its oxygen, inside the residue, because `PDBFile` groups
    atoms into residues by file order and an atom appended after the chain would
    become a residue of its own.

    The new bond is also written as a CONECT record. OpenMM does not need it -
    `residues.xml` gives ASP an HD2-OD2 bond and GLU an HE2-OE2 bond already, so
    the proton is bonded by name once the residue has been renamed back to ASP
    or GLU on load - but it costs one line, duplicate bonds are discarded, and
    it means the file states its own bond rather than relying on a data file
    elsewhere happening to carry an atom the charged residue never has. The
    failure it guards against is silent in the PDB and loud hours later: an
    unbonded proton matches no template and dies in `createSystem`.
    """
    if len(args) != 2:
        raise SystemExit("protonate needs <pdb> <chain:resid> --to=ASH")
    src, selector = args
    variant = (OPTS.get("to") or "").upper()
    if variant not in PROTON_VARIANTS:
        raise SystemExit(
            "--to= wants one of %s, not %r"
            % (", ".join(sorted(PROTON_VARIANTS)), OPTS.get("to")))
    spec = PROTON_VARIANTS[variant]
    if "unsupported" in spec:
        raise SystemExit("--to=%s is not implemented: %s"
                         % (variant, spec["unsupported"]))
    parts = selector.split(":")
    if len(parts) != 2 or not all(parts):
        raise SystemExit("the residue wants CHAIN:RESID, not %r" % selector)
    chain, resid = parts

    out = OPTS.get("out") or "%s_%s%s_%s.pdb" % (
        src.rsplit(".", 1)[0], chain, resid, variant)
    if os.path.abspath(out) == os.path.abspath(src):
        raise SystemExit(
            "--out is the input file. `protonate` writes a second structure so "
            "the two protonation states can be compared; it does not edit one "
            "in place. Nothing was written.")

    lines = open(src).read().splitlines(keepends=True)
    atoms = _pdb_atoms(lines)
    if not atoms:
        raise SystemExit("%s holds no ATOM or HETATM records" % src)
    residue = [a for a in atoms if a[3] == chain and a[4] == resid]
    if not residue:
        raise SystemExit(
            "no residue %s:%s in %s. Nothing was written." % (chain, resid, src))
    names = {a[1]: a for a in residue}
    resnames = {a[2] for a in residue}
    if resnames != {spec["parent"]}:
        raise SystemExit(
            "%s:%s is %s, and %s is built from %s. Nothing was written."
            % (chain, resid, "/".join(sorted(resnames)), variant,
               spec["parent"]))

    proton = spec["add"]
    if proton in names:
        raise SystemExit(
            "%s:%s already carries %s, so it is already the protonated form. "
            "Nothing was written." % (chain, resid, proton))
    needed = [spec["carbon"], spec["other"], spec["oxygen"]]
    absent = [n for n in needed if n not in names]
    if absent:
        raise SystemExit(
            "%s:%s is missing %s, which %s is built from. It holds: %s. "
            "Nothing was written."
            % (chain, resid, ", ".join(absent), proton, " ".join(sorted(names))))
    try:
        xyz = place_carboxyl_proton(*[names[n][5] for n in needed])
    except ValueError as exc:
        raise SystemExit("cannot place %s on %s:%s: %s. Nothing was written."
                         % (proton, chain, resid, exc))
    anchor = names[spec["oxygen"]]

    serials = [int(lines[a[0]][6:11]) for a in atoms
               if lines[a[0]][6:11].strip()]
    serial = max(serials or [len(atoms)]) + 1
    if serial > 99999:
        raise SystemExit(
            "the highest serial in %s is %d, so a new atom cannot be numbered "
            "inside the five columns a PDB gives it. This structure needs a "
            "format that can hold it (mmCIF) rather than a proton squeezed "
            "into a field that has run out. Nothing was written."
            % (src, serial - 1))

    residue_rows = {a[0] for a in residue}
    # Measured before anything is written: a clash is a reason not to have a
    # structure, not something to find out after optimising it.
    contacts = proton_contacts(atoms, xyz, residue_rows)
    if contacts and contacts[0][0] < CLASH_A:
        d, other = contacts[0]
        raise SystemExit(
            "%s on %s:%s lands %.2f A from %s:%s:%s, under %.1f A, so it is "
            "inside that atom rather than beside it. The placement is geometric "
            "and does no search, so this frame cannot carry the proton in the "
            "anti position; pick another frame or relax the neighbour first. "
            "Nothing was written."
            % (proton, chain, resid, d, other[3], other[4], other[1], CLASH_A))

    body = []
    for index, line in enumerate(lines):
        if index in residue_rows:
            line = _rename(line, variant)
        body.append(line)
        if index == anchor[0]:
            new = _hydrogen_line(line, proton, xyz)
            body.append(new[:6] + "%5d" % serial + new[11:])
    # Before END: readers stop there, so a CONECT appended after it is dropped
    # without complaint.
    anchor_serial = int(lines[anchor[0]][6:11])
    end = next((i for i, l in enumerate(body) if l.startswith("END")),
               len(body))
    body[end:end] = ["CONECT%5d%5d\n" % (anchor_serial, serial)]

    with open(out, "w") as fh:
        fh.write("".join(body))

    print("%s:%s %s -> %s, written to %s"
          % (chain, resid, spec["parent"], variant, out))
    o, c, other = (names[spec["oxygen"]][5], names[spec["carbon"]][5],
                   names[spec["other"]][5])
    n = _unit(_cross(_sub(c, o), _sub(other, o)))
    print("  %s at %.3f %.3f %.3f, serial %d, CONECT to %s"
          % (proton, xyz[0], xyz[1], xyz[2], serial, spec["oxygen"]))
    print("  %s-%s %.4f A, %s-%s-%s %.2f deg, out of plane %.4f A"
          % (spec["oxygen"], proton, _norm(_sub(xyz, o)),
             spec["carbon"], spec["oxygen"], proton,
             math.degrees(math.acos(max(-1.0, min(1.0, _dot(
                 _unit(_sub(c, o)), _unit(_sub(xyz, o))))))),
             abs(_dot(_sub(xyz, o), n))))
    print("  %s-%s %.3f A, anti (syn would be %.3f A)"
          % (proton, spec["other"], _norm(_sub(xyz, other)),
             _norm(_sub(_syn(c, other, o), other))))
    print("  nearest outside the residue: "
          + ", ".join("%.2f A to %s:%s:%s" % (d, a[3], a[4], a[1])
                      for d, a in contacts))
    print("  the residue's atom set is now exactly ff14SB's %s template; the "
          "residue is neutral where %s was -1, so the system's total charge "
          "rises by 1" % (variant, spec["parent"]))


def _syn(carbon, other_o, oxygen):
    """The syn position, for reporting only: the proton turned to face `other_o`.

    `protonate` never writes this. It is printed beside the anti distance so the
    line says which conformer was built in a way that can be checked, rather
    than asserting it.
    """
    return place_carboxyl_proton(carbon, other_o, oxygen,
                                 angle_deg=-COH_ANGLE)


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
    """Accept or reject a candidate structure against one named hypothesis.

    This mode exists because of D123: the pose everything ran on from October
    had the substrate bound the other way round, with a substituent oxygen on
    the metal and the carbonyl 5.43 A away, and the number sat in the decision
    log for two weeks without being read as a defect. Nothing here is new
    knowledge - every criterion was already measured and recorded. What was
    missing was a command that says yes or no, so a start structure is not
    accepted by eye. Exits non-zero if anything required at this stage fails.

    `--hypothesis=<name>` picks one set from `michaelis_hypotheses`, because a
    criterion is only ever a criterion WITHIN a mechanism. D125 wrote "the
    closest substrate oxygen must be the beta-lactam carbonyl" as though it
    were a fact about Michaelis complexes; it is a fact about the
    carbonyl-on-metal mechanism, which is the thing under test, and D126 found
    a candidate complex that violates it while satisfying the attack geometry.
    A gate that rejects a competing hypothesis by assuming the first one is not
    a gate, it is a way of never testing anything. So each hypothesis carries
    its own criteria and is named in the output, and a structure is rejected
    only against a stated mechanism.
    """
    stage = OPTS.get("at", "start")
    sets = site.get("michaelis_hypotheses") or {}
    name = OPTS.get("hypothesis")
    if name is not None:
        if name not in sets:
            raise SystemExit("no hypothesis %r; the site file has: %s"
                             % (name, ", ".join(sorted(sets)) or "none"))
        criteria = sets[name]
    elif sets:
        name = sorted(sets)[0] if len(sets) == 1 else None
        if name is None:
            raise SystemExit("several hypotheses (%s): say which with "
                             "--hypothesis=<name>" % ", ".join(sorted(sets)))
        criteria = sets[name]
    else:
        name, criteria = "michaelis_criteria", site.get("michaelis_criteria") or []
    if not criteria:
        raise SystemExit("no criteria in the site file: nothing to check")
    worst = 0
    for path in args:
        at = load_pdb(path)
        print("==", path, "   hypothesis '%s', required at: %s" % (name, stage))
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
        print("   %s under '%s'"
              % ("REJECTED" if worst else "accepted at stage '%s'" % stage, name))
    raise SystemExit(worst)


MODES = {"sphere": mode_sphere, "nucleophiles": mode_nucleophiles,
         "acceptors": mode_acceptors, "joint": mode_joint, "poses": mode_poses,
         "graft": mode_graft, "shell": mode_shell,
         "contacts": mode_contacts, "check": mode_check,
         "protonate": mode_protonate}

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
    if mode in ("sphere", "nucleophiles", "poses", "graft", "check",
                "protonate"):
        MODES[mode](site, rest)
    else:
        if len(rest) != 2:
            raise SystemExit("%s needs <trajectory> <topology>" % mode)
        MODES[mode](site, rest[0], rest[1])
