"""Active-site geometry analysis, driven by a site description.

    python site.py sphere       <pdb> [pdb ...]     metal sphere and protonation
    python site.py nucleophiles <pdb> [pdb ...]     attack geometry + proton acceptors
    python site.py acceptors    <traj> <topology>   what CAN take the proton, counted
    python site.py joint        <traj> <topology>   attack AND donate, same molecule

Every atom, every candidate base and every criterion comes from `site.yaml`
(override with --site), so this script carries no knowledge of a particular
enzyme and a new system needs a new YAML file rather than new code.

The four modes answer four different questions, in rough order of cost:

  sphere        what is coordinated, and where are the hydrogens?
  nucleophiles  for one structure: which waters could attack, and for each of
                their protons, what is within hydrogen-bond reach?
  acceptors     over a trajectory, during the frames where a water is genuinely
                on the Burgi-Dunitz line, what is near its protons at all? Use
                this when you do not want to assume which base is involved.
  joint         over a trajectory, does any single water both reach the attack
                geometry and donate to a named base, directly or through one
                bridging water? A general-base mechanism needs both at once.
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


# --- entry point -------------------------------------------------------------

MODES = {"sphere": mode_sphere, "nucleophiles": mode_nucleophiles,
         "acceptors": mode_acceptors, "joint": mode_joint}

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    site_file = next((a.split("=", 1)[1] for a in sys.argv[1:]
                      if a.startswith("--site=")), "site.yaml")
    if len(args) < 2 or args[0] not in MODES:
        raise SystemExit(__doc__)
    site = yaml.safe_load(open(site_file))["site"]
    mode, rest = args[0], args[1:]
    if mode in ("sphere", "nucleophiles"):
        MODES[mode](site, rest)
    else:
        if len(rest) != 2:
            raise SystemExit("%s needs <trajectory> <topology>" % mode)
        MODES[mode](site, rest[0], rest[1])
