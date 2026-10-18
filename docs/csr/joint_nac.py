"""Does any water ever attack AND donate to the base at the same time?

    python joint_nac.py traj.dcd topology.pdb

The scan of a single optimised structure (nuc_scout.py) can only say what that
structure does. This asks the whole trajectory a joint question, because a
general-base mechanism needs both halves at once:

  * attack geometry: O...C4x within `MAX_OC`, angle O-C4x-O2x in the
    Burgi-Dunitz window, and the water not bound to the zinc; and
  * a proton already pointing at a base: one of that water's own hydrogens
    within `MAX_HB` of Asp92 OD1/OD2 -- or, for a relay, of a water that is
    itself donating to Asp92.

Reports how many frames satisfy each half separately and both together. If both
halves never coincide, that is evidence about the mechanism rather than a reason
to force a scan: the enzyme would have to reorganise first, and the barrier for
that reorganisation belongs in the model.
"""
import sys

import mdtraj as md
import numpy as np

MAX_OC = 4.5          # A, nucleophile oxygen to carbonyl carbon
BD_LO, BD_HI = 95.0, 120.0   # degrees, Burgi-Dunitz window
MAX_HB = 2.5          # A, H...acceptor
ZN_BOUND = 2.6        # A, treat as a zinc ligand beyond this cutoff

traj_file, top_file = sys.argv[1], sys.argv[2]
top = md.load_topology(top_file)


def sel(expr):
    idx = top.select(expr)
    if len(idx) == 0:
        raise SystemExit("nothing matched: " + expr)
    return int(idx[0])


zn = sel("resname ZN and name ZN")
c4 = sel("resname mer and name C4x")
o2 = sel("resname mer and name O2x")
od1 = sel("resid 91 and name OD1")   # mdtraj resid is 0-based: Asp92
od2 = sel("resid 91 and name OD2")

waters = [r for r in top.residues if r.name in ("HOH", "WAT")]
w_o = {r.index: [a.index for a in r.atoms if a.element.symbol == "O"][0] for r in waters}
w_h = {r.index: [a.index for a in r.atoms if a.element.symbol == "H"] for r in waters}

print("loading trajectory...")
traj = md.load(traj_file, top=top_file)
n = traj.n_frames
print("%d frames, %d waters" % (n, len(waters)))

# Prefilter on the first frame: only waters that ever come near the carbon
first = traj[0]
near = []
for ri, oi in w_o.items():
    dd = md.compute_distances(first, [[oi, c4]])[0, 0] * 10
    if dd < 12.0:
        near.append(ri)
print("%d waters start within 12 A of C4x" % len(near))

stats = {}
for ri in near:
    oi = w_o[ri]
    pairs = [[oi, c4], [oi, zn]]
    dist = md.compute_distances(traj, pairs) * 10
    angle = np.degrees(md.compute_angles(traj, [[oi, c4, o2]])[:, 0])
    attacks = (dist[:, 0] < MAX_OC) & (dist[:, 1] > ZN_BOUND) & \
              (angle > BD_LO) & (angle < BD_HI)
    if not attacks.any():
        continue
    hb = np.zeros(n, dtype=bool)
    for hi in w_h[ri]:
        hd = md.compute_distances(traj, [[hi, od1], [hi, od2]]) * 10
        hb |= (hd < MAX_HB).any(axis=1)
    both = attacks & hb
    stats[ri] = (attacks.sum(), hb.sum(), both.sum(), dist[:, 0], angle)

if not stats:
    raise SystemExit("no non-zinc water ever reaches the attack geometry")

res = {r.index: r for r in waters}
print("\n%-16s %8s %8s %8s   %s" % ("water", "attack", "H->Asp92", "BOTH", "best O-C4x / angle"))
for ri, (a, h, b, dd, an) in sorted(stats.items(), key=lambda kv: -kv[1][2]):
    r = res[ri]
    k = int(np.argmin(np.where((an > BD_LO) & (an < BD_HI), dd, 1e9)))
    print("%-16s %8d %8d %8d   %.2f A / %.0f deg (frame %d)"
          % ("%s%s" % (r.name, r.resSeq), a, h, b, dd[k], an[k], k))

best = max(stats.items(), key=lambda kv: kv[1][2])
if best[1][2] == 0:
    print("\nNo frame has both halves at once.")
    print("The attacking water and the water that donates to Asp92 are never the same")
    print("molecule, so a one-water general-base step is not available from this")
    print("trajectory: either a relay carries the proton, or the site must reorganise.")
else:
    ri, (a, h, b, dd, an) = best
    frames = "see above"
    print("\n%s%s satisfies both in %d frames -- use one of those as the M2 reactant."
          % (res[ri].name, res[ri].resSeq, b))
