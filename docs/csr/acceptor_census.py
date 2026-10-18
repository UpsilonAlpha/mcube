"""What can accept the attacking water's proton, counted over the frames where it attacks?

    python acceptor_census.py traj.dcd topology.pdb [resSeq]

joint_nac.py asks a yes/no question about one named base. This asks the open
question instead: during the frames where a water is genuinely on the
Burgi-Dunitz line, what is within hydrogen-bond reach of its protons at all?

That matters because a mechanism does not require the proton transfer to be
concerted with the attack. If the nearest acceptor is another water, the proton
leaves along a solvent chain and those waters belong in the QM region; if it is
a carboxylate or a neutral lysine, that residue is the base. Either way the
answer should come from the structures rather than from the hypothesis.
"""
import collections
import sys

import mdtraj as md
import numpy as np

MAX_OC = 4.5
BD_LO, BD_HI = 95.0, 120.0
MAX_HB = 2.6
ZN_BOUND = 2.6

traj_file, top_file = sys.argv[1], sys.argv[2]
want = sys.argv[3] if len(sys.argv) > 3 else None
top = md.load_topology(top_file)


def one(expr):
    idx = top.select(expr)
    if len(idx) == 0:
        raise SystemExit("nothing matched: " + expr)
    return int(idx[0])


zn, c4, o2 = one("resname ZN and name ZN"), one("resname mer and name C4x"), one("resname mer and name O2x")
traj = md.load(traj_file, top=top_file)
n = traj.n_frames
print("%d frames" % n)

waters = [r for r in top.residues if r.name in ("HOH", "WAT")]
w_o = {r.index: [a.index for a in r.atoms if a.element.symbol == "O"][0] for r in waters}
w_h = {r.index: [a.index for a in r.atoms if a.element.symbol == "H"] for r in waters}
res_of = {r.index: r for r in waters}

# Pick the water that attacks most often, unless one was named
best, best_n, best_mask = None, 0, None
for ri, oi in w_o.items():
    if md.compute_distances(traj[0], [[oi, c4]])[0, 0] * 10 > 12:
        continue
    if want is not None and str(res_of[ri].resSeq) != want:
        continue
    d = md.compute_distances(traj, [[oi, c4], [oi, zn]]) * 10
    a = np.degrees(md.compute_angles(traj, [[oi, c4, o2]])[:, 0])
    mask = (d[:, 0] < MAX_OC) & (d[:, 1] > ZN_BOUND) & (a > BD_LO) & (a < BD_HI)
    if mask.sum() > best_n:
        best, best_n, best_mask = ri, int(mask.sum()), mask
if best is None:
    raise SystemExit("no attacking water found")
r = res_of[best]
print("attacking water: %s%s, on the Burgi-Dunitz line in %d of %d frames\n"
      % (r.name, r.resSeq, best_n, n))

# Every heavy atom that is ever close to either proton during those frames
sub = traj[best_mask]
heavy = [a.index for a in top.atoms
         if a.element.symbol in ("O", "N") and a.residue.index != best]
counts = collections.Counter()
closest = {}
for hi in w_h[best]:
    pairs = [[hi, k] for k in heavy]
    d = md.compute_distances(sub, pairs) * 10
    near = d < MAX_HB
    for col in np.flatnonzero(near.any(axis=0)):
        k = heavy[col]
        a = top.atom(k)
        label = "%s%s %s" % (a.residue.name, a.residue.resSeq, a.name)
        counts[label] += int(near[:, col].sum())
        closest[label] = min(closest.get(label, 9.9), float(d[:, col].min()))

if not counts:
    print("Nothing comes within %.1f A of either proton while it is attacking." % MAX_HB)
else:
    print("%-22s %8s %10s" % ("acceptor", "frames", "closest"))
    for label, c in counts.most_common(12):
        print("%-22s %8d %9.2f A" % (label, c, closest[label]))
    print("\n(out of %d attacking frames)" % best_n)
