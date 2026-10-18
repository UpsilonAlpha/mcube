"""D42: proton-transfer check in a clean four-coordinate sphere.

Both states start from the pre-QM/MM-MD NVT end point, which already has the
canonical reactant: Zn + His88 + His164 + OH-(O01, one H) + H2O(HOH485), with
Lys226 as NH3+ and Asp92 OD2 at its crystallographic 2.79 A.

The same two sphere constraints are applied to BOTH states, so the energy
difference reflects only where the proton sits:
  Zn-Asp92 OD2 = 2.80 A   (crystallographic; stops the migration seen in the MD)
  Zn-mer O2x   = 2.90 A   (carbonyl polarised, not inner-sphere; the preprocess
                           metal restraint had pinned it at 1.94 A, shorter than
                           the Zn-OH- bond)
Each state then pins the shared proton on one side, so neither can collapse into
the other:
  lysH: NZ(Lys226)-HZ1 = 1.03 A   -> proton on the lysine, O01 stays hydroxide
  ohH:  O01-HZ1        = 0.98 A   -> proton on O01, O01 becomes water
Forcing both sides is a change from D30, where ohH was left unconstrained and
relied on the transfer happening by itself; in a clean sphere it may not, which
would have made the two states identical.
"""
import os
from pathlib import Path

import yaml

from mcube.preprocess import _locked

SPHERE = [["B:1:ZN", "A:92:OD2", 2.80], ["B:1:ZN", "C:1:O2x", 2.90]]
states = {
    "qmmm_cleanR_lysH": SPHERE + [["A:226:NZ", "A:226:HZ1", 1.03]],
    "qmmm_cleanR_ohH":  SPHERE + [["C:2:O01", "A:226:HZ1", 0.98]],
}
c = Path("config.yaml")
with _locked(c):
    d = yaml.safe_load(open(c))
    for name, cons in states.items():
        d[name] = {"inputs": dict(skip_equilibrate=True, singlepoint=True,
                                  opt_maxiter=600, opt_constraints=cons,
                                  sp_cores=16)}
    yaml.safe_dump(d, open("config.yaml.edit.tmp", "w"), sort_keys=False, default_flow_style=False)
    os.replace("config.yaml.edit.tmp", c)

# every selector must resolve in the starting structure
at = {}
start = "md_results/NVTtrajectory_lastframe.pdb"
for i, l in enumerate(x for x in open(start) if x.startswith(("ATOM", "HETATM"))):
    at[f"{l[21]}:{l[22:26].strip()}:{l[12:16].strip()}"] = i
bad = []
for name, cons in states.items():
    sels = [a for p in cons for a in p[:2]]
    for s in sels:
        if s not in at:
            bad.append((name, s))
    print(name)
    for p in cons:
        print(f"   hold {p[0]}({at.get(p[0],'?')}) - {p[1]}({at.get(p[1],'?')}) at {p[2]} A")
print("MISSING SELECTORS:", bad if bad else "none")
