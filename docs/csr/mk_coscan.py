"""D58: can the site open by releasing the substrate carbonyl, and does that help?

The Zn-Asp92 scan (D57) showed Asp92 is a real ligand - 9.1 kcal/mol at DFT to
pull it to 2.80 A - and, more importantly, that removing it does **not** bring
the nucleophile closer: O01-C4x went 3.62 to 3.70 A across the whole scan. So
opening that particular coordination site does not make the site competent.

The carbonyl is the other removable ligand, at Zn-O2x 2.13 A in the relaxed
structure. Releasing it is different in kind from releasing Asp92, because it
frees the substrate to rotate: the two Zn-bound waters sit 34-36 degrees off the
Burgi-Dunitz line precisely because they are cis to a carbonyl pinned on the
same metal. The question is whether the attack geometry improves as the carbonyl
leaves, and at what energetic cost - against the obvious trade-off that leaving
also gives up the Lewis-acid activation of the electrophile.

Monitored along the way: both Zn-water distances to C4x, the rest of the first
shell, and the lactam C4x-N1x bond. Angles are not available through `monitor`,
which only does distances, so run docs/csr/zn_sphere.py over the scan_*.pdb
files afterwards to get the attack angles.

Same start structure as the Asp92 scan (protdft/zncheck/relax/QM_MM_optimized.pdb)
so the two coordinates are directly comparable.
"""
import os
from pathlib import Path

import yaml

from mcube.preprocess import _locked

MON = [["B:1:ZN", "A:92:OD2"], ["B:1:ZN", "C:2:O01"], ["B:1:ZN", "D:1:O"],
       ["B:1:ZN", "A:88:NE2"], ["B:1:ZN", "A:164:NE2"],
       ["C:2:O01", "C:1:C4x"], ["D:1:O", "C:1:C4x"], ["C:1:C4x", "C:1:N1x"]]
SECTION = "scan_znco"
inputs = dict(atoms=["B:1:ZN", "C:1:O2x"], start=2.15, end=4.00, points=9,
              opt_maxiter=300, monitor=MON, resume=True)

c = Path("config.yaml")
with _locked(c):
    d = yaml.safe_load(open(c))
    d[SECTION] = {"inputs": inputs}
    yaml.safe_dump(d, open("t.tmp", "w"), sort_keys=False, default_flow_style=False)
    os.replace("t.tmp", c)

start = "protdft/zncheck/relax/QM_MM_optimized.pdb"
at = {}
for i, l in enumerate(x for x in open(start) if x.startswith(("ATOM", "HETATM"))):
    at[f"{l[21]}:{l[22:26].strip()}:{l[12:16].strip()}"] = i
sels = inputs["atoms"] + [a for p in MON for a in p]
missing = [s for s in sels if s not in at]
print(SECTION, inputs["atoms"], "%.2f -> %.2f A in %d points"
      % (inputs["start"], inputs["end"], inputs["points"]))
print("start:", start, "exists =", Path(start).exists())
print("MISSING SELECTORS:", missing if missing else "none")
