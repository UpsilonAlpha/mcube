"""D53: is the zinc really six-coordinate, or is that an artefact?

Everything in D52 rests on the zinc holding both Asp92 and the substrate
carbonyl in addition to two waters. Yang et al., J Chem Inf Model 2021,
61:5658 (PMID 34748329) found that classical MM MD "overfavors octahedral
geometries" by pulling extra waters into the zinc shell, and recommended
DFTB3 QM/MM MD followed by DFT QM/MM optimisation to correct it. Our pipeline
equilibrated with MM first, with preprocess metal restraints on everything
within 3.0 A of the zinc, so this is exactly the documented failure mode.

The test: relax the structure, walk Zn-Asp92 OD2 out from 2.05 to 3.05 A at
xTB/MM, then recompute a few points at DFT/MM. If xTB's minimum sits near 2.1
and DFT's near 2.8, the six-coordinate site is an artefact and the mechanism
built on it has to be revisited. Zn-O2x, Zn-O01 and Zn-HOH485 are monitored so
the whole first shell is visible along the coordinate.

Sections created:
  qmmm_znrelax   - plain xTB/MM optimisation of the md2 end point
  scan_znasp     - the 9-point relaxed scan
  qmmm_sp_dft    - a single point only, reused at several scan geometries
"""
import os
from pathlib import Path

import yaml

from mcube.preprocess import _locked

MON = [["B:1:ZN", "C:1:O2x"], ["B:1:ZN", "C:2:O01"], ["B:1:ZN", "D:1:O"],
       ["B:1:ZN", "A:88:NE2"], ["B:1:ZN", "A:164:NE2"], ["C:2:O01", "C:1:C4x"]]
sections = {
    "qmmm_znrelax": dict(skip_equilibrate=True, singlepoint=False, opt_maxiter=400),
    "scan_znasp": dict(atoms=["B:1:ZN", "A:92:OD2"], start=2.05, end=3.05, points=9,
                       opt_maxiter=300, monitor=MON, resume=True),
}
# One section per scan point. They must be separate, because mcube writes each
# run's results back into its own section and a shared one would be overwritten
# three times over, leaving only the last energy.
for point in ("00", "03", "06", "08"):
    sections[f"qmmm_sp_dft_{point}"] = dict(
        # skip both the MD and the optimisation: score the geometry it is given
        skip_equilibrate=True, skip_optimize=True, singlepoint=True, sp_cores=16)
c = Path("config.yaml")
with _locked(c):
    d = yaml.safe_load(open(c))
    for name, inp in sections.items():
        d[name] = {"inputs": inp}
    yaml.safe_dump(d, open("t.tmp", "w"), sort_keys=False, default_flow_style=False)
    os.replace("t.tmp", c)

at = {}
start = "protdft/md2_ohH/QM_MM_lastframe.pdb"
for i, l in enumerate(x for x in open(start) if x.startswith(("ATOM", "HETATM"))):
    at[f"{l[21]}:{l[22:26].strip()}:{l[12:16].strip()}"] = i
sels = sections["scan_znasp"]["atoms"] + [a for p in MON for a in p]
missing = [s for s in sels if s not in at]
for name in sections:
    print(name)
print("scan coordinate:", sels[0], at.get(sels[0]), "-", sels[1], at.get(sels[1]))
print("MISSING SELECTORS:", missing if missing else "none")
