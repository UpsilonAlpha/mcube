"""D50: 25 ps QM/MM MD from the corrected resting-state reactant.

D49 found that in the original 50 ps trajectory no water ever both reaches the
Burgi-Dunitz window and donates to Asp92. But that trajectory was generated in
the degraded sphere, with Asp92 coordinated to the zinc at 2.07 A and Lys226
already deprotonated, so it may simply be the wrong state to ask.

This re-runs the sampling from `protdft/cleanR/ohH/QM_MM_optimized.pdb`, the
resting state established in D48: two His, O01 as water, HOH485, Asp92 free at
2.80 A and Lys226 neutral.

No constraints. The QM/MM MD applies no MM restraints, so the metal site is
described by xTB alone and is free to find its own arrangement; whether Asp92
returns to the zinc is then a result rather than an assumption. (In CphA,
Asp120 is both a zinc ligand and the general base, so it coordinating is not
by itself disqualifying.)

Afterwards, rerun `joint_nac.py` on the new trajectory.
"""
import os
from pathlib import Path

import yaml

from mcube.preprocess import _locked

SECTION = "qmmm_md2_ohH"
inputs = dict(
    skip_equilibrate=False,   # this stage is the MD
    skip_optimize=True,       # no optimisation afterwards: frames are what we want
    singlepoint=False,
    time=25,                  # ps; ~9.5 h at the measured 1.37 s/step
    timestep=0.001,
    traj_frequency=50,        # 500 frames
    temp=300,
    cores=4,                  # the benchmarked optimum with xTB thread binding
)
c = Path("config.yaml")
with _locked(c):
    d = yaml.safe_load(open(c))
    d[SECTION] = {"inputs": inputs}
    yaml.safe_dump(d, open("t.tmp", "w"), sort_keys=False, default_flow_style=False)
    os.replace("t.tmp", c)

start = Path("protdft/cleanR/ohH/QM_MM_optimized.pdb")
print(SECTION, "->", inputs["time"], "ps from", start)
print("start structure exists:", start.exists())
