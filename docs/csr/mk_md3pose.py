"""D59: is the substrate pose robust, or is it inherited from the MM constraints?

The MM pipeline held Zn-carbonyl with a rigid holonomic constraint at 1.94 A
through all of NPT and NVT, so every structure since descends from a forced
starting point. Releasing the carbonyl by 0.8 A and watching it relax back
(md2) proves little, because that is a local nudge within the same pose.

The carbonyl scan produced something better: relaxed structures with the
substrate genuinely displaced, out to Zn-O2x 4.00 A. Starting unconstrained
QM/MM MD from the far end asks the real question.

  carbonyl returns to the zinc  -> the pose is robust, and the six-coordinate
                                   Michaelis complex is a genuine conclusion
  it settles somewhere else     -> we have a different Michaelis complex, and
                                   the CphA/EcoRV question deserves re-asking
                                   on the better geometry (nuc_scout.py and
                                   joint_nac.py on the new trajectory)

Same MD settings as md2 so the two trajectories are directly comparable.
"""
import os
from pathlib import Path

import yaml

from mcube.preprocess import _locked

SECTION = "qmmm_md3_pose"
inputs = dict(skip_equilibrate=False, skip_optimize=True, singlepoint=False,
              time=25, timestep=0.001, traj_frequency=50, temp=300, cores=4)

c = Path("config.yaml")
with _locked(c):
    d = yaml.safe_load(open(c))
    d[SECTION] = {"inputs": inputs}
    yaml.safe_dump(d, open("t.tmp", "w"), sort_keys=False, default_flow_style=False)
    os.replace("t.tmp", c)

start = Path("protdft/coscan/scan_08.pdb")
print(SECTION, "%d ps from %s" % (inputs["time"], start))
print("start exists:", start.exists())
if not start.exists():
    print("WARNING: the carbonyl scan has not produced its last point yet")
