"""D45: time a DFT/MM gradient on a GPU, against the measured CPU cost.

Three sections, all starting from the same NVT end point as D42, all doing a few
constrained QM/MM optimisation steps with no single point, so the per-step wall
time in optimization_energies.log IS the DFT/MM gradient cost:

  gpu_wb97x  wB97X/def2-SVP on GPU   - directly comparable to the CPU measurement
                                       (~4-5 h per energy on 16 cores, job 29997628)
  gpu_r2scan r2SCAN/def2-SVP on GPU  - the proposed geometry level (meta-GGA, no
                                       exact exchange, so it should be much cheaper)
  gpu_r2scan_pc as above with GPU_pcgrad, to see what moving the MM point-charge
                                       gradient onto the device is worth

The same two sphere constraints as D42 are kept so the geometry is the reactant we
actually care about, not an arbitrary structure.
"""
import os
from pathlib import Path

import yaml

from mcube.preprocess import _locked

SPHERE = [["B:1:ZN", "A:92:OD2", 2.80], ["B:1:ZN", "C:1:O2x", 2.90]]
BASE = dict(skip_equilibrate=True, singlepoint=False, opt_maxiter=3,
            opt_constraints=SPHERE, qm_engine="dft", dft_platform="GPU",
            dft_basis="def2-svp", dft_gridlevel=3, dft_densityfit=True)
sections = {
    "qmmm_gpu_wb97x": dict(BASE, dft_functional="wb97x"),
    "qmmm_gpu_r2scan": dict(BASE, dft_functional="r2scan", dft_scf_dispersion="D4"),
    "qmmm_gpu_r2scan_pc": dict(BASE, dft_functional="r2scan", dft_scf_dispersion="D4",
                               dft_gpu_pcgrad=True),
}
c = Path("config.yaml")
with _locked(c):
    d = yaml.safe_load(open(c))
    for name, inp in sections.items():
        d[name] = {"inputs": inp}
    yaml.safe_dump(d, open("t.tmp", "w"), sort_keys=False, default_flow_style=False)
    os.replace("t.tmp", c)
for name, inp in sections.items():
    print(name, inp["dft_functional"], "pcgrad", inp.get("dft_gpu_pcgrad", False))
