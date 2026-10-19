"""D55/D56: DFT/MM geometry optimisation of the zinc site, with no xTB in the path.

Cost at this size has surprised us twice, so this creates a short timing probe
alongside the two real optimisations. Run the probe first and read the per-step
wall time out of optimization_energies.log before launching the others.

Level: r2SCAN-D4/def2-TZVP. A meta-GGA has no exact exchange, so it is several
times cheaper than wB97X at the same basis, and it is what the wB97X-3c authors
recommend for geometries; wB97X-3c itself is weaker for 3d transition-metal
geometries and its vDZP basis is in neither PySCF nor the Basis Set Exchange
(D55). Dispersion goes through ASH's own keyword so it is in the gradient.

Two starting coordinations, because a local optimisation cannot leave its basin
and starting from one would prove nothing:

  dftopt_fromMM   restrained MM end point - carbonyl 1.94 A, Asp92 2.79 A
  dftopt_from4c   the four-coordinate state - Asp92 2.80 A, carbonyl 2.90 A

No constraints in either: where Asp92 and the carbonyl end up is the result.
"""
import os
from pathlib import Path

import yaml

from mcube.preprocess import _locked

LEVEL = dict(qm_engine="dft", dft_functional="r2scan", dft_basis="def2-tzvp",
             dft_scf_dispersion="D4", dft_gridlevel=3, dft_densityfit=True,
             dft_cores=16, skip_equilibrate=True, singlepoint=False)
sections = {
    # three steps only: this exists to measure the gradient cost
    "dftopt_probe": dict(LEVEL, opt_maxiter=3),
    "dftopt_fromMM": dict(LEVEL, opt_maxiter=250),
    "dftopt_from4c": dict(LEVEL, opt_maxiter=250),
}
c = Path("config.yaml")
with _locked(c):
    d = yaml.safe_load(open(c))
    for name, inp in sections.items():
        d[name] = {"inputs": inp}
    yaml.safe_dump(d, open("t.tmp", "w"), sort_keys=False, default_flow_style=False)
    os.replace("t.tmp", c)

starts = {"dftopt_fromMM": "md_results/NVTtrajectory_lastframe.pdb",
          "dftopt_from4c": "protdft/cleanR/ohH/QM_MM_optimized.pdb"}
for name in sections:
    print(name, sections[name]["dft_functional"], sections[name]["dft_basis"],
          "maxiter", sections[name]["opt_maxiter"])
for name, path in starts.items():
    print("%-16s <- %s  exists=%s" % (name, path, Path(path).exists()))
