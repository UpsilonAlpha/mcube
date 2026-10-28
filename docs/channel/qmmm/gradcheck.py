#!/usr/bin/env python3
"""
Finite-difference check of the QM/MM gradient of one exchange.py job (decision log C22).

The DFT/MM optimisations climbed 30-45 kcal/mol after an early minimum while geomeTRIC rated the
steps as good, which is what an energy and a gradient that do not belong together look like.
This builds the QMMMTheory exactly as `mcube qmmm` does from the job's config.yaml (same MM
system, QM region, metal parameters and QM level), computes the analytic gradient once, and
compares chosen Cartesian components with central differences of the energy. Each energy is also
split into ASH's QM and MM parts, so the term that disagrees can be named.

Atoms (0-based) default to: the ion, the nearest QM water O, the nearest MM water O, an MM water O
about 6 A away (inside the 8 A active region) and one about 12 A away (outside it).

  python gradcheck.py --job $S/channel/qmmm/water/LA_f0024/LA --h 0.005 --out gradcheck.json

Run it in the job directory's environment (exchange.slurm's modules); cost is 1 + 6 per atom QM/MM
energies at the job's optimisation level.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

BOHR = 0.52917721092          # A per Bohr
HARTREE_KCAL = 627.509474


def build(job: Path, cores: int):
    from ash import Fragment
    from mcube import qmmm as Q
    config = str(job / "config.yaml")
    preprocess = Q.load_section(config, "preprocess")
    md = Q.load_section(config, "md")
    inputs = Q.stage_inputs(config, "qmmm", "qmmm")
    pdb_file, _ = Q.resolve_start_structure(config, "qmmm", inputs, [])
    forcefields, implicit = Q.mm_settings(inputs, md, config)
    mm = Q.make_mm(forcefields, pdb_file, implicit, cores, inputs.get("platform", "CPU"),
                   inputs.get("properties"))
    inputs = Q.resolve_stage_atoms(inputs, mm)
    fragment = Fragment(pdbfile=pdb_file)
    qm_atoms, boundary = Q.resolve_qm_region(inputs, preprocess, pdb_file)
    Q.apply_stage_metal_parameters(inputs, preprocess, md, mm, forcefields, qm_atoms)
    theory = Q.make_qmmm(inputs, mm, fragment, qm_atoms, boundary, cores)
    return theory, fragment, inputs, qm_atoms, pdb_file


def evaluate(theory, fragment, coords, inputs, grad=False):
    out = theory.run(current_coords=coords, elems=fragment.elems, Grad=grad,
                     charge=inputs.get("charge"), mult=inputs.get("mult"))
    energy, gradient = (out if grad else (out, None))
    parts = {k: float(getattr(theory, k)) for k in ("QMenergy", "MMenergy")
             if isinstance(getattr(theory, k, None), (int, float, np.floating))}
    return float(energy), gradient, parts


def default_atoms(fragment, qm_atoms, pdb_file):
    from openmm.app import PDBFile
    top = PDBFile(pdb_file).topology
    xyz = np.asarray(fragment.coords)
    qm = set(int(a) for a in qm_atoms)
    ion = [a.index for a in top.atoms() if a.residue.name in ("CA", "LA", "LU")
           and len(list(a.residue.atoms())) == 1][0]
    water_o = [a.index for a in top.atoms() if a.residue.name in ("HOH", "WAT") and a.element.symbol == "O"]
    d = {i: float(np.linalg.norm(xyz[i] - xyz[ion])) for i in water_o}
    qm_w = min((i for i in water_o if i in qm), key=d.get)
    mm_w = sorted((i for i in water_o if i not in qm), key=d.get)
    near = mm_w[0]
    mid = min(mm_w, key=lambda i: abs(d[i] - 6.0))
    far = min(mm_w, key=lambda i: abs(d[i] - 12.0))
    return {"ion": ion, "qm_water_O": qm_w, "mm_water_O_near": near, "mm_water_O_6A": mid,
            "mm_water_O_12A": far}, d


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--job", required=True)
    p.add_argument("--h", type=float, default=0.005, help="displacement in A")
    p.add_argument("--atoms", help="comma-separated 0-based indices (default: see above)")
    p.add_argument("--cores", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    p.add_argument("--out", default="gradcheck.json")
    args = p.parse_args()

    job = Path(args.job).resolve()
    work = job / "gradcheck"           # keep ASH/PySCF scratch files out of the job directory
    work.mkdir(exist_ok=True)
    os.chdir(work)
    theory, fragment, inputs, qm_atoms, pdb_file = build(job, args.cores)
    x0 = np.array(fragment.coords, dtype=float)
    if args.atoms:
        atoms = {f"atom{i}": int(i) for i in args.atoms.split(",")}
        dist = {}
    else:
        atoms, dist = default_atoms(fragment, qm_atoms, pdb_file)
    e0, g0, parts0 = evaluate(theory, fragment, x0, inputs, grad=True)
    g0 = np.asarray(g0)
    report = {"job": str(job), "pdb_file": pdb_file, "h_A": args.h, "E0": e0, "parts0": parts0,
              "n_qm": len(qm_atoms), "atoms": {}}
    print(f"E0 {e0:.10f} Eh  parts {parts0}", flush=True)
    for label, i in atoms.items():
        rows = []
        for c in range(3):
            e = {}
            for s in (+1, -1):
                x = x0.copy(); x[i, c] += s * args.h
                e[s] = evaluate(theory, fragment, x, inputs)
            fd = (e[+1][0] - e[-1][0]) / (2 * args.h)                     # Eh/A
            an = float(g0[i, c]) / BOHR                                     # Eh/Bohr -> Eh/A
            fd_parts = {k: (e[+1][2][k] - e[-1][2][k]) / (2 * args.h) for k in e[+1][2] if k in e[-1][2]}
            rows.append({"comp": "xyz"[c], "analytic": an, "fd": fd, "diff": an - fd,
                         "fd_parts": fd_parts})
            print(f"{label} ({i}, {'QM' if i in set(qm_atoms) else 'MM'}) {'xyz'[c]}: analytic {an:+.6f}  "
                  f"fd {fd:+.6f}  diff {an - fd:+.2e} Eh/A  fd_parts {fd_parts}", flush=True)
        report["atoms"][label] = {"index": i, "qm": i in set(qm_atoms),
                                  "dist_to_ion_A": round(dist.get(i, float('nan')), 2), "rows": rows}
    worst = max(abs(r["diff"]) for a in report["atoms"].values() for r in a["rows"])
    report["max_abs_diff_Eh_per_A"] = worst
    report["max_abs_diff_kcal_per_A"] = worst * HARTREE_KCAL
    (job / args.out).write_text(json.dumps(report, indent=1))
    print(f"max |analytic - fd| = {worst:.2e} Eh/A ({worst * HARTREE_KCAL:.3f} kcal/mol/A)")


if __name__ == "__main__":
    main()
