#!/usr/bin/env python3
"""
Relaxed QM/MM distance scan (xTB/MM, same model as the qmmm stage).

The distance between two atoms is stepped from `start` to `end`; at each point
the active region is optimised with that distance fixed, starting from the
previous point's geometry. Used to build intermediates and products inside the
reactant's box (identical atoms for NEB) and to get a first energy profile.

Config (the scan is defined in the section named by --section, default `scan`;
anything not given there is taken from `scan.inputs`):

    scan_attack:
      inputs:
        atoms: ["C:2:O01", "C:1:C4x"]   # i, j: "chain:resid:atom" or 0-based indices
        start: 3.0            # Å; default: the current distance
        end: 1.45             # Å
        points: 16            # including both ends
        constraints: []       # extra bonds held during the scan, [[i, j]] or [[i, j, d]]
        monitor: []           # distances [[i, j], ...] reported at every point
        opt_maxiter: 200

Charge, multiplicity, force fields, xTB settings, cores and the active region
come from qmmm.inputs. Outputs (in the working directory): scan_XX.pdb per
point and scan_energies.csv; results go to the section.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Any

import numpy as np

from .preprocess import save_result, load_section
from .qmmm import (
    resolve_qm_region,
    stage_inputs,
    mm_settings,
    make_mm,
    make_qmmm,
    resolve_actatoms,
    distance,
    check_in_active,
    water_constraints,
    bond_constraints,
    clean_geometric_files,
    write_pdb,
    resolve_stage_atoms,
    optimize,
)
from .threads import resolve_cores

from ash import *

HARTREE_TO_KCAL = 627.509474


def run(args) -> Dict[str, Any]:
    section = getattr(args, "section", None) or "scan"
    preprocess = load_section(args.config, "preprocess")
    md = load_section(args.config, "md")
    qmmm_inputs = load_section(args.config, "qmmm").get("inputs") or {}
    inputs = stage_inputs(args.config, "scan", section)

    if "atoms" not in inputs or "end" not in inputs:
        raise KeyError(f"{section}.inputs needs `atoms: [i, j]` and `end` (and usually `points`)")
    npoints = int(inputs.get("points", 10))
    if npoints < 2:
        raise ValueError("points must be at least 2")

    pdb_file = getattr(args, "pdb", None) or inputs.get("pdb_file") or \
        (load_section(args.config, "qmmm").get("results") or {}).get("optimized_pdb")
    if pdb_file is None:
        raise ValueError("No starting structure: give --pdb (e.g. the reactant's QM_MM_optimized.pdb)")

    charge, mult = qmmm_inputs.get("charge"), qmmm_inputs.get("mult")
    cores = resolve_cores(inputs.get("cores", qmmm_inputs.get("cores")))
    forcefields, implicit = mm_settings(qmmm_inputs, md)

    print("[STEP] Building MM system...")
    mm = make_mm(forcefields, pdb_file, implicit, cores)
    fragment = Fragment(pdbfile=pdb_file)
    inputs = resolve_stage_atoms(inputs, mm)
    qmmm_inputs = resolve_stage_atoms(qmmm_inputs, mm)
    i, j = (int(a) for a in inputs["atoms"])

    print("[STEP] Building QM region...")
    qm_atoms, boundary_excluded = resolve_qm_region({**qmmm_inputs, **{k: inputs[k] for k in ('qm_residues', 'extra_qm_residues') if k in inputs}}, preprocess, pdb_file)
    qmmm_theory = make_qmmm(qmmm_inputs, mm, fragment, qm_atoms, boundary_excluded, cores)

    act_inputs = {**qmmm_inputs, **{k: inputs[k] for k in ("actradius", "actregion_origin") if k in inputs}}
    actatoms = resolve_actatoms(getattr(args, "actatoms", None), act_inputs, mm, fragment,
                                qm_atoms, boundary_excluded)
    extra = inputs.get("constraints") or []
    check_in_active([[i, j]] + list(extra), actatoms, "scan")
    monitor = [[int(a), int(b)] for a, b in (inputs.get("monitor") or [])]

    d0 = distance(fragment.coords, i, j)
    start = float(inputs.get("start", d0))
    end = float(inputs["end"])
    targets = [float(x) for x in np.linspace(start, end, npoints)]
    print(f"[INFO] Scan {i}-{j}: current {d0:.3f} Å; {npoints} points {start:.3f} -> {end:.3f} Å "
          f"(step {targets[1] - targets[0]:+.3f} Å); active region {len(actatoms)} atoms")

    # Held at their starting values throughout (waters rigid, extra pairs as given)
    held = bond_constraints(water_constraints(mm, actatoms), fragment.coords) + \
        bond_constraints(extra, fragment.coords)

    result: Dict[str, Any] = {
        "pdb_file": Path(pdb_file).resolve().as_posix(),
        "atoms": [i, j],
        "constraints": bond_constraints(extra, fragment.coords),
        "points": [],
    }
    csvfile = Path("scan_energies.csv")
    with open(csvfile, "w", newline="") as f:
        csv.writer(f).writerow(["point", "target_A", "distance_A", "energy_Eh", "rel_kcal_mol", "converged"]
                               + [f"d_{a}_{b}" for a, b in monitor] + ["pdb"])

    e_first = None
    for n, d in enumerate(targets):
        print(f"\n[STEP] Scan point {n}/{npoints - 1}: {i}-{j} = {d:.3f} Å")
        clean_geometric_files()
        # Fresh lists each time: ASH rewrites the indices to active-region indices in place
        constraints = {'bond': [list(c) for c in held] + [[i, j, d]]}
        # A point that runs out of steps keeps its last geometry (converged: false) and the scan goes on
        energy, converged = optimize(fragment, qmmm_theory, actatoms, {"charge": charge, "mult": mult},
                                     constraints, True, inputs.get("opt_maxiter", 200))
        energy = float(energy)
        if e_first is None:
            e_first = energy
        pdb = write_pdb(mm, fragment.coords, f"scan_{n:02d}.pdb")
        point = {
            "target": round(d, 4),
            "distance": round(distance(fragment.coords, i, j), 4),
            "energy": energy,
            "rel_kcal_mol": round((energy - e_first) * HARTREE_TO_KCAL, 3),
            "converged": converged,
            "monitor": {f"{a}-{b}": round(distance(fragment.coords, a, b), 3) for a, b in monitor},
            "pdb": pdb,
        }
        result["points"].append(point)
        with open(csvfile, "a", newline="") as f:
            csv.writer(f).writerow([n, point["target"], point["distance"], f"{energy:.8f}",
                                    point["rel_kcal_mol"], converged]
                                   + list(point["monitor"].values()) + [pdb])
        print(f"[INFO] Point {n}: {point['distance']:.3f} Å, {point['rel_kcal_mol']:+.2f} kcal/mol"
              f"{'' if converged else ' (NOT converged)'}; {point['monitor']}")
        # Saved after every point so a stopped job keeps what it finished
        result["energies_csv"] = csvfile.resolve().as_posix()
        save_result(result, filename=args.config, section=section)

    rel = [p["rel_kcal_mol"] for p in result["points"]]
    top = int(np.argmax(rel))
    print(f"\n[INFO] Scan maximum at point {top} ({targets[top]:.3f} Å): {rel[top]:.2f} kcal/mol "
          f"relative to point 0; last point {rel[-1]:+.2f} kcal/mol")
    result["max_point"] = top
    save_result(result, filename=args.config, section=section)
    return result
