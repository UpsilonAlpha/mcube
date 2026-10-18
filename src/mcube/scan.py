#!/usr/bin/env python3
"""
Relaxed QM/MM distance scan, 1D or 2D (xTB/MM, same model as the qmmm stage).

1D: the distance between two atoms is stepped from `start` to `end`; at each
point the active region is optimised with that distance fixed, starting from
the previous point's geometry. Used to build intermediates and products inside
the reactant's box (identical atoms for NEB) and to get a first energy profile.

2D: a second distance (`atoms2`, `start2`, `end2`, `points2`) makes a grid,
e.g. nucleophile–carbon distance x proton–acceptor distance, so that a proton
transfer coupled to bond formation is mapped instead of assumed. Rows run over
the second coordinate; within a row the first coordinate is swept in
alternating directions ("snake" order), so every point starts from an adjacent
grid point.

Config (the scan is defined in the section named by --section, default `scan`;
anything not given there is taken from `scan.inputs`):

    scan_attack:
      inputs:
        atoms: ["C:2:O01", "C:1:C4x"]   # i, j: "chain:resid:atom" or 0-based indices
        start: 3.0            # Å; default: the current distance
        end: 1.45             # Å
        points: 16            # including both ends
        atoms2: ["D:1:H1", "A:92:OD1"]  # optional second coordinate (2D grid)
        start2: 1.8
        end2: 1.0
        points2: 6
        constraints: []       # extra bonds held during the scan, [[i, j]] or [[i, j, d]]
        monitor: []           # distances [[i, j], ...] reported at every point
        opt_maxiter: 200
        resume: true          # continue from scan_energies.csv if it exists

Charge, multiplicity, force fields, xTB settings, cores and the active region
come from qmmm.inputs (`extra_qm_residues`/`qm_residues` may be set per
section). Outputs (in the working directory): one PDB per point
(scan_XX.pdb, or scan_KK_XX.pdb for a grid) and scan_energies.csv; results go
to the section after every point.
"""

from __future__ import annotations

import csv
import gc
from pathlib import Path
from typing import Dict, Any, List, Tuple

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


def grid_order(n1: int, n2: int) -> List[Tuple[int, int]]:
    """(row, column) visiting order: rows over coordinate 2, columns snake back and forth."""
    order = []
    for k in range(n2):
        cols = range(n1) if k % 2 == 0 else reversed(range(n1))
        order.extend((k, n) for n in cols)
    return order


def targets_for(inputs: Dict[str, Any], suffix: str, current: float) -> List[float]:
    npoints = int(inputs.get(f"points{suffix}", 10))
    if npoints < 2:
        raise ValueError(f"points{suffix} must be at least 2")
    start = float(inputs.get(f"start{suffix}", current))
    end = float(inputs[f"end{suffix}"])
    return [float(x) for x in np.linspace(start, end, npoints)]


def read_done(csvfile: Path) -> List[Dict[str, str]]:
    if not csvfile.exists():
        return []
    with open(csvfile, newline="") as f:
        return list(csv.DictReader(f))


def run(args) -> Dict[str, Any]:
    section = getattr(args, "section", None) or "scan"
    preprocess = load_section(args.config, "preprocess")
    md = load_section(args.config, "md")
    qmmm_inputs = load_section(args.config, "qmmm").get("inputs") or {}
    inputs = stage_inputs(args.config, "scan", section)

    if "atoms" not in inputs or "end" not in inputs:
        raise KeyError(f"{section}.inputs needs `atoms: [i, j]` and `end` (and usually `points`)")
    two_d = inputs.get("atoms2") is not None
    if two_d and "end2" not in inputs:
        raise KeyError(f"{section}.inputs: a 2D scan needs `atoms2` and `end2` (and usually `points2`)")

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
    i2, j2 = (int(a) for a in inputs["atoms2"]) if two_d else (None, None)

    print("[STEP] Building QM region...")
    qm_atoms, boundary_excluded = resolve_qm_region(
        {**qmmm_inputs, **{k: inputs[k] for k in ('qm_residues', 'extra_qm_residues') if k in inputs}},
        preprocess, pdb_file)
    qmmm_theory = make_qmmm(qmmm_inputs, mm, fragment, qm_atoms, boundary_excluded, cores)

    act_inputs = {**qmmm_inputs, **{k: inputs[k] for k in ("actradius", "actregion_origin") if k in inputs}}
    actatoms = resolve_actatoms(getattr(args, "actatoms", None), act_inputs, mm, fragment,
                                qm_atoms, boundary_excluded)
    extra = inputs.get("constraints") or []
    check_in_active([[i, j]] + ([[i2, j2]] if two_d else []) + list(extra), actatoms, "scan")
    monitor = [[int(a), int(b)] for a, b in (inputs.get("monitor") or [])]

    targets1 = targets_for(inputs, "", distance(fragment.coords, i, j))
    targets2 = targets_for(inputs, "2", distance(fragment.coords, i2, j2)) if two_d else [None]
    order = grid_order(len(targets1), len(targets2))
    print(f"[INFO] Scan {i}-{j}: {targets1[0]:.3f} -> {targets1[-1]:.3f} Å in {len(targets1)} points"
          + (f"; x {i2}-{j2}: {targets2[0]:.3f} -> {targets2[-1]:.3f} Å in {len(targets2)} points "
             f"({len(order)} grid points)" if two_d else "")
          + f"; active region {len(actatoms)} atoms")

    # Held at their starting values throughout (MM waters rigid, extra pairs as given)
    held = bond_constraints(water_constraints(mm, actatoms, qm_atoms), fragment.coords) + \
        bond_constraints(extra, fragment.coords)

    result: Dict[str, Any] = {
        "pdb_file": Path(pdb_file).resolve().as_posix(),
        "atoms": [i, j],
        "constraints": bond_constraints(extra, fragment.coords),
        "points": [],
    }
    if two_d:
        result["atoms2"] = [i2, j2]

    csvfile = Path("scan_energies.csv")
    if two_d:
        header = ["point", "i1", "i2", "target1_A", "distance1_A", "target2_A", "distance2_A"]
    else:
        header = ["point", "target_A", "distance_A"]
    header += ["energy_Eh", "rel_kcal_mol", "converged"] + [f"d_{a}_{b}" for a, b in monitor] + ["pdb"]

    # Resume: skip points already in the CSV and continue from the last structure written
    done = read_done(csvfile) if inputs.get("resume", True) else []
    e_first = None
    if done:
        e_first = float(done[0]["energy_Eh"])
        last_pdb = done[-1]["pdb"]
        fragment.coords = np.array(Fragment(pdbfile=last_pdb).coords)
        print(f"[INFO] Resuming: {len(done)} of {len(order)} points already done; continuing from {last_pdb}")
        result["points"] = [{k: v for k, v in r.items()} for r in done]
    else:
        with open(csvfile, "w", newline="") as f:
            csv.writer(f).writerow(header)

    for n_done, (k, n) in enumerate(order):
        if n_done < len(done):
            continue
        d1, d2 = targets1[n], targets2[k]
        label = f"{i}-{j} = {d1:.3f} Å" + (f", {i2}-{j2} = {d2:.3f} Å" if two_d else "")
        print(f"\n[STEP] Scan point {n_done}/{len(order) - 1}: {label}")
        clean_geometric_files()
        # Fresh lists each time: ASH rewrites the indices to active-region indices in place
        bonds = [list(c) for c in held] + [[i, j, d1]] + ([[i2, j2, d2]] if two_d else [])
        # A point that runs out of steps keeps its last geometry (converged: false) and the scan goes on
        energy, converged = optimize(fragment, qmmm_theory, actatoms, {"charge": charge, "mult": mult},
                                     {'bond': bonds}, True, inputs.get("opt_maxiter", 200))
        energy = float(energy)
        if e_first is None:
            e_first = energy
        name = f"scan_{k:02d}_{n:02d}.pdb" if two_d else f"scan_{n:02d}.pdb"
        pdb = write_pdb(mm, fragment.coords, name)
        point = {
            "target": round(d1, 4),
            "distance": round(distance(fragment.coords, i, j), 4),
            "energy": energy,
            "rel_kcal_mol": round((energy - e_first) * HARTREE_TO_KCAL, 3),
            "converged": converged,
            "monitor": {f"{a}-{b}": round(distance(fragment.coords, a, b), 3) for a, b in monitor},
            "pdb": pdb,
        }
        if two_d:
            point.update({"i1": n, "i2": k, "target2": round(d2, 4),
                          "distance2": round(distance(fragment.coords, i2, j2), 4)})
            row = [n_done, n, k, point["target"], point["distance"], point["target2"], point["distance2"]]
        else:
            row = [n, point["target"], point["distance"]]
        row += [f"{energy:.8f}", point["rel_kcal_mol"], converged] + list(point["monitor"].values()) + [pdb]
        result["points"].append(point)
        with open(csvfile, "a", newline="") as f:
            csv.writer(f).writerow(row)
        print(f"[INFO] Point {n_done}: {label}: {point['rel_kcal_mol']:+.2f} kcal/mol"
              f"{'' if converged else ' (NOT converged)'}; {point['monitor']}")
        # Saved after every point so a stopped job keeps what it finished
        result["energies_csv"] = csvfile.resolve().as_posix()
        save_result(result, filename=args.config, section=section)
        gc.collect()

    rel = [float(p["rel_kcal_mol"]) for p in result["points"]]
    top = int(np.argmax(rel))
    print(f"\n[INFO] Highest point {top}: {rel[top]:.2f} kcal/mol relative to the first point; "
          f"last point {rel[-1]:+.2f} kcal/mol")
    result["max_point"] = top
    save_result(result, filename=args.config, section=section)
    return result
