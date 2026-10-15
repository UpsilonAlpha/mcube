#!/usr/bin/env python3
"""
Transition-state refinement, frequency check and IRC (xTB/MM, same model as qmmm).

Starting from a TS guess (usually the NEB climbing image, saddlepoint.pdb):
1. TS optimisation (geomeTRIC, TSOpt) with an exact partial Hessian over
   `hessian_atoms` (the rest of the active region from a model Hessian), as
   ASH's NEBTS does.
2. Numerical frequencies over `hessian_atoms`: a TS has exactly one imaginary
   frequency, and its mode should be the reaction coordinate.
3. IRC in both directions:
   - `irc_method: mode` (default): displace along the imaginary mode by
     ±`irc_step` Å (largest atom displacement) and minimise each side. A cheap
     quasi-IRC that shows which minima the TS connects.
   - `irc_method: sella`: Sella's mass-weighted IRC on the active region
     (needs `sella`; expensive for a large active region).

Config (section named by --section, default `ts`; anything not given there is
taken from `ts.inputs`):

    ts_attack:
      inputs:
        hessian_atoms: [4203, 4208, 4209, 4210, 4255, 4256]  # default: the QM atoms
        monitor: [[4255, 4208], [4208, 4210]]   # distances to report at TS and endpoints
        tsopt: true
        tsopt_maxiter: 100
        convergence_setting: null    # geomeTRIC preset, e.g. GAU
        modelhessian: unit           # rest of the Hessian: unit, Almloef, Lindh, Schlegel, xtb
        hessian_npoint: 1            # 1: one-sided (3N gradients), 2: central (6N)
        freq: true
        freq_npoint: 2
        freq_workers: 1              # >1 runs NumFreq displacements in parallel
        imag_threshold: 50           # |cm-1| below which imaginary modes are reported as noise
        irc: true
        irc_method: mode
        irc_step: 0.15
        opt_maxiter: 300
        rigid_water: true            # hold active-region waters rigid, as in qmmm minimisations

Charge, multiplicity, force fields, xTB settings, cores and the active region
come from qmmm.inputs; use --actatoms with the reactant's active_atoms.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Dict, Any, List

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
    optimizer_constraints,
    clean_geometric_files,
    write_pdb,
)

from ash import *

HARTREE_TO_KCAL = 627.509474


def monitor_distances(coords, pairs) -> Dict[str, float]:
    return {f"{int(i)}-{int(j)}": round(distance(coords, int(i), int(j)), 3) for i, j in (p[:2] for p in pairs)}


def run(args) -> Dict[str, Any]:
    section = getattr(args, "section", None) or "ts"
    preprocess = load_section(args.config, "preprocess")
    md = load_section(args.config, "md")
    qmmm_inputs = load_section(args.config, "qmmm").get("inputs") or {}
    inputs = stage_inputs(args.config, "ts", section)

    pdb_file = getattr(args, "pdb", None) or inputs.get("pdb_file") or \
        (load_section(args.config, "neb").get("results") or {}).get("saddlepoint_pdb")
    if pdb_file is None:
        raise ValueError("No TS guess: give --pdb (e.g. the NEB's saddlepoint.pdb)")

    charge, mult = qmmm_inputs.get("charge"), qmmm_inputs.get("mult")
    cores = int(inputs.get("cores", qmmm_inputs.get("cores", 1)))
    workers = max(1, min(int(inputs.get("freq_workers", 1)), cores))
    forcefields, implicit = mm_settings(qmmm_inputs, md)

    print("[STEP] Building MM system...")
    mm = make_mm(forcefields, pdb_file, implicit, cores)
    fragment = Fragment(pdbfile=pdb_file, charge=charge, mult=mult)

    print("[STEP] Building QM region...")
    qm_atoms, boundary_excluded = resolve_qm_region(qmmm_inputs, preprocess, pdb_file)
    # Parallel NumFreq runs `workers` gradients at once, so give each fewer xTB threads
    qmmm_theory = make_qmmm(qmmm_inputs, mm, fragment, qm_atoms, boundary_excluded,
                            cores // workers if workers > 1 else cores,
                            runmode='library' if workers > 1 else None)

    act_inputs = {**qmmm_inputs, **{k: inputs[k] for k in ("actradius", "actregion_origin") if k in inputs}}
    actatoms = resolve_actatoms(getattr(args, "actatoms", None), act_inputs, mm, fragment,
                                qm_atoms, boundary_excluded)
    hessatoms: List[int] = [int(a) for a in (inputs.get("hessian_atoms") or qm_atoms)]
    monitor = inputs.get("monitor") or []
    check_in_active([[a, a] for a in hessatoms], actatoms, "hessian_atoms")
    print(f"[INFO] Active region {len(actatoms)} atoms; Hessian over {len(hessatoms)} atoms "
          f"(~{3 * len(hessatoms) * int(inputs.get('freq_npoint', 2))} gradients per frequency job)")

    water = water_constraints(mm, actatoms) if inputs.get("rigid_water", True) else []

    def constraints():
        # Fresh lists each call: ASH rewrites the indices to active-region indices in place
        return optimizer_constraints(water, [], fragment.coords)[0] if water else None

    def numfreq(frag, npoint):
        qmmm_theory.set_numcores(cores // workers if workers > 1 else cores)
        return NumFreq(theory=qmmm_theory, fragment=frag, hessatoms=hessatoms, npoint=npoint,
                       charge=charge, mult=mult, numcores=workers,
                       runmode='parallel' if workers > 1 else 'serial', printlevel=0)

    result: Dict[str, Any] = {
        "pdb_file": Path(pdb_file).resolve().as_posix(),
        "hessian_atoms": hessatoms,
        "guess_distances": monitor_distances(fragment.coords, monitor),
    }

    # ---------- 1. TS optimisation ----------
    ts_energy = None
    if inputs.get("tsopt", True):
        print("[STEP] Partial Hessian for the TS optimisation...")
        freq0 = numfreq(fragment, int(inputs.get("hessian_npoint", 1)))
        guess_freqs = [float(f) for f in freq0.frequencies]
        result["guess_imaginary"] = [round(f, 1) for f in guess_freqs if f < 0]
        print(f"[INFO] TS guess imaginary frequencies: {result['guess_imaginary']}")
        full_hessian = approximate_full_Hessian_from_smaller(
            fragment, freq0.hessian, hessatoms, large_atomindices=actatoms,
            restHessian=inputs.get("modelhessian", "unit"), charge=charge, mult=mult)
        write_hessian(full_hessian, hessfile="Hessian_partial_TS")

        print("[STEP] TS optimisation...")
        qmmm_theory.set_numcores(cores)
        clean_geometric_files()
        opt = Optimizer(fragment=fragment, theory=qmmm_theory, TSOpt=True,
                        hessian="file:" + str(Path("Hessian_partial_TS").resolve()),
                        ActiveRegion=True, actatoms=actatoms, maxiter=inputs.get("tsopt_maxiter", 100),
                        convergence_setting=inputs.get("convergence_setting"),
                        constraints=constraints(), charge=charge, mult=mult)
        ts_energy = float(opt.energy)
        result["ts_pdb"] = write_pdb(mm, fragment.coords, "ts_optimized.pdb")
        result["tsopt_traj"] = Path("geometric_OPTtraj_Full.xyz").resolve().as_posix()
    if ts_energy is None:
        qmmm_theory.set_numcores(cores)
        ts_energy = float(Singlepoint(theory=qmmm_theory, fragment=fragment, charge=charge, mult=mult).energy)
    result["ts_energy"] = ts_energy
    result["ts_distances"] = monitor_distances(fragment.coords, monitor)
    save_result(result, filename=args.config, section=section)

    # ---------- 2. Frequencies ----------
    mode = None
    if inputs.get("freq", True):
        print("[STEP] Frequencies at the TS...")
        freq = numfreq(fragment, int(inputs.get("freq_npoint", 2)))
        freqs = [float(f) for f in freq.frequencies]
        threshold = float(inputs.get("imag_threshold", 50))
        imag = [f for f in freqs if f < -threshold]
        result["imaginary_frequencies"] = [round(f, 1) for f in imag]
        result["small_imaginary"] = [round(f, 1) for f in freqs if -threshold <= f < 0]
        result["lowest_frequencies"] = [round(f, 1) for f in sorted(freqs)[:8]]
        k = int(np.argmin(freqs))
        mode = np.asarray(freq.normal_modes)[k].reshape(-1, 3)
        amp = np.linalg.norm(mode, axis=1)
        top = np.argsort(amp)[::-1][:8]
        result["ts_mode_atoms"] = {int(hessatoms[a]): round(float(amp[a] / amp.max()), 3) for a in top}
        if len(imag) == 1:
            print(f"[INFO] One imaginary frequency, {imag[0]:.1f} cm-1: a first-order saddle point")
        else:
            print(f"[WARNING] {len(imag)} imaginary frequencies below -{threshold:g} cm-1: {imag} "
                  "(a TS has exactly one)")
        print(f"[INFO] Atoms moving most in the lowest mode ({freqs[k]:.1f} cm-1): {result['ts_mode_atoms']}")
        save_result(result, filename=args.config, section=section)

    # ---------- 3. IRC ----------
    if inputs.get("irc", True):
        method = str(inputs.get("irc_method", "mode")).lower()
        qmmm_theory.set_numcores(cores)
        if method == "sella":
            print("[STEP] Sella IRC (forward and reverse)...")
            SellaIRC(theory=qmmm_theory, fragment=copy.deepcopy(fragment), charge=charge, mult=mult,
                     actatoms=actatoms, maxiter=inputs.get("opt_maxiter", 300),
                     IRC_dx=inputs.get("irc_dx", 0.1))
            result["irc_traj"] = Path("sella_optim_full.xyz" if actatoms else "sella_IRC.xyz").resolve().as_posix()
        elif method == "mode":
            if mode is None:
                raise ValueError("irc_method: mode needs the TS mode: set freq: true")
            step = float(inputs.get("irc_step", 0.15))
            disp = mode / np.linalg.norm(mode, axis=1).max() * step
            result["irc"] = {}
            for label, sign in (("forward", 1.0), ("reverse", -1.0)):
                print(f"[STEP] Quasi-IRC {label}: {sign * step:+.2f} Å along the TS mode, then minimise...")
                frag = copy.deepcopy(fragment)
                coords = np.array(frag.coords, dtype=float)
                coords[hessatoms] += sign * disp
                frag.coords = coords
                clean_geometric_files()
                opt = Optimizer(fragment=frag, theory=qmmm_theory, ActiveRegion=True, actatoms=actatoms,
                                maxiter=inputs.get("opt_maxiter", 300), constraints=constraints(),
                                charge=charge, mult=mult)
                e = float(opt.energy)
                result["irc"][label] = {
                    "pdb": write_pdb(mm, frag.coords, f"irc_{label}.pdb"),
                    "energy": e,
                    "barrier_kcal_mol": round((ts_energy - e) * HARTREE_TO_KCAL, 2),
                    "distances": monitor_distances(frag.coords, monitor),
                }
                if Path("geometric_OPTtraj_Full.xyz").exists():
                    Path("geometric_OPTtraj_Full.xyz").replace(f"irc_{label}_opttraj.xyz")
                print(f"[INFO] IRC {label}: E(TS) - E = {result['irc'][label]['barrier_kcal_mol']} kcal/mol; "
                      f"distances {result['irc'][label]['distances']}")
                save_result(result, filename=args.config, section=section)
        else:
            raise ValueError(f"irc_method must be 'mode' or 'sella', not {method!r}")

    save_result(result, filename=args.config, section=section)
    return result
