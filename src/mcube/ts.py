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
        hessian_atoms: ["B:1:ZN", "C:1:C4x", "C:1:O2x", "C:2:O01"]  # default: the QM atoms
        monitor: [["C:2:O01", "C:1:C4x"]]   # distances to report at TS and endpoints
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
        rigid_water: false           # off by default: a held water distance becomes a
                                     # geomeTRIC `$set` convergence criterion that rigid
                                     # TIP3P cannot satisfy (D137). QM waters are never
                                     # constrained either way.

Charge, multiplicity, force fields, xTB settings, cores and the active region
come from qmmm.inputs; set `actatoms_file` to the reactant's active_atoms.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Dict, Any, List

import numpy as np

from .preprocess import resolve_cores, save_result, save_input, load_section
from .qmmm import (
    apply_stage_metal_parameters,
    resolve_qm_region,
    stage_inputs,
    mm_settings,
    make_mm,
    make_qmmm,
    ACT_KEYS,
    actatoms_record,
    resolve_actatoms,
    resolve_start_structure,
    distance,
    check_in_active,
    rigid_water_pairs,
    optimizer_constraints,
    clean_geometric_files,
    write_pdb,
    resolve_stage_atoms,
    qm_level_inputs,
    qm_level_label,
)

from ash import *

HARTREE_TO_KCAL = 627.509474


def monitor_distances(coords, pairs) -> Dict[str, float]:
    return {f"{int(i)}-{int(j)}": round(distance(coords, int(i), int(j)), 3) for i, j in (p[:2] for p in pairs)}


def check_ts_mode(amplitudes, hessatoms, monitor, floor: float = 0.3):
    """Whether the imaginary mode is carried by the atoms the step is about.

    One imaginary frequency makes a structure a first-order saddle. It does not
    make it the transition state of *this* reaction: a saddle whose mode is a
    water libration or a methyl rotation passes the count and means nothing.
    `monitor`'s pairs are the bonds the step makes or breaks, so their atoms
    should carry a substantial share of the displacement.

    Returns (ok, share, names_of_reacting_atoms_found). `ok` is None when there
    is nothing to check against - no monitor pairs, or none of their atoms in
    the Hessian - which is not a pass.
    """
    index = {int(a): i for i, a in enumerate(hessatoms)}
    wanted = sorted({int(a) for pair in (monitor or []) for a in pair[:2] if int(a) in index})
    if not wanted or amplitudes.max() <= 0:
        return None, None, wanted
    share = float(max(amplitudes[index[a]] for a in wanted) / amplitudes.max())
    return share >= floor, share, wanted


def run(args) -> Dict[str, Any]:
    section = getattr(args, "section", None) or "ts"
    preprocess = load_section(args.config, "preprocess")
    md = load_section(args.config, "md")
    qmmm_inputs = load_section(args.config, "qmmm").get("inputs") or {}
    inputs = stage_inputs(args.config, "ts", section)
    result["qm_level"] = qm_level_label(qmmm_inputs)
    qmmm_inputs = qm_level_inputs(qmmm_inputs, inputs)  # qm_engine/dft_* may be set per section

    pdb_file, pdb_source = resolve_start_structure(
        args.config, section, inputs,
        [((load_section(args.config, "neb").get("results") or {}).get("saddlepoint_pdb"),
          "neb.results.saddlepoint_pdb")], stage="ts", what="TS guess")
    print(f"[INFO] TS guess: {pdb_file} (from {pdb_source})")
    if not inputs.get("pdb_file"):
        save_input(args.config, section, "pdb_file", pdb_file)

    charge, mult = qmmm_inputs.get("charge"), qmmm_inputs.get("mult")
    cores = resolve_cores(inputs.get("cores", qmmm_inputs.get("cores")))
    workers = max(1, min(int(inputs.get("freq_workers", 1)), cores))
    if workers > 1:
        print(f"[INFO] {workers} NumFreq workers run in parallel: run this stage "
              "with OMP_PROC_BIND/OMP_PLACES unset (the xtb-parallel recipe)")
    forcefields, implicit = mm_settings(qmmm_inputs, md, args.config)

    print("[STEP] Building MM system...")
    mm = make_mm(forcefields, pdb_file, implicit, cores)
    fragment = Fragment(pdbfile=pdb_file, charge=charge, mult=mult)
    inputs = resolve_stage_atoms(inputs, mm)
    qmmm_inputs = resolve_stage_atoms(qmmm_inputs, mm)

    print("[STEP] Building QM region...")
    qm_atoms, boundary_excluded = resolve_qm_region({**qmmm_inputs, **{k: inputs[k] for k in ('qm_residues', 'extra_qm_residues') if k in inputs}}, preprocess, pdb_file)
    # The same metal parameters as md, before the QMMMTheory is built
    apply_stage_metal_parameters(qmmm_inputs, preprocess, md, mm, forcefields, qm_atoms)
    # Parallel NumFreq runs `workers` gradients at once, so give each fewer xTB threads
    qmmm_theory = make_qmmm(qmmm_inputs, mm, fragment, qm_atoms, boundary_excluded,
                            cores // workers if workers > 1 else cores,
                            )   # not 'library': no QM/MM point-charge gradient there

    act_inputs = {**qmmm_inputs, **{k: inputs[k] for k in ACT_KEYS if k in inputs}}
    actatoms = resolve_actatoms(args.config, act_inputs, mm, fragment,
                                qm_atoms, boundary_excluded, pdb_file)
    result.update(actatoms_record(mm, pdb_file, actatoms))
    hessatoms: List[int] = [int(a) for a in (inputs.get("hessian_atoms") or qm_atoms)]
    monitor = inputs.get("monitor") or []
    check_in_active([[a, a] for a in hessatoms], actatoms, "hessian_atoms")
    print(f"[INFO] Active region {len(actatoms)} atoms; Hessian over {len(hessatoms)} atoms "
          f"(~{3 * len(hessatoms) * int(inputs.get('freq_npoint', 2))} gradients per frequency job)")

    water = rigid_water_pairs(mm, actatoms, inputs, qm_atoms)

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
        ok, share, wanted = check_ts_mode(amp, hessatoms, inputs.get("monitor"))
        result["ts_mode_is_reaction_coordinate"] = ok
        if share is not None:
            result["ts_mode_reacting_share"] = round(share, 3)
        if ok is None:
            print("[WARNING] no `monitor` pairs among the Hessian atoms, so whether this "
                  "mode is the reaction coordinate was not checked - only that the "
                  "structure is a saddle point")
        elif ok:
            print(f"[INFO] The mode is carried by the reacting atoms {wanted} "
                  f"({share:.0%} of the largest displacement)")
        else:
            print(f"[WARNING] the reacting atoms {wanted} carry only {share:.0%} of the "
                  "largest displacement in this mode: it is a saddle point, but likely "
                  "not the transition state of this step. Compare `ts_mode_atoms` against "
                  "the bonds the step should be making or breaking")
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
