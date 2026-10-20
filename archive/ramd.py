#!/usr/bin/env python3
"""
τRAMD: relative residence times from random acceleration MD.

Protocol (Kokh et al. 2018, 2020): several replicas, each started from its own
unbiased equilibration (new velocities, optionally a different snapshot); from
each, many RAMD trajectories with new velocities and random force directions,
each run until the binder leaves (COM distance > max_distance). τ per replica
is the bootstrapped time at which half of its trajectories have dissociated;
τ is the mean over replicas (mcube.ramd_analysis).

The classical system is exactly the md stage's (force fields, HBond
constraints, rigid water, metal-site constraints), so metal sites stay held as
in every other MM stage. md restraints are NOT applied unless given in
ramd.inputs.restraints (a restraint tying the binder to the receptor would
prevent unbinding, and is refused).

Work is restartable and splittable: each replica lives in replica_RR/, stop
times are appended to ramd_times.csv (one row per trajectory, file-locked), and
finished trajectories are skipped. `--replicas 2` runs only replica 2 (one SLURM
array task per replica); `--analyse` only analyses ramd_times.csv.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Any, List

import numpy as np
import openmm
from openmm import unit
from openmm.app import PDBFile

from .preprocess import save_result, load_section, resolve_atoms, _locked
from .qmmm import stage_inputs
from .md import build_md_theory, md_constraints
from .ramd_engine import RAMD, com_and_radius
from .ramd_analysis import tau_ramd, format_report
from .threads import resolve_cores

WATER_IONS = {"HOH", "WAT", "SOL", "TIP3", "NA", "CL", "K", "Na+", "Cl-", "K+", "NA+", "CL-"}
CSV_FIELDS = ["replica", "trajectory", "seed", "dissociated", "time_ps", "steps", "direction_changes",
              "start_distance_A", "final_distance_A", "wall_s"]


def residue_atoms(topology, spec) -> List[int]:
    """
    Atoms of the residues in `spec`: {chain: [resid, ...]} or {chain: "all"} (the
    preprocess residue format), or a list of atoms (indices / "chain:resid:atom").
    """
    if isinstance(spec, list):
        return [int(i) for i in resolve_atoms(spec, topology)]
    if not isinstance(spec, dict):
        raise ValueError(f"binder/receptor must be {{chain: [resids] or 'all'}} or a list of atoms, not {spec!r}")
    want = {str(ch): (None if str(v).lower() == "all" else {str(r) for r in v}) for ch, v in spec.items()}
    idx = [a.index for a in topology.atoms()
           if a.residue.chain.id in want and (want[a.residue.chain.id] is None
                                              or a.residue.id.strip() in want[a.residue.chain.id])]
    if not idx:
        raise ValueError(f"No atoms match {spec}")
    return idx


def default_receptor(topology, binder) -> List[int]:
    """Everything that is not the binder, water or a monatomic ion."""
    b = set(binder)
    return [a.index for a in topology.atoms() if a.index not in b and a.residue.name not in WATER_IONS]


def check_tethers(pairs, binder, label):
    """A constraint/restraint between a binder atom and a non-binder atom would hold the binder in place."""
    b = set(binder)
    bad = [p for p in pairs or [] if (int(p[0]) in b) != (int(p[1]) in b)]
    if bad:
        raise ValueError(f"{label} {bad} tie the binder to the rest of the system and would prevent "
                         "unbinding; drop them in ramd.inputs")


def read_times(csvfile: Path):
    if not csvfile.exists():
        return []
    with open(csvfile) as f:
        return list(csv.DictReader(f))


def append_row(csvfile: Path, row: Dict[str, Any]):
    with _locked(csvfile):
        new = not csvfile.exists()
        with open(csvfile, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
            if new:
                w.writeheader()
            w.writerow({k: row[k] for k in CSV_FIELDS})


def analyse(rows, inputs) -> Dict[str, Any]:
    times, cens = {}, {}
    for r in rows:
        times.setdefault(r["replica"], []).append(float(r["time_ps"]))
        cens.setdefault(r["replica"], []).append(str(r["dissociated"]) not in ("True", "true", "1"))
    return tau_ramd(times, cens, rounds=int(inputs.get("bootstrap_rounds", 5000)),
                    fraction=float(inputs.get("bootstrap_fraction", 0.8)),
                    seed=inputs.get("seed", 0), min_trajectories=int(inputs.get("min_trajectories", 10)))


def run(args) -> Dict[str, Any]:
    section = getattr(args, "section", None) or "ramd"
    inputs = stage_inputs(args.config, "ramd", section)
    preprocess = load_section(args.config, "preprocess")
    md = load_section(args.config, "md")
    md_inputs = md.get("inputs") or {}
    md_results = md.get("results") or {}
    csvfile = Path(inputs.get("times_csv", "ramd_times.csv")).resolve()

    if getattr(args, "analyse", False):
        res = analyse(read_times(csvfile), inputs)
        print(format_report(res))
        res["times_csv"] = csvfile.as_posix()
        save_result(res, filename=args.config, section=section)
        return res

    if inputs.get("binder") is None:
        raise KeyError(f"{section}.inputs.binder is required, e.g. {{C: all}} or {{B: [1]}}")
    pdb_file = inputs.get("pdb_file") or md_results.get("lastframe")
    if pdb_file is None:
        raise ValueError("No start structure: set ramd.inputs.pdb_file or run md first")
    replica_pdbs = inputs.get("replica_pdbs")  # optional: a different snapshot per replica

    # ----- the md stage's classical system -----
    if "constraints" in inputs:
        constraints = resolve_atoms(inputs["constraints"], pdb_file)
    else:
        constraints = md_constraints(md_inputs, preprocess, pdb_file)
    restraints = resolve_atoms(inputs.get("restraints") or None, pdb_file)
    cores = resolve_cores(inputs.get("cores", md_inputs.get("cores")))
    mm, forcefields = build_md_theory(md_inputs, preprocess, pdb_file, constraints, restraints, cores=cores)
    topology = mm.topology
    implicit = bool(md_inputs.get("implicit"))

    binder = residue_atoms(topology, inputs["binder"])
    receptor = residue_atoms(topology, inputs["receptor"]) if inputs.get("receptor") else \
        default_receptor(topology, binder)
    check_tethers(constraints, binder, "constraints")
    check_tethers(restraints, binder, "restraints")

    # ----- geometry checks: exit distance and box size -----
    pdb = PDBFile(pdb_file)
    pos = np.array(pdb.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    masses = np.array([mm.system.getParticleMass(i).value_in_unit(unit.dalton)
                       for i in range(mm.system.getNumParticles())])
    com_b, r_b = com_and_radius(pos, masses, binder)
    com_r, r_r = com_and_radius(pos, masses, receptor)
    bound = 10 * float(np.linalg.norm(com_b - com_r))
    margin = float(inputs.get("exit_margin", 10.0))
    max_distance = float(inputs.get("max_distance") or (10 * (r_r + r_b) + margin))
    print(f"[INFO] Binder {len(binder)} atoms (radius {10 * r_b:.1f} Å), receptor {len(receptor)} atoms "
          f"(radius {10 * r_r:.1f} Å); bound COM distance {bound:.1f} Å; exit at {max_distance:.1f} Å")
    if max_distance <= bound:
        raise ValueError(f"max_distance {max_distance:.1f} Å is not beyond the bound COM distance {bound:.1f} Å")
    box = pdb.topology.getPeriodicBoxVectors()
    if not implicit:
        if box is None:
            raise ValueError(f"{pdb_file} has no CRYST1 box")
        lmin = 10 * min(box[0][0], box[1][1], box[2][2]).value_in_unit(unit.nanometer)
        cutoff = 10.0
        need = max_distance + 10 * (r_r + r_b) + cutoff
        if max_distance >= lmin / 2:
            raise ValueError(f"max_distance {max_distance:.1f} Å must be below half the box ({lmin / 2:.1f} Å): "
                             "the minimum-image COM distance cannot grow past it")
        if lmin < need:
            print(f"[WARNING] Box {lmin:.1f} Å < {need:.1f} Å (exit distance + receptor and binder radii + "
                  f"{cutoff:.0f} Å cutoff): near the exit the binder can meet the receptor's periodic image. "
                  "Increase preprocess.solvent_padding.")

    timestep = float(inputs.get("timestep", md_inputs.get("nvt_timestep") or 0.002))
    eval_steps = max(1, int(round(float(inputs.get("eval_interval", 0.1)) / timestep)))
    settings = {
        "force": float(inputs.get("force", 14.0)),
        "eval_steps": eval_steps,
        "rmin": float(inputs.get("rmin", 0.025)),
        "max_distance": max_distance,
        "timestep": timestep,
        "temperature": float(inputs.get("temp", md_inputs.get("temp") or 300.0)),
        "friction": float(inputs.get("friction", 1.0)),
        "pressure": None if implicit else inputs.get("pressure", 1.0),
        "periodic": not implicit,
        "remove_cm_motion": bool(inputs.get("remove_cm_motion", True)),
        "platform": inputs.get("platform"),
        "properties": inputs.get("properties"),
    }
    if settings["platform"] == "CPU" and not settings["properties"]:
        settings["properties"] = {"Threads": str(cores)}
    print(f"[INFO] RAMD settings: {settings}")
    engine = RAMD(mm.system, topology, binder, receptor, **settings)

    n_rep = int(inputs.get("replicas", 4))
    n_traj = int(inputs.get("trajectories", 15))
    replicas = getattr(args, "replicas", None) or list(range(n_rep))
    base_seed = int(inputs.get("seed", 2026))
    max_time = float(inputs.get("max_time", 10000.0))
    done = {(int(r["replica"]), int(r["trajectory"])) for r in read_times(csvfile)}

    result: Dict[str, Any] = {
        "binder_atoms": len(binder), "receptor_atoms": len(receptor),
        "bound_distance_A": round(bound, 2), "max_distance_A": round(max_distance, 2),
        "settings": {k: v for k, v in settings.items() if k != "properties"},
        "forcefields": forcefields, "times_csv": csvfile.as_posix(),
    }

    for r in replicas:
        rdir = Path(f"replica_{r:02d}")
        rdir.mkdir(exist_ok=True)
        state_file = rdir / "equilibrated.xml"
        if state_file.exists():
            state = openmm.XmlSerializer.deserialize(state_file.read_text())
            print(f"[INFO] Replica {r}: reusing {state_file}")
        else:
            start = replica_pdbs[r] if replica_pdbs else pdb_file
            spdb = PDBFile(start)
            engine.set_coordinates(spdb.getPositions(), spdb.topology.getPeriodicBoxVectors()
                                   if not implicit else None)
            eq = float(inputs.get("equil_time", 100.0))
            print(f"[STEP] Replica {r}: {eq:g} ps unbiased equilibration from {start}")
            state = engine.equilibrate(eq, seed=base_seed + 1000 * r + 1,
                                       report=str(rdir / "equilibration.dcd"))
            state_file.write_text(openmm.XmlSerializer.serialize(state))
        for t in range(n_traj):
            if (r, t) in done:
                continue
            seed = base_seed + 1000 * r + 10 * t + 7
            print(f"[STEP] Replica {r}, trajectory {t} (seed {seed})")
            out = engine.run_trajectory(
                state, seed, max_time,
                dcd=str(rdir / f"traj_{t:02d}.dcd") if inputs.get("traj_interval", 0) else None,
                dcd_interval_ps=float(inputs.get("traj_interval", 0)),
                exit_pdb=str(rdir / f"traj_{t:02d}_end.pdb") if inputs.get("exit_pdb", True) else None)
            row = {"replica": r, "trajectory": t, "seed": seed, **out}
            append_row(csvfile, row)
            print(f"[INFO] {'dissociated' if out['dissociated'] else 'NOT dissociated'} after "
                  f"{out['time_ps']:.1f} ps ({out['direction_changes']} direction changes, "
                  f"{out['wall_s']:.0f} s wall)")

    res = analyse(read_times(csvfile), inputs)
    print(format_report(res))
    result.update(res)
    save_result(result, filename=args.config, section=section)
    return result
