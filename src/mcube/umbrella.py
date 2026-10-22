#!/usr/bin/env python3
"""
Umbrella sampling of one atom (an ion) along the membrane normal, on the MM or
the QM/MM model, and the free-energy analysis that goes with it.

The collective variable is the z displacement of the ion from the mass-weighted
centre of a reference group (the channel), z = z_ion - z_COM; a flat-bottomed
cylinder on the xy distance keeps the ion over the pore in bulk. Both are one
CustomCVForce added to the OpenMM system after everything else (after the
QMMMTheory, whose setup zeroes bonded terms in forces already present), with
uniquely named parameters, in force group 31 so energies can leave it out.

Every window is its own job, and the stages of the method are actions:

    mcube umbrella config.yaml pull                          steered MM run that seeds every window
    mcube umbrella config.yaml window --window 37            one MM window (OpenMM, GPU)
    mcube umbrella config.yaml window --window 7 --section umbrella_qmmm    one QM/MM window (ASH)
    mcube umbrella config.yaml energies --window 7 --level xtb [--section S]  per-frame energies
    mcube umbrella config.yaml wham [--section S] [--reweight dft]           PMF
    mcube umbrella config.yaml stitch                        MM <-> QM/MM join at the constriction ends

Inputs are `umbrella.inputs`, overridden by the `--section` inputs (the same
layering as scan/neb/ts), so the QM/MM leg states only what differs.

MM windows run on OpenMM directly, with the system the md stage builds
(membrane barostat, HBond constraints, 1.5 amu hydrogens): ASH's pure-MM MD has
no per-step hook, so it could not log the CV without writing the whole system
every few steps. QM/MM windows run through ASH (`OpenMM_MD` on a QMMMTheory) and
log the CV from a trajectory of the ion and reference atoms only.

The metal's hydration shell, as preprocess would detect it in each window's own
starting frame, is restrained during the MM equilibration only and released for
sampling; QM/MM is never restrained. The QM region of a QM/MM window is the
ion, `qm_waters` nearest waters (default: the detected shell) and the fixed
`qm_residues`, rebuilt on each window's starting frame; its charge is the
sum of the MM charges it replaces, unless `charge` is given.

Results are written per window under `<section>.results.windows.<NNN>` with a
locked read-modify-write, so window jobs can share one config.
"""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import yaml

from .preprocess import (
    COORDINATION_FORCE,
    _coords_and_cell,
    _locked,
    build_qm_region,
    compare_coordination,
    coordination_shell,
    load_section,
    qm_region_settings,
    resolve_atoms,
    resolve_cores,
)

KB = 0.0019872041          # kcal/mol/K
HARTREE_TO_KCAL = 627.509474
TO_KJ_NM2 = 4.184 * 100.0   # kcal/mol/A^2 -> kJ/mol/nm^2
BIAS_GROUP = 31
BIAS_NAME = "mcube_umbrella"
COORDINATION_NAME = "mcube_coordination"

DEFAULTS = {
    "k": 10.0,                # kcal/mol/A^2
    "engine": "mm",
    "pull_speed": 10.0,       # A/ns
    "equilibrate": 100.0,     # ps
    "time": 5000.0,           # ps
    "timestep": 0.004,        # ps
    "cv_frequency": 25,       # steps
    "traj_frequency": 25000,  # steps
    "platform": "CPU",
}


# ==============================
# Config
# ==============================

def umbrella_inputs(config: str, section: Optional[str] = None) -> Dict[str, Any]:
    """`umbrella.inputs`, overridden by `section`'s inputs, over the defaults above."""
    with open(config) as fh:
        data = yaml.safe_load(fh) or {}
    inputs = dict(DEFAULTS)
    inputs.update((data.get("umbrella") or {}).get("inputs") or {})
    if section and section != "umbrella":
        inputs.update((data.get(section) or {}).get("inputs") or {})
    return inputs


def temperature_of(inputs: Dict[str, Any], md: Dict[str, Any]) -> float:
    """K: the section's `temp`, else md's (the windows sample md's ensemble)."""
    return float(inputs.get("temp", (md.get("inputs") or {}).get("temp", 300.0)))


def window_centres(inputs: Dict[str, Any]) -> List[float]:
    """Window centres in A, from z_start to z_end inclusive, `spacing` apart."""
    start, end, spacing = (float(inputs[k]) for k in ("z_start", "z_end", "spacing"))
    if spacing <= 0 or end <= start:
        raise ValueError("umbrella needs z_start < z_end and spacing > 0")
    n = int(round((end - start) / spacing)) + 1
    return [round(start + i * spacing, 6) for i in range(n)]


def window_key(index: int) -> str:
    return f"{int(index):03d}"


def section_windows(config: str, section: str) -> Dict[str, Dict[str, Any]]:
    return ((load_section(config, section).get("results") or {}).get("windows")) or {}


def find_window(config: str, section: str, z0: float) -> Optional[Dict[str, Any]]:
    """The window entry of `section` centred on z0, if there is one."""
    for entry in section_windows(config, section).values():
        if entry.get("z0") is not None and abs(float(entry["z0"]) - z0) < 1e-6:
            return entry
    return None


def save_nested(config: str, section: str, keys: List[str], value) -> None:
    """Set config[section]['results'][keys...] = value (dicts merge), under the config lock.

    `save_result` replaces a section's whole `results`, which would let one
    window job erase another's; this changes only its own entry.
    """
    path = Path(config)
    with _locked(path):
        with open(path) as fh:
            data = yaml.safe_load(fh) or {}
        node = data.setdefault(section, {})
        if not isinstance(node.get("results"), dict):
            node["results"] = {}
        node = node["results"]
        for k in keys[:-1]:
            if not isinstance(node.get(k), dict):
                node[k] = {}
            node = node[k]
        if isinstance(value, dict) and isinstance(node.get(keys[-1]), dict):
            node[keys[-1]].update(value)
        else:
            node[keys[-1]] = value
        tmp = path.with_name(f"{path.name}.{os.getpid()}.tmp")
        with open(tmp, "w") as fh:
            yaml.safe_dump(data, fh, sort_keys=False, default_flow_style=False)
        tmp.replace(path)
    print(f"[INFO] Updated -> {path} [{section}.results.{'.'.join(keys)}]")


# ==============================
# Geometry and the bias
# ==============================

def reference_atoms(topology, selection) -> List[int]:
    """Reference group: an mdtraj selection string (PDB numbers are `resSeq`), or a
    list of indices / `chain:resid:atom` selectors."""
    if isinstance(selection, str):
        import mdtraj
        atoms = [int(a) for a in mdtraj.Topology.from_openmm(topology).select(selection)]
    else:
        atoms = [int(a) for a in resolve_atoms(list(selection), topology)]
    if not atoms:
        raise ValueError(f"umbrella reference {selection!r} selects no atoms")
    return atoms


def ion_atom(topology, selector) -> int:
    return int(resolve_atoms(selector, topology))


def atom_masses(topology) -> np.ndarray:
    return np.array([a.element.mass._value if a.element else 0.0 for a in topology.atoms()])


def displacement(coords, ion: int, reference, masses, box) -> np.ndarray:
    """Ion minus reference COM (A), per frame, minimum image on an orthorhombic box.

    `coords` is (frames, atoms, 3) or (atoms, 3), `box` (3,) or (frames, 3), in A.
    The reference group is unwrapped about its first atom first, so a group that
    a trajectory writer split across the cell still has the right centre.
    """
    coords = np.asarray(coords, dtype=float)
    single = coords.ndim == 2
    if single:
        coords = coords[None]
    box = np.broadcast_to(np.asarray(box, dtype=float), (coords.shape[0], 3))[:, None, :]
    ref = coords[:, reference, :]
    ref = ref - box * np.round((ref - ref[:, :1, :]) / box)
    m = np.asarray(masses, dtype=float)[reference][None, :, None]
    com = (ref * m).sum(axis=1) / m.sum(axis=1)
    d = coords[:, ion, :] - com
    d = d - box[:, 0, :] * np.round(d / box[:, 0, :])
    return d[0] if single else d


def cv_values(coords, ion, reference, masses, box):
    """(z, r_xy) in A for each frame."""
    d = np.atleast_2d(displacement(coords, ion, reference, masses, box))
    return d[:, 2], np.hypot(d[:, 0], d[:, 1])


def image_ion(coords, ion: int, reference, masses, box) -> np.ndarray:
    """Copy of coords with the ion moved to the image nearest the reference COM.

    The bias works on raw positions (no periodic wrapping in the force), so the
    ion and the channel must start in the same image; OpenMM never wraps the
    positions it integrates, so they stay that way for the run.
    """
    coords = np.array(coords, dtype=float)
    box = np.asarray(box, dtype=float)
    ref = coords[reference]
    if np.any(np.ptp(ref, axis=0) > box / 2):
        raise ValueError("the umbrella reference group spans more than half the box: it is "
                         "split across the cell; give a structure with the channel whole")
    m = np.asarray(masses)[reference]
    com = (ref * m[:, None]).sum(axis=0) / m.sum()
    d = coords[ion] - com
    coords[ion] -= box * np.round(d / box)
    return coords


def bias_force(ion: int, reference, k: float, z0: float, cylinder=None, masses=None):
    """The umbrella: 0.5 k (z - z0)^2 + 0.5 k_cyl max(0, r - R)^2 (k in kcal/mol/A^2, A).

    The reference centre is weighted by element masses (`masses`, per atom), not by
    the System's, which carry repartitioned hydrogen masses: the analysis computes
    the same centre from the trajectory, so the two must agree.
    """
    import openmm
    z = openmm.CustomCentroidBondForce(2, "z2-z1")
    r = openmm.CustomCentroidBondForce(2, "sqrt((x2-x1)^2+(y2-y1)^2)")
    for f in (z, r):
        if masses is None:
            f.addGroup([int(a) for a in reference])
        else:
            f.addGroup([int(a) for a in reference], [float(masses[int(a)]) for a in reference])
        f.addGroup([int(ion)])
        f.addBond([0, 1], [])
    cyl = cylinder or {}
    force = openmm.CustomCVForce("0.5*mcube_us_k*(zcv-mcube_us_z0)^2"
                                 " + 0.5*mcube_cyl_k*max(0, rcv-mcube_cyl_r)^2")
    force.addCollectiveVariable("zcv", z)
    force.addCollectiveVariable("rcv", r)
    force.addGlobalParameter("mcube_us_k", float(k) * TO_KJ_NM2)
    force.addGlobalParameter("mcube_us_z0", float(z0) / 10.0)
    force.addGlobalParameter("mcube_cyl_k", float(cyl.get("k", 0.0)) * TO_KJ_NM2)
    force.addGlobalParameter("mcube_cyl_r", float(cyl.get("radius", 1000.0)) / 10.0)
    force.setName(BIAS_NAME)
    force.setForceGroup(BIAS_GROUP)
    return force


def coordination_force(pdb_file: str, preprocess_inputs: Dict[str, Any]):
    """Harmonic restraints on the shell preprocess would detect in this frame, at the
    distances it has. Returns (force, shell)."""
    import openmm
    shell = coordination_shell(pdb_file, preprocess_inputs)
    force = openmm.CustomBondForce("0.5*k_coord*(r-r0_coord)^2")
    force.addPerBondParameter("k_coord")
    force.addPerBondParameter("r0_coord")
    for row in shell:
        force.addBond(int(row["metal_index"]), int(row["donor_index"]),
                      [COORDINATION_FORCE * TO_KJ_NM2, float(row["distance"]) / 10.0])
    force.setUsesPeriodicBoundaryConditions(True)
    force.setName(COORDINATION_NAME)
    return force, shell


def release(force, context) -> None:
    """Switch the coordination restraints off for sampling."""
    for i in range(force.getNumBonds()):
        a, b, (_, r0) = force.getBondParameters(i)
        force.setBondParameters(i, a, b, [0.0, r0])
    force.updateParametersInContext(context)


# ==============================
# The MM model
# ==============================

def pdb_box(pdb) -> np.ndarray:
    """Orthorhombic box (A) of a PDBFile."""
    cell = pdb.topology.getUnitCellDimensions()
    if cell is None:
        raise ValueError("umbrella needs a periodic structure (CRYST1)")
    return np.array([cell[i]._value for i in range(3)]) * 10.0


def mm_system(pdb, inputs: Dict[str, Any], preprocess: Dict[str, Any], md: Dict[str, Any]):
    """The md stage's classical system for this structure, with the same metal parameters.
    Returns (system, forcefields)."""
    from openmm import unit
    from openmm.app import ForceField, PME, HBonds
    from . import ions
    from . import membrane as membrane_mod
    from .md import md_forcefields, md_membrane

    md_in = md.get("inputs") or {}
    if md_in.get("implicit"):
        raise NotImplementedError("umbrella sampling needs explicit solvent")
    forcefields = (md.get("results") or {}).get("forcefields") or md_forcefields(md_in, preprocess)
    temp = temperature_of(inputs, md)
    cutoff = float(inputs.get("nonbonded_cutoff", md_in.get("nonbonded_cutoff", membrane_mod.NONBONDED_CUTOFF)))
    if md_membrane(md_in, preprocess):
        system = membrane_mod.membrane_system(
            pdb.topology, forcefields, temperature=temp, nonbonded_cutoff=cutoff,
            pressure=md_in.get("pressure", membrane_mod.PRESSURE),
            surface_tension=md_in.get("surface_tension", membrane_mod.SURFACE_TENSION),
            barostat_frequency=md_in.get("barostat_frequency", membrane_mod.BAROSTAT_FREQUENCY))
    else:
        # md's explicit production is NVT, so the windows are too
        system = ForceField(*forcefields).createSystem(
            pdb.topology, nonbondedMethod=PME, nonbondedCutoff=cutoff / 10.0 * unit.nanometer,
            constraints=HBonds, rigidWater=True, hydrogenMass=1.5 * unit.amu)
    metals = (preprocess.get("inputs") or {}).get("metals")
    which = ions.metal_parameter_set(inputs, md)
    if metals and which != "none":
        ions.apply_metal_parameters(system, pdb.topology, metals, which, forcefields=forcefields,
                                    polarizabilities=inputs.get("polarizabilities",
                                                                md_in.get("polarizabilities")))
    return system, forcefields


def make_simulation(pdb, system, positions_A, inputs, md):
    import openmm
    from openmm import unit
    from openmm.app import Simulation
    temp = temperature_of(inputs, md)
    integrator = openmm.LangevinMiddleIntegrator(temp * unit.kelvin, 1.0 / unit.picosecond,
                                                 float(inputs["timestep"]) * unit.picoseconds)
    platform = openmm.Platform.getPlatformByName(inputs.get("platform", "CPU"))
    props = dict(inputs.get("properties") or {})
    if inputs.get("platform", "CPU") == "CPU":
        # OpenMM's CPU platform otherwise starts one thread per core of the whole
        # node: on a 2-core allocation of a 128-core node a window ran ~2000 threads
        props.setdefault("Threads", str(resolve_cores(inputs.get("cores"))))
    simulation = Simulation(pdb.topology, system, integrator, platform, props)
    simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())
    simulation.context.setPositions(np.asarray(positions_A) / 10.0 * unit.nanometer)
    return simulation, temp


def box_vectors(vectors_nm):
    """Three openmm Vec3 box vectors (nm) from a 3x3 array."""
    from openmm import Vec3, unit
    return [Vec3(*(float(x) for x in v)) * unit.nanometer for v in vectors_nm]


def write_state_pdb(simulation, topology, filename: str) -> str:
    from openmm.app import PDBFile
    state = simulation.context.getState(getPositions=True)
    topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
    with open(filename, "w") as fh:
        PDBFile.writeFile(topology, state.getPositions(), fh, keepIds=True)
    return Path(filename).resolve().as_posix()


def write_coords_pdb(topology, coords_A, filename: str) -> str:
    from openmm import unit
    from openmm.app import PDBFile
    with open(filename, "w") as fh:
        PDBFile.writeFile(topology, np.asarray(coords_A) / 10.0 * unit.nanometer, fh, keepIds=True)
    return Path(filename).resolve().as_posix()


def cv_now(force, context):
    z, r = force.getCollectiveVariableValues(context)
    return 10.0 * z, 10.0 * r


# ==============================
# Actions
# ==============================

def run_pull(config: str, section: str, inputs, preprocess, md) -> Dict[str, Any]:
    """Steered MM run along z that saves a seed structure as z0 passes each centre."""
    from openmm.app import PDBFile
    pdb_file = inputs.get("pdb_file") or (md.get("results") or {}).get("lastframe")
    if not pdb_file:
        raise ValueError("pull needs a starting structure: umbrella.inputs.pdb_file or md.results.lastframe")
    pdb_file = Path(pdb_file).resolve().as_posix()
    work = Path("pull")
    work.mkdir(exist_ok=True)
    os.chdir(work)

    pdb = PDBFile(pdb_file)
    box = pdb_box(pdb)
    ion = ion_atom(pdb.topology, inputs["ion"])
    reference = reference_atoms(pdb.topology, inputs["reference"])
    masses = atom_masses(pdb.topology)
    coords = image_ion(np.array(pdb.getPositions()._value) * 10.0, ion, reference, masses, box)

    system, _ = mm_system(pdb, inputs, preprocess, md)
    z_start = float(cv_values(coords, ion, reference, masses, box)[0][0])
    k = float(inputs.get("pull_k", inputs["k"]))
    bias = bias_force(ion, reference, k, z_start, inputs.get("cylinder"), masses)
    system.addForce(bias)
    simulation, temp = make_simulation(pdb, system, coords, inputs, md)
    print("[STEP] Minimisation...")
    simulation.minimizeEnergy(maxIterations=500)
    simulation.context.setVelocitiesToTemperature(temp)

    centres = window_centres(inputs)
    near, far = sorted((centres[0], centres[-1]), key=lambda c: abs(c - z_start))
    legs = [(z_start, near), (near, far)]
    chunk = int(inputs["cv_frequency"])
    dt = float(inputs["timestep"])
    step_dz = float(inputs["pull_speed"]) / 1000.0 * dt * chunk   # A per chunk
    print(f"[INFO] Pull: CV starts at {z_start:.2f} A; {near:.2f} -> {far:.2f} A at "
          f"{inputs['pull_speed']} A/ns, k {k} kcal/mol/A^2; {len(centres)} window centres")

    seeds, z0, t = {}, z_start, 0.0
    with open("pull_cv.csv", "w", newline="") as fh:
        out = csv.writer(fh)
        out.writerow(["time_ps", "z0_A", "z_A", "r_A"])
        for a, b in legs:
            direction = 1.0 if b >= a else -1.0
            n = int(math.ceil(abs(b - a) / step_dz)) if b != a else 0
            for i in range(1, n + 1):
                previous = z0
                z0 = a + direction * min(i * step_dz, abs(b - a))
                simulation.context.setParameter("mcube_us_z0", z0 / 10.0)
                simulation.step(chunk)
                t += chunk * dt
                z, r = cv_now(bias, simulation.context)
                out.writerow([round(t, 4), round(z0, 4), round(z, 4), round(r, 4)])
                lo, hi = sorted((previous, z0))
                for index, c in enumerate(centres):
                    key = window_key(index)
                    if key not in seeds and lo - 1e-9 <= c <= hi + 1e-9:
                        seed = write_state_pdb(simulation, pdb.topology, f"seed_{key}.pdb")
                        seeds[key] = {"z0": c, "seed": seed, "seed_cv": round(z, 3)}
                        print(f"[INFO] seed {key}: z0 {c:.2f} A, CV {z:.2f} A")

    result = {"start": Path(pdb_file).resolve().as_posix(), "z_start": round(z_start, 3),
              "cv": Path("pull_cv.csv").resolve().as_posix(), "k": k,
              "speed_A_per_ns": float(inputs["pull_speed"]), "seeds": len(seeds)}
    save_nested(config, section, ["pull"], result)
    for key, entry in seeds.items():
        save_nested(config, section, ["windows", key], entry)
    return result


def window_seed(config: str, section: str, inputs, z0: float, key: str) -> str:
    """This window's starting structure: its own seed, else (QM/MM) the MM window at
    the same z0, else the pull seed for that z0."""
    own = section_windows(config, section).get(key) or {}
    if own.get("seed"):
        return own["seed"]
    if section != "umbrella":
        mm = find_window(config, "umbrella", z0)
        if mm and mm.get("lastframe"):
            return mm["lastframe"]
        if mm and mm.get("seed"):
            return mm["seed"]
    raise ValueError(f"no starting structure for window {key} (z0 {z0} A): run `pull` first, "
                     "or give `seed` in its results")


def run_mm_window(config, section, inputs, preprocess, md, index: int) -> Dict[str, Any]:
    from openmm.app import PDBFile, DCDReporter, StateDataReporter
    centres = window_centres(inputs)
    z0, key = centres[index], window_key(index)
    seed = Path(window_seed(config, section, inputs, z0, key)).resolve().as_posix()
    work = Path(f"{section}_{key}")
    work.mkdir(exist_ok=True)
    os.chdir(work)

    pdb = PDBFile(seed)
    box = pdb_box(pdb)
    ion = ion_atom(pdb.topology, inputs["ion"])
    reference = reference_atoms(pdb.topology, inputs["reference"])
    masses = atom_masses(pdb.topology)
    coords = image_ion(np.array(pdb.getPositions()._value) * 10.0, ion, reference, masses, box)
    start = write_coords_pdb(pdb.topology, coords, "start.pdb")

    system, _ = mm_system(pdb, inputs, preprocess, md)
    bias = bias_force(ion, reference, float(inputs["k"]), z0, inputs.get("cylinder"), masses)
    system.addForce(bias)
    prep_in = preprocess.get("inputs") or {}
    restraint, shell = (coordination_force(start, prep_in) if prep_in.get("metals") else (None, []))
    if restraint is not None:
        system.addForce(restraint)
    simulation, temp = make_simulation(pdb, system, coords, inputs, md)
    print(f"[INFO] Window {key}: z0 {z0} A, k {inputs['k']} kcal/mol/A^2, "
          f"{len(shell)} shell contacts restrained for equilibration")
    simulation.minimizeEnergy(maxIterations=500)
    simulation.context.setVelocitiesToTemperature(temp)

    dt, chunk = float(inputs["timestep"]), int(inputs["cv_frequency"])
    n_eq = int(round(float(inputs["equilibrate"]) / dt))
    if n_eq:
        print(f"[STEP] Equilibration, {inputs['equilibrate']} ps, shell restrained...")
        simulation.step(n_eq)
    if restraint is not None:
        release(restraint, simulation.context)

    n_prod = int(round(float(inputs["time"]) / dt))
    freq = int(inputs["traj_frequency"])
    simulation.reporters.append(DCDReporter("window.dcd", freq))
    simulation.reporters.append(StateDataReporter("window.csv", freq, step=True, time=True,
                                                  potentialEnergy=True, temperature=True,
                                                  volume=True))
    print(f"[STEP] Sampling, {inputs['time']} ps, CV every {chunk} steps...")
    with open("cv.csv", "w", newline="") as fh:
        out = csv.writer(fh)
        out.writerow(["time_ps", "z_A", "r_A"])
        for i in range(n_prod // chunk):
            simulation.step(chunk)
            z, r = cv_now(bias, simulation.context)
            out.writerow([round((i + 1) * chunk * dt, 4), round(z, 4), round(r, 4)])
    last = write_state_pdb(simulation, pdb.topology, "lastframe.pdb")

    after = coordination_shell(last, prep_in) if prep_in.get("metals") else []
    entry = {"z0": z0, "k": float(inputs["k"]), "engine": "mm", "seed": seed, "start": start,
             "cv": Path("cv.csv").resolve().as_posix(), "dcd": Path("window.dcd").resolve().as_posix(),
             "lastframe": last, "equilibrate_ps": float(inputs["equilibrate"]),
             "time_ps": float(inputs["time"]), "timestep": dt,
             "traj_frequency": freq, "shell_start": len(shell), "shell_end": len(after),
             "done": True}
    save_nested(config, section, ["windows", key], entry)
    return entry


def umbrella_qm_region(pdb_file: str, inputs, preprocess, ion: int):
    """QM atoms for a QM/MM window: the ion, `qm_waters` nearest waters (default: the
    detected shell's residues) and the fixed `qm_residues`. Returns (atoms, boundary, residues)."""
    from openmm.app import PDBFile
    prep_in = preprocess.get("inputs") or {}
    topology = PDBFile(pdb_file).topology
    atoms = list(topology.atoms())
    residues: Dict[str, Dict[int, str]] = {}

    def add(res):
        residues.setdefault(res.chain.id, {})[int(res.id)] = res.name

    add(atoms[ion].residue)
    n_waters = inputs.get("qm_waters")
    if n_waters is None:
        for row in coordination_shell(pdb_file, prep_in):
            if row["metal_index"] == ion:
                add(atoms[row["donor_index"]].residue)
    else:
        coords, cell = _coords_and_cell(pdb_file)
        coords = np.asarray(coords)
        oxygens = [a.index for a in atoms if a.residue.name == "HOH" and a.element is not None
                   and a.element.symbol == "O"]
        d = coords[oxygens] - coords[ion]
        if cell:
            d -= np.asarray(cell) * np.round(d / np.asarray(cell))
        for i in np.argsort(np.linalg.norm(d, axis=1))[:int(n_waters)]:
            add(atoms[oxygens[i]].residue)
    for chain, resids in (inputs.get("qm_residues") or {}).items():
        found = {int(r.id): r for r in topology.residues() if r.chain.id == str(chain)}
        missing = [r for r in resids if int(r) not in found]
        if missing:
            raise ValueError(f"qm_residues: chain {chain} residues {missing} not in {pdb_file}")
        for r in resids:
            add(found[int(r)])
    qm_atoms, boundary = build_qm_region(pdb_file, residues, **qm_region_settings(prep_in))
    return qm_atoms, boundary, residues


def qm_charge(system, qm_atoms) -> tuple:
    """(integer charge, raw sum) of the MM charges the QM region replaces."""
    import openmm
    nb = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    total = sum(nb.getParticleParameters(int(i))[0]._value for i in qm_atoms)
    return int(round(total)), round(total, 4)


def qmmm_level_inputs(config: str, inputs) -> Dict[str, Any]:
    """qmmm.inputs (charge, mult, xTB, cores ...) with this section's QM keys on top."""
    from .qmmm import qm_level_inputs, stage_inputs
    qin = qm_level_inputs(stage_inputs(config, "qmmm"), inputs)
    for k in ("charge", "mult", "cores", "metal_parameters", "polarizabilities"):
        if k in inputs:
            qin[k] = inputs[k]
    return qin


def read_special_xyz(path: str) -> np.ndarray:
    """Frames of an XYZ file (A), as (frames, atoms, 3)."""
    frames, lines = [], open(path).read().split("\n")
    i = 0
    while i < len(lines) and lines[i].strip():
        n = int(lines[i])
        frames.append([[float(x) for x in l.split()[1:4]] for l in lines[i + 2:i + 2 + n]])
        i += n + 2
    return np.array(frames)


def run_qmmm_window(config, section, inputs, preprocess, md, index: int) -> Dict[str, Any]:
    from ash import Fragment, OpenMM_MD
    from openmm.app import PDBFile
    from .qmmm import apply_stage_metal_parameters, make_mm, make_qmmm, mm_settings

    centres = window_centres(inputs)
    z0, key = centres[index], window_key(index)
    seed = Path(window_seed(config, section, inputs, z0, key)).resolve().as_posix()
    work = Path(f"{section}_{key}")
    work.mkdir(exist_ok=True)
    os.chdir(work)

    pdb = PDBFile(seed)
    box = pdb_box(pdb)
    ion = ion_atom(pdb.topology, inputs["ion"])
    reference = reference_atoms(pdb.topology, inputs["reference"])
    masses = atom_masses(pdb.topology)
    start = write_coords_pdb(pdb.topology,
                             image_ion(np.array(pdb.getPositions()._value) * 10.0, ion, reference,
                                       masses, box), "start.pdb")

    qin = qmmm_level_inputs(config, inputs)
    forcefields, implicit = mm_settings(qin, md)
    cores = resolve_cores(qin.get("cores"))
    mm = make_mm(forcefields, start, implicit, cores, inputs.get("platform", "CPU"),
                 inputs.get("properties"))
    qm_atoms, boundary, residues = umbrella_qm_region(start, inputs, preprocess, ion)
    charge, raw = qm_charge(mm.system, qm_atoms)
    if "charge" in inputs:
        charge = int(inputs["charge"])
    qin["charge"] = charge
    print(f"[INFO] QM region: {len(qm_atoms)} atoms, charge {charge} (MM charges sum {raw}), "
          f"mult {qin.get('mult')}")
    apply_stage_metal_parameters(qin, preprocess, md, mm, forcefields, qm_atoms)
    fragment = Fragment(pdbfile=start)
    theory = make_qmmm(qin, mm, fragment, qm_atoms, boundary, cores)
    # After the QMMMTheory: its setup edits the forces already in the system
    mm.system.addForce(bias_force(ion, reference, float(inputs["k"]), z0, inputs.get("cylinder"),
                                  masses))

    special = [ion] + list(reference)
    dt, chunk = float(inputs["timestep"]), int(inputs["cv_frequency"])
    Path("wrapped_special_traj.xyz").unlink(missing_ok=True)
    # One OpenMM_MD call: each call on a QMMMTheory adds another external force to the
    # system, so equilibration is the discarded start of this run, not a second run.
    OpenMM_MD(fragment=fragment, theory=theory, timestep=dt,
              simulation_time=float(inputs["equilibrate"]) + float(inputs["time"]),
              traj_frequency=int(inputs["traj_frequency"]), temperature=temperature_of(inputs, md),
              integrator="LangevinMiddleIntegrator", coupling_frequency=1, trajfilename="window",
              charge=charge, mult=qin.get("mult"),
              specialatoms=special, specialtraj_frequency=chunk)

    frames = read_special_xyz("wrapped_special_traj.xyz")
    sub_masses = np.concatenate([[masses[ion]], masses[reference]])
    z, r = cv_values(frames, 0, list(range(1, len(special))), sub_masses, box)
    with open("cv.csv", "w", newline="") as fh:
        out = csv.writer(fh)
        out.writerow(["time_ps", "z_A", "r_A"])
        for i, (zi, ri) in enumerate(zip(z, r)):
            out.writerow([round(i * chunk * dt, 4), round(float(zi), 4), round(float(ri), 4)])

    last = Path("window_lastframe.pdb")
    lastframe = last.resolve().as_posix() if last.exists() else None
    prep_in = preprocess.get("inputs") or {}
    changes = (compare_coordination(coordination_shell(start, prep_in),
                                    coordination_shell(lastframe, prep_in))
               if lastframe and prep_in.get("metals") else [])
    entry = {"z0": z0, "k": float(inputs["k"]), "engine": "qmmm", "seed": seed, "start": start,
             "cv": Path("cv.csv").resolve().as_posix(),
             "dcd": Path("window.dcd").resolve().as_posix(), "lastframe": lastframe,
             "equilibrate_ps": float(inputs["equilibrate"]), "time_ps": float(inputs["time"]),
             "timestep": dt, "traj_frequency": int(inputs["traj_frequency"]),
             "qm_atoms": len(qm_atoms), "qm_residues": residues, "charge": charge,
             "mult": qin.get("mult"), "shell_changes": changes, "done": True}
    save_nested(config, section, ["windows", key], entry)
    return entry


def run_energies(config, section, inputs, preprocess, md, index: int, level: str,
                 stride: int, qm_section: str) -> Dict[str, Any]:
    """Unbiased potential energy (kcal/mol) of every `stride`-th saved frame of a window,
    at `mm`, `xtb` (xTB/MM) or `dft` (DFT/MM). The QM region of each frame follows
    `qm_section`'s rule, so MM and QM/MM frames are scored by the same Hamiltonian."""
    import mdtraj
    from openmm import unit
    from openmm.app import PDBFile

    key = window_key(index)
    entry = section_windows(config, section).get(key) or {}
    if not entry.get("done"):
        raise ValueError(f"window {key} of {section} has not finished")
    traj = mdtraj.load(entry["dcd"], top=entry["lastframe"])[::max(1, int(stride))]
    pdb = PDBFile(entry["lastframe"])
    ion = ion_atom(pdb.topology, inputs["ion"])
    reference = reference_atoms(pdb.topology, inputs["reference"])
    # Writers keep molecules whole but image each one on its own, so the ion's
    # waters can sit a box length away from it. Re-image every molecule about the
    # channel chains, as one QM region needs (PME energies do not change).
    chains = {traj.topology.atom(i).residue.chain.index for i in reference}
    anchors = [set(c.atoms) for c in traj.topology.chains if c.index in chains]
    traj.image_molecules(inplace=True, anchor_molecules=anchors)
    masses = atom_masses(pdb.topology)
    xyz = traj.xyz * 10.0
    boxes = traj.unitcell_lengths * 10.0
    z = cv_values(xyz, ion, reference, masses, boxes)[0]
    times = traj.time
    os.chdir(Path(entry["dcd"]).parent)

    energies = []
    if level == "mm":
        system, _ = mm_system(pdb, inputs, preprocess, md)
        simulation, _ = make_simulation(pdb, system, xyz[0], dict(inputs, timestep=0.001), md)
        groups = {g for g in range(32) if g != BIAS_GROUP}
        for f in range(traj.n_frames):
            simulation.context.setPeriodicBoxVectors(*box_vectors(traj.unitcell_vectors[f]))
            simulation.context.setPositions(traj.xyz[f] * unit.nanometer)
            e = simulation.context.getState(getEnergy=True, groups=groups).getPotentialEnergy()
            energies.append(e.value_in_unit(unit.kilocalorie_per_mole))
    elif level in ("xtb", "dft"):
        from ash import Fragment, QMMMTheory, Singlepoint
        from .qmmm import (apply_stage_metal_parameters, make_dft, make_mm, make_qmmm,
                           mm_settings)
        qm_rule = umbrella_inputs(config, qm_section)
        qin = qmmm_level_inputs(config, qm_rule)
        forcefields, implicit = mm_settings(qin, md)
        cores = resolve_cores(qin.get("dft_cores" if level == "dft" else "cores", qin.get("cores")))
        for f in range(traj.n_frames):
            frame_pdb = f"energy_frame_{f:04d}.pdb"
            pdb.topology.setPeriodicBoxVectors(box_vectors(traj.unitcell_vectors[f]))
            write_coords_pdb(pdb.topology, xyz[f], frame_pdb)
            mm = make_mm(forcefields, frame_pdb, implicit, cores)
            qm_atoms, boundary, _ = umbrella_qm_region(frame_pdb, qm_rule, preprocess, ion)
            charge = int(qm_rule["charge"]) if "charge" in qm_rule else qm_charge(mm.system, qm_atoms)[0]
            apply_stage_metal_parameters(qin, preprocess, md, mm, forcefields, qm_atoms)
            fragment = Fragment(pdbfile=frame_pdb)
            if level == "xtb":
                theory = make_qmmm(dict(qin, charge=charge), mm, fragment, qm_atoms, boundary, cores)
            else:
                theory = QMMMTheory(qm_theory=make_dft(qin, cores, qin.get("mult"), "sp"),
                                    mm_theory=mm, fragment=fragment, qmatoms=qm_atoms,
                                    excludeboundaryatomlist=boundary, embedding="electrostatic",
                                    printlevel=1)
            res = Singlepoint(theory=theory, fragment=fragment, charge=charge, mult=qin.get("mult"))
            energies.append(float(res.energy) * HARTREE_TO_KCAL)
            Path(frame_pdb).unlink(missing_ok=True)
    else:
        raise ValueError("level must be mm, xtb or dft")

    out_csv = Path(f"energies_{level}.csv").resolve()
    with open(out_csv, "w", newline="") as fh:
        out = csv.writer(fh)
        out.writerow(["frame", "time_ps", "z_A", "energy_kcal"])
        for f, (t, zi, e) in enumerate(zip(times, z, energies)):
            out.writerow([f * max(1, int(stride)), round(float(t), 4), round(float(zi), 4), e])
    save_nested(config, section, ["windows", key], {f"energies_{level}": out_csv.as_posix(),
                                                    f"energies_{level}_stride": int(stride)})
    return {"energies": out_csv.as_posix(), "frames": len(energies)}


# ==============================
# Free energy
# ==============================

def wham(windows, temperature: float, bin_width: float, tol: float = 1e-8,
         max_iter: int = 200000):
    """
    WHAM (Kumar et al. 1992) on a 1D grid. `windows` is a list of
    (z0, k, samples, weights-or-None), k in kcal/mol/A^2 and z in A. Per-sample
    weights (e.g. exp(-beta dE) for reweighting to another level) are normalised
    to each window's sample count, so a window still counts by its samples.
    Returns (bin centres, F in kcal/mol with min 0, NaN where unsampled, f_i).
    """
    beta = 1.0 / (KB * temperature)
    allz = np.concatenate([np.asarray(w[2], dtype=float) for w in windows])
    edges = np.arange(allz.min(), allz.max() + bin_width, bin_width)
    if len(edges) < 2:
        edges = np.array([allz.min() - bin_width / 2, allz.min() + bin_width / 2])
    centres = 0.5 * (edges[1:] + edges[:-1])
    counts, n = [], []
    for z0, k, samples, weights in windows:
        samples = np.asarray(samples, dtype=float)
        w = np.ones_like(samples) if weights is None else np.asarray(weights, dtype=float)
        w = w * len(samples) / w.sum()
        counts.append(np.histogram(samples, bins=edges, weights=w)[0])
        n.append(len(samples))
    counts, n = np.array(counts), np.array(n, dtype=float)
    bias = np.array([0.5 * float(k) * (centres - float(z0)) ** 2 for z0, k, _, _ in windows])
    numerator = counts.sum(axis=0)
    f = np.zeros(len(windows))
    for _ in range(max_iter):
        # log P_k = log sum_i n_ik - log sum_i N_i exp(beta (f_i - U_ik))
        a = np.log(n)[:, None] + beta * (f[:, None] - bias)
        amax = a.max(axis=0)
        log_den = amax + np.log(np.exp(a - amax).sum(axis=0))
        with np.errstate(divide="ignore"):
            log_p = np.log(numerator) - log_den
        b = log_p[None, :] - beta * bias
        b = np.where(np.isfinite(b), b, -np.inf)
        bmax = b.max(axis=1)
        f_new = -(bmax + np.log(np.exp(b - bmax[:, None]).sum(axis=1))) / beta
        f_new -= f_new[0]
        if np.max(np.abs(f_new - f)) < tol:
            f = f_new
            break
        f = f_new
    F = np.where(np.isfinite(log_p), -log_p / beta, np.nan)
    return centres, F - np.nanmin(F), f


def block_bootstrap(windows, temperature, bin_width, n_boot: int = 50, n_blocks: int = 10,
                    seed: int = 0):
    """Standard deviation of F per bin over block-bootstrap resamples of each window."""
    rng = np.random.default_rng(seed)
    centres, F, _ = wham(windows, temperature, bin_width)
    samples = []
    for _ in range(n_boot):
        resampled = []
        for z0, k, z, w in windows:
            z = np.asarray(z)
            w = None if w is None else np.asarray(w)
            blocks = np.array_split(np.arange(len(z)), min(n_blocks, len(z)))
            pick = np.concatenate([blocks[i] for i in rng.integers(0, len(blocks), len(blocks))])
            resampled.append((z0, k, z[pick], None if w is None else w[pick]))
        c, Fb, _ = wham(resampled, temperature, bin_width)
        samples.append(np.interp(centres, c, Fb, left=np.nan, right=np.nan))
    return centres, F, np.nanstd(np.array(samples), axis=0)


def overlaps(windows, bin_width: float) -> List[float]:
    """Histogram overlap sum(min(p_i, p_j)) between neighbouring windows (by z0)."""
    ordered = sorted(windows, key=lambda w: float(w[0]))
    allz = np.concatenate([np.asarray(w[2]) for w in ordered])
    edges = np.arange(allz.min(), allz.max() + bin_width, bin_width)
    hists = [np.histogram(w[2], bins=edges)[0] / len(w[2]) for w in ordered]
    return [float(np.minimum(a, b).sum()) for a, b in zip(hists, hists[1:])]


def bar(w_forward, w_reverse, temperature: float) -> float:
    """Bennett acceptance ratio free energy A1 - A0 (kcal/mol) from forward work
    U1 - U0 on state-0 samples and reverse work U0 - U1 on state-1 samples."""
    beta = 1.0 / (KB * temperature)
    wf, wr = np.asarray(w_forward, dtype=float), np.asarray(w_reverse, dtype=float)
    shift = np.median(np.concatenate([wf, -wr]))      # large constant offsets (QM energies)
    wf, wr = wf - shift, wr + shift
    m = math.log(len(wf) / len(wr))

    def g(d):
        return (np.sum(1.0 / (1.0 + np.exp(np.clip(m + beta * (wf - d), -700, 700))))
                - np.sum(1.0 / (1.0 + np.exp(np.clip(-m + beta * (wr + d), -700, 700)))))

    lo = min(wf.min(), -wr.max()) - 50.0 / beta
    hi = max(wf.max(), -wr.min()) + 50.0 / beta
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(mid) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi) + shift


def read_csv(path: str) -> Dict[str, np.ndarray]:
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]} if rows else {}


def window_samples(entry, reweight: Optional[str], temperature: float):
    """(samples, weights) for WHAM from one window entry."""
    if reweight:
        hi, lo = read_csv(entry[f"energies_{reweight}"]), read_csv(entry["energies_xtb"])
        if not np.array_equal(hi["frame"], lo["frame"]):
            raise ValueError(f"window z0 {entry['z0']}: energies_{reweight} and energies_xtb "
                             "are not on the same frames")
        keep = hi["time_ps"] >= (entry.get("equilibrate_ps", 0.0) if entry.get("engine") == "qmmm" else 0.0)
        de = (hi["energy_kcal"] - lo["energy_kcal"])[keep]
        beta = 1.0 / (KB * temperature)
        return hi["z_A"][keep], np.exp(-beta * (de - de.min()))
    cv = read_csv(entry["cv"])
    skip = entry.get("equilibrate_ps", 0.0) if entry.get("engine") == "qmmm" else 0.0
    keep = cv["time_ps"] >= skip
    return cv["z_A"][keep], None


def run_wham(config, section, inputs, md, reweight: Optional[str], n_boot: int) -> Dict[str, Any]:
    temperature = temperature_of(inputs, md)
    bin_width = float(inputs.get("bin_width", float(inputs["spacing"]) / 5.0))
    windows = []
    for key, entry in sorted(section_windows(config, section).items()):
        if not entry.get("done"):
            continue
        z, w = window_samples(entry, reweight, temperature)
        if len(z):
            windows.append((float(entry["z0"]), float(entry["k"]), z, w))
    if len(windows) < 2:
        raise ValueError(f"{section}: fewer than two finished windows")
    centres, F, err = block_bootstrap(windows, temperature, bin_width, n_boot=n_boot)
    ov = overlaps(windows, bin_width)
    poor = [round(float(sorted(w[0] for w in windows)[i]), 3) for i, o in enumerate(ov) if o < 0.05]
    if poor:
        print(f"[WARNING] window overlap below 0.05 after z0 = {poor}: add windows or soften k there")
    name = f"{section}_pmf{'_' + reweight if reweight else ''}.csv"
    with open(name, "w", newline="") as fh:
        out = csv.writer(fh)
        out.writerow(["z_A", "F_kcal", "err_kcal"])
        for row in zip(centres, F, err):
            out.writerow([round(float(v), 5) for v in row])
    result = {"csv": Path(name).resolve().as_posix(), "windows": len(windows),
              "temperature": temperature, "bin_width": bin_width, "bootstrap": n_boot,
              "reweight": reweight, "min_overlap": round(min(ov), 4) if ov else None,
              "poor_overlap_after": poor}
    save_nested(config, section, ["pmf" + (f"_{reweight}" if reweight else "")], result)
    return result


def stitch_offset(mm_entry, qm_entry, temperature: float) -> Dict[str, Any]:
    """MM -> xTB/MM free energy at one window centre by BAR."""
    mm_mm, mm_qm = read_csv(mm_entry["energies_mm"]), read_csv(mm_entry["energies_xtb"])
    qm_mm, qm_qm = read_csv(qm_entry["energies_mm"]), read_csv(qm_entry["energies_xtb"])
    wf = mm_qm["energy_kcal"] - mm_mm["energy_kcal"]
    wr = qm_mm["energy_kcal"] - qm_qm["energy_kcal"]
    return {"z0": float(mm_entry["z0"]), "dA_kcal": bar(wf, wr, temperature),
            "n_mm": len(wf), "n_qmmm": len(wr)}


def combine_pmfs(z_mm, F_mm, z_qm, F_qm, z_lo: float, z_hi: float, dA_lo: float, dA_hi: float):
    """
    One profile: MM outside [z_lo, z_hi], the directly sampled QM/MM profile inside
    (aligned to MM at z_lo), and MM beyond z_hi shifted by dA_hi - dA_lo. The
    thermodynamic cycle says the QM/MM rise across the constriction equals the MM
    rise plus dA_hi - dA_lo; `closure` is how far the sampled profiles miss that,
    and shows as the step at z_hi.
    """
    z_mm, F_mm, z_qm, F_qm = (np.asarray(a, dtype=float) for a in (z_mm, F_mm, z_qm, F_qm))
    fm = lambda z: float(np.interp(z, z_mm, F_mm))
    fq = lambda z: float(np.interp(z, z_qm, F_qm))
    closure = (fq(z_hi) - fq(z_lo)) - (fm(z_hi) - fm(z_lo)) - (dA_hi - dA_lo)
    grid = np.union1d(z_mm, z_qm)
    out = []
    for z in grid:
        if z < z_lo:
            out.append(fm(z))
        elif z <= z_hi:
            out.append(fq(z) - fq(z_lo) + fm(z_lo))
        else:
            out.append(fm(z) + dA_hi - dA_lo)
    out = np.array(out)
    return grid, out - np.nanmin(out), closure


def run_stitch(config, qm_section: str, inputs_mm, inputs_qm, md) -> Dict[str, Any]:
    temperature = temperature_of(inputs_qm, md)
    centres = window_centres(inputs_qm)
    ends = []
    for z in (centres[0], centres[-1]):
        mm_entry, qm_entry = find_window(config, "umbrella", z), find_window(config, qm_section, z)
        if not mm_entry or not qm_entry:
            raise ValueError(f"stitch needs an MM and a QM/MM window at z0 = {z} A")
        for e, name in ((mm_entry, "umbrella"), (qm_entry, qm_section)):
            for level in ("mm", "xtb"):
                if not e.get(f"energies_{level}"):
                    raise ValueError(f"{name} window at z0 {z}: run `energies --level {level}` first")
        ends.append(stitch_offset(mm_entry, qm_entry, temperature))
    pmf_mm = (load_section(config, "umbrella").get("results") or {}).get("pmf")
    pmf_qm = (load_section(config, qm_section).get("results") or {}).get("pmf")
    if not pmf_mm or not pmf_qm:
        raise ValueError("run `wham` for both the MM and the QM/MM sections first")
    a, b = read_csv(pmf_mm["csv"]), read_csv(pmf_qm["csv"])
    ok_a, ok_b = np.isfinite(a["F_kcal"]), np.isfinite(b["F_kcal"])
    grid, F, closure = combine_pmfs(a["z_A"][ok_a], a["F_kcal"][ok_a], b["z_A"][ok_b],
                                    b["F_kcal"][ok_b], centres[0], centres[-1],
                                    ends[0]["dA_kcal"], ends[1]["dA_kcal"])
    name = f"stitched_{qm_section}.csv"
    with open(name, "w", newline="") as fh:
        out = csv.writer(fh)
        out.writerow(["z_A", "F_kcal"])
        for z, f in zip(grid, F):
            out.writerow([round(float(z), 5), round(float(f), 5)])
    result = {"csv": Path(name).resolve().as_posix(), "ends": ends,
              "closure_kcal": round(closure, 3)}
    print(f"[INFO] MM -> xTB/MM: dA {ends[0]['dA_kcal']:.2f} (z {ends[0]['z0']}) and "
          f"{ends[1]['dA_kcal']:.2f} kcal/mol (z {ends[1]['z0']}); cycle closure {closure:.2f} kcal/mol")
    save_nested(config, qm_section, ["stitch"], result)
    return result


# ==============================
# Entry point
# ==============================

def run(args) -> Dict[str, Any]:
    config = str(Path(args.config).resolve())
    section = getattr(args, "section", None) or "umbrella"
    inputs = umbrella_inputs(config, section)
    preprocess = load_section(config, "preprocess")
    md = load_section(config, "md")
    action = args.action
    index = getattr(args, "window", None)
    if action in ("window", "energies") and index is None:
        raise SystemExit(f"umbrella {action} needs --window N")

    if action == "pull":
        return run_pull(config, section, inputs, preprocess, md)
    if action == "window":
        if str(inputs["engine"]).lower() == "qmmm":
            return run_qmmm_window(config, section, inputs, preprocess, md, index)
        return run_mm_window(config, section, inputs, preprocess, md, index)
    if action == "energies":
        return run_energies(config, section, inputs, preprocess, md, index,
                            str(args.level).lower(), int(args.stride), args.qm_section)
    if action == "wham":
        return run_wham(config, section, inputs, md, getattr(args, "reweight", None),
                        int(getattr(args, "bootstrap", 50)))
    if action == "stitch":
        qm_section = section if section != "umbrella" else args.qm_section
        return run_stitch(config, qm_section, umbrella_inputs(config),
                          umbrella_inputs(config, qm_section), md)
    raise SystemExit(f"unknown umbrella action {action!r}")
