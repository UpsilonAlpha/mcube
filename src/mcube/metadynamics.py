#!/usr/bin/env python3
"""
Well-tempered metadynamics of one ion (an MM model, the umbrella stage's system) on its
carboxylate-oxygen coordination number n_COO, optionally together with z.

n_COO is sum_j 1/(1 + (r_ij/r0)^12) over the Asp OD1/OD2 and Glu OE1/OE2 oxygens of
charged Asp and Glu (ASH/GLH and ASPP/GLUP are the neutral acids and are left out), with
r0 = `cn_cutoff`: the same switching function as umbrella's logged coordination numbers.
z is umbrella's: the ion minus the mass-weighted centre of `reference`, in A. Two modes:

    cvs: [n_coo]       n_COO biased; z held at `z0` by umbrella's harmonic restraint (`k`)
                       and the ion over the pore by its cylinder: metadynamics inside an
                       umbrella window, for a carboxylate that never lets go of the ion
    cvs: [z, n_coo]    both biased in 2D on `grid.z`; the cylinder still applies, the z
                       restraint is off, and `z_walls` adds optional flat-bottomed walls

    mcube metadynamics config.yaml run [--section S]    the run (OpenMM, GPU)
    mcube metadynamics config.yaml fes [--section S]    free energies from the hills

Inputs are `metadynamics.inputs` over the keys it shares with `umbrella.inputs` (ion,
reference, cylinder, metal parameters ...: see SHARED), overridden by the `--section`
inputs, over DEFAULTS. A run writes `<section>/`: cv.csv (time, z, r, n_COO, cn_water,
cn_protein, bias energy), hills.csv (every Gaussian: time, centre, height), OpenMM's
bias grid in bias/, lastframe.pdb and, with `traj_frequency` > 0, run.dcd. Running `run`
again extends the run: it starts from lastframe, OpenMM reloads the bias from bias/, and
cv.csv and hills.csv are appended to (velocities are redrawn).

The run is openmm.app.Metadynamics on the umbrella system, not ASH's OpenMM_metadynamics.
ASH's native back end is a thin wrapper around the same class: its coordination-number CV
is a different switching function (1/(1+x^6)) and its "custom" CV takes a ready-made
BiasVariable, so it would add nothing but an ASH OpenMMTheory in place of the md stage's
system (membrane barostat, HMR, 12-6-4 / C4 pairs) and one step() call for the whole run,
with no way to log the CVs along it. PLUMED is not needed.

`fes` rebuilds the bias from hills.csv on the grid (the same Gaussians OpenMM adds), and
the well-tempered free energy F = -(gamma/(gamma-1)) V, projected onto each biased CV
(and F(z) in 2D). A smooth n_COO near the cutoff is not a state, so n_COO is reported as
integer states: F_n = -kT ln of the integral of exp(-F/kT) over [n-0.5, n+0.5]. As a
convergence check the state free energies are also given from the hills of the first
50 % and 75 % of the run.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import yaml

from .preprocess import load_section
from .umbrella import (
    KB,
    TO_KJ_NM2,
    atom_masses,
    bias_force,
    carboxylate_oxygens,
    centroid_cv,
    cn_atom_groups,
    cn_cutoff,
    cn_cv,
    cv_now,
    find_window,
    image_ion,
    ion_atom,
    load_structure,
    make_simulation,
    mm_system,
    pdb_box,
    read_csv,
    reference_atoms,
    save_nested,
    temperature_of,
    write_coords_pdb,
    write_state_pdb,
)

KJ_TO_KCAL = 1.0 / 4.184
WALL_NAME = "mcube_mtd_walls"
CVS = ("n_coo", "z")

DEFAULTS = {
    "cvs": ["n_coo"],
    "k": 10.0,                # kcal/mol/A^2, the z restraint when z is not biased
    "height": 1.0,            # kJ/mol, initial Gaussian height
    "sigma": {"n_coo": 0.1, "z": 0.3},              # Gaussian widths (z in A)
    "grid": {"n_coo": {"min": 0.0, "max": 6.0}},    # bins: OpenMM's default 5 (max-min)/sigma
    "biasfactor": 10.0,
    "frequency": 500,         # steps between hills
    "time": 10000.0,          # ps
    "timestep": 0.004,        # ps
    "cv_frequency": 250,      # steps
    "traj_frequency": 25000,  # steps; 0 writes no DCD
    "platform": "CPU",
}

# umbrella.inputs keys that mean the same here, so one config states them once
SHARED = ("ion", "reference", "cylinder", "k", "timestep", "platform", "properties", "temp",
          "cores", "metal_parameters", "polarizabilities", "c4_set", "c4_pairs", "cn_cutoff",
          "nonbonded_cutoff")


# ==============================
# Config
# ==============================

def metadynamics_inputs(config: str, section: Optional[str] = None) -> Dict[str, Any]:
    """DEFAULTS < umbrella.inputs (SHARED keys) < metadynamics.inputs < `section`'s inputs.
    `sigma` and `grid` merge per CV, so a section can change one CV's settings only."""
    with open(config) as fh:
        data = yaml.safe_load(fh) or {}
    layers = [{k: v for k, v in ((data.get("umbrella") or {}).get("inputs") or {}).items()
               if k in SHARED},
              (data.get("metadynamics") or {}).get("inputs") or {}]
    if section and section != "metadynamics":
        layers.append((data.get(section) or {}).get("inputs") or {})
    inputs = {k: (dict(v) if isinstance(v, dict) else v) for k, v in DEFAULTS.items()}
    inputs["grid"] = {k: dict(v) for k, v in DEFAULTS["grid"].items()}
    for layer in layers:
        for key, value in layer.items():
            if key == "sigma" and isinstance(value, dict):
                inputs["sigma"].update(value)
            elif key == "grid" and isinstance(value, dict):
                for cv, g in value.items():
                    inputs["grid"].setdefault(cv, {}).update(g or {})
            else:
                inputs[key] = value
    return inputs


def cv_list(inputs) -> List[str]:
    cvs = [str(c).lower() for c in inputs["cvs"]]
    if not cvs or len(set(cvs)) != len(cvs) or any(c not in CVS for c in cvs):
        raise ValueError(f"metadynamics cvs must be distinct names from {list(CVS)}, got {inputs['cvs']}")
    # `coordinate` is umbrella's key and it is a shared input, so a section that sets it
    # to `distance` would otherwise reach here and be biased along z without comment.
    if str(inputs.get("coordinate", "z")).lower() != "z":
        raise ValueError("metadynamics biases the CVs in `cvs` and only supports the channel "
                         f"coordinate; `coordinate: {inputs['coordinate']}` is an umbrella key. "
                         "For a distance CV use umbrella windows (`coordinate: distance`).")
    # `z0`, not `cv0`: this stage is z-only, and `run_seed` and the restraint centre read
    # `z0`. Accepting `cv0` here let a section pass validation and then silently hold z at
    # 0.0 A and seed from md's lastframe, which is worse than being told to rename a key.
    if "z" not in cvs and inputs.get("z0") is None:
        raise ValueError("metadynamics without z biased holds z at z0: give `z0`")
    for c in cvs:
        g = (inputs.get("grid") or {}).get(c) or {}
        if g.get("min") is None or g.get("max") is None or float(g["max"]) <= float(g["min"]):
            raise ValueError(f"metadynamics needs grid.{c}.min < grid.{c}.max")
        if float((inputs.get("sigma") or {}).get(c) or 0) <= 0:
            raise ValueError(f"metadynamics needs sigma.{c} > 0")
    return cvs


def grid_axis(inputs, cv: str) -> np.ndarray:
    """Grid points of one CV in its config units (A for z), as OpenMM tabulates the bias:
    `bins` points from min to max inclusive, by default ceil(5 (max - min) / sigma)."""
    g = inputs["grid"][cv]
    lo, hi = float(g["min"]), float(g["max"])
    bins = g.get("bins") or int(np.ceil(5 * (hi - lo) / float(inputs["sigma"][cv])))
    return np.linspace(lo, hi, int(bins))


def unit_scale(cv: str) -> float:
    """Config units per OpenMM unit (z: A per nm; n_COO is a number)."""
    return 10.0 if cv == "z" else 1.0


def run_dir(section: str, inputs) -> Path:
    tag = inputs.get("tag")
    return Path(f"{section}_{tag}" if tag else section)


# ==============================
# Forces
# ==============================

def logged_metadynamics(*args, **kwargs):
    """openmm.app.Metadynamics that also records every Gaussian it adds in `.hills`
    as (time from `.clock()`, CV values in OpenMM units, height in kJ/mol), so the bias
    at any point of the run can be rebuilt. It hooks the private _addGaussian, the one
    place the hill's centre and well-tempered height are both known."""
    from openmm import unit
    from openmm.app import Metadynamics

    class Logged(Metadynamics):
        def _addGaussian(self, position, height, context):
            self.hills.append((self.clock(), [float(p) for p in position],
                               height.value_in_unit(unit.kilojoules_per_mole)))
            super()._addGaussian(position, height, context)

    meta = Logged(*args, **kwargs)
    meta.hills, meta.clock = [], lambda: 0.0
    return meta


def bias_variables(inputs, cvs, ion: int, reference, masses, coo: List[int], cutoff: float):
    """One BiasVariable per CV, on its own CV force, with the grid of grid_axis."""
    from openmm.app import BiasVariable
    out = []
    for cv in cvs:
        force = (centroid_cv("z2-z1", ion, reference, masses) if cv == "z"
                 else cn_cv(ion, coo, cutoff))
        axis = grid_axis(inputs, cv) / unit_scale(cv)
        out.append(BiasVariable(force, float(axis[0]), float(axis[-1]),
                                float(inputs["sigma"][cv]) / unit_scale(cv), False,
                                gridWidth=len(axis)))
    return out


def wall_force(ion: int, reference, masses, walls):
    """Flat-bottomed walls on z: 0.5 k (max(0, z - max)^2 + max(0, min - z)^2), A and kcal/mol/A^2."""
    import openmm
    force = openmm.CustomCVForce("0.5*mcube_wall_k*(max(0, zw-mcube_wall_hi)^2"
                                 " + max(0, mcube_wall_lo-zw)^2)")
    force.addCollectiveVariable("zw", centroid_cv("z2-z1", ion, reference, masses))
    force.addGlobalParameter("mcube_wall_k", float(walls.get("k", 10.0)) * TO_KJ_NM2)
    force.addGlobalParameter("mcube_wall_lo", float(walls["min"]) / 10.0)
    force.addGlobalParameter("mcube_wall_hi", float(walls["max"]) / 10.0)
    force.setName(WALL_NAME)
    return force


# ==============================
# Run
# ==============================

def run_seed(config: str, inputs, md) -> str:
    """`seed`, else the lastframe of the umbrella window at z0 (section `seed_section`,
    default umbrella), else md's lastframe."""
    if inputs.get("seed"):
        return inputs["seed"]
    if inputs.get("z0") is not None:
        window = find_window(config, inputs.get("seed_section", "umbrella"), float(inputs["z0"]))
        if window and window.get("lastframe"):
            return window["lastframe"]
    last = (md.get("results") or {}).get("lastframe")
    if not last:
        raise ValueError("metadynamics needs a starting structure: `seed`, an umbrella window at z0, "
                         "or md.results.lastframe")
    return last


def run_metadynamics(config, section, inputs, preprocess, md) -> Dict[str, Any]:
    from openmm.app import DCDReporter
    cvs = cv_list(inputs)
    previous = (load_section(config, section).get("results") or {}).get("run") or {}
    work = run_dir(section, inputs)
    hills_file = (work / "hills.csv").resolve()
    extend = (bool(previous.get("lastframe")) and hills_file.exists()
              and previous.get("hills") == hills_file.as_posix())
    seed = previous["lastframe"] if extend else run_seed(config, inputs, md)
    seed = Path(seed).resolve().as_posix()
    (work / "bias").mkdir(parents=True, exist_ok=True)
    os.chdir(work)

    pdb = load_structure(seed, preprocess)
    box = pdb_box(pdb)
    ion = ion_atom(pdb.topology, inputs["ion"])
    reference = reference_atoms(pdb.topology, inputs["reference"])
    masses = atom_masses(pdb.topology)
    coords = image_ion(np.array(pdb.getPositions()._value) * 10.0, ion, reference, masses, box)
    start = write_coords_pdb(pdb.topology, coords, "start.pdb")
    coo = carboxylate_oxygens(pdb.topology, ion)
    if "n_coo" in cvs and not coo:
        raise ValueError("n_coo: the structure has no charged Asp/Glu carboxylate oxygens")
    cutoff = cn_cutoff(inputs, preprocess)

    system, _ = mm_system(pdb, inputs, preprocess, md, config)
    # The umbrella force holds z (off when z is biased) and the cylinder, and logs the CVs
    groups = dict(cn_atom_groups(pdb.topology, ion), coo=coo)
    k = 0.0 if "z" in cvs else float(inputs["k"])
    z0 = float(inputs.get("z0") or 0.0)
    logger = bias_force(ion, reference, k, z0, inputs.get("cylinder"), masses, groups, cutoff)
    system.addForce(logger)
    if "z" in cvs and inputs.get("z_walls"):
        system.addForce(wall_force(ion, reference, masses, inputs["z_walls"]))
    temp = temperature_of(inputs, md)
    freq = int(inputs["frequency"])
    save = int(inputs.get("save_frequency") or 10 * freq)
    meta = logged_metadynamics(system, bias_variables(inputs, cvs, ion, reference, masses, coo, cutoff),
                        temp, float(inputs["biasfactor"]), float(inputs["height"]), freq,
                        saveFrequency=save, biasDir=str(Path("bias").resolve()))
    simulation, _ = make_simulation(pdb, system, coords, inputs, md)
    if not extend:
        simulation.minimizeEnergy(maxIterations=500)
    simulation.context.setVelocitiesToTemperature(temp)

    dt, chunk = float(inputs["timestep"]), int(inputs["cv_frequency"])
    t0 = float(previous.get("time_ps", 0.0)) if extend else 0.0
    meta.clock = lambda: t0 + simulation.currentStep * dt
    n_prod = int(round(float(inputs["time"]) / dt))
    traj = int(inputs.get("traj_frequency") or 0)
    if traj:
        simulation.reporters.append(DCDReporter("run.dcd", traj, append=extend and Path("run.dcd").exists()))
    group = meta._force.getForceGroup()
    print(f"[INFO] Metadynamics {section}: {'extending from ' + str(t0) + ' ps, ' if extend else ''}"
          f"CVs {cvs}, {inputs['time']} ps, hills of {inputs['height']} kJ/mol every {freq} steps, "
          f"bias factor {inputs['biasfactor']}, "
          + (f"z free (walls {inputs.get('z_walls')})" if "z" in cvs else f"z held at {z0} A (k {k})")
          + f"; {len(coo)} carboxylate O")
    names = ["z_A", "r_A", "cn_water", "cn_protein", "n_coo"]
    mode = "a" if extend else "w"
    with open("cv.csv", mode, newline="") as fh, open("hills.csv", mode, newline="") as fh_h:
        out, hills = csv.writer(fh), csv.writer(fh_h)
        if not extend:
            out.writerow(["time_ps"] + names + ["bias_kJ"])
            hills.writerow(["time_ps"] + [f"{c}_A" if c == "z" else c for c in cvs] + ["height_kJ"])
        written = 0
        for _ in range(n_prod // chunk):
            meta.step(simulation, chunk)
            values = cv_now(logger, simulation.context)
            energy = simulation.context.getState(getEnergy=True, groups={group}).getPotentialEnergy()
            out.writerow([round(meta.clock(), 4)] + [round(float(v), 4) for v in values]
                         + [round(energy._value, 4)])
            for t, pos, h in meta.hills[written:]:
                hills.writerow([round(t, 4)] + [round(p * unit_scale(c), 5) for p, c in zip(pos, cvs)]
                               + [h])
            written = len(meta.hills)
    last = write_state_pdb(simulation, pdb.topology, "lastframe.pdb")

    cv = read_csv("cv.csv")
    entry = {"cvs": cvs, "seed": seed, "start": start, "z0": inputs.get("z0"),
             "k": k, "cylinder": inputs.get("cylinder"), "z_walls": inputs.get("z_walls"),
             "height": float(inputs["height"]), "sigma": {c: float(inputs["sigma"][c]) for c in cvs},
             "grid": {c: {"min": float(grid_axis(inputs, c)[0]), "max": float(grid_axis(inputs, c)[-1]),
                          "bins": len(grid_axis(inputs, c))} for c in cvs},
             "biasfactor": float(inputs["biasfactor"]), "frequency": freq, "temperature": temp,
             "cn_cutoff": cutoff, "carboxylate_oxygens": len(coo),
             "cv": Path("cv.csv").resolve().as_posix(), "hills": Path("hills.csv").resolve().as_posix(),
             "bias_dir": Path("bias").resolve().as_posix(),
             "dcd": Path("run.dcd").resolve().as_posix() if traj else None,
             "lastframe": last, "timestep": dt, "time_ps": round(t0 + n_prod // chunk * chunk * dt, 4),
             "runs": int(previous.get("runs", 0)) + 1 if extend else 1,
             "n_coo_mean": round(float(cv["n_coo"].mean()), 3) if len(cv.get("n_coo", [])) else None,
             "done": True}
    save_nested(config, section, ["run"], entry)
    return entry


# ==============================
# Free energy
# ==============================

def bias_on_grid(centres, heights, sigmas, axes) -> np.ndarray:
    """Sum of Gaussians h exp(-0.5 sum_d ((x_d - c_d)/s_d)^2) on the grid `axes` (one array
    per CV, in the order of the CVs), as OpenMM adds them. Shape: one axis per CV."""
    centres = np.atleast_2d(np.asarray(centres, dtype=float))
    heights = np.asarray(heights, dtype=float)
    V = np.zeros(tuple(len(a) for a in axes))
    if not len(heights):
        return V
    gauss = [np.exp(-0.5 * ((np.asarray(a)[None, :] - centres[:, d][:, None]) / s) ** 2)
             for d, (a, s) in enumerate(zip(axes, sigmas))]
    if len(axes) == 1:
        return heights @ gauss[0]
    return np.einsum("h,hi,hj->ij", heights, gauss[0], gauss[1])


def well_tempered_fes(V, biasfactor: float) -> np.ndarray:
    """F = -(gamma/(gamma-1)) V, shifted to min 0 (same units as V)."""
    F = -(float(biasfactor) / (float(biasfactor) - 1.0)) * np.asarray(V, dtype=float)
    return F - F.min()


def project(F, keep: int, kT: float) -> np.ndarray:
    """F along axis `keep`: -kT ln sum over the other axes of exp(-F/kT), min 0."""
    F = np.asarray(F, dtype=float)
    if F.ndim == 1:
        return F - F.min()
    other = tuple(i for i in range(F.ndim) if i != keep)
    m = F.min()
    P = np.exp(-(F - m) / kT).sum(axis=other)
    out = -kT * np.log(P)
    return out - out.min()


def integer_states(x, F, kT: float) -> Dict[int, float]:
    """F_n = -kT ln sum over grid points with x in [n - 0.5, n + 0.5) of exp(-F/kT) dx,
    for each integer n on the grid, relative to the lowest state."""
    x, F = np.asarray(x, dtype=float), np.asarray(F, dtype=float)
    dx = np.gradient(x) if len(x) > 1 else np.ones_like(x)
    w = np.exp(-(F - F.min()) / kT) * dx
    states = {}
    for n in range(int(np.ceil(x.min() - 0.5)), int(np.floor(x.max() + 0.5)) + 1):
        inside = (x >= n - 0.5) & (x < n + 0.5)
        if inside.any() and w[inside].sum() > 0:
            states[n] = -kT * np.log(w[inside].sum())
    low = min(states.values())
    return {n: round(float(f - low), 3) for n, f in states.items()}


def fes_from_hills(hills, cvs, inputs, temperature: float, fraction: float = 1.0):
    """(axes, F in kcal/mol on the grid) from the hills of the first `fraction` of the run."""
    t = hills["time_ps"]
    keep = t <= t.max() * fraction + 1e-9 if len(t) else np.zeros(0, bool)
    cols = [f"{c}_A" if c == "z" else c for c in cvs]
    centres = np.column_stack([hills[c][keep] for c in cols]) if keep.any() else np.zeros((0, len(cvs)))
    axes = [grid_axis(inputs, c) for c in cvs]
    V = bias_on_grid(centres, hills["height_kJ"][keep] * KJ_TO_KCAL,
                     [float(inputs["sigma"][c]) for c in cvs], axes)
    return axes, well_tempered_fes(V, float(inputs["biasfactor"]))


def run_fes(config, section, inputs, md) -> Dict[str, Any]:
    run = (load_section(config, section).get("results") or {}).get("run")
    if not run or not run.get("hills"):
        raise ValueError(f"{section}: run `metadynamics run` first")
    cvs = run["cvs"]
    # the grid and widths the run used, so a later config edit cannot change the analysis
    inputs = dict(inputs, cvs=cvs, sigma=run["sigma"], grid=run["grid"], biasfactor=run["biasfactor"])
    temperature = float(run.get("temperature") or temperature_of(inputs, md))
    kT = KB * temperature
    hills = read_csv(run["hills"])
    if not hills:
        raise ValueError(f"{section}: no hills in {run['hills']}")
    os.chdir(Path(run["hills"]).parent)
    result: Dict[str, Any] = {"temperature": temperature, "hills": len(hills["time_ps"]),
                              "time_ps": run.get("time_ps")}
    convergence = {}
    for fraction in (0.5, 0.75, 1.0):
        axes, F = fes_from_hills(hills, cvs, inputs, temperature, fraction)
        if fraction == 1.0:
            for d, cv in enumerate(cvs):
                name = f"fes_{cv}.csv"
                with open(name, "w", newline="") as fh:
                    out = csv.writer(fh)
                    out.writerow([f"{cv}_A" if cv == "z" else cv, "F_kcal"])
                    for x, f in zip(axes[d], project(F, d, kT)):
                        out.writerow([round(float(x), 5), round(float(f), 5)])
                result[f"fes_{cv}"] = Path(name).resolve().as_posix()
            if len(cvs) == 2:
                with open("fes_2d.csv", "w", newline="") as fh:
                    out = csv.writer(fh)
                    out.writerow([f"{c}_A" if c == "z" else c for c in cvs] + ["F_kcal"])
                    for i, a in enumerate(axes[0]):
                        for j, b in enumerate(axes[1]):
                            out.writerow([round(float(a), 5), round(float(b), 5), round(float(F[i, j]), 5)])
                result["fes_2d"] = Path("fes_2d.csv").resolve().as_posix()
        if "n_coo" in cvs:
            d = cvs.index("n_coo")
            convergence[f"{int(fraction * 100)}%"] = integer_states(axes[d], project(F, d, kT), kT)
    if convergence:
        states = convergence["100%"]
        result["n_coo_states_kcal"] = states
        result["convergence"] = convergence
        drift = [abs(convergence["75%"].get(n, np.nan) - f) for n, f in states.items()]
        drift = [x for x in drift if np.isfinite(x)]
        result["max_state_change_75_100_kcal"] = round(max(drift), 3) if drift else None
        with open("n_coo_states.csv", "w", newline="") as fh:
            out = csv.writer(fh)
            out.writerow(["n_coo", "F_50_kcal", "F_75_kcal", "F_kcal"])
            for n in sorted(states):
                out.writerow([n] + [convergence[p].get(n, "") for p in ("50%", "75%", "100%")])
        result["n_coo_states"] = Path("n_coo_states.csv").resolve().as_posix()
        print(f"[INFO] n_COO states (kcal/mol, from all hills): {states}; "
              f"largest change since 75 %: {result['max_state_change_75_100_kcal']}")
    save_nested(config, section, ["fes"], result)
    return result


# ==============================
# Entry point
# ==============================

def run(args) -> Dict[str, Any]:
    config = str(Path(args.config).resolve())
    section = getattr(args, "section", None) or "metadynamics"
    inputs = metadynamics_inputs(config, section)
    md = load_section(config, "md")
    if args.action == "run":
        return run_metadynamics(config, section, inputs, load_section(config, "preprocess"), md)
    if args.action == "fes":
        return run_fes(config, section, inputs, md)
    raise SystemExit(f"unknown metadynamics action {args.action!r}")
