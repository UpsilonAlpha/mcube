#!/usr/bin/env python3
"""
Bulk-water check of an ion model, before any channel window: does the model give the
first hydration shell the literature measures? Desolvation is the step the channel PMF
turns on, and a model that holds its first shell too tightly or too loosely in bulk
will get the dehydration barrier wrong in the channel too.

Modes:

  run       MM in a cubic TIP3P box (NPT, 300 K) with mcube's metal parameters, then the
            shell analysis. Writes system.pdb (the final box, ion at the centre, residues
            numbered 1..N per chain, ready for `mcube preprocess` with skip_solvate), and with
            --dcd-every full-system frames.dcd (for the QM/MM swaps of docs/channel/qmmm).
  qmmm-config  an mcube config (preprocess + qmmm, QM/MM MD only) for that box, the QM
            region being the ion and the shell preprocess detects inside the MM first minimum.
  analyse   the same shell analysis on any trajectory, e.g. the xTB/MM `QM_MM.dcd` of
            `mcube qmmm` run on that box.

  python bulk_ion.py run --ion LA --set iod --ns 4 --out LA_iod
  python bulk_ion.py analyse --top LA_iod/system.pdb --traj qmmm/QM_MM.dcd --ion LA --out LA_xtb

Observables (report.json): ion-O RDF first peak and first minimum; coordination number
(waters with O inside the first minimum) as mean, spread and histogram; first-shell
exchange events and residence times, counting a water as having left only after it has
been outside for `--tstar` ps (2 ps, the usual Impey-Madden-McDonald tolerance).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

ION = {"CA": ("Ca", 2), "LA": ("La", 3), "LU": ("Lu", 3)}
R_MAX = 6.0      # A, RDF range
DR = 0.02        # A, RDF bin


# ==============================
# Shell analysis
# ==============================

def shell_distances(frames, box, ion: int, oxygens):
    """Ion-O distances (A) per frame, minimum image, shape (frames, waters)."""
    d = frames[:, oxygens, :] - frames[:, ion:ion + 1, :]
    d -= box[:, None, :] * np.round(d / box[:, None, :])
    return np.linalg.norm(d, axis=2)


def rdf(dist, volumes):
    """g(r) of ion-O from per-frame distances, normalised by the mean O density."""
    edges = np.arange(0.0, R_MAX + DR, DR)
    hist = np.histogram(dist.ravel(), bins=edges)[0].astype(float)
    r = 0.5 * (edges[1:] + edges[:-1])
    shell = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    density = dist.shape[1] / np.mean(volumes)
    return r, hist / (dist.shape[0] * density * shell)


def first_shell(r, g):
    """
    (peak, minimum) of the first shell in A.

    For a tightly bound ion g(r) is flat at (or near) zero between the first and second
    shells, so the "minimum" is a plateau, and its first bin sits at the tail of the peak.
    Taking that bin as the shell boundary (and so as the QM-region cutoff) can cut off
    first-shell waters, so the minimum is the middle of the first plateau within 5% of
    the peak height of the lowest value.
    """
    window = (r > 1.8) & (r < 3.4)
    peak = r[window][np.argmax(g[window])]
    g_peak = float(g[window].max())
    after = np.where((r > peak) & (r < peak + 1.6))[0]
    low = g[after].min()
    flat = after[g[after] <= low + 0.05 * g_peak]
    run = [flat[0]]
    for i in flat[1:]:
        if i != run[-1] + 1:
            break
        run.append(i)
    minimum = float(np.mean(r[run]))
    return float(peak), minimum


def residence(inside, dt_ps: float, tstar_ps: float):
    """
    First-shell residence: a water leaves only after `tstar` ps outside. Returns the
    number of completed stays (enter and leave inside the run), their mean length in ps,
    and the number of exchange events (departures).
    """
    nstar = max(1, int(round(tstar_ps / dt_ps)))
    stays, departures = [], 0
    for w in range(inside.shape[1]):
        col = inside[:, w]
        start, out_run, t = None, 0, 0
        for t, flag in enumerate(col):
            if flag:
                out_run = 0
                if start is None:
                    start = t
            elif start is not None:
                out_run += 1
                if out_run >= nstar:
                    end = t - out_run + 1
                    departures += 1
                    if start > 0:            # entered during the run: a complete stay
                        stays.append((end - start) * dt_ps)
                    start, out_run = None, 0
    return len(stays), (float(np.mean(stays)) if stays else None), departures


def analyse(frames, box, ion: int, oxygens, dt_ps: float, tstar_ps: float = 2.0) -> dict:
    if len(frames) < 2:
        raise ValueError(f"only {len(frames)} frame(s) to analyse: the trajectory is shorter than "
                         "the time discarded with --skip")
    dist = shell_distances(frames, box, ion, oxygens)
    r, g = rdf(dist, np.prod(box, axis=1))
    peak, minimum = first_shell(r, g)
    inside = dist < minimum
    cn = inside.sum(axis=1)
    n_stays, mean_stay, departures = residence(inside, dt_ps, tstar_ps)
    t_total = frames.shape[0] * dt_ps
    return {
        "frames": int(frames.shape[0]), "dt_ps": dt_ps, "time_ps": round(t_total, 2),
        "rdf_peak_A": round(peak, 3), "rdf_min_A": round(minimum, 3),
        "g_peak": round(float(g[np.argmin(np.abs(r - peak))]), 2),
        "cn_mean": round(float(cn.mean()), 3), "cn_std": round(float(cn.std()), 3),
        "cn_histogram": {int(k): int(v) for k, v in zip(*np.unique(cn, return_counts=True))},
        "exchanges": departures,
        "exchange_rate_per_ns": round(departures / (t_total / 1000.0), 3),
        "completed_stays": n_stays,
        "mean_residence_ps": round(mean_stay, 2) if mean_stay else None,
        "tstar_ps": tstar_ps,
        "rdf": {"r_A": [round(float(x), 3) for x in r], "g": [round(float(x), 4) for x in g]},
    }


# ==============================
# MM run
# ==============================

def threads() -> int:
    """mcube's rule (preprocess.resolve_cores) without importing ASH: $OMP_NUM_THREADS,
    else $SLURM_CPUS_PER_TASK, else 1. OpenMM's own default is every core on the node."""
    import os
    for var in ("OMP_NUM_THREADS", "SLURM_CPUS_PER_TASK"):
        value = (os.environ.get(var) or "").strip()
        if value.isdigit() and int(value) > 0:
            return int(value)
    return 1


def build_box(ion: str, box_nm: float, waters=None):
    import openmm
    from openmm import unit
    from openmm.app import ForceField, Modeller, Topology, element
    top = Topology()
    chain = top.addChain("A")
    top.addAtom(ion, element.Element.getBySymbol(ION[ion][0]), top.addResidue(ion, chain))
    centre = openmm.Vec3(box_nm / 2, box_nm / 2, box_nm / 2)
    modeller = Modeller(top, [centre] * unit.nanometer)
    ff = ForceField("amber14/tip3p.xml")
    if waters:
        # A fixed number of waters, as in a published protocol (cubic box chosen by Modeller)
        modeller.addSolvent(ff, model="tip3p", numAdded=int(waters), boxShape="cube",
                            positiveIon="Na+", negativeIon="Cl-", neutralize=True)
    else:
        modeller.addSolvent(ff, model="tip3p",
                            boxSize=openmm.Vec3(box_nm, box_nm, box_nm) * unit.nanometer,
                            positiveIon="Na+", negativeIon="Cl-", neutralize=True)
    return modeller, ff


def centred_positions(topology, pos_nm, box_nm, ion: int):
    """Positions with the ion at the box centre and every residue wrapped whole into the box."""
    pos = np.array(pos_nm) - np.array(pos_nm[ion]) + box_nm / 2
    for res in topology.residues():
        idx = [a.index for a in res.atoms()]
        shift = np.floor(pos[idx[0]] / box_nm) * box_nm
        pos[idx] -= shift
    return pos


def run(args):
    import openmm
    from openmm import unit
    from openmm.app import PME, HBonds, PDBFile, Simulation
    from mcube import ions

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    modeller, ff = build_box(args.ion, args.box, args.waters)
    top = modeller.topology
    system = ff.createSystem(top, nonbondedMethod=PME, nonbondedCutoff=1.0 * unit.nanometer,
                             constraints=HBonds, rigidWater=True)
    system.addForce(openmm.MonteCarloBarostat(1.0 * unit.bar, args.temp * unit.kelvin))
    changed = ions.apply_metal_parameters(system, top, [args.ion], args.set,
                                          forcefields=["amber14/tip3p.xml"])
    ion = changed[0]
    nb = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    q, sigma, eps = nb.getParticleParameters(ion)

    integrator = openmm.LangevinMiddleIntegrator(args.temp * unit.kelvin,
                                                 args.friction / unit.picosecond,
                                                 args.timestep * unit.picoseconds)
    import time
    platform = openmm.Platform.getPlatformByName(args.platform)
    props = {"Threads": str(threads())} if args.platform == "CPU" else {"Precision": "mixed"}
    sim = Simulation(top, system, integrator, platform, props)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=1000)
    sim.context.setVelocitiesToTemperature(args.temp * unit.kelvin)
    t0 = time.time()
    n_eq = int(args.equil / args.timestep)
    sim.step(n_eq)
    # Printed so a run's speed (and whether it will fit its time limit) shows in the log
    print(f"[INFO] {args.platform}: equilibration {args.equil} ps in {time.time() - t0:.0f} s, "
          f"{n_eq * args.timestep * 1e-3 / max(time.time() - t0, 1e-9) * 86400:.1f} ns/day", flush=True)

    oxygens = [a.index for a in top.atoms() if a.residue.name == "HOH" and a.name == "O"]
    stride = int(round(args.dt / args.timestep))
    n = int(round(args.ns * 1000.0 / args.dt))
    frames = np.empty((n, len(oxygens) + 1, 3), dtype=np.float32)
    boxes = np.empty((n, 3))
    keep = [ion] + oxygens
    if args.dcd_every:
        # Full-system frames for QM/MM vertical swaps (docs/channel/qmmm/exchange.py), from the
        # start of sampling; the topology is the system.pdb written at the end (same atom order)
        from openmm.app import DCDReporter
        sim.reporters.append(DCDReporter(str(out / "frames.dcd"),
                                         int(round(args.dcd_every / args.timestep))))
    t0 = time.time()
    for i in range(n):
        if i and i % max(1, n // 10) == 0:
            done = i * args.dt / 1000.0
            print(f"[INFO] {done:.2f}/{args.ns} ns, "
                  f"{done / max(time.time() - t0, 1e-9) * 86400:.1f} ns/day", flush=True)
        sim.step(stride)
        state = sim.context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        frames[i] = pos[keep]
        v = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom)
        boxes[i] = [v[0][0], v[1][1], v[2][2]]

    report = analyse(frames, boxes, 0, list(range(1, len(keep))), args.dt, args.tstar)
    report.update({"mode": "mm", "ion": args.ion, "set": args.set, "temp_K": args.temp,
                   "timestep_ps": args.timestep, "friction_per_ps": args.friction,
                   "dcd_every_ps": args.dcd_every,
                   "atoms": top.getNumAtoms(), "waters": len(oxygens),
                   "sigma_nm": round(sigma.value_in_unit(unit.nanometer), 5),
                   "epsilon_kJ": round(eps.value_in_unit(unit.kilojoule_per_mole), 5),
                   "c4_kcal_A4": ions.c4_value(ION[args.ion][0], ION[args.ion][1])
                   if str(args.set) == "1264" else None})

    state = sim.context.getState(getPositions=True)
    box_nm = np.array([state.getPeriodicBoxVectors()[k][k].value_in_unit(unit.nanometer)
                       for k in range(3)])
    pos = centred_positions(top, state.getPositions(asNumpy=True).value_in_unit(unit.nanometer),
                            box_nm, ion)
    top.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
    with open(out / "system.pdb", "w") as fh:
        PDBFile.writeFile(top, pos * unit.nanometer, fh)     # residues renumbered 1..N per chain
    write_report(report, out)


def analyse_traj(args):
    import mdtraj
    # several files are joined in the order given (a resumed qmmm run: QM_MM_part1.dcd, ..., QM_MM.dcd)
    traj = mdtraj.load(args.traj, top=args.top)
    t = traj.topology
    ion = next(a.index for a in t.atoms if a.residue.name == args.ion)
    oxygens = [a.index for a in t.atoms if a.residue.name == "HOH" and a.name == "O"]
    frames = traj.xyz[:, [ion] + oxygens, :] * 10.0
    boxes = traj.unitcell_lengths * 10.0
    dt = float(args.dt) if args.dt else float(traj.timestep)
    report = analyse(frames, boxes, 0, list(range(1, len(oxygens) + 1)), dt, args.tstar)
    report.update({"mode": "trajectory", "ion": args.ion,
                   "traj": [str(Path(p).resolve()) for p in args.traj],
                   "skip_ps": args.skip})
    if args.skip:
        k = int(round(args.skip / dt))
        report = {**analyse(frames[k:], boxes[k:], 0, list(range(1, len(oxygens) + 1)), dt,
                            args.tstar), **{k2: v for k2, v in report.items()
                                            if k2 in ("mode", "ion", "traj", "skip_ps")}}
    write_report(report, Path(args.out))


def qmmm_config(args):
    """
    An mcube config for xTB/MM MD of the MM box: preprocess (skip_solvate) detects the
    hydration shell inside the MM first minimum and makes ion + shell the QM region; qmmm
    runs QM/MM MD only. The same metal parameters as the MM run reach the QM/MM system.
    """
    import yaml
    run_dir = Path(args.run).resolve()
    mm = json.loads((run_dir / "report.json").read_text())
    pdb = (run_dir / "system.pdb").as_posix()
    config = {
        "preprocess": {"inputs": {
            "protein_pdb": pdb, "skip_solvate": True, "metals": [args.ion],
            "cutoff": round(mm["rdf_min_A"], 2),
            "default_forcefields": ["amber14/tip3p.xml"]}},
        "md": {"inputs": {"metal_parameters": str(mm["set"])}},
        "qmmm": {"inputs": {
            "pdb_file": pdb, "forcefields": ["amber14/tip3p.xml"],
            "charge": ION[args.ion][1], "mult": 1, "xtb_method": "GFN2",
            "cores": args.cores, "temp": mm["temp_K"], "timestep": 0.001, "time": args.ps,
            "traj_frequency": args.every, "skip_optimize": True,
            # a timed-out QM/MM MD continues from ASH's state file on resubmission (time is
            # the total); with no state file it starts fresh, so it is safe on the first run
            "resume": True,
            "metal_parameters": str(mm["set"])}},
    }
    for key, value in (("xtb_etemp", args.etemp), ("xtb_pc_hardness", args.pc_hardness)):
        if value is not None:
            config["qmmm"]["inputs"][key] = value
    if args.shell2:
        # Second shell in the QM region too: preprocess takes every water inside this cutoff
        config["preprocess"]["inputs"]["cutoff"] = args.shell2
    out = run_dir / args.dir / "config.yaml"
    out.parent.mkdir(exist_ok=True)
    out.write_text(yaml.safe_dump(config, sort_keys=False))
    print(out)


def write_report(report: dict, out: Path):
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "report.json", "w") as fh:
        json.dump(report, fh, indent=1)
    brief = {k: v for k, v in report.items() if k != "rdf"}
    print(json.dumps(brief, indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="mode", required=True)
    r = sub.add_parser("run")
    r.add_argument("--ion", required=True, choices=sorted(ION))
    r.add_argument("--set", default="iod")
    r.add_argument("--ns", type=float, default=4.0, help="production length")
    r.add_argument("--equil", type=float, default=200.0, help="ps of NPT before sampling")
    r.add_argument("--dt", type=float, default=0.5, help="ps between samples")
    r.add_argument("--box", type=float, default=3.5, help="nm, cubic")
    r.add_argument("--temp", type=float, default=300.0)
    r.add_argument("--tstar", type=float, default=2.0)
    r.add_argument("--platform", default="CPU", help="OpenMM platform (CPU, or CUDA on gpu_cuda)")
    r.add_argument("--timestep", type=float, default=0.002, help="ps")
    r.add_argument("--friction", type=float, default=1.0, help="Langevin friction, 1/ps")
    r.add_argument("--waters", type=int, help="solvate with this many waters instead of --box")
    r.add_argument("--dcd-every", type=float, help="ps between full-system frames in frames.dcd "
                   "(for docs/channel/qmmm/exchange.py); none by default")
    r.add_argument("--out", required=True)
    a = sub.add_parser("analyse")
    a.add_argument("--top", required=True)
    a.add_argument("--traj", required=True, nargs="+", help="one or more DCDs, joined in order")
    a.add_argument("--ion", required=True, choices=sorted(ION))
    a.add_argument("--dt", type=float, help="ps between frames (default: from the file)")
    a.add_argument("--skip", type=float, default=0.0, help="ps to discard at the start")
    a.add_argument("--tstar", type=float, default=2.0)
    a.add_argument("--out", required=True)
    q = sub.add_parser("qmmm-config")
    q.add_argument("--run", required=True, help="a `run` output directory")
    q.add_argument("--ion", required=True, choices=sorted(ION))
    q.add_argument("--ps", type=float, default=10.0, help="QM/MM MD length")
    q.add_argument("--every", type=int, default=20, help="steps between saved frames")
    q.add_argument("--cores", type=int, default=4)
    q.add_argument("--etemp", type=float, help="xTB electronic temperature (K)")
    q.add_argument("--pc-hardness", type=lambda v: v if v == "elements" else float(v),
                   help="xTB hardness of the MM point charges: a number (ASH default 1000, bare "
                        "charges) or 'elements' (xTB's per-element hardness, smeared charges)")
    q.add_argument("--shell2", type=float, help="QM cutoff (A) taking in the second shell")
    q.add_argument("--dir", default="qmmm", help="subdirectory of --run for this QM/MM run")
    args = p.parse_args()
    {"run": run, "analyse": analyse_traj, "qmmm-config": qmmm_config}[args.mode](args)


if __name__ == "__main__":
    main()
