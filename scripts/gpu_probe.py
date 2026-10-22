#!/usr/bin/env python3
"""Can this environment actually run OpenMM on a GPU, and how fast?

    python gpu_probe.py                              platforms, context on each
    python gpu_probe.py system.pdb [extra.xml ...]   the above, plus a benchmark

Any force field files after the structure are added to the amber14 pair, which
is what a system holding a ligand or a non-standard residue needs - without
them `createSystem` fails on a missing template and the benchmark says nothing
about the GPU.

Run it *on a GPU node*, inside a job. Listing platforms on a login or
interactive node only says what OpenMM was compiled with; whether a platform
initialises depends on the driver being present, which is a different question
and the one that matters.

Note what the platform list does and does not tell you. OpenMM reaches AMD
cards through HIP or OpenCL, so a build reporting only Reference/CPU/CUDA cannot
use a ROCm card at all, whatever the queue offers. A build reporting CUDA can
still fail to create a context - no driver, a MIG slice not visible, or a card
already busy - which is why this creates one rather than trusting the list.

The benchmark uses the same settings as mcube's `md` stage (PME, 12 A cutoff,
HBond constraints, rigid water, 1.5 amu hydrogens, 4 fs) so the ns/day it
reports is comparable with a real run rather than with a toy.
"""

from __future__ import annotations

import sys
import time

import openmm
import openmm.app as app
import openmm.unit as unit

TIMESTEP_FS = 4.0
BENCH_STEPS = 2000
FORCEFIELDS = ("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")


def platforms() -> list:
    return [openmm.Platform.getPlatform(i).getName()
            for i in range(openmm.Platform.getNumPlatforms())]


def context_on(name: str):
    """(device_description, None) if a context opens on `name`, else (None, why)."""
    system = openmm.System()
    for _ in range(2):
        system.addParticle(1.0 * unit.amu)
    force = openmm.NonbondedForce()
    for _ in range(2):
        force.addParticle(0.0, 1.0 * unit.nanometer, 0.0)
    system.addForce(force)
    try:
        platform = openmm.Platform.getPlatformByName(name)
        context = openmm.Context(system, openmm.VerletIntegrator(0.001), platform)
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}"
    detail = []
    for prop in ("DeviceName", "DeviceIndex", "Precision", "UseBlockingSync"):
        try:
            detail.append(f"{prop}={platform.getPropertyValue(context, prop)}")
        except Exception:
            pass
    del context
    return (", ".join(detail) or "context opened"), None


def benchmark(pdb_file: str, name: str, extra=()) -> str:
    """ns/day for `pdb_file` on platform `name`, at mcube's md settings."""
    pdb = app.PDBFile(pdb_file)
    field = app.ForceField(*FORCEFIELDS, *extra)
    system = field.createSystem(pdb.topology, nonbondedMethod=app.PME,
                                nonbondedCutoff=1.2 * unit.nanometer,
                                constraints=app.HBonds, rigidWater=True,
                                hydrogenMass=1.5 * unit.amu)
    integrator = openmm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, TIMESTEP_FS * unit.femtosecond)
    sim = app.Simulation(pdb.topology, system, integrator,
                         openmm.Platform.getPlatformByName(name))
    sim.context.setPositions(pdb.positions)
    sim.minimizeEnergy(maxIterations=50)
    sim.step(100)                                  # warm the kernels, then time
    start = time.perf_counter()
    sim.step(BENCH_STEPS)
    seconds = time.perf_counter() - start
    ns = BENCH_STEPS * TIMESTEP_FS * 1e-6
    per_day = ns * 86400 / seconds
    return (f"{system.getNumParticles()} atoms, {BENCH_STEPS} steps in "
            f"{seconds:.1f} s -> {per_day:.1f} ns/day")


def main(argv) -> int:
    print("openmm", openmm.version.version)
    names = platforms()
    print("platforms:", names)
    for absent in ("HIP", "OpenCL"):
        if absent not in names:
            print(f"[INFO] no {absent} platform: AMD/ROCm cards are unusable by this build")

    usable = []
    for name in names:
        if name in ("Reference", "CPU"):
            continue
        detail, why = context_on(name)
        print(f"{name}: {detail or 'UNAVAILABLE'}" + (f" ({why})" if why else ""))
        if detail:
            usable.append(name)
    if not usable:
        print("[WARNING] no GPU platform opened a context: this is a CPU-only node or "
              "the driver is missing. Platform list alone does not establish GPU support.")

    if len(argv) > 1:
        for name in usable + ["CPU"]:
            try:
                print(f"{name}: {benchmark(argv[1], name, argv[2:])}")
            except Exception as exc:
                print(f"{name}: benchmark failed - {type(exc).__name__}: {exc}")
    return 0 if usable else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
