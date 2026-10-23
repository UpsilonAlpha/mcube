#!/usr/bin/env python3
"""
Is the Lu3+ 12-6-4 coordination number (7.6 in OpenMM, 9.0 in Li, Song & Merz 2015, SI.4)
an OpenMM/mcube problem, or does Amber's own engine give it too?

check_1264.py showed mcube's energies equal those of the Amber-built prmtop, but it evaluated
both in OpenMM, so a fault in how OpenMM evaluates the C4 term would pass it. Here the same
Amber-built box (tleap + ParmEd add12_6_4, as in check_1264.py) is simulated with no mcube
code at all, under the paper's protocol (1 fs, Langevin 5 ps^-1, 1 atm isotropic, 10 Å
cutoff, PME, SHAKE on H, rigid water, 300 K), in either engine:

  openmm   OpenMM on the prmtop (AmberPrmtopFile reads LENNARD_JONES_CCOEF)
  sander   AmberTools' sander with lj1264 = 1, run in 100 ps chunks so a job cut short by
           the time limit still leaves usable frames

and the shell is analysed with bulk_ion.analyse, the code behind every CN in the log.

  python amber_crosscheck.py --pdb LU_1264_lsm/system.pdb --ion LU --engine sander --out xc_LU_sander
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bulk_ion import analyse                      # noqa: E402
from check_1264 import TLEAP, amber_pdb, box_of   # noqa: E402

SANDER_IN = """{title}
 &cntrl
  imin=0, irest={irest}, ntx={ntx}, nstlim={nstlim}, dt=0.001,
  ntc=2, ntf=2, cut=10.0, ntb=2, ntp=1, barostat=1, taup=1.0, pres0=1.0,
  ntt=3, gamma_ln=5.0, temp0=300.0, tempi=300.0, ig=-1,
  ntpr=1000, ntwx=500, ntwr=10000, ioutfm=1, iwrap=1,
  lj1264=1,
 /
"""


def build(pdb_src, ion):
    a, b, c = box_of(str(pdb_src))
    amber_pdb(str(pdb_src), "amber.pdb")
    Path("leap.in").write_text(TLEAP.format(pdb="amber.pdb", a=a, b=b, c=c))
    subprocess.run(["tleap", "-f", "leap.in"], check=True, stdout=open("tleap.log", "w"),
                   stderr=subprocess.STDOUT)
    pol = Path(os.environ["CONDA_PREFIX"]) / "dat" / "leap" / "parm" / "lj_1264_pol.dat"
    Path("parmed.in").write_text(f"add12_6_4 :{ion} watermodel TIP3P polfile {pol}\n"
                                 "outparm sys1264.prmtop\nquit\n")
    subprocess.run(["parmed", "-i", "parmed.in", "sys.prmtop"], check=True,
                   stdout=open("parmed.log", "w"), stderr=subprocess.STDOUT)


def run_openmm(equil_ps, prod_ps, cores):
    import openmm
    from openmm import unit
    from openmm.app import AmberInpcrdFile, AmberPrmtopFile, DCDReporter, HBonds, PME, Simulation
    prmtop = AmberPrmtopFile("sys1264.prmtop")
    inpcrd = AmberInpcrdFile("sys.inpcrd")
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0 * unit.nanometer,
                                 constraints=HBonds, rigidWater=True)
    system.addForce(openmm.MonteCarloBarostat(1.0 * unit.atmosphere, 300 * unit.kelvin))
    integ = openmm.LangevinMiddleIntegrator(300 * unit.kelvin, 5.0 / unit.picosecond,
                                            0.001 * unit.picosecond)
    sim = Simulation(prmtop.topology, system, integ, openmm.Platform.getPlatformByName("CPU"),
                     {"Threads": str(cores)})
    sim.context.setPositions(inpcrd.positions)
    sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    sim.minimizeEnergy()
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin)
    sim.step(int(equil_ps * 1000))
    sim.reporters.append(DCDReporter("prod.dcd", 500))
    t0 = time.time()
    sim.step(int(prod_ps * 1000))
    return ["prod.dcd"], round(prod_ps / 1000 / ((time.time() - t0) / 86400), 1)


def run_sander(equil_ps, prod_ps, chunk_ps=100):
    def go(name, inp, rst_in, nstlim, irest):
        Path(f"{name}.in").write_text(SANDER_IN.format(title=name, irest=irest, ntx=5 if irest else 1,
                                                       nstlim=nstlim))
        subprocess.run(["sander", "-O", "-i", f"{name}.in", "-o", f"{name}.out", "-p",
                        "sys1264.prmtop", "-c", rst_in, "-r", f"{name}.rst", "-x", f"{name}.nc",
                        "-inf", f"{name}.info"], check=True)
    # minimise first so the tleap box starts clean
    Path("min.in").write_text("min\n &cntrl\n  imin=1, maxcyc=2000, ncyc=1000, cut=10.0, ntb=1,"
                              " lj1264=1,\n /\n")
    subprocess.run(["sander", "-O", "-i", "min.in", "-o", "min.out", "-p", "sys1264.prmtop",
                    "-c", "sys.inpcrd", "-r", "min.rst"], check=True)
    go("equil", None, "min.rst", int(equil_ps * 1000), 0)
    prev, trajs, t0 = "equil.rst", [], time.time()
    done = 0.0
    for k in range(int(np.ceil(prod_ps / chunk_ps))):
        length = min(chunk_ps, prod_ps - done)
        go(f"prod{k:02d}", None, prev, int(round(length * 1000)), 1)
        done += length
        prev = f"prod{k:02d}.rst"
        trajs.append(f"prod{k:02d}.nc")
        speed = round(done / 1000 / ((time.time() - t0) / 86400), 1)
        shell_report(trajs, speed)      # usable even if the job is cut short later
    return trajs, speed


def shell_report(trajs, speed, engine="sander"):
    import mdtraj
    traj = mdtraj.load(trajs, top="sys1264.prmtop")
    top = traj.topology
    ion = next(a.index for a in top.atoms if a.residue.name == ARGS.ion)
    oxygens = [a.index for a in top.atoms if a.residue.name in ("WAT", "HOH") and a.element.symbol == "O"]
    frames = traj.xyz[:, [ion] + oxygens, :] * 10.0
    boxes = traj.unitcell_lengths * 10.0
    report = analyse(frames, boxes, 0, list(range(1, len(oxygens) + 1)), 0.5, 2.0)
    report = {k: v for k, v in report.items() if k != "rdf"}
    report.update({"engine": engine, "ion": ARGS.ion, "frames_ns": round(len(traj) * 0.5 / 1000, 2),
                   "ns_per_day": speed, "waters": len(oxygens)})
    Path("report.json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report), flush=True)


def main():
    global ARGS
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pdb", required=True)
    p.add_argument("--ion", required=True)
    p.add_argument("--engine", choices=["openmm", "sander"], required=True)
    p.add_argument("--equil", type=float, default=500.0, help="NPT ps before sampling (paper: 500)")
    p.add_argument("--prod", type=float, default=2000.0, help="sampled ps (paper: 2000)")
    p.add_argument("--cores", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    p.add_argument("--out", required=True)
    ARGS = p.parse_args()
    pdb = Path(ARGS.pdb).resolve()
    out = Path(ARGS.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    os.chdir(out)
    build(pdb, ARGS.ion)
    if ARGS.engine == "openmm":
        trajs, speed = run_openmm(ARGS.equil, ARGS.prod, ARGS.cores)
    else:
        trajs, speed = run_sander(ARGS.equil, ARGS.prod)
    shell_report(trajs, speed, ARGS.engine)


if __name__ == "__main__":
    main()
