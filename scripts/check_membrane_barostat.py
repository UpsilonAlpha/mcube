#!/usr/bin/env python3
"""
Check that mcube's membrane mode gives a bilayer a sane area per lipid and a stable
box, and that an isotropic barostat (the only kind ASH adds) cannot.

Runs short NPT simulations on OpenMM's own pre-equilibrated lipid patch
(`openmm/app/data/POPC.pdb`, 128 lipids and 5120 waters), built through
`mcube.membrane.membrane_system`, i.e. the same code path the md stage uses. No ASH
and no input structure needed, so it also serves as a check that the installed OpenMM
has usable lipid parameters.

Three cases:
    equilibrated   the patch as it comes, membrane barostat -> area per lipid and box
                   must stay put
    stretched      the box scaled by --stretch in x and y, membrane barostat -> the
                   area per lipid must relax back down at constant thickness
    isotropic      the same stretched box with an isotropic MonteCarloBarostat ->
                   control: it restores the density by thinning the bilayer instead,
                   and leaves the area per lipid too high

Per case it writes <case>.csv (time_ps, box_x, box_y, box_z, area_per_lipid,
p_p_separation, upper_leaflet, lower_leaflet) and prints the starting values and the
mean over the second half of the run. The P-P separation (mean z of the phosphates of
one leaflet minus the other) is the bilayer thickness, and the leaflet counts should
stay equal: a bilayer that has started to distort loses both.

POPC is roughly 60-68 Å² per lipid and 36-39 Å P-P at 300 K.

Usage:
    python check_membrane_barostat.py [--ps 40] [--lipid POPC] [--stretch 1.08]
        [--temp 300] [--timestep 0.004] [--platform CPU] [--outdir .]
        [--cases equilibrated stretched isotropic]
"""

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
import openmm
import openmm.app as app
from openmm import unit

from mcube import membrane

FORCEFIELDS = ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]  # mcube's defaults
CASES = ("equilibrated", "stretched", "isotropic")


def patch(lipid):
    """OpenMM's pre-equilibrated patch for this lipid type."""
    name = membrane.check_lipid_type(lipid) + ".pdb"
    return app.PDBFile(str(Path(app.forcefield.__file__).parent / "data" / name))


def phosphorus(topology):
    return [a.index for a in topology.atoms()
            if a.element is not None and a.element.symbol == "P"
            and a.residue.name.upper() in membrane.LIPID_RESNAMES]


def geometry(positions, box, p_atoms):
    """Area per lipid, bilayer thickness and lipids per leaflet, from one frame."""
    lx, ly, lz = (box[i][i].value_in_unit(unit.angstrom) for i in range(3))
    z = np.array([positions[i][2].value_in_unit(unit.angstrom) for i in p_atoms])
    upper, lower = z[z > z.mean()], z[z <= z.mean()]
    return {
        "box_x": lx, "box_y": ly, "box_z": lz,
        "area_per_lipid": membrane.area_per_lipid(lx, ly, len(p_atoms)),
        "p_p_separation": upper.mean() - lower.mean(),
        "upper_leaflet": len(upper), "lower_leaflet": len(lower),
    }


def stretch_laterally(positions, box, factor):
    """Scale the box and the coordinates in x and y, i.e. put the bilayer under tension."""
    box = [openmm.Vec3(box[0][0] * factor, box[0][1], box[0][2]),
           openmm.Vec3(box[1][0], box[1][1] * factor, box[1][2]),
           box[2]]
    positions = unit.Quantity(
        [openmm.Vec3(p[0] * factor, p[1] * factor, p[2])
         for p in positions.value_in_unit(unit.nanometer)], unit.nanometer)
    return positions, box


def build(case, args):
    pdb = patch(args.lipid)
    topology, positions = pdb.topology, pdb.positions
    box = topology.getPeriodicBoxVectors()
    if case != "equilibrated":
        positions, box = stretch_laterally(positions, box, args.stretch)
        topology.setPeriodicBoxVectors(box)

    system = membrane.membrane_system(topology, membrane.lipid_forcefields(FORCEFIELDS),
                                      temperature=args.temp)
    if case == "isotropic":
        # What ASH would give a membrane system: one isotropic MonteCarloBarostat
        names = [f.__class__.__name__ for f in system.getForces()]
        system.removeForce(names.index(membrane.MEMBRANE_BAROSTAT))
        system.addForce(openmm.MonteCarloBarostat(1 * unit.bar, args.temp * unit.kelvin,
                                                 membrane.BAROSTAT_FREQUENCY))
    return topology, positions, box, system


def run(case, args):
    topology, positions, box, system = build(case, args)
    print(f"[{case}] barostat: {', '.join(membrane.barostat_names(system))}")

    integrator = openmm.LangevinMiddleIntegrator(args.temp * unit.kelvin, 1 / unit.picosecond,
                                                 args.timestep * unit.picoseconds)
    simulation = app.Simulation(topology, system, integrator,
                                openmm.Platform.getPlatformByName(args.platform))
    simulation.context.setPositions(positions)
    simulation.context.setPeriodicBoxVectors(*box)
    simulation.minimizeEnergy(maxIterations=500)
    simulation.context.setVelocitiesToTemperature(args.temp * unit.kelvin)

    p_atoms = phosphorus(topology)
    steps = int(args.ps / args.timestep)
    every = max(1, steps // 40)
    rows, start = [], time.time()
    for done in range(0, steps + every, every):
        state = simulation.context.getState(getPositions=True)
        row = {"time_ps": done * args.timestep}
        row.update(geometry(state.getPositions(), state.getPeriodicBoxVectors(), p_atoms))
        rows.append(row)
        if done < steps:
            simulation.step(every)

    out = Path(args.outdir) / f"{case}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    half = rows[len(rows) // 2:]
    mean = {k: float(np.mean([r[k] for r in half]))
            for k in ("box_x", "box_y", "box_z", "area_per_lipid", "p_p_separation")}
    spread = float(np.std([r["area_per_lipid"] for r in half]))
    first = rows[0]
    print(f"[{case}] {args.ps:g} ps in {(time.time() - start) / 60:.1f} min -> {out}")
    print(f"[{case}] start: area per lipid {first['area_per_lipid']:.1f} Å², box "
          f"{first['box_x']:.1f} x {first['box_y']:.1f} x {first['box_z']:.1f} Å, "
          f"P-P {first['p_p_separation']:.1f} Å")
    print(f"[{case}] second half: area per lipid {mean['area_per_lipid']:.1f} ± {spread:.2f} Å², "
          f"box {mean['box_x']:.1f} x {mean['box_y']:.1f} x {mean['box_z']:.1f} Å, "
          f"P-P {mean['p_p_separation']:.1f} Å, leaflets "
          f"{rows[-1]['upper_leaflet']}/{rows[-1]['lower_leaflet']}")
    sys.stdout.flush()
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ps", type=float, default=40, help="length of each run in ps")
    parser.add_argument("--lipid", default="POPC")
    parser.add_argument("--stretch", type=float, default=1.08,
                        help="lateral strain applied in the stretched cases")
    parser.add_argument("--temp", type=float, default=300)
    parser.add_argument("--timestep", type=float, default=0.004,
                        help="ps; 4 fs as mcube's production default (HBonds + 1.5 amu H)")
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--cases", nargs="+", choices=CASES, default=list(CASES))
    args = parser.parse_args()
    for case in args.cases:
        run(case, args)


if __name__ == "__main__":
    main()
