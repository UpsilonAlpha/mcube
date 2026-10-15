"""
Random acceleration MD (RAMD) on an OpenMM System.

A constant-magnitude force F acts on the binder's centre of mass along a random
direction. Every `eval_steps` the binder COM displacement relative to the
receptor COM is checked: if it moved less than `rmin`, a new random direction
is drawn (Lüdemann et al. 2000; Kokh et al. 2018/2020). A trajectory ends when
the binder–receptor COM distance reaches `max_distance` (dissociated) or after
`max_time`.

Monitored quantities are collective variables of a zero-energy CustomCVForce,
so they are read with getCollectiveVariableValues() instead of transferring all
positions every check. Pure OpenMM: no ASH dependency.
"""

from __future__ import annotations

import time
from typing import Dict, Any, Optional, Sequence

import numpy as np
import openmm
from openmm import unit
from openmm.app import Simulation, DCDReporter, PDBFile

KCAL_A_TO_KJ_NM = 41.84  # 1 kcal/mol/Å = 41.84 kJ/mol/nm


def random_direction(rng: np.random.Generator) -> np.ndarray:
    v = rng.normal(size=3)
    return v / np.linalg.norm(v)


class RAMD:
    """
    Parameters (lengths in Å, times in ps, force in kcal/mol/Å, as in mcube configs):
        system, topology   OpenMM System (copied, not modified) and Topology
        binder, receptor   0-based atom indices
        force              RAMD force magnitude (14 kcal/mol/Å = 585.76 kJ/mol/nm)
        eval_steps         MD steps between direction checks
        rmin               minimum binder COM displacement (relative to receptor COM) per check
        max_distance       binder–receptor COM distance that ends a trajectory
        periodic           minimum-image COM distance (explicit solvent)
        pressure           bar; adds a MonteCarloBarostat (periodic only); None for NVT
        remove_cm_motion   keep the system's CMMotionRemover (False removes it)
    """

    def __init__(self, system, topology, binder: Sequence[int], receptor: Sequence[int], *,
                 force: float = 14.0, eval_steps: int = 50, rmin: float = 0.025,
                 max_distance: float = 30.0, timestep: float = 0.002, temperature: float = 300.0,
                 friction: float = 1.0, pressure: Optional[float] = 1.0, periodic: bool = True,
                 remove_cm_motion: bool = True, platform: Optional[str] = None,
                 properties: Optional[Dict[str, str]] = None):
        if set(binder) & set(receptor):
            raise ValueError("binder and receptor share atoms")
        self.topology = topology
        self.force_kj_nm = force * KCAL_A_TO_KJ_NM
        self.eval_steps = int(eval_steps)
        self.rmin_nm = rmin / 10.0
        self.max_distance_nm = max_distance / 10.0
        self.timestep = timestep
        self.temperature = temperature
        self.periodic = periodic

        system = openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(system))
        if not remove_cm_motion:
            for k in reversed(range(system.getNumForces())):
                if isinstance(system.getForce(k), openmm.CMMotionRemover):
                    system.removeForce(k)

        # RAMD force on the binder COM: E = -F·r_com, so the force is +F along the direction.
        # Linear in the coordinates, so periodic wrapping does not matter.
        ramd = openmm.CustomCentroidBondForce(1, "-(fx*x1 + fy*y1 + fz*z1)")
        for name in ("fx", "fy", "fz"):
            ramd.addGlobalParameter(name, 0.0)
        ramd.addGroup([int(i) for i in binder])  # weights omitted: mass-weighted COM
        ramd.addBond([0], [])
        system.addForce(ramd)

        # Monitored quantities as CVs of a zero-energy force (no effect on the dynamics)
        def centroid_cv(expression, pbc):
            f = openmm.CustomCentroidBondForce(2, expression)
            f.addGroup([int(i) for i in binder])
            f.addGroup([int(i) for i in receptor])
            f.addBond([0, 1], [])
            f.setUsesPeriodicBoundaryConditions(pbc)
            return f
        self.cv = openmm.CustomCVForce("0")
        self.cv.addCollectiveVariable("sep", centroid_cv("distance(g1,g2)", periodic))
        for axis in "xyz":
            self.cv.addCollectiveVariable(f"d{axis}", centroid_cv(f"{axis}1-{axis}2", False))
        system.addForce(self.cv)

        if pressure is not None and periodic:
            system.addForce(openmm.MonteCarloBarostat(pressure * unit.bar, temperature * unit.kelvin))
        self.system = system

        integrator = openmm.LangevinMiddleIntegrator(temperature * unit.kelvin, friction / unit.picosecond,
                                                     timestep * unit.picosecond)
        plat = openmm.Platform.getPlatformByName(platform) if platform else None
        if plat is None:
            self.simulation = Simulation(topology, system, integrator)
        else:
            self.simulation = Simulation(topology, system, integrator, plat, properties or {})
        self.context = self.simulation.context

    # ---------- helpers ----------
    def set_direction(self, d: np.ndarray):
        for name, val in zip(("fx", "fy", "fz"), self.force_kj_nm * np.asarray(d)):
            self.context.setParameter(name, float(val))

    def monitor(self):
        """(binder–receptor COM distance in nm, binder-minus-receptor COM vector in nm)."""
        v = self.cv.getCollectiveVariableValues(self.context)
        return v[0], np.array(v[1:4])

    def box_min_nm(self) -> float:
        a, b, c = self.context.getState().getPeriodicBoxVectors()
        return min(a[0], b[1], c[2]).value_in_unit(unit.nanometer)

    def set_coordinates(self, positions, box_vectors=None):
        if box_vectors is not None:
            self.context.setPeriodicBoxVectors(*box_vectors)
        self.context.setPositions(positions)
        self.context.applyConstraints(1e-6)

    # ---------- stages ----------
    def equilibrate(self, time_ps: float, seed: int, report: Optional[str] = None, interval_ps: float = 10.0):
        """Unbiased MD (RAMD force off) from the current coordinates with new velocities."""
        self.set_direction(np.zeros(3))
        self.context.setVelocitiesToTemperature(self.temperature * unit.kelvin, int(seed))
        nsteps = int(round(time_ps / self.timestep))
        if report:
            self.simulation.reporters.append(DCDReporter(report, max(1, int(round(interval_ps / self.timestep)))))
        self.simulation.step(nsteps)
        self.simulation.reporters.clear()
        return self.context.getState(getPositions=True, getVelocities=True, getParameters=True,
                                     enforcePeriodicBox=False)

    def run_trajectory(self, start_state, seed: int, max_time: float, dcd: Optional[str] = None,
                       dcd_interval_ps: float = 0.0, exit_pdb: Optional[str] = None) -> Dict[str, Any]:
        """One RAMD trajectory from `start_state` with new velocities. Times in ps."""
        rng = np.random.default_rng(seed)
        self.context.setState(start_state)
        self.context.setTime(0.0)
        self.context.setVelocitiesToTemperature(self.temperature * unit.kelvin, int(seed) % 2**31)
        if dcd and dcd_interval_ps > 0:
            self.simulation.reporters.append(DCDReporter(dcd, max(1, int(round(dcd_interval_ps / self.timestep)))))

        self.set_direction(random_direction(rng))
        sep0, prev = self.monitor()
        if sep0 >= self.max_distance_nm:
            raise ValueError(f"Start COM distance {10 * sep0:.1f} Å is already beyond max_distance "
                             f"{10 * self.max_distance_nm:.1f} Å")
        half_box = self.box_min_nm() / 2 if self.periodic else np.inf
        max_steps = int(round(max_time / self.timestep))
        steps, changes, dissociated = 0, 0, False
        wall = time.time()
        sep = sep0
        while steps < max_steps:
            self.simulation.step(self.eval_steps)
            steps += self.eval_steps
            sep, cur = self.monitor()
            if sep >= self.max_distance_nm:
                dissociated = True
                break
            moved = np.linalg.norm(cur - prev)
            # A jump of half a box is a periodic wrap of the raw COMs, not slow motion
            if moved < self.rmin_nm and moved < half_box:
                self.set_direction(random_direction(rng))
                changes += 1
            prev = cur
        self.simulation.reporters.clear()
        if exit_pdb:
            pos = self.context.getState(getPositions=True).getPositions()
            with open(exit_pdb, "w") as f:
                PDBFile.writeFile(self.topology, pos, f, keepIds=True)
        return {
            "dissociated": dissociated,
            "time_ps": round(steps * self.timestep, 4),
            "steps": steps,
            "direction_changes": changes,
            "start_distance_A": round(10 * sep0, 3),
            "final_distance_A": round(10 * sep, 3),
            "wall_s": round(time.time() - wall, 1),
        }


def com_and_radius(positions_nm: np.ndarray, masses: np.ndarray, idx: Sequence[int]):
    """Mass-weighted COM and the largest atom distance from it (both nm)."""
    x = positions_nm[list(idx)]
    m = masses[list(idx)]
    com = (x * m[:, None]).sum(0) / m.sum()
    return com, float(np.max(np.linalg.norm(x - com, axis=1)))
