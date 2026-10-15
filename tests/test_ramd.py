"""
Tests for the RAMD engine and τ analysis on a toy OpenMM system (no ASH needed):

    pytest tests/test_ramd.py

The toy "receptor" is a ring of LJ particles held by position restraints; the
"binder" is a bonded 3-particle cluster sitting in the ring's attractive well,
in a periodic box with no solvent.
"""

import numpy as np
import openmm
import pytest
from openmm import unit
from openmm.app import Topology, Element

from mcube.ramd_engine import RAMD, KCAL_A_TO_KJ_NM
from mcube.ramd_analysis import tau_ramd, ks_exponential, replica_tau

BOX = 6.0  # nm


def toy(n_rec=12):
    system = openmm.System()
    system.setDefaultPeriodicBoxVectors(openmm.Vec3(BOX, 0, 0), openmm.Vec3(0, BOX, 0), openmm.Vec3(0, 0, BOX))
    top = Topology()
    chain_r, chain_b = top.addChain("A"), top.addChain("B")
    pos = []
    centre = np.array([BOX / 2] * 3)
    restr = openmm.CustomExternalForce("500*periodicdistance(x,y,z,x0,y0,z0)^2")
    for p in ("x0", "y0", "z0"):
        restr.addPerParticleParameter(p)
    for k in range(n_rec):
        a = 2 * np.pi * k / n_rec
        x = centre + 0.55 * np.array([np.cos(a), np.sin(a), 0.0])
        i = system.addParticle(40.0)
        restr.addParticle(i, list(x))
        top.addAtom("AR", Element.getBySymbol("Ar"), top.addResidue("REC", chain_r, str(k + 1)))
        pos.append(x)
    system.addForce(restr)
    res_b = top.addResidue("LIG", chain_b, "1")
    bonds = openmm.HarmonicBondForce()
    for k in range(3):
        system.addParticle(20.0)
        top.addAtom(f"C{k}", Element.getBySymbol("C"), res_b)
        pos.append(centre + np.array([0.1 * (k - 1), 0.0, 0.0]))
    bonds.addBond(n_rec, n_rec + 1, 0.1, 5000.0)
    bonds.addBond(n_rec + 1, n_rec + 2, 0.1, 5000.0)
    system.addForce(bonds)
    nb = openmm.NonbondedForce()
    nb.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
    nb.setCutoffDistance(1.2)
    for i in range(n_rec + 3):
        nb.addParticle(0.0, 0.35, 2.0 if i >= n_rec else 0.0)  # binder–ring attraction only
    for i in range(n_rec):
        for j in range(i + 1, n_rec):
            nb.addException(i, j, 0.0, 0.1, 0.0)
    nb.addException(n_rec, n_rec + 1, 0.0, 0.1, 0.0)
    nb.addException(n_rec + 1, n_rec + 2, 0.0, 0.1, 0.0)
    nb.addException(n_rec, n_rec + 2, 0.0, 0.1, 0.0)
    system.addForce(nb)
    system.addForce(openmm.CMMotionRemover())
    receptor = list(range(n_rec))
    binder = [n_rec, n_rec + 1, n_rec + 2]
    return system, top, np.array(pos), binder, receptor


def engine(force=14.0, **kw):
    system, top, pos, binder, receptor = toy()
    opts = dict(force=force, eval_steps=50, rmin=0.025, max_distance=20.0, timestep=0.002,
                temperature=300.0, friction=1.0, pressure=None, platform="Reference")
    opts.update(kw)
    e = RAMD(system, top, binder, receptor, **opts)
    e.set_coordinates(pos * unit.nanometer)
    return e, pos, binder, receptor


def test_net_force_on_binder_is_F_along_direction():
    e, pos, binder, _ = engine(force=14.0)
    d = np.array([1.0, 2.0, -2.0]) / 3.0

    def binder_force():
        f = e.context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer)
        return f[binder].sum(0)
    e.set_direction(np.zeros(3))
    f0 = binder_force()
    e.set_direction(d)
    f1 = binder_force()
    assert np.allclose(f1 - f0, 14.0 * KCAL_A_TO_KJ_NM * d, atol=1e-3)


def test_monitor_matches_positions():
    e, pos, binder, receptor = engine()
    m = np.array([e.system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(len(pos))])
    cb = (pos[binder] * m[binder, None]).sum(0) / m[binder].sum()
    cr = (pos[receptor] * m[receptor, None]).sum(0) / m[receptor].sum()
    sep, vec = e.monitor()
    assert sep == pytest.approx(np.linalg.norm(cb - cr), abs=1e-5)
    assert np.allclose(vec, cb - cr, atol=1e-5)


def test_strong_force_dissociates_and_weak_does_not():
    e, *_ = engine(force=40.0)
    state = e.equilibrate(2.0, seed=1)
    out = e.run_trajectory(state, seed=5, max_time=200.0)
    assert out["dissociated"] and 0 < out["time_ps"] < 200.0
    assert out["final_distance_A"] >= 20.0

    e0, *_ = engine(force=0.0)
    state0 = e0.equilibrate(2.0, seed=1)
    out0 = e0.run_trajectory(state0, seed=5, max_time=10.0)
    assert not out0["dissociated"] and out0["time_ps"] == pytest.approx(10.0)
    # No force means no binder motion beyond thermal: directions get redrawn sometimes or not, but never escape
    assert out0["final_distance_A"] < 20.0


def test_trajectories_restart_from_the_same_state():
    e, *_ = engine(force=40.0)
    state = e.equilibrate(2.0, seed=1)
    a = e.run_trajectory(state, seed=11, max_time=200.0)
    b = e.run_trajectory(state, seed=11, max_time=200.0)
    assert a["start_distance_A"] == b["start_distance_A"]


def test_cm_motion_remover_can_be_dropped():
    e, *_ = engine(remove_cm_motion=False)
    assert not any(isinstance(e.system.getForce(k), openmm.CMMotionRemover) for k in range(e.system.getNumForces()))
    e2, *_ = engine()
    assert any(isinstance(e2.system.getForce(k), openmm.CMMotionRemover) for k in range(e2.system.getNumForces()))


def test_tau_of_exponential_times():
    rng = np.random.default_rng(3)
    tau_exp = 1000.0
    reps = {r: rng.exponential(tau_exp, size=30) for r in range(4)}
    res = tau_ramd(reps, rounds=2000, seed=1)
    # τRAMD's τ is the 50% time: ln2 * the exponential time constant
    assert res["tau_ps"] == pytest.approx(np.log(2) * tau_exp, rel=0.35)
    assert all(v["ks"] < 0.45 for v in res["replicas"].values())
    assert res["n_replicas"] == 4


def test_censored_majority_is_flagged():
    times = [100, 200, 5000, 5000, 5000, 5000]
    res = replica_tau(times, censored=[False, False, True, True, True, True], rounds=200)
    assert res["lower_bound"]


def test_ks_perfect_fit_is_small():
    tau = 100.0
    t = -tau * np.log(1 - (np.arange(1, 201) - 0.5) / 200)  # exact quantiles
    assert ks_exponential(t, tau) < 0.01
