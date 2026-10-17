"""
Tests for the explicit membrane mode on OpenMM's own pre-equilibrated POPC patch
(no ASH needed):

    pytest tests/test_membrane.py

The patch (openmm/app/data/POPC.pdb) is the same one `Modeller.addMembrane` tiles to
build a bilayer, so it is a real lipid/water system with a known area per lipid.
"""

import sys
from pathlib import Path

import openmm
import pytest
from openmm import unit
from openmm.app import ForceField, PDBFile
from openmm.app import forcefield as openmm_forcefield

from mcube import membrane

FORCEFIELDS = ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]  # mcube's defaults


@pytest.fixture(scope="module")
def patch():
    """OpenMM's pre-equilibrated POPC patch: 128 lipids and 5120 waters in a box."""
    data = Path(openmm_forcefield.__file__).parent / "data" / "POPC.pdb"
    return PDBFile(str(data))


@pytest.fixture(scope="module")
def system(patch):
    return membrane.membrane_system(patch.topology, membrane.lipid_forcefields(FORCEFIELDS),
                                    temperature=300)


def forces(system):
    return [f.__class__.__name__ for f in system.getForces()]


# ==============================
# Force field
# ==============================

def test_lipid_forcefield_ships_with_openmm():
    assert membrane.check_lipid_forcefield() == membrane.LIPID_XML


def test_missing_lipid_forcefield_is_an_error():
    with pytest.raises(FileNotFoundError, match="does not provide"):
        membrane.check_lipid_forcefield("amber14/no-such-lipid.xml")


def test_lipid_forcefield_is_appended_once():
    once = membrane.lipid_forcefields(FORCEFIELDS)
    assert once == FORCEFIELDS + [membrane.LIPID_XML]
    assert membrane.lipid_forcefields(once) == once  # idempotent
    assert FORCEFIELDS == ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]  # not mutated


def test_lipid_type_is_checked_against_openmms_patches():
    assert membrane.check_lipid_type("popc") == "POPC"
    assert membrane.check_lipid_type(None) == "POPC"
    with pytest.raises(ValueError, match="POPG"):
        membrane.check_lipid_type("POPG")


def test_lipids_can_be_typed(patch):
    """Without the lipid XML the defaults cannot type POP, which is the whole point."""
    with pytest.raises(ValueError, match="[Nn]o template"):
        ForceField(*FORCEFIELDS).createSystem(patch.topology)


# ==============================
# Barostat
# ==============================

def test_membrane_barostat_lets_area_and_thickness_move_separately():
    barostat = membrane.membrane_barostat(300)
    assert barostat.getXYMode() == openmm.MonteCarloMembraneBarostat.XYIsotropic
    assert barostat.getZMode() == openmm.MonteCarloMembraneBarostat.ZFree
    assert barostat.getDefaultSurfaceTension().value_in_unit(unit.bar * unit.nanometer) == 0.0
    assert barostat.getDefaultTemperature().value_in_unit(unit.kelvin) == 300
    assert barostat.getDefaultPressure().value_in_unit(unit.bar) == 1.0
    assert barostat.getFrequency() == membrane.BAROSTAT_FREQUENCY


def test_membrane_barostat_takes_the_settings_given():
    barostat = membrane.membrane_barostat(310, pressure=1.5, surface_tension=20.0, frequency=50)
    assert barostat.getDefaultTemperature().value_in_unit(unit.kelvin) == 310
    assert barostat.getDefaultPressure().value_in_unit(unit.bar) == 1.5
    assert barostat.getDefaultSurfaceTension().value_in_unit(unit.bar * unit.nanometer) == 20.0
    assert barostat.getFrequency() == 50


def test_barostat_needs_a_temperature():
    with pytest.raises(ValueError, match="temp"):
        membrane.membrane_barostat(None)


def test_one_barostat_is_required():
    system = openmm.System()
    with pytest.raises(RuntimeError, match="found none"):
        membrane.check_single_barostat(system)

    system.addForce(membrane.membrane_barostat(300))
    membrane.check_single_barostat(system)

    # The failure this guards against: ASH adding its isotropic barostat on top
    system.addForce(openmm.MonteCarloBarostat(1 * unit.bar, 300 * unit.kelvin))
    with pytest.raises(RuntimeError, match="MonteCarloBarostat"):
        membrane.check_single_barostat(system)


def test_an_isotropic_barostat_alone_is_rejected():
    system = openmm.System()
    system.addForce(openmm.MonteCarloBarostat(1 * unit.bar, 300 * unit.kelvin))
    with pytest.raises(RuntimeError, match="MonteCarloMembraneBarostat"):
        membrane.check_single_barostat(system)


def test_ash_would_not_notice_the_membrane_barostat():
    """
    Why the md stage must not pass `barostat=` for a membrane system: ASH's guard is
    `if "MonteCarloBarostat" not in [f.__class__.__name__ for f in system.getForces()]`,
    which does not match MonteCarloMembraneBarostat, so ASH would add a second barostat.
    """
    system = openmm.System()
    system.addForce(membrane.membrane_barostat(300))
    assert "MonteCarloBarostat" not in forces(system)


# ==============================
# System
# ==============================

def test_membrane_system_is_periodic_pme_with_constrained_hydrogens(system, patch):
    assert membrane.barostat_names(system) == [membrane.MEMBRANE_BAROSTAT]
    assert system.usesPeriodicBoundaryConditions()
    nonbonded = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    assert nonbonded.getNonbondedMethod() == openmm.NonbondedForce.PME
    assert nonbonded.getCutoffDistance().value_in_unit(unit.angstrom) == pytest.approx(
        membrane.NONBONDED_CUTOFF)
    # HBonds + rigid water: every bond to a hydrogen is constrained
    assert system.getNumConstraints() > 0
    hydrogens = [a.index for a in patch.topology.atoms() if a.element.symbol == "H"]
    assert system.getParticleMass(hydrogens[0]).value_in_unit(unit.amu) == pytest.approx(
        membrane.HYDROGEN_MASS)


def test_membrane_system_needs_a_periodic_box(patch):
    topology = PDBFile(str(Path(openmm_forcefield.__file__).parent / "data" / "POPC.pdb")).topology
    topology.setPeriodicBoxVectors(None)
    with pytest.raises(ValueError, match="CRYST1"):
        membrane.membrane_system(topology, membrane.lipid_forcefields(FORCEFIELDS), temperature=300)


def test_membrane_system_runs(system, patch):
    """The system is usable: a context can be built and gives a finite energy."""
    integrator = openmm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                                 0.002 * unit.picoseconds)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("CPU"))
    context.setPositions(patch.positions)
    energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole)
    assert energy == energy and abs(energy) < 1e9  # finite, not NaN


# ==============================
# Reporting
# ==============================

def test_lipids_are_counted(patch):
    assert membrane.count_lipids(patch.topology) == 128


def test_area_per_lipid_of_the_equilibrated_patch(patch):
    """OpenMM's patch is pre-equilibrated, so it should already sit at a sane POPC area."""
    lx, ly, _ = patch.topology.getUnitCellDimensions().value_in_unit(unit.angstrom)
    apl = membrane.area_per_lipid(lx, ly, membrane.count_lipids(patch.topology))
    assert 60.0 < apl < 68.0


def test_area_per_lipid_counts_both_leaflets():
    # 64 lipids per leaflet over a 64 x 64 Å box
    assert membrane.area_per_lipid(64.0, 64.0, 128) == pytest.approx(64.0)
    with pytest.raises(ValueError):
        membrane.area_per_lipid(64.0, 64.0, 1)


def test_membrane_report_writes_a_csv_and_yaml_safe_results(patch, tmp_path, monkeypatch):
    """
    mdtraj is stubbed out so the report can be checked without reading a trajectory:
    what matters here is the CSV, and that the numbers are plain floats rather than
    numpy scalars, since they are written into the config file with `yaml.safe_dump`.
    """
    import types

    import numpy as np
    import yaml

    monkeypatch.chdir(tmp_path)
    trajectory = types.SimpleNamespace(
        unitcell_lengths=np.array([[6.37, 6.69, 7.31], [6.40, 6.72, 7.28]], dtype=np.float32),
        topology=types.SimpleNamespace(to_openmm=lambda: patch.topology))
    monkeypatch.setitem(sys.modules, "mdtraj", types.SimpleNamespace(
        load=lambda trajfile, top: trajectory))

    report = membrane.membrane_report("lastframe.pdb", "traj.dcd")

    assert report["lipids"] == 128
    assert 60.0 < report["area_per_lipid"] < 68.0
    # The report goes into the config file, so it has to survive yaml.safe_dump
    assert yaml.safe_load(yaml.safe_dump(report))["lipids"] == 128
    rows = Path("membrane_box.csv").read_text().splitlines()
    assert rows[0].startswith("frame,box_x,box_y,box_z,area_per_lipid")
    assert len(rows) == 3  # header + two frames


def test_box_table_has_a_row_per_frame():
    rows = membrane.box_table([(60.0, 60.0, 90.0), (62.0, 62.0, 88.0)], 128)
    assert [r["frame"] for r in rows] == [0, 1]
    assert rows[1]["area_per_lipid"] == pytest.approx(62.0 * 62.0 / 64)
    assert rows[1]["box_z"] == 88.0
