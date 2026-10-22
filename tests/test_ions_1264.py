"""
The 12-6-4 metal model and the metal-parameter set shared by md and QM/MM:

    pytest tests/test_ions_1264.py

Needs OpenMM (for real force fields and energies), not ASH.
"""

import math

import pytest

pytest.importorskip("openmm", reason="the 12-6-4 force is an OpenMM force")

import openmm  # noqa: E402
from openmm import unit  # noqa: E402
from openmm.app import ForceField, NoCutoff, PME, Topology, element  # noqa: E402

from mcube import ions  # noqa: E402

WATER_XML = ["amber14/tip3p.xml"]


def ion_and_water(metal="La", resname="LA", r=2.5, box=None):
    """One ion and one TIP3P water, the oxygen r A from the ion along x."""
    top = Topology()
    chain = top.addChain()
    res = top.addResidue(resname, chain)
    top.addAtom(resname, element.Element.getBySymbol(metal), res)
    water = top.addResidue("HOH", chain)
    o = top.addAtom("O", element.oxygen, water)
    h1 = top.addAtom("H1", element.hydrogen, water)
    h2 = top.addAtom("H2", element.hydrogen, water)
    top.addBond(o, h1)
    top.addBond(o, h2)
    if box:
        top.setPeriodicBoxVectors([openmm.Vec3(box, 0, 0), openmm.Vec3(0, box, 0),
                                   openmm.Vec3(0, 0, box)] * unit.nanometer)
    positions = [openmm.Vec3(0, 0, 0), openmm.Vec3(r / 10, 0, 0),
                 openmm.Vec3(r / 10 + 0.0757, 0.0586, 0), openmm.Vec3(r / 10 - 0.0757, 0.0586, 0)]
    return top, positions * unit.nanometer


def system_for(top, periodic=False):
    ff = ForceField(*WATER_XML)
    if periodic:
        return ff.createSystem(top, nonbondedMethod=PME, nonbondedCutoff=1.0 * unit.nanometer,
                               rigidWater=False)
    return ff.createSystem(top, nonbondedMethod=NoCutoff, rigidWater=False)


def c4_energy(system, positions):
    """Energy (kcal/mol) of the C4 force alone, on CPU: unlike Reference it enforces
    that every nonbonded force has the same exclusions."""
    force = next(f for f in system.getForces() if f.getName() == "mcube_1264_C4")
    force.setForceGroup(7)
    context = openmm.Context(system, openmm.VerletIntegrator(0.001),
                             openmm.Platform.getPlatformByName("CPU"))
    context.setPositions(positions)
    state = context.getState(getEnergy=True, groups={7})
    return state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)


class TestTables:
    def test_1264_lookup_is_the_li_song_merz_row(self):
        sigma, eps = ions.lookup("1264", "LA", 3)
        assert sigma == pytest.approx(2 * 1.758 / 2 ** (1 / 6) / 10)
        assert eps == pytest.approx(0.17997960 * 4.184)

    def test_iod_lookup_is_unchanged(self):
        sigma, _ = ions.lookup("iod", "LA", 3)
        assert sigma == pytest.approx(2 * 1.718 / 2 ** (1 / 6) / 10)

    def test_c4_values(self):
        assert ions.c4_value("LA", 3) == 152.0
        assert ions.c4_value("LU", 3) == 249.0
        assert ions.c4_value("CA", 2) is None   # +2 rows not transcribed

    def test_numeric_set_name_from_yaml(self):
        """`metal_parameters: 1264` arrives from YAML as an int."""
        assert ions.lookup(1264, "LU", 3) is not None

    def test_polarisability_file(self):
        table = ions.read_polarizabilities()
        assert table["OW"] == ions.WATER_POL == 1.444
        assert table["CT"] == 1.061 and table["Na+"] == 0.2495

    @pytest.mark.parametrize("atom_class, expected", [
        ("tip3p-O", 1.444), ("tip3p-H", 0.0), ("tip3p_standard-Cl-", 3.235),
        ("CT", 1.061), ("O2", 0.569), ("cA", 1.061), ("oP", 0.569), ("pA", 1.538),
    ])
    def test_class_translation(self, atom_class, expected):
        assert ions.type_polarizability(atom_class, ions.read_polarizabilities()) == expected

    def test_override_wins(self):
        assert ions.type_polarizability("CT", ions.read_polarizabilities(), {"CT": 2.0}) == 2.0

    def test_every_protein_lipid_and_water_class_is_covered(self):
        table = ions.read_polarizabilities()
        for xml in ("amber14/protein.ff14SB.xml", "amber14/lipid17.xml"):
            ff = ForceField(xml)
            missing = {t.atomClass for t in ff._atomTypes.values()
                       if ions.type_polarizability(t.atomClass, table) is None}
            assert not missing, (xml, missing)
        for c in ("tip3p-O", "tip3p-H", "tip3p_standard-Na+", "tip3p_standard-Cl-"):
            assert ions.type_polarizability(c, table) is not None


class TestApply:
    def test_c4_energy_is_minus_c4_over_r4(self):
        top, pos = ion_and_water(r=2.5)
        system = system_for(top)
        changed = ions.apply_metal_parameters(system, top, ["LA"], "1264", forcefields=WATER_XML)
        assert changed == [0]
        # Water hydrogens have zero polarisability, so only La-O counts, at C4(La) exactly
        assert c4_energy(system, pos) == pytest.approx(-152.0 / 2.5 ** 4, rel=1e-5)

    def test_lj_terms_are_the_1264_set(self):
        top, _ = ion_and_water()
        system = system_for(top)
        ions.apply_metal_parameters(system, top, ["LA"], "1264", forcefields=WATER_XML)
        nb = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
        charge, sigma, eps = nb.getParticleParameters(0)
        assert charge.value_in_unit(unit.elementary_charge) == pytest.approx(3.0)
        assert sigma.value_in_unit(unit.nanometer) == pytest.approx(ions.lookup("1264", "LA", 3)[0])

    def test_qm_ion_has_no_c4(self):
        """In QM/MM the embedding polarises the QM region; C4 would count it twice."""
        top, pos = ion_and_water()
        system = system_for(top)
        ions.apply_metal_parameters(system, top, ["LA"], "1264", forcefields=WATER_XML,
                                    qm_atoms=[0])
        assert c4_energy(system, pos) == 0.0

    def test_applying_twice_adds_one_force(self):
        top, _ = ion_and_water()
        system = system_for(top)
        for _ in range(2):
            ions.apply_metal_parameters(system, top, ["LA"], "1264", forcefields=WATER_XML)
        assert sum(f.getName() == "mcube_1264_C4" for f in system.getForces()) == 1

    def test_periodic_system_takes_the_nonbonded_cutoff(self):
        top, _ = ion_and_water(box=3.0)
        system = system_for(top, periodic=True)
        ions.apply_metal_parameters(system, top, ["LA"], "1264", forcefields=WATER_XML)
        force = next(f for f in system.getForces() if f.getName() == "mcube_1264_C4")
        assert force.getNonbondedMethod() == openmm.CustomNonbondedForce.CutoffPeriodic
        assert force.getCutoffDistance().value_in_unit(unit.nanometer) == pytest.approx(1.0)

    def test_untabulated_ion_raises_under_1264(self):
        top, _ = ion_and_water(metal="Ca", resname="CA")
        with pytest.raises(ValueError, match="no 1264 parameters for CA"):
            ions.apply_metal_parameters(system_for(top), top, ["CA"], "1264",
                                        forcefields=WATER_XML)

    def test_iod_adds_no_c4(self):
        top, _ = ion_and_water(metal="Ca", resname="CA")
        system = system_for(top)
        assert ions.apply_metal_parameters(system, top, ["CA"], "iod") == [0]
        assert not any(f.getName() == "mcube_1264_C4" for f in system.getForces())

    def test_missing_polarisability_names_the_type(self, monkeypatch):
        top, _ = ion_and_water()
        monkeypatch.setattr(ions, "read_polarizabilities", lambda: {"OW": 1.444})
        with pytest.raises(ValueError, match="tip3p-H"):
            ions.apply_metal_parameters(system_for(top), top, ["LA"], "1264",
                                        forcefields=WATER_XML)


class TestStageSet:
    def test_stage_input_wins(self):
        assert ions.metal_parameter_set({"metal_parameters": "1264"},
                                        {"results": {"metal_parameters": "iod"}}) == "1264"

    def test_falls_back_to_what_md_used(self):
        """F1: QM/MM uses md's set unless told otherwise."""
        assert ions.metal_parameter_set({}, {"results": {"metal_parameters": "1264"}}) == "1264"

    def test_then_md_inputs_then_iod(self):
        assert ions.metal_parameter_set({}, {"inputs": {"metal_parameters": "none"}}) == "none"
        assert ions.metal_parameter_set({}, {}) == "iod"
