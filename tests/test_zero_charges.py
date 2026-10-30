"""
`zero_charges`: a residue's MM point charges deleted for a per-residue
electrostatic decomposition:

    pytest tests/test_zero_charges.py

Runs without ASH. `charge_deletion` is pure (charges in, record out),
`resolve_residues` only needs an OpenMM topology, and `apply_zero_charges` only
needs an object carrying an OpenMM system, a topology and ASH's own charge list -
so `mcube.qmmm` and `mcube.sp`, which both `from ash import *` at module scope,
are imported against a stub `ash` module, the approach tests/test_restraints.py
and tests/test_umbrella.py take. Needs OpenMM (the `mcube-test` env). The QM/MM
energies themselves are run on Bunya.
"""

import sys
import types

import numpy as np
import pytest

pytest.importorskip("openmm", reason="point charges live in an OpenMM NonbondedForce")

import openmm  # noqa: E402
from openmm import unit  # noqa: E402
from openmm.app import Topology, element  # noqa: E402

# Every ASH name qmmm.py and preprocess.py pull in at import; none is called here.
ASH_NAMES = ("Fragment", "OpenMMTheory", "QMMMTheory", "OpenMM_MD", "OpenMM_Modeller",
             "Singlepoint", "xTBTheory", "PySCFTheory", "small_molecule_parameterizer",
             "merge_pdb_files", "actregiondefine", "read_intlist_from_file",
             "Optimizer", "write_pdbfile")


@pytest.fixture
def qmmm(monkeypatch):
    """`mcube.qmmm`, imported against a stub `ash` (nothing here calls into it)."""
    stub = types.ModuleType("ash")
    for name in ASH_NAMES:
        setattr(stub, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", stub)
    for name in [m for m in sys.modules if m.startswith("mcube")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    from mcube import qmmm as module
    return module


@pytest.fixture
def preprocess(qmmm):
    from mcube import preprocess as module
    return module


@pytest.fixture
def sp(qmmm):
    from mcube import sp as module
    return module


# Arg229 (3 atoms), Lys226 (2 atoms) and a zinc, in one chain.
ATOMS = [("A", "229", "N", "N"), ("A", "229", "CA", "C"), ("A", "229", "NH1", "N"),
         ("A", "226", "N", "N"), ("A", "226", "NZ", "N"),
         ("B", "1", "ZN", "Zn")]
CHARGES = [-0.4, 0.1, -0.8, -0.4, -0.3, 2.0]


def topology():
    top = Topology()
    chains, residues = {}, {}
    for chain_id, resid, name, symbol in ATOMS:
        chain = chains.setdefault(chain_id, top.addChain(chain_id))
        key = (chain_id, resid)
        if key not in residues:
            residues[key] = top.addResidue("RES", chain, id=resid)
        top.addAtom(name, element.Element.getBySymbol(symbol), residues[key])
    return top


def system_with_charges(charges=CHARGES):
    """A system whose only force is a NonbondedForce carrying `charges`."""
    system = openmm.System()
    nonbonded = openmm.NonbondedForce()
    for q in charges:
        system.addParticle(12.0 * unit.amu)
        nonbonded.addParticle(q, 0.3, 0.5)
    system.addForce(nonbonded)
    return system, nonbonded


class FakeMM:
    """Enough of ASH's OpenMMTheory for `apply_zero_charges`: the system, the
    topology, and ASH's own copy of the charges, which QMMMTheory reads."""

    def __init__(self, charges=CHARGES):
        self.system, self.nonbonded = system_with_charges(charges)
        self.topology = topology()
        self.charges = list(charges)
        self.charges_original = list(charges)

    def particle_charges(self):
        return [self.nonbonded.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
                for i in range(self.nonbonded.getNumParticles())]


# ---- resolve_residues: the selector, one field shorter than an atom selector ----

class TestResolveResidues:
    def test_a_residue_selector_gives_every_atom_of_that_residue(self, preprocess):
        assert preprocess.resolve_residues(["A:229", "A:226"], topology()) == \
            {"A:229": [0, 1, 2], "A:226": [3, 4]}

    def test_the_chain_is_part_of_the_selector(self, preprocess):
        assert preprocess.resolve_residues(["B:1"], topology()) == {"B:1": [5]}

    def test_an_atom_selector_is_refused(self, preprocess):
        with pytest.raises(ValueError, match="must be 'chain:resid'"):
            preprocess.resolve_residues(["A:229:NH1"], topology())

    def test_an_index_is_refused(self, preprocess):
        with pytest.raises(ValueError, match="must be 'chain:resid'"):
            preprocess.resolve_residues([229], topology())

    def test_a_residue_not_in_the_structure_is_refused(self, preprocess):
        with pytest.raises(ValueError, match="No residue 'A:90'"):
            preprocess.resolve_residues(["A:90"], topology())

    def test_a_number_used_twice_in_one_chain_is_refused(self, preprocess):
        """A selector that matches two residues does not say which, and zeroing
        both would act on atoms nobody named."""
        top = Topology()
        chain = top.addChain("A")
        for _ in range(2):           # two residue objects, both numbered 229
            res = top.addResidue("RES", chain, id="229")
            top.addAtom("N", element.Element.getBySymbol("N"), res)
        with pytest.raises(ValueError, match="matches 2 residues"):
            preprocess.resolve_residues(["A:229"], top)


# ---- charge_deletion: pure, so the bookkeeping is testable on its own ----

class TestChargeDeletion:
    def test_it_counts_the_atoms_and_the_charge_removed(self, qmmm):
        record = qmmm.charge_deletion(CHARGES, {"A:229": [0, 1, 2], "A:226": [3, 4]}, [5])
        assert record == {"key": "zero_charges", "scope": "qm_mm_and_mm_mm",
                          "residues": ["A:229", "A:226"], "atoms": 5,
                          "charge_removed": -1.8, "atom_indices": [0, 1, 2, 3, 4]}

    def test_the_record_says_which_of_the_two_keys_it_was(self, qmmm):
        record = qmmm.charge_deletion(CHARGES, {"A:229": [0, 1, 2]}, [5],
                                      "zero_charges_qm_only")
        assert record["key"] == "zero_charges_qm_only"
        assert record["scope"] == "qm_mm"

    def test_the_refusal_names_the_key_that_asked(self, qmmm):
        with pytest.raises(ValueError, match="zero_charges_qm_only: B:1"):
            qmmm.charge_deletion(CHARGES, {"B:1": [5]}, [5], "zero_charges_qm_only")

    def test_the_residues_are_recorded_in_the_order_the_config_gave(self, qmmm):
        record = qmmm.charge_deletion(CHARGES, {"A:226": [3, 4], "A:229": [0, 1, 2]}, [])
        assert record["residues"] == ["A:226", "A:229"]

    def test_a_qm_residue_raises_and_is_named(self, qmmm):
        with pytest.raises(ValueError) as excinfo:
            qmmm.charge_deletion(CHARGES, {"A:229": [0, 1, 2], "B:1": [5]}, [2, 5])
        message = str(excinfo.value)
        # Both offenders named: A:229 overlaps the QM region through one atom
        assert "A:229" in message and "B:1" in message
        assert "electrostatic embedding has already removed" in message

    def test_a_residue_outside_the_qm_region_is_fine(self, qmmm):
        record = qmmm.charge_deletion(CHARGES, {"A:226": [3, 4]}, [0, 1, 2, 5])
        assert record["atoms"] == 2 and record["charge_removed"] == -0.7


# ---- apply_zero_charges: the key absent, and the charges it zeroes ----

class TestApplyZeroCharges:
    def test_no_key_changes_nothing(self, qmmm):
        mm = FakeMM()
        assert qmmm.apply_zero_charges({}, mm, [5]) is None
        assert mm.particle_charges() == pytest.approx(CHARGES)
        assert mm.charges == CHARGES

    def test_an_empty_list_changes_nothing(self, qmmm):
        mm = FakeMM()
        assert qmmm.apply_zero_charges({"zero_charges": []}, mm, [5]) is None
        assert mm.particle_charges() == pytest.approx(CHARGES)

    def test_every_atom_of_a_named_residue_is_zeroed(self, qmmm):
        mm = FakeMM()
        record = qmmm.apply_zero_charges({"zero_charges": ["A:229"]}, mm, [5])
        assert mm.particle_charges() == pytest.approx([0.0, 0.0, 0.0, -0.4, -0.3, 2.0])
        assert record == {"key": "zero_charges", "scope": "qm_mm_and_mm_mm",
                          "residues": ["A:229"], "atoms": 3, "charge_removed": -1.1}

    def test_ashs_own_charge_list_is_zeroed_too(self, qmmm):
        # QMMMTheory hands mm.charges to the QM code as the point charges, so the
        # OpenMM system alone is not enough
        mm = FakeMM()
        qmmm.apply_zero_charges({"zero_charges": ["A:226"]}, mm, [5])
        assert mm.charges == [-0.4, 0.1, -0.8, 0.0, 0.0, 2.0]
        assert mm.charges_original == [-0.4, 0.1, -0.8, 0.0, 0.0, 2.0]

    def test_the_vdw_terms_are_untouched(self, qmmm):
        mm = FakeMM()
        qmmm.apply_zero_charges({"zero_charges": ["A:229"]}, mm, [5])
        for i in range(3):
            _, sigma, epsilon = mm.nonbonded.getParticleParameters(i)
            assert sigma.value_in_unit(unit.nanometer) == pytest.approx(0.3)
            assert epsilon.value_in_unit(unit.kilojoule_per_mole) == pytest.approx(0.5)

    def test_a_residue_in_the_qm_region_raises_before_anything_is_changed(self, qmmm):
        mm = FakeMM()
        with pytest.raises(ValueError, match="B:1"):
            qmmm.apply_zero_charges({"zero_charges": ["A:226", "B:1"]}, mm, [5])
        assert mm.particle_charges() == pytest.approx(CHARGES)

    def test_the_info_line_names_the_residues_and_the_charge(self, qmmm, capsys):
        qmmm.apply_zero_charges({"zero_charges": ["A:229", "A:226"]}, FakeMM(), [5])
        out = capsys.readouterr().out
        assert "[INFO] MM charges zeroed on A:229, A:226" in out
        assert "5 atoms, -1.8000 e removed" in out


# ---- the record reaches the stage's results ----

def test_the_stage_records_what_it_zeroed(qmmm, monkeypatch):
    """`results.zeroed_charges`, because an energy from a modified Hamiltonian
    has to say so itself. The stage is not run here - only the record it keeps."""
    result = {}
    zeroed = qmmm.apply_zero_charges({"zero_charges": ["A:229"]}, FakeMM(), [5])
    if zeroed is not None:
        result["zeroed_charges"] = zeroed
    assert result["zeroed_charges"]["residues"] == ["A:229"]
    assert result["zeroed_charges"]["atoms"] == 3
    assert "atom_indices" not in result["zeroed_charges"]


# ---- zero_charges_qm_only: the QM/MM interaction form ----

class FakeTheory:
    """Enough of ASH's QMMMTheory for the detection to have something to find:
    the full-system charge lists, the MM-subset list of point charges and the MM
    atom indices that index it, plus attributes that are not charges at all."""

    def __init__(self, mm, qm_atoms=(5,), subset=True):
        natoms = len(mm.charges)
        self.qmatoms = list(qm_atoms)
        self.mmatoms = [i for i in range(natoms) if i not in set(qm_atoms)]
        self.charges = mm.charges                       # ASH shares the MM theory's list
        self.charges_qmregionzeroed = [0.0 if i in set(qm_atoms) else q
                                       for i, q in enumerate(mm.charges)]
        if subset:
            self.pointcharges = [mm.charges[i] for i in self.mmatoms]
        # Decoys: per-atom, numeric, and not charges
        self.masses = [12.0] * natoms
        self.nothing = [0.0] * natoms
        self.coords = np.zeros((natoms, 3))
        self.printlevel = 1


class TestApplyZeroChargesQmOnly:
    def test_no_key_changes_nothing(self, qmmm):
        mm = FakeMM()
        assert qmmm.apply_zero_charges_qm_only({}, FakeTheory(mm), mm, [5]) == (None, [])

    def test_the_other_key_is_left_to_the_mm_system_path(self, qmmm):
        mm = FakeMM()
        record, found = qmmm.apply_zero_charges_qm_only(
            {"zero_charges": ["A:229"]}, FakeTheory(mm), mm, [5])
        assert (record, found) == (None, [])

    def test_the_mm_system_keeps_its_own_electrostatics(self, qmmm):
        """The whole point of this key: the OpenMM charges are untouched, so the
        residue's MM-MM interaction is still there and cannot contaminate a
        difference taken between two structures."""
        mm = FakeMM()
        qmmm.apply_zero_charges_qm_only({"zero_charges_qm_only": ["A:229"]},
                                        FakeTheory(mm), mm, [5])
        assert mm.particle_charges() == pytest.approx(CHARGES)

    def test_the_qm_visible_lists_are_zeroed(self, qmmm):
        mm = FakeMM()
        theory = FakeTheory(mm)
        record, found = qmmm.apply_zero_charges_qm_only(
            {"zero_charges_qm_only": ["A:229"]}, theory, mm, [5])
        assert theory.charges == [0.0, 0.0, 0.0, -0.4, -0.3, 2.0]
        assert theory.charges_qmregionzeroed == [0.0, 0.0, 0.0, -0.4, -0.3, 0.0]
        assert record["key"] == "zero_charges_qm_only" and record["scope"] == "qm_mm"
        assert record["atoms"] == 3 and record["charge_removed"] == -1.1
        assert len(found) >= 2

    def test_an_mm_subset_list_is_mapped_through_the_mm_indices(self, qmmm):
        """`pointcharges` is indexed by MM atom, not by system atom, so writing a
        system index into it would zero the wrong atom."""
        mm = FakeMM()
        theory = FakeTheory(mm)
        qmmm.apply_zero_charges_qm_only({"zero_charges_qm_only": ["A:226"]},
                                        theory, mm, [5])
        # A:226 is atoms 3 and 4, which are positions 3 and 4 of the five MM atoms
        assert theory.pointcharges == [-0.4, 0.1, -0.8, 0.0, 0.0]

    def test_what_is_not_a_charge_list_is_left_alone(self, qmmm):
        mm = FakeMM()
        theory = FakeTheory(mm)
        record, _ = qmmm.apply_zero_charges_qm_only(
            {"zero_charges_qm_only": ["A:229"]}, theory, mm, [5])
        assert theory.masses == [12.0] * 6
        assert theory.nothing == [0.0] * 6
        assert "theory.masses" not in record["pointcharge_lists"]
        assert "theory.nothing" not in record["pointcharge_lists"]

    def test_the_lists_it_wrote_are_recorded(self, qmmm):
        mm = FakeMM()
        record, _ = qmmm.apply_zero_charges_qm_only(
            {"zero_charges_qm_only": ["A:229"]}, FakeTheory(mm), mm, [5])
        assert "theory.charges_qmregionzeroed" in record["pointcharge_lists"]
        assert "theory.pointcharges" in record["pointcharge_lists"]

    def test_a_residue_in_the_qm_region_raises(self, qmmm):
        mm = FakeMM()
        theory = FakeTheory(mm)
        with pytest.raises(ValueError, match="zero_charges_qm_only: B:1"):
            qmmm.apply_zero_charges_qm_only({"zero_charges_qm_only": ["B:1"]},
                                            theory, mm, [5])
        assert theory.charges == CHARGES

    def test_finding_no_charge_list_raises_rather_than_deleting_nothing(self, qmmm):
        """A silent no-op would be written down as a contribution of zero."""
        mm = FakeMM()
        del mm.charges, mm.charges_original   # an ASH that keeps them somewhere else
        bare = types.SimpleNamespace(qmatoms=[5], masses=[12.0] * 6, printlevel=1)
        with pytest.raises(ValueError, match="none of this QMMMTheory's attributes"):
            qmmm.apply_zero_charges_qm_only({"zero_charges_qm_only": ["A:229"]},
                                            bare, mm, [5])

    def test_a_list_that_refuses_the_write_raises(self, qmmm):
        class Stubborn(list):
            def __setitem__(self, *a):      # a copy-on-write or read-only sequence
                pass

        mm = FakeMM()
        theory = types.SimpleNamespace(qmatoms=[5], charges=Stubborn(CHARGES))
        mm.charges = Stubborn(CHARGES)
        mm.charges_original = Stubborn(CHARGES)
        with pytest.raises(ValueError, match="did not take the zeros"):
            qmmm.apply_zero_charges_qm_only({"zero_charges_qm_only": ["A:229"]},
                                            theory, mm, [5])

    def test_the_info_line_says_what_it_did(self, qmmm, capsys):
        mm = FakeMM()
        qmmm.apply_zero_charges_qm_only({"zero_charges_qm_only": ["A:229", "A:226"]},
                                        FakeTheory(mm), mm, [5])
        out = capsys.readouterr().out
        assert "[INFO] QM-visible point charges zeroed on A:229, A:226" in out
        assert "5 atoms, -1.8000 e removed" in out
        assert "MM system keeps its own electrostatics" in out


class TestMutuallyExclusive:
    def test_both_keys_raise(self, qmmm):
        with pytest.raises(ValueError, match="are both set"):
            qmmm.zero_charge_request({"zero_charges": ["A:229"],
                                      "zero_charges_qm_only": ["A:229"]})

    def test_both_keys_raise_from_the_mm_system_path_too(self, qmmm):
        mm = FakeMM()
        with pytest.raises(ValueError, match="are both set"):
            qmmm.apply_zero_charges({"zero_charges": ["A:229"],
                                     "zero_charges_qm_only": ["A:226"]}, mm, [5])
        assert mm.particle_charges() == pytest.approx(CHARGES)

    def test_both_keys_raise_from_the_qm_only_path_too(self, qmmm):
        mm = FakeMM()
        with pytest.raises(ValueError, match="are both set"):
            qmmm.apply_zero_charges_qm_only({"zero_charges": ["A:229"],
                                             "zero_charges_qm_only": ["A:226"]},
                                            FakeTheory(mm), mm, [5])

    def test_one_of_them_empty_is_not_asking_for_both(self, qmmm):
        assert qmmm.zero_charge_request(
            {"zero_charges": [], "zero_charges_qm_only": ["A:229"]}) == \
            ("zero_charges_qm_only", ["A:229"])

    def test_neither_key_is_no_request(self, qmmm):
        assert qmmm.zero_charge_request({"charge": -1}) == (None, [])


class TestVerifiedAfterRun:
    def test_zeros_still_in_place_are_recorded_as_verified(self, qmmm):
        mm = FakeMM()
        record, found = qmmm.apply_zero_charges_qm_only(
            {"zero_charges_qm_only": ["A:229"]}, FakeTheory(mm), mm, [5])
        assert qmmm.confirm_zero_charges(found) == []
        qmmm.zero_charges_verified(record, found)
        assert record["verified_after_run"] is True
        assert "restored_lists" not in record

    def test_charges_put_back_during_the_run_are_recorded_and_warned_about(self, qmmm, capsys):
        """ASH re-deriving the point charges off the OpenMM system mid-run is the
        one failure this cannot rule out in advance, so it is read back."""
        mm = FakeMM()
        theory = FakeTheory(mm)
        record, found = qmmm.apply_zero_charges_qm_only(
            {"zero_charges_qm_only": ["A:229"]}, theory, mm, [5])
        capsys.readouterr()
        theory.charges_qmregionzeroed[1] = 0.1          # ASH put it back
        assert qmmm.confirm_zero_charges(found) == ["theory.charges_qmregionzeroed"]
        qmmm.zero_charges_verified(record, found)
        assert record["verified_after_run"] is False
        assert record["restored_lists"] == ["theory.charges_qmregionzeroed"]
        assert "[WARNING]" in capsys.readouterr().out

    def test_nothing_to_verify_leaves_the_record_alone(self, qmmm):
        record = {"key": "zero_charges", "scope": "qm_mm_and_mm_mm"}
        assert qmmm.zero_charges_verified(record, []) is record
        assert "verified_after_run" not in record


# ---- the key reaches the stages that read it ----

def test_sp_inherits_zero_charges_from_qmmm_inputs(sp, tmp_path):
    """The decomposition is run as `mcube sp` single points, so `zero_charges` on
    `qmmm.inputs` has to reach the sp stage like every other QM/MM input."""
    config = tmp_path / "config.yaml"
    config.write_text("qmmm:\n  inputs:\n    zero_charges: [A:229, A:226]\n"
                      "sp:\n  inputs:\n    sp_basis: def2-svp\n")
    assert sp.sp_inputs(str(config), "sp")["zero_charges"] == ["A:229", "A:226"]


def test_the_section_can_override_the_residues(sp, tmp_path):
    """One config, one section per residue: each writes its own `zeroed_charges`."""
    config = tmp_path / "config.yaml"
    config.write_text("qmmm:\n  inputs:\n    zero_charges: [A:229]\n"
                      "sp_arg90:\n  inputs:\n    zero_charges: [A:90]\n")
    assert sp.sp_inputs(str(config), "sp_arg90")["zero_charges"] == ["A:90"]
