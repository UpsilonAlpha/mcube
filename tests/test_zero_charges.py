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
        assert record == {"residues": ["A:229", "A:226"], "atoms": 5,
                          "charge_removed": -1.8, "atom_indices": [0, 1, 2, 3, 4]}

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
        assert record == {"residues": ["A:229"], "atoms": 3, "charge_removed": -1.1}

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
