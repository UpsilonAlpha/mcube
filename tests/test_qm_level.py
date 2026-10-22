"""Every stage records which QM method produced its energies.

A barrier assembled from opt_energy, mep_energies or saddlepoint_energy carries
no provenance unless the level is recorded, and GFN2-xTB energies are not
barriers. Runs without ASH via the stub used in tests/test_umbrella.py.
"""

import sys
import types

import pytest

ASH_NAMES = ("Fragment", "OpenMMTheory", "QMMMTheory", "xTBTheory", "PySCFTheory",
             "Singlepoint", "Optimizer", "NumFreq", "NEB", "NEBTS", "OpenMM_Opt",
             "OpenMM_MD", "OpenMM_Modeller", "OpenMM_box_equilibration",
             "Gentle_warm_up_MD", "small_molecule_parameterizer", "merge_pdb_files")


@pytest.fixture
def qmmm(monkeypatch):
    stub = types.ModuleType("ash")
    for name in ASH_NAMES:
        setattr(stub, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", stub)
    for name in [m for m in sys.modules if m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    from mcube import qmmm
    return qmmm


class TestQmLevelLabel:
    def test_default_is_gfn2_xtb(self, qmmm):
        assert qmmm.qm_level_label({}) == "GFN2-xTB"

    def test_xtb_method_is_honoured(self, qmmm):
        assert qmmm.qm_level_label({"xtb_method": "GFN1"}) == "GFN1-xTB"

    def test_dft_reports_functional_and_basis(self, qmmm):
        assert qmmm.qm_level_label(
            {"qm_engine": "dft", "dft_functional": "b3lyp", "dft_basis": "def2-tzvp"}
        ) == "b3lyp/def2-tzvp"

    def test_dft_falls_back_to_the_sp_keys(self, qmmm):
        """dft_* falls back to sp_*, so a level set once is not silently dropped."""
        assert qmmm.qm_level_label(
            {"qm_engine": "dft", "sp_functional": "wb97m-v", "sp_basis": "def2-tzvp"}
        ) == "wb97m-v/def2-tzvp"

    def test_dft_defaults_are_named_not_blank(self, qmmm):
        assert qmmm.qm_level_label({"qm_engine": "dft"}) == "wb97x-v/def2-svp"

    def test_pyscf_is_treated_as_dft(self, qmmm):
        assert qmmm.qm_level_label({"qm_engine": "pyscf"}) == "wb97x-v/def2-svp"

    def test_an_xtb_label_is_never_mistaken_for_a_dft_one(self, qmmm):
        """The point of the key: the two are distinguishable in the results."""
        assert "/" not in qmmm.qm_level_label({})
        assert "/" in qmmm.qm_level_label({"qm_engine": "dft"})
