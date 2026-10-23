"""xTB's library runmode is refused: pytest tests/test_xtb_runmode.py

mcube only ever runs xTB as the QM region of a QMMMTheory, and ASH's library
path does not implement QM/MM point charges. Choosing it would drop the
electrostatic embedding rather than fail, so it must be refused outright.
Runs without ASH via the stub used in tests/test_umbrella.py.
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


def test_library_from_the_config_is_refused(qmmm):
    with pytest.raises(ValueError, match="point charges"):
        qmmm.make_qm({"xtb_runmode": "library"}, 4)


def test_library_passed_by_a_stage_is_refused(qmmm):
    """neb and ts used to force this; the guard catches it wherever it comes from."""
    with pytest.raises(ValueError, match="point charges"):
        qmmm.make_qm({}, 4, runmode="library")


def test_case_does_not_evade_the_guard(qmmm):
    with pytest.raises(ValueError, match="point charges"):
        qmmm.make_qm({"xtb_runmode": "Library"}, 4)


def test_the_message_says_what_to_use_instead(qmmm):
    with pytest.raises(ValueError) as exc:
        qmmm.make_qm({}, 4, runmode="library")
    assert "inputfile" in str(exc.value)


def test_inputfile_is_accepted(qmmm):
    """The stub xTBTheory returns None; reaching it at all means no refusal."""
    assert qmmm.make_qm({"xtb_runmode": "inputfile"}, 4) is None
    assert qmmm.make_qm({}, 4) is None            # the default


def test_dft_is_unaffected(qmmm):
    """The guard is about xTB; it must not interfere with the DFT engine."""
    assert qmmm.qm_level_label({"qm_engine": "dft"}) == "wb97x-v/def2-svp"
