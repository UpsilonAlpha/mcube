"""qmmm resolves its starting structure visibly: pytest tests/test_qmmm_structure.py

The fallback to `md.results.lastframe` is a convenience for the usual
preprocess -> md -> qmmm chain, but it silently reached into another run's
directory once, and QM/MM started on the wrong structure converges and looks
fine. Runs without ASH via the stub used in tests/test_umbrella.py.
"""

import sys
import types

import pytest

pytest.importorskip("openmm", reason="qmmm reads structures through OpenMM")

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


@pytest.fixture
def home(tmp_path):
    """A config directory with a frame of its own, and a frame in a sibling run."""
    here = tmp_path / "run_a"
    here.mkdir()
    (here / "config.yaml").write_text("qmmm:\n  inputs: {}\n")
    (here / "mine.pdb").write_text("END\n")
    other = tmp_path / "run_b"
    other.mkdir()
    (other / "theirs.pdb").write_text("END\n")
    return types.SimpleNamespace(
        config=str(here / "config.yaml"),
        mine=str(here / "mine.pdb"),
        theirs=str(other / "theirs.pdb"),
    )


def resolve(qmmm, home, inputs=None, md=None, pdb=None):
    return qmmm.resolve_start_structure(home.config, "qmmm", inputs or {}, md or {}, pdb)


class TestPrecedence:
    def test_pdb_argument_wins(self, qmmm, home):
        path, source = resolve(qmmm, home, inputs={"pdb_file": home.mine},
                               md={"lastframe": home.mine}, pdb=home.mine)
        assert source == "--pdb" and path.endswith("mine.pdb")

    def test_section_input_beats_the_md_fallback(self, qmmm, home):
        _, source = resolve(qmmm, home, inputs={"pdb_file": home.mine},
                            md={"lastframe": home.theirs})
        assert source == "qmmm.inputs.pdb_file"

    def test_md_lastframe_is_used_inside_the_config_directory(self, qmmm, home):
        _, source = resolve(qmmm, home, md={"lastframe": home.mine})
        assert source == "md.results.lastframe"

    def test_the_path_comes_back_absolute(self, qmmm, home):
        path, _ = resolve(qmmm, home, md={"lastframe": home.mine})
        assert path == str(home.mine).replace("\\", "/")


class TestFailuresAreClear:
    def test_md_lastframe_outside_the_config_directory_is_refused(self, qmmm, home):
        """The real bug: a frame from a different run, taken silently."""
        with pytest.raises(ValueError) as exc:
            resolve(qmmm, home, md={"lastframe": home.theirs})
        text = str(exc.value)
        assert "outside this config's directory" in text
        assert "theirs.pdb" in text              # names the offending frame
        assert "`qmmm.inputs`" in text           # and how to fix it

    def test_an_explicit_path_outside_the_directory_is_allowed(self, qmmm, home):
        """Only the silent fallback is refused; naming it is a deliberate choice."""
        _, source = resolve(qmmm, home, inputs={"pdb_file": home.theirs})
        assert source == "qmmm.inputs.pdb_file"

    def test_nothing_at_all_names_the_three_ways(self, qmmm, home):
        with pytest.raises(KeyError) as exc:
            resolve(qmmm, home)
        text = str(exc.value)
        assert "pdb_file" in text and "--pdb" in text and "md" in text

    def test_a_named_but_missing_file_says_which_input_named_it(self, qmmm, home):
        with pytest.raises(FileNotFoundError, match="qmmm.inputs.pdb_file"):
            resolve(qmmm, home, inputs={"pdb_file": home.mine + ".nope"})
