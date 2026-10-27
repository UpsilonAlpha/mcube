"""qmmm resolves its starting structure visibly: pytest tests/test_qmmm_structure.py

The fallback to `md.results.lastframe` is a convenience for the usual
preprocess -> md -> qmmm chain, but it silently reached into another run's
directory once, and QM/MM started on the wrong structure converges and looks
fine. There is no `--pdb` to override it any more: a structure that is not in
the config is a result nothing records. Runs without ASH via the stub used in
tests/test_umbrella.py.
"""

import sys
import types

import pytest

pytest.importorskip("openmm", reason="qmmm reads structures through OpenMM")

ASH_NAMES = ("Fragment", "OpenMMTheory", "QMMMTheory", "xTBTheory", "PySCFTheory",
             "Singlepoint", "Optimizer", "NumFreq", "NEB", "NEBTS", "OpenMM_Opt",
             "OpenMM_MD", "OpenMM_Modeller", "OpenMM_box_equilibration",
             "Gentle_warm_up_MD", "small_molecule_parameterizer", "merge_pdb_files",
             "read_intlist_from_file", "actregiondefine", "getwaterconstraintslist")


@pytest.fixture
def qmmm(monkeypatch):
    stub = types.ModuleType("ash")
    for name in ASH_NAMES:
        setattr(stub, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", stub)
    # `mcube` itself too: `from mcube import qmmm` would otherwise hand back the
    # stale attribute on the cached package, imported under another file's stub.
    for name in [m for m in sys.modules if m == "mcube" or m.startswith("mcube.")]:
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


def resolve(qmmm, home, inputs=None, md=None):
    return qmmm.resolve_start_structure(
        home.config, "qmmm", inputs or {},
        [((md or {}).get("lastframe"), "md.results.lastframe")])


class TestPrecedence:
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

    def test_nothing_at_all_names_both_ways(self, qmmm, home):
        with pytest.raises(KeyError) as exc:
            resolve(qmmm, home)
        text = str(exc.value)
        assert "`qmmm.inputs`" in text               # the key to set
        assert "md.results.lastframe" in text        # and the stage that writes one

    def test_a_relative_pdb_file_is_read_against_the_config(self, qmmm, home, monkeypatch, tmp_path):
        """The trap: the job script runs each stage from its own output directory,
        so a relative `pdb_file` used to mean a different file per stage."""
        monkeypatch.chdir(tmp_path / "run_b")
        path, source = resolve(qmmm, home, inputs={"pdb_file": "mine.pdb"})
        assert source == "qmmm.inputs.pdb_file"
        assert path == str(home.mine).replace("\\", "/")

    def test_a_named_but_missing_file_says_which_input_named_it(self, qmmm, home):
        with pytest.raises(FileNotFoundError, match="qmmm.inputs.pdb_file"):
            resolve(qmmm, home, inputs={"pdb_file": home.mine + ".nope"})


class TestActiveRegionContainsTheQMRegion:
    """A QM atom outside the active region is frozen and still treated quantum
    mechanically, and nothing in the output says so. The union is unconditional."""

    QM = [10, 11, 12, 900]          # 900 is the far end of a QM residue

    def actatoms(self, qmmm, monkeypatch, written, inputs=None, config="config.yaml"):
        monkeypatch.setattr(qmmm, "read_intlist_from_file", lambda *a: list(written))
        monkeypatch.setattr(qmmm, "actregiondefine", lambda **k: None)
        return qmmm.resolve_actatoms(config, inputs or {}, None, None, self.QM, [])

    def test_a_qm_atom_outside_the_radius_is_added(self, qmmm, monkeypatch):
        got = self.actatoms(qmmm, monkeypatch, [10, 11, 12, 13], {"actradius": 8, "actregion_origin": "B:1:ZN"})
        assert got == [10, 11, 12, 13, 900]

    def test_it_is_a_no_op_when_the_radius_already_covers_the_qm_region(self, qmmm, monkeypatch):
        got = self.actatoms(qmmm, monkeypatch, [10, 11, 12, 13, 900], {"actradius": 8, "actregion_origin": "B:1:ZN"})
        assert got == [10, 11, 12, 13, 900]

    def test_a_region_given_in_the_config_is_completed_too(self, qmmm, monkeypatch):
        """A reactant's region reused for the product must still cover QM."""
        got = self.actatoms(qmmm, monkeypatch, [], {"actatoms": "10-12"})
        assert got == [10, 11, 12, 900]

    def test_it_says_so_when_it_adds_any(self, qmmm, monkeypatch, capsys):
        self.actatoms(qmmm, monkeypatch, [10, 11, 12], {"actradius": 8, "actregion_origin": "B:1:ZN"})
        out = capsys.readouterr().out
        assert "fell outside the active region" in out and "900" in out

    def test_it_stays_quiet_when_it_adds_none(self, qmmm, monkeypatch, capsys):
        self.actatoms(qmmm, monkeypatch, list(self.QM), {"actradius": 8, "actregion_origin": "B:1:ZN"})
        assert "fell outside" not in capsys.readouterr().out

    def test_no_radius_and_no_file_is_still_just_the_qm_atoms(self, qmmm, monkeypatch):
        assert self.actatoms(qmmm, monkeypatch, [1, 2, 3]) == self.QM

    def test_actradius_without_an_origin_or_boundary_residue_is_refused(self, qmmm, monkeypatch):
        monkeypatch.setattr(qmmm, "read_intlist_from_file", lambda *a: [10])
        monkeypatch.setattr(qmmm, "actregiondefine", lambda **k: None)
        with pytest.raises(ValueError, match="no actregion_origin"):
            qmmm.resolve_actatoms("config.yaml", {"actradius": 8}, None, None, self.QM, [])
