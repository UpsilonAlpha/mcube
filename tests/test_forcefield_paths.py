"""Force-field paths read against the config: pytest tests/test_forcefield_paths.py

A run directory is meant to be copyable. Absolute force-field paths are not:
three of the CSR configs were reading `openff_MER.xml` and `hydroxide.xml` out
of a *different* run's directory, inherited when the config was copied, and
nothing complained because the files existed. Relative-to-the-stage paths were
no better - `../hydroxide.xml` meant different files from `md_x_results/` and
`protdft/x/`, which is what put a substrate carbon where the zinc should have
been (D115's sibling).

So everything resolves against the config's directory. What these tests pin is
the three things that rule must not break: OpenMM's bundled names, absolute
paths someone wrote deliberately, and the `..` entries that would otherwise be
silently repointed at the scratch root.
"""

import sys
import types
from pathlib import Path

import pytest

from tests.test_qmmm_structure import ASH_NAMES


@pytest.fixture
def resolve_forcefields(monkeypatch):
    """`preprocess` imports ASH at module scope, which is not installed here."""
    stub = types.ModuleType("ash")
    for name in ASH_NAMES:
        setattr(stub, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", stub)
    for name in [m for m in sys.modules if m == "mcube" or m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    from mcube.preprocess import resolve_forcefields as fn
    return fn


@pytest.fixture
def home(tmp_path):
    """A config directory with two force-field files beside the config."""
    (tmp_path / "hydroxide.xml").write_text("<ForceField/>")
    (tmp_path / "water_nuc.xml").write_text("<ForceField/>")
    (tmp_path / "preprocess_results").mkdir()
    (tmp_path / "preprocess_results" / "openff_MER.xml").write_text("<ForceField/>")
    config = tmp_path / "config.yaml"
    config.write_text("preprocess: {}\n")
    return config


def test_a_file_beside_the_config_is_resolved(resolve_forcefields, home):
    out = resolve_forcefields(["hydroxide.xml"], home)
    assert out == [(home.parent / "hydroxide.xml").resolve().as_posix()]


def test_a_file_in_a_subdirectory_is_resolved(resolve_forcefields, home):
    out = resolve_forcefields(["preprocess_results/openff_MER.xml"], home)
    assert out[0].endswith("preprocess_results/openff_MER.xml")
    assert out[0].startswith(home.parent.resolve().as_posix())


def test_openmms_bundled_names_are_left_alone(resolve_forcefields, home):
    """`amber14/protein.ff14SB.xml` is not a file beside the config; OpenMM
    resolves it from its own data directory and must keep seeing it verbatim."""
    names = ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]
    assert resolve_forcefields(list(names), home) == names


def test_an_absolute_path_is_left_alone(resolve_forcefields, home):
    """It already says what it means. Rewriting it would change which file is
    used, which is a decision, not a path fix."""
    other = "/scratch/user/someone/CSR_crystal/hydroxide.xml"
    assert resolve_forcefields([other], home) == [other]


def test_a_dotdot_entry_is_left_alone_and_reported(resolve_forcefields, home, capsys):
    """The load-bearing exclusion. `<config dir>/../water_nuc.xml` is the
    scratch root, where spare copies accumulate: resolving it would silently
    prefer a stranger over the file the stage reads today. It is reported
    instead, because under this rule it no longer means what it did."""
    (home.parent.parent / "water_nuc.xml").write_text("<ForceField/>")
    out = resolve_forcefields(["../water_nuc.xml"], home)
    assert out == ["../water_nuc.xml"], "a .. entry was rewritten"
    printed = capsys.readouterr().out
    assert "../water_nuc.xml" in printed and "drop the '../'" in printed


def test_a_missing_file_is_passed_through_not_invented(resolve_forcefields, home):
    """Anything unrecognised goes to OpenMM, which raises its own error naming
    the file. Inventing a path here would bury that."""
    assert resolve_forcefields(["nosuch.xml"], home) == ["nosuch.xml"]


def test_a_mixed_list_keeps_its_order(resolve_forcefields, home):
    """The order of a force-field list is meaningful: later files override
    earlier ones, which is how water_nuc.xml replaces the hydroxide."""
    out = resolve_forcefields(
        ["amber14/protein.ff14SB.xml", "hydroxide.xml", "/abs/x.xml"], home)
    assert out[0] == "amber14/protein.ff14SB.xml"
    assert out[1].endswith("/hydroxide.xml")
    assert Path(out[1]).is_absolute(), out[1]
    assert out[2] == "/abs/x.xml"


def test_empty_and_none_survive(resolve_forcefields, home):
    assert resolve_forcefields([], home) == []
    assert resolve_forcefields(None, home) is None


def test_resolution_does_not_depend_on_the_working_directory(resolve_forcefields, home, tmp_path, monkeypatch):
    """The whole point. The same config must name the same files whichever
    directory a stage happens to run in - that is what `../` could not do."""
    deep = home.parent / "protdft" / "qmmm_x"
    deep.mkdir(parents=True)
    monkeypatch.chdir(deep)
    from_deep = resolve_forcefields(["hydroxide.xml"], home)
    monkeypatch.chdir(tmp_path)
    from_top = resolve_forcefields(["hydroxide.xml"], home)
    assert from_deep == from_top
