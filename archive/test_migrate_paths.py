"""Path migration: pytest tests/test_migrate_paths.py

What this must get right is not the rewriting, which is mechanical, but the
one case it must refuse: an absolute path into *another* run. Three CSR
configs are in that state, inherited by copying, and rewriting one to a bare
name would silently change which force field the run uses. That is a decision
about the science, so the script reports it and leaves it alone.
"""

import importlib.util

import pytest

from pathlib import Path  # noqa: E402

SCRIPT = Path(__file__).resolve().parents[1] / "docs" / "csr" / "migrate_paths.py"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("migrate_paths", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def run(tmp_path):
    """A run directory with its own files, beside a second run."""
    home = tmp_path / "CSR_nac"
    (home / "preprocess_results").mkdir(parents=True)
    (home / "hydroxide.xml").write_text("<ForceField/>")
    (home / "preprocess_results" / "openff_MER.xml").write_text("<ForceField/>")
    other = tmp_path / "CSR_crystal"
    other.mkdir()
    (other / "hydroxide.xml").write_text("<ForceField/>")
    return home, other


def write(home, entries, key="forcefields"):
    body = ["md:\n", "  inputs:\n", f"    {key}:\n"]
    body += [f"    - {e}\n" for e in entries]
    body += ["    cores: 4\n"]
    config = home / "config.yaml"
    config.write_text("".join(body))
    return config


def entries_of(config):
    return [l.strip()[2:] for l in config.read_text().splitlines()
            if l.strip().startswith("- ")]


def test_a_dotdot_entry_into_this_run_is_rewritten(mod, run):
    """`../hydroxide.xml` only ever worked because the stage ran one level
    down; under the new rule it points outside the run entirely."""
    home, _ = run
    config = write(home, ["amber14/tip3p.xml", "../hydroxide.xml"])
    # A decoy where home/".."/value literally points. The entry never meant
    # this file, and resolving it there is the hazard the rule avoids.
    (home.parent / "hydroxide.xml").write_text("<DECOY/>")
    mod.migrate(config, apply=True)
    assert entries_of(config) == ["amber14/tip3p.xml", "hydroxide.xml"]
    assert (home / "hydroxide.xml").read_text() == "<ForceField/>"


def test_an_absolute_path_into_this_run_becomes_relative(mod, run):
    home, _ = run
    config = write(home, [str(home / "hydroxide.xml")])
    mod.migrate(config, apply=True)
    assert entries_of(config) == ["hydroxide.xml"]


def test_an_absolute_path_into_another_run_is_left_alone(mod, run, capsys):
    """The case that matters. Rewriting it to a bare name would switch the run
    to a different file - which is what copying a config quietly did."""
    home, other = run
    foreign = str(other / "hydroxide.xml")
    config = write(home, [foreign])
    mod.migrate(config, apply=True)
    assert entries_of(config) == [foreign], "a foreign path was rewritten"
    out = capsys.readouterr().out
    assert "foreign" in out and "ANOTHER run" in out


def test_openmm_names_are_untouched(mod, run):
    home, _ = run
    names = ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]
    config = write(home, list(names))
    mod.migrate(config, apply=True)
    assert entries_of(config) == names


def test_an_already_relative_entry_is_untouched(mod, run):
    home, _ = run
    config = write(home, ["preprocess_results/openff_MER.xml"])
    mod.migrate(config, apply=True)
    assert entries_of(config) == ["preprocess_results/openff_MER.xml"]


def test_a_dry_run_writes_nothing(mod, run):
    home, _ = run
    config = write(home, [str(home / "hydroxide.xml")])
    before = config.read_text()
    mod.migrate(config, apply=False)
    assert config.read_text() == before
    assert not config.with_suffix(".yaml.bak_paths").exists()


def test_applying_leaves_a_backup(mod, run):
    home, _ = run
    config = write(home, [str(home / "hydroxide.xml")])
    before = config.read_text()
    mod.migrate(config, apply=True)
    assert config.with_suffix(".yaml.bak_paths").read_text() == before


def test_xml_outside_a_forcefield_list_is_not_touched(mod, run):
    """Only force-field lists. A `.xml` under some other key is not ours."""
    home, _ = run
    config = home / "config.yaml"
    config.write_text(
        "md:\n  inputs:\n    ligand_xml:\n    - %s\n" % (home / "hydroxide.xml"))
    mod.migrate(config, apply=True)
    assert entries_of(config) == [str(home / "hydroxide.xml")]


class TestRepointing:
    """--repoint sends a config at its own copy of a force field, but only
    when that is a no-op. Five CSR configs read force fields out of another
    run, inherited by copying; for one of them the local copy has different
    charges, so repointing it would change the force field the run used while
    looking like a path fix."""

    def test_an_identical_local_copy_is_used(self, mod, run, tmp_path):
        home, other = run
        mod.REPOINT = True
        try:
            config = write(home, [str(other / "hydroxide.xml")])
            mod.migrate(config, apply=True)
            assert entries_of(config) == ["hydroxide.xml"]
        finally:
            mod.REPOINT = False

    def test_a_copy_differing_only_by_date_is_used(self, mod, run, tmp_path):
        """OpenFF stamps its output with the generation date, so two runs of
        the same preprocess differ in that one line and nothing else."""
        home, other = run
        (home / "preprocess_results" / "openff.xml").write_text(
            "<ForceField>\n  <DateGenerated>2026-10-21</DateGenerated>\n  <x/>\n</ForceField>")
        (other / "openff.xml").write_text(
            "<ForceField>\n  <DateGenerated>2026-10-15</DateGenerated>\n  <x/>\n</ForceField>")
        mod.REPOINT = True
        try:
            config = write(home, [str(other / "openff.xml")])
            mod.migrate(config, apply=True)
            assert entries_of(config) == ["preprocess_results/openff.xml"]
        finally:
            mod.REPOINT = False

    def test_a_copy_with_different_parameters_is_refused(self, mod, run, capsys):
        """The case that matters: CSR_unrestrained2's own openff_MER.xml has
        different atom charges from the CSR_crystal one it reads. Repointing
        would change the force field, not the path."""
        home, other = run
        (home / "preprocess_results" / "openff.xml").write_text(
            '<ForceField>\n  <Atom name="C1x" charge="-0.092"/>\n</ForceField>')
        (other / "openff.xml").write_text(
            '<ForceField>\n  <Atom name="C1x" charge="-0.144"/>\n</ForceField>')
        mod.REPOINT = True
        try:
            foreign = str(other / "openff.xml")
            config = write(home, [foreign])
            mod.migrate(config, apply=True)
            assert entries_of(config) == [foreign], "parameters were changed"
            assert "DIFFERENT contents" in capsys.readouterr().out
        finally:
            mod.REPOINT = False

    def test_a_run_without_its_own_copy_is_refused(self, mod, run, capsys):
        home, other = run
        (other / "openff.xml").write_text("<ForceField/>")
        mod.REPOINT = True
        try:
            foreign = str(other / "openff.xml")
            config = write(home, [foreign])
            mod.migrate(config, apply=True)
            assert entries_of(config) == [foreign]
            assert "no openff.xml in this run" in capsys.readouterr().out
        finally:
            mod.REPOINT = False

    def test_without_the_flag_nothing_is_repointed(self, mod, run):
        home, other = run
        foreign = str(other / "hydroxide.xml")
        config = write(home, [foreign])
        mod.migrate(config, apply=True)
        assert entries_of(config) == [foreign]
