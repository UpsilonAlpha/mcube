"""sections.py checks the structure a section names: pytest tests/test_sections_structures.py

The selector check (a scan whose coordinate starts beyond hydrogen-bond distance
is not a proton transfer) needs a structure to check against, and that structure
used to be a separate `start` key that mcube never read - so the thing checked
here and the thing computed on could differ. Both are `inputs.pdb_file` now, and
the two NEB endpoints are `reactant_pdb`/`product_pdb`.

Needs only PyYAML.
"""

import importlib.util
import sys
import types

import pytest

yaml = pytest.importorskip("yaml")
pytest.importorskip("openmm", reason="preprocess reads structures through OpenMM")

# sections.py uses mcube.preprocess for the config lock, and that imports ASH at
# module scope; nothing here calls into it. Same stub as tests/test_umbrella.py.
ASH_NAMES = ("Fragment", "OpenMM_Modeller", "small_molecule_parameterizer",
             "merge_pdb_files")

from pathlib import Path  # noqa: E402

SCRIPT = Path(__file__).resolve().parents[1] / "docs" / "csr" / "sections.py"

# One zinc and one histidine nitrogen: enough for a selector to resolve or not.
PDB = ("HETATM    1 ZN   ZN  B   1       0.000   0.000   0.000\n"
       "ATOM      2  NE2 HIS A  88       0.000   0.000   2.100\n")


@pytest.fixture
def mod(monkeypatch):
    stub = types.ModuleType("ash")
    for name in ASH_NAMES:
        setattr(stub, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", stub)
    for name in [m for m in sys.modules if m == "mcube" or m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    spec = importlib.util.spec_from_file_location("sections", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def home(tmp_path, monkeypatch):
    """A config directory with one structure, run from somewhere else entirely."""
    (tmp_path / "md_results").mkdir()
    (tmp_path / "md_results" / "last.pdb").write_text(PDB)
    (tmp_path / "config.yaml").write_text("preprocess:\n  inputs: {}\n")
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(tmp_path)        # sections.py writes ./config.yaml
    return tmp_path


def site(home, body):
    path = home / "site.yaml"
    path.write_text(yaml.safe_dump({"sections": body}, sort_keys=False))
    return str(path)


def written(home, section):
    data = yaml.safe_load(open(home / "config.yaml"))
    return (data.get(section) or {}).get("inputs")


class TestTheStructureIsAnInput:
    def test_a_section_is_written_when_its_selectors_resolve(self, mod, home):
        file = site(home, {"znhis": {"stage": "scan", "inputs": {
            "pdb_file": "md_results/last.pdb", "atoms": ["B:1:ZN", "A:88:NE2"]}}})
        mod.main(file, ["znhis"])
        assert written(home, "scan_znhis")["atoms"] == ["B:1:ZN", "A:88:NE2"]

    def test_a_selector_that_is_not_in_the_structure_stops_everything(self, mod, home):
        """Nothing is written, so a bad section cannot be half-added."""
        file = site(home, {"znhis": {"stage": "scan", "inputs": {
            "pdb_file": "md_results/last.pdb", "atoms": ["B:1:ZN", "A:226:NZ"]}}})
        with pytest.raises(SystemExit):
            mod.main(file, ["znhis"])
        assert written(home, "scan_znhis") is None

    def test_a_missing_structure_is_named_with_its_key(self, mod, home, capsys):
        file = site(home, {"znhis": {"stage": "scan",
                                     "inputs": {"pdb_file": "md_results/gone.pdb"}}})
        with pytest.raises(SystemExit):
            mod.main(file, ["znhis"])
        out = capsys.readouterr().out
        assert "pdb_file missing" in out and "gone.pdb" in out

    def test_a_path_is_read_against_the_config_not_the_shell(self, mod, home, monkeypatch):
        """mcube reads it that way, so the check has to look in the same place."""
        monkeypatch.chdir(home / "elsewhere")
        monkeypatch.setattr(mod, "Path", mod.Path)       # no-op; keeps the import visible
        file = site(home, {"znhis": {"stage": "scan", "inputs": {
            "pdb_file": "md_results/last.pdb", "atoms": ["B:1:ZN"]}}})
        monkeypatch.chdir(home)
        mod.main(file, ["znhis"])
        assert written(home, "scan_znhis")["pdb_file"] == "md_results/last.pdb"


class TestRefusals:
    def test_a_leftover_start_key_says_where_it_went(self, mod, home, capsys):
        """`start` described a structure mcube never read; leaving it silently
        unused would check one structure and compute on another."""
        file = site(home, {"znhis": {"stage": "scan",
                                     "start": "md_results/last.pdb",
                                     "inputs": {"atoms": ["B:1:ZN"]}}})
        with pytest.raises(SystemExit):
            mod.main(file, ["znhis"])
        out = capsys.readouterr().out
        assert "`start` is gone" in out and "inputs.pdb_file" in out

    def test_a_neb_section_needs_both_endpoints(self, mod, home, capsys):
        file = site(home, {"step1": {"stage": "neb", "inputs": {
            "reactant_pdb": "md_results/last.pdb"}}})
        with pytest.raises(SystemExit):
            mod.main(file, ["step1"])
        out = capsys.readouterr().out
        assert "reactant_pdb" in out and "product_pdb" in out

    def test_a_neb_with_both_endpoints_checks_the_reactant(self, mod, home):
        """Endpoints are one system in two geometries, so the reactant's atom
        order is the one the selectors have to resolve in."""
        file = site(home, {"step1": {"stage": "neb", "inputs": {
            "reactant_pdb": "md_results/last.pdb",
            "product_pdb": "md_results/last.pdb",
            "hessian_atoms": ["B:1:ZN"]}}})
        mod.main(file, ["step1"])
        assert written(home, "neb_step1")["product_pdb"] == "md_results/last.pdb"


class TestTheBaseIsTheConfigNotTheSiteFile:
    """site.yaml lives in the repo; the structures live in the run directory.

    The existing path test could not tell the two apart, because its helper
    writes site.yaml into the same directory as config.yaml. In reality they
    are far apart - `python $M/docs/csr/sections.py $M/docs/csr/site.yaml ohHmin`
    run from /scratch/.../CSR_nac - and resolving against site.yaml's directory
    looked for the run's frames inside the git checkout, refusing every section
    with "pdb_file missing" for a file that was plainly there.
    """

    def test_a_site_file_outside_the_run_still_finds_the_structure(
            self, mod, home, tmp_path):
        repo = tmp_path.parent / "repo_docs_csr"
        repo.mkdir(exist_ok=True)
        path = repo / "site.yaml"
        path.write_text(yaml.safe_dump({"sections": {"znhis": {
            "stage": "scan",
            "inputs": {"pdb_file": "md_results/last.pdb",
                       "atoms": ["B:1:ZN"]}}}}, sort_keys=False))
        mod.main(str(path), ["znhis"])
        assert written(home, "scan_znhis")["pdb_file"] == "md_results/last.pdb"

    def test_a_structure_beside_the_site_file_does_not_count(
            self, mod, home, tmp_path, capsys):
        """The mirror image: a file next to site.yaml is in the repo, not in
        the run, and must not satisfy the check."""
        repo = tmp_path.parent / "repo_decoy"
        (repo / "md_results").mkdir(parents=True, exist_ok=True)
        (repo / "md_results" / "decoy.pdb").write_text(PDB)
        path = repo / "site.yaml"
        path.write_text(yaml.safe_dump({"sections": {"znhis": {
            "stage": "scan",
            "inputs": {"pdb_file": "md_results/decoy.pdb"}}}}, sort_keys=False))
        with pytest.raises(SystemExit):
            mod.main(str(path), ["znhis"])
        assert "decoy.pdb" in capsys.readouterr().out


class TestRestraintSelectorsAreChecked:
    """D128: `restraints` was not in SELECTOR_KEYS, so an `md` section whose only
    selectors were restraint pairs reported `0 selectors checked` and was written
    without any of them being resolved. `md.md_restraints` puts that same list
    through `resolve_atoms`, so a name that does not resolve would surface hours
    into the run instead of at the point the section is added."""

    def test_a_restraint_pair_that_resolves_is_counted(self, mod, home, capsys):
        file = site(home, {"sub": {"stage": "md", "inputs": {
            "pdb_file": "md_results/last.pdb",
            "restraints": [["B:1:ZN", "A:88:NE2", 2.2]]}}})
        mod.main(file, ["sub"])
        assert written(home, "md_sub")["restraints"] == [["B:1:ZN", "A:88:NE2", 2.2]]
        assert "0 selectors checked" not in capsys.readouterr().out

    def test_a_restraint_selector_that_is_absent_stops_everything(self, mod, home):
        file = site(home, {"sub": {"stage": "md", "inputs": {
            "pdb_file": "md_results/last.pdb",
            "restraints": [["B:1:ZN", "A:226:NZ", 2.2]]}}})
        with pytest.raises(SystemExit):
            mod.main(file, ["sub"])
        assert written(home, "md_sub") is None

    def test_the_distance_is_not_mistaken_for_a_selector(self, mod, home, capsys):
        """The third element is a number; only strings with two colons count."""
        file = site(home, {"sub": {"stage": "md", "inputs": {
            "pdb_file": "md_results/last.pdb",
            "restraints": [["B:1:ZN", "A:88:NE2", 2.2]]}}})
        mod.main(file, ["sub"])
        assert "2 selectors checked" in capsys.readouterr().out

class TestADoubledStagePrefix:
    """`sp_qmmm_ref` under `stage: sp` is written as `sp_sp_qmmm_ref`, while
    `run.slurm sp sp_qmmm_ref` looks up `sp_qmmm_ref`. Jobs 31176144/31176145
    died in 8 s on that (D152), and nothing before the submission caught it.
    """

    def test_a_name_starting_with_its_stage_is_refused(self, mod, home, capsys):
        file = site(home, {"sp_qmmm_ref": {"stage": "sp", "inputs": {
            "pdb_file": "md_results/last.pdb"}}})
        with pytest.raises(SystemExit):
            mod.main(file, ["sp_qmmm_ref"])
        out = capsys.readouterr().out
        assert "sp_sp_qmmm_ref" in out and "qmmm_ref" in out

    def test_nothing_is_written_under_either_key(self, mod, home):
        file = site(home, {"sp_qmmm_ref": {"stage": "sp", "inputs": {
            "pdb_file": "md_results/last.pdb"}}})
        with pytest.raises(SystemExit):
            mod.main(file, ["sp_qmmm_ref"])
        assert written(home, "sp_sp_qmmm_ref") is None
        assert written(home, "sp_qmmm_ref") is None

    def test_a_name_equal_to_its_stage_is_refused(self, mod, home):
        file = site(home, {"sp": {"stage": "sp", "inputs": {
            "pdb_file": "md_results/last.pdb"}}})
        with pytest.raises(SystemExit):
            mod.main(file, ["sp"])

    def test_the_undoubled_name_is_written_as_the_key_mcube_is_given(self, mod, home):
        file = site(home, {"qmmm_ref": {"stage": "sp", "inputs": {
            "pdb_file": "md_results/last.pdb"}}})
        mod.main(file, ["qmmm_ref"])
        assert written(home, "sp_qmmm_ref")["pdb_file"] == "md_results/last.pdb"

    def test_a_stage_name_merely_contained_in_the_name_is_fine(self, mod, home):
        """`aspout` is not an `sp` section doubled; only a leading `<stage>_` is."""
        file = site(home, {"aspout": {"stage": "scan", "inputs": {
            "pdb_file": "md_results/last.pdb", "atoms": ["B:1:ZN"]}}})
        mod.main(file, ["aspout"])
        assert written(home, "scan_aspout") is not None
