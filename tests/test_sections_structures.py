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
