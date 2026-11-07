"""sections.py's declarative gates: pytest tests/test_sections_require.py

Every `docs/csr/*.sh` launcher was the same shell script with a different
hand-written preflight block in the middle, and they diverged: `aspqm2.sh`
*warned* when something was restrained where `aspdft.sh` *refused*. A gate
written per run is a gate that is wrong somewhere. These are the same checks as
config, so a new arm is a section in site.yaml and nothing else.

The one that earns its keep is `same_as`: a two-arm comparison with a second
variable changed by accident decides nothing, and nothing in the output of
either run would say so.

Needs only PyYAML.
"""

import importlib.util
import sys
import types
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")
pytest.importorskip("openmm", reason="sections.py imports mcube.preprocess")

ASH_NAMES = ("Fragment", "OpenMM_Modeller", "small_molecule_parameterizer",
             "merge_pdb_files")

SCRIPT = Path(__file__).resolve().parents[1] / "docs" / "csr" / "sections.py"

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
    (tmp_path / "md_results").mkdir()
    (tmp_path / "md_results" / "last.pdb").write_text(PDB)
    (tmp_path / "config.yaml").write_text("preprocess:\n  inputs: {}\n")
    monkeypatch.chdir(tmp_path)
    return tmp_path


def site(home, body):
    path = home / "site.yaml"
    path.write_text(yaml.safe_dump({"sections": body}, sort_keys=False))
    return str(path)


def written(home, section):
    data = yaml.safe_load(open(home / "config.yaml"))
    return (data.get(section) or {}).get("inputs")


BASE = {"pdb_file": "md_results/last.pdb", "charge": 2, "rigid_water": False}


class TestUnsetAndEquals:
    def test_a_section_passes_when_its_requirements_hold(self, mod, home):
        file = site(home, {"arm": {
            "stage": "qmmm",
            "require": {"unset": ["opt_constraints"], "equals": {"charge": 2}},
            "inputs": dict(BASE)}})
        mod.main(file, ["arm"])
        assert written(home, "qmmm_arm")["charge"] == 2

    def test_a_key_that_must_be_unset_stops_everything(self, mod, home):
        """D134: opt_constraints on a qmmm run makes convergence unreachable."""
        file = site(home, {"arm": {
            "stage": "qmmm",
            "require": {"unset": ["opt_constraints"]},
            "inputs": dict(BASE, opt_constraints=[["B:1:ZN", "A:88:NE2", 2.8]])}})
        with pytest.raises(SystemExit):
            mod.main(file, ["arm"])
        assert written(home, "qmmm_arm") is None

    def test_a_wrong_value_stops_everything(self, mod, home):
        file = site(home, {"arm": {
            "stage": "qmmm",
            "require": {"equals": {"charge": 2}},
            "inputs": dict(BASE, charge=0)}})
        with pytest.raises(SystemExit):
            mod.main(file, ["arm"])
        assert written(home, "qmmm_arm") is None

    def test_false_is_a_value_and_not_an_absence(self, mod, home):
        """`rigid_water: false` must satisfy equals, not read as unset."""
        file = site(home, {"arm": {
            "stage": "qmmm",
            "require": {"equals": {"rigid_water": False}},
            "inputs": dict(BASE)}})
        mod.main(file, ["arm"])
        assert written(home, "qmmm_arm")["rigid_water"] is False

    def test_an_unknown_require_key_is_refused_rather_than_ignored(self, mod, home):
        """A typo in a gate is a gate that does not run, which is the worst case."""
        file = site(home, {"arm": {
            "stage": "qmmm",
            "require": {"equal": {"charge": 2}},
            "inputs": dict(BASE)}})
        with pytest.raises(SystemExit):
            mod.main(file, ["arm"])


class TestSameAs:
    def ref_written(self, home):
        data = yaml.safe_load(open(home / "config.yaml"))
        data["qmmm_ref"] = {"inputs": dict(BASE, qm_engine="dft", opt_maxiter=60)}
        (home / "config.yaml").write_text(yaml.safe_dump(data, sort_keys=False))

    def test_identical_except_the_named_keys_passes(self, mod, home):
        self.ref_written(home)
        file = site(home, {"arm": {
            "stage": "qmmm",
            "require": {"same_as": {"section": "ref",
                                    "except": ["qm_engine", "opt_maxiter"]}},
            "inputs": dict(BASE, opt_maxiter=300)}})
        mod.main(file, ["arm"])
        assert written(home, "qmmm_arm")["opt_maxiter"] == 300

    def test_a_second_difference_stops_everything_and_is_named(self, mod, home):
        """The whole point: a comparison with two variables changed is not one."""
        self.ref_written(home)
        file = site(home, {"arm": {
            "stage": "qmmm",
            "require": {"same_as": {"section": "ref", "except": ["qm_engine"]}},
            "inputs": dict(BASE, charge=1, opt_maxiter=60)}})
        with pytest.raises(SystemExit):
            mod.main(file, ["arm"])
        assert written(home, "qmmm_arm") is None

    def test_the_reference_is_read_from_config_not_the_site_file(self, mod, home):
        """config.yaml is what the job reads, so it is what must match."""
        self.ref_written(home)
        file = site(home, {
            "ref": {"stage": "qmmm", "inputs": dict(BASE, qm_engine="dft",
                                                    opt_maxiter=60, charge=7)},
            "arm": {"stage": "qmmm",
                    "require": {"same_as": {"section": "ref",
                                            "except": ["qm_engine", "opt_maxiter"]}},
                    "inputs": dict(BASE, opt_maxiter=300)}})
        # site.yaml's `ref` has charge 7 and config.yaml's has 2; the arm has 2,
        # so this passes only if the written section is the reference.
        mod.main(file, ["arm"])
        assert written(home, "qmmm_arm")["charge"] == 2

    def test_a_missing_reference_is_an_error_not_a_pass(self, mod, home):
        file = site(home, {"arm": {
            "stage": "qmmm",
            "require": {"same_as": {"section": "nosuch"}},
            "inputs": dict(BASE)}})
        with pytest.raises(SystemExit):
            mod.main(file, ["arm"])

    def test_a_bare_string_names_the_reference(self, mod, home):
        self.ref_written(home)
        file = site(home, {"arm": {
            "stage": "qmmm", "require": {"same_as": "ref"},
            "inputs": dict(BASE, qm_engine="dft", opt_maxiter=60)}})
        mod.main(file, ["arm"])
        assert written(home, "qmmm_arm")["qm_engine"] == "dft"


class TestSubmit:
    def test_the_sbatch_line_comes_from_config(self, mod):
        cmd = mod.sbatch_command("arm", "qmmm_arm", "qmmm",
                                 {"cores": 8, "mem": "16G", "time": "24:00:00",
                                  "partition": "general"}, "/r/run.slurm")
        assert cmd == ["sbatch", "-J", "arm", "-p", "general",
                       "--cpus-per-task=8", "--mem=16G", "-t", "24:00:00",
                       "/r/run.slurm", "qmmm", "qmmm_arm"]

    def test_extra_flags_pass_through(self, mod):
        """A GPU submission needs `-q gpu`, which run.slurm hardcodes wrongly."""
        cmd = mod.sbatch_command("arm", "qmmm_arm", "qmmm",
                                 {"extra": ["-q", "gpu"]}, "/r/run.slurm")
        assert cmd[-5:-3] == ["-q", "gpu"]

    def test_a_dry_run_checks_without_writing(self, mod, home, capsys):
        file = site(home, {"arm": {
            "stage": "qmmm", "submit": {"cores": 8},
            "require": {"equals": {"charge": 2}},
            "inputs": dict(BASE)}})
        mod.main(file, ["arm"], dry_run=True)
        assert written(home, "qmmm_arm") is None
        assert "sbatch" in capsys.readouterr().out

    def test_a_dry_run_still_refuses_a_failing_gate(self, mod, home):
        file = site(home, {"arm": {
            "stage": "qmmm", "submit": {"cores": 8},
            "require": {"equals": {"charge": 2}},
            "inputs": dict(BASE, charge=0)}})
        with pytest.raises(SystemExit):
            mod.main(file, ["arm"], dry_run=True)
