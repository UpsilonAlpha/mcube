"""Thread placement derived per stage: pytest tests/test_env.py

Needs only PyYAML: cli.py imports no stage, which is the point -- a job script
calls `mcube env` before it has decided how to run one.
"""

import textwrap

import pytest

from mcube.cli import env_exports, thread_policy


@pytest.fixture
def config(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(textwrap.dedent("""
        qmmm: {inputs: {cores: 4}}
        ts: {inputs: {freq_workers: 1}}
        ts_parallel: {inputs: {freq_workers: 8}}
        qmmm_dft: {inputs: {qm_engine: dft}}
        sp: {inputs: {}}
    """))
    return str(path)


@pytest.mark.parametrize("stage, inputs, threads, bind", [
    ("qmmm", {}, 8, True),                       # one xTB process: binding was worth 2-9x
    ("scan", {}, 8, True),
    ("ts", {"freq_workers": 1}, 8, True),
    ("ts", {"freq_workers": 8}, 8, False),       # parallel xTB: would stack on one core
    ("neb", {}, 8, False),                       # images always run in parallel
    ("md", {}, 8, False),                        # OpenMM is in-process
    ("sp", {}, 1, False),                        # PySCF sets its own count
    ("qmmm", {"qm_engine": "DFT"}, 1, False),    # so the xTB rule is wrong here
    ("preprocess", {}, 1, False),
])
def test_policy(stage, inputs, threads, bind):
    policy = thread_policy(stage, inputs, 16)
    assert (policy["OMP_NUM_THREADS"], policy["bind"]) == (threads, bind)
    assert policy["why"]


def test_one_thread_minimum():
    assert thread_policy("qmmm", {}, 1)["OMP_NUM_THREADS"] == 1


def test_binding_stage_exports_all_three(config, monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
    out = env_exports(config, "qmmm")
    assert "export OMP_NUM_THREADS=8" in out
    assert "export OMP_PROC_BIND=close" in out and "export OMP_PLACES=cores" in out


def test_non_binding_stage_clears_stale_settings(config, monkeypatch):
    """They may already be set in the environment the job inherited."""
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
    out = env_exports(config, "sp")
    assert "unset OMP_PROC_BIND OMP_PLACES" in out and "export OMP_PROC_BIND" not in out


@pytest.mark.parametrize("section, binds", [(None, True), ("ts_parallel", False)])
def test_section_inputs_override_the_base(config, monkeypatch, section, binds):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
    assert ("export OMP_PROC_BIND=close" in env_exports(config, "ts", section)) is binds


def test_warns_only_when_the_allocation_is_wasted(config, monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "64")
    assert "WARNING" in env_exports(config, "sp")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
    assert "WARNING" not in env_exports(config, "sp")


def test_output_is_shell_safe(config, monkeypatch):
    """It is eval'd by the job script, so every line must be quoted."""
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "64")
    out = env_exports(config, "sp")
    assert all(l.startswith(("export ", "unset ", "echo '")) for l in out.splitlines())
    assert out.endswith("\n")


def test_missing_sections_are_not_an_error(tmp_path, monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    path = tmp_path / "bare.yaml"
    path.write_text("preprocess: {inputs: {}}\n")
    assert "export OMP_NUM_THREADS=" in env_exports(str(path), "qmmm")
