"""Thread placement derived per stage: pytest tests/test_env.py

Needs nothing but PyYAML: cli.py imports no stage, which is the point -- a job
script calls `mcube env` before it has decided how to run anything.
"""

import textwrap

import pytest

from mcube.cli import USEFUL_CPUS, env_exports, thread_policy


@pytest.fixture
def config(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(textwrap.dedent("""
        preprocess: {inputs: {}}
        md: {inputs: {}}
        qmmm: {inputs: {cores: 4}}
        scan: {inputs: {}}
        ts: {inputs: {freq_workers: 1}}
        ts_parallel: {inputs: {freq_workers: 8}}
        qmmm_dft: {inputs: {qm_engine: dft}}
        sp: {inputs: {}}
    """))
    return str(path)


class TestThreadPolicy:
    @pytest.mark.parametrize("stage", ["qmmm", "scan", "ts"])
    def test_single_xtb_stages_bind(self, stage):
        """One xTB process at a time: binding was worth 2-9x on SMT nodes."""
        assert thread_policy(stage, {}, 16)["bind"] is True

    def test_half_the_allocation_is_one_thread_per_physical_core(self):
        assert thread_policy("qmmm", {}, 16)["OMP_NUM_THREADS"] == 8

    def test_one_thread_minimum(self):
        assert thread_policy("qmmm", {}, 1)["OMP_NUM_THREADS"] == 1

    def test_sp_leaves_the_count_to_pyscf(self):
        policy = thread_policy("sp", {}, 64)
        assert policy["OMP_NUM_THREADS"] == 1 and policy["bind"] is False

    def test_a_dft_engine_switches_any_stage_to_the_pyscf_policy(self):
        """qm_engine: dft means PySCF in-process, so the xTB rule is wrong."""
        policy = thread_policy("qmmm", {"qm_engine": "dft"}, 16)
        assert policy["OMP_NUM_THREADS"] == 1 and policy["bind"] is False

    def test_dft_engine_is_case_insensitive(self):
        assert thread_policy("scan", {"qm_engine": "DFT"}, 16)["bind"] is False

    def test_neb_never_binds(self):
        """Images run as parallel xTB processes; binding stacks them on one core."""
        assert thread_policy("neb", {}, 16)["bind"] is False

    def test_ts_stops_binding_once_there_are_workers(self):
        assert thread_policy("ts", {"freq_workers": 1}, 16)["bind"] is True
        assert thread_policy("ts", {"freq_workers": 8}, 16)["bind"] is False

    def test_md_does_not_bind_because_openmm_is_in_process(self):
        assert thread_policy("md", {}, 16)["bind"] is False

    def test_every_policy_explains_itself(self):
        for stage in ("preprocess", "md", "qmmm", "neb", "sp"):
            assert thread_policy(stage, {}, 8)["why"]


class TestEnvExports:
    def test_binding_stage_exports_all_three(self, config, monkeypatch):
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
        out = env_exports(config, "qmmm")
        assert "export OMP_NUM_THREADS=8" in out
        assert "export OMP_PROC_BIND=close" in out and "export OMP_PLACES=cores" in out

    def test_non_binding_stage_clears_stale_settings(self, config, monkeypatch):
        """The variables may be set in the user's environment already."""
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
        out = env_exports(config, "sp")
        assert "unset OMP_PROC_BIND OMP_PLACES" in out
        assert "export OMP_PROC_BIND" not in out

    def test_section_inputs_override_the_base(self, config, monkeypatch):
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
        assert "export OMP_PROC_BIND=close" in env_exports(config, "ts")
        assert "export OMP_PROC_BIND=close" not in env_exports(config, "ts", "ts_parallel")

    def test_a_dft_section_switches_policy(self, config, monkeypatch):
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
        assert "export OMP_NUM_THREADS=1" in env_exports(config, "qmmm", "qmmm_dft")

    def test_warns_when_the_allocation_is_wasted(self, config, monkeypatch):
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "64")
        out = env_exports(config, "sp")
        assert "WARNING" in out and str(USEFUL_CPUS["sp"]) in out

    def test_quiet_when_the_allocation_is_sensible(self, config, monkeypatch):
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
        assert "WARNING" not in env_exports(config, "sp")

    def test_output_is_shell_safe(self, config, monkeypatch):
        """It is eval'd by the job script, so the notes must be quoted."""
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "64")
        for line in env_exports(config, "sp").splitlines():
            assert line.startswith(("export ", "unset ", "echo '"))
        assert env_exports(config, "sp").endswith("\n")

    def test_missing_sections_are_not_an_error(self, tmp_path, monkeypatch):
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
        path = tmp_path / "bare.yaml"
        path.write_text("preprocess: {inputs: {}}\n")
        assert "export OMP_NUM_THREADS=" in env_exports(str(path), "qmmm")
