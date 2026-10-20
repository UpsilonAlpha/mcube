"""Thread-count resolution: pytest tests/test_cores.py

Replaces the `resolve_cores` tests from the archived threads.py. Thread
*placement* is no longer mcube's business (the job script sets `OMP_PROC_BIND`
and friends), so all that is left to test is how a stage decides a count.

Needs the ASH environment, because mcube.preprocess imports ASH.
"""

import pytest

pytest.importorskip("ash", reason="mcube.preprocess imports ASH")

from mcube.preprocess import resolve_cores  # noqa: E402


@pytest.fixture(autouse=True)
def clean_env(monkeypatch):
    for var in ("OMP_NUM_THREADS", "SLURM_CPUS_PER_TASK"):
        monkeypatch.delenv(var, raising=False)


class TestResolveCores:
    def test_config_wins(self, monkeypatch):
        monkeypatch.setenv("OMP_NUM_THREADS", "4")
        assert resolve_cores(16) == 16

    def test_strings_from_yaml_are_accepted(self):
        assert resolve_cores("8") == 8

    def test_falls_back_to_omp_num_threads(self, monkeypatch):
        monkeypatch.setenv("OMP_NUM_THREADS", "4")
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
        assert resolve_cores(None) == 4

    def test_then_to_the_slurm_allocation(self, monkeypatch):
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
        assert resolve_cores(None) == 8

    def test_one_when_nothing_says_otherwise(self):
        """Better to be slow than to oversubscribe a shared node."""
        assert resolve_cores(None) == 1

    @pytest.mark.parametrize("value", ["", "   ", "0", "all", "-4", "4.5"])
    def test_unusable_env_values_are_ignored(self, monkeypatch, value):
        monkeypatch.setenv("OMP_NUM_THREADS", value)
        assert resolve_cores(None) == 1

    def test_label_appears_in_the_message(self, monkeypatch, capsys):
        """The log has to say which setting was missing: a stage resolves several."""
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
        resolve_cores(None, "sp_cores")
        assert "sp_cores" in capsys.readouterr().out
