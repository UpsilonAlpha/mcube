"""Tests for mcube.threads (thread-binding defaults): pytest tests/test_threads.py"""

import pytest

from mcube import threads


@pytest.fixture(autouse=True)
def clean_env(monkeypatch):
    for k in ("OMP_PROC_BIND", "OMP_PLACES", "MCUBE_THREAD_BINDING", threads._SET_BY_MCUBE):
        monkeypatch.delenv(k, raising=False)


def test_sets_binding_when_unset(monkeypatch):
    n = threads.configure_threads(verbose=False)
    import os
    assert os.environ["OMP_PROC_BIND"] == "close" and os.environ["OMP_PLACES"] == "cores"
    assert n >= 1


def test_respects_user_settings(monkeypatch):
    import os
    monkeypatch.setenv("OMP_PROC_BIND", "spread")
    threads.configure_threads(verbose=False)
    assert os.environ["OMP_PROC_BIND"] == "spread" and "OMP_PLACES" not in os.environ
    threads.unbind_threads("test")  # not ours: left alone
    assert os.environ["OMP_PROC_BIND"] == "spread"


def test_opt_out(monkeypatch):
    import os
    monkeypatch.setenv("MCUBE_THREAD_BINDING", "0")
    threads.configure_threads(verbose=False)
    assert "OMP_PROC_BIND" not in os.environ


def test_unbind_undoes_mcube_defaults():
    import os
    threads.configure_threads(verbose=False)
    threads.unbind_threads("parallel workers")
    assert os.environ["OMP_PROC_BIND"] == "false" and "OMP_PLACES" not in os.environ


def test_resolve_cores(monkeypatch, capsys):
    monkeypatch.setattr(threads, "physical_cores", lambda: 4)
    assert threads.resolve_cores(None) == 4
    assert threads.resolve_cores(2) == 2
    assert threads.resolve_cores(8) == 8
    assert "exceeds the 4 physical cores" in capsys.readouterr().out


def test_physical_cores_counts_smt_siblings_once(monkeypatch, tmp_path):
    # CPUs 0,1 are siblings on core 0; CPU 2 is core 1
    topo = {0: ("0", "0"), 1: ("0", "0"), 2: ("0", "1")}
    monkeypatch.setattr(threads, "cpu_set", lambda: [0, 1, 2])
    real_open = open

    def fake_open(path, *a, **k):
        if path.startswith("/sys/devices/system/cpu/cpu"):
            cpu = int(path.split("/cpu")[2].split("/")[0])
            return real_open(_write(tmp_path, path, topo[cpu][0 if path.endswith("package_id") else 1]))
        return real_open(path, *a, **k)
    monkeypatch.setattr("builtins.open", fake_open)
    assert threads.physical_cores() == 2


def _write(tmp_path, path, value):
    f = tmp_path / path.replace("/", "_")
    f.write_text(value)
    return str(f)
