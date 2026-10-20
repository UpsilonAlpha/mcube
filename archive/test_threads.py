"""Tests for mcube.threads (thread-binding defaults): pytest tests/test_threads.py"""

import os

import pytest

from mcube import threads


@pytest.fixture(autouse=True)
def clean_env(monkeypatch):
    for k in ("OMP_PROC_BIND", "OMP_PLACES", "MCUBE_THREAD_BINDING", threads._BIND_XTB):
        monkeypatch.delenv(k, raising=False)


def test_python_process_is_never_bound():
    threads.configure_threads(verbose=False)
    assert "OMP_PROC_BIND" not in os.environ and "OMP_PLACES" not in os.environ
    assert os.environ[threads._BIND_XTB] == "1"


def test_respects_user_settings(monkeypatch):
    monkeypatch.setenv("OMP_PROC_BIND", "spread")
    threads.configure_threads(verbose=False)
    assert os.environ["OMP_PROC_BIND"] == "spread"
    assert threads.xtb_dir() is None


def test_opt_out(monkeypatch):
    monkeypatch.setenv("MCUBE_THREAD_BINDING", "0")
    threads.configure_threads(verbose=False)
    assert threads.xtb_dir() is None


def test_unbind_stops_xtb_binding():
    threads.configure_threads(verbose=False)
    threads.unbind_threads("parallel workers")
    assert threads.xtb_dir() is None


@pytest.mark.skipif(os.name == "nt", reason="wrapper is a POSIX shell script")
def test_xtb_wrapper(monkeypatch, tmp_path):
    fake = tmp_path / "bin" / "xtb"
    fake.parent.mkdir()
    fake.write_text("#!/bin/sh\n")
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", str(fake.parent) + os.pathsep + os.environ.get("PATH", ""))
    threads.configure_threads(verbose=False)
    d = threads.xtb_dir(str(tmp_path))
    text = (tmp_path / ".mcube_xtb" / "xtb").read_text()
    assert d.endswith(".mcube_xtb") and "OMP_PROC_BIND=close" in text and str(fake) in text
    assert "XTBHOME=" in text and os.access(tmp_path / ".mcube_xtb" / "xtb", os.X_OK)


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
        if str(path).startswith("/sys/devices/system/cpu/cpu"):
            cpu = int(str(path).split("/cpu")[2].split("/")[0])
            f = tmp_path / str(path).replace("/", "_")
            f.write_text(topo[cpu][0 if str(path).endswith("package_id") else 1])
            return real_open(f)
        return real_open(path, *a, **k)
    monkeypatch.setattr("builtins.open", fake_open)
    assert threads.physical_cores() == 2
