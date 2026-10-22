"""The TS imaginary mode must belong to the reaction: pytest tests/test_ts_mode.py

Runs without ASH by importing mcube against a stub, as tests/test_umbrella.py
does. Needs numpy only; nothing here reads a structure.
"""

import sys
import types

import numpy as np
import pytest

ASH_NAMES = ("Fragment", "OpenMMTheory", "QMMMTheory", "xTBTheory", "NumFreq", "Optimizer",
             "SellaIRC", "Singlepoint", "OpenMM_Opt", "OpenMM_MD", "OpenMM_Modeller",
             "OpenMM_box_equilibration", "Gentle_warm_up_MD", "NEB", "NEBTS",
             "small_molecule_parameterizer", "merge_pdb_files")

# Hessian over five atoms; indices 10 and 11 are the bond the step breaks.
HESSATOMS = [10, 11, 20, 21, 22]
MONITOR = [[10, 11]]


@pytest.fixture
def ts(monkeypatch):
    stub = types.ModuleType("ash")
    for name in ASH_NAMES:
        setattr(stub, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", stub)
    for name in [m for m in sys.modules if m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    from mcube import ts
    return ts


def test_mode_on_the_reacting_atoms_passes(ts):
    amp = np.array([1.0, 0.8, 0.1, 0.05, 0.05])
    ok, share, wanted = ts.check_ts_mode(amp, HESSATOMS, MONITOR)
    assert ok is True and share == 1.0 and wanted == [10, 11]


def test_mode_on_a_spectator_is_flagged(ts):
    """A water libration: one imaginary frequency, wrong mode."""
    amp = np.array([0.05, 0.03, 1.0, 0.2, 0.1])
    ok, share, _ = ts.check_ts_mode(amp, HESSATOMS, MONITOR)
    assert ok is False and share == pytest.approx(0.05)


def test_borderline_respects_the_floor(ts):
    amp = np.array([0.31, 0.1, 1.0, 0.0, 0.0])
    assert ts.check_ts_mode(amp, HESSATOMS, MONITOR)[0] is True
    assert ts.check_ts_mode(amp, HESSATOMS, MONITOR, floor=0.5)[0] is False


def test_no_monitor_pairs_is_unknown_not_a_pass(ts):
    amp = np.array([1.0, 0.8, 0.1, 0.05, 0.05])
    ok, share, wanted = ts.check_ts_mode(amp, HESSATOMS, None)
    assert ok is None and share is None and wanted == []


def test_monitor_atoms_outside_the_hessian_are_unknown(ts):
    """Monitoring an atom the Hessian never covered cannot vouch for the mode."""
    amp = np.array([1.0, 0.8, 0.1, 0.05, 0.05])
    assert ts.check_ts_mode(amp, HESSATOMS, [[90, 91]])[0] is None


def test_selector_pairs_with_trailing_numbers_are_handled(ts):
    """`monitor` rows may carry a distance; only the first two entries are atoms."""
    amp = np.array([1.0, 0.8, 0.1, 0.05, 0.05])
    ok, _, wanted = ts.check_ts_mode(amp, HESSATOMS, [[10, 11, 2.2]])
    assert ok is True and wanted == [10, 11]
