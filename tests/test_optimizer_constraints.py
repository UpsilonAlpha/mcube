"""optimizer_constraints and rigid_water_pairs: the constraints handed to ASH's Optimizer."""

import pytest

pytest.importorskip("ash", reason="mcube.qmmm imports ASH")

from mcube import qmmm as qmmm_mod  # noqa: E402
from mcube.qmmm import optimizer_constraints, rigid_water_pairs  # noqa: E402

COORDS = [[0.0, 0.0, 0.0], [0.9572, 0.0, 0.0], [-0.24, 0.927, 0.0], [3.0, 0.0, 0.0]]


def test_no_waters_and_no_user_pairs_means_no_constraints():
    # rigid_water: false with no opt_constraints: geomeTRIC must see none (cart refuses any)
    assert optimizer_constraints([], [], COORDS) == (None, False)


def test_rigid_waters_alone_are_frozen():
    constraints, constrainvalue = optimizer_constraints([[0, 1], [0, 2]], [], COORDS)
    assert constraints == {"bond": [[0, 1], [0, 2]]}
    assert constrainvalue is False


def test_user_pairs_without_waters_get_their_distance():
    constraints, constrainvalue = optimizer_constraints([], [[0, 3]], COORDS)
    assert constrainvalue is True
    assert constraints == {"bond": [[0, 3, 3.0]]}


def test_waters_are_copied_not_shared():
    water = [[0, 1]]
    constraints, _ = optimizer_constraints(water, [], COORDS)
    constraints["bond"][0][0] = 99        # ASH rewrites indices in place
    assert water == [[0, 1]]

# ---- rigid_water_pairs: off by default, because a `$set` water distance gates
# convergence and rigid TIP3P cannot satisfy it (D137). ----

def test_rigid_water_is_off_by_default():
    # mm is never touched when the option is off, so None is a valid stand-in
    assert rigid_water_pairs(None, [0, 1, 2], {}) == []


def test_rigid_water_false_holds_nothing():
    assert rigid_water_pairs(None, [0, 1, 2], {"rigid_water": False}) == []


def test_rigid_water_true_delegates_to_water_constraints(monkeypatch):
    seen = {}

    def fake(mm, actatoms, qm_atoms=()):
        seen.update(mm=mm, actatoms=list(actatoms), qm_atoms=tuple(qm_atoms))
        return [[0, 1], [0, 2]]

    monkeypatch.setattr(qmmm_mod, "water_constraints", fake)
    pairs = rigid_water_pairs("MM", [0, 1, 2], {"rigid_water": True}, (5,))
    assert pairs == [[0, 1], [0, 2]]
    assert seen == {"mm": "MM", "actatoms": [0, 1, 2], "qm_atoms": (5,)}


def test_default_scan_point_sets_only_its_own_bond():
    # What a scan now hands geomeTRIC: one `$set` bond, not that plus every water
    constraints, constrainvalue = optimizer_constraints(
        rigid_water_pairs(None, [0, 1, 2, 3], {}), [[0, 3, 2.5]], COORDS)
    assert constrainvalue is True
    assert constraints == {"bond": [[0, 3, 2.5]]}
