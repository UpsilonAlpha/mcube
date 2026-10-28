"""optimizer_constraints: the constraint set handed to ASH's Optimizer (qmmm `rigid_water`)."""

import pytest

pytest.importorskip("ash", reason="mcube.qmmm imports ASH")

from mcube.qmmm import optimizer_constraints  # noqa: E402

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
