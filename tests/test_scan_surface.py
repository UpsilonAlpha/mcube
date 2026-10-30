"""write_ash_surface: the profile ASH's surface tools read back.

A 2D scan's points carry `target1_A`/`target2_A`, a 1D scan's carry `target_A`,
and asking only for the 1D spelling killed every 2D scan at the very end - after
the last point had converged (D149).
"""

import pytest

pytest.importorskip("ash", reason="mcube.scan imports ASH")

from mcube.scan import write_ash_surface  # noqa: E402


def read(path):
    return [line.split() for line in open(path).read().splitlines()]


def test_1d_points_from_the_csv(tmp_path):
    points = [{"target_A": "1.90", "energy_Eh": "-447.3"},
              {"target_A": "2.10", "energy_Eh": "-447.4"}]
    out = write_ash_surface(points, False, str(tmp_path / "s.txt"))
    assert read(out) == [["1.90000000", "-447.30000000"],
                         ["2.10000000", "-447.40000000"]]


def test_2d_points_from_the_csv_use_target1(tmp_path):
    # The spelling job 31155441 died on
    points = [{"target1_A": "4.70", "target2_A": "3.60", "energy_Eh": "-447.31"},
              {"target1_A": "2.10", "target2_A": "1.90", "energy_Eh": "-447.28"}]
    out = write_ash_surface(points, True, str(tmp_path / "s.txt"))
    assert read(out) == [["4.70000000", "3.60000000", "-447.31000000"],
                         ["2.10000000", "1.90000000", "-447.28000000"]]


def test_2d_points_made_in_this_run(tmp_path):
    points = [{"target1": 4.7, "target2": 3.6, "energy": -447.31}]
    out = write_ash_surface(points, True, str(tmp_path / "s.txt"))
    assert read(out) == [["4.70000000", "3.60000000", "-447.31000000"]]


def test_a_point_missing_its_coordinate_still_raises(tmp_path):
    with pytest.raises(KeyError):
        write_ash_surface([{"energy_Eh": "-447.3"}], False, str(tmp_path / "s.txt"))
