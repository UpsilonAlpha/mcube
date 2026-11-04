"""Monitored angles over a trajectory: pytest tests/test_angle_occupancy.py

`contact_occupancy` answers how close two atoms come, which is not enough to
say whether a nucleophile is attacking a carbon or lying across it. The CSR-1
work measured that angle with one-off scripts for three frames (D198) and so
had no record of it in any config, which is the same failure D168 records for
contacts.

The behaviours pinned here are the ones that carry history. A triple that does
not resolve is REPORTED, not skipped (D120). The selectors go through
`resolve_atoms`, so a config means one thing in every stage. The MIDDLE
selector is the vertex, because that is the order a chemist writes an angle in.
And with no `window` there is no criterion to score, so occupancy is None
rather than a number that looks like a result.

Needs MDAnalysis, which comes with the ASH environment.
"""

import sys
import types

import pytest

pytest.importorskip("MDAnalysis", reason="postprocess reads trajectories through MDAnalysis")
pytest.importorskip("openmm", reason="preprocess imports OpenMM at module scope")

from tests.test_qmmm_structure import ASH_NAMES


@pytest.fixture
def postprocess(monkeypatch):
    stub = types.ModuleType("ash")
    for name in ASH_NAMES:
        setattr(stub, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", stub)
    for name in [m for m in sys.modules if m == "mcube" or m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    from mcube import postprocess
    return postprocess


def atom(serial, name, resname, chain, seq, xyz, element):
    return ("HETATM%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
            % (serial, name, resname, chain, seq, xyz[0], xyz[1], xyz[2], element))


@pytest.fixture
def system(tmp_path):
    """A carbonyl carbon as the vertex, its oxygen fixed, and a mobile nucleophile.

    The legs are 1.5 A each so the angle is the only thing the frames change.
    """
    pdb = tmp_path / "top.pdb"
    pdb.write_text(
        "CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1           1\n"
        + atom(1, "C4x", "mer", "C", 1, (0.0, 0.0, 0.0), "C")
        + atom(2, "O2x", "mer", "C", 1, (1.5, 0.0, 0.0), "O")
        + atom(3, "O01", "OH", "C", 2, (0.0, 1.5, 0.0), "O")
        + "END\n")
    return pdb


def trajectory(tmp_path, system, degrees, radius=1.5):
    """Frames putting O01 at each named angle from the C4x-O2x leg."""
    import math
    lines = system.read_text().splitlines(keepends=True)
    out = tmp_path / "traj.pdb"
    with out.open("w") as fh:
        for n, deg in enumerate(degrees, start=1):
            fh.write("MODEL %8d\n" % n)
            for line in lines:
                if line.startswith("HETATM") and line[12:16].strip() == "O01":
                    x = radius * math.cos(math.radians(deg))
                    y = radius * math.sin(math.radians(deg))
                    line = line[:30] + "%8.3f%8.3f" % (x, y) + line[46:]
                if line.startswith(("ATOM", "HETATM")):
                    fh.write(line)
            fh.write("ENDMDL\n")
    return out


MONITOR = {"O01-C4x-O2x": ["C:2:O01", "C:1:C4x", "C:1:O2x"]}


class TestAngles:
    def test_the_middle_selector_is_the_vertex(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [60.0])
        rep = postprocess.angle_occupancy(str(system), str(traj), MONITOR)
        assert rep["frames"] == 1
        assert rep["triples"][0]["mean"] == pytest.approx(60.0, abs=0.2)

    def test_the_spread_is_the_range_visited(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [60.0, 95.0, 130.0])
        row = postprocess.angle_occupancy(str(system), str(traj), MONITOR)["triples"][0]
        assert row["min"] == pytest.approx(60.0, abs=0.2)
        assert row["max"] == pytest.approx(130.0, abs=0.2)
        assert row["mean"] == pytest.approx(95.0, abs=0.2)

    def test_a_window_scores_the_fraction_inside_it(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [60.0, 100.0, 110.0, 150.0])
        rep = postprocess.angle_occupancy(str(system), str(traj), MONITOR,
                                          window=[95, 120])
        assert rep["window"] == [95.0, 120.0]
        assert rep["triples"][0]["occupancy"] == 50.0

    def test_no_window_means_no_occupancy_rather_than_zero(self, postprocess,
                                                           tmp_path, system):
        """A None is a missing criterion; a 0.0 would read as a failed one."""
        traj = trajectory(tmp_path, system, [60.0, 100.0])
        rep = postprocess.angle_occupancy(str(system), str(traj), MONITOR)
        assert rep["window"] is None
        assert rep["triples"][0]["occupancy"] is None

    def test_stride_skips_frames(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [60.0, 150.0, 60.0, 150.0])
        rep = postprocess.angle_occupancy(str(system), str(traj), MONITOR, stride=2)
        assert rep["frames"] == 2
        assert rep["triples"][0]["max"] == pytest.approx(60.0, abs=0.2)


class TestMinimumImage:
    def test_a_leg_across_the_periodic_boundary_is_wrapped(self, postprocess,
                                                           tmp_path):
        """Without the minimum image this angle comes out of a 60 A leg."""
        pdb = tmp_path / "top.pdb"
        pdb.write_text(
            "CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1           1\n"
            + atom(1, "C4x", "mer", "C", 1, (0.0, 0.0, 0.0), "C")
            + atom(2, "O2x", "mer", "C", 1, (1.5, 0.0, 0.0), "O")
            + atom(3, "O01", "OH", "C", 2, (0.0, 58.5, 0.0), "O")
            + "END\n")
        traj = tmp_path / "traj.pdb"
        traj.write_text("MODEL        1\n"
                        + pdb.read_text().replace("CRYST1   60.000   60.000"
                                                  "   60.000  90.00  90.00"
                                                  "  90.00 P 1           1\n", "")
                        .replace("END\n", "")
                        + "ENDMDL\nEND\n")
        rep = postprocess.angle_occupancy(str(pdb), str(traj), MONITOR)
        # O01 is 1.5 A away across the boundary, at -y, so the angle is 90 deg.
        assert rep["triples"][0]["mean"] == pytest.approx(90.0, abs=0.2)


class TestMissingIsReported:
    def test_an_unresolvable_selector_is_reported_not_skipped(self, postprocess,
                                                              tmp_path, system):
        traj = trajectory(tmp_path, system, [60.0])
        monitor = dict(MONITOR)
        monitor["O01-C4x-NZ"] = ["C:2:O01", "C:1:C4x", "A:226:NZ"]
        rep = postprocess.angle_occupancy(str(system), str(traj), monitor)
        assert [r["triple"] for r in rep["triples"]] == ["O01-C4x-O2x"]
        assert len(rep["missing"]) == 1
        assert rep["missing"][0]["triple"] == "O01-C4x-NZ"
        assert "A:226:NZ" in rep["missing"][0]["why"]

    def test_a_pair_given_where_a_triple_is_wanted_is_reported(self, postprocess,
                                                               tmp_path, system):
        traj = trajectory(tmp_path, system, [60.0])
        monitor = {"two atoms": ["C:2:O01", "C:1:C4x"]}
        rep = postprocess.angle_occupancy(str(system), str(traj), monitor)
        assert rep["triples"] == []
        assert "three selectors" in rep["missing"][0]["why"]
