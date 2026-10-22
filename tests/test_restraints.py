"""Detected coordination becomes a distance restraint: pytest tests/test_restraints.py

Needs the ASH environment, because mcube.md and mcube.preprocess import ASH.
"""

import pytest

pytest.importorskip("ash", reason="mcube.md imports ASH")

from mcube.md import md_restraints  # noqa: E402
from mcube.preprocess import COORDINATION_FORCE, pair_distances  # noqa: E402

# Zn, one donor at 2.10 A, and a pair that is only close across a box face.
ATOMS = [("B", "1", "ZN", 0.0, 0.0, 0.0),
         ("A", "88", "NE2", 0.0, 0.0, 2.1),
         ("E", "500", "O", 0.0, 0.0, 0.5),
         ("E", "501", "O", 0.0, 0.0, 9.5)]


def write_pdb(path, cell=None):
    with open(path, "w") as fh:
        if cell:
            fh.write("CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1\n"
                     % (cell, cell, cell))
        for n, (chain, resid, name, x, y, z) in enumerate(ATOMS, 1):
            fh.write("HETATM%5d %-4s %-3s %1s%4s    %8.3f%8.3f%8.3f\n"
                     % (n, name, "LIG", chain, resid, x, y, z))
    return str(path)


@pytest.fixture
def pdb(tmp_path):
    return write_pdb(tmp_path / "system.pdb")


@pytest.fixture
def periodic_pdb(tmp_path):
    return write_pdb(tmp_path / "periodic.pdb", cell=10.0)


class TestPairDistances:
    def test_measures_the_current_distance(self, pdb):
        assert pair_distances(pdb, [[0, 1]]) == [2.1]

    def test_minimum_image_when_the_file_carries_a_cell(self, periodic_pdb):
        # 9.0 A apart down the z axis, 1.0 A apart across the face
        assert pair_distances(periodic_pdb, [[2, 3]]) == [1.0]

    def test_no_cell_means_no_wrapping(self, pdb):
        assert pair_distances(pdb, [[2, 3]]) == [9.0]

    def test_trailing_entries_are_ignored(self, pdb):
        assert pair_distances(pdb, [[0, 1, 2.2, 20.0]]) == [2.1]


class TestMdRestraints:
    def test_rows_preprocess_wrote_are_used_as_given(self, pdb):
        """preprocess writes [i, j, d, k], so an edited k in the config is what runs."""
        restraints, detected = md_restraints(
            {}, {"results": {"coordination": [[0, 1, 2.1, 12.5]]}}, pdb)
        assert restraints == [[0, 1, 2.1, 12.5]]
        # Recorded apart from chosen restraints: it asserts nothing new (D68).
        assert detected == restraints

    def test_a_bare_pair_is_measured_at_this_run_s_geometry(self, pdb):
        restraints, _ = md_restraints(
            {}, {"results": {"coordination": [[0, 1]]}}, pdb)
        assert restraints == [[0, 1, 2.1, COORDINATION_FORCE]]

    def test_nothing_detected_gives_no_restraints(self, pdb):
        assert md_restraints({}, {"results": {}}, pdb) == (None, [])

    def test_an_explicit_restraint_on_the_same_pair_wins(self, pdb):
        restraints, detected = md_restraints(
            {"restraints": [[0, 1, 2.2, 20.0]]},
            {"results": {"coordination": [[0, 1, 2.1, 50.0]]}}, pdb)
        assert restraints == [[0, 1, 2.2, 20.0]]
        assert detected == []

    def test_chosen_and_detected_restraints_coexist(self, pdb):
        restraints, detected = md_restraints(
            {"restraints": [[0, 3, 3.0, 20.0]]},
            {"results": {"coordination": [[0, 1, 2.1, 50.0]]}}, pdb)
        assert restraints == [[0, 3, 3.0, 20.0], [0, 1, 2.1, 50.0]]
        assert detected == [[0, 1, 2.1, 50.0]]

    def test_selectors_resolve(self, pdb):
        restraints, _ = md_restraints(
            {}, {"results": {"coordination": [["B:1:ZN", "A:88:NE2"]]}}, pdb)
        assert restraints == [[0, 1, 2.1, COORDINATION_FORCE]]
