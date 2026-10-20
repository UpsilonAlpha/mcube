"""Restraint provenance: pytest tests/test_provenance.py

A harmonic restraint says "we do not know this distance and are imposing one".
A constraint on an auto-detected contact says "this distance is structural".
Chaining runs used to turn the first into the second silently, because the
detector only sees the number in the file. These are the tests for the stamp
that makes it visible and the check that refuses it.

Needs the ASH environment, because mcube.preprocess imports ASH.
"""

import pytest

pytest.importorskip("ash", reason="mcube.preprocess imports ASH")

from mcube.preprocess import (  # noqa: E402
    atom_selectors,
    check_inherited_restraints,
    read_restraint_stamps,
    stamp_restraints,
)

ATOMS = [("B", "1", "ZN"), ("C", "1", "O2x"), ("A", "88", "NE2")]


def write_pdb(path, atoms=ATOMS):
    with open(path, "w") as fh:
        fh.write("CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1\n")
        for n, (chain, resid, name) in enumerate(atoms, 1):
            fh.write("HETATM%5d %-4s %-3s %1s%4s    %8.3f%8.3f%8.3f\n"
                     % (n, name, "LIG", chain, resid, 0.0, 0.0, float(n)))
        fh.write("END\n")
    return str(path)


@pytest.fixture
def structure(tmp_path):
    return write_pdb(tmp_path / "out.pdb")


class TestStamping:
    def test_selectors_are_read_from_pdb_columns(self, structure):
        assert atom_selectors(structure) == ["B:1:ZN", "C:1:O2x", "A:88:NE2"]

    def test_round_trip(self, structure):
        """Recorded as names, not indices, so resolvation cannot shift them."""
        stamp_restraints([structure], [[0, 1, 2.2, 20.0]], structure)
        assert read_restraint_stamps(structure) == [("B:1:ZN", "C:1:O2x", 2.2, 20.0)]

    def test_coordinates_survive(self, structure):
        stamp_restraints([structure], [[0, 1, 2.2, 20.0]], structure)
        atoms = [l for l in open(structure) if l.startswith(("ATOM", "HETATM"))]
        assert len(atoms) == len(ATOMS)

    def test_stamping_twice_does_not_duplicate(self, structure):
        for _ in range(2):
            stamp_restraints([structure], [[0, 1, 2.2, 20.0]], structure)
        assert len(read_restraint_stamps(structure)) == 1

    def test_no_restraints_writes_nothing(self, structure):
        before = open(structure).read()
        stamp_restraints([structure], [], structure)
        assert open(structure).read() == before

    def test_missing_files_are_skipped(self, structure, tmp_path):
        stamp_restraints([str(tmp_path / "nope.pdb"), None], [[0, 1, 2.2, 20.0]], structure)

    def test_unstamped_structure_reads_empty(self, structure):
        assert read_restraint_stamps(structure) == []
        assert read_restraint_stamps(None) == []


class TestRefusingToLaunder:
    @pytest.fixture
    def inherited(self, structure):
        stamp_restraints([structure], [[0, 1, 2.2, 20.0]], structure)
        return read_restraint_stamps(structure)

    def test_refuses_and_names_the_pair(self, structure, inherited):
        with pytest.raises(ValueError) as excinfo:
            check_inherited_restraints([[0, 1], [0, 2]], inherited, structure, {})
        message = str(excinfo.value)
        assert "B:1:ZN" in message and "C:1:O2x" in message and "2.20" in message
        assert "acknowledge_restrained" in message
        assert "constraint_exclude_residues" in message

    def test_pair_order_does_not_matter(self, structure, inherited):
        with pytest.raises(ValueError):
            check_inherited_restraints([[1, 0]], inherited, structure, {})

    def test_acknowledging_it_allows_the_constraint(self, structure, inherited):
        check_inherited_restraints([[0, 1]], inherited, structure,
                                   {"acknowledge_restrained": ["B:1:ZN C:1:O2x"]})

    def test_acknowledging_accepts_a_pair_as_a_list(self, structure, inherited):
        check_inherited_restraints([[0, 1]], inherited, structure,
                                   {"acknowledge_restrained": [["B:1:ZN", "C:1:O2x"]]})

    def test_acknowledging_one_pair_does_not_allow_another(self, structure):
        stamp_restraints([structure], [[0, 1, 2.2, 20.0], [0, 2, 2.0, 20.0]], structure)
        with pytest.raises(ValueError):
            check_inherited_restraints([[0, 1], [0, 2]], read_restraint_stamps(structure),
                                       structure, {"acknowledge_restrained": ["B:1:ZN C:1:O2x"]})

    def test_contacts_that_were_never_restrained_pass(self, structure, inherited):
        """Real coordination -- histidines, waters -- must not be obstructed."""
        check_inherited_restraints([[0, 2]], inherited, structure, {})

    def test_no_stamps_is_a_no_op(self, structure):
        check_inherited_restraints([[0, 1]], [], structure, {})

    def test_no_constraints_is_a_no_op(self, structure, inherited):
        check_inherited_restraints(None, inherited, structure, {})
        check_inherited_restraints([], inherited, structure, {})
