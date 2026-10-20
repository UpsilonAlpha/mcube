"""preprocess refuses to turn a restraint into a constraint: pytest tests/test_provenance.py

Needs the ASH environment, because mcube.preprocess imports ASH.
"""

import pytest

pytest.importorskip("ash", reason="mcube.preprocess imports ASH")

from mcube.preprocess import atom_names, check_inherited_restraints, to_selectors  # noqa: E402

ATOMS = [("B", "1", "ZN"), ("C", "1", "O2x"), ("A", "88", "NE2")]
INDEXED = [[0, 1, 2.2, 20.0]]                           # how older configs stored it
SELECTED = [["B:1:ZN", "C:1:O2x", 2.2, 20.0]]           # how md records it now
DETECTED = [[0, 1], [0, 2]]                             # contacts preprocess just found


@pytest.fixture
def pdb(tmp_path):
    path = tmp_path / "system.pdb"
    with open(path, "w") as fh:
        for n, (chain, resid, name) in enumerate(ATOMS, 1):
            fh.write("HETATM%5d %-4s %-3s %1s%4s    %8.3f%8.3f%8.3f\n"
                     % (n, name, "LIG", chain, resid, 0.0, 0.0, float(n)))
    return str(path)


class TestSelectors:
    def test_read_from_pdb_columns(self, pdb):
        assert atom_names(pdb) == ["B:1:ZN", "C:1:O2x", "A:88:NE2"]

    @pytest.mark.parametrize("entries, expected", [
        ([[0, 1, 2.2, 20.0]], [["B:1:ZN", "C:1:O2x", 2.2, 20.0]]),   # indices
        ([["B:1:ZN", "C:1:O2x"]], [["B:1:ZN", "C:1:O2x"]]),          # already selectors
        ([[0, "C:1:O2x", 2.2]], [["B:1:ZN", "C:1:O2x", 2.2]]),       # mixed
        (["B:1:ZN C:1:O2x"], [["B:1:ZN", "C:1:O2x"]]),               # one string
        (None, []),
        ([], []),
    ])
    def test_normalisation(self, pdb, entries, expected):
        assert to_selectors(entries, pdb) == expected

    def test_trailing_numbers_are_left_alone(self, pdb):
        assert to_selectors([[0, 1, 2.2, 20.0]], pdb)[0][2:] == [2.2, 20.0]


class TestRefusingToLaunder:
    def test_refuses_and_says_which_pair_and_how_to_resolve_it(self, pdb):
        with pytest.raises(ValueError) as excinfo:
            check_inherited_restraints(DETECTED, SELECTED, pdb, {})
        message = str(excinfo.value)
        assert "B:1:ZN" in message and "C:1:O2x" in message and "2.2" in message
        assert "constraint_exclude_residues" in message and "acknowledge_restrained" in message

    @pytest.mark.parametrize("restraints", [INDEXED, SELECTED])
    def test_either_stored_form_is_caught(self, pdb, restraints):
        """Configs written before restraints were recorded as selectors still work."""
        with pytest.raises(ValueError):
            check_inherited_restraints(DETECTED, restraints, pdb, {})

    def test_pair_order_does_not_matter(self, pdb):
        with pytest.raises(ValueError):
            check_inherited_restraints([[1, 0]], SELECTED, pdb, {})

    @pytest.mark.parametrize("acknowledged", [[[0, 1]], [["B:1:ZN", "C:1:O2x"]],
                                              ["B:1:ZN C:1:O2x"]])
    def test_acknowledging_the_pair_allows_it(self, pdb, acknowledged):
        check_inherited_restraints(DETECTED, SELECTED, pdb,
                                   {"acknowledge_restrained": acknowledged})

    @pytest.mark.parametrize("constraints, restraints", [
        ([[0, 2]], SELECTED),   # a contact that was never restrained
        (DETECTED, []),         # md imposed nothing
        (DETECTED, None),       # md has not run
        (None, SELECTED),       # no metals, so nothing detected
    ])
    def test_leaves_everything_else_alone(self, pdb, constraints, restraints):
        check_inherited_restraints(constraints, restraints, pdb, {})
