"""`preprocess.inputs.remove_residues`: pytest tests/test_remove_residues.py

A substrate that coordinates a metal takes a position something else occupied,
and in this project's site that is one of two crystallographic waters (D120,
D196). There was no config route to that at all - the only way was to edit the
input PDB by hand - so what is pinned here is the behaviour that makes the
config route trustworthy:

  * the residue NAME is checked, so a stale selector fails instead of deleting
    whatever has since taken that number;
  * a selector that matches nothing raises, rather than succeeding quietly on
    a structure it was never meant for;
  * CONECT records naming a deleted atom go with it, because a dangling one
    surfaces several steps later as an error about a bond.
"""

import sys
import types

import pytest


@pytest.fixture(scope="module", autouse=True)
def _ash():
    """`preprocess` imports ASH at import time, and these two functions are
    pure text processing that has nothing to do with it - so it is stubbed
    rather than skipped, and the guards below are checked everywhere."""
    if "ash" not in sys.modules:
        stub = types.ModuleType("ash")
        for name in ("Fragment", "OpenMM_Modeller",
                     "small_molecule_parameterizer", "merge_pdb_files"):
            setattr(stub, name, lambda *a, **k: None)
        sys.modules["ash"] = stub


@pytest.fixture(scope="module")
def pre(_ash):
    from mcube import preprocess
    return preprocess


WATER = """\
ATOM      1  N   HIS A  88      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  NE2 HIS A  88      11.000  10.000  10.000  1.00  0.00           N
HETATM    3 ZN    ZN B   1      12.000  10.000  10.000  1.00  0.00          ZN
HETATM    4  O   HOH D   1      13.900  10.000  10.000  1.00  0.00           O
HETATM    5  H1  HOH D   1      14.400  10.700  10.000  1.00  0.00           H
HETATM    6  H2  HOH D   1      14.400   9.300  10.000  1.00  0.00           H
HETATM    7  O   HOH D   2      20.000  10.000  10.000  1.00  0.00           O
CONECT    4    5    6
CONECT    3    4
CONECT    1    2
END
"""


@pytest.fixture
def pdb(tmp_path):
    path = tmp_path / "site.pdb"
    path.write_text(WATER)
    return str(path)


def atoms(path):
    return [line for line in open(path)
            if line.startswith(("ATOM", "HETATM"))]


def conects(path):
    return [line.rstrip() for line in open(path) if line.startswith("CONECT")]


class TestParsing:
    def test_a_string_selector(self, pre):
        assert pre.parse_residue_selectors(["D:1:HOH"]) == [("D", 1, "HOH")]

    def test_a_mapping_is_accepted_like_every_other_selector_input(self, pre):
        assert pre.parse_residue_selectors({"D": {1: "HOH"}}) == [("D", 1, "HOH")]

    def test_the_name_is_lowercased_in_the_config_and_upper_here(self, pre):
        assert pre.parse_residue_selectors(["D:1:hoh"]) == [("D", 1, "HOH")]

    def test_a_selector_without_a_name_is_refused(self, pre):
        """Chain and number alone would delete whatever now sits there."""
        with pytest.raises(ValueError) as exc:
            pre.parse_residue_selectors(["D:1"])
        assert "RESNAME" in str(exc.value)

    def test_a_non_numeric_residue_is_refused(self, pre):
        with pytest.raises(ValueError) as exc:
            pre.parse_residue_selectors(["D:first:HOH"])
        assert "residue number" in str(exc.value)

    def test_nothing_is_nothing(self, pre):
        assert pre.parse_residue_selectors(None) == []


class TestRemoval:
    def test_the_whole_residue_goes_hydrogens_included(self, pre, pdb):
        out, rows = pre.remove_residues(pdb, ["D:1:HOH"])
        assert len(atoms(out)) == 4
        assert not [line for line in atoms(out) if line[21] == "D"
                    and line[22:26].strip() == "1"]
        assert rows == ["D:1:HOH (3 atoms)"]

    def test_the_other_water_stays(self, pre, pdb):
        out, _ = pre.remove_residues(pdb, ["D:1:HOH"])
        assert [line for line in atoms(out)
                if line[21] == "D" and line[22:26].strip() == "2"]

    def test_conect_records_naming_a_removed_atom_go_too(self, pre, pdb):
        """A CONECT pointing at a serial that no longer exists fails later, in
        the reader, as an error about the bond rather than the deletion."""
        out, _ = pre.remove_residues(pdb, ["D:1:HOH"])
        assert conects(out) == ["CONECT    1    2"]

    def test_the_input_file_is_not_touched(self, pre, pdb):
        pre.remove_residues(pdb, ["D:1:HOH"])
        assert len(atoms(pdb)) == 7

    def test_nothing_requested_returns_the_file_unchanged(self, pre, pdb):
        out, rows = pre.remove_residues(pdb, [])
        assert rows == [] and out.endswith("site.pdb")

    def test_several_at_once(self, pre, pdb):
        out, rows = pre.remove_residues(pdb, ["D:1:HOH", "D:2:HOH"])
        assert len(atoms(out)) == 3
        assert rows == ["D:1:HOH (3 atoms)", "D:2:HOH (1 atom)"]


class TestRefusals:
    def test_a_missing_residue_number_raises(self, pre, pdb):
        with pytest.raises(ValueError) as exc:
            pre.remove_residues(pdb, ["D:9:HOH"])
        assert "no residue 9" in str(exc.value)

    def test_the_wrong_name_at_the_right_number_raises(self, pre, pdb):
        """The case this guard exists for: solvating renumbers everything after
        the solute, so a selector written against an earlier build can point at
        a real residue that is not the intended one."""
        with pytest.raises(ValueError) as exc:
            pre.remove_residues(pdb, ["B:1:HOH"])
        assert "is ZN, not HOH" in str(exc.value)

    def test_a_refusal_removes_nothing_at_all(self, pre, pdb, tmp_path):
        """Partial success would leave a half-built system that still runs."""
        with pytest.raises(ValueError):
            pre.remove_residues(pdb, ["D:1:HOH", "D:9:HOH"])
        assert not list(tmp_path.glob("*_removed.pdb"))
