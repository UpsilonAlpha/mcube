"""What detection restrains: pytest tests/test_coordination_detection.py

`coordination_from_rows` has said since it was written that the QM region wants
whole residues and restraints want individual atoms. `find_coordinators` - the
one caller that builds restraints - took the per-residue default anyway, so
only the closest donor of each residue was ever held and a group with two
donors inside the cutoff had one pinned and the other free.

That is the geometry D120 and D130 both record failing: Asp92's OD2 held at
2.80 while OD1 pivoted onto the metal at 1.99. It bit again on the product,
whose ligand has O5x at 2.92 and O4x at 2.94 - per residue, O4x is left free
0.02 A further out than the oxygen beside it (D175).

Nothing tested this, which is why it survived two entries about it. These do.

Needs OpenMM, which preprocess imports at module scope.
"""

import sys
import types

import pytest

pytest.importorskip("openmm", reason="preprocess imports OpenMM at module scope")

from tests.test_qmmm_structure import ASH_NAMES


@pytest.fixture
def preprocess(monkeypatch):
    stub = types.ModuleType("ash")
    for name in ASH_NAMES:
        setattr(stub, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", stub)
    for name in [m for m in sys.modules if m == "mcube" or m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    from mcube import preprocess
    return preprocess


def atom(serial, name, resname, chain, seq, xyz, element=None):
    return ("HETATM%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
            % (serial, name, resname, chain, seq, xyz[0], xyz[1], xyz[2],
               element or name[0]))


@pytest.fixture
def carboxylate(tmp_path):
    """A zinc, a histidine, and a carboxylate with BOTH oxygens inside 3.0 A.

    The product's own geometry in miniature: O5x at 2.92 and O4x at 2.94 are
    0.02 A apart, and per-residue detection keeps exactly one of them.
    """
    pdb = tmp_path / "s.pdb"
    pdb.write_text(
        "CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1           1\n"
        + atom(1, "ZN", "ZN", "B", 1, (0.0, 0.0, 0.0), "ZN")
        + atom(2, "NE2", "HIS", "A", 88, (2.09, 0.0, 0.0), "N")
        + atom(3, "OD2", "ASP", "A", 92, (0.0, 2.80, 0.0), "O")
        + atom(4, "OD1", "ASP", "A", 92, (0.0, 0.0, 2.94), "O")
        + "END\n")
    return str(pdb)


INPUTS = {"metals": ["ZN"], "cutoff": 3.0, "donor_elements": ["N", "O", "S"]}


def pairs(coordination):
    return {(row[0], row[1]): row[2] for row in coordination}


class TestRestraintsArePerAtom:
    def test_both_oxygens_of_a_carboxylate_are_restrained(self, preprocess,
                                                          carboxylate):
        """The fix. Per residue this list has three rows and OD1 is free."""
        _, coordination = preprocess.find_coordinators(carboxylate, INPUTS)
        got = pairs(coordination)
        assert ("B:1:ZN", "A:92:OD2") in got
        assert ("B:1:ZN", "A:92:OD1") in got
        assert len(coordination) == 3

    def test_each_is_held_at_its_own_distance(self, preprocess, carboxylate):
        """Not at the closer oxygen's distance: the restraint preserves the
        geometry the file has, which is the whole premise."""
        got = pairs(preprocess.find_coordinators(carboxylate, INPUTS)[1])
        assert got[("B:1:ZN", "A:92:OD2")] == pytest.approx(2.80, abs=0.01)
        assert got[("B:1:ZN", "A:92:OD1")] == pytest.approx(2.94, abs=0.01)

    def test_a_donor_outside_the_cutoff_is_not_restrained(self, preprocess,
                                                          tmp_path):
        pdb = tmp_path / "far.pdb"
        pdb.write_text(
            "CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1           1\n"
            + atom(1, "ZN", "ZN", "B", 1, (0.0, 0.0, 0.0), "ZN")
            + atom(2, "NE2", "HIS", "A", 88, (2.09, 0.0, 0.0), "N")
            + atom(3, "OD1", "ASP", "A", 92, (0.0, 3.81, 0.0), "O")
            + "END\n")
        got = pairs(preprocess.find_coordinators(str(pdb), INPUTS)[1])
        assert ("B:1:ZN", "A:92:OD1") not in got

    def test_the_cutoff_decides_and_three_includes_what_two_point_five_omits(
            self, preprocess, carboxylate):
        """The 2.5 A cutoff left everything from 2.5 to 3.0 A unheld rather
        than protected, which is how Asp92 OD2 at 2.80 came to collapse to
        2.03 under MM (D173, D174)."""
        tight = pairs(preprocess.find_coordinators(
            carboxylate, dict(INPUTS, cutoff=2.5))[1])
        assert ("B:1:ZN", "A:92:OD2") not in tight
        wide = pairs(preprocess.find_coordinators(carboxylate, INPUTS)[1])
        assert ("B:1:ZN", "A:92:OD2") in wide


class TestTheQmRegionIsStillPerResidue:
    def test_a_residue_with_two_donors_appears_once(self, preprocess, carboxylate):
        """Naming the same residue twice sets the same dict entry twice, so
        per-atom restraints cannot double a QM region."""
        qm, _ = preprocess.find_coordinators(carboxylate, INPUTS)
        assert qm["A"] == {88: "HIS", 92: "ASP"}

    def test_the_metal_is_in_the_qm_region(self, preprocess, carboxylate):
        qm, _ = preprocess.find_coordinators(carboxylate, INPUTS)
        assert qm["B"] == {1: "ZN"}


class TestExclusions:
    def test_an_excluded_residue_gets_no_restraint(self, preprocess, carboxylate):
        """`constraint_exclude_residues` is right for a substrate that must be
        free to react, and it is what left the product's zinc contacts unheld
        during an equilibration (D173)."""
        inputs = dict(INPUTS, constraint_exclude_residues={"A": [92]})
        got = pairs(preprocess.find_coordinators(carboxylate, inputs)[1])
        assert ("B:1:ZN", "A:92:OD2") not in got
        assert ("B:1:ZN", "A:92:OD1") not in got
        assert ("B:1:ZN", "A:88:NE2") in got

    def test_an_excluded_residue_is_still_restrained_when_only_qm_excluded(
            self, preprocess, carboxylate):
        """`exclude_residues` keeps a residue out of the QM region and keeps
        its restraint - the two exclusions are not the same switch."""
        inputs = dict(INPUTS, exclude_residues={"A": [92]})
        qm, coordination = preprocess.find_coordinators(carboxylate, inputs)
        assert 92 not in qm.get("A", {})
        assert ("B:1:ZN", "A:92:OD2") in pairs(coordination)
