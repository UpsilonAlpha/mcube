"""How the sphere changed, and whether anything was holding it:
pytest tests/test_compare_coordination.py

`compare_coordination` has reported APPEARED/LEFT/MOVED since D120, where a Gln
amide oxygen went from 5.54 to 1.79 A unwatched. What it did not report was
whether a restraint was on the pair, and the advice printed beside it in `md`
assumed none were: "a contact that appeared was not constrained, and a
non-bonded metal has no Pauli repulsion". Under the standing rule that
everything inside the coordination radius is restrained at its crystal distance
(D186) that is wrong about most lines it prints - a held pair that moved was
moved by the restraint, which calls for loosening the restraint rather than for
distrusting the force field.
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


def row(donor, distance, metal="B:1:ZN", resname="ASP"):
    return {"metal": metal, "donor": donor, "resname": resname,
            "distance": distance, "metal_index": 0, "donor_index": 1}


BEFORE = [row("A:92:OD2", 2.80), row("A:88:NE2", 2.22)]


class TestTheReportItself:

    def test_an_unchanged_sphere_gives_no_lines(self, preprocess):
        assert preprocess.compare_coordination(BEFORE, BEFORE) == []

    def test_appeared_left_and_moved(self, preprocess):
        after = [row("A:92:OD2", 2.00), row("A:93:OE1", 1.79)]
        lines = preprocess.compare_coordination(BEFORE, after)
        assert any(l.startswith("  APPEARED") and "A:93:OE1" in l for l in lines)
        assert any(l.startswith("  LEFT") and "A:88:NE2" in l for l in lines)
        assert any(l.startswith("  MOVED") and "A:92:OD2" in l for l in lines)

    def test_a_move_under_the_threshold_is_not_a_change(self, preprocess):
        after = [row("A:92:OD2", 2.90), row("A:88:NE2", 2.22)]
        assert preprocess.compare_coordination(BEFORE, after) == []

    def test_no_restraint_list_tags_nothing(self, preprocess):
        """The default, and every caller that predates the argument."""
        after = [row("A:92:OD2", 2.00), row("A:88:NE2", 2.22)]
        lines = preprocess.compare_coordination(BEFORE, after)
        assert lines and not any("[restrained]" in l for l in lines)


class TestTheRestraintTag:

    def test_a_held_pair_is_tagged(self, preprocess):
        after = [row("A:92:OD2", 2.00), row("A:88:NE2", 2.22)]
        lines = preprocess.compare_coordination(
            BEFORE, after, [["B:1:ZN", "A:92:OD2", 2.80]])
        moved = [l for l in lines if "A:92:OD2" in l][0]
        assert "[restrained]" in moved

    def test_a_free_pair_beside_a_held_one_is_not(self, preprocess):
        after = [row("A:92:OD2", 2.00), row("A:93:OE1", 1.79)]
        lines = preprocess.compare_coordination(
            BEFORE, after, [["B:1:ZN", "A:92:OD2", 2.80]])
        appeared = [l for l in lines if "A:93:OE1" in l][0]
        assert "[restrained]" not in appeared
        assert "[restrained]" in [l for l in lines if "A:92:OD2" in l][0]

    def test_the_order_of_the_pair_does_not_matter(self, preprocess):
        """`to_selectors` writes metal-first, the shell row is metal-first too,
        but a hand-written `md.inputs.restraints` need not be."""
        after = [row("A:92:OD2", 2.00), row("A:88:NE2", 2.22)]
        lines = preprocess.compare_coordination(
            BEFORE, after, [["A:92:OD2", "B:1:ZN"]])
        assert "[restrained]" in [l for l in lines if "A:92:OD2" in l][0]

    def test_a_string_entry_is_accepted(self, preprocess):
        """Restraints may be written as one string per pair."""
        after = [row("A:92:OD2", 2.00), row("A:88:NE2", 2.22)]
        lines = preprocess.compare_coordination(
            BEFORE, after, ["B:1:ZN A:92:OD2 2.80"])
        assert "[restrained]" in [l for l in lines if "A:92:OD2" in l][0]

    def test_an_unresolved_entry_is_ignored_not_fatal(self, preprocess):
        """Indices rather than selectors, which `to_selectors` leaves alone when
        it cannot read the structure. Nothing is tagged; nothing raises."""
        after = [row("A:92:OD2", 2.00), row("A:88:NE2", 2.22)]
        lines = preprocess.compare_coordination(BEFORE, after, [[0, 1, 2.80]])
        assert lines and not any("[restrained]" in l for l in lines)

    def test_a_tag_does_not_disturb_the_numbers(self, preprocess):
        after = [row("A:92:OD2", 2.00), row("A:88:NE2", 2.22)]
        lines = preprocess.compare_coordination(
            BEFORE, after, [["B:1:ZN", "A:92:OD2", 2.80]])
        moved = [l for l in lines if "A:92:OD2" in l][0]
        assert "2.80 -> 2.00 A" in moved
