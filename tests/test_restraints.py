"""Detected coordination becomes a distance restraint: pytest tests/test_restraints.py

Runs without ASH. Nothing tested here touches it, but `md.py` and
`preprocess.py` both `from ash import *` at module scope, so mcube is imported
against a stub module instead - the same approach as tests/test_umbrella.py.
Needs OpenMM (the `mcube-test` env), because preprocess reads PDBs through it.
"""

import sys
import types

import pytest

pytest.importorskip("openmm", reason="preprocess reads structures through OpenMM")

# Zn, one donor at 2.10 A, and a pair that is only close across a box face.
ATOMS = [("B", "1", "ZN", 0.0, 0.0, 0.0),
         ("A", "88", "NE2", 0.0, 0.0, 2.1),
         ("E", "500", "O", 0.0, 0.0, 0.5),
         ("E", "501", "O", 0.0, 0.0, 9.5)]

# Every ASH name md.py and preprocess.py pull in at import; none is called here.
ASH_NAMES = ("Fragment", "OpenMMTheory", "OpenMM_Opt", "OpenMM_MD",
             "OpenMM_box_equilibration", "Gentle_warm_up_MD", "OpenMM_Modeller",
             "small_molecule_parameterizer", "merge_pdb_files", "Singlepoint",
             "QMMMTheory", "xTBTheory")


@pytest.fixture
def mc(monkeypatch):
    """`mcube.md` and `mcube.preprocess`, imported against a stub `ash`."""
    stub = types.ModuleType("ash")
    for name in ASH_NAMES:
        setattr(stub, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", stub)
    for name in [m for m in sys.modules if m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    from mcube import md, preprocess
    return types.SimpleNamespace(md=md, preprocess=preprocess)


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
    def test_measures_the_current_distance(self, mc, pdb):
        assert mc.preprocess.pair_distances(pdb, [[0, 1]]) == [2.1]

    def test_minimum_image_when_the_file_carries_a_cell(self, mc, periodic_pdb):
        # 9.0 A apart down the z axis, 1.0 A apart across the face
        assert mc.preprocess.pair_distances(periodic_pdb, [[2, 3]]) == [1.0]

    def test_no_cell_means_no_wrapping(self, mc, pdb):
        assert mc.preprocess.pair_distances(pdb, [[2, 3]]) == [9.0]

    def test_trailing_entries_are_ignored(self, mc, pdb):
        assert mc.preprocess.pair_distances(pdb, [[0, 1, 2.2, 20.0]]) == [2.1]


class TestMdRestraints:
    def test_rows_preprocess_wrote_are_used_as_given(self, mc, pdb):
        """preprocess writes [i, j, d, k], so an edited k in the config is what runs."""
        restraints, detected, _ = mc.md.md_restraints(
            {}, {"results": {"coordination": [[0, 1, 2.1, 12.5]]}}, pdb)
        assert restraints == [[0, 1, 2.1, 12.5]]
        # Recorded apart from chosen restraints: it asserts nothing new.
        assert detected == restraints

    def test_a_bare_pair_is_measured_at_this_runs_geometry(self, mc, pdb):
        restraints, _, _ = mc.md.md_restraints(
            {}, {"results": {"coordination": [[0, 1]]}}, pdb)
        assert restraints == [[0, 1, 2.1, mc.preprocess.COORDINATION_FORCE]]

    def test_nothing_detected_gives_no_restraints(self, mc, pdb):
        assert mc.md.md_restraints({}, {"results": {}}, pdb)[:2] == (None, [])

    def test_an_explicit_restraint_on_the_same_pair_wins(self, mc, pdb):
        restraints, detected, _ = mc.md.md_restraints(
            {"restraints": [[0, 1, 2.2, 20.0]]},
            {"results": {"coordination": [[0, 1, 2.1, 50.0]]}}, pdb)
        assert restraints == [[0, 1, 2.2, 20.0]]
        assert detected == []

    def test_chosen_and_detected_restraints_coexist(self, mc, pdb):
        restraints, detected, _ = mc.md.md_restraints(
            {"restraints": [[0, 3, 3.0, 20.0]]},
            {"results": {"coordination": [[0, 1, 2.1, 50.0]]}}, pdb)
        assert restraints == [[0, 3, 3.0, 20.0], [0, 1, 2.1, 50.0]]
        assert detected == [[0, 1, 2.1, 50.0]]

    def test_an_asked_row_inherits_the_default_force_constant(self, mc, pdb):
        """One rule for every pair: a named restraint is filled like a detected one."""
        restraints, detected, _ = mc.md.md_restraints(
            {"restraints": [[0, 1, 2.2]]}, {"results": {}}, pdb)
        assert restraints == [[0, 1, 2.2, mc.preprocess.COORDINATION_FORCE]]
        assert detected == []

    def test_an_asked_row_with_no_distance_is_held_where_it_is(self, mc, pdb):
        restraints, _, _ = mc.md.md_restraints(
            {"restraints": [[0, 1]]}, {"results": {}}, pdb)
        assert restraints == [[0, 1, 2.1, mc.preprocess.COORDINATION_FORCE]]

    def test_coordination_force_is_read_from_the_config(self, mc, pdb):
        """The literature-based default is a config input, not a buried constant."""
        restraints, _, _ = mc.md.md_restraints(
            {"restraints": [[0, 1]]},
            {"inputs": {"coordination_force": 30.0}, "results": {}}, pdb)
        assert restraints == [[0, 1, 2.1, 30.0]]

    def test_selectors_resolve(self, mc, pdb):
        restraints, _, _ = mc.md.md_restraints(
            {}, {"results": {"coordination": [["B:1:ZN", "A:88:NE2"]]}}, pdb)
        assert restraints == [[0, 1, 2.1, mc.preprocess.COORDINATION_FORCE]]


# The fixture above writes no element column, which is all the distance tests
# need. Element identity is what the D115 checks turn on, so these write it.
def write_site(path, atoms):
    """atoms: (chain, resid, name, element, x, y, z), in file order."""
    with open(path, "w") as fh:
        for n, (chain, resid, name, element, x, y, z) in enumerate(atoms, 1):
            fh.write("HETATM%5d %-4s %-3s %1s%4s    %8.3f%8.3f%8.3f"
                     "                      %2s\n"
                     % (n, name, "LIG", chain, resid, x, y, z, element))
    return str(path)


ZN_NE2 = [("B", "1", "ZN", "ZN", 0.0, 0.0, 0.0),
          ("A", "88", "NE2", "N", 0.0, 0.0, 2.1)]


class TestCoordinationIsAboutThisStructure:
    """D115: `preprocess.results.coordination` held raw atom indices, `md` applied
    them to whatever structure it was given, and a frame built by moving one
    proton between residues had renumbered every atom between the two. The
    restraints meant for the zinc held meropenem's C1x, the zinc ran unrestrained
    through an MM minimisation and lost an anchoring histidine, and nothing in the
    output said so - the run simply looked normal."""

    SECTIONS = {"inputs": {"metals": ["ZN"], "donor_elements": ["N", "O", "S"],
                           "cutoff": 2.8},
                "results": {}}

    def sections(self, coordination=None):
        out = {"inputs": dict(self.SECTIONS["inputs"]), "results": {}}
        if coordination is not None:
            out["results"]["coordination"] = coordination
        return out

    def test_selectors_survive_renumbering(self, mc, tmp_path):
        """What preprocess writes now. The same two rows resolve to different
        indices in two structures, which is the whole point."""
        shifted = write_site(tmp_path / "shifted.pdb",
                             [("Z", "9", "XX", "C", 0.0, 0.0, 9.9)] + ZN_NE2)
        rows = [["B:1:ZN", "A:88:NE2", 2.1, 50.0]]
        here, _, _ = mc.md.md_restraints(
            {}, self.sections(rows), write_site(tmp_path / "here.pdb", ZN_NE2))
        there, _, _ = mc.md.md_restraints({}, self.sections(rows), shifted)
        assert here == [[0, 1, 2.1, 50.0]]
        assert there == [[1, 2, 2.1, 50.0]]      # one atom inserted in front

    def test_indices_that_no_longer_name_a_metal_are_refused(self, mc, tmp_path):
        pdb = write_site(tmp_path / "site.pdb", ZN_NE2)
        """The D115 row shape: the pair still resolves, and names the wrong atoms."""
        with pytest.raises(ValueError) as exc:
            mc.md.md_restraints({}, self.sections([[1, 2, 2.1, 50.0]]), pdb)
        text = str(exc.value)
        assert "index 1 should be a metal" in text and "ZN" in text   # expected
        assert "N HIS" in text or "N" in text                         # and found
        assert "coordination_redetect" in text                        # and the way out

    def test_a_donor_index_that_names_a_hydrogen_is_refused(self, mc, tmp_path):
        """One of D115's recorded donors was the nucleophile's hydrogen."""
        with_h = write_site(tmp_path / "withh.pdb",
                            ZN_NE2 + [("E", "502", "H03", "H", 0.0, 0.0, 1.0)])
        with pytest.raises(ValueError, match="should be a donor"):
            mc.md.md_restraints({}, self.sections([[0, 2, 2.1, 50.0]]), with_h)

    def test_an_index_past_the_end_is_refused(self, mc, tmp_path):
        pdb = write_site(tmp_path / "site.pdb", ZN_NE2)
        with pytest.raises(ValueError, match="past the end"):
            mc.md.md_restraints({}, self.sections([[0, 99, 2.1, 50.0]]), pdb)

    def test_indices_that_still_name_a_metal_and_a_donor_are_used(self, mc, tmp_path):
        """The usual case - the same structure preprocess ran on - is untouched."""
        pdb = write_site(tmp_path / "site.pdb", ZN_NE2)
        restraints, _, source = mc.md.md_restraints(
            {}, self.sections([[0, 1, 2.1, 50.0]]), pdb)
        assert restraints == [[0, 1, 2.1, 50.0]]
        assert source == "preprocess.results.coordination"


class TestRedetectOnThisStructure:
    """The atom order of an edited frame cannot be preserved: the proton has to
    leave one residue and join another, and a PDB residue's atoms must be
    contiguous. So re-detecting on the structure in front of the stage is the only
    way to run one at all - with the same rule and the same inputs preprocess
    used."""

    def sections(self):
        return {"inputs": {"metals": ["ZN"], "donor_elements": ["N", "O", "S"],
                           "cutoff": 2.8, "coordination_force": 30.0},
                "results": {"coordination": [[99, 98, 9.9, 1.0]]}}

    def test_it_detects_on_the_structure_given_not_what_was_saved(self, mc, tmp_path):
        """The saved row names atoms that are not there; re-detection ignores it."""
        here = write_site(tmp_path / "system.pdb", ZN_NE2)
        restraints, detected, source = mc.md.md_restraints(
            {"coordination_redetect": True}, self.sections(), here)
        assert restraints == [[0, 1, 2.1, 30.0]]     # Zn - NE2, as measured here
        assert detected == restraints
        assert "redetected" in source and "system.pdb" in source

    def test_the_distances_are_the_ones_this_structure_has(self, mc, tmp_path):
        """Which is the point, and why the source is recorded: a re-detected run
        is not comparable with one that used the saved numbers."""
        moved = write_site(tmp_path / "moved.pdb",
                           [("B", "1", "ZN", "ZN", 0.0, 0.0, 0.0),
                            ("A", "88", "NE2", "N", 0.0, 0.0, 2.4)])
        restraints, _, _ = mc.md.md_restraints(
            {"coordination_redetect": True}, self.sections(), moved)
        assert restraints == [[0, 1, 2.4, 30.0]]

    def test_the_cutoff_is_the_one_preprocess_used(self, mc, tmp_path):
        """Same rule, not a looser one: a contact beyond `cutoff` is not a contact."""
        far = write_site(tmp_path / "far.pdb",
                         [("B", "1", "ZN", "ZN", 0.0, 0.0, 0.0),
                          ("A", "88", "NE2", "N", 0.0, 0.0, 3.5)])
        restraints, _, _ = mc.md.md_restraints(
            {"coordination_redetect": True}, self.sections(), far)
        assert restraints is None

    def test_it_says_so_when_there_are_no_metals_to_detect(self, mc, tmp_path):
        pdb = write_site(tmp_path / "site.pdb", ZN_NE2)
        with pytest.raises(KeyError, match="metals"):
            mc.md.md_restraints({"coordination_redetect": True},
                                {"inputs": {}, "results": {}}, pdb)


class TestWriteAshSurface:
    """scan.py's ASH-format output, including the resumed-CSV key spellings."""

    @pytest.fixture
    def scan(self, monkeypatch, mc):
        from mcube import scan
        return scan

    def test_1d_rows_are_coordinate_then_energy(self, scan, tmp_path):
        out = scan.write_ash_surface(
            [{"target": 3.0, "energy": -100.5}, {"target": 2.5, "energy": -100.4}],
            False, str(tmp_path / "s.txt"))
        assert [l.split() for l in open(out)] == \
            [["3.00000000", "-100.50000000"], ["2.50000000", "-100.40000000"]]

    def test_2d_rows_carry_both_coordinates(self, scan, tmp_path):
        out = scan.write_ash_surface(
            [{"target": 3.0, "target2": 1.8, "energy": -100.5}], True, str(tmp_path / "s.txt"))
        assert open(out).read().split() == ["3.00000000", "1.80000000", "-100.50000000"]

    def test_points_resumed_from_the_csv_use_its_headers(self, scan, tmp_path):
        out = scan.write_ash_surface(
            [{"target_A": "3.0", "target2_A": "1.8", "energy_Eh": "-100.5"}],
            True, str(tmp_path / "s.txt"))
        assert open(out).read().split() == ["3.00000000", "1.80000000", "-100.50000000"]

    def test_a_point_missing_both_spellings_raises(self, scan, tmp_path):
        with pytest.raises(KeyError):
            scan.write_ash_surface([{"nope": 1}], False, str(tmp_path / "s.txt"))

class TestCoordinationPresent:
    """A metal system with no detected coordination must not run unrestrained."""

    def test_metals_with_coordination_passes(self, mc):
        mc.md.check_coordination_present({"metals": ["ZN"]}, {"coordination": [[0, 1, 2.1, 50.0]]})

    def test_no_metals_is_fine(self, mc):
        """Protein-protein: nothing to hold, so nothing to complain about."""
        mc.md.check_coordination_present({"metals": []}, {})
        mc.md.check_coordination_present({}, {})

    def test_metals_without_coordination_is_refused(self, mc):
        with pytest.raises(KeyError, match="would run unrestrained"):
            mc.md.check_coordination_present({"metals": ["ZN"]}, {})

    def test_the_old_key_is_named_so_the_cause_is_obvious(self, mc):
        """A config predating the rename is the likely cause; say so."""
        with pytest.raises(KeyError) as exc:
            mc.md.check_coordination_present({"metals": ["ZN"]}, {"constraints": [[0, 1]]})
        assert "predates the rename" in str(exc.value)
        assert "Re-run `mcube preprocess`" in str(exc.value)

    def test_md_restraints_refuses_rather_than_silently_dropping_the_site(self, mc, pdb):
        """The whole point: this is reached through the normal path, not only directly."""
        with pytest.raises(KeyError, match="would run unrestrained"):
            mc.md.md_restraints({"restraints": [[0, 1, 2.2]]},
                                {"inputs": {"metals": ["ZN"]},
                                 "results": {"constraints": [[0, 1]]}}, pdb)

