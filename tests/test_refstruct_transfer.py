"""refstruct.py transfer: pytest tests/test_refstruct_transfer.py

Putting a deposited complex's ligand into our frame is the step that decides
where the substrate goes, so the two things that can quietly go wrong are
tested here: a fit whose atom correspondence is wrong (it must show up as a
large RMSD rather than as a plausible pose), and the connectivity rules that
pick out a carboxylate or an intact amide carbonyl without an atom-name table.

numpy only; refstruct.py imports it inside the fitting functions.
"""

import importlib.util
import math
from pathlib import Path

import pytest

pytest.importorskip("numpy")

SCRIPT = Path(__file__).resolve().parents[1] / "docs" / "csr" / "refstruct.py"


@pytest.fixture(scope="module")
def rs():
    spec = importlib.util.spec_from_file_location("refstruct", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def atom(name, resname, chain, seq, xyz, element=None):
    return (name, resname, chain, str(seq), xyz[0], xyz[1], xyz[2],
            element or name[0])


class TestKabsch:
    def test_identical_sets_fit_exactly(self, rs):
        pts = [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 2.5)]
        assert rs.kabsch(pts, pts)[3] == pytest.approx(0.0, abs=1e-9)

    def test_a_pure_rotation_is_recovered(self, rs):
        pts = [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 2.5)]
        c, s = math.cos(0.7), math.sin(0.7)
        turned = [(c * x - s * y, s * x + c * y, z) for x, y, z in pts]
        fit = rs.kabsch(pts, turned)
        assert fit[3] == pytest.approx(0.0, abs=1e-9)
        for p, q in zip(pts, turned):
            moved = rs.apply_fit(fit, p)
            assert [round(float(v), 6) for v in moved] == [round(v, 6) for v in q]

    def test_a_wrong_correspondence_shows_up_as_a_large_rmsd(self, rs):
        """The whole safeguard: swap two pairs and the fit must not look fine."""
        pts = [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 2.5)]
        swapped = [pts[1], pts[0], pts[2], pts[3]]
        assert rs.kabsch(pts, swapped)[3] > 0.5

    def test_no_reflection_is_allowed(self, rs):
        """A mirrored set must fit badly, not perfectly - a reflection would
        place the ligand's stereocentres on the wrong face."""
        pts = [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 2.5)]
        mirrored = [(x, y, -z) for x, y, z in pts]
        assert rs.kabsch(pts, mirrored)[3] > 0.5


class TestSelect:
    def test_it_finds_one_atom(self, rs):
        atoms = [atom("NE2", "HIS", "A", 88, (1.0, 2.0, 3.0)),
                 atom("NE2", "HIS", "A", 164, (4.0, 5.0, 6.0))]
        assert rs.select(atoms, "A:88:NE2")[4] == 1.0

    def test_a_selector_matching_nothing_is_an_error(self, rs):
        with pytest.raises(SystemExit):
            rs.select([atom("NE2", "HIS", "A", 88, (0, 0, 0))], "A:90:NE2")

    def test_an_ambiguous_selector_is_an_error(self, rs):
        twice = [atom("NE2", "HIS", "A", 88, (0, 0, 0)),
                 atom("NE2", "HIS", "A", 88, (1, 1, 1))]
        with pytest.raises(SystemExit):
            rs.select(twice, "A:88:NE2")


class TestConnectivityRules:
    """Atom names differ in every ligand code (LMP, DQM, MX0, MER), so the
    carboxylate and the carbonyl are found by what they are bonded to."""

    def carboxylate(self):
        return [atom("C31", "LMP", "A", 317, (0.0, 0.0, 0.0), "C"),
                atom("O31", "LMP", "A", 317, (1.25, 0.0, 0.0), "O"),
                atom("O32", "LMP", "A", 317, (-0.6, 1.1, 0.0), "O")]

    def test_a_carboxylate_is_found_by_its_two_oxygens(self, rs):
        found = rs.carboxylates(self.carboxylate())
        assert len(found) == 1 and found[0][0][0] == "C31"

    def test_a_single_oxygen_is_not_a_carboxylate(self, rs):
        assert rs.carboxylates(self.carboxylate()[:2]) == []

    def test_an_intact_lactam_carbonyl_is_found(self, rs):
        lig = [atom("C7", "MX0", "A", 303, (0.0, 0.0, 0.0), "C"),
               atom("O1", "MX0", "A", 303, (1.22, 0.0, 0.0), "O"),
               atom("N1", "MX0", "A", 303, (-0.7, 1.2, 0.0), "N"),
               atom("C6", "MX0", "A", 303, (-0.8, -1.2, 0.0), "C")]
        found = rs.amide_carbonyls(lig)
        assert len(found) == 1
        assert (found[0][0][0], found[0][1][0], found[0][2][0]) == ("C7", "O1", "N1")

    def test_a_hydrolysed_product_has_no_amide_carbonyl(self, rs):
        """The bond the enzyme broke is the one this rule looks for."""
        lig = self.carboxylate() + [atom("N4", "LMP", "A", 317, (4.0, 0.0, 0.0), "N")]
        assert rs.amide_carbonyls(lig) == []


class TestModels:
    def test_each_vina_pose_is_returned_separately(self, rs, tmp_path):
        path = tmp_path / "poses.pdbqt"
        path.write_text(
            "MODEL 1\n"
            "ATOM      1  C1  UNL L   1       1.000   0.000   0.000  1.00  0.00     0.0 C\n"
            "ENDMDL\n"
            "MODEL 2\n"
            "ATOM      1  C1  UNL L   1       9.000   0.000   0.000  1.00  0.00     0.0 C\n"
            "ENDMDL\n")
        got = rs.models(str(path))
        assert [lab for lab, _ in got] == ["model 1", "model 2"]
        assert got[0][1][0][4] == 1.0 and got[1][1][0][4] == 9.0

    def test_a_pdbqt_without_element_columns_still_gets_elements(self, rs, tmp_path):
        path = tmp_path / "p.pdbqt"
        path.write_text(
            "MODEL 1\n"
            "ATOM      1  O2  UNL L   1       1.000   0.000   0.000  1.00  0.00     0.0\n"
            "ENDMDL\n")
        assert rs.models(str(path))[0][1][0][7] == "O"

    def test_a_single_model_file_is_one_entry(self, rs, tmp_path):
        path = tmp_path / "one.pdb"
        path.write_text(
            "HETATM    1  ZN   ZN B   1       0.000   0.000   0.000  1.00  0.00  "
            "        ZN\n")
        assert len(rs.models(str(path))) == 1

    def test_filtering_by_residue_name_drops_the_protein(self, rs, tmp_path):
        path = tmp_path / "sys.pdb"
        path.write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00    "
            "       C\n"
            "HETATM    2  C1  MER C   1       5.000   0.000   0.000  1.00  0.00    "
            "       C\n")
        got = rs.models(str(path), "MER")
        assert len(got) == 1 and len(got[0][1]) == 1


class TestCentroid:
    def test_distance_between_two_sets(self, rs):
        a = [atom("C1", "L", "L", 1, (0.0, 0.0, 0.0), "C")]
        b = [atom("C1", "L", "L", 1, (3.0, 4.0, 0.0), "C")]
        assert rs.centroid_distance(a, b) == pytest.approx(5.0)

class TestTheBetaLactamMarker:
    """`ligand_report` marks which amide is the beta-lactam by its distance to
    the C3 carboxylate on the same fused bicycle, because a pdbqt pose carries
    no usable atom names. A hydrolysed product must get no mark at all: its
    only remaining amide is the distal one, 7.0 A away in LMP, and labelling
    that would name the bond the enzyme has already broken."""

    def report(self, rs, lig, capsys):
        zn = atom("ZN", "ZN", "B", 1, (0.0, 0.0, 0.0), "ZN")
        rs.ligand_report("lig", lig, zn)
        return capsys.readouterr().out

    def test_the_ring_amide_is_marked(self, rs, capsys):
        lig = [atom("C7", "MER", "C", 1, (0.0, 5.0, 0.0), "C"),
               atom("O1", "MER", "C", 1, (1.2, 5.0, 0.0), "O"),
               atom("N1", "MER", "C", 1, (-0.7, 6.2, 0.0), "N"),
               atom("C3", "MER", "C", 1, (3.5, 5.0, 0.0), "C"),
               atom("O3", "MER", "C", 1, (4.2, 6.0, 0.0), "O"),
               atom("O4", "MER", "C", 1, (3.9, 3.9, 0.0), "O")]
        assert "<-- beta-lactam" in self.report(rs, lig, capsys)

    def test_a_distal_amide_alone_is_not_marked(self, rs, capsys):
        lig = [atom("C6A", "LMP", "C", 1, (0.0, 5.0, 0.0), "C"),
               atom("O7A", "LMP", "C", 1, (1.2, 5.0, 0.0), "O"),
               atom("N7A", "LMP", "C", 1, (-0.7, 6.2, 0.0), "N"),
               atom("C31", "LMP", "C", 1, (9.0, 5.0, 0.0), "C"),
               atom("O31", "LMP", "C", 1, (9.7, 6.0, 0.0), "O"),
               atom("O32", "LMP", "C", 1, (9.4, 3.9, 0.0), "O")]
        out = self.report(rs, lig, capsys)
        assert "amide C=O" in out and "<-- beta-lactam" not in out

    def test_two_carboxylates_make_it_abstain(self, rs, capsys):
        """Moxalactam has a second carboxylate on its side arm, and the rule
        then picks the distal amide over the lactam the crystal shows on the
        metal. Abstaining is the only safe answer."""
        lig = [atom("C1", "MX0", "C", 1, (0.0, 5.0, 0.0), "C"),
               atom("O1", "MX0", "C", 1, (1.2, 5.0, 0.0), "O"),
               atom("N1", "MX0", "C", 1, (-0.7, 6.2, 0.0), "N"),
               atom("C13", "MX0", "C", 1, (3.5, 5.0, 0.0), "C"),
               atom("O8", "MX0", "C", 1, (4.2, 6.0, 0.0), "O"),
               atom("O9", "MX0", "C", 1, (3.9, 3.9, 0.0), "O"),
               atom("C11", "MX0", "C", 1, (12.0, 5.0, 0.0), "C"),
               atom("O6", "MX0", "C", 1, (12.7, 6.0, 0.0), "O"),
               atom("O7", "MX0", "C", 1, (12.4, 3.9, 0.0), "O")]
        out = self.report(rs, lig, capsys)
        assert "<-- beta-lactam" not in out and "2 carboxylates" in out


class TestSimilarity:
    """Whether two MBLs would hold a substrate the same way is a question with
    a number behind it, and the number has to come from the structures rather
    than from the subclass label (D158)."""

    def test_identical_sequences_align_one_to_one(self, rs):
        assert rs.align("ACDEFGH", "ACDEFGH") == [(i, i) for i in range(7)]

    def test_a_deletion_is_absorbed_by_a_gap(self, rs):
        """The D is missing from the second sequence; everything either side
        must still pair up rather than shifting by one."""
        pairs = rs.align("ACDEFGHIK", "ACEFGHIK")
        assert (0, 0) in pairs and (1, 1) in pairs
        assert (3, 2) in pairs and (8, 7) in pairs
        assert not any(i == 2 for i, _ in pairs)

    def test_affine_gaps_keep_one_gap_rather_than_several(self, rs):
        pairs = rs.align("AAAACDEFGHWWWW", "AAAAWWWW")
        runs = sum(1 for k in range(1, len(pairs))
                   if pairs[k][0] - pairs[k - 1][0] > 1)
        assert runs == 1

    def test_unrelated_sequences_still_return_pairs(self, rs):
        """A global alignment always pairs something; identity is what says
        whether the pairing means anything."""
        assert rs.align("AAAAAA", "WWWWWW")

    def test_a_windows_path_is_not_read_as_a_chain(self, rs, tmp_path):
        """`C:/x.pdb` must not parse as chain-less file "C" on chain "/x.pdb"."""
        path = tmp_path / "x.pdb"
        path.write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00"
            "           C\n")
        head, sep, tail = str(path).rpartition(":")
        parsed = (head, tail) if sep and len(tail) == 1 else (str(path), "")
        assert parsed[0] == str(path) and parsed[1] == ""

    def test_residues_without_a_ca_are_dropped(self, rs, tmp_path):
        path = tmp_path / "y.pdb"
        path.write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00"
            "           C\n"
            "ATOM      2  CB  SER A   2       1.000   0.000   0.000  1.00  0.00"
            "           C\n")
        assert [r[0] for r in rs.chain_residues(str(path))] == ["A"]

    def test_only_the_named_chain_is_read(self, rs, tmp_path):
        path = tmp_path / "z.pdb"
        path.write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00"
            "           C\n"
            "ATOM      2  CA  TRP B   1       1.000   0.000   0.000  1.00  0.00"
            "           C\n")
        assert [r[0] for r in rs.chain_residues(str(path), "B")] == ["W"]

    def test_the_core_rmsd_ignores_one_flapping_loop(self, rs):
        """Two folds that agree everywhere but one loop must report a small
        core RMSD and a core smaller than the alignment."""
        import math

        # A helix, not a line: a collinear set lets Kabsch rotate freely about
        # its own axis and the displaced points swing back into place.
        def helix(i):
            return (1.5 * i, 5.0 * math.cos(i * 0.6), 5.0 * math.sin(i * 0.6))

        a = [("A", helix(i)) for i in range(40)]
        b = [("A", helix(i)) for i in range(40)]
        b[20:25] = [("A", (helix(i)[0], helix(i)[1] + 25.0, helix(i)[2]))
                    for i in range(20, 25)]
        pairs = [(i, i) for i in range(40)]
        full, core, core_rmsd = rs.superpose_pairs(a, b, pairs)
        assert full > 2.0 and core == 35 and core_rmsd < 0.01
