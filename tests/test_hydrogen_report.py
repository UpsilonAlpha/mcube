"""Hydrogens, read off a structure: pytest tests/test_hydrogen_report.py

Hydrogens are the one part of a system no crystal structure supplies - every
one was placed by pdbfixer, Modeller, OpenFF or a graft - and nothing
downstream refuses a wrongly placed one. A missing amide H changes a backbone,
an extra carboxylate H changes a charge, and an H sitting on another atom is a
NaN waiting for the integrator to find it: the product arm's 100 ns MD died at
34.14 ns that way, with every template matched and no warning anywhere (D172).

What these tests pin is that the check needs no force field and no
connectivity, so it runs on a built system, on a trajectory frame, or on
someone else's PDB; and that it distinguishes the states a name alone cannot -
HID from HIE, ASP from ASH.

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


def atom(serial, name, resname, chain, seq, xyz):
    element = name[0] if name[0] != "H" else "H"
    return ("ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
            % (serial, name, resname, chain, seq, xyz[0], xyz[1], xyz[2], element))


def write(tmp_path, atoms, name="s.pdb"):
    path = tmp_path / name
    path.write_text(
        "CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1           1\n"
        + "".join(atom(n + 1, *a) for n, a in enumerate(atoms)) + "END\n")
    return str(path)


def residue(resname, chain, seq, extra=(), x=0.0):
    """A minimal residue: N with its amide H, CA, and whatever else is given."""
    out = [("N", resname, chain, seq, (x, 0.0, 0.0)),
           ("H", resname, chain, seq, (x, 1.01, 0.0)),
           ("CA", resname, chain, seq, (x + 1.5, 0.0, 0.0))]
    out += [(n, resname, chain, seq, p) for n, p in extra]
    return out


class TestProtonationStates:
    """A residue name does not say which tautomer it is. `residue_variants`
    asks for one; this reports what the file actually has."""

    def test_a_delta_histidine_reads_as_hid(self, preprocess, tmp_path):
        atoms = residue("HIS", "A", 88, [("ND1", (3.0, 0.0, 0.0)),
                                         ("HD1", (3.0, 1.01, 0.0)),
                                         ("NE2", (4.5, 0.0, 0.0))])
        rep = preprocess.hydrogen_report(write(tmp_path, atoms))
        his = [s for s in rep["protonation"] if s["residue"].startswith("HIS")][0]
        assert his["present"] == ["HD1"] and his["absent"] == ["HE2"]
        assert "HID" in his["reads_as"]

    def test_both_hydrogens_read_as_doubly_protonated(self, preprocess, tmp_path):
        atoms = residue("HIS", "A", 88, [("ND1", (3.0, 0.0, 0.0)),
                                         ("HD1", (3.0, 1.01, 0.0)),
                                         ("NE2", (4.5, 0.0, 0.0)),
                                         ("HE2", (4.5, 1.01, 0.0))])
        rep = preprocess.hydrogen_report(write(tmp_path, atoms))
        his = [s for s in rep["protonation"] if s["residue"].startswith("HIS")][0]
        assert his["present"] == ["HD1", "HE2"] and his["absent"] == []

    def test_a_bare_aspartate_is_not_reported_as_protonated(self, preprocess, tmp_path):
        """Asp92's charge is the whole question on this site: a stray HD2 would
        remove the electrostatic reason for the contact being measured."""
        atoms = residue("ASP", "A", 92, [("OD1", (3.0, 0.0, 0.0)),
                                         ("OD2", (4.5, 0.0, 0.0))])
        rep = preprocess.hydrogen_report(write(tmp_path, atoms))
        asp = [s for s in rep["protonation"] if s["residue"].startswith("ASP")][0]
        assert asp["present"] == [] and asp["absent"] == ["HD2"]

    def test_a_protonated_aspartate_says_so(self, preprocess, tmp_path):
        atoms = residue("ASP", "A", 92, [("OD1", (3.0, 0.0, 0.0)),
                                         ("OD2", (4.5, 0.0, 0.0)),
                                         ("HD2", (4.5, 1.01, 0.0))])
        rep = preprocess.hydrogen_report(write(tmp_path, atoms))
        asp = [s for s in rep["protonation"] if s["residue"].startswith("ASP")][0]
        assert "ASH" in asp["reads_as"]


class TestMissingHydrogens:
    def test_an_amide_nitrogen_without_its_h_is_reported(self, preprocess, tmp_path):
        atoms = [("N", "ALA", "A", 10, (0.0, 0.0, 0.0)),
                 ("CA", "ALA", "A", 10, (1.5, 0.0, 0.0)),
                 ("HA", "ALA", "A", 10, (1.5, 1.01, 0.0))]
        rep = preprocess.hydrogen_report(write(tmp_path, atoms))
        assert rep["amide_n_without_h"] == ["ALA A:10"]

    def test_proline_has_no_amide_h_and_is_not_reported(self, preprocess, tmp_path):
        atoms = [("N", "PRO", "A", 11, (0.0, 0.0, 0.0)),
                 ("CA", "PRO", "A", 11, (1.5, 0.0, 0.0)),
                 ("HA", "PRO", "A", 11, (1.5, 1.01, 0.0))]
        rep = preprocess.hydrogen_report(write(tmp_path, atoms))
        assert rep["amide_n_without_h"] == []

    def test_a_terminal_h1_counts_as_the_amide_h(self, preprocess, tmp_path):
        atoms = [("N", "ALA", "A", 1, (0.0, 0.0, 0.0)),
                 ("H1", "ALA", "A", 1, (0.0, 1.01, 0.0)),
                 ("CA", "ALA", "A", 1, (1.5, 0.0, 0.0))]
        rep = preprocess.hydrogen_report(write(tmp_path, atoms))
        assert rep["amide_n_without_h"] == []

    def test_a_residue_with_no_hydrogens_at_all_is_reported(self, preprocess, tmp_path):
        """A PDB straight from the bank is nothing but these: a residue whose
        hydrogens were never added, which every later stage will type as though
        they were."""
        atoms = [("N", "ALA", "A", 10, (0.0, 0.0, 0.0)),
                 ("CA", "ALA", "A", 10, (1.5, 0.0, 0.0)),
                 ("C", "ALA", "A", 10, (3.0, 0.0, 0.0))]
        rep = preprocess.hydrogen_report(write(tmp_path, atoms))
        assert rep["residues_without_hydrogens"] == ["ALA A:10"]


class TestGeometry:
    def test_an_unattached_hydrogen_is_reported(self, preprocess, tmp_path):
        """Nothing within bonding distance: a placement that went into solvent."""
        atoms = residue("ALA", "A", 10) + [("HB", "ALA", "A", 10, (20.0, 20.0, 20.0))]
        rep = preprocess.hydrogen_report(write(tmp_path, atoms))
        assert [o["hydrogen"] for o in rep["unattached_hydrogens"]] == ["ALA A:10:HB"]

    def test_a_hydrogen_between_two_heavy_atoms_is_reported(self, preprocess, tmp_path):
        """This is the NaN mechanism: an H close enough to two heavy atoms is
        either a clash or an ambiguous assignment, and the integrator finds it
        long after the minimiser declared success."""
        atoms = [("N", "ALA", "A", 10, (0.0, 0.0, 0.0)),
                 ("H", "ALA", "A", 10, (1.0, 0.0, 0.0)),
                 ("CA", "ALA", "A", 10, (2.0, 0.0, 0.0))]
        rep = preprocess.hydrogen_report(write(tmp_path, atoms))
        assert rep["hydrogen_clashes"]
        assert rep["hydrogen_clashes"][0]["hydrogen"] == "ALA A:10:H"
        assert len(rep["hydrogen_clashes"][0]["near"]) == 2

    def test_a_normal_hydrogen_is_not_reported(self, preprocess, tmp_path):
        rep = preprocess.hydrogen_report(write(tmp_path, residue("ALA", "A", 10)))
        assert rep["unattached_hydrogens"] == [] and rep["hydrogen_clashes"] == []

    def test_the_minimum_image_is_applied(self, preprocess, tmp_path):
        """An H at one box edge and its heavy atom at the other are bonded, and
        a check that said otherwise would report every edge atom."""
        atoms = [("N", "ALA", "A", 10, (0.5, 0.0, 0.0)),
                 ("H", "ALA", "A", 10, (59.9, 0.0, 0.0)),
                 ("CA", "ALA", "A", 10, (2.0, 0.0, 0.0))]
        rep = preprocess.hydrogen_report(write(tmp_path, atoms))
        assert rep["unattached_hydrogens"] == []


class TestTheLigandFormula:
    """The check that catches a ligand hydrogenated to the wrong state, which
    no template match refuses - OpenFF types whatever it is given."""

    LIG = [("C1", "LIG", "D", 1, (0.0, 0.0, 0.0)),
           ("O1", "LIG", "D", 1, (1.3, 0.0, 0.0)),
           ("H1", "LIG", "D", 1, (1.3, 1.01, 0.0)),
           ("H2", "LIG", "D", 1, (0.0, 1.01, 0.0)),
           ("H3", "LIG", "D", 1, (0.0, -1.01, 0.0)),
           ("H4", "LIG", "D", 1, (-1.01, 0.0, 0.0))]

    def test_a_matching_smiles_agrees(self, preprocess, tmp_path):
        pytest.importorskip("rdkit")
        rep = preprocess.hydrogen_report(write(tmp_path, self.LIG),
                                         ligand_resname="LIG", ligand_smiles="CO")
        assert rep["ligand_formula_check"]["agrees"]
        assert rep["ligand_formula_check"]["from_smiles"] == "C1H4O1"

    def test_a_missing_hydrogen_disagrees(self, preprocess, tmp_path):
        """Methanol less one H is methoxide: same atoms to a template matcher,
        a different charge to the SCF."""
        pytest.importorskip("rdkit")
        rep = preprocess.hydrogen_report(write(tmp_path, self.LIG[:-1]),
                                         ligand_resname="LIG", ligand_smiles="CO")
        chk = rep["ligand_formula_check"]
        assert not chk["agrees"] and chk["copies"] == ["C1H3O1"]

    def test_the_formula_is_reported_without_a_smiles(self, preprocess, tmp_path):
        rep = preprocess.hydrogen_report(write(tmp_path, self.LIG))
        assert rep["nonstandard_residues"][0]["formula"] == "C1H4O1"
        assert "ligand_formula_check" not in rep

    def test_a_named_ligand_that_is_absent_is_reported(self, preprocess, tmp_path):
        pytest.importorskip("rdkit")
        rep = preprocess.hydrogen_report(write(tmp_path, self.LIG),
                                         ligand_resname="XXX", ligand_smiles="CO")
        assert rep["ligand_formula_check"]["agrees"] is False
        assert "no residue named XXX" in rep["ligand_formula_check"]["status"]


def test_the_report_prints(preprocess, tmp_path, capsys):
    """`report_hydrogens` is what a run actually shows, so it has to survive a
    report with every list empty as well as a full one."""
    rep = preprocess.hydrogen_report(write(tmp_path, residue("ALA", "A", 10)))
    preprocess.report_hydrogens(rep)
    out = capsys.readouterr().out
    assert "Every amide N has its H" in out


class TestMonatomicIons:
    """The product system has a zinc and 56 sodiums. Reported as residues
    without hydrogens they were 57 lines burying the one that mattered - and a
    monatomic ion has no hydrogens to be missing (D173)."""

    def test_an_ion_is_not_a_residue_without_hydrogens(self, preprocess, tmp_path):
        atoms = [("ZN", "ZN", "B", 1, (0.0, 0.0, 0.0)),
                 ("NA", "NA", "F", 1, (9.0, 0.0, 0.0))]
        rep = preprocess.hydrogen_report(write(tmp_path, atoms))
        assert rep["residues_without_hydrogens"] == []
        assert rep["monatomic_ions"] == {"ZN": 1, "NA": 1}

    def test_they_are_counted_not_listed(self, preprocess, tmp_path):
        atoms = [("NA", "NA", "F", n, (9.0 * n, 0.0, 0.0)) for n in range(1, 6)]
        rep = preprocess.hydrogen_report(write(tmp_path, atoms))
        assert rep["monatomic_ions"] == {"NA": 5}
        assert rep["nonstandard_residues"] == []

    def test_a_lone_water_oxygen_is_an_ion_by_this_rule(self, preprocess, tmp_path):
        """The honest limit of counting atoms: one oxygen named HOH is
        indistinguishable from an ion here, so it is counted rather than
        warned about. Water missing its hydrogens is caught at the template
        match instead, where tip3p has no one-atom residue to match."""
        atoms = [("O", "HOH", "C", 1, (0.0, 0.0, 0.0))]
        rep = preprocess.hydrogen_report(write(tmp_path, atoms))
        assert rep["monatomic_ions"] == {"HOH": 1}
