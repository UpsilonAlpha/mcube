"""A missing preprocess input says where it looked: pytest tests/test_ligand_inputs.py

RDKit answers a missing file with `OSError: Bad input file <name>`, which reads
like a malformed structure and sends you to check the PDB columns. It is nearly
always the directory instead: `run.slurm` does `cd "$OUT"` before calling mcube,
so `preprocess.inputs.protein_pdb` and `.ligand` are read from the stage's
working directory while every other stage reads `pdb_file` against the config's.
Job 31299717 died on that twice, six seconds each (D161).

Needs rdkit, which comes with the ASH environment; skipped elsewhere.
"""

import sys
import types

import pytest

pytest.importorskip("rdkit", reason="preprocess reads ligands through RDKit")
pytest.importorskip("openmm", reason="preprocess imports OpenMM")

ASH_NAMES = ("Fragment", "OpenMM_Modeller", "small_molecule_parameterizer",
             "merge_pdb_files")


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


class TestAMissingLigand:
    def test_it_raises_FileNotFoundError_not_OSError(self, preprocess, tmp_path,
                                                     monkeypatch):
        monkeypatch.chdir(tmp_path)
        with pytest.raises(FileNotFoundError):
            preprocess.ligand_mol_from_pdb("./gone.pdb", 0)

    def test_the_message_names_the_working_directory(self, preprocess, tmp_path,
                                                     monkeypatch):
        """The directory is the answer, so it has to be in the message."""
        monkeypatch.chdir(tmp_path)
        with pytest.raises(FileNotFoundError) as exc:
            preprocess.ligand_mol_from_pdb("./gone.pdb", 0)
        assert str(tmp_path) in str(exc.value)

    def test_the_message_contrasts_it_with_pdb_file(self, preprocess, tmp_path,
                                                    monkeypatch):
        monkeypatch.chdir(tmp_path)
        with pytest.raises(FileNotFoundError) as exc:
            preprocess.ligand_mol_from_pdb("./gone.pdb", 0)
        assert "pdb_file" in str(exc.value)


class TestALigandWithoutBonds:
    def test_no_CONECT_is_refused_by_name(self, preprocess, tmp_path, monkeypatch):
        """The other way this file goes wrong, and the reason proximity bonding
        is off: distance-derived bonds put a spurious O-O into this very
        molecule once (D21)."""
        monkeypatch.chdir(tmp_path)
        (tmp_path / "bare.pdb").write_text(
            "HETATM    1  C1  LIG L   1       0.000   0.000   0.000  1.00  0.00"
            "           C\n"
            "HETATM    2  O1  LIG L   1       1.250   0.000   0.000  1.00  0.00"
            "           O\n"
            "END\n")
        with pytest.raises(ValueError) as exc:
            preprocess.ligand_mol_from_pdb("./bare.pdb", 0)
        assert "CONECT" in str(exc.value)

    def test_CONECT_is_enough_to_get_a_molecule(self, preprocess, tmp_path,
                                                monkeypatch):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "ok.pdb").write_text(
            "HETATM    1  C1  LIG L   1       0.000   0.000   0.000  1.00  0.00"
            "           C\n"
            "HETATM    2  O1  LIG L   1       1.250   0.000   0.000  1.00  0.00"
            "           O\n"
            "CONECT    1    2\n"
            "END\n")
        mol = preprocess.ligand_mol_from_pdb("./ok.pdb", 0)
        assert mol.GetNumAtoms() == 2 and mol.GetNumBonds() == 1
