"""The ligand's atom names after parameterising: pytest tests/test_ligand_atom_map.py

OpenFF re-derives ligand atom names from the molecular graph, so the names in
the ligand file a config was written against are not the names in the system
that runs. The hydrolysed meropenem product changed 15 of 53, and three of
those were a permutation among names BOTH files use - `O01`->`O6x`,
`O5x`->`O4x`, `O6x`->`O5x`. A selector written against the input therefore
named a real but different atom: nothing was missing, nothing raised, and a
distance between the wrong pair of oxygens was reported under the right label
(D169).

What these tests pin is the permutation case, since that is the one no other
check can catch, and that a map which cannot be built is reported rather than
guessed at.

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


def atom(serial, name, xyz):
    return ("HETATM%5d %-4s LMP D 401    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
            % (serial, name, xyz[0], xyz[1], xyz[2], name[:1]))


def write(tmp_path, stem, named):
    """A ligand PDB from (name, xyz) pairs."""
    path = tmp_path / f"{stem}.pdb"
    path.write_text("".join(atom(n + 1, name, xyz)
                            for n, (name, xyz) in enumerate(named)) + "END\n")
    return str(path)


# The three oxygens of the product's hydrolysed carboxylate end, at the
# coordinates D159 measured, with the names each file gives them.
INPUT = [("C4x", (0.0, 0.0, 0.0)),
         ("O01", (1.0, 0.0, 0.0)),
         ("O5x", (2.0, 0.0, 0.0)),
         ("O6x", (3.0, 0.0, 0.0))]
PARAMETERISED = [("C4x", (0.0, 0.0, 0.0)),
                 ("O6x", (1.0, 0.0, 0.0)),
                 ("O4x", (2.0, 0.0, 0.0)),
                 ("O5x", (3.0, 0.0, 0.0))]


class TestThePermutation:
    """The failure mode that nothing else catches: every name in the config
    exists in the system, on a different atom."""

    def test_only_the_changed_names_are_reported(self, preprocess, tmp_path):
        got = preprocess.ligand_atom_map(write(tmp_path, "in", INPUT),
                                         write(tmp_path, "par", PARAMETERISED))
        assert got == {"O01": "O6x", "O5x": "O4x", "O6x": "O5x"}

    def test_an_unchanged_name_is_absent(self, preprocess, tmp_path):
        """`C4x` and `N1x` survived, which is why the bond being broken could be
        named the same way in both. An empty entry would imply a change."""
        got = preprocess.ligand_atom_map(write(tmp_path, "in", INPUT),
                                         write(tmp_path, "par", PARAMETERISED))
        assert "C4x" not in got

    def test_the_map_is_by_coordinate_not_by_order(self, preprocess, tmp_path):
        """Order is not a correspondence a renaming cannot fool, and the
        parameteriser is free to reorder. Writing the same atoms backwards must
        give the same map."""
        got = preprocess.ligand_atom_map(
            write(tmp_path, "in", INPUT),
            write(tmp_path, "par", list(reversed(PARAMETERISED))))
        assert got == {"O01": "O6x", "O5x": "O4x", "O6x": "O5x"}


class TestNoRenaming:
    def test_identical_names_give_an_empty_map(self, preprocess, tmp_path):
        """An empty map is the statement that the config's own names are the
        system's names, so it must not be the same value as `could not tell`."""
        got = preprocess.ligand_atom_map(write(tmp_path, "in", INPUT),
                                         write(tmp_path, "par", INPUT))
        assert got == {}

    def test_a_small_coordinate_shift_is_still_the_same_atom(self, preprocess,
                                                             tmp_path):
        moved = [(n, (x + 0.02, y, z)) for n, (x, y, z) in PARAMETERISED]
        got = preprocess.ligand_atom_map(write(tmp_path, "in", INPUT),
                                         write(tmp_path, "par", moved))
        assert got == {"O01": "O6x", "O5x": "O4x", "O6x": "O5x"}


class TestAMapThatCannotBeBuilt:
    """Reported and empty, never guessed. The parameterisation has already
    succeeded by this point, so this is not worth failing a run over - but a
    silent empty map would read as `no renaming happened`."""

    def test_a_different_atom_count_is_reported(self, preprocess, tmp_path, capsys):
        got = preprocess.ligand_atom_map(write(tmp_path, "in", INPUT),
                                         write(tmp_path, "par", PARAMETERISED[:3]))
        assert got == {}
        assert "Cannot map ligand atom names" in capsys.readouterr().out

    def test_an_atom_with_no_match_is_reported(self, preprocess, tmp_path, capsys):
        far = [(n, (x, y, z + 9.0)) for n, (x, y, z) in PARAMETERISED]
        assert preprocess.ligand_atom_map(write(tmp_path, "in", INPUT),
                                          write(tmp_path, "par", far)) == {}
        out = capsys.readouterr().out
        assert "matches 0 atoms" in out

    def test_two_atoms_on_top_of_each_other_are_reported(self, preprocess,
                                                          tmp_path, capsys):
        """An ambiguous match is not a 50/50 guess about which name a selector
        means."""
        doubled = [PARAMETERISED[0], PARAMETERISED[1],
                   ("O4x", (1.0, 0.0, 0.0)), PARAMETERISED[3]]
        assert preprocess.ligand_atom_map(write(tmp_path, "in", INPUT),
                                          write(tmp_path, "par", doubled)) == {}
        assert "matches 2 atoms" in capsys.readouterr().out

    def test_a_missing_file_is_reported(self, preprocess, tmp_path, capsys):
        """`{resname}.pdb` is written by the parameteriser into the working
        directory, so its absence means something earlier went wrong."""
        assert preprocess.ligand_atom_map(write(tmp_path, "in", INPUT),
                                          str(tmp_path / "absent.pdb")) == {}
        assert "no file at" in capsys.readouterr().out

    def test_two_empty_files_are_reported(self, preprocess, tmp_path, capsys):
        """Equal counts, so the count check passes; nothing is mapped, which
        would otherwise look like agreement."""
        empty = tmp_path / "empty.pdb"
        empty.write_text("END\n")
        assert preprocess.ligand_atom_map(str(empty), str(empty)) == {}
        assert "Cannot map ligand atom names" in capsys.readouterr().out


class TestTheTolerance:
    def test_it_is_configurable(self, preprocess, tmp_path, capsys):
        moved = [(n, (x + 0.2, y, z)) for n, (x, y, z) in PARAMETERISED]
        assert preprocess.ligand_atom_map(write(tmp_path, "in", INPUT),
                                          write(tmp_path, "par", moved)) == {}
        got = preprocess.ligand_atom_map(write(tmp_path, "in", INPUT),
                                         write(tmp_path, "par", moved),
                                         tolerance=0.5)
        assert got == {"O01": "O6x", "O5x": "O4x", "O6x": "O5x"}

    def test_the_default_is_tight_enough_to_separate_bonded_atoms(self, preprocess):
        """A bond is ~1.0 A at the shortest, so the default must be far below
        it or a hydrogen could match its own heavy atom."""
        assert preprocess.NAME_MAP_TOLERANCE < 0.5
