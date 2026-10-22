"""NEB refuses endpoints that are not one system: pytest tests/test_neb_endpoints.py

Runs without ASH, which neb.py imports at module scope, by importing mcube
against a stub - the same approach as tests/test_umbrella.py. Needs OpenMM,
because atom_names reads the structure through it.
"""

import sys
import types

import pytest

pytest.importorskip("openmm", reason="atom_names reads structures through OpenMM")

ASH_NAMES = ("Fragment", "OpenMMTheory", "QMMMTheory", "xTBTheory", "NEB", "NEBTS",
             "Singlepoint", "OpenMM_Opt", "OpenMM_MD", "OpenMM_Modeller",
             "OpenMM_box_equilibration", "Gentle_warm_up_MD",
             "small_molecule_parameterizer", "merge_pdb_files", "Optimizer", "NumFreq")

# Two waters and a zinc: enough for an order mismatch that preserves the count.
REACTANT = [("B", "1", "ZN", "ZN"), ("E", "1", "HOH", "O"), ("E", "2", "HOH", "O")]
SWAPPED = [("B", "1", "ZN", "ZN"), ("E", "2", "HOH", "O"), ("E", "1", "HOH", "O")]


@pytest.fixture
def neb(monkeypatch):
    stub = types.ModuleType("ash")
    for name in ASH_NAMES:
        setattr(stub, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", stub)
    for name in [m for m in sys.modules if m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    from mcube import neb
    return neb


class Frag:
    """Stands in for an ASH Fragment, which only numatoms is read from."""

    def __init__(self, n):
        self.numatoms = n


def write(path, atoms):
    with open(path, "w") as fh:
        for n, (chain, resid, resname, name) in enumerate(atoms, 1):
            fh.write("HETATM%5d %-4s %-3s %1s%4s    %8.3f%8.3f%8.3f\n"
                     % (n, name, resname, chain, resid, 0.0, 0.0, float(n)))
    return str(path)


def test_same_system_passes(neb, tmp_path):
    a = write(tmp_path / "r.pdb", REACTANT)
    b = write(tmp_path / "p.pdb", REACTANT)
    neb.check_same_system(a, b, Frag(3), Frag(3))


def test_different_atom_count_is_refused(neb, tmp_path):
    a = write(tmp_path / "r.pdb", REACTANT)
    b = write(tmp_path / "p.pdb", REACTANT[:2])
    with pytest.raises(ValueError, match="NEB needs one system"):
        neb.check_same_system(a, b, Frag(3), Frag(2))


def test_same_count_but_reordered_is_refused(neb, tmp_path):
    """The case an atom-count check cannot see, and the one that matters."""
    a = write(tmp_path / "r.pdb", REACTANT)
    b = write(tmp_path / "p.pdb", SWAPPED)
    with pytest.raises(ValueError, match="different order"):
        neb.check_same_system(a, b, Frag(3), Frag(3))


def test_the_message_names_the_first_offending_atoms(neb, tmp_path):
    a = write(tmp_path / "r.pdb", REACTANT)
    b = write(tmp_path / "p.pdb", SWAPPED)
    with pytest.raises(ValueError) as exc:
        neb.check_same_system(a, b, Frag(3), Frag(3))
    text = str(exc.value)
    assert "atom 1: reactant E:1:O, product E:2:O" in text
    assert "2 positions differ" in text
    assert "scan" in text          # points at how to build a product that passes
