"""Residue numbers past 9999 in OpenMM-written PDBs (membrane boxes have ~20k waters)."""

import io

import pytest

pytest.importorskip("ash")
from openmm import Vec3, unit                     # noqa: E402
from openmm.app import Element, PDBFile, Topology  # noqa: E402

from mcube.preprocess import atom_names, pdb_resseq  # noqa: E402


def test_decodes_what_openmm_writes(tmp_path):
    top = Topology()
    chain = top.addChain()
    pos = []
    for i in range(12000):
        res = top.addResidue("HOH", chain)
        top.addAtom("O", Element.getBySymbol("O"), res)
        pos.append(Vec3(0.1 * (i % 50), 0.1 * ((i // 50) % 50), 0.1 * (i // 2500)))
    path = tmp_path / "box.pdb"
    with open(path, "w") as fh:
        PDBFile.writeFile(top, pos * unit.nanometer, fh)
    lines = [l for l in open(path) if l.startswith(("ATOM", "HETATM"))]
    ours = [pdb_resseq(l[22:26]) for l in lines]
    theirs = [int(r.id) for r in PDBFile(str(path)).topology.residues()]
    assert ours == theirs
    assert "A1F5" in lines[10500] and ours[10500] == 10501
    assert atom_names(str(path))[10500].split(":")[1] == "10501"


def test_plain_numbers():
    assert pdb_resseq("   7") == 7
    assert pdb_resseq("9999") == 9999
