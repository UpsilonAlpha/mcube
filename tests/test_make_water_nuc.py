"""ohH input structure: pytest tests/test_make_water_nuc.py

The whole protonation test rests on the two arms differing by one hydrogen and
nothing else, so the geometry this places and the refusals it makes are the
things worth pinning. A hydrogen dropped into a clash, or a residue left named
OH, would both survive to the cluster and fail there hours later.
"""

import importlib.util
import math
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "docs" / "csr" / "make_water_nuc.py"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("make_water_nuc", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def pdb_line(serial, name, resname, chain, resseq, xyz):
    return (f"HETATM{serial:5d} {name:<4s}{resname:>4s} {chain}{resseq:>4s}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00\n")


@pytest.fixture
def lysH(tmp_path):
    """A minimal lysH input: a zinc, a hydroxide, and one distant spectator.

    With CONECT records, as the real input has: the hydroxide is not a standard
    residue, so PDBFile cannot infer its bond from atom names."""
    lines = [
        pdb_line(1, "ZN", "ZN", "B", "1", (0.0, 0.0, 0.0)),
        pdb_line(2, "O01", "OH", "C", "2", (1.9, 0.0, 0.0)),
        pdb_line(3, "H03", "OH", "C", "2", (2.4, 0.83, 0.0)),
        pdb_line(4, "CA", "GLY", "A", "9", (12.0, 0.0, 0.0)),
        "CONECT    2    3\n",
        "CONECT    3    2\n",
        "END\n",
    ]
    path = tmp_path / "in.pdb"
    path.write_text("".join(lines))
    return path


# ---- the geometry -----------------------------------------------------------

def test_h04_has_tip3p_geometry(mod):
    zn, o01, h03 = (0.0, 0.0, 0.0), (1.9, 0.0, 0.0), (2.4, 0.83, 0.0)
    h04 = mod.place_h04(zn, o01, h03)
    assert mod.norm(mod.sub(h04, o01)) == pytest.approx(mod.OH_LEN, abs=1e-6)
    a = mod.unit(mod.sub(h03, o01))
    b = mod.unit(mod.sub(h04, o01))
    angle = math.degrees(math.acos(mod.dot(a, b)))
    assert angle == pytest.approx(mod.HOH_ANGLE, abs=1e-4)


def test_h04_is_in_the_zn_o_h_plane(mod):
    zn, o01, h03 = (0.0, 0.0, 0.0), (1.9, 0.0, 0.0), (2.4, 0.83, 0.3)
    h04 = mod.place_h04(zn, o01, h03)
    normal = mod.unit(mod.cross(mod.sub(h03, o01), mod.sub(zn, o01)))
    assert mod.dot(normal, mod.sub(h04, o01)) == pytest.approx(0.0, abs=1e-6)


def test_h04_points_away_from_the_metal(mod):
    """Both hydrogens of a metal-bound water point away from the cation. The
    rejected alternative solution folds H04 back towards the zinc."""
    zn, o01, h03 = (0.0, 0.0, 0.0), (1.9, 0.0, 0.0), (2.4, 0.83, 0.0)
    h04 = mod.place_h04(zn, o01, h03)
    assert mod.norm(mod.sub(h04, zn)) > mod.norm(mod.sub(o01, zn))


def test_collinear_geometry_is_refused(mod):
    """With H03 along the Zn-O axis there is no plane and no 'away' side."""
    with pytest.raises(ValueError, match="collinear"):
        mod.place_h04((0.0, 0.0, 0.0), (1.9, 0.0, 0.0), (2.9, 0.0, 0.0))


# ---- the file it writes -----------------------------------------------------

def test_writes_h04_and_renames_the_residue(mod, lysH, tmp_path, capsys):
    out = tmp_path / "out.pdb"
    mod.main([str(lysH), str(out)])
    text = out.read_text()
    assert "H04" in text
    assert "NUW" in text and " OH " not in text
    atoms = [l for l in text.splitlines() if l.startswith(("ATOM", "HETATM"))]
    assert len(atoms) == 5, "one hydrogen added to the four input atoms"
    assert "H04" in capsys.readouterr().out


def test_serials_are_renumbered_without_duplicates(mod, lysH, tmp_path):
    out = tmp_path / "out.pdb"
    mod.main([str(lysH), str(out)])
    serials = [int(l[6:11]) for l in out.read_text().splitlines()
               if l.startswith(("ATOM", "HETATM"))]
    assert serials == list(range(1, len(serials) + 1))


def test_the_new_bond_is_written_as_conect_before_end(mod, lysH, tmp_path):
    """PDBFile stops at END, so a CONECT after it is dropped in silence and the
    new hydrogen arrives unbonded - which the force field then rejects."""
    out = tmp_path / "out.pdb"
    mod.main([str(lysH), str(out)])
    lines = out.read_text().splitlines()
    conect = [i for i, l in enumerate(lines) if l.startswith("CONECT")]
    end = next(i for i, l in enumerate(lines) if l.startswith("END"))
    assert conect and max(conect) < end, "a CONECT record sits after END"
    o01 = next(int(l[6:11]) for l in lines
               if l.startswith(("ATOM", "HETATM")) and l[12:16].strip() == "O01")
    h04 = next(int(l[6:11]) for l in lines
               if l.startswith(("ATOM", "HETATM")) and l[12:16].strip() == "H04")
    bonded = set()
    for l in lines:
        if not l.startswith("CONECT"):
            continue
        fields = [l[p:p + 5].strip() for p in range(6, 31, 5)]
        fields = [int(f) for f in fields if f]
        for other in fields[1:]:
            bonded.add(frozenset((fields[0], other)))
    assert frozenset((o01, h04)) in bonded, "O01-H04 was not recorded"


def test_input_without_conect_is_refused(mod, tmp_path):
    lines = [
        pdb_line(1, "ZN", "ZN", "B", "1", (0.0, 0.0, 0.0)),
        pdb_line(2, "O01", "OH", "C", "2", (1.9, 0.0, 0.0)),
        pdb_line(3, "H03", "OH", "C", "2", (2.4, 0.83, 0.0)),
        "END\n",
    ]
    src = tmp_path / "noconect.pdb"
    src.write_text("".join(lines))
    out = tmp_path / "out.pdb"
    with pytest.raises(SystemExit, match="CONECT"):
        mod.main([str(src), str(out)])
    assert not out.exists()


def test_a_clash_is_refused_and_nothing_is_written(mod, tmp_path):
    """A spectator atom sitting where H04 belongs must stop the run. The
    placement looks fine in a PDB and only fails later, in a minimisation."""
    lines = [
        pdb_line(1, "ZN", "ZN", "B", "1", (0.0, 0.0, 0.0)),
        pdb_line(2, "O01", "OH", "C", "2", (1.9, 0.0, 0.0)),
        pdb_line(3, "H03", "OH", "C", "2", (2.4, 0.83, 0.0)),
        pdb_line(4, "NE2", "HIS", "A", "88", (2.4, -0.83, 0.0)),
        "END\n",
    ]
    src = tmp_path / "clash.pdb"
    src.write_text("".join(lines))
    out = tmp_path / "out.pdb"
    with pytest.raises(SystemExit, match="clash cutoff"):
        mod.main([str(src), str(out)])
    assert not out.exists()


def test_missing_atoms_are_named_in_the_refusal(mod, tmp_path):
    src = tmp_path / "bare.pdb"
    src.write_text(pdb_line(1, "CA", "GLY", "A", "9", (0.0, 0.0, 0.0)))
    out = tmp_path / "out.pdb"
    with pytest.raises(SystemExit, match="Zn.*O01.*H03|could not find"):
        mod.main([str(src), str(out)])
    assert not out.exists()


def test_running_twice_is_refused(mod, lysH, tmp_path):
    out = tmp_path / "out.pdb"
    mod.main([str(lysH), str(out)])
    again = tmp_path / "again.pdb"
    with pytest.raises(SystemExit, match="already"):
        mod.main([str(out), str(again)])
    assert not again.exists()
