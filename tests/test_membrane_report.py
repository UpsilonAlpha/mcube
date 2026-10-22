"""
Tests for scripts/membrane_report.py on OpenMM's pre-equilibrated POPC patch (no ASH):

    pytest tests/test_membrane_report.py
"""

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("openmm", reason="the membrane report reads OpenMM's POPC patch")
mdtraj = pytest.importorskip("mdtraj", reason="the membrane report uses mdtraj")

from openmm.app import forcefield as openmm_forcefield  # noqa: E402

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "membrane_report.py"


@pytest.fixture(scope="module")
def patch_files(tmp_path_factory):
    """The POPC patch as a PDB plus a two-frame DCD (the same frame twice)."""
    data = Path(openmm_forcefield.__file__).parent / "data" / "POPC.pdb"
    out = tmp_path_factory.mktemp("patch")
    traj = mdtraj.load(str(data))
    traj = mdtraj.join([traj, traj])
    traj[0].save_pdb(str(out / "top.pdb"))
    traj.save_dcd(str(out / "traj.dcd"))
    return out


def run(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], check=True,
                          capture_output=True, text=True)


def test_membrane_only_report(patch_files):
    run(patch_files / "top.pdb", patch_files / "traj.dcd", "--out", patch_files / "r", "--frame-ps", 10)
    rows = (patch_files / "r.tsv").read_text().splitlines()
    header, last = rows[0].split("\t"), dict(zip(rows[0].split("\t"), rows[-1].split("\t")))
    assert len(rows) == 3  # header + 2 frames
    assert "bb_rmsd_core" not in header  # no protein, so no receptor columns
    assert 60 < float(last["apl_gross"]) < 70  # POPC patch: 66.6 A^2 per lipid
    assert 35 < float(last["thickness_PP"]) < 40
    assert int(float(last["leaflet_imbalance"])) == 0
    assert (patch_files / "r.txt").read_text().startswith("frames=2")


def test_script_loads_without_side_effects():
    spec = importlib.util.spec_from_file_location("membrane_report", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.residue_range("290-330") == (290, 330)
    assert module.residue_range("12") == (12, 12)
    assert module.named_residues("gate=319,322", 2) == ("gate", [319, 322])
    with pytest.raises(Exception):
        module.named_residues("gate=319", 2)
