"""Coordination over a trajectory, in mcube: pytest tests/test_coordination_occupancy.py

`md` reports the sphere before and after a run - two frames. Occupancy over the
whole trajectory lived only in `docs/csr/site.py`, so the numbers in D134's
table came from outside the package and outside the config (D167). What these
tests pin is that the in-framework version measures by the SAME rule that built
the restraints, and keeps the two distinctions a weaker version loses: a
bidentate carboxylate counted as two ligands rather than one residue, and the
coordination number as a distribution rather than a mean.

Needs MDAnalysis, which comes with the ASH environment.
"""

import sys
import types

import pytest

pytest.importorskip("MDAnalysis", reason="postprocess reads trajectories through MDAnalysis")
pytest.importorskip("openmm", reason="preprocess imports OpenMM at module scope")

from tests.test_qmmm_structure import ASH_NAMES


@pytest.fixture
def postprocess(monkeypatch):
    stub = types.ModuleType("ash")
    for name in ASH_NAMES:
        setattr(stub, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", stub)
    for name in [m for m in sys.modules if m == "mcube" or m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    from mcube import postprocess
    return postprocess


def atom(serial, name, resname, chain, seq, xyz, element):
    return ("HETATM%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
            % (serial, name, resname, chain, seq, xyz[0], xyz[1], xyz[2], element))


ZN = (0.0, 0.0, 0.0)
INPUTS = {"metals": ["ZN"], "cutoff": 2.6, "donor_elements": ["N", "O", "S"]}


@pytest.fixture
def system(tmp_path):
    """A zinc, a histidine NE2, and a carboxylate whose two oxygens both reach.

    OD1 starts outside the cutoff and can be moved in, which is what a
    bidentate carboxylate looks like forming.
    """
    pdb = tmp_path / "top.pdb"
    pdb.write_text(
        "CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1           1\n"
        + atom(1, "ZN", "ZN", "B", 1, ZN, "ZN")
        + atom(2, "NE2", "HIS", "A", 88, (2.10, 0.0, 0.0), "N")
        + atom(3, "OD2", "ASP", "A", 92, (0.0, 2.05, 0.0), "O")
        + atom(4, "OD1", "ASP", "A", 92, (0.0, 0.0, 3.00), "O")
        + "END\n")
    return pdb


def trajectory(tmp_path, system, od1_z):
    """A multi-frame PDB trajectory moving OD1 to each given z."""
    lines = system.read_text().splitlines(keepends=True)
    out = tmp_path / "traj.pdb"
    with out.open("w") as fh:
        for n, z in enumerate(od1_z, start=1):
            fh.write("MODEL %8d\n" % n)
            for line in lines:
                if line.startswith("HETATM") and line[12:16].strip() == "OD1":
                    line = line[:46] + "%8.3f" % z + line[54:]
                if line.startswith(("ATOM", "HETATM")):
                    fh.write(line)
            fh.write("ENDMDL\n")
    return out


class TestTheRuleIsTheRestraintRule:
    def test_the_cutoff_comes_from_preprocess_inputs(self, postprocess, tmp_path, system):
        """Not from a CLI flag. `--cutoff 2.6` never entered the config, so the
        config could not say what had been measured."""
        traj = trajectory(tmp_path, system, [3.0, 3.0])
        rep = postprocess.coordination_occupancy(str(system), str(traj), INPUTS)
        assert rep["cutoff"] == 2.6

    def test_an_explicit_radius_overrides_it(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [3.0, 3.0])
        rep = postprocess.coordination_occupancy(str(system), str(traj), INPUTS,
                                                 radius=2.0)
        assert rep["cutoff"] == 2.0

    def test_a_donor_inside_the_cutoff_is_fully_occupied(self, postprocess, tmp_path,
                                                         system):
        traj = trajectory(tmp_path, system, [3.0, 3.0])
        rep = postprocess.coordination_occupancy(str(system), str(traj), INPUTS)
        his = [r for r in rep["ligands"] if r["resname"] == "HIS"][0]
        assert his["occupancy"] == 100.0
        assert his["min"] == his["max"] == 2.10

    def test_frames_are_counted(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [3.0, 3.0, 3.0])
        rep = postprocess.coordination_occupancy(str(system), str(traj), INPUTS)
        assert rep["frames"] == 3

    def test_a_stride_subsamples(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [3.0] * 6)
        rep = postprocess.coordination_occupancy(str(system), str(traj), INPUTS,
                                                 stride=3)
        assert rep["frames"] == 2


class TestTheBidentateDistinction:
    """D166 measured Asp92 monodentate at 2.02 in one run and bidentate
    (2.06/2.12) in another. Keyed by residue those are the same edge."""

    def test_both_oxygens_are_reported_separately(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [2.10, 2.10])
        rep = postprocess.coordination_occupancy(str(system), str(traj), INPUTS)
        asp = [r["donor"] for r in rep["ligands"] if r["resname"] == "ASP"]
        assert len(asp) == 2
        assert any(d.endswith("OD1") for d in asp)
        assert any(d.endswith("OD2") for d in asp)

    def test_per_residue_collapses_them(self, postprocess, tmp_path, system):
        """Available, but not the default."""
        traj = trajectory(tmp_path, system, [2.10, 2.10])
        rep = postprocess.coordination_occupancy(str(system), str(traj), INPUTS,
                                                 per_residue=True)
        asp = [r for r in rep["ligands"] if r["resname"] == "ASP"]
        assert len(asp) == 1

    def test_a_partial_occupancy_is_a_fraction_not_a_presence(self, postprocess,
                                                              tmp_path, system):
        """OD1 inside for one of four frames."""
        traj = trajectory(tmp_path, system, [2.10, 3.00, 3.00, 3.00])
        rep = postprocess.coordination_occupancy(str(system), str(traj), INPUTS)
        od1 = [r for r in rep["ligands"] if r["donor"].endswith("OD1")][0]
        assert od1["occupancy"] == 25.0
        assert od1["min"] == 2.10 and od1["max"] == 3.00


class TestTheCoordinationNumber:
    def test_it_is_a_distribution_not_a_mean(self, postprocess, tmp_path, system):
        """A site flipping between 2 and 3 is not a site sitting at 2.5."""
        traj = trajectory(tmp_path, system, [2.10, 3.00, 2.10, 3.00])
        rep = postprocess.coordination_occupancy(str(system), str(traj), INPUTS)
        assert rep["coordination_number"] == {"2": 2, "3": 2}

    def test_a_steady_site_reports_one_number(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [3.00, 3.00, 3.00])
        rep = postprocess.coordination_occupancy(str(system), str(traj), INPUTS)
        assert rep["coordination_number"] == {"2": 3}


class TestTheMargin:
    def test_a_donor_just_outside_is_listed_at_zero(self, postprocess, tmp_path,
                                                    system):
        """2.7 against a 2.6 cutoff is a borderline call about the site, not an
        absence, and hiding it makes the cutoff look like a measurement."""
        traj = trajectory(tmp_path, system, [2.70, 2.70])
        rep = postprocess.coordination_occupancy(str(system), str(traj), INPUTS)
        od1 = [r for r in rep["ligands"] if r["donor"].endswith("OD1")][0]
        assert od1["occupancy"] == 0.0 and od1["mean"] == 2.70

    def test_far_donors_are_not_listed(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [9.00, 9.00])
        rep = postprocess.coordination_occupancy(str(system), str(traj), INPUTS)
        assert not [r for r in rep["ligands"] if r["donor"].endswith("OD1")]


class TestTheTopologyMustMatch:
    def test_a_mismatched_atom_count_is_refused(self, postprocess, tmp_path, system):
        """Silently zipping the wrong topology onto a trajectory would report
        one atom's distances under another atom's name."""
        traj = trajectory(tmp_path, system, [2.10])
        short = tmp_path / "short.pdb"
        short.write_text(atom(1, "ZN", "ZN", "B", 1, ZN, "ZN") + "END\n")
        with pytest.raises(ValueError) as exc:
            postprocess.coordination_occupancy(str(short), str(traj), INPUTS)
        assert "atoms" in str(exc.value)
