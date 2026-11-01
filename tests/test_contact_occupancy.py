"""Monitored contacts over a trajectory: pytest tests/test_contact_occupancy.py

`coordination_occupancy` answers what is on the metal; this answers everything
else, for the pairs the config names. It was `docs/csr/site.py contacts` - a
script, so a result from it was not recorded in any config (D168).

Two behaviours carry history and are pinned here. A pair that does not resolve
is REPORTED, not skipped: D120 is the entry where the oxygen that migrated onto
the metal was in no monitor at all, and a selector that silently resolves to
nothing is the same failure with extra steps. And the selectors go through
`resolve_atoms`, the same parser restraints use, so a config means one thing
everywhere.

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


@pytest.fixture
def system(tmp_path):
    """A zinc, a carboxylate oxygen, and a far-off amide oxygen."""
    pdb = tmp_path / "top.pdb"
    pdb.write_text(
        "CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1           1\n"
        + atom(1, "ZN", "ZN", "B", 1, (0.0, 0.0, 0.0), "ZN")
        + atom(2, "OD2", "ASP", "A", 92, (2.00, 0.0, 0.0), "O")
        + atom(3, "OE1", "GLN", "A", 93, (0.0, 8.00, 0.0), "O")
        + "END\n")
    return pdb


def trajectory(tmp_path, system, oe1_y):
    """Frames moving the Gln oxygen to each given y."""
    lines = system.read_text().splitlines(keepends=True)
    out = tmp_path / "traj.pdb"
    with out.open("w") as fh:
        for n, y in enumerate(oe1_y, start=1):
            fh.write("MODEL %8d\n" % n)
            for line in lines:
                if line.startswith("HETATM") and line[12:16].strip() == "OE1":
                    line = line[:38] + "%8.3f" % y + line[46:]
                if line.startswith(("ATOM", "HETATM")):
                    fh.write(line)
            fh.write("ENDMDL\n")
    return out


MONITOR = {"Zn-Asp92 OD2": ["B:1:ZN", "A:92:OD2"],
           "Zn-Gln93 OE1": ["B:1:ZN", "A:93:OE1"]}


class TestOccupancy:
    def test_a_held_contact_is_fully_occupied(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [8.0, 8.0])
        rep = postprocess.contact_occupancy(str(system), str(traj), MONITOR)
        asp = [r for r in rep["pairs"] if r["pair"] == "Zn-Asp92 OD2"][0]
        assert asp["occupancy"] == 100.0 and asp["mean"] == 2.00

    def test_a_far_contact_is_never_occupied(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [8.0, 8.0])
        rep = postprocess.contact_occupancy(str(system), str(traj), MONITOR)
        gln = [r for r in rep["pairs"] if r["pair"] == "Zn-Gln93 OE1"][0]
        assert gln["occupancy"] == 0.0 and gln["min"] == 8.00

    def test_occupancy_is_a_fraction_of_frames(self, postprocess, tmp_path, system):
        """Half the frames inside 3.5 A. A contact that is there half the time
        is not a contact that is there."""
        traj = trajectory(tmp_path, system, [3.0, 3.0, 8.0, 8.0])
        rep = postprocess.contact_occupancy(str(system), str(traj), MONITOR)
        gln = [r for r in rep["pairs"] if r["pair"] == "Zn-Gln93 OE1"][0]
        assert gln["occupancy"] == 50.0
        assert gln["min"] == 3.00 and gln["max"] == 8.00

    def test_the_cutoff_is_configurable(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [3.0, 3.0])
        rep = postprocess.contact_occupancy(str(system), str(traj), MONITOR,
                                            cutoff=2.5)
        gln = [r for r in rep["pairs"] if r["pair"] == "Zn-Gln93 OE1"][0]
        assert gln["occupancy"] == 0.0

    def test_pairs_are_sorted_by_occupancy(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [8.0, 8.0])
        rep = postprocess.contact_occupancy(str(system), str(traj), MONITOR)
        assert [r["pair"] for r in rep["pairs"]][0] == "Zn-Asp92 OD2"

    def test_a_stride_subsamples(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [8.0] * 6)
        rep = postprocess.contact_occupancy(str(system), str(traj), MONITOR,
                                            stride=3)
        assert rep["frames"] == 2


class TestAPairThatDoesNotResolve:
    """Reported, not skipped. A monitor that silently matches nothing is how
    the atom that mattered went unwatched for weeks (D120)."""

    def test_a_missing_atom_is_listed_with_a_reason(self, postprocess, tmp_path,
                                                    system):
        traj = trajectory(tmp_path, system, [8.0])
        rep = postprocess.contact_occupancy(
            str(system), str(traj), {"Zn-nothing": ["B:1:ZN", "A:999:OD1"]})
        assert [r["pair"] for r in rep["missing"]] == ["Zn-nothing"]
        assert "A:999:OD1" in rep["missing"][0]["why"]

    def test_it_does_not_stop_the_others(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [8.0])
        monitor = dict(MONITOR, **{"Zn-nothing": ["B:1:ZN", "A:999:OD1"]})
        rep = postprocess.contact_occupancy(str(system), str(traj), monitor)
        assert len(rep["pairs"]) == 2 and len(rep["missing"]) == 1

    def test_a_malformed_entry_is_also_reported(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [8.0])
        rep = postprocess.contact_occupancy(
            str(system), str(traj), {"lonely": ["B:1:ZN"]})
        assert rep["missing"][0]["pair"] == "lonely"
        assert "two selectors" in rep["missing"][0]["why"]


class TestTheProfile:
    def test_it_reports_a_minimum_and_the_counts(self, postprocess, tmp_path,
                                                  system):
        """-RT ln P over the range the trajectory visited. The per-bin count is
        part of the answer: a bin with a handful of frames is noise."""
        traj = trajectory(tmp_path, system, [3.0] * 8 + [4.0] * 2)
        rep = postprocess.contact_occupancy(str(system), str(traj), MONITOR)
        prof = postprocess.contact_profile(rep, "Gln93", width=0.5)
        assert len(prof) == 1
        assert prof[0]["minimum"] == pytest.approx(3.25, abs=0.3)
        assert prof[0]["sampled_frames"] == 10
        assert sum(b["frames"] for b in prof[0]["bins"]) == 10

    def test_the_minimum_has_zero_free_energy(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [3.0] * 8 + [4.0] * 2)
        rep = postprocess.contact_occupancy(str(system), str(traj), MONITOR)
        prof = postprocess.contact_profile(rep, "Gln93", width=0.5)
        assert min(b["dG"] for b in prof[0]["bins"]) == 0.0

    def test_it_selects_by_substring(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [3.0, 4.0])
        rep = postprocess.contact_occupancy(str(system), str(traj), MONITOR)
        assert postprocess.contact_profile(rep, "nothing-matches") == []

    def test_temperature_scales_the_energies(self, postprocess, tmp_path, system):
        traj = trajectory(tmp_path, system, [3.0] * 8 + [4.0] * 2)
        rep = postprocess.contact_occupancy(str(system), str(traj), MONITOR)
        cold = postprocess.contact_profile(rep, "Gln93", width=0.5, temperature=150.0)
        warm = postprocess.contact_profile(rep, "Gln93", width=0.5, temperature=300.0)
        assert max(b["dG"] for b in warm[0]["bins"]) > \
            max(b["dG"] for b in cold[0]["bins"])


class TestTheDistancesAreNotKeptInTheResult:
    def test_the_report_carries_them_for_the_profile(self, postprocess, tmp_path,
                                                     system):
        """Needed to build a profile, but megabytes per pair - `run` drops them
        before saving, since the config is read by hand."""
        traj = trajectory(tmp_path, system, [8.0, 8.0])
        rep = postprocess.contact_occupancy(str(system), str(traj), MONITOR)
        assert len(rep["_distances"]["Zn-Asp92 OD2"]) == 2
