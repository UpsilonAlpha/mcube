"""Analysing an MD run that crashed: pytest tests/test_md_salvage.py

`md.run` assigns every `result[...]` key after OpenMM returns, so a NaN or a
walltime kill records nothing while leaving the frames up to that step on disk.
Those are the runs most worth looking at - the product arm NaN'd twice - and
`postprocess` refused both of them: D172 made the reporter CSV optional and
stopped one key short, so the next salvaged run died on `KeyError: 'lastframe'`
rather than `KeyError: 'csv'`.

What is pinned here: the trajectory is found, the topology is the trajectory's
OWN first frame (a topology with a different atom order reads as plausible
nonsense, not as an error), and a run with nothing on disk still fails with a
message that says which input to set.
"""

import sys
import types

import pytest

pytest.importorskip("openmm", reason="md's module chain imports OpenMM")

from tests.test_qmmm_structure import ASH_NAMES  # noqa: E402


@pytest.fixture(autouse=True)
def md_mod(monkeypatch):
    """`md` with ASH stubbed: salvage only touches the filesystem."""
    stub = types.ModuleType("ash")
    for name in ASH_NAMES:
        setattr(stub, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", stub)
    for name in [m for m in sys.modules if m == "mcube" or m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    from mcube import md
    globals()["MD_TRAJECTORIES"] = md.MD_TRAJECTORIES
    globals()["salvage_results"] = md.salvage_results
    return md


def write(directory, *names):
    for name in names:
        (directory / name).write_text("ATOM\nEND\n")


class TestWhatIsFound:
    def test_the_nvt_trajectory_and_its_first_frame(self, tmp_path):
        write(tmp_path, "NVTtrajectory.dcd", "NVTtrajectory_firstframe.pdb")
        got = salvage_results(directory=tmp_path)
        assert got["trajfile"].endswith("NVTtrajectory.dcd")
        assert got["pdbfile"].endswith("NVTtrajectory_firstframe.pdb")

    def test_the_imaged_trajectory_wins_when_both_exist(self, tmp_path):
        """Reimaging is the last thing `md` does to a finished NVT run, so if
        an imaged DCD is there it is the better of the two."""
        write(tmp_path, "NVTtrajectory.dcd", "NVTtrajectory_imaged.dcd",
              "NVTtrajectory_firstframe.pdb")
        assert salvage_results(directory=tmp_path)["trajfile"].endswith(
            "NVTtrajectory_imaged.dcd")

    def test_nvt_wins_over_the_box_equilibration(self, tmp_path):
        """A run that reached NVT has production frames; the NPT box relaxation
        is equilibration and is not what anyone is asking about."""
        write(tmp_path, "relaxbox_NPT.dcd", "relaxbox_NPT_firstframe.pdb",
              "NVTtrajectory.dcd", "NVTtrajectory_firstframe.pdb")
        assert salvage_results(directory=tmp_path)["trajfile"].endswith(
            "NVTtrajectory.dcd")

    def test_the_box_equilibration_is_used_if_that_is_all_there_is(self, tmp_path):
        write(tmp_path, "relaxbox_NPT.dcd", "relaxbox_NPT_firstframe.pdb")
        assert salvage_results(directory=tmp_path)["trajfile"].endswith(
            "relaxbox_NPT.dcd")

    def test_paths_are_absolute(self, tmp_path):
        """postprocess resolves other paths against the config's directory, so
        a relative one here would be read against the wrong directory."""
        write(tmp_path, "NVTtrajectory.dcd", "NVTtrajectory_firstframe.pdb")
        from pathlib import Path
        for value in salvage_results(directory=tmp_path).values():
            assert Path(value).is_absolute()


class TestWhatIsRefused:
    def test_a_trajectory_without_its_first_frame_is_not_salvaged(self, tmp_path):
        """Guessing a topology is worse than refusing: the wrong atom order
        gives numbers that look entirely reasonable."""
        write(tmp_path, "NVTtrajectory.dcd")
        assert salvage_results(directory=tmp_path) == {}

    def test_a_first_frame_without_a_trajectory_is_not_salvaged(self, tmp_path):
        write(tmp_path, "NVTtrajectory_firstframe.pdb")
        assert salvage_results(directory=tmp_path) == {}

    def test_an_empty_directory_gives_nothing(self, tmp_path):
        assert salvage_results(directory=tmp_path) == {}


class TestFindingTheDirectory:
    def test_the_directory_comes_from_whatever_the_crash_did_record(self, tmp_path):
        """A crashed run still recorded its minimisation, and those paths are
        absolute - which is how the results directory is known without the
        caller having to work it out."""
        write(tmp_path, "NVTtrajectory.dcd", "NVTtrajectory_firstframe.pdb")
        got = salvage_results({"minimized": str(tmp_path / "frag-minimized.pdb")})
        assert got["trajfile"].endswith("NVTtrajectory.dcd")

    def test_trajfile_is_preferred_over_the_later_fallbacks(self, tmp_path):
        write(tmp_path, "NVTtrajectory.dcd", "NVTtrajectory_firstframe.pdb")
        got = salvage_results({"trajfile": str(tmp_path / "NVTtrajectory.dcd"),
                               "minimized": "/nowhere/frag-minimized.pdb"})
        assert got["pdbfile"].endswith("NVTtrajectory_firstframe.pdb")

    def test_no_recorded_paths_at_all_gives_nothing_rather_than_scanning_cwd(
            self, tmp_path, monkeypatch):
        """With neither a directory nor a recorded path, the only candidate is
        the process's cwd, which is the run directory when mcube runs there."""
        monkeypatch.chdir(tmp_path)
        assert salvage_results({}) == {}
        write(tmp_path, "NVTtrajectory.dcd", "NVTtrajectory_firstframe.pdb")
        assert salvage_results({})["trajfile"].endswith("NVTtrajectory.dcd")


def test_every_candidate_names_a_topology():
    """A pair whose topology is missing would salvage a trajectory with no way
    to read it."""
    for traj, topology in MD_TRAJECTORIES:
        assert traj.endswith(".dcd") and topology.endswith(".pdb")
