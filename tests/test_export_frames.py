"""Writing a trajectory frame back out: pytest tests/test_export_frames.py

A crashed run has to be restarted from somewhere, and the only thing a stage
can be pointed at is a file. The product arm's MD died at 34.14 ns with
`Particle coordinate is NaN`, leaving a good trajectory and no `lastframe`, so
the restart structure had to come out of the DCD - and the frames nearest the
crash are the ones whose geometry is least trustworthy, which is why negative
indices matter (D172).

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


def atom(serial, name, resname, chain, seq, xyz):
    return ("ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
            % (serial, name, resname, chain, seq, xyz[0], xyz[1], xyz[2], name[:1]))


@pytest.fixture
def system(tmp_path):
    pdb = tmp_path / "top.pdb"
    pdb.write_text(
        "CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1           1\n"
        + atom(1, "N", "ALA", "A", 1, (0.0, 0.0, 0.0))
        + atom(2, "CA", "ALA", "A", 1, (1.5, 0.0, 0.0))
        + "END\n")
    return pdb


@pytest.fixture
def traj(tmp_path, system):
    """Five frames, the CA walking along x so each frame is identifiable."""
    lines = system.read_text().splitlines(keepends=True)
    out = tmp_path / "traj.pdb"
    with out.open("w") as fh:
        for n, x in enumerate([1.5, 2.5, 3.5, 4.5, 5.5], start=1):
            fh.write("MODEL %8d\n" % n)
            for line in lines:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    line = line[:30] + "%8.3f" % x + line[38:]
                if line.startswith("ATOM"):
                    fh.write(line)
            fh.write("ENDMDL\n")
    return out


def ca_x(path):
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
            return round(float(line[30:38]), 1)
    raise AssertionError("no CA in " + str(path))


class TestWhichFrame:
    def test_a_positive_index_is_that_frame(self, postprocess, tmp_path, system,
                                            traj, monkeypatch):
        monkeypatch.chdir(tmp_path)
        got = postprocess.export_frames(str(system), str(traj), [2])
        assert len(got) == 1 and ca_x(got[0]) == 3.5

    def test_minus_one_is_the_last_frame(self, postprocess, tmp_path, system,
                                        traj, monkeypatch):
        monkeypatch.chdir(tmp_path)
        got = postprocess.export_frames(str(system), str(traj), [-1])
        assert ca_x(got[0]) == 5.5

    def test_a_couple_of_frames_back(self, postprocess, tmp_path, system, traj,
                                     monkeypatch):
        """What a restart after a NaN asks for: not the last frame, which is the
        one nearest whatever went wrong."""
        monkeypatch.chdir(tmp_path)
        got = postprocess.export_frames(str(system), str(traj), [-3])
        assert ca_x(got[0]) == 3.5

    def test_several_frames_at_once(self, postprocess, tmp_path, system, traj,
                                    monkeypatch):
        monkeypatch.chdir(tmp_path)
        got = postprocess.export_frames(str(system), str(traj), [0, -1])
        assert [ca_x(g) for g in got] == [1.5, 5.5]

    def test_the_index_is_in_the_filename(self, postprocess, tmp_path, system,
                                          traj, monkeypatch):
        """Resolved, not as given: `frame_-3.pdb` would say nothing about which
        geometry a later run was started from."""
        monkeypatch.chdir(tmp_path)
        got = postprocess.export_frames(str(system), str(traj), [-3])
        assert got[0].endswith("frame_2.pdb")

    def test_the_prefix_is_configurable(self, postprocess, tmp_path, system,
                                        traj, monkeypatch):
        monkeypatch.chdir(tmp_path)
        got = postprocess.export_frames(str(system), str(traj), [0], out="restart")
        assert got[0].endswith("restart_0.pdb")


class TestRefusals:
    def test_a_frame_past_the_end_is_refused(self, postprocess, tmp_path, system,
                                             traj, monkeypatch):
        """Silently clamping to the last frame would restart a run from a
        geometry nobody chose."""
        monkeypatch.chdir(tmp_path)
        with pytest.raises(ValueError) as exc:
            postprocess.export_frames(str(system), str(traj), [99])
        assert "5 frames" in str(exc.value)

    def test_a_negative_index_past_the_start_is_refused(self, postprocess, tmp_path,
                                                        system, traj, monkeypatch):
        monkeypatch.chdir(tmp_path)
        with pytest.raises(ValueError):
            postprocess.export_frames(str(system), str(traj), [-99])


def test_the_written_frame_keeps_the_atom_names(postprocess, tmp_path, system,
                                                traj, monkeypatch):
    """The point of writing through the topology: the result is the same system
    in another geometry, which is what `md` resumes from and what `neb` demands
    of its endpoints."""
    monkeypatch.chdir(tmp_path)
    got = postprocess.export_frames(str(system), str(traj), [1])
    names = [line[12:16].strip() for line in open(got[0])
             if line.startswith(("ATOM", "HETATM"))]
    assert names == ["N", "CA"]
