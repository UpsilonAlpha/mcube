"""Classical MD resume: pytest tests/test_md_resume.py

`qmmm` has had `resume: true` since D109. `md` had not, so a 100 ns
equilibration that hit its walltime restarted from zero - and the files it
needed were there the whole time, because ASH writes `OpenMM_MD_state.xml` as
it runs whichever stage calls `OpenMM_MD` (D182).

Two things are specific to this stage rather than inherited from `qmmm`. The
time and timestep live under `nvt_time` and `nvt_timestep`, not `time` and
`timestep`, so the same helper has to be told which keys to read - passing the
wrong one silently resumes against `None` and refuses, or worse, runs the full
time again. And `md` has phases *before* the production MD: minimisation and
NPT. A resumed run must skip them, because the state file carries positions and
velocities, so re-equilibrating the box spends the walltime that resuming
exists to save and then throws the coordinates away.

The state files here are written by OpenMM itself, so the header being parsed
is the real format.
"""

import sys
import types

import pytest

pytest.importorskip("openmm", reason="resume reads OpenMM state files")

from tests.test_qmmm_structure import ASH_NAMES  # noqa: E402
from tests.test_qmmm_resume import write_state  # noqa: E402


@pytest.fixture
def md(monkeypatch):
    stub = types.ModuleType("ash")
    for name in ASH_NAMES:
        setattr(stub, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", stub)
    for name in [m for m in sys.modules if m == "mcube" or m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    from mcube import md
    return md


class TestTheKeysAreThisStagesKeys:
    def test_the_remainder_comes_from_nvt_time(self, md, tmp_path, monkeypatch):
        """100 ns asked for, 34 ns done, so 66 ns left - the arithmetic that
        the product arm's NaN crash would have needed."""
        monkeypatch.chdir(tmp_path)
        write_state(tmp_path / "OpenMM_MD_state.xml", 1000, timestep_ps=0.004)
        state, remaining = md.resolve_resume(
            {"resume": True, "nvt_time": 100.0, "nvt_timestep": 0.004},
            time_key="nvt_time", timestep_key="nvt_timestep")
        assert state is not None
        assert remaining == pytest.approx(96.0, abs=1e-6)

    def test_the_qmmm_keys_are_still_the_default(self, md, tmp_path, monkeypatch):
        """`qmmm` reaches the same helper through a re-export and passes no
        keys, so the defaults must stay what that stage uses."""
        monkeypatch.chdir(tmp_path)
        write_state(tmp_path / "OpenMM_MD_state.xml", 1000, timestep_ps=0.001)
        state, remaining = md.resolve_resume({"resume": True, "time": 50.0})
        assert state is not None and remaining == pytest.approx(49.0, abs=1e-6)

    def test_reading_the_wrong_key_is_refused_not_guessed(self, md, tmp_path,
                                                          monkeypatch):
        """A config with `nvt_time` read as `time` would find None. Refusing
        names the key, so the message is actionable."""
        monkeypatch.chdir(tmp_path)
        write_state(tmp_path / "OpenMM_MD_state.xml", 1000)
        with pytest.raises(ValueError) as exc:
            md.resolve_resume({"resume": True, "nvt_time": 100.0},
                              time_key="time")
        assert "`time` is not set" in str(exc.value)

    def test_the_timestep_key_decides_what_counts_as_finished(self, md, tmp_path,
                                                              monkeypatch):
        """4 fs, not 1 fs: a run inside one timestep of its target is done, and
        the threshold has to be this stage's own step."""
        monkeypatch.chdir(tmp_path)
        write_state(tmp_path / "OpenMM_MD_state.xml", 1000, timestep_ps=0.004)
        _, remaining = md.resolve_resume(
            {"resume": True, "nvt_time": 4.001, "nvt_timestep": 0.004},
            time_key="nvt_time", timestep_key="nvt_timestep")
        assert remaining == 0.0


class TestWithoutResume:
    def test_the_requested_time_passes_through(self, md, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        write_state(tmp_path / "OpenMM_MD_state.xml", 1000, timestep_ps=0.004)
        state, remaining = md.resolve_resume(
            {"nvt_time": 100.0}, time_key="nvt_time", timestep_key="nvt_timestep")
        assert state is None and remaining == 100.0

    def test_resume_with_no_state_file_starts_fresh(self, md, tmp_path,
                                                    monkeypatch, capsys):
        monkeypatch.chdir(tmp_path)
        state, remaining = md.resolve_resume(
            {"resume": True, "nvt_time": 100.0, "nvt_timestep": 0.004},
            time_key="nvt_time", timestep_key="nvt_timestep")
        assert state is None and remaining == 100.0
        assert "nothing to resume from" in capsys.readouterr().out


class TestTheTrajectorySegments:
    def test_the_previous_trajectory_is_moved_aside(self, md, tmp_path, monkeypatch):
        """ASH builds its DCDReporter with append=False, so a restart truncates
        the trajectory it resumed from unless that file is moved first."""
        monkeypatch.chdir(tmp_path)
        (tmp_path / "NVTtrajectory.dcd").write_bytes(b"frames")
        kept = md.preserve_trajectory("NVTtrajectory.dcd")
        assert kept.name == "NVTtrajectory_part1.dcd"
        assert not (tmp_path / "NVTtrajectory.dcd").exists()
        assert kept.read_bytes() == b"frames"

    def test_segments_accumulate(self, md, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        for expected in ("NVTtrajectory_part1.dcd", "NVTtrajectory_part2.dcd"):
            (tmp_path / "NVTtrajectory.dcd").write_bytes(b"x")
            assert md.preserve_trajectory("NVTtrajectory.dcd").name == expected

    def test_the_imaged_copy_is_kept_too(self, md, tmp_path, monkeypatch):
        """`md` writes both a raw and an imaged DCD, and a restart would
        overwrite either one it is pointed at."""
        monkeypatch.chdir(tmp_path)
        (tmp_path / "NVTtrajectory_imaged.dcd").write_bytes(b"y")
        kept = md.preserve_trajectory("NVTtrajectory_imaged.dcd")
        assert kept.name == "NVTtrajectory_imaged_part1.dcd"

    def test_nothing_to_preserve_is_not_an_error(self, md, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        assert md.preserve_trajectory("NVTtrajectory.dcd") is None


class TestTheLastFrameWarning:
    def test_it_names_this_stages_own_last_frame(self, md, tmp_path,
                                                 monkeypatch, capsys):
        """Starting from the file this run is about to overwrite is survivable -
        coordinates come from the state - but the recorded `pdb_file` then no
        longer holds what the run started from, which is worth saying."""
        monkeypatch.chdir(tmp_path)
        write_state(tmp_path / "OpenMM_MD_state.xml", 1000, timestep_ps=0.004)
        (tmp_path / "NVTtrajectory_lastframe.pdb").write_text("END\n")
        md.resolve_resume(
            {"resume": True, "nvt_time": 100.0, "nvt_timestep": 0.004},
            pdb_file=str(tmp_path / "NVTtrajectory_lastframe.pdb"),
            time_key="nvt_time", timestep_key="nvt_timestep",
            lastframe="NVTtrajectory_lastframe.pdb")
        assert "will overwrite" in capsys.readouterr().out

    def test_a_missing_last_frame_is_not_a_warning(self, md, tmp_path,
                                                   monkeypatch, capsys):
        """A fresh directory has no last frame, and resolving its path must not
        be what fails."""
        monkeypatch.chdir(tmp_path)
        write_state(tmp_path / "OpenMM_MD_state.xml", 1000, timestep_ps=0.004)
        md.resolve_resume(
            {"resume": True, "nvt_time": 100.0, "nvt_timestep": 0.004},
            pdb_file=str(tmp_path / "start.pdb"),
            time_key="nvt_time", timestep_key="nvt_timestep",
            lastframe="NVTtrajectory_lastframe.pdb")
        assert "will overwrite" not in capsys.readouterr().out


def test_qmmm_still_exports_the_helpers(md, monkeypatch):
    """They moved to `md` and `qmmm` re-exports them; every existing caller and
    test in the tree reaches them through `qmmm`."""
    from mcube import qmmm
    assert qmmm.resolve_resume is md.resolve_resume
    assert qmmm.preserve_trajectory is md.preserve_trajectory
    assert qmmm.read_state_time is md.read_state_time
    assert qmmm.RESTART_STATEFILES == md.RESTART_STATEFILES


class TestItSaysWhenItCouldHaveResumed:
    """The cost of the missing wiring was a run restarting from zero; leaving
    `resume` unset reproduces that by hand, so the run says what it is about to
    discard (D183). `resume` stays opt-in - for a trajectory, continuing and
    replacing are different things, unlike `scan`, where skipping a finished
    point cannot change a result."""

    def test_an_unset_resume_with_a_state_file_present_warns(self, md, tmp_path,
                                                             monkeypatch, capsys):
        monkeypatch.chdir(tmp_path)
        write_state(tmp_path / "OpenMM_MD_state.xml", 1000, timestep_ps=0.004)
        state, remaining = md.resolve_resume(
            {"nvt_time": 100.0, "nvt_timestep": 0.004},
            time_key="nvt_time", timestep_key="nvt_timestep")
        out = capsys.readouterr().out
        assert state is None and remaining == 100.0
        assert "4.000 ps of a previous run" in out
        assert "resume: true" in out and "nvt_time" in out

    def test_a_clean_directory_says_nothing(self, md, tmp_path, monkeypatch, capsys):
        monkeypatch.chdir(tmp_path)
        md.resolve_resume({"nvt_time": 100.0}, time_key="nvt_time")
        assert "WARNING" not in capsys.readouterr().out

    def test_a_truncated_state_file_still_warns(self, md, tmp_path, monkeypatch,
                                                capsys):
        """An unreadable state file is the case most worth mentioning, not
        least worth: something wrote it and then stopped."""
        monkeypatch.chdir(tmp_path)
        (tmp_path / "OpenMM_MD_state.xml").write_text("<State time=")
        md.resolve_resume({"nvt_time": 100.0}, time_key="nvt_time")
        assert "an unreadable amount" in capsys.readouterr().out
