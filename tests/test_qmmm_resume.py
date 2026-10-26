"""QM/MM MD resume: pytest tests/test_qmmm_resume.py

`resume: true` picks a run up where it stopped. Two things make that delicate.
ASH restores positions and velocities from the state file and then steps
`simulation_time` further, so the remaining time is mcube's arithmetic; getting
it wrong silently doubles or truncates an expensive run. And ASH only builds its
DCDReporter with `append=False`, so a restart overwrites the trajectory it
resumed from unless that file is moved aside first.

The state files here are written by OpenMM itself rather than hand-rolled, so
the header these tests parse is the real format.
"""

import os
import sys
import types

import pytest

pytest.importorskip("openmm", reason="resume reads OpenMM state files")

from tests.test_qmmm_structure import ASH_NAMES  # noqa: E402


@pytest.fixture
def qmmm(monkeypatch):
    stub = types.ModuleType("ash")
    for name in ASH_NAMES:
        setattr(stub, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", stub)
    for name in [m for m in sys.modules if m == "mcube" or m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    from mcube import qmmm
    return qmmm


def write_state(path, steps, timestep_ps=0.001):
    """A genuine OpenMM state file advanced `steps` steps."""
    import openmm as mm
    import openmm.unit as u

    system = mm.System()
    for _ in range(3):
        system.addParticle(12.0)
    integrator = mm.LangevinMiddleIntegrator(
        300 * u.kelvin, 1 / u.picosecond, timestep_ps * u.picoseconds)
    context = mm.Context(system, integrator,
                         mm.Platform.getPlatformByName("Reference"))
    context.setPositions([(0, 0, 0), (0.1, 0, 0), (0, 0.1, 0)])
    if steps:
        integrator.step(steps)
    state = context.getState(getPositions=True, getVelocities=True)
    path.write_text(mm.XmlSerializer.serialize(state))
    return path


# ---- reading the elapsed time ------------------------------------------------

def test_reads_time_and_steps_from_a_real_state_file(qmmm, tmp_path):
    state = write_state(tmp_path / "OpenMM_MD_state.xml", 21000)
    elapsed, steps = qmmm.read_state_time(state)
    assert steps == 21000
    assert elapsed == pytest.approx(21.0, abs=1e-6)


def test_header_only_so_the_size_of_the_file_does_not_matter(qmmm, tmp_path):
    """Parsing as XML would read every position and velocity to reach two root
    attributes. Truncating the file past the header must not matter."""
    state = write_state(tmp_path / "OpenMM_MD_state.xml", 5000)
    clipped = tmp_path / "clipped.xml"
    clipped.write_text(state.read_text()[:400])
    assert qmmm.read_state_time(clipped)[0] == pytest.approx(5.0, abs=1e-6)


def test_a_state_file_truncated_inside_its_header_is_refused(qmmm, tmp_path):
    """An over-quota write leaves a file that looks present but is not whole.
    Guessing a time from it would silently shorten the resumed run."""
    state = write_state(tmp_path / "OpenMM_MD_state.xml", 5000)
    text = state.read_text()
    broken = tmp_path / "broken.xml"
    broken.write_text(text[:text.index("time=")])
    with pytest.raises(ValueError, match="truncated|complete OpenMM"):
        qmmm.read_state_time(broken)


def test_a_file_that_is_not_a_state_file_is_refused(qmmm, tmp_path):
    other = tmp_path / "OpenMM_MD_state.xml"
    other.write_text('<?xml version="1.0" ?>\n<System openmmVersion="8.4"></System>')
    with pytest.raises(ValueError, match="complete OpenMM"):
        qmmm.read_state_time(other)


# ---- deciding what to run ----------------------------------------------------

def test_without_resume_the_requested_time_is_passed_through(qmmm, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    write_state(tmp_path / "OpenMM_MD_state.xml", 21000)
    statefile, run_time = qmmm.resolve_resume({"time": 50})
    assert statefile is None and run_time == 50


def test_resume_runs_only_the_remainder(qmmm, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    write_state(tmp_path / "OpenMM_MD_state.xml", 21000)
    statefile, run_time = qmmm.resolve_resume({"time": 50, "resume": True})
    assert statefile.name == "OpenMM_MD_state.xml"
    assert run_time == pytest.approx(29.0, abs=1e-5)


def test_resume_with_nothing_to_resume_from_starts_fresh(qmmm, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    statefile, run_time = qmmm.resolve_resume({"time": 50, "resume": True})
    assert statefile is None and run_time == 50
    assert "nothing to resume from" in capsys.readouterr().out


def test_the_newest_state_file_wins(qmmm, tmp_path, monkeypatch):
    """ASH writes both OpenMM_MD_state.xml and OpenMM_MD_final_state.xml; the
    periodic one is newer than the final one of an earlier, shorter segment."""
    monkeypatch.chdir(tmp_path)
    old = write_state(tmp_path / "OpenMM_MD_final_state.xml", 10000)
    write_state(tmp_path / "OpenMM_MD_state.xml", 30000)
    os.utime(old, (1, 1))
    statefile, run_time = qmmm.resolve_resume({"time": 50, "resume": True})
    assert statefile.name == "OpenMM_MD_state.xml"
    assert run_time == pytest.approx(20.0, abs=1e-5)


def test_a_finished_run_is_not_rerun(qmmm, tmp_path, monkeypatch):
    """A run that reached its target reports 50 ps + 3e-11, not 50 ps, because
    OpenMM accumulates the timestep. Comparing against zero would queue an
    expensive job to integrate a few femtoseconds."""
    monkeypatch.chdir(tmp_path)
    write_state(tmp_path / "OpenMM_MD_state.xml", 50000)
    elapsed, _ = qmmm.read_state_time(tmp_path / "OpenMM_MD_state.xml")
    assert elapsed != 50.0, "the drift this guards against has to be real"
    statefile, run_time = qmmm.resolve_resume({"time": 50, "resume": True})
    assert statefile is not None, "reported so the caller can tell skip from fresh"
    assert run_time == 0.0


def test_a_run_a_single_step_short_is_still_resumed(qmmm, tmp_path, monkeypatch):
    """The cutoff is one timestep, so it must not swallow real remaining time."""
    monkeypatch.chdir(tmp_path)
    write_state(tmp_path / "OpenMM_MD_state.xml", 49990)
    _, run_time = qmmm.resolve_resume(
        {"time": 50, "resume": True, "timestep": 0.001})
    assert run_time == pytest.approx(0.01, abs=1e-5)


def test_resume_without_a_target_time_is_refused(qmmm, tmp_path, monkeypatch):
    """`time` is the total wanted, not the remainder. With no total there is
    nothing to subtract from, and ASH would fall back to its own default."""
    monkeypatch.chdir(tmp_path)
    write_state(tmp_path / "OpenMM_MD_state.xml", 21000)
    with pytest.raises(ValueError, match="total simulation time"):
        qmmm.resolve_resume({"resume": True})


def test_warns_when_the_start_structure_is_the_file_the_md_overwrites(
        qmmm, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    write_state(tmp_path / "OpenMM_MD_state.xml", 21000)
    frame = tmp_path / "QM_MM_lastframe.pdb"
    frame.write_text("END\n")
    qmmm.resolve_resume({"time": 50, "resume": True}, str(frame))
    assert "will overwrite" in capsys.readouterr().out


# ---- not losing the trajectory ----------------------------------------------

def test_the_previous_trajectory_is_moved_aside(qmmm, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "QM_MM.dcd").write_bytes(b"first")
    assert qmmm.preserve_trajectory().name == "QM_MM_part1.dcd"
    assert (tmp_path / "QM_MM_part1.dcd").read_bytes() == b"first"
    assert not (tmp_path / "QM_MM.dcd").exists()


def test_segments_accumulate_rather_than_overwrite_each_other(qmmm, tmp_path, monkeypatch):
    """Two resumes must leave two readable pieces, not one."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "QM_MM.dcd").write_bytes(b"first")
    qmmm.preserve_trajectory()
    (tmp_path / "QM_MM.dcd").write_bytes(b"second")
    qmmm.preserve_trajectory()
    assert (tmp_path / "QM_MM_part1.dcd").read_bytes() == b"first"
    assert (tmp_path / "QM_MM_part2.dcd").read_bytes() == b"second"


def test_nothing_to_preserve_is_not_an_error(qmmm, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert qmmm.preserve_trajectory() is None
