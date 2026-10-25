"""docs/csr/traj_frame.py: pytest tests/test_traj_frame.py

Needs nothing but the standard library: the script writes a PDB template's
coordinates from an XYZ frame, so there is no ASH or OpenMM in the path.
"""

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "docs" / "csr" / "traj_frame.py"


@pytest.fixture(scope="module")
def tf():
    spec = importlib.util.spec_from_file_location("traj_frame", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# Two frames of the same three atoms, the second displaced by +1 in x.
XYZ = """3
Iteration 0 Energy -422.3655127128
O     0.000   0.000   0.000
H     0.957   0.000   0.000
H    -0.240   0.927   0.000
3
Iteration 1 Energy -423.0021856719
O     1.000   0.000   0.000
H     1.957   0.000   0.000
H     0.760   0.927   0.000
"""

PDB = ("HETATM    1  O   HOH E 500       0.000   0.000   0.000  1.00  0.00           O\n"
       "HETATM    2  H1  HOH E 500       0.957   0.000   0.000  1.00  0.00           H\n"
       "HETATM    3  H2  HOH E 500      -0.240   0.927   0.000  1.00  0.00           H\n"
       "END\n")


@pytest.fixture
def files(tmp_path):
    (tmp_path / "t.xyz").write_text(XYZ)
    (tmp_path / "tpl.pdb").write_text(PDB)
    return tmp_path


def run(tf, files, index=None):
    out = files / "out.pdb"
    argv = ["traj_frame.py", str(files / "t.xyz"), str(files / "tpl.pdb"), str(out)]
    if index is not None:
        argv.append(str(index))
    assert tf.main(argv) == 0
    return out.read_text()


class TestFrameChoice:
    def test_the_last_frame_is_the_default(self, tf, files):
        text = run(tf, files)
        assert "   1.000   0.000   0.000" in text      # the displaced O

    def test_an_explicit_index_picks_that_frame(self, tf, files):
        text = run(tf, files, 0)
        assert "   0.000   0.000   0.000" in text

    def test_a_negative_index_counts_from_the_end(self, tf, files):
        assert run(tf, files, -2) == run(tf, files, 0)

    def test_an_index_past_the_end_is_refused(self, tf, files):
        with pytest.raises(IndexError):
            run(tf, files, 7)

    def test_the_energy_is_read_from_the_comment_line(self, tf, files):
        (frame, _) = tf.pick(str(files / "t.xyz"), -1)
        assert tf.energy_of(frame[0]) == pytest.approx(-423.0021856719)


class TestTemplateIsRespected:
    def test_everything_but_the_coordinates_is_kept(self, tf, files):
        text = run(tf, files)
        assert text.rstrip().endswith("END")
        assert "HOH E 500" in text and "1.00  0.00" in text

    def test_a_different_atom_count_is_refused(self, tf, files):
        (files / "tpl.pdb").write_text(PDB.replace(
            "HETATM    3  H2  HOH E 500      -0.240   0.927   0.000  1.00  0.00           H\n", ""))
        with pytest.raises(ValueError, match="has 2 atoms, the frame has 3"):
            run(tf, files)

    def test_a_different_atom_order_is_refused(self, tf, files):
        """The guard that matters: same count, wrong order."""
        (files / "t.xyz").write_text(XYZ.replace("O     1.000", "N     1.000"))
        with pytest.raises(ValueError, match="not the same system"):
            run(tf, files)
