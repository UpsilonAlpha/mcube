"""docs/csr/site.py check mode: pytest tests/test_site_check.py

These are regression tests for two mistakes that were made, not hypotheticals.

D123: everything from October ran on a pose bound the other way round - a
substituent oxygen 3.13 A from the zinc and the carbonyl that has to be attacked
5.43 A away. Both numbers were recorded on the day and neither was read as a
defect, because nothing turned them into a verdict. The first test here is that
pose: it must be rejected.

D120: the zinc went five-coordinate during a 100 ns MM run when an unrestrained
carboxylate oxygen rotated onto it. That one is not caught by checking the
input, only by checking the output, which is why the mode takes a stage and why
there is a test for each end of a run.

Needs pyyaml only: no mdtraj, no ASH, no OpenMM.
"""

import contextlib
import importlib.util
import io
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "docs" / "csr" / "site.py"


@pytest.fixture(scope="module")
def site_mod():
    spec = importlib.util.spec_from_file_location("site_script", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CRITERIA = [
    {"name": "closest substrate O", "closest_to": "B:1:ZN",
     "among": ["C:1:O2x", "C:1:O1x"], "must_be": "C:1:O2x",
     "required_at": ["start", "reactant"]},
    {"name": "Zn - carbonyl O2x", "atoms": ["B:1:ZN", "C:1:O2x"],
     "range": [2.4, 3.0], "required_at": ["start", "reactant"]},
    {"name": "Zn coordination number", "coordination": "B:1:ZN",
     "range": [4, 4], "required_at": ["reactant"]},
    {"name": "O01 - C4x - O2x angle", "atoms": ["C:2:O01", "C:1:C4x", "C:1:O2x"],
     "range": [95, 120], "required_at": ["reactant"]},
]

SITE = {
    "metal": "B:1:ZN",
    "electrophile": "C:1:C4x",
    "carbonyl_o": "C:1:O2x",
    "burgi_dunitz": {"metal_bound": 2.6, "angle": [95, 120], "max_distance": 4.5,
                     "hbond": 2.6},
    "shell": {"cutoff": 2.6, "margin": 1.0, "donor_elements": ["O", "N", "S"]},
    "michaelis_criteria": CRITERIA,
}


def atom(serial, name, resname, chain, resseq, xyz, element):
    return ("HETATM%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s"
            % (serial, name, resname, chain, resseq, xyz[0], xyz[1], xyz[2], element))


def write_pdb(tmp_path, name, rows):
    """rows: (name, resname, chain, resseq, (x, y, z), element), Zn at the origin."""
    lines = [atom(i + 1, *r) for i, r in enumerate(rows)]
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\nEND\n")
    return str(p)


ZN = ("ZN", "ZN", "B", 1, (0.0, 0.0, 0.0), "ZN")
HIS88 = ("NE2", "HIS", "A", 88, (0.0, 2.05, 0.0), "N")
HIS164 = ("NE2", "HIS", "A", 164, (0.0, 0.0, 2.05), "N")


def good_rows(carbonyl=2.55, extra=()):
    """A pose with the beta-lactam carbonyl nearest the metal, as L1's is."""
    return [ZN, HIS88, HIS164,
            ("O01", "OH", "C", 2, (-2.00, 0.0, 0.0), "O"),
            ("O2x", "mer", "C", 1, (carbonyl, 0.0, 0.0), "O"),
            ("O1x", "mer", "C", 1, (13.27, 0.0, 0.0), "O"),
            ("C4x", "mer", "C", 1, (3.90, 1.20, 0.0), "C")] + list(extra)


def run(site_mod, site, path, stage):
    """(exit code, printed lines). The mode exits non-zero on any FAIL."""
    site_mod.OPTS = {"at": stage}
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        with pytest.raises(SystemExit) as e:
            site_mod.mode_check(site, [path])
    return e.value.code, buf.getvalue().splitlines()


def verdicts(lines):
    """{criterion name: PASS|FAIL|WARN} from the printed table."""
    out = {}
    for ln in lines:
        head = ln.split()
        if head and head[0] in ("PASS", "FAIL", "WARN"):
            for c in CRITERIA:
                if c["name"] in ln:
                    out[c["name"]] = head[0]
    return out


class TestTheOctoberPose:
    """D123's pose must be rejected, and for the right reason."""

    def test_flipped_pose_is_rejected(self, site_mod, tmp_path):
        rows = [ZN, HIS88, HIS164,
                ("O01", "OH", "C", 2, (-2.00, 0.0, 0.0), "O"),
                ("O1x", "mer", "C", 1, (3.13, 0.0, 0.0), "O"),
                ("O2x", "mer", "C", 1, (5.43, 0.0, 0.0), "O"),
                ("C4x", "mer", "C", 1, (5.00, 1.50, 0.0), "C")]
        code, lines = run(site_mod, SITE,
                          write_pdb(tmp_path, "oct.pdb", rows), "start")
        assert code == 1
        assert "REJECTED" in "\n".join(lines)
        v = verdicts(lines)
        assert v["closest substrate O"] == "FAIL"
        assert v["Zn - carbonyl O2x"] == "FAIL"
        assert "must be C:1:O2x" in "\n".join(lines)

    def test_pose_2_is_accepted_at_start(self, site_mod, tmp_path):
        code, lines = run(site_mod, SITE,
                          write_pdb(tmp_path, "p2.pdb", good_rows()), "start")
        assert code == 0
        v = verdicts(lines)
        assert v["closest substrate O"] == "PASS"
        assert v["Zn - carbonyl O2x"] == "PASS"

    def test_a_tightly_held_carbonyl_also_fails(self, site_mod, tmp_path):
        """D34 held it at 2.16-2.32; the surviving criterion is LOOSE (D124)."""
        code, _ = run(site_mod, SITE,
                      write_pdb(tmp_path, "tight.pdb", good_rows(carbonyl=2.20)),
                      "start")
        assert code == 1


class TestStages:
    """A criterion that describes the reactant must not reject a start structure."""

    def test_angle_is_advisory_at_start(self, site_mod, tmp_path):
        # good_rows puts O01-C4x-O2x near 50 deg, which is what pose 2 gives.
        path = write_pdb(tmp_path, "ang.pdb", good_rows())
        code, lines = run(site_mod, SITE, path, "start")
        assert code == 0
        assert verdicts(lines)["O01 - C4x - O2x angle"] == "WARN"

    def test_the_same_angle_fails_at_reactant(self, site_mod, tmp_path):
        path = write_pdb(tmp_path, "ang.pdb", good_rows())
        code, lines = run(site_mod, SITE, path, "reactant")
        assert code == 1
        assert verdicts(lines)["O01 - C4x - O2x angle"] == "FAIL"


class TestCoordination:
    """D120: the fifth ligand arrived during the MM, so the output is what to check."""

    def test_four_coordinate_passes(self, site_mod, tmp_path):
        code, lines = run(site_mod, SITE,
                          write_pdb(tmp_path, "four.pdb", good_rows()), "reactant")
        assert verdicts(lines)["Zn coordination number"] == "PASS"

    def test_the_D120_fifth_ligand_fails(self, site_mod, tmp_path):
        """Asp92 OD1 at 1.99 A, the oxygen nothing was watching."""
        od1 = ("OD1", "ASP", "A", 92, (0.0, -1.99, 0.0), "O")
        code, lines = run(site_mod, SITE,
                          write_pdb(tmp_path, "five.pdb", good_rows(extra=[od1])),
                          "reactant")
        assert code == 1
        assert verdicts(lines)["Zn coordination number"] == "FAIL"
        # the detail line names what is on the metal, so the fifth is visible
        assert "A:92:ASP:OD1 1.99" in "\n".join(lines)

    def test_the_fifth_ligand_only_warns_at_start(self, site_mod, tmp_path):
        """The grafted structure inherits that sphere; the restrained MM fixes it."""
        od1 = ("OD1", "ASP", "A", 92, (0.0, -1.99, 0.0), "O")
        code, lines = run(site_mod, SITE,
                          write_pdb(tmp_path, "five.pdb", good_rows(extra=[od1])),
                          "start")
        assert code == 0
        assert verdicts(lines)["Zn coordination number"] == "WARN"


class TestMissingAtoms:
    def test_an_absent_atom_fails_rather_than_crashing(self, site_mod, tmp_path):
        rows = [r for r in good_rows() if r[0] != "O2x"]
        code, lines = run(site_mod, SITE,
                          write_pdb(tmp_path, "no2x.pdb", rows), "start")
        assert code == 1
        assert "not in the structure: C:1:O2x" in "\n".join(lines)

    def test_a_site_with_no_criteria_is_an_error(self, site_mod, tmp_path):
        site = dict(SITE, michaelis_criteria=[])
        site_mod.OPTS = {"at": "start"}
        with pytest.raises(SystemExit) as e:
            site_mod.mode_check(site, [write_pdb(tmp_path, "x.pdb", good_rows())])
        assert "michaelis_criteria" in str(e.value)
