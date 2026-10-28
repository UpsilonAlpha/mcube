"""docs/csr/site.py shell mode: pytest tests/test_site_shell.py

The quantity under test is the coordination *number*, which is why this mode
exists at all (D114: semiempirical methods misjudge 4- vs 5-coordinate zinc by
more than this mechanism's barriers differ). So the fixture puts donors at
distances chosen either side of the cutoff, and the tests check that the count
is the one the geometry implies - including that a donor just outside the cutoff
is reported rather than dropped, since dropping it would make the cutoff look
like a measurement.

Needs mdtraj and pyyaml, nothing else: no ASH, no OpenMM.
"""

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "docs" / "csr" / "site.py"


@pytest.fixture(scope="module")
def site_mod():
    spec = importlib.util.spec_from_file_location("site_script", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SITE = {
    "metal": "B:1:ZN",
    "electrophile": "C:1:C4x",
    "carbonyl_o": "C:1:O2x",
    "burgi_dunitz": {"metal_bound": 2.6, "angle": [95, 120], "max_distance": 4.5,
                     "hbond": 2.6},
    "shell": {"cutoff": 2.6, "margin": 1.0, "donor_elements": ["O", "N", "S"]},
}


def atom(serial, name, resname, chain, resseq, x, element):
    """One ATOM record on the x axis, so the distance to a Zn at 0 is x."""
    return ("HETATM%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s"
            % (serial, name, resname, chain, resseq, x, 0.0, 0.0, element))


# Three histidine nitrogens and one hydroxide oxygen inside 2.6 A, one water at
# 2.70 just outside it, and one at 5.00 that is not in the shell on any reading.
ROWS = [
    atom(1, "ZN", "ZN", "B", 1, 0.000, "ZN"),
    atom(2, "NE2", "HIS", "A", 88, 2.100, "N"),
    atom(3, "NE2", "HIS", "A", 164, 2.150, "N"),
    atom(4, "NE2", "HIS", "A", 230, 2.200, "N"),
    atom(5, "O01", "OH", "C", 2, 1.900, "O"),
    atom(6, "O", "HOH", "D", 1, 2.700, "O"),
    atom(7, "O", "HOH", "D", 2, 5.000, "O"),
]


@pytest.fixture
def pdb(tmp_path):
    path = tmp_path / "frame.pdb"
    path.write_text("\n".join(ROWS) + "\nEND\n")
    return str(path)


def run(site_mod, pdb, site=None, capsys=None):
    site_mod.mode_shell(site or SITE, pdb, pdb)
    return capsys.readouterr().out


class TestCoordinationNumber:
    def test_the_count_is_the_donors_inside_the_cutoff(self, site_mod, pdb, capsys):
        out = run(site_mod, pdb, capsys=capsys)
        # Three His N at 2.10-2.20 and O01 at 1.90; the 2.70 water is outside.
        assert "coordination number: mean 4.00" in out
        assert "4: 1 frames (100%)" in out

    def test_a_wider_cutoff_takes_in_the_borderline_water(self, site_mod, pdb,
                                                          capsys):
        site = dict(SITE, shell=dict(SITE["shell"], cutoff=2.8))
        out = run(site_mod, pdb, site, capsys)
        assert "coordination number: mean 5.00" in out

    def test_the_cutoff_defaults_to_metal_bound(self, site_mod, pdb, capsys):
        site = {k: v for k, v in SITE.items() if k != "shell"}
        out = run(site_mod, pdb, site, capsys)
        assert "cutoff 2.60 A" in out
        assert "coordination number: mean 4.00" in out


class TestWhatIsReported:
    def test_every_bound_ligand_is_named_with_its_occupancy(self, site_mod, pdb,
                                                            capsys):
        out = run(site_mod, pdb, capsys=capsys)
        for ligand in ("HIS88:NE2", "HIS164:NE2", "HIS230:NE2", "OH2:O01"):
            assert ligand in out
        assert out.count("100.0") >= 4

    def test_a_donor_just_outside_the_cutoff_is_still_listed(self, site_mod, pdb,
                                                             capsys):
        # 2.70 against a 2.60 cutoff is a borderline call about the site, not an
        # absence, so it is printed - with 0% occupancy and said to be outside.
        out = run(site_mod, pdb, capsys=capsys)
        assert "HOH1:O" in out
        assert "outside the cutoff in every frame" in out

    def test_a_donor_well_beyond_the_margin_is_not_listed(self, site_mod, pdb,
                                                          capsys):
        out = run(site_mod, pdb, capsys=capsys)
        assert "HOH2:O" not in out

    def test_the_metal_itself_is_not_counted_as_its_own_ligand(self, site_mod, pdb,
                                                               capsys):
        out = run(site_mod, pdb, capsys=capsys)
        assert "ZN1:ZN" not in out
