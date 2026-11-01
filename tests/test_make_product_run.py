"""make_product_run.py: pytest tests/test_make_product_run.py

The product run is the reactant run with one molecule swapped, and the two are
only comparable if nothing else changed. The things that can go wrong quietly:
the hydroxide being deleted rather than turned into a water, which leaves a
three-coordinate zinc no crystal supports; and a results block copied across,
which would point the new run at the old run's files.

Needs pyyaml only.
"""

import importlib.util
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")

SCRIPT = Path(__file__).resolve().parents[1] / "docs" / "csr" / "make_product_run.py"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("make_product_run", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def line(serial, name, resname, chain, seq, x=0.0):
    return ("HETATM%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00\n"
            % (serial, name, resname, chain, seq, x, 0.0, 0.0))


@pytest.fixture
def start(tmp_path):
    """A protein atom, the zinc, a crystal water, the hydroxide and the substrate."""
    path = tmp_path / "start.pdb"
    path.write_text(
        "CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1\n"
        + line(1, "NE2", "HIS", "A", 88)
        + line(2, "ZN", "ZN", "B", 1, 2.2)
        + line(3, "O", "HOH", "W", 1, 1.9)
        + line(4, "O01", "OH", "C", 2, 2.14)
        + line(5, "H03", "OH", "C", 2, 3.0)
        + line(6, "C4x", "mer", "C", 1, 5.0)
        + line(7, "O2x", "mer", "C", 1, 6.0)
        + "END\n")
    return path


class TestTheProteinFile:
    def test_the_substrate_is_removed(self, mod, start, tmp_path):
        out = tmp_path / "p.pdb"
        mod.product_protein(start, out)
        assert "mer" not in out.read_text()

    def test_the_hydroxide_becomes_a_water(self, mod, start, tmp_path):
        """Not deleted: the position is a crystal water in 6NC5 and deleting it
        leaves a three-coordinate zinc."""
        out = tmp_path / "p.pdb"
        counts = mod.product_protein(start, out)
        text = out.read_text()
        assert counts["converted"] == 1
        assert " OH " not in text
        assert text.count("HOH") == 2          # the crystal water and this one

    def test_the_hydroxide_keeps_its_coordinates(self, mod, start, tmp_path):
        out = tmp_path / "p.pdb"
        mod.product_protein(start, out)
        kept = [l for l in out.read_text().splitlines()
                if l.startswith("HETATM") and l[12:16].strip() == "O01"]
        assert len(kept) == 1 and kept[0][30:38].strip() == "2.140"

    def test_the_hydroxide_proton_is_dropped(self, mod, start, tmp_path):
        """preprocess adds the water's two hydrogens; keeping one of them would
        give the residue three."""
        out = tmp_path / "p.pdb"
        mod.product_protein(start, out)
        assert "H03" not in out.read_text()

    def test_the_protein_the_zinc_and_the_crystal_water_survive(self, mod, start,
                                                                tmp_path):
        out = tmp_path / "p.pdb"
        mod.product_protein(start, out)
        text = out.read_text()
        for keep in ("HIS", " ZN ", "CRYST1"):
            assert keep in text

    def test_it_ends_with_one_END(self, mod, start, tmp_path):
        out = tmp_path / "p.pdb"
        mod.product_protein(start, out)
        assert out.read_text().count("END") == 1


class TestTheConfig:
    def reactant(self):
        return {"preprocess": {"inputs": {"protein_pdb": "./old.pdb",
                                          "ligand": "./mer.pdb",
                                          "ligand_resname": "MER",
                                          "ligand_smiles": "CCO",
                                          "charge": 0, "cutoff": 2.5,
                                          "residue_variants": {"A": {88: "HID"}}},
                               "results": {"processed_pdb": "/other/run.pdb"}},
                "md": {"inputs": {"nvt_time": 1000},
                       "results": {"lastframe": "/other/last.pdb"}},
                "qmmm": {"inputs": {"charge": 2}},
                "scan_aspin": {"inputs": {}, "results": {"energies_csv": "x"}}}

    def test_the_products_own_inputs_are_set(self, mod):
        out = mod.product_config(self.reactant(), "./p.pdb", "./h.pdb", "hyd")
        got = out["preprocess"]["inputs"]
        assert got["protein_pdb"] == "./p.pdb" and got["ligand"] == "./h.pdb"
        assert got["ligand_resname"] == "hyd" and got["charge"] == -1
        assert got["ligand_smiles"] == mod.PRODUCT_SMILES

    def test_the_rest_of_preprocess_is_inherited(self, mod):
        """Same cutoff, same histidines: the two states are comparable only if
        nothing else differs."""
        got = mod.product_config(self.reactant(), "a", "b", "hyd")["preprocess"]["inputs"]
        assert got["cutoff"] == 2.5 and got["residue_variants"] == {"A": {88: "HID"}}

    def test_no_results_are_carried_over(self, mod):
        """A copied results block points the new run at the old run's files."""
        out = mod.product_config(self.reactant(), "a", "b", "hyd")
        assert all("results" not in body for body in out.values())

    def test_md_and_qmmm_inputs_come_across(self, mod):
        out = mod.product_config(self.reactant(), "a", "b", "hyd")
        assert out["md"]["inputs"]["nvt_time"] == 1000
        assert out["qmmm"]["inputs"]["charge"] == 2

    def test_other_sections_are_not_copied(self, mod):
        """`scan_aspin` belongs to the reactant's coordinate, not to this run."""
        assert "scan_aspin" not in mod.product_config(self.reactant(), "a", "b", "hyd")

    def test_the_stale_ligand_xml_is_cleared(self, mod):
        cfg = self.reactant()
        cfg["preprocess"]["inputs"]["ligand_xmlfile"] = "/other/openff_MER.xml"
        got = mod.product_config(cfg, "a", "b", "hyd")["preprocess"]["inputs"]
        assert got["ligand_xmlfile"] is None


class TestTheSmiles:
    def test_it_is_the_minus_one_zwitterion(self, mod):
        """+1 on the pyrrolidine, -1 on each carboxylate: D13's figure."""
        assert mod.PRODUCT_SMILES.count("[O-]") == 2
        assert mod.PRODUCT_SMILES.count("[NH2+]") == 1

    def test_the_leaving_nitrogen_keeps_its_proton(self, mod):
        """The ring N is written bare, so SMILES gives it the implicit H the
        hydrolysis delivers - hyd_zw.pdb's N1x-H03."""
        assert "1NC(" in mod.PRODUCT_SMILES


class TestWhereTheInputsLand:
    """`run.slurm` does `cd "$OUT"` before calling mcube, so preprocess reads
    `protein_pdb` and `ligand` relative to preprocess_results/, not to the
    config. Writing them one level up gives RDKit's opaque "Bad input file"
    six seconds in (D161)."""

    def build(self, mod, tmp_path, start):
        src = tmp_path / "react"
        src.mkdir()
        (src / "config.yaml").write_text(
            "preprocess:@  inputs:@    protein_pdb: ./start.pdb@".replace("@", "\n"))
        (src / "start.pdb").write_text(start.read_text())
        lig = tmp_path / "lig.pdb"
        lig.write_text(
            "HETATM    1  C1x mer C   1       0.000   0.000   0.000\n")
        dst = tmp_path / "prod"
        mod.main(str(src), str(dst), str(lig))
        return dst

    def test_both_inputs_go_into_preprocess_results(self, mod, tmp_path, start):
        dst = self.build(mod, tmp_path, start)
        assert (dst / "preprocess_results" / "product_protein.pdb").is_file()
        assert (dst / "preprocess_results" / "lig.pdb").is_file()
        assert not (dst / "product_protein.pdb").exists()

    def test_the_config_names_them_without_a_directory(self, mod, tmp_path, start):
        """They are read from the stage's own working directory, so the path in
        the config must not carry one."""
        dst = self.build(mod, tmp_path, start)
        inputs = yaml.safe_load(
            (dst / "config.yaml").read_text())["preprocess"]["inputs"]
        assert inputs["protein_pdb"] == "./product_protein.pdb"
        assert inputs["ligand"] == "./lig.pdb"
