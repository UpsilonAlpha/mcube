"""ohH config derivation: pytest tests/test_make_ohH_config.py

The comparison's validity rests on the two configs differing in three keys, so
what is worth pinning is that it changes those three and carries everything else
across untouched - and that it drops the lysH run's results, since a stale
`results.lastframe` is how a stage silently starts from another run's structure.
"""

import importlib.util

import pytest

yaml = pytest.importorskip("yaml")

from pathlib import Path  # noqa: E402

SCRIPT = Path(__file__).resolve().parents[1] / "docs" / "csr" / "make_ohH_config.py"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("make_ohH_config", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


LYSH = {
    "preprocess": {
        "inputs": {
            "protein_pdb": "./csr_crystal_start.pdb",
            "ligand": "./mer_crystal.pdb",
            "forcefields": ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml",
                            "../hydroxide.xml"],
            "residue_variants": {"A": {53: "HIE", 88: "HID", 164: "HID"}},
            "cutoff": 3.0,
            "require_his_variants": True,
        },
        "timestamp": "2026-10-23T08:21:52",
        "results": {"processed_pdb": "/scratch/.../finalsystem.pdb"},
    },
    "md": {
        "inputs": {"temp": 300, "timestep": 0.004},
        "results": {"lastframe": "/scratch/CSR_nac/md/last.pdb"},
    },
    "qmmm": {"inputs": {"charge": 2, "cores": 4}},
}


@pytest.fixture
def src(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(LYSH, sort_keys=False))
    return path


def written(mod, src, tmp_path, name="CSR_ohH/config.yaml"):
    dest = tmp_path / name
    mod.main([str(src), str(dest)])
    return yaml.safe_load(dest.read_text())


# ---- the three changes ------------------------------------------------------

def test_changes_exactly_the_three_preprocess_keys(mod, src, tmp_path):
    out = written(mod, src, tmp_path)
    pre = out["preprocess"]["inputs"]
    assert pre["protein_pdb"] == "./csr_crystal_start_ohH.pdb"
    assert "../water_nuc.xml" in pre["forcefields"]
    assert "../hydroxide.xml" not in pre["forcefields"]
    assert pre["residue_variants"]["A"][226] == "LYN"


def test_his_variants_survive(mod, src, tmp_path):
    """require_his_variants would reject the run without them, and changing a
    histidine tautomer between the arms would confound the comparison."""
    out = written(mod, src, tmp_path)
    chain_a = out["preprocess"]["inputs"]["residue_variants"]["A"]
    assert {53: "HIE", 88: "HID", 164: "HID"}.items() <= chain_a.items()


def test_everything_else_is_carried_across_untouched(mod, src, tmp_path):
    out = written(mod, src, tmp_path)
    pre = out["preprocess"]["inputs"]
    assert pre["ligand"] == "./mer_crystal.pdb"
    assert pre["cutoff"] == 3.0
    assert pre["require_his_variants"] is True
    assert out["md"]["inputs"] == {"temp": 300, "timestep": 0.004}
    assert out["qmmm"]["inputs"] == {"charge": 2, "cores": 4}
    # Only the force field list changed, not its other entries or their order.
    assert out["preprocess"]["inputs"]["forcefields"][:2] == [
        "amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]


def test_ph_is_recorded(mod, src, tmp_path):
    """It was previously left to preprocess's default, so the assumed pH lived
    in the code rather than in the config."""
    assert written(mod, src, tmp_path)["preprocess"]["inputs"]["pH"] == 8.0


# ---- what it drops ---------------------------------------------------------

def test_results_and_timestamps_are_dropped_at_every_depth(mod, src, tmp_path):
    """A stale results.lastframe is how a stage silently starts from another
    run's structure - the failure mode resolve_start_structure exists for."""
    out = written(mod, src, tmp_path)
    assert "results" not in out["preprocess"]
    assert "timestamp" not in out["preprocess"]
    assert "results" not in out["md"]
    assert yaml.safe_dump(out).find("lastframe") == -1


def test_the_stages_themselves_are_kept(mod, src, tmp_path):
    out = written(mod, src, tmp_path)
    assert set(out) == {"preprocess", "md", "qmmm"}


# ---- refusals --------------------------------------------------------------

def test_refuses_to_overwrite(mod, src, tmp_path):
    dest = tmp_path / "exists.yaml"
    dest.write_text("x: 1")
    with pytest.raises(SystemExit, match="already exists"):
        mod.main([str(src), str(dest)])
    assert dest.read_text() == "x: 1"


def test_refuses_a_config_that_is_not_the_lysh_arm(mod, tmp_path):
    other = dict(LYSH)
    other["preprocess"] = {"inputs": dict(LYSH["preprocess"]["inputs"],
                                          forcefields=["amber14/tip3p.xml"])}
    path = tmp_path / "other.yaml"
    path.write_text(yaml.safe_dump(other))
    with pytest.raises(SystemExit, match="hydroxide.xml"):
        mod.main([str(path), str(tmp_path / "out.yaml")])


def test_refuses_if_226_already_has_a_variant(mod, tmp_path):
    other = {"preprocess": {"inputs": dict(
        LYSH["preprocess"]["inputs"],
        residue_variants={"A": {88: "HID", 226: "LYS"}})}}
    path = tmp_path / "other.yaml"
    path.write_text(yaml.safe_dump(other))
    with pytest.raises(SystemExit, match="226 already"):
        mod.main([str(path), str(tmp_path / "out.yaml")])


def test_refuses_without_chain_a_variants(mod, tmp_path):
    other = {"preprocess": {"inputs": dict(LYSH["preprocess"]["inputs"],
                                           residue_variants={})}}
    path = tmp_path / "other.yaml"
    path.write_text(yaml.safe_dump(other))
    with pytest.raises(SystemExit, match="chain A"):
        mod.main([str(path), str(tmp_path / "out.yaml")])
