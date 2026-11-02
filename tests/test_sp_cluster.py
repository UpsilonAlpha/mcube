"""`sp` on a gas-phase cluster: pytest tests/test_sp_cluster.py

`sp` is a QM/MM single point on a solvated system, which cannot answer the
question a cluster model is for: whether the semi-empirical method is right
about this chemistry. Imhof's own benchmark - the paper this project's method
follows - found AM1/d made 4→5 coordination 7.55 kcal/mol too favourable and
bidentate→monodentate carboxylate 7 too costly, so the only way to know
whether GFN2 does the same here is to score one small model both ways
(D187, plan item 7, open since D9).

So two things are pinned here. `cluster: true` must build **no MM system** -
a cluster has no point charges, and silently embedding one would make the
benchmark a QM/MM calculation wearing a cluster's name. And `sp_engine` must
reach the same xTB construction every other stage uses, not a second copy of
it, or the thing being benchmarked is not the thing being run.

Needs OpenMM, which the module chain imports at import time.
"""

import sys
import types

import pytest

pytest.importorskip("openmm", reason="sp's module chain imports OpenMM")

from tests.test_qmmm_structure import ASH_NAMES


@pytest.fixture
def sp_mod(monkeypatch):
    """`sp` with ASH stubbed, recording what it was asked to build."""
    calls = {"mm": 0, "qmmm": 0, "singlepoint": [], "fragment": []}

    stub = types.ModuleType("ash")
    for name in ASH_NAMES:
        setattr(stub, name, lambda *a, **k: None)

    class Frag:
        def __init__(self, *a, **k):
            calls["fragment"].append(k)
            self.numatoms = 37

    class Theory:
        def __init__(self, *a, **k):
            self.kind = k.get("kind", "theory")
            self.functional = k.get("functional", "wb97x-v")

    def singlepoint(**k):
        calls["singlepoint"].append(k)
        return types.SimpleNamespace(energy=-1234.5)

    def qmmm(**k):
        calls["qmmm"] += 1
        return Theory()

    stub.Fragment = Frag
    stub.Singlepoint = singlepoint
    stub.QMMMTheory = qmmm
    stub.xTBTheory = lambda **k: Theory(kind="xtb", **{})
    monkeypatch.setitem(sys.modules, "ash", stub)
    for name in [m for m in sys.modules if m == "mcube" or m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)

    from mcube import sp as module
    monkeypatch.setattr(module, "make_mm",
                        lambda *a, **k: calls.__setitem__("mm", calls["mm"] + 1))
    monkeypatch.setattr(module, "make_dft",
                        lambda *a, **k: Theory(kind="dft"))
    monkeypatch.setattr(module, "make_qm",
                        lambda inputs, cores: Theory(kind="xtb:" + str(
                            inputs.get("qm_engine"))))
    module._calls = calls
    return module


def config(tmp_path, **inputs):
    """A config whose `sp` section is a cluster, and a one-atom PDB."""
    import yaml
    pdb = tmp_path / "cluster.pdb"
    pdb.write_text("HETATM    1 ZN    ZN B   1       0.000   0.000   0.000"
                   "  1.00  0.00          ZN\nEND\n")
    cfg = tmp_path / "config.yaml"
    section = {"cluster": True, "pdb_file": str(pdb), "charge": -1, "mult": 1}
    section.update(inputs)
    cfg.write_text(yaml.safe_dump(
        {"preprocess": {"inputs": {}}, "md": {"inputs": {}}, "sp": {"inputs": section}}))
    return types.SimpleNamespace(config=str(cfg), section="sp")


class TestNoMmSystem:
    def test_a_cluster_builds_no_mm_and_no_qmmm_theory(self, sp_mod, tmp_path):
        """The whole point: a cluster has no point charges, and embedding one
        silently would make this a QM/MM number wearing a cluster's name."""
        sp_mod.run(config(tmp_path))
        assert sp_mod._calls["mm"] == 0
        assert sp_mod._calls["qmmm"] == 0

    def test_the_single_point_gets_the_fragment_and_the_charge(self, sp_mod, tmp_path):
        sp_mod.run(config(tmp_path))
        call = sp_mod._calls["singlepoint"][-1]
        assert call["charge"] == -1 and call["mult"] == 1
        assert call["fragment"].numatoms == 37

    def test_the_result_records_that_it_was_a_cluster(self, sp_mod, tmp_path):
        """A QM/MM energy and a gas-phase one are not comparable, so the config
        has to say which this was."""
        out = sp_mod.run(config(tmp_path))
        assert out["cluster"] is True and out["cluster_atoms"] == 37
        assert out["charge"] == -1


class TestTheEngineSwitch:
    def test_dft_is_the_default(self, sp_mod, tmp_path):
        out = sp_mod.run(config(tmp_path))
        assert out["sp_level"] == "wb97x-v/def2-svp"

    def test_xtb_goes_through_make_qm(self, sp_mod, tmp_path):
        """Not a second xTB construction here: the method, electronic
        temperature and SCC options must be the ones every other stage uses, or
        the thing benchmarked is not the thing run."""
        out = sp_mod.run(config(tmp_path, sp_engine="xtb"))
        assert out["sp_level"] == "GFN2-xTB"
        assert sp_mod._calls["singlepoint"][-1]["theory"].kind == "xtb:xtb"

    def test_the_xtb_method_is_reported_not_assumed(self, sp_mod, tmp_path):
        out = sp_mod.run(config(tmp_path, sp_engine="xtb", xtb_method="GFN1"))
        assert out["sp_level"] == "GFN1-xTB"

    def test_an_unknown_engine_is_refused(self, sp_mod, tmp_path):
        with pytest.raises(ValueError) as exc:
            sp_mod.run(config(tmp_path, sp_engine="dftb"))
        assert "sp_engine" in str(exc.value) and "dftb" in str(exc.value)

    def test_the_two_engines_give_two_comparable_records(self, sp_mod, tmp_path):
        """The benchmark is the difference between them on ONE geometry, so both
        must record the same `pdb_file` and differ only in level."""
        a = sp_mod.run(config(tmp_path))
        b = sp_mod.run(config(tmp_path, sp_engine="xtb"))
        assert a["pdb_file"] == b["pdb_file"]
        assert a["sp_level"] != b["sp_level"]
