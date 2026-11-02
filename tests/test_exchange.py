"""
The exchange stage (analysis side of the QM/MM vertical swaps): config layering and checks,
source-balanced pooling against the frame-weighted mean, the periodic re-referencing, and that
every number lands in `exchange.results`. Synthetic run directories with known energies, so
no QM is run:

    pytest tests/test_exchange.py

ASH is replaced by an empty module (mcube.preprocess imports it), as in test_umbrella.
"""

import json
import sys
import textwrap
import types

import pytest
import yaml

HARTREE_KCAL = 627.509474
EV_KCAL = 23.060548


@pytest.fixture
def ex(monkeypatch):
    module = types.ModuleType("ash")
    for name in ("Fragment", "OpenMM_Modeller", "small_molecule_parameterizer", "merge_pdb_files"):
        setattr(module, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", module)
    for name in [m for m in sys.modules if m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    from mcube import exchange
    return exchange


def write_run(root, frames):
    """frames: [(key, kind, source, source_ion, {ion: energy_kcal})] -> run dir with frames.json
    and one job config.yaml per (frame, ion), energies in Hartree as mcube saves them."""
    meta = {}
    for key, kind, source, source_ion, energies in frames:
        jobs = {}
        for ion, e in energies.items():
            job = root / kind / key.split("/")[-1] / ion
            job.mkdir(parents=True)
            (job / "config.yaml").write_text(yaml.safe_dump({
                "qmmm": {"inputs": {"qm_residues": {"A": [1]}},
                         "results": {"opt_energy": e / HARTREE_KCAL, "converged": True}},
                "sp": {"results": {"sp_energy": e / HARTREE_KCAL}}}))
            jobs[ion] = job.as_posix()
        meta[key] = {"kind": kind, "source": source, "source_ion": source_ion, "time_ns": 1.0,
                     "qm_residues": {"A": [1]}, "jobs": jobs}
    (root / "frames.json").write_text(json.dumps(meta))
    return root


def write_config(tmp_path, body):
    path = tmp_path / "config.yaml"
    path.write_text(textwrap.dedent(body))
    return str(path)


def unbalanced_run(root):
    """La -> Lu, two La-sampled site sources and one Lu, water 2 + 2 frames.
    Per source ddE: LA 1, LAchbe 9, LU 2 (balanced 4, sd sqrt(19)); frame-weighted 5."""
    frames = []
    for src, ion, dE in (("LA", "LA", -5.0), ("LAchbe", "LA", 3.0), ("LU", "LU", -10.0)):
        for i in range(2):
            frames.append((f"site/{src}_f{i}", "site", src, ion, {"LA": 0.0, "LU": dE}))
    for src, dE in (("LA", -6.0), ("LU", -12.0)):
        for i in range(2):
            frames.append((f"water/{src}_f{i}", "water", src, src, {"LA": 0.0, "LU": dE}))
    return write_run(root, frames)


def run_stage(ex, config, action, section=None):
    return ex.run(types.SimpleNamespace(config=config, action=action, section=section))


class TestInputs:
    def test_defaults_and_layering(self, ex, tmp_path):
        cfg = write_config(tmp_path, """
            exchange:
              inputs:
                analyse: {electrostatics: periodic}
            inbox:
              inputs:
                analyse: {lra: {lasite: [CAla, LA]}}
                potential: {components: true}
        """)
        inp = ex.exchange_inputs(cfg, "inbox")
        assert inp["analyse"]["electrostatics"] == ["periodic"]
        assert inp["analyse"]["lra"] == {"lasite": ["CAla", "LA"]}
        assert inp["analyse"]["pooling"] == "balanced"          # default kept beside the override
        assert inp["potential"] == {"radii": [8, 10, 15, 20, 30], "pme": True, "components": True}

    @pytest.mark.parametrize("body, match", [
        ("exchange: {inputs: {analyse: {electrostatic: periodic}}}", "unknown keys"),
        ("exchange: {inputs: {analyse: {electrostatics: [R15]}}}", "embed_radii"),
        ("exchange: {inputs: {analyse: {pooling: frames}}}", "pooling"),
        ("exchange: {inputs: {analyse: {lra: {x: [LA]}}}}", "two site sources"),
        ("exchange: {inputs: {outdir: x}}", "unknown keys"),
    ])
    def test_mistyped_inputs_raise(self, ex, tmp_path, body, match):
        with pytest.raises(ValueError, match=match):
            ex.exchange_inputs(write_config(tmp_path, body))

    def test_missing_frames_json(self, ex, tmp_path):
        cfg = write_config(tmp_path, "exchange: {inputs: {run_dir: nowhere}}")
        with pytest.raises(FileNotFoundError, match="frames.json"):
            run_stage(ex, cfg, "analyse")


class TestAnalyse:
    def test_balanced_pooling_and_results_in_config(self, ex, tmp_path):
        unbalanced_run(tmp_path / "run")
        cfg = write_config(tmp_path, "exchange: {inputs: {run_dir: run}}")
        run_stage(ex, cfg, "analyse")
        res = yaml.safe_load(open(cfg))["exchange"]["results"]["analyse"]
        direct = res["treatments"]["direct"]
        lalu = direct["ddE"]["LA->LU"]
        assert lalu["ddE"]["ddE"] == pytest.approx(4.0)
        assert lalu["ddE"]["err"] == pytest.approx(19 ** 0.5, abs=0.01)
        assert lalu["ddE"]["n_sources"] == 3
        assert {k: v["ddE"] for k, v in lalu["by_source"].items()} == pytest.approx(
            {"LA": 1.0, "LAchbe": 9.0, "LU": 2.0})
        report = json.loads(open(direct["report"]).read())
        assert report["levels"]["sp"]["ddE"]["LA->LU"]["all"]["ddE"] == pytest.approx(5.0)
        assert res["pooling"] == "balanced" and res["level"] == "sp"
        assert "CA->LA" not in direct["ddE"]                       # no Ca jobs, no Ca row

    def test_frame_weighted_pooling_on_request(self, ex, tmp_path):
        unbalanced_run(tmp_path / "run")
        cfg = write_config(tmp_path, "exchange: {inputs: {run_dir: run, analyse: {pooling: all}}}")
        res = run_stage(ex, cfg, "analyse")
        assert res["treatments"]["direct"]["ddE"]["LA->LU"]["ddE"]["ddE"] == pytest.approx(5.0)

    def test_periodic_shifts_by_charge_change_times_potential(self, ex, tmp_path):
        run = write_run(tmp_path / "run", [
            ("site/CAx_f0", "site", "CAx", "CA", {"CA": 0.0, "LA": 10.0}),
            ("water/CA_f0", "water", "CA", "CA", {"CA": 0.0, "LA": 4.0})])
        (run / "potential.json").write_text(json.dumps({"frames": {
            "site/CAx_f0": {"phi_direct_V": 1.0, "phi_pme_V": 0.5, "phi_V": {}},
            "water/CA_f0": {"phi_direct_V": 2.0, "phi_pme_V": 2.0, "phi_V": {}}}}))
        cfg = write_config(tmp_path, """
            exchange:
              inputs:
                run_dir: run
                analyse: {electrostatics: [direct, periodic]}
        """)
        res = run_stage(ex, cfg, "analyse")
        t = res["treatments"]
        assert t["direct"]["ddE"]["CA->LA"]["ddE"]["ddE"] == pytest.approx(6.0)
        # dq = +1 e; the site frame's dE moves by -dq * (direct - pme) * 23.06
        assert t["periodic"]["ddE"]["CA->LA"]["ddE"]["ddE"] == pytest.approx(6.0 - 0.5 * EV_KCAL, abs=0.01)
        assert t["periodic"]["report"].endswith("report_periodic.json")

    def test_rerun_replaces_dropped_treatments(self, ex, tmp_path):
        unbalanced_run(tmp_path / "run")
        cfg = write_config(tmp_path, "exchange: {inputs: {run_dir: run}}")
        run_stage(ex, cfg, "analyse")
        data = yaml.safe_load(open(cfg))
        data["exchange"]["results"]["analyse"]["treatments"]["stale"] = {"report": "old"}
        open(cfg, "w").write(yaml.safe_dump(data))
        run_stage(ex, cfg, "analyse")
        assert set(yaml.safe_load(open(cfg))["exchange"]["results"]["analyse"]["treatments"]) == {"direct"}

    def test_section_writes_its_own_results(self, ex, tmp_path):
        unbalanced_run(tmp_path / "run")
        cfg = write_config(tmp_path, """
            exchange: {inputs: {run_dir: run}}
            exchange_all: {inputs: {analyse: {pooling: all}}}
        """)
        run_stage(ex, cfg, "analyse", section="exchange_all")
        data = yaml.safe_load(open(cfg))
        assert "results" not in data["exchange"]
        assert data["exchange_all"]["results"]["analyse"]["pooling"] == "all"


def test_cli_has_the_stage():
    from mcube.cli import build_parser
    args = build_parser().parse_args(["exchange", "config.yaml", "analyse", "--section", "inbox"])
    assert (args.command, args.action, args.section) == ("exchange", "analyse", "inbox")
