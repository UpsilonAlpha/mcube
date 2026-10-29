"""
The metadynamics stage: the n_COO atom selection, config layering, the hills -> bias ->
well-tempered free energy -> integer-state math, the CLI wiring, and a short run of
both modes on a La3+ water box (CPU):

    pytest tests/test_metadynamics.py

ASH is replaced by an empty module (mcube.preprocess imports it), as in test_umbrella.
"""

import sys
import textwrap
import types

import numpy as np
import pytest
import yaml

pytest.importorskip("openmm", reason="the stage runs on OpenMM")

import openmm  # noqa: E402
from openmm import unit  # noqa: E402

KB = 0.0019872041


@pytest.fixture
def mtd(monkeypatch):
    module = types.ModuleType("ash")
    for name in ("Fragment", "OpenMM_Modeller", "small_molecule_parameterizer", "merge_pdb_files"):
        setattr(module, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", module)
    for name in [m for m in sys.modules if m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    from mcube import metadynamics as module_
    return module_


def write(tmp_path, body):
    path = tmp_path / "config.yaml"
    path.write_text(textwrap.dedent(body))
    return str(path)


class TestCarboxylates:
    def topology(self):
        from openmm.app import Element, Topology
        top = Topology()
        chain = top.addChain("A")
        O, C = Element.getBySymbol("O"), Element.getBySymbol("C")
        spec = [("ASP", ["CG", "OD1", "OD2", "O"]), ("GLU", ["CD", "OE1", "OE2"]),
                ("ASH", ["OD1", "OD2"]), ("GLH", ["OE1", "OE2"]), ("SER", ["OG"]),
                ("GLU", ["OE1", "OE2", "OXT"]), ("HOH", ["O"])]
        for name, atoms in spec:
            res = top.addResidue(name, chain)
            for a in atoms:
                top.addAtom(a, O if a.startswith("O") else C, res)
        return top

    def test_charged_asp_glu_side_chain_oxygens_only(self, mtd):
        from mcube.umbrella import carboxylate_oxygens
        top = self.topology()
        names = [(a.residue.name, a.name) for a in top.atoms()]
        picked = [names[i] for i in carboxylate_oxygens(top)]
        # no backbone O, no neutral ASH/GLH, no Ser OG, no C-terminal OXT, no water
        assert picked == [("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2"),
                          ("GLU", "OE1"), ("GLU", "OE2")]

    def test_ion_is_never_its_own_ligand(self, mtd):
        from mcube.umbrella import carboxylate_oxygens
        top = self.topology()
        assert 1 not in carboxylate_oxygens(top, ion=1)


class TestConfig:
    def test_layering(self, mtd, tmp_path):
        config = write(tmp_path, """
            umbrella:
              inputs: {ion: "C:1:LA", reference: "chain A and name CA", k: 5.0,
                       cylinder: {radius: 10.0, k: 10.0}, time: 5000, cv_frequency: 25,
                       c4_set: kantakevicius2022_chbe_tip3p}
            metadynamics:
              inputs: {z0: 0.0, sigma: {n_coo: 0.15}}
            metadynamics_lu:
              inputs: {ion: "C:1:LU", grid: {n_coo: {max: 4.0}}}
            """)
        inputs = mtd.metadynamics_inputs(config, "metadynamics_lu")
        # shared umbrella keys come through; umbrella's run-length keys do not
        assert inputs["reference"] == "chain A and name CA" and inputs["k"] == 5.0
        assert inputs["c4_set"] == "kantakevicius2022_chbe_tip3p"
        assert inputs["time"] == mtd.DEFAULTS["time"]
        assert inputs["cv_frequency"] == mtd.DEFAULTS["cv_frequency"]
        assert inputs["ion"] == "C:1:LU" and inputs["z0"] == 0.0
        # sigma and grid merge per CV
        assert inputs["sigma"] == {"n_coo": 0.15, "z": 0.3}
        assert inputs["grid"]["n_coo"] == {"min": 0.0, "max": 4.0}
        assert mtd.DEFAULTS["grid"]["n_coo"]["max"] == 6.0
        assert mtd.cv_list(inputs) == ["n_coo"]

    def test_n_coo_alone_needs_z0(self, mtd):
        with pytest.raises(ValueError, match="z0"):
            mtd.cv_list(dict(mtd.DEFAULTS))

    def test_z_needs_a_grid(self, mtd):
        with pytest.raises(ValueError, match="grid.z"):
            mtd.cv_list(dict(mtd.DEFAULTS, cvs=["z", "n_coo"]))

    def test_unknown_cv(self, mtd):
        with pytest.raises(ValueError):
            mtd.cv_list(dict(mtd.DEFAULTS, cvs=["rmsd"], z0=0.0))

    def test_grid_axis_matches_openmm_default(self, mtd):
        inputs = dict(mtd.DEFAULTS)
        axis = mtd.grid_axis(inputs, "n_coo")
        assert len(axis) == int(np.ceil(5 * 6.0 / 0.1)) and axis[0] == 0.0 and axis[-1] == 6.0

    def test_seed_prefers_the_umbrella_window_at_z0(self, mtd, tmp_path):
        config = write(tmp_path, """
            umbrella:
              results: {windows: {"070": {z0: 0.0, lastframe: w.pdb}}}
            """)
        md = {"results": {"lastframe": "md.pdb"}}
        assert mtd.run_seed(config, {"z0": 0.0}, md) == "w.pdb"
        assert mtd.run_seed(config, {"z0": 3.0}, md) == "md.pdb"
        assert mtd.run_seed(config, {"z0": 0.0, "seed": "s.pdb"}, md) == "s.pdb"


class TestFes:
    def test_one_hill_on_the_grid(self, mtd):
        x = np.linspace(0, 4, 81)
        V = mtd.bias_on_grid([[2.0]], [3.0], [0.25], [x])
        assert V.max() == pytest.approx(3.0) and x[V.argmax()] == pytest.approx(2.0)
        assert V[np.isclose(x, 2.25)][0] == pytest.approx(3.0 * np.exp(-0.5))

    def test_2d_is_a_product(self, mtd):
        a, b = np.linspace(-1, 1, 21), np.linspace(0, 3, 31)
        V = mtd.bias_on_grid([[0.0, 1.0]], [2.0], [0.5, 0.2], [a, b])
        assert V.shape == (21, 31)
        expect = 2.0 * np.exp(-0.5 * (a[:, None] / 0.5) ** 2) * np.exp(-0.5 * ((b[None, :] - 1) / 0.2) ** 2)
        assert np.allclose(V, expect)

    def test_well_tempered_factor(self, mtd):
        V = np.array([0.0, 1.0, 4.0])
        # gamma = 5: F = -(5/4) V, shifted to min 0
        assert np.allclose(mtd.well_tempered_fes(V, 5.0), [5.0, 3.75, 0.0])

    def test_integer_states_from_two_wells(self, mtd):
        """F with wells at n = 1 and n = 2, the second 2 kcal/mol higher: the state free
        energies follow from the integrals, not from the grid values at 1 and 2."""
        T = 300.0
        kT = KB * T
        x = np.linspace(0, 4, 401)
        F = np.minimum(20 * (x - 1) ** 2, 20 * (x - 2) ** 2 + 2.0)
        states = mtd.integer_states(x, F, kT)
        assert states[1] == 0.0
        assert states[2] == pytest.approx(2.0, abs=0.05)
        assert states[0] > 4 and states[3] > 4
        # a broader well holds more probability: the same depth, twice the width, so it
        # is the lower state by kT ln 2
        F2 = np.minimum(20 * (x - 1) ** 2, 5 * (x - 2) ** 2)
        states = mtd.integer_states(x, F2, kT)
        assert states[2] == 0.0 and states[1] == pytest.approx(kT * np.log(2.0), abs=0.05)

    def test_projection(self, mtd):
        kT = 0.6
        F = np.array([[0.0, 0.0], [1.0, 5.0]])
        assert np.allclose(mtd.project(F, 0, kT), [0.0, -kT * np.log((np.exp(-1 / kT) + np.exp(-5 / kT)) / 2)])

    def test_hills_recover_a_double_well(self, mtd):
        """Converged well-tempered hills, V = -(1 - 1/gamma) F, give F back by both routes."""
        x_true = np.linspace(0, 4, 401)
        F_true = np.minimum(20 * (x_true - 1) ** 2, 20 * (x_true - 2) ** 2 + 2.0)
        gamma, sigma = 8.0, 0.02
        # Gaussians every grid point, heights proportional to the target bias
        V_target = -(1 - 1 / gamma) * F_true
        V_target -= V_target.min()
        spacing = x_true[1] - x_true[0]
        heights = V_target * spacing / (sigma * np.sqrt(2 * np.pi))
        inputs = dict(mtd.DEFAULTS, cvs=["n_coo"], biasfactor=gamma, sigma={"n_coo": sigma},
                      grid={"n_coo": {"min": 0.0, "max": 4.0, "bins": 401}})
        hills = {"time_ps": np.arange(401.0), "n_coo": x_true, "height_kJ": heights / mtd.KJ_TO_KCAL}
        axes, F = mtd.fes_from_hills(hills, ["n_coo"], inputs, 300.0)
        # away from the edges, where the Gaussians are cut off, and the cusp between the
        # wells, which they round off
        inner = (x_true > 0.5) & (x_true < 3.0) & (np.abs(x_true - 1.55) > 0.1)
        assert np.allclose(F[inner], (F_true - F_true.min())[inner], atol=0.05)
        states = mtd.integer_states(axes[0], F, KB * 300.0)
        assert states[2] == pytest.approx(2.0, abs=0.1)
        # half the hills: the first half of the run only
        _, F_half = mtd.fes_from_hills(hills, ["n_coo"], inputs, 300.0, fraction=0.5)
        assert not np.allclose(F_half, F)


class TestHillLog:
    def test_logged_hills_rebuild_openmm_bias(self, mtd):
        """The bias rebuilt from the logged hills is the one OpenMM applied."""
        from openmm.app import Element, Simulation, Topology
        from mcube.umbrella import cn_cv
        top = Topology()
        chain = top.addChain("A")
        for name, el in (("LA", "La"), ("O", "O"), ("O", "O")):
            top.addAtom(name, Element.getBySymbol(el), top.addResidue(name, chain))
        system = openmm.System()
        for m in (138.9, 16.0, 16.0):
            system.addParticle(m)
        system.setDefaultPeriodicBoxVectors(*(openmm.Vec3(*v) for v in np.eye(3) * 3.0))
        inputs = dict(mtd.DEFAULTS, cvs=["n_coo"], z0=0.0, sigma={"n_coo": 0.1},
                      grid={"n_coo": {"min": 0.0, "max": 3.0}})
        var = mtd.bias_variables(inputs, ["n_coo"], 0, [1], None, [1, 2], 3.0)
        meta = mtd.logged_metadynamics(system, var, 300.0, 5.0, 1.0, 5)
        sim = Simulation(top, system, openmm.LangevinMiddleIntegrator(300, 1, 0.002),
                         openmm.Platform.getPlatformByName("Reference"))
        sim.context.setPositions(np.array([[0, 0, 0], [0.25, 0, 0], [0, 0.4, 0]]) * unit.nanometer)
        sim.context.setVelocitiesToTemperature(300)
        meta.clock = lambda: sim.currentStep * 0.002
        meta.step(sim, 50)
        assert len(meta.hills) == 10
        centres = [h[1] for h in meta.hills]
        heights = [h[2] for h in meta.hills]
        V = mtd.bias_on_grid(centres, heights, [0.1], [mtd.grid_axis(inputs, "n_coo")])
        assert np.allclose(V, meta._totalBias, atol=1e-9)
        # well tempered: later hills are lower
        assert heights[0] == pytest.approx(1.0) and heights[-1] < heights[0]


class TestCli:
    def test_parser(self):
        from mcube.cli import build_parser
        args = build_parser().parse_args(["metadynamics", "c.yaml", "fes", "--section", "mtd_lu"])
        assert (args.command, args.action, args.section) == ("metadynamics", "fes", "mtd_lu")
        assert build_parser().parse_args(["metadynamics", "c.yaml", "run"]).section == "metadynamics"

    def test_thread_policy(self):
        from mcube.cli import thread_policy
        assert thread_policy("metadynamics", {}, 8)["bind"] is False


class TestEndToEnd:
    """Both modes on a small La3+ water box, a few water oxygens standing in for carboxylates."""

    def build(self, tmp_path):
        from openmm.app import ForceField, Modeller, PDBFile, Topology, element
        top = Topology()
        chain = top.addChain("A")
        top.addAtom("LA", element.Element.getBySymbol("La"), top.addResidue("LA", chain))
        modeller = Modeller(top, [openmm.Vec3(1.2, 1.2, 1.2)] * unit.nanometer)
        modeller.addSolvent(ForceField("amber14/tip3p.xml"), boxSize=openmm.Vec3(2.4, 2.4, 2.4) * unit.nanometer,
                            positiveIon="Na+", negativeIon="Cl-", neutralize=True)
        pdb = tmp_path / "system.pdb"
        with open(pdb, "w") as fh:
            PDBFile.writeFile(modeller.topology, modeller.positions, fh, keepIds=True)
        from openmm.app import PDBFile as P
        p = P(str(pdb))
        xyz = np.array(p.getPositions()._value) * 10.0
        oxygens = [a.index for a in p.topology.atoms() if a.residue.name == "HOH" and a.name == "O"]
        order = np.argsort(np.linalg.norm(xyz[oxygens] - xyz[0], axis=1))
        return str(pdb), [oxygens[i] for i in order[10:15]], [oxygens[i] for i in order[:4]]

    @pytest.mark.parametrize("mode", ["n_coo", "2d"])
    def test_run_extend_and_fes(self, mtd, tmp_path, monkeypatch, mode):
        pdb, reference, ligands = self.build(tmp_path)
        monkeypatch.setattr(mtd, "carboxylate_oxygens", lambda top, ion=None: ligands)
        extra = ("cvs: [n_coo]\n                z0: 0.0" if mode == "n_coo" else
                 "cvs: [z, n_coo]\n                grid: {z: {min: -3.0, max: 3.0}}\n"
                 "                z_walls: {min: -2.5, max: 2.5, k: 10.0}")
        config = write(tmp_path, f"""
            preprocess:
              inputs: {{metals: [LA], cutoff: 3.0, default_forcefields: [amber14/tip3p.xml]}}
            md:
              inputs: {{temp: 300}}
              results: {{forcefields: [amber14/tip3p.xml], lastframe: {pdb}}}
            umbrella:
              inputs:
                ion: "A:1:LA"
                reference: {reference}
                k: 10.0
                cylinder: {{radius: 5.0, k: 10.0}}
                metal_parameters: "1264"
            metadynamics:
              inputs:
                {extra}
                height: 2.0
                frequency: 5
                time: 0.1
                timestep: 0.002
                cv_frequency: 5
                traj_frequency: 25
            """)
        monkeypatch.chdir(tmp_path)
        args = types.SimpleNamespace(config=config, section="metadynamics", action="run")
        entry = mtd.run(args)
        assert entry["done"] and entry["runs"] == 1 and entry["time_ps"] == pytest.approx(0.1)
        cv = mtd.read_csv(entry["cv"])
        assert list(cv) == ["time_ps", "z_A", "r_A", "cn_water", "cn_protein", "n_coo", "bias_kJ"]
        assert len(cv["time_ps"]) == 10 and cv["n_coo"].max() > 0.5
        hills = mtd.read_csv(entry["hills"])
        assert len(hills["time_ps"]) == 10 and hills["height_kJ"][0] == pytest.approx(2.0)
        assert ("z_A" in hills) is (mode == "2d")

        monkeypatch.chdir(tmp_path)
        entry = mtd.run(args)                       # again: extends the same run
        assert entry["runs"] == 2 and entry["time_ps"] == pytest.approx(0.2)
        cv = mtd.read_csv(entry["cv"])
        assert len(cv["time_ps"]) == 20 and cv["time_ps"][-1] == pytest.approx(0.2)
        assert len(mtd.read_csv(entry["hills"])["time_ps"]) == 20

        monkeypatch.chdir(tmp_path)
        args.action = "fes"
        result = mtd.run(args)
        assert result["hills"] == 20
        assert min(result["n_coo_states_kcal"].values()) == 0.0
        assert set(result["convergence"]) == {"50%", "75%", "100%"}
        assert ("fes_2d" in result) is (mode == "2d") and ("fes_z" in result) is (mode == "2d")
        data = yaml.safe_load(open(config))
        assert data["metadynamics"]["results"]["fes"]["hills"] == 20
        assert data["metadynamics"]["inputs"]["height"] == 2.0


class TestKeySpellings:
    def test_cv0_does_not_stand_in_for_z0(self, mtd):
        """umbrella reads `cv0`; metadynamics is z-only and reads `z0` in `run_seed` and
        for the restraint centre. Accepting `cv0` here would pass validation and then
        silently hold z at 0.0 A (ion-channel session's review of 2bfa012)."""
        inputs = {"cvs": ["n_coo"], "cv0": 2.5,
                  "grid": {"n_coo": {"min": 0.0, "max": 8.0}}, "sigma": {"n_coo": 0.2}}
        with pytest.raises(ValueError, match="z0"):
            mtd.cv_list(inputs)

    def test_an_umbrella_distance_section_is_refused(self, mtd):
        """`coordinate` is a shared input, so without this a distance section would be
        biased along z without comment."""
        inputs = {"cvs": ["z"], "coordinate": "distance",
                  "grid": {"z": {"min": -5.0, "max": 5.0}}, "sigma": {"z": 0.5}}
        with pytest.raises(ValueError, match="umbrella key"):
            mtd.cv_list(inputs)
