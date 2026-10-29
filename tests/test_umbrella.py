"""
The umbrella stage: CV geometry, the bias force, WHAM, BAR and the stitch, and the
per-window config bookkeeping:

    pytest tests/test_umbrella.py

ASH is replaced by an empty module (mcube.preprocess imports it), so these run
with OpenMM and mdtraj only; the MD itself is exercised on Bunya.
"""

import math
import sys
import textwrap
import types

import numpy as np
import pytest
import yaml

pytest.importorskip("openmm", reason="the bias is an OpenMM force")

import openmm  # noqa: E402
from openmm import unit  # noqa: E402

KB = 0.0019872041


@pytest.fixture
def umbrella(monkeypatch):
    module = types.ModuleType("ash")
    for name in ("Fragment", "OpenMM_Modeller", "small_molecule_parameterizer", "merge_pdb_files"):
        setattr(module, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", module)
    for name in [m for m in sys.modules if m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    from mcube import umbrella as module_
    return module_


def write(tmp_path, body):
    path = tmp_path / "config.yaml"
    path.write_text(textwrap.dedent(body))
    return str(path)


class TestConfig:
    def test_centres_include_both_ends(self, umbrella):
        assert umbrella.window_centres({"z_start": -1.0, "z_end": 1.0, "spacing": 0.5}) == \
            [-1.0, -0.5, 0.0, 0.5, 1.0]

    def test_bad_range_raises(self, umbrella):
        with pytest.raises(ValueError):
            umbrella.window_centres({"z_start": 1.0, "z_end": -1.0, "spacing": 0.5})

    def test_section_overrides_base_over_defaults(self, umbrella, tmp_path):
        config = write(tmp_path, """
            umbrella:
              inputs: {ion: "C:1:LA", k: 10.0, z_start: -35, z_end: 35, spacing: 0.5}
            umbrella_qmmm:
              inputs: {engine: qmmm, z_start: -8, z_end: 8, timestep: 0.001}
            """)
        inputs = umbrella.umbrella_inputs(config, "umbrella_qmmm")
        assert inputs["engine"] == "qmmm" and inputs["z_start"] == -8
        assert inputs["ion"] == "C:1:LA" and inputs["k"] == 10.0
        assert inputs["cv_frequency"] == umbrella.DEFAULTS["cv_frequency"]

    def test_window_entries_merge_and_inputs_survive(self, umbrella, tmp_path):
        """Window jobs share the config: each must change only its own entry."""
        config = write(tmp_path, "umbrella:\n  inputs: {k: 10.0}\n")
        umbrella.save_nested(config, "umbrella", ["windows", "000"], {"z0": -1.0, "seed": "a.pdb"})
        umbrella.save_nested(config, "umbrella", ["windows", "001"], {"z0": -0.5})
        umbrella.save_nested(config, "umbrella", ["windows", "000"], {"done": True})
        data = yaml.safe_load(open(config))
        assert data["umbrella"]["inputs"] == {"k": 10.0}
        assert data["umbrella"]["results"]["windows"]["000"] == \
            {"z0": -1.0, "seed": "a.pdb", "done": True}
        assert umbrella.find_window(config, "umbrella", -0.5) == {"z0": -0.5}


class TestFrames:
    def test_qmmm_windows_always_keep_frames(self, umbrella, tmp_path):
        config = write(tmp_path, "umbrella:\n  inputs: {}\n")
        assert umbrella.frames_wanted(config, "umbrella_qmmm", {"engine": "qmmm"}, 0.0)

    def test_auto_picks_the_ends_of_the_qmmm_range(self, umbrella, tmp_path):
        config = write(tmp_path, """
            umbrella:
              inputs: {z_start: -35, z_end: 35, spacing: 0.5}
            umbrella_qmmm:
              inputs: {engine: qmmm, z_start: -8, z_end: 8}
            """)
        inputs = umbrella.umbrella_inputs(config)
        wanted = [z for z in umbrella.window_centres(inputs)
                  if umbrella.frames_wanted(config, "umbrella", inputs, z)]
        assert wanted == [-8.0, 8.0]

    def test_auto_without_a_qmmm_section_keeps_none(self, umbrella, tmp_path):
        config = write(tmp_path, "umbrella:\n  inputs: {z_start: -1, z_end: 1, spacing: 0.5}\n")
        assert not umbrella.frames_wanted(config, "umbrella", umbrella.umbrella_inputs(config), 0.0)

    def test_explicit_list(self, umbrella, tmp_path):
        config = write(tmp_path, "umbrella:\n  inputs: {}\n")
        inputs = {"frames_at": [-0.5, 2.0]}
        assert umbrella.frames_wanted(config, "umbrella", inputs, 2.0)
        assert not umbrella.frames_wanted(config, "umbrella", inputs, 0.0)


class TestGeometry:
    MASSES = np.array([12.0, 12.0, 138.9])

    def test_z_is_ion_minus_reference_com(self, umbrella):
        coords = np.array([[0, 0, 0], [0, 0, 2.0], [3.0, 4.0, 11.0]])
        z, r = umbrella.cv_values(coords, 2, [0, 1], self.MASSES, [100, 100, 100])
        assert z[0] == pytest.approx(10.0) and r[0] == pytest.approx(5.0)

    def test_minimum_image(self, umbrella):
        coords = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 95.0]])
        z, _ = umbrella.cv_values(coords, 2, [0, 1], self.MASSES, [100, 100, 100])
        assert z[0] == pytest.approx(-5.0)

    def test_reference_split_across_the_cell(self, umbrella):
        """A trajectory writer may put half the channel a box length away."""
        coords = np.array([[0, 0, 1.0], [0, 0, 99.0], [0, 0, 10.0]])
        z, _ = umbrella.cv_values(coords, 2, [0, 1], self.MASSES, [100, 100, 100])
        assert z[0] == pytest.approx(10.0)

    def test_image_ion_moves_only_the_ion(self, umbrella):
        coords = np.array([[0, 0, 0], [0, 0, 2.0], [0, 0, 91.0]])
        out = umbrella.image_ion(coords, 2, [0, 1], self.MASSES, np.array([100.0] * 3))
        assert out[2, 2] == pytest.approx(-9.0)
        assert np.array_equal(out[:2], coords[:2])

    def test_split_reference_is_refused(self, umbrella):
        coords = np.array([[0, 0, 1.0], [0, 0, 99.0], [0, 0, 10.0]])
        with pytest.raises(ValueError, match="split across the cell"):
            umbrella.image_ion(coords, 2, [0, 1], self.MASSES, np.array([100.0] * 3))

    def test_special_xyz(self, umbrella, tmp_path):
        path = tmp_path / "s.xyz"
        path.write_text("2\nframe\nLa 0 0 1\nC 0 0 0\n2\nframe\nLa 0 0 2\nC 0 0 0\n")
        frames = umbrella.read_special_xyz(str(path))
        assert frames.shape == (2, 2, 3) and frames[1, 0, 2] == 2.0


class TestBias:
    def context(self, umbrella, positions, z0=1.0, k=10.0, cylinder=None, masses=None):
        system = openmm.System()
        for m in (12.0, 12.0, 138.9):
            system.addParticle(m)
        force = umbrella.bias_force(2, [0, 1], k, z0, cylinder, masses)
        system.addForce(force)
        ctx = openmm.Context(system, openmm.VerletIntegrator(0.001),
                             openmm.Platform.getPlatformByName("Reference"))
        ctx.setPositions(np.asarray(positions) / 10.0 * unit.nanometer)
        return ctx, force

    def energy(self, ctx):
        state = ctx.getState(getEnergy=True, groups={31})
        return state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)

    def test_harmonic_in_z(self, umbrella):
        ctx, force = self.context(umbrella, [[0, 0, 0], [0, 0, 2.0], [0, 0, 4.0]], z0=1.0)
        assert self.energy(ctx) == pytest.approx(0.5 * 10.0 * 2.0 ** 2, rel=1e-6)
        z, r = umbrella.cv_now(force, ctx)
        assert z == pytest.approx(3.0) and r == pytest.approx(0.0)

    def test_cylinder_is_flat_inside(self, umbrella):
        ctx, _ = self.context(umbrella, [[0, 0, 0], [0, 0, 2.0], [3.0, 0, 1.0]],
                              z0=0.0, cylinder={"radius": 4.0, "k": 5.0})
        assert self.energy(ctx) == pytest.approx(0.0, abs=1e-9)

    def test_cylinder_outside(self, umbrella):
        ctx, _ = self.context(umbrella, [[0, 0, 0], [0, 0, 2.0], [6.0, 0, 1.0]],
                              z0=0.0, cylinder={"radius": 4.0, "k": 5.0})
        assert self.energy(ctx) == pytest.approx(0.5 * 5.0 * 2.0 ** 2, rel=1e-6)

    def test_window_centre_is_a_parameter(self, umbrella):
        """The pull moves the centre through the context without rebuilding anything."""
        ctx, _ = self.context(umbrella, [[0, 0, 0], [0, 0, 2.0], [0, 0, 4.0]], z0=1.0)
        ctx.setParameter("mcube_us_cv0", 0.3)
        assert self.energy(ctx) == pytest.approx(0.0, abs=1e-9)

    def test_explicit_weights_override_system_masses(self, umbrella):
        """With repartitioned hydrogen masses the System's weights are not the elements'."""
        ctx, force = self.context(umbrella, [[0, 0, 0], [0, 0, 3.0], [0, 0, 4.0]],
                                  masses=[1.0, 2.0, 138.9])
        assert umbrella.cv_now(force, ctx)[0] == pytest.approx(2.0)

    def test_coordination_numbers_are_logged_not_biased(self, umbrella):
        """CN = sum 1/(1+(r/r0)^12): ~1 inside r0, ~0 outside; it adds no energy."""
        system = openmm.System()
        for m in (12.0, 12.0, 138.9, 16.0, 16.0, 16.0):
            system.addParticle(m)
        system.setDefaultPeriodicBoxVectors(*(openmm.Vec3(*v) for v in np.eye(3) * 5.0))
        force = umbrella.bias_force(2, [0, 1], 10.0, 0.0, None, None,
                                    {"water": [3, 4], "protein": [5]}, cn_cutoff=3.0)
        system.addForce(force)
        ctx = openmm.Context(system, openmm.VerletIntegrator(0.001),
                             openmm.Platform.getPlatformByName("Reference"))
        # ion at z = 0 above the reference centre; waters at 2 and 5 A, a protein O at 3 A
        pos = [[0, 0, -1.0], [0, 0, 1.0], [0, 0, 0], [2.0, 0, 0], [0, 5.0, 0], [0, 0, 3.0]]
        ctx.setPositions(np.asarray(pos) / 10.0 * unit.nanometer)
        z, r, cn_water, cn_protein = umbrella.cv_now(force, ctx)
        s = lambda d: 1 / (1 + (d / 3.0) ** 12)
        assert cn_water == pytest.approx(s(2.0) + s(5.0), rel=1e-6)
        assert cn_protein == pytest.approx(0.5, rel=1e-6)
        assert self.energy(ctx) == pytest.approx(0.0, abs=1e-9)

    def test_coordination_release(self, umbrella):
        system = openmm.System()
        system.addParticle(1.0), system.addParticle(1.0)
        force = openmm.CustomBondForce("0.5*k_coord*(r-r0_coord)^2")
        force.addPerBondParameter("k_coord"), force.addPerBondParameter("r0_coord")
        force.addBond(0, 1, [100.0, 0.2])
        system.addForce(force)
        ctx = openmm.Context(system, openmm.VerletIntegrator(0.001),
                             openmm.Platform.getPlatformByName("Reference"))
        ctx.setPositions([openmm.Vec3(0, 0, 0), openmm.Vec3(0.3, 0, 0)] * unit.nanometer)
        umbrella.release(force, ctx)
        assert ctx.getState(getEnergy=True).getPotentialEnergy()._value == pytest.approx(0.0)


def sample_windows(F, z0s, k, T, n, rng):
    """Exact samples from each biased distribution exp(-beta (F + w)) on a fine grid."""
    beta = 1.0 / (KB * T)
    grid = np.linspace(-6, 6, 24001)
    out = []
    for z0 in z0s:
        p = np.exp(-beta * (F(grid) + 0.5 * k * (grid - z0) ** 2))
        cdf = np.cumsum(p)
        cdf /= cdf[-1]
        out.append((z0, k, np.interp(rng.random(n), cdf, grid), None))
    return out


class TestWham:
    T = 300.0

    def test_recovers_a_known_profile(self, umbrella):
        rng = np.random.default_rng(1)
        F = lambda z: 2.0 * np.exp(-z ** 2)          # a 2 kcal/mol barrier
        windows = sample_windows(F, np.arange(-4, 4.01, 0.5), 10.0, self.T, 20000, rng)
        centres, pmf, _ = umbrella.wham(windows, self.T, 0.1)
        keep = np.abs(centres) < 3.5
        ref = F(centres[keep]) - F(centres[keep]).min()
        got = pmf[keep] - np.nanmin(pmf[keep])
        assert np.nanmax(np.abs(got - ref)) < 0.15

    def test_weights_reweight_to_another_profile(self, umbrella):
        """Per-sample exp(-beta dE) weights shift the PMF by dE(z): weighted perturbation."""
        rng = np.random.default_rng(2)
        flat = lambda z: 0.0 * z
        windows = sample_windows(flat, np.arange(-3, 3.01, 0.5), 10.0, self.T, 20000, rng)
        beta = 1.0 / (KB * self.T)
        dE = lambda z: 0.3 * z                       # a linear correction
        weighted = [(z0, k, z, np.exp(-beta * dE(z))) for z0, k, z, _ in windows]
        centres, pmf, _ = umbrella.wham(weighted, self.T, 0.1)
        keep = np.abs(centres) < 2.5
        ref = dE(centres[keep]) - dE(centres[keep]).min()
        assert np.nanmax(np.abs(pmf[keep] - np.nanmin(pmf[keep]) - ref)) < 0.15

    def test_bootstrap_errors_are_small_and_finite(self, umbrella):
        rng = np.random.default_rng(3)
        windows = sample_windows(lambda z: 0 * z, np.arange(-2, 2.01, 0.5), 10.0, self.T, 5000, rng)
        centres, _, err = umbrella.block_bootstrap(windows, self.T, 0.1, n_boot=10)
        assert np.nanmax(err[np.abs(centres) < 1.5]) < 0.2

    def test_overlap_flags_a_gap(self, umbrella):
        rng = np.random.default_rng(4)
        windows = sample_windows(lambda z: 0 * z, [-3.0, -2.5, 2.5], 10.0, self.T, 5000, rng)
        ov = umbrella.overlaps(windows, 0.1)
        assert ov[0] > 0.2 and ov[1] < 0.01


class TestHydration:
    def test_one_state_is_not_bimodal(self, umbrella):
        rng = np.random.default_rng(6)
        assert not umbrella.bimodal(rng.normal(9.0, 0.4, 4000))

    def test_two_separated_states_are(self, umbrella):
        rng = np.random.default_rng(7)
        v = np.concatenate([rng.normal(9.0, 0.3, 2000), rng.normal(6.0, 0.3, 2000)])
        assert umbrella.bimodal(v)

    def test_neighbouring_states_are_one_broad_state(self, umbrella):
        """8 <-> 9 hopping is ordinary exchange, not a frozen second population."""
        rng = np.random.default_rng(8)
        assert not umbrella.bimodal(rng.choice([8.0, 9.0], 4000))

    def test_profile_flags_the_frozen_window(self, umbrella, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        rng = np.random.default_rng(9)
        config = write(tmp_path, "umbrella:\n  inputs: {}\n")
        for key, z0, values in (("000", -1.0, rng.normal(9.0, 0.3, 500)),
                                ("001", 0.0, np.concatenate([rng.normal(9.0, 0.3, 250),
                                                             rng.normal(5.0, 0.3, 250)]))):
            path = tmp_path / f"cv_{key}.csv"
            with open(path, "w") as fh:
                fh.write("time_ps,z_A,r_A,cn_water,cn_protein\n")
                for i, v in enumerate(values):
                    fh.write(f"{i * 0.1},{z0},0.0,{v},0.0\n")
            umbrella.save_nested(config, "umbrella", ["windows", key],
                                 {"z0": z0, "cv": str(path), "done": True, "engine": "mm"})
        result = umbrella.hydration_profile(config, "umbrella")
        assert result["bimodal_at"] == [0.0]
        rows = umbrella.read_csv(result["csv"])
        assert rows["cn_water_mean"][0] == pytest.approx(9.0, abs=0.1)


class TestBar:
    def test_gaussian_work(self, umbrella):
        """For Gaussian work dA = mu - beta sigma^2 / 2, whatever the constant offset."""
        T, mu, sigma, offset = 300.0, 1.5, 0.6, -1.0e5
        beta = 1.0 / (KB * T)
        rng = np.random.default_rng(5)
        wf = rng.normal(mu, sigma, 20000) + offset
        wr = rng.normal(-mu + beta * sigma ** 2, sigma, 20000) - offset
        expected = mu - beta * sigma ** 2 / 2 + offset
        assert umbrella.bar(wf, wr, T) == pytest.approx(expected, abs=0.03)


class TestStitch:
    def test_consistent_profiles_close(self, umbrella):
        z = np.linspace(-10, 10, 201)
        F_mm = 0.1 * z ** 2
        F_qm = F_mm + 0.2 * (z + 10)                   # correction rises 4 kcal/mol across
        grid, F, closure = umbrella.combine_pmfs(z, F_mm, z, F_qm, -5.0, 5.0, 1.0, 3.0)
        assert closure == pytest.approx(0.0, abs=1e-9)
        # MM left of the constriction, QM/MM inside, MM + dA_hi - dA_lo beyond
        at = lambda x: F[np.argmin(np.abs(grid - x))]
        assert at(8.0) - at(-8.0) == pytest.approx(2.0, abs=1e-6)

    def test_inconsistent_profiles_report_closure(self, umbrella):
        z = np.linspace(-10, 10, 201)
        _, _, closure = umbrella.combine_pmfs(z, 0 * z, z, 0 * z, -5.0, 5.0, 0.0, 1.5)
        assert closure == pytest.approx(-1.5)


class TestCli:
    def test_parser(self):
        from mcube.cli import build_parser
        args = build_parser().parse_args(["umbrella", "c.yaml", "window", "--window", "7",
                                          "--section", "umbrella_qmmm"])
        assert (args.action, args.window, args.section) == ("window", 7, "umbrella_qmmm")

    def test_thread_policy(self):
        from mcube.cli import thread_policy
        assert thread_policy("umbrella", {"engine": "mm"}, 8)["bind"] is False
        assert thread_policy("umbrella", {"engine": "qmmm"}, 8)["bind"] is True


class TestEndToEnd:
    """pull -> two MM windows -> MM energies -> WHAM on a small La3+ water box (CPU)."""

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
        return str(pdb)

    def test_pipeline(self, umbrella, tmp_path, monkeypatch):
        pdb = self.build(tmp_path)
        # A compact "channel": the five water oxygens nearest the box centre
        from openmm.app import PDBFile
        system_pdb = PDBFile(pdb)
        xyz = np.array(system_pdb.getPositions()._value) * 10.0
        oxygens = [a.index for a in system_pdb.topology.atoms() if a.residue.name == "HOH" and a.name == "O"]
        reference = [oxygens[i] for i in np.argsort(np.linalg.norm(xyz[oxygens] - 12.0, axis=1))[1:6]]
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
                z_start: -0.5
                z_end: 0.5
                spacing: 0.5
                k: 10.0
                pull_speed: 20000.0
                equilibrate: 0.02
                time: 0.1
                timestep: 0.002
                cv_frequency: 5
                traj_frequency: 10
                frames_at: [-0.5]
                metal_parameters: "1264"
            """)
        monkeypatch.chdir(tmp_path)
        args = types.SimpleNamespace(config=config, section="umbrella", window=None, level="mm",
                                     stride=1, qm_section="umbrella_qmmm", reweight=None, bootstrap=5)
        args.action = "pull"
        umbrella.run(args)
        windows = umbrella.section_windows(config, "umbrella")
        assert set(windows) == {"000", "001", "002"}
        for index in (0, 1):
            monkeypatch.chdir(tmp_path)
            args.action, args.window = "window", index
            entry = umbrella.run(args)
            assert entry["done"] and entry["shell_start"] > 0
            assert entry["frames"] is (index == 0)
            cv = umbrella.read_csv(entry["cv"])
            assert len(cv["z_A"]) == 10
            # A first shell is logged (the box is barely equilibrated, so only a range)
            assert 3.0 < entry["cn"]["water"]["mean"] < 12.0
            assert entry["cn"]["protein"]["mean"] == 0.0
        monkeypatch.chdir(tmp_path)
        args.action, args.window = "energies", 0
        out = umbrella.run(args)
        assert out["frames"] == 5
        monkeypatch.chdir(tmp_path)
        args.window = 1
        with pytest.raises(ValueError, match="saved no full-system frames"):
            umbrella.run(args)
        monkeypatch.chdir(tmp_path)
        args.action = "wham"
        result = umbrella.run(args)
        assert result["windows"] == 2
        data = yaml.safe_load(open(config))
        assert data["umbrella"]["results"]["pull"]["seeds"] == 3
        assert "energies_mm" in data["umbrella"]["results"]["windows"]["000"]


class TestReferenceSelection:
    """`chain A` in a reference selection (mdtraj itself only knows `chainid N`)."""

    def _top(self):
        from openmm.app import Element, Topology
        top = Topology()
        for letter, n in (("A", 3), ("B", 2)):
            chain = top.addChain(letter)
            for i in range(n):
                res = top.addResidue("ALA", chain)
                top.addAtom("CA", Element.getBySymbol("C"), res)
                top.addAtom("CB", Element.getBySymbol("C"), res)
        return top

    def test_chain_letter(self, umbrella):
        pytest.importorskip("mdtraj")
        top = self._top()
        assert umbrella.reference_atoms(top, "chain A and name CA") == [0, 2, 4]
        assert umbrella.reference_atoms(top, "chain B and name CA") == [6, 8]
        assert umbrella.reference_atoms(top, "chainid 1 and name CA") == [6, 8]

    def test_unknown_chain(self, umbrella):
        pytest.importorskip("mdtraj")
        with pytest.raises(ValueError, match="no chain 'Z'"):
            umbrella.reference_atoms(self._top(), "chain Z and name CA")


class TestLoadStructure:
    """md's imaged lastframe has no CONECT (lipids lose their bonds) and splits molecules."""

    def _files(self, tmp_path):
        from openmm import Vec3, unit
        from openmm.app import Element, PDBFile, Topology
        top = Topology()
        chain = top.addChain("A")
        res = top.addResidue("LIG", chain)
        a = top.addAtom("C1", Element.getBySymbol("C"), res)
        b = top.addAtom("C2", Element.getBySymbol("C"), res)
        top.addBond(a, b)
        top.setPeriodicBoxVectors([Vec3(3, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 3)] * unit.nanometer)
        processed = tmp_path / "processed.pdb"
        with open(processed, "w") as fh:
            PDBFile.writeFile(top, [Vec3(0.1, 1, 1), Vec3(0.25, 1, 1)] * unit.nanometer, fh)
        assert any(l.startswith("CONECT") for l in open(processed))
        # the same atoms, C2 wrapped to the far side of the box, no CONECT
        frame = tmp_path / "lastframe.pdb"
        with open(frame, "w") as fh:
            PDBFile.writeFile(top, [Vec3(0.1, 1, 1), Vec3(2.95, 1, 1)] * unit.nanometer, fh)
        frame.write_text("".join(l for l in open(frame) if not l.startswith("CONECT")))
        return str(processed), str(frame)

    def test_bonds_and_whole(self, umbrella, tmp_path):
        processed, frame = self._files(tmp_path)
        pdb = umbrella.load_structure(frame, {"results": {"processed_pdb": processed}})
        assert pdb.topology.getNumBonds() == 1
        xyz = np.array(pdb.getPositions()._value) * 10.0
        assert abs(np.linalg.norm(xyz[1] - xyz[0]) - 1.5) < 1e-3
        assert np.allclose(umbrella.pdb_box(pdb), [30, 30, 30])

    def test_without_processed_pdb(self, umbrella, tmp_path):
        _, frame = self._files(tmp_path)
        pdb = umbrella.load_structure(frame, {})
        assert pdb.topology.getNumBonds() == 0


class TestTallReference:
    """A channel taller than half the box (here 55 A of Cα in an 83 A box) is whole."""

    def test_tall_whole_reference(self, umbrella):
        z = np.arange(0.0, 56.0, 3.8)                      # a chain of Cα up the z axis
        ref = np.stack([np.zeros_like(z), np.zeros_like(z), z + 15.0], axis=1)
        ion = np.array([[0.0, 0.0, 10.0]])
        coords = np.vstack([ref, ion])
        masses = [12.0] * len(ref) + [139.0]
        box = np.array([120.0, 120.0, 83.0])
        out = umbrella.image_ion(coords, len(ref), list(range(len(ref))), masses, box)
        com = ref[:, 2].mean()
        zcv, _ = umbrella.cv_values(out, len(ref), list(range(len(ref))), masses, box)
        assert zcv[0] == pytest.approx(10.0 - com)

    def test_tall_reference_split_by_a_writer(self, umbrella):
        z = np.arange(0.0, 56.0, 3.8)
        ref = np.stack([np.zeros_like(z), np.zeros_like(z), z + 40.0], axis=1)
        wrapped = ref.copy()
        wrapped[:, 2] = wrapped[:, 2] % 83.0               # the top of the chain wraps to the bottom
        coords = np.vstack([wrapped, [[0.0, 0.0, 30.0]]])
        masses = [12.0] * len(ref) + [139.0]
        box = np.array([120.0, 120.0, 83.0])
        zcv, _ = umbrella.cv_values(coords, len(ref), list(range(len(ref))), masses, box)
        assert zcv[0] == pytest.approx(30.0 - ref[:, 2].mean())


class TestDistanceCoordinate:
    """
    `coordinate: distance` - one umbrella stage for a channel's z and for a
    metal-ligand separation (CSR-1's Zn-Asp92 OD1, D133). z stays the default, so
    every existing channel config is unaffected. What is covered here is what a
    wrong unit or a half-converted config would otherwise break silently.
    """

    def test_z_is_the_default(self, umbrella):
        assert umbrella.coordinate_of({}) == "z"
        assert umbrella.coordinate_of({"coordinate": "DISTANCE"}) == "distance"

    def test_an_unknown_coordinate_is_refused(self, umbrella):
        with pytest.raises(ValueError, match="coordinate"):
            umbrella.coordinate_of({"coordinate": "angle"})

    def test_one_place_converts_units(self, umbrella):
        """k in kcal/mol/A^2 -> kJ/mol/nm^2 and the centre in A -> nm, together.

        A factor of 10 missing on the centre, or 418.4 on k, gives a PMF that looks
        plausible and is wrong, so both conversions live in one helper.
        """
        k, cv0 = umbrella.to_openmm(10.0, 2.5)
        assert k == pytest.approx(10.0 * 4.184 * 100.0)
        assert cv0 == pytest.approx(0.25)

    def bias_context(self, umbrella, positions, cv0, k=10.0, pair=(0, 2), **kw):
        system = openmm.System()
        for _ in positions:
            system.addParticle(12.0)
        force = umbrella.bias_force(None, None, k, cv0, coordinate="distance", pair=list(pair), **kw)
        system.addForce(force)
        ctx = openmm.Context(system, openmm.VerletIntegrator(0.001),
                             openmm.Platform.getPlatformByName("Reference"))
        ctx.setPositions(np.asarray(positions, dtype=float) / 10.0 * unit.nanometer)
        return ctx, force

    def energy(self, ctx):
        return ctx.getState(getEnergy=True, groups={31}).getPotentialEnergy().value_in_unit(
            unit.kilocalorie_per_mole)

    def test_harmonic_in_the_pair_separation(self, umbrella):
        """0.5 k (r - cv0)^2 on the two named atoms, and nothing else."""
        ctx, force = self.bias_context(
            umbrella, [[0, 0, 0], [0, 0, 9.0], [0, 0, 4.0]], cv0=2.0)
        assert self.energy(ctx) == pytest.approx(0.5 * 10.0 * 2.0 ** 2, rel=1e-6)
        assert umbrella.cv_now(force, ctx)[0] == pytest.approx(4.0)

    def test_zero_at_the_window_centre(self, umbrella):
        ctx, _ = self.bias_context(umbrella, [[0, 0, 0], [0, 0, 9.0], [3.0, 0, 0]], cv0=3.0)
        assert self.energy(ctx) == pytest.approx(0.0, abs=1e-9)

    def test_the_centre_is_a_settable_parameter(self, umbrella):
        ctx, _ = self.bias_context(umbrella, [[0, 0, 0], [0, 0, 9.0], [0, 0, 4.0]], cv0=2.0)
        ctx.setParameter("mcube_us_cv0", 0.4)
        assert self.energy(ctx) == pytest.approx(0.0, abs=1e-9)

    def test_a_cylinder_is_refused(self, umbrella):
        """There is no axis to confine around, so a leftover cylinder is an error."""
        with pytest.raises(ValueError, match="cylinder"):
            self.bias_context(umbrella, [[0, 0, 0], [0, 0, 9.0], [0, 0, 4.0]], cv0=2.0,
                              cylinder={"radius": 4.0, "k": 5.0})

    def test_the_pair_must_be_exactly_two_atoms(self, umbrella):
        with pytest.raises(ValueError, match="exactly two"):
            umbrella.cv_pair(None, {"atoms": ["B:1:ZN"]})
        with pytest.raises(ValueError, match="exactly two"):
            umbrella.cv_pair(None, {"atoms": ["B:1:ZN", "A:92:OD1", "A:92:OD2"]})

    def test_channel_keys_left_in_the_config_are_refused(self, umbrella):
        """A half-converted config: ignoring `reference` would bias the wrong thing."""
        with pytest.raises(ValueError, match="reference"):
            umbrella.cv_pair(None, {"atoms": ["B:1:ZN", "A:92:OD1"],
                                    "reference": "chain A and name CA"})

    def test_window_bounds_take_either_spelling(self, umbrella):
        neutral = {"cv_start": 1.9, "cv_end": 3.6, "spacing": 0.1}
        channel = {"z_start": 1.9, "z_end": 3.6, "spacing": 0.1}
        assert umbrella.window_centres(neutral) == umbrella.window_centres(channel)
        assert umbrella.window_centres(neutral)[0] == 1.9
        assert umbrella.window_centres(neutral)[-1] == 3.6

    def test_both_spellings_at_once_is_refused(self, umbrella):
        with pytest.raises(ValueError, match="not both"):
            umbrella.window_centres({"cv_start": 1.9, "z_start": 2.0, "cv_end": 3.6,
                                     "spacing": 0.1})

    def test_saved_windows_are_read_by_either_key(self, umbrella):
        """`cv0` now, `z0` in the channel runs already on disk."""
        assert umbrella.entry_cv0({"cv0": 2.5}) == 2.5
        assert umbrella.entry_cv0({"z0": -4.0}) == -4.0
        assert umbrella.entry_cv0({}) is None

    def test_the_cv_column_is_found_by_name(self, umbrella):
        assert umbrella.cv_column({"time_ps": 1, "distance_A": 2}) == 2
        assert umbrella.cv_column({"time_ps": 1, "z_A": 3, "r_A": 4}) == 3
        with pytest.raises(KeyError):
            umbrella.cv_column({"time_ps": 1})

    def test_distance_values_use_the_minimum_image(self, umbrella):
        coords = np.array([[[0.0, 0.0, 0.0], [0.0, 0.0, 9.0]]])
        assert umbrella.distance_values(coords, (0, 1), [10.0, 10.0, 10.0])[0] == pytest.approx(1.0)
        assert umbrella.distance_values(coords, (0, 1), [100.0, 100.0, 100.0])[0] == pytest.approx(9.0)


class TestDistanceCoordinateLimits:
    """
    What the distance coordinate deliberately does NOT do, and the provenance it
    records. From the ion-channel session's review of 2bfa012.
    """

    def test_rescoring_refuses_a_distance_section(self, umbrella, tmp_path):
        """`energies` images frames about the reference chains and writes a z column.

        A distance section cannot even set `ion`, so this used to KeyError partway
        through a rescoring job. It is refused on the config instead, before the window
        state is looked at.
        """
        config = write(tmp_path, """
            umbrella: {inputs: {coordinate: distance, atoms: [B:1:ZN, A:92:OD1]}}
            """)
        with pytest.raises(ValueError, match="coordinate: z"):
            umbrella.run_energies(config, "umbrella",
                                  {"coordinate": "distance", "atoms": ["B:1:ZN", "A:92:OD1"]},
                                  {}, {}, 0, "mm", 1, "umbrella_qmmm")

    def test_recorded_qm_residues_are_the_ones_resolved(self, umbrella, tmp_path):
        """A section's own `extra_qm_residues` wins inside `resolve_qm_region`, so the
        entry has to record the resolved atoms' residues, not preprocess's saved list."""
        from openmm.app import PDBFile, Topology, element
        top = Topology()
        chain = top.addChain("A")
        for resid, name, symbol in ((92, "ASP", "O"), (93, "GLN", "O")):
            residue = top.addResidue(name, chain, id=str(resid))
            top.addAtom("O", element.Element.getBySymbol(symbol), residue)
        pdb = tmp_path / "region.pdb"
        with open(pdb, "w") as fh:
            PDBFile.writeFile(top, [openmm.Vec3(0, 0, 0), openmm.Vec3(0, 0, 0.3)] * unit.nanometer,
                              fh, keepIds=True)
        assert umbrella.residues_of(str(pdb), [0]) == {"A": {92: "ASP"}}
        assert umbrella.residues_of(str(pdb), [0, 1]) == {"A": {92: "ASP", 93: "GLN"}}
