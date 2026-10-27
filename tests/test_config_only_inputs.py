"""Nothing that changes a result comes from the command line:
pytest tests/test_config_only_inputs.py

A flag that carries data means `config.yaml` no longer says what was computed,
which is what reproducing a published number depends on. The concrete failure
that motivated this: `md` had always read its structure from
`md.inputs.pdb_file` and silently ignored `--pdb`, so passing a frame on the
command line equilibrated `preprocess.results.processed_pdb` instead and the
result looked entirely normal. The same shape of bug was available on every
other stage.

Two halves: the parser refuses the removed flags, and each stage takes the same
information from its section instead - including the active region, which used
to be a file path beside the run.

`cli` needs only PyYAML; the rest is imported against a stub `ash`, as in
tests/test_umbrella.py.
"""

import sys
import textwrap
import types

import pytest

from mcube.cli import build_parser

ASH_NAMES = ("Fragment", "OpenMMTheory", "QMMMTheory", "xTBTheory", "PySCFTheory",
             "Singlepoint", "Optimizer", "NumFreq", "NEB", "NEBTS", "OpenMM_Opt",
             "OpenMM_MD", "OpenMM_Modeller", "OpenMM_box_equilibration",
             "Gentle_warm_up_MD", "small_molecule_parameterizer", "merge_pdb_files",
             "read_intlist_from_file", "actregiondefine", "getwaterconstraintslist")


@pytest.fixture
def mc(monkeypatch):
    """mcube's stage modules, imported against a stub `ash`."""
    stub = types.ModuleType("ash")
    for name in ASH_NAMES:
        setattr(stub, name, lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "ash", stub)
    for name in [m for m in sys.modules if m == "mcube" or m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    from mcube import md, neb, preprocess, qmmm
    return types.SimpleNamespace(md=md, neb=neb, preprocess=preprocess, qmmm=qmmm)


def parse(argv):
    return build_parser().parse_args(argv)


# =========================================================
# The parser
# =========================================================

class TestRemovedFlags:
    """argparse exits 2 on an unknown option, so a stale job script fails at
    submit time instead of running with the flag silently dropped."""

    @pytest.mark.parametrize("argv", [
        ["qmmm", "config.yaml", "--pdb", "start.pdb"],
        ["sp", "config.yaml", "--pdb", "start.pdb"],
        ["scan", "config.yaml", "--pdb", "start.pdb"],
        ["ts", "config.yaml", "--pdb", "start.pdb"],
        ["qmmm", "config.yaml", "--actatoms", "active_atoms"],
        ["scan", "config.yaml", "--actatoms", "active_atoms"],
        ["ts", "config.yaml", "--actatoms", "active_atoms"],
        ["neb", "config.yaml", "--actatoms", "active_atoms"],
        ["neb", "config.yaml", "--reactant", "r.pdb", "--product", "p.pdb"],
        ["postprocess", "config.yaml", "--trajectory", "traj.dcd"],
    ])
    def test_a_data_carrying_flag_is_refused(self, argv):
        with pytest.raises(SystemExit) as exc:
            parse(argv)
        assert exc.value.code == 2

    def test_neb_runs_from_the_config_alone(self):
        """It could not before: --product and --reactant were required, so the
        one stage whose inputs were only ever flags was unrunnable from a config."""
        args = parse(["neb", "config.yaml", "--section", "neb_step1"])
        assert args.section == "neb_step1"
        assert not hasattr(args, "reactant") and not hasattr(args, "product")


class TestFlagsThatStay:
    """`--section`, `--qm-section`, `--window` and the umbrella/metadynamics
    actions address the config or the work; they carry no data of their own."""

    @pytest.mark.parametrize("argv, attr, value", [
        (["qmmm", "config.yaml", "--section", "qmmm_product"], "section", "qmmm_product"),
        (["md", "config.yaml", "--section", "md_nac"], "section", "md_nac"),
        (["umbrella", "config.yaml", "window", "--window", "3"], "window", 3),
        (["umbrella", "config.yaml", "energies", "--qm-section", "u_qmmm"],
         "qm_section", "u_qmmm"),
        (["metadynamics", "config.yaml", "fes"], "action", "fes"),
        (["neb", "config.yaml", "--cores", "16"], "cores", 16),
    ])
    def test_still_accepted(self, argv, attr, value):
        assert getattr(parse(argv), attr) == value


# =========================================================
# Paths in `inputs`
# =========================================================

class TestConfigRelativePaths:
    """Each stage is launched from its own output directory, so a relative path
    in `inputs` used to mean a different file depending on which stage read it.
    Reading them against the config makes one path mean one file."""

    def test_a_relative_path_resolves_against_the_config(self, mc, tmp_path):
        config = tmp_path / "run" / "config.yaml"
        config.parent.mkdir()
        got = mc.preprocess.config_relative(str(config), "md_results/last.pdb")
        assert got == (tmp_path / "run" / "md_results" / "last.pdb").as_posix()

    def test_an_absolute_path_is_left_alone(self, mc, tmp_path):
        absolute = (tmp_path / "elsewhere" / "last.pdb").as_posix()
        assert mc.preprocess.config_relative("run/config.yaml", absolute) == absolute

    def test_nothing_stays_nothing(self, mc):
        assert mc.preprocess.config_relative("config.yaml", None) is None


# =========================================================
# The active region
# =========================================================

class TestAtomRanges:
    """A region is hundreds to thousands of indices. As a YAML list that is one
    line per atom, which is why it lived in a file; as ranges it is one line."""

    def test_runs_collapse(self, mc):
        assert mc.preprocess.format_atom_ranges([4203, 4204, 4205, 5001]) == "4203-4205,5001"

    def test_a_round_trip_keeps_every_atom(self, mc):
        atoms = [1, 2, 3, 7, 10, 11]
        text = mc.preprocess.format_atom_ranges(atoms)
        assert mc.preprocess.parse_atom_ranges(text) == atoms

    @pytest.mark.parametrize("value, expected", [
        ("10-12,900", [10, 11, 12, 900]),
        ([10, "12-14"], [10, 12, 13, 14]),
        ([10, 11], [10, 11]),
        ("7", [7]),
    ])
    def test_every_accepted_form(self, mc, value, expected):
        assert mc.preprocess.parse_atom_ranges(value) == expected

    def test_a_selector_says_what_to_use_instead(self, mc):
        """Selectors resolve per structure; an index list does not, which is the
        point of sharing one region between two geometries."""
        with pytest.raises(ValueError, match="actregion_origin"):
            mc.preprocess.parse_atom_ranges("B:1:ZN", "actatoms")

    def test_a_backwards_range_is_refused(self, mc):
        with pytest.raises(ValueError, match="counts backwards"):
            mc.preprocess.parse_atom_ranges("9-1", "actatoms")


class TestActiveRegionFromTheConfig:
    QM = [10, 11, 900]

    def config(self, tmp_path, body=""):
        path = tmp_path / "config.yaml"
        path.write_text(body or "qmmm:\n  inputs: {}\n")
        return str(path)

    def resolve(self, mc, config, inputs):
        return mc.qmmm.resolve_actatoms(config, inputs, None, None, self.QM, [])

    def test_a_region_in_the_config_is_used(self, mc, tmp_path):
        got = self.resolve(mc, self.config(tmp_path), {"actatoms": "10-12"})
        assert got == [10, 11, 12, 900]

    def test_another_section_s_recorded_region_is_reused(self, mc, tmp_path):
        """The reactant, the product and the NEB between them must share one
        region: `actradius` is measured per structure, so the same radius gives
        slightly different regions on two endpoints, and two energies computed
        over different movable atoms are not comparable."""
        config = self.config(tmp_path, "qmmm_reactant:\n  results:\n    actatoms: 10-12\n")
        got = self.resolve(mc, config, {"actatoms_from": "qmmm_reactant"})
        assert got == [10, 11, 12, 900]

    def test_actatoms_wins_over_actatoms_from(self, mc, tmp_path):
        config = self.config(tmp_path, "qmmm_reactant:\n  results:\n    actatoms: 500-501\n")
        got = self.resolve(mc, config, {"actatoms": "10-12",
                                        "actatoms_from": "qmmm_reactant"})
        assert got == [10, 11, 12, 900]

    def test_a_section_with_no_recorded_region_names_it(self, mc, tmp_path):
        """Pointing at a section that has not run yet is the likely mistake, and
        an empty region would silently freeze everything but the QM atoms."""
        with pytest.raises(KeyError) as exc:
            self.resolve(mc, self.config(tmp_path), {"actatoms_from": "qmmm_reactant"})
        text = str(exc.value)
        assert "qmmm_reactant.results.actatoms" in text   # what is missing
        assert "actatoms" in text                         # and what to set instead


# =========================================================
# NEB endpoints
# =========================================================

class TestNebEndpoints:
    """Both minima are required and neither has a default: the two structures
    are the step the barrier belongs to, so a path that is not in the config
    leaves nothing saying what was computed."""

    @pytest.fixture
    def home(self, tmp_path):
        (tmp_path / "config.yaml").write_text("neb:\n  inputs: {}\n")
        for name in ("r.pdb", "p.pdb"):
            (tmp_path / name).write_text("END\n")
        return tmp_path

    def test_both_endpoints_come_from_the_section(self, mc, home):
        react, prod = mc.neb.resolve_endpoints(
            str(home / "config.yaml"), "neb_step1",
            {"reactant_pdb": "r.pdb", "product_pdb": "p.pdb"})
        assert react == (home / "r.pdb").resolve().as_posix()
        assert prod == (home / "p.pdb").resolve().as_posix()

    def test_a_resolved_endpoint_is_pinned_into_the_config(self, mc, home):
        """So a rerun uses the same two structures even if the relative path
        would now resolve elsewhere."""
        config = str(home / "config.yaml")
        mc.neb.resolve_endpoints(config, "neb_step1",
                                 {"reactant_pdb": "r.pdb", "product_pdb": "p.pdb"})
        pinned = mc.preprocess.load_section(config, "neb_step1")["inputs"]
        assert pinned["reactant_pdb"] == (home / "r.pdb").resolve().as_posix()

    def test_a_missing_key_is_named(self, mc, home):
        with pytest.raises(KeyError) as exc:
            mc.neb.resolve_endpoints(str(home / "config.yaml"), "neb_step1",
                                     {"reactant_pdb": "r.pdb"})
        text = str(exc.value)
        assert "product_pdb" in text and "neb_step1.inputs" in text

    def test_a_missing_file_is_named_with_its_key(self, mc, home):
        with pytest.raises(KeyError) as exc:
            mc.neb.resolve_endpoints(str(home / "config.yaml"), "neb_step1",
                                     {"reactant_pdb": "r.pdb", "product_pdb": "gone.pdb"})
        text = str(exc.value)
        assert "product_pdb" in text and "gone.pdb" in text


# =========================================================
# md, the stage the flag actually bit
# =========================================================

class TestMdTakesItsStructureFromTheConfig:
    """`md` used to accept `--pdb` and ignore it, so a frame passed on the
    command line equilibrated `preprocess.results.processed_pdb` instead and
    nothing in the output said so. There is no flag to ignore now, and the
    remaining way to get it wrong - a path relative to the wrong directory -
    fails loudly and names the key."""

    def config(self, tmp_path, pdb_file):
        path = tmp_path / "config.yaml"
        path.write_text(textwrap.dedent("""
            md:
              inputs:
                pdb_file: %s
            preprocess:
              inputs: {}
              results:
                processed_pdb: processed.pdb
            """) % pdb_file)
        return str(path)

    def test_a_missing_structure_names_the_key_and_the_rule(self, mc, tmp_path):
        config = self.config(tmp_path, "md_results/gone.pdb")
        with pytest.raises(FileNotFoundError) as exc:
            mc.md.run(types.SimpleNamespace(config=config))
        text = str(exc.value)
        assert "md.inputs.pdb_file" in text                      # the key at fault
        assert (tmp_path / "md_results").as_posix() in text      # resolved from the config
        assert "config's directory" in text                      # and why it looked there

    def test_a_relative_structure_beside_the_config_is_found(self, mc, tmp_path, monkeypatch):
        """Resolution does not depend on where the job script ran the stage from."""
        (tmp_path / "start.pdb").write_text("END\n")
        config = self.config(tmp_path, "start.pdb")
        (tmp_path / "elsewhere").mkdir()
        monkeypatch.chdir(tmp_path / "elsewhere")
        # Past the structure check; it then fails on the stub, not on the path.
        with pytest.raises(Exception) as exc:
            mc.md.run(types.SimpleNamespace(config=config))
        assert not isinstance(exc.value, FileNotFoundError)


# =========================================================
# The region's provenance
# =========================================================

class Atom:
    """Stands in for an OpenMM topology atom: element, residue and name."""

    def __init__(self, element, resname, name):
        self.element = types.SimpleNamespace(symbol=element)
        self.residue = types.SimpleNamespace(name=resname)
        self.name = name


class Topology:
    def __init__(self, atoms):
        self._atoms = atoms

    def atoms(self):
        return iter(self._atoms)


def topology(swap=False):
    """Three atoms; `swap` moves the proton onto the next residue, as a real edit
    does - the atom count is unchanged, which is why counting cannot catch it."""
    atoms = [Atom("Zn", "ZN", "ZN"), Atom("H", "LYS", "HZ1"), Atom("O", "HOH", "O")]
    if swap:
        atoms[1] = Atom("H", "HOH", "H04")
    return Topology(atoms)


class TestRegionProvenance:
    """`actatoms_from` is a one-line config edit, so the region it copies has to
    carry the structure it was derived from. D115: saved indices stayed
    syntactically valid while the structure under them changed, and the run
    restrained a substrate carbon as the metal."""

    def test_the_fingerprint_changes_when_an_atom_changes_residue(self, mc):
        atoms = [0, 1, 2]
        before = mc.preprocess.atom_fingerprint(topology(), atoms)
        after = mc.preprocess.atom_fingerprint(topology(swap=True), atoms)
        assert before != after

    def test_the_same_structure_gives_the_same_fingerprint(self, mc):
        assert (mc.preprocess.atom_fingerprint(topology(), [0, 2])
                == mc.preprocess.atom_fingerprint(topology(), [0, 2]))

    def test_an_index_past_the_end_says_the_structure_is_a_different_one(self, mc):
        with pytest.raises(ValueError, match="different structure"):
            mc.preprocess.atom_fingerprint(topology(), [0, 9])

    def recorded(self, mc, tmp_path, swap=False):
        """A config in which `qmmm_reactant` has run and recorded its region."""
        mm = types.SimpleNamespace(topology=topology(swap=swap))
        record = mc.qmmm.actatoms_record(mm, str(tmp_path / "reactant.pdb"), [0, 1, 2])
        path = tmp_path / "config.yaml"
        path.write_text("qmmm_reactant:\n  results:\n"
                        + "".join("    %s: %s\n" % (k, v) for k, v in record.items()))
        return str(path)

    def test_reuse_on_the_same_structure_is_allowed(self, mc, tmp_path):
        config = self.recorded(mc, tmp_path)
        mm = types.SimpleNamespace(topology=topology())
        got = mc.qmmm.resolve_actatoms(config, {"actatoms_from": "qmmm_reactant"},
                                       mm, None, [0], [], "product.pdb")
        assert got == [0, 1, 2]

    def test_reuse_on_a_renumbered_structure_is_refused(self, mc, tmp_path):
        """The ohH comparison: two sections of one system whose atom order differs
        by one hydrogen. Freezing the wrong atoms would show up nowhere."""
        config = self.recorded(mc, tmp_path)
        mm = types.SimpleNamespace(topology=topology(swap=True))
        with pytest.raises(ValueError) as exc:
            mc.qmmm.resolve_actatoms(config, {"actatoms_from": "qmmm_reactant"},
                                     mm, None, [0], [], "product.pdb")
        text = str(exc.value)
        assert "qmmm_reactant" in text                 # which section's region
        assert "reactant.pdb" in text                   # what it was built on
        assert "product.pdb" in text                    # what it was applied to
        assert "actregion_origin" in text               # and what to do instead

    def test_a_region_recorded_without_a_fingerprint_is_still_usable(self, mc, tmp_path):
        """Configs from before this check keep working; they just get no check."""
        path = tmp_path / "config.yaml"
        path.write_text("qmmm_reactant:\n  results:\n    actatoms: 0-2\n")
        mm = types.SimpleNamespace(topology=topology(swap=True))
        got = mc.qmmm.resolve_actatoms(str(path), {"actatoms_from": "qmmm_reactant"},
                                       mm, None, [0], [], "product.pdb")
        assert got == [0, 1, 2]
