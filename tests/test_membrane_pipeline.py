"""
Tests that the preprocess and md stages route a membrane system the right way, and
that the explicit and implicit modes are untouched when `membrane` is absent:

    pytest tests/test_membrane_pipeline.py

ASH is replaced by a module that only records the keyword arguments mcube passes it,
so these tests run without ASH installed and check mcube's side of the interface.
"""

import sys
import types
from pathlib import Path

import pytest

# ASH is stubbed below, but mcube's own imports still reach OpenMM.
pytest.importorskip("openmm", reason="mcube.preprocess imports OpenMM")

import openmm
import pytest
import yaml
from openmm import XmlSerializer, unit
from openmm.app import Topology

# A two-residue protein, numbered 1..N per chain as check_residue_numbering requires
PDB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O
ATOM      5  N   ALA A   2       3.332   1.550   0.000  1.00  0.00           N
ATOM      6  CA  ALA A   2       4.000   2.840   0.000  1.00  0.00           C
ATOM      7  C   ALA A   2       5.500   2.700   0.000  1.00  0.00           C
ATOM      8  O   ALA A   2       6.000   1.580   0.000  1.00  0.00           O
TER
END
"""


class Recorder:
    """Stands in for an ASH function, keeping the keyword arguments it was given."""

    def __init__(self, result=None, side_effect=None):
        self.calls = []
        self.result = result
        self.side_effect = side_effect

    def __call__(self, *args, **kwargs):
        self.calls.append(kwargs)
        if self.side_effect is not None:
            self.side_effect(**kwargs)
        return self.result

    @property
    def kwargs(self):
        assert len(self.calls) == 1, f"expected one call, got {len(self.calls)}"
        return self.calls[0]


class Theory:
    """An OpenMMTheory stand-in that also deserialises the system XML mcube wrote."""

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.system = None
        if kwargs.get("xmlsystemfile"):
            with open(kwargs["xmlsystemfile"]) as f:
                self.system = XmlSerializer.deserialize(f.read())


@pytest.fixture
def ash(monkeypatch):
    """A fake `ash` module, in place before mcube.preprocess / mcube.md are imported."""
    module = types.ModuleType("ash")
    module.OpenMMTheory = Theory
    module.Fragment = Recorder()
    module.OpenMM_Modeller = Recorder(result=(None, None))
    module.small_molecule_parameterizer = Recorder()
    module.merge_pdb_files = Recorder()
    for name in ("OpenMM_Opt", "Gentle_warm_up_MD", "OpenMM_box_equilibration", "OpenMM_MD",
                 "MDtraj_imagetraj", "Singlepoint"):
        setattr(module, name, Recorder())
    monkeypatch.setitem(sys.modules, "ash", module)
    for name in [m for m in sys.modules if m.startswith("mcube.")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    return module


@pytest.fixture
def md(ash):
    from mcube import md as md_module
    return md_module


class FakeForceField:
    """
    Returns an empty periodic System, so the routing can be tested without typing a real
    bilayer (mcube.membrane is tested against OpenMM's POPC patch in test_membrane.py).
    """

    instances = []

    def __init__(self, *xmlfiles):
        self.xmlfiles = list(xmlfiles)
        FakeForceField.instances.append(self)

    def createSystem(self, topology, **kwargs):
        self.createSystem_kwargs = kwargs
        system = openmm.System()
        system.setDefaultPeriodicBoxVectors(*topology.getPeriodicBoxVectors())
        return system


@pytest.fixture
def forcefield(monkeypatch, md):
    FakeForceField.instances = []
    monkeypatch.setattr(md.membrane_mod, "ForceField", FakeForceField)
    monkeypatch.setattr(md, "ForceField", FakeForceField)
    return FakeForceField


def boxed_topology():
    topology = Topology()
    topology.setPeriodicBoxVectors([openmm.Vec3(6.4, 0, 0), openmm.Vec3(0, 6.4, 0),
                                    openmm.Vec3(0, 0, 9.0)] * unit.nanometer)
    return topology


@pytest.fixture
def pdb(monkeypatch, md, tmp_path):
    """A box-carrying topology in place of reading a solvated PDB file."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(md, "PDBFile", lambda path: types.SimpleNamespace(topology=boxed_topology()))


def sections(**md_inputs):
    preprocess = {"inputs": {"default_forcefields": ["amber14/protein.ff14SB.xml",
                                                     "amber14/tip3p.xml"]}}
    return md_inputs, preprocess


# ==============================
# md: membrane
# ==============================

def test_membrane_system_is_passed_to_ash_as_a_system_xml(md, forcefield, pdb):
    inputs, preprocess = sections(membrane=True, temp=300)
    theory, forcefields = md.build_md_theory(inputs, preprocess, "finalsystem.pdb", None, None)

    assert theory.kwargs["xmlsystemfile"] == "system.xml"
    assert theory.kwargs["periodic"] is True
    assert "xmlfiles" not in theory.kwargs
    assert Path("system.xml").is_file()
    # The lipid force field is added, and recorded for the stages that follow
    assert forcefields[-1] == md.membrane_mod.LIPID_XML
    assert forcefield.instances[-1].xmlfiles[-1] == md.membrane_mod.LIPID_XML


def test_the_serialised_membrane_system_has_exactly_one_barostat(md, forcefield, pdb):
    inputs, preprocess = sections(membrane=True, temp=310, pressure=1.0, surface_tension=0.0)
    theory, _ = md.build_md_theory(inputs, preprocess, "finalsystem.pdb", None, None)

    assert md.membrane_mod.barostat_names(theory.system) == ["MonteCarloMembraneBarostat"]
    barostat = theory.system.getForce(0)
    assert barostat.getDefaultTemperature().value_in_unit(unit.kelvin) == 310
    assert barostat.getXYMode() == openmm.MonteCarloMembraneBarostat.XYIsotropic
    assert barostat.getZMode() == openmm.MonteCarloMembraneBarostat.ZFree


def test_membrane_system_uses_pme_and_hbond_constraints(md, forcefield, pdb):
    from openmm.app import HBonds, PME

    inputs, preprocess = sections(membrane=True, temp=300)
    md.build_md_theory(inputs, preprocess, "finalsystem.pdb", None, None)

    kwargs = forcefield.instances[-1].createSystem_kwargs
    assert kwargs["nonbondedMethod"] is PME
    assert kwargs["constraints"] is HBonds
    assert kwargs["rigidWater"] is True
    assert kwargs["hydrogenMass"].value_in_unit(unit.amu) == pytest.approx(1.5)
    assert kwargs["nonbondedCutoff"].value_in_unit(unit.angstrom) == pytest.approx(12.0)


def test_membrane_is_inherited_from_preprocess(md):
    _, preprocess = sections()
    preprocess["inputs"]["membrane"] = True
    assert md.md_membrane({}, preprocess) is True
    # An explicit md.inputs setting wins, in either direction
    assert md.md_membrane({"membrane": False}, preprocess) is False
    assert md.md_membrane({"membrane": True}, {"inputs": {}}) is True


def test_membrane_needs_a_temperature(md, forcefield, pdb):
    inputs, preprocess = sections(membrane=True)
    with pytest.raises(ValueError, match="temp"):
        md.build_md_theory(inputs, preprocess, "finalsystem.pdb", None, None)


# ==============================
# md: the existing modes are unchanged
# ==============================

def test_explicit_solvent_is_unchanged_without_membrane(md, forcefield, pdb):
    for inputs, preprocess in (sections(), sections(membrane=False)):
        theory, forcefields = md.build_md_theory(inputs, preprocess, "finalsystem.pdb", None, None)

        assert theory.kwargs["xmlfiles"] == ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]
        assert theory.kwargs["periodic"] is True
        assert theory.kwargs["rigidwater"] is True
        assert "xmlsystemfile" not in theory.kwargs
        assert forcefields == ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]
        assert not Path("system.xml").exists()  # no system is built for the explicit path


def test_implicit_solvent_is_unchanged_without_membrane(md, forcefield, pdb):
    from openmm.app import CutoffNonPeriodic

    inputs, preprocess = sections(implicit=True, implicit_eps=20)
    theory, forcefields = md.build_md_theory(inputs, preprocess, "finalsystem.pdb", None, None)

    assert theory.kwargs["xmlsystemfile"] == "system.xml"
    assert theory.kwargs["periodic"] is False
    assert theory.kwargs["nonbondedMethod_noPBC"] == "CutoffNonPeriodic"
    assert forcefields[-1] == "implicit/gbn2.xml"
    assert md.membrane_mod.LIPID_XML not in forcefields
    kwargs = forcefield.instances[-1].createSystem_kwargs
    assert kwargs["nonbondedMethod"] is CutoffNonPeriodic
    assert kwargs["solventDielectric"] == 20
    assert md.membrane_mod.barostat_names(theory.system) == []  # implicit MD has no barostat


def test_membrane_and_implicit_together_are_rejected(md, forcefield, pdb):
    inputs, preprocess = sections(membrane=True, implicit=True, temp=300)
    theory, _ = md.build_md_theory(inputs, preprocess, "finalsystem.pdb", None, None)
    # build_md_theory keeps the implicit path (it is tested first); run() refuses the pair
    assert theory.kwargs["periodic"] is False

    config = Path("config.yaml")
    config.write_text(yaml.safe_dump({"md": {"inputs": dict(inputs, pdb_file="finalsystem.pdb")},
                                      "preprocess": preprocess}))
    with pytest.raises(ValueError, match="cannot both be set"):
        md.run(types.SimpleNamespace(config=str(config)))


# ==============================
# preprocess
# ==============================

@pytest.fixture
def preprocess_run(ash, tmp_path, monkeypatch):
    """`mcube preprocess` on a two-residue PDB, with ASH's modeller recorded."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "protein.pdb").write_text(PDB)

    def write_finalsystem(**kwargs):
        (tmp_path / "finalsystem.pdb").write_text(PDB)

    ash.OpenMM_Modeller = Recorder(result=(None, None), side_effect=write_finalsystem)
    from mcube import preprocess as preprocess_module

    def run(**extra_inputs):
        config = tmp_path / "config.yaml"
        inputs = {"protein_pdb": "protein.pdb",
                  "default_forcefields": ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]}
        inputs.update(extra_inputs)
        config.write_text(yaml.safe_dump({"preprocess": {"inputs": inputs}}))
        preprocess_module.OpenMM_Modeller = ash.OpenMM_Modeller
        preprocess_module.run(types.SimpleNamespace(config=str(config)))
        return ash.OpenMM_Modeller.kwargs

    return run


def test_preprocess_builds_a_membrane(preprocess_run):
    kwargs = preprocess_run(membrane=True, membrane_lipidtype="popc", membrane_padding=12.0,
                            membrane_center_z=0.0, ionic_strength=0.15)
    assert kwargs["membrane"] is True
    assert kwargs["membrane_lipidtype"] == "POPC"
    assert kwargs["membrane_padding"] == 12.0
    assert kwargs["membraneCenterZ"] == 0.0
    assert kwargs["watermodel"] == "tip3p"
    assert kwargs["ionicstrength"] == 0.15
    assert "solvent_padding" not in kwargs  # addMembrane sizes the water box itself


def test_preprocess_membrane_defaults(preprocess_run):
    kwargs = preprocess_run(membrane=True)
    assert kwargs["membrane_lipidtype"] == "POPC"
    assert kwargs["membrane_padding"] == 10.0
    assert kwargs["membraneCenterZ"] == 0.0


def test_preprocess_explicit_solvation_is_unchanged(preprocess_run):
    kwargs = preprocess_run()
    assert "membrane" not in kwargs
    assert kwargs["solvent_padding"] == 10.0
    assert kwargs["pH"] == 8.0


def test_preprocess_implicit_solvation_is_unchanged(preprocess_run):
    kwargs = preprocess_run(implicit=True)
    assert kwargs["implicit"] is True
    assert "membrane" not in kwargs
    assert kwargs["pH"] == 7.0


def test_preprocess_rejects_membrane_with_implicit(preprocess_run):
    with pytest.raises(ValueError, match="cannot both be set"):
        preprocess_run(membrane=True, implicit=True)


def test_preprocess_rejects_an_unsupported_lipid(preprocess_run):
    with pytest.raises(ValueError, match="POPG"):
        preprocess_run(membrane=True, membrane_lipidtype="POPG")
