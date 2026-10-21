#!/usr/bin/env python3

from __future__ import annotations
import argparse
import os
import json
from pathlib import Path
from typing import Any, Dict

import openmm
from openmm import XmlSerializer, unit
from openmm.app import PDBFile, ForceField, CutoffNonPeriodic, HBonds

from ash import *

from .preprocess import (
    atom_names,
    compare_coordination,
    coordination_shell,
    report_coordination,
    resolve_cores,
    to_selectors,
    save_result,
    load_section,
    resolve_atoms,
)
from . import membrane as membrane_mod

# ==============================
# Unified run() function
# ==============================

def md_forcefields(inputs: Dict[str, Any], preprocess: Dict[str, Any]) -> list:
    """Force-field XML files: md.inputs overrides, else preprocess defaults + ligand XML."""
    ligand_forcefields = inputs.get("ligand_xml")
    if ligand_forcefields is None:
        ligand_forcefields = (preprocess.get("results") or {}).get("ligand_xml")

    default_forcefields = inputs.get("default_forcefields")
    if default_forcefields is None:
        default_forcefields = preprocess["inputs"]["default_forcefields"]
    default_forcefields = list(default_forcefields)

    if ligand_forcefields:
        default_forcefields.append(ligand_forcefields)
    return default_forcefields


def md_membrane(inputs: Dict[str, Any], preprocess: Dict[str, Any]) -> bool:
    """
    Whether this is a membrane system: `md.inputs.membrane`, else what preprocess built.

    It defaults to the preprocess setting because the two have to agree: a bilayer built
    by preprocess but run without the membrane barostat would be squeezed isotropically.
    """
    value = inputs.get("membrane")
    if value is None:
        value = (preprocess.get("inputs") or {}).get("membrane")
    return bool(value)


def restrain_positions(mm, atom_indices, pdb_file: str, force_constant: float):
    """Harmonic positional restraints, referenced to `pdb_file`'s coordinates.

    Holds the metal site in place while the solvent relaxes around it, which is
    what a non-bonded metal cannot do for itself: over tens of ps at 300 K a
    second-shell donor migrates into the first shell and stays at a distance no
    real metal-ligand bond reaches. Guven et al., J. Phys. Chem. B 2024, 128,
    12416 report the same for VIM-2 with two zinc parameter sets and hold the
    site for the whole classical stage instead of equilibrating it.

    A restraint rather than a constraint, deliberately. Constraints on a metal
    site fail two ways: many of them sharing one centre over-determine the
    constraint network, and OpenMM refuses any constraint touching a frozen
    (massless) atom, which rules out freezing a donor that carries a hydrogen.
    A force has neither problem. It also carries no element-specific
    parameters, so this works for any metal without parameterisation - unlike a
    bonded, dummy-cation or 12-6-4 model.

    The restraint references fixed Cartesian coordinates, so it is in mild
    tension with a barostat rescaling the box around it - more so with the
    anisotropic membrane barostat than an isotropic one. That is the usual
    trade in a restrained equilibration and the reason restraints are stepped
    down or released; nothing here holds the site during the QM/MM stage.
    """
    positions = PDBFile(pdb_file).positions
    force = openmm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter(
        "k", force_constant * unit.kilocalorie_per_mole / unit.angstrom ** 2)
    for name in ("x0", "y0", "z0"):
        force.addPerParticleParameter(name)
    for i in atom_indices:
        xyz = positions[int(i)].value_in_unit(unit.nanometer)
        force.addParticle(int(i), [xyz[0], xyz[1], xyz[2]])
    mm.add_force(force)
    return force


def build_md_theory(inputs: Dict[str, Any], preprocess: Dict[str, Any], pdb_file: str,
                    constraints, restraints, cores=None):
    """
    The classical system used by the md stage: force fields,
    HBond constraints, rigid water, 1.5 amu hydrogens, metal-site bond constraints
    and optional harmonic restraints. Returns (OpenMMTheory, force-field list).
    """
    default_forcefields = md_forcefields(inputs, preprocess)

    if inputs.get("implicit"):
        default_forcefields.append("implicit/gbn2.xml")
        pdb = PDBFile(pdb_file)
        topology = pdb.topology
        forcefield = ForceField(*default_forcefields)
        # ASH cannot change constraints or the nonbonded method of a system read from
        # xmlsystemfile, so they have to be set here.
        system = forcefield.createSystem(
            topology,
            nonbondedMethod=CutoffNonPeriodic,
            nonbondedCutoff=2.0 * unit.nanometer,
            constraints=HBonds,
            hydrogenMass=1.5 * unit.amu,
            soluteDielectric=1.0,
            solventDielectric=inputs.get("implicit_eps"),
        )

        # Assuming 'system' is your OpenMM System object
        with open('system.xml', 'w') as output:
            output.write(XmlSerializer.serialize(system))

        mm = OpenMMTheory(
            xmlsystemfile="system.xml",
            pdbfile=pdb_file,
            periodic=False,
            numcores=cores,
            platform=inputs.get("platform", "CPU"),  # e.g. CUDA on a GPU node
            properties=inputs.get("properties"),
            autoconstraints='HBonds',
            constraints=constraints,
            restraints=restraints,
            nonbondedMethod_noPBC='CutoffNonPeriodic',
            nonbonded_cutoff_noPBC=20,
        )

    elif md_membrane(inputs, preprocess):
        default_forcefields = membrane_mod.lipid_forcefields(
            default_forcefields, inputs.get("membrane_lipid_xml") or membrane_mod.LIPID_XML)
        topology = PDBFile(pdb_file).topology
        # Same reason as the implicit branch above: ASH cannot set constraints or the
        # nonbonded method on a system read from xmlsystemfile, so they are set here. On
        # top of that, ASH can only add an isotropic MonteCarloBarostat, which would tie
        # the bilayer area to its thickness, so the membrane barostat goes in here too.
        system = membrane_mod.membrane_system(
            topology,
            default_forcefields,
            temperature=inputs.get("temp"),
            nonbonded_cutoff=inputs.get("nonbonded_cutoff", membrane_mod.NONBONDED_CUTOFF),
            pressure=inputs.get("pressure", membrane_mod.PRESSURE),
            surface_tension=inputs.get("surface_tension", membrane_mod.SURFACE_TENSION),
            barostat_frequency=inputs.get("barostat_frequency", membrane_mod.BAROSTAT_FREQUENCY),
        )

        with open('system.xml', 'w') as output:
            output.write(XmlSerializer.serialize(system))

        mm = OpenMMTheory(
            xmlsystemfile="system.xml",
            pdbfile=pdb_file,
            periodic=True,
            numcores=cores,
            platform=inputs.get("platform", "CPU"),  # e.g. CUDA on a GPU node
            properties=inputs.get("properties"),
            autoconstraints='HBonds',
            constraints=constraints,
            restraints=restraints,
            rigidwater=True,
        )

    else:
        mm = OpenMMTheory(
            xmlfiles=default_forcefields,
            pdbfile=pdb_file,
            periodic=True,
            numcores=cores,
            platform=inputs.get("platform", "CPU"),  # e.g. CUDA on a GPU node
            properties=inputs.get("properties"),
            autoconstraints='HBonds',
            constraints=constraints,
            restraints=restraints,
            rigidwater=True,
        )
    return mm, default_forcefields


def md_constraints(inputs: Dict[str, Any], preprocess: Dict[str, Any], pdb_file: str):
    """Metal-site bond constraints: md.inputs.constraints, else those found by preprocess."""
    constraints = inputs.get("constraints")
    if constraints is None:
        constraints = (preprocess.get("results") or {}).get("constraints")
    # Atoms may be named as "chain:resid:atom" instead of indices; none (e.g. no metals) -> None
    return resolve_atoms(constraints, pdb_file) or None


def run(args) -> Dict[str, Any]:
    result: Dict[str, Any] = {}

    preprocess = load_section(args.config, "preprocess")
    md = load_section(args.config, "md")
    inputs = md.get("inputs")
    if inputs is None:
        raise KeyError(f"No md.inputs section in {args.config}")

    pdb_file = inputs.get("pdb_file")
    if pdb_file is None:
        pdb_file = preprocess["results"]["processed_pdb"]

    membrane = md_membrane(inputs, preprocess)
    if membrane and inputs.get("implicit"):
        raise ValueError("`membrane` and `implicit` cannot both be set in the md inputs: a bilayer "
                         "needs an explicit, periodic box.")

    constraints = md_constraints(inputs, preprocess, pdb_file)
    # Harmonic distance restraints [[i, j, distance_A, k_kcal_per_mol_A2], ...], e.g. to pull a
    # substrate into place during equilibration. Constraints stay as [i, j] pairs (current distance).
    restraints = resolve_atoms(inputs.get("restraints") or None, pdb_file)
    if restraints:
        print(f"[INFO] {len(restraints)} distance restraints: {restraints}")

    # The coordination sphere this run starts from. Reported again at the end,
    # because a contact that was not there at setup was never constrained, and a
    # non-bonded metal will happily pull one in.
    preprocess_inputs = preprocess.get("inputs") or {}
    shell_before = coordination_shell(pdb_file, preprocess_inputs)
    report_coordination(shell_before, "before MD")
    result["coordination_start"] = shell_before

    # `restrain_coordination` holds the detected metal site with positional
    # restraints instead of constraining it. The same atoms either way, from the
    # same detection; see restrain_positions for why a restraint and not a
    # constraint. Applied after the theory is built, so the constraints are not.
    restrained_atoms = []
    if inputs.get("restrain_coordination"):
        restrained_atoms = sorted({int(a) for c in (constraints or []) for a in c[:2]})
        constraints = None

    mm, default_forcefields = build_md_theory(inputs, preprocess, pdb_file, constraints, restraints,
                                              cores=resolve_cores(inputs.get("cores")))

    if restrained_atoms:
        force_constant = float(inputs.get("coordination_force") or 10.0)
        restrain_positions(mm, restrained_atoms, pdb_file, force_constant)
        print(f"[INFO] restrain_coordination: {len(restrained_atoms)} atoms held by "
              f"harmonic positional restraint at {force_constant} kcal/mol/A^2 "
              "(the metal and the donors within constraint_cutoff); the metal bond "
              "constraints are not applied")
        names = atom_names(pdb_file)
        result["restrained_atoms"] = [names[i] for i in restrained_atoms]

    fragment = Fragment(pdbfile=pdb_file)
    result["forcefields"] = default_forcefields
    # Recorded as selectors, not the indices OpenMM was given: indices only mean
    # something against the structure they were resolved on.
    result["restraints"] = to_selectors(restraints, pdb_file)

    print("[STEP] Minimization...")
    OpenMM_Opt(fragment=fragment, theory=mm, maxiter=500, tolerance=1)
    result["minimized"] = Path("frag-minimized.pdb").resolve().as_posix()
    save_result(result, filename=args.config, section="md")

    if inputs.get("gentle"):
        Gentle_warm_up_MD(theory=mm, fragment=fragment, time_steps=[0.0005,0.001,0.004],
                    steps=[10,50,10000], temperatures=[1,10,300])
        result["gentle"] = Path("warmup_MD_cycle2_lastframe.pdb").resolve().as_posix()
        fragment = Fragment(pdbfile=result["gentle"])
        save_result(result, filename=args.config, section="md")

    if inputs.get("implicit") and not inputs.get("skip_npt"):
        print("[INFO] Implicit solvent has no periodic box: skipping NPT equilibration")

    if membrane and not inputs.get("skip_npt"):
        print("[STEP] NPT Equilibration (membrane barostat)...")
        # Not OpenMM_box_equilibration: it always asks ASH for a MonteCarloBarostat, which
        # is isotropic, and ASH's "already present" check looks for that exact class name,
        # so it would not see the MonteCarloMembraneBarostat and the system would end up
        # with two barostats. OpenMM_MD without `barostat=` leaves the membrane barostat
        # that is already in the serialised system as the only one; ASH's own removal pass
        # in that case drops only MonteCarloBarostat/AndersenThermostat, not the membrane
        # barostat. Checked either side of the run.
        membrane_mod.check_single_barostat(mm.system)
        OpenMM_MD(
            fragment=fragment,
            theory=mm,
            timestep=inputs.get("npt_timestep"),
            simulation_time=inputs.get("npt_time", 1000),
            traj_frequency=inputs.get("traj_frequency", 1000),
            temperature=inputs.get("temp"),
            integrator='LangevinMiddleIntegrator',
            coupling_frequency=1,
            trajfilename='relaxbox_NPT',
            trajectory_file_option='DCD',
            datafilename="nptsim.csv",
        )
        membrane_mod.check_single_barostat(mm.system)
        result["npt_csv"] = Path("nptsim.csv").resolve().as_posix()
        # Box size and area per lipid over the run. Diagnostic only, so a failure here
        # must not throw away an equilibrated box.
        try:
            result["membrane"] = membrane_mod.membrane_report(
                Path("relaxbox_NPT_lastframe.pdb").resolve().as_posix(), "relaxbox_NPT.dcd")
        except Exception as exc:
            print(f"[WARNING] Could not report the membrane box size: {exc}")
        save_result(result, filename=args.config, section="md")

    if not inputs.get("skip_npt") and not inputs.get("implicit") and not membrane:
        print("[STEP] NPT Equilibration...")
        OpenMM_box_equilibration(
            fragment=fragment,
            theory=mm,
            datafilename="nptsim.csv",
            numsteps_per_NPT=10000,
            volume_threshold=1.0,
            density_threshold=0.001,
            temperature=inputs.get("temp"),
            timestep=inputs.get("npt_timestep"),
            traj_frequency=100,
            trajfilename='relaxbox_NPT',
            trajectory_file_option='DCD',
            coupling_frequency=1,
        )
    
    if not inputs.get("skip_nvt"):
        print("[STEP] Production MD...")
        if membrane:
            # The membrane barostat stays in the system, so production runs at constant
            # pressure and zero surface tension rather than at constant volume. That is
            # the right ensemble for a bilayer: it keeps the area per lipid free.
            print("[INFO] Membrane system: production MD keeps the membrane barostat (NPT, "
                  "tensionless), so the box area and thickness stay free")
            membrane_mod.check_single_barostat(mm.system)
        OpenMM_MD(
            fragment=fragment,
            theory=mm,
            timestep=inputs.get("nvt_timestep"),
            simulation_time=inputs.get("nvt_time"),
            traj_frequency=inputs.get("traj_frequency", 1000),
            temperature=inputs.get("temp"),
            integrator='LangevinMiddleIntegrator',
            coupling_frequency=1,
            trajfilename='NVTtrajectory',
            trajectory_file_option='DCD',
            datafilename="nvtsim.csv",
        )
    
        result["csv"] = Path("nvtsim.csv").resolve().as_posix()

        if inputs.get("implicit"):
            # No unit cell to reimage into
            result["trajfile"] = Path("NVTtrajectory.dcd").resolve().as_posix()
            result["lastframe"] = Path("NVTtrajectory_lastframe.pdb").resolve().as_posix()
        else:
            print("[STEP] Reimaging trajectory...")
            MDtraj_imagetraj(
                "NVTtrajectory.dcd",
                "NVTtrajectory_lastframe.pdb",
                format='DCD',
            )
            result["trajfile"] = Path("NVTtrajectory_imaged.dcd").resolve().as_posix()
            result["lastframe"] = Path("NVTtrajectory_lastframe_imaged.pdb").resolve().as_posix()

        if membrane:
            membrane_mod.check_single_barostat(mm.system)
            try:
                result["membrane_production"] = membrane_mod.membrane_report(
                    Path("NVTtrajectory_lastframe.pdb").resolve().as_posix(), "NVTtrajectory.dcd",
                    csvfile="membrane_box_production.csv")
            except Exception as exc:
                print(f"[WARNING] Could not report the membrane box size: {exc}")

    final = result.get("lastframe") or result.get("gentle") or result.get("minimized")
    if final and Path(final).exists():
        shell_after = coordination_shell(final, preprocess_inputs)
        report_coordination(shell_after, "after MD")
        result["coordination_end"] = shell_after
        changes = compare_coordination(shell_before, shell_after)
        if changes:
            print()
            print("The coordination sphere changed during MD:")
            for line in changes:
                print(line)
            print("Check these are physical before building on this structure: a "
                  "contact that appeared was not constrained, and a non-bonded "
                  "metal has no Pauli repulsion to stop a donor collapsing onto it.")
        else:
            print()
            print("The coordination sphere is unchanged.")

    print("[DONE] Workflow complete")
    save_result(result,  filename=args.config, section="md")
    return result
