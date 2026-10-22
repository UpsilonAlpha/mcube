#!/usr/bin/env python3

from __future__ import annotations
import argparse
import os
import json
from pathlib import Path
from typing import Any, Dict, Optional

import openmm
from openmm import XmlSerializer, unit
from openmm.app import PDBFile, ForceField, CutoffNonPeriodic, HBonds

from ash import *

from .preprocess import (
    COORDINATION_FORCE,
    compare_coordination,
    coordination_shell,
    pair_distances,
    report_coordination,
    resolve_cores,
    to_selectors,
    save_result,
    load_section,
    resolve_atoms,
)
from . import ions
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



def apply_metal_parameters(mm, pdb_file: str, metals, which: str) -> list:
    """Override the metals' Lennard-Jones terms with a published ion set.

    Amber's stock parameters put the metal-oxygen pair minimum too close - about
    1.69 A for Zn(II) against a real bond of 2.0-2.1 A - so the metal pulls
    donors in harder than it should. See ions.py for the set and its sources.

    The ion is identified from the structure's element and the charge the force
    field itself assigned, so nothing has to be specified: a metal in the PDB is
    looked up and replaced. An element with no entry is left alone and said so.
    """
    nonbonded = next((f for f in mm.system.getForces()
                      if isinstance(f, openmm.NonbondedForce)), None)
    if nonbonded is None:
        print("[WARNING] no NonbondedForce: metal parameters unchanged")
        return []

    wanted = {m.upper() for m in (metals or [])}
    atoms = list(PDBFile(pdb_file).topology.atoms())
    changed = []
    for atom in atoms:
        element = (atom.element.symbol if atom.element else "").upper()
        if element not in wanted:
            continue
        charge, sigma, epsilon = nonbonded.getParticleParameters(atom.index)
        entry = ions.lookup(which, element, charge.value_in_unit(unit.elementary_charge))
        if entry is None:
            print(f"[WARNING] no {which} parameters for {element} with charge "
                  f"{charge}: left at the force field's own values")
            continue
        new_sigma, new_epsilon = entry
        nonbonded.setParticleParameters(atom.index, charge,
                                        new_sigma * unit.nanometer,
                                        new_epsilon * unit.kilojoule_per_mole)
        changed.append(atom.index)
        print(f"[INFO] {which.upper()} parameters for {element} "
              f"{atom.residue.chain.id}:{atom.residue.id}: sigma "
              f"{sigma.value_in_unit(unit.nanometer):.4f} -> {new_sigma:.4f} nm, "
              f"epsilon {epsilon.value_in_unit(unit.kilojoule_per_mole):.4f} -> "
              f"{new_epsilon:.4f} kJ/mol")
    return changed


def build_md_theory(inputs: Dict[str, Any], preprocess: Dict[str, Any], pdb_file: str,
                    restraints, cores=None):
    """
    The classical system used by the md stage: force fields, HBond constraints,
    rigid water, 1.5 amu hydrogens and the harmonic distance restraints that
    hold the metal site. Returns (OpenMMTheory, force-field list).
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
            restraints=restraints,
            rigidwater=True,
        )
    return mm, default_forcefields


def md_inputs(config: str, section: Optional[str] = None) -> Dict[str, Any]:
    """`md.inputs`, overridden by `section`'s inputs.

    Layered so a second MD run - releasing a restraint, extending the time,
    starting from the first's last frame - inherits the first's whole setup and
    states only what differs. Both runs then read one
    `preprocess.results.coordination`, instead of a duplicated config directory
    whose copy can drift from the original - which is how a distance one run
    imposed came to be read as detected coordination by the next.
    """
    inputs: Dict[str, Any] = dict(load_section(config, "md").get("inputs") or {})
    if section and section != "md":
        inputs.update(load_section(config, section).get("inputs") or {})
    return inputs


def md_restraints(inputs: Dict[str, Any], preprocess: Dict[str, Any], pdb_file: str):
    """Harmonic distance restraints: the detected metal coordination, plus any asked for.

    The coordination a force field has no bonded parameters for is held by
    restraint rather than rigid constraint. ASH, whose workflow this pipeline
    follows, sanctions the substitution directly - "we could also have added
    harmonic bond restraints instead of rigid constraints" - and it avoids both
    ways a constrained metal site fails: many constraints sharing one centre
    over-determine the constraint network, and OpenMM refuses a constraint that
    touches a massless atom. A restraint also carries no element-specific
    parameters, so it works for any metal without parameterisation.

    `preprocess.results.coordination` already carries the distance and force
    constant for every contact it detected, printed into the config, so editing
    a row there is what changes the run. A row given as a bare `[i, j]` pair is
    restrained at whatever distance this run starts from, at the default
    constant.

    Returns (restraints, detected) with `detected` the auto-generated subset, so
    the results can record which restraints were a measurement and which were a
    choice. Only a distance someone chose can launder an assumption into a later
    run: holding a contact where it already was asserts nothing new.
    """
    asked = [list(r) for r in (resolve_atoms(inputs.get("restraints") or None, pdb_file) or [])]
    rows = resolve_atoms((preprocess.get("results") or {}).get("coordination"), pdb_file) or []
    # An explicit restraint on a pair wins: the user named a distance for it.
    named = {frozenset(r[:2]) for r in asked}
    rows = [list(r) for r in rows if frozenset(r[:2]) not in named]
    bare = [r for r in rows if len(r) < 3]
    measured = dict(zip((frozenset(r[:2]) for r in bare), pair_distances(pdb_file, bare)))
    detected = [r if len(r) >= 4 else
                [r[0], r[1], r[2] if len(r) > 2 else measured[frozenset(r[:2])],
                 COORDINATION_FORCE]
                for r in rows]
    return (asked + detected) or None, detected


def run(args) -> Dict[str, Any]:
    result: Dict[str, Any] = {}

    preprocess = load_section(args.config, "preprocess")
    section = getattr(args, "section", None) or "md"
    inputs = md_inputs(args.config, section)
    if not inputs:
        raise KeyError(f"No md.inputs section in {args.config}")
    if section != "md":
        print(f"[INFO] md.inputs overridden by {section}.inputs; results -> {section}")

    pdb_file = inputs.get("pdb_file")
    if pdb_file is None:
        pdb_file = preprocess["results"]["processed_pdb"]

    membrane = md_membrane(inputs, preprocess)
    if membrane and inputs.get("implicit"):
        raise ValueError("`membrane` and `implicit` cannot both be set in the md inputs: a bilayer "
                         "needs an explicit, periodic box.")

    restraints, detected = md_restraints(inputs, preprocess, pdb_file)
    if restraints:
        print(f"[INFO] {len(restraints)} harmonic distance restraints "
              f"({len(detected)} detected metal coordination): {restraints}")

    # The coordination sphere this run starts from. Reported again at the end,
    # because a contact that was not there at setup was never restrained, and a
    # non-bonded metal will happily pull one in.
    preprocess_inputs = preprocess.get("inputs") or {}
    shell_before = coordination_shell(pdb_file, preprocess_inputs)
    report_coordination(shell_before, "before MD")
    result["coordination_start"] = shell_before

    mm, default_forcefields = build_md_theory(inputs, preprocess, pdb_file, restraints,
                                              cores=resolve_cores(inputs.get("cores")))

    # IOD by default: Amber's stock metal radii put the metal-oxygen minimum too
    # close (~1.69 A for Zn(II) against 2.0-2.1 A measured). `none` keeps whatever
    # the force field supplies. Recorded in the results, so a run says which it used.
    metal_set = inputs.get("metal_parameters", "iod")
    if metal_set and str(metal_set).lower() != "none":
        apply_metal_parameters(mm, pdb_file, preprocess_inputs.get("metals"), metal_set)
        result["metal_parameters"] = str(metal_set).lower()

    fragment = Fragment(pdbfile=pdb_file)
    result["forcefields"] = default_forcefields
    # Recorded as selectors, not the indices OpenMM was given: indices only mean
    # something against the structure they were resolved on. The two lists are
    # kept apart because only a distance this run chose can turn into apparent
    # structure in a later one - a detected pair restrained at the distance it
    # already had asserts nothing new. `check_inherited_restraints` enforces it.
    chosen = [r for r in (restraints or []) if r not in detected]
    result["restraints"] = to_selectors(chosen, pdb_file)
    result["coordination_restraints"] = to_selectors(detected, pdb_file)

    print("[STEP] Minimization...")
    OpenMM_Opt(fragment=fragment, theory=mm, maxiter=500, tolerance=1)
    result["minimized"] = Path("frag-minimized.pdb").resolve().as_posix()
    save_result(result, filename=args.config, section=section)

    if inputs.get("gentle"):
        Gentle_warm_up_MD(theory=mm, fragment=fragment, time_steps=[0.0005,0.001,0.004],
                    steps=[10,50,10000], temperatures=[1,10,300])
        result["gentle"] = Path("warmup_MD_cycle2_lastframe.pdb").resolve().as_posix()
        fragment = Fragment(pdbfile=result["gentle"])
        save_result(result, filename=args.config, section=section)

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
        save_result(result, filename=args.config, section=section)

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
    save_result(result,  filename=args.config, section=section)
    return result
