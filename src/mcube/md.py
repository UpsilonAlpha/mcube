#!/usr/bin/env python3

from __future__ import annotations
import argparse
import os
import json
from pathlib import Path
from typing import Any, Dict

from openmm import XmlSerializer, unit
from openmm.app import PDBFile, ForceField, CutoffNonPeriodic, HBonds

from ash import *

from .preprocess import (
    save_result,
    load_section,
)

# ==============================
# Unified run() function
# ==============================

def run(args) -> Dict[str, Any]:
    result: Dict[str, Any] = {}

    preprocess = load_section(args.config, "preprocess")
    md = load_section(args.config, "md")
    inputs = md.get("inputs")
    if inputs is None:
        raise KeyError(f"No md.inputs section in {args.config}")

    constraints = inputs.get("constraints")
    if constraints is None:
        constraints = preprocess["results"]["constraints"]

    pdb_file = inputs.get("pdb_file")
    if pdb_file is None:
        pdb_file = preprocess["results"]["processed_pdb"]
    
    ligand_forcefields = inputs.get("ligand_xml")
    if ligand_forcefields is None:
        ligand_forcefields = preprocess["results"]["ligand_xml"]
    
    default_forcefields = inputs.get("default_forcefields")
    if default_forcefields is None:
        default_forcefields = preprocess["inputs"]["default_forcefields"]
    default_forcefields = list(default_forcefields)


    if ligand_forcefields:
        default_forcefields.append(ligand_forcefields)

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
            numcores=inputs.get("cores"),
            autoconstraints='HBonds',
            constraints=constraints,
            nonbondedMethod_noPBC='CutoffNonPeriodic',
            nonbonded_cutoff_noPBC=20,
        )

    else:
        mm = OpenMMTheory(
            xmlfiles=default_forcefields,
            pdbfile=pdb_file,
            periodic=True,
            numcores=inputs.get("cores"),
            autoconstraints='HBonds',
            constraints=constraints,
            rigidwater=True,
        )

    fragment = Fragment(pdbfile=pdb_file)
    result["forcefields"] = default_forcefields

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

    if not inputs.get("skip_npt") and not inputs.get("implicit"):
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

    print("[DONE] Workflow complete")
    save_result(result,  filename=args.config, section="md")
    return result
