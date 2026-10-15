#!/usr/bin/env python3
"""
QM/MM workflow (xTB) with automatic QM region selection.

Pipeline:
1. Take the last MD frame (or a user-supplied PDB)
2. Select QM atoms (metal coordination, from preprocess)
3. Build QM/MM system
4. Run short QM/MM MD
5. QM/MM geometry optimisation (optionally with an active region)
6. Optional PySCF QM/MM single point with density, MEP and HOMO/LUMO cube files
"""

from __future__ import annotations

from typing import Dict, Any
from pathlib import Path

from openmm import unit
from openmm.app import PDBFile

from .preprocess import (
    save_result,
    load_section,
    build_qm_region,
    qm_region_settings,
)

from ash import *


def resolve_qm_region(inputs: Dict[str, Any], preprocess: Dict[str, Any], pdb_file: str):
    """
    QM atoms and boundary-excluded atoms, in order of preference:
    1. `qm_residues` given in this stage's inputs (rebuilt on pdb_file)
    2. `qm_atoms` / `boundary_excluded` indices saved by preprocess
    3. `qm_residues` saved by preprocess (configs from before indices were saved)
    """
    settings = qm_region_settings(preprocess.get("inputs") or {})
    results = preprocess.get("results") or {}

    if inputs.get("qm_residues") is not None:
        return build_qm_region(pdb_file, inputs["qm_residues"], **settings)

    if results.get("qm_atoms") is not None:
        print(f"[INFO] Using {len(results['qm_atoms'])} QM atoms saved by preprocess")
        return results["qm_atoms"], results.get("boundary_excluded", [])

    return build_qm_region(pdb_file, results["qm_residues"], **settings)


def run(args) -> Dict[str, Any]:

    preprocess = load_section(args.config, "preprocess")
    md = load_section(args.config, "md")
    qmmm = load_section(args.config, "qmmm")
    inputs = qmmm.get("inputs")
    if inputs is None:
        raise KeyError(f"No qmmm.inputs section in {args.config}")
    # Results go to this section, so several structures (e.g. NEB endpoints) can share one config
    section = getattr(args, "section", None) or "qmmm"
    md_results = md.get("results") or {}

    pdb_file = getattr(args, "pdb", None) or inputs.get("pdb_file") or md_results.get("lastframe")

    forcefields = inputs.get("forcefields")
    if forcefields is None:
        forcefields = md_results.get("forcefields")

    implicit = inputs.get("implicit", (md.get("inputs") or {}).get("implicit", False))

    print("[STEP] Building MM system...")
    mm = OpenMMTheory(
        xmlfiles=forcefields,
        pdbfile=pdb_file,
        autoconstraints=None,
        periodic=not implicit,
        rigidwater=False,
    )

    fragment = Fragment(pdbfile=pdb_file)

    print("[STEP] Building QM region...")
    qm_atoms, boundary_excluded = resolve_qm_region(inputs, preprocess, pdb_file)
    result: Dict[str, Any] = {
        "qm_atoms": qm_atoms,
        "boundary_excluded": boundary_excluded,
    }
    save_result(result, filename=args.config, section=section)

    xtb = xTBTheory(xtbmethod=inputs.get("xtb_method", "GFN2"), numcores=inputs.get("cores"),
                    runmode=inputs.get("xtb_runmode", "inputfile"))

    qmmm_theory = QMMMTheory(
        qm_theory=xtb,
        mm_theory=mm,
        fragment=fragment,
        qmatoms=qm_atoms,
        excludeboundaryatomlist=boundary_excluded,
        embedding="electrostatic",
        printlevel=1,
    )


    if not inputs.get("skip_equilibrate"):
        print("[STEP] Running QM/MM MD...")
        OpenMM_MD(
            fragment=fragment,
            theory=qmmm_theory,
            timestep=inputs.get("timestep"),
            simulation_time=inputs.get("time"),
            traj_frequency=inputs.get("traj_frequency", 1000),
            temperature=inputs.get("temp"),
            integrator='LangevinMiddleIntegrator',
            trajfilename="QM_MM",
            coupling_frequency=1,
            charge=inputs.get("charge"),
            mult=inputs.get("mult"),
        )

        result["trajfile"] = Path("QM_MM.dcd").resolve().as_posix()
        result["lastframe"] = Path("QM_MM_lastframe.pdb").resolve().as_posix()
        save_result(result, filename=args.config, section=section)

    if not inputs.get("skip_optimize"):

        # Start from the QM/MM MD end point if there is one, else from pdb_file
        if Path("QM_MM_lastframe.pdb").exists():
            fragment = Fragment(pdbfile="QM_MM_lastframe.pdb")

        if getattr(args, "actatoms", None):
            # Reuse another run's active region (e.g. the reactant's for the product)
            actatoms = read_intlist_from_file(args.actatoms)
        elif inputs.get("actradius"):
            originatom = inputs.get("actregion_origin")
            if originatom is None:
                if not boundary_excluded:
                    raise ValueError("actradius is set but there is no actregion_origin and no "
                                     "non-standard QM residue to centre the active region on")
                originatom = boundary_excluded[0]
            actregiondefine(mmtheory=mm, fragment=fragment, radius=inputs.get("actradius"), originatom=originatom)
            actatoms = read_intlist_from_file("active_atoms")
        else:
            actatoms = qm_atoms

        waterconlist = getwaterconstraintslist(openmmtheoryobject=mm, atomlist=actatoms, watermodel='tip3p')
        waterconstraints = {'bond': waterconlist}

        # geomeTRIC appends to its trajectory files, so remove any left from earlier runs
        for f in ("geometric_OPTtraj_Full.xyz", "geometric_OPTtraj_QMregion.xyz"):
            Path(f).unlink(missing_ok=True)

        print("[STEP] Running QM/MM Optimisation...")
        Optimizer(fragment=fragment, theory=qmmm_theory, ActiveRegion=True, actatoms=actatoms, maxiter=200,
            constraints=waterconstraints, charge=inputs.get("charge") , mult=inputs.get("mult") )

        # Full optimised system with the MM topology, e.g. as an NEB endpoint
        with open("QM_MM_optimized.pdb", "w") as f:
            PDBFile.writeFile(mm.topology, fragment.coords * unit.angstrom, f, keepIds=True)
        result["opttraj"] = Path("geometric_OPTtraj_Full.xyz").resolve().as_posix()
        result["optimized_pdb"] = Path("QM_MM_optimized.pdb").resolve().as_posix()
        save_result(result, filename=args.config, section=section)

    if inputs.get("singlepoint", True):
        # QM/MM single point on the current (optimised) geometry: CA-CB cuts get link atoms and
        # the MM environment enters as point charges. No continuum solvation: ASH's PySCF
        # interface ignores it when point charges are present.
        mult = inputs.get("mult")
        pyscf = PySCFTheory(scf_type="RKS" if mult == 1 else "UKS",
                            functional=inputs.get("sp_functional", "wb97x-v"),
                            basis=inputs.get("sp_basis", "def2-svp"),
                            numcores=inputs.get("cores"))
        qmmm_sp = QMMMTheory(
            qm_theory=pyscf,
            mm_theory=mm,
            fragment=fragment,
            qmatoms=qm_atoms,
            excludeboundaryatomlist=boundary_excluded,
            embedding="electrostatic",
            printlevel=1,
        )
        print("[STEP] Running PySCF QM/MM single point...")
        Singlepoint(theory=qmmm_sp, fragment=fragment, charge=inputs.get("charge"), mult=mult)

        mol = pyscf.mol
        mf = pyscf.mf
        dm = mf.make_rdm1()
        if mult != 1:
            dm = dm[0] + dm[1]  # total density from alpha + beta
        pyscf.cubegen_density(mol, "density.cube", dm)
        pyscf.cubegen_mep(mol, "mep.cube", dm)

        # Alpha orbitals for UKS
        mo_coeff = mf.mo_coeff if mult == 1 else mf.mo_coeff[0]
        homo = mol.nelec[0] - 1
        lumo = mol.nelec[0]
        pyscf.cubegen_orbital(mol, "HOMO.cube", mo_coeff[:, homo])
        pyscf.cubegen_orbital(mol, "LUMO.cube", mo_coeff[:, lumo])

        result["cubes"] = [Path(f).resolve().as_posix() for f in
                           ("density.cube", "mep.cube", "HOMO.cube", "LUMO.cube")]
        save_result(result, filename=args.config, section=section)

    return result
