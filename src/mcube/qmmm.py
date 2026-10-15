#!/usr/bin/env python3
"""
QM/MM workflow (xTB) with automatic QM region selection.

Pipeline:
1. Take the last MD frame (or a user-supplied PDB)
2. Select QM atoms (metal coordination, from preprocess)
3. Build QM/MM system
4. Run short QM/MM MD
5. QM/MM geometry optimisation (optionally with an active region and with
   user bond constraints, `opt_constraints`)
6. Optional PySCF QM/MM single point with density, MEP and HOMO/LUMO cube files
"""

from __future__ import annotations

import copy
from typing import Dict, Any, List, Optional
from pathlib import Path

import numpy as np

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


# ==============================
# Helpers shared by the QM/MM stages (qmmm, scan, neb, ts)
# ==============================

def stage_inputs(config: str, base: str, section: Optional[str] = None) -> Dict[str, Any]:
    """Inputs of the `base` section, overridden by those of `section` if it has any
    (e.g. scan_attack.inputs over scan.inputs)."""
    inputs = dict(load_section(config, base).get("inputs") or {})
    if section and section != base:
        inputs.update(load_section(config, section).get("inputs") or {})
    return inputs


def mm_settings(qmmm_inputs: Dict[str, Any], md: Dict[str, Any]):
    """Force-field files and implicit-solvent flag: qmmm.inputs, else what md used."""
    forcefields = qmmm_inputs.get("forcefields")
    if forcefields is None:
        forcefields = (md.get("results") or {}).get("forcefields")
    implicit = qmmm_inputs.get("implicit", (md.get("inputs") or {}).get("implicit", False))
    return forcefields, implicit


def make_mm(forcefields, pdb_file: str, implicit: bool, cores: int):
    return OpenMMTheory(
        xmlfiles=forcefields,
        pdbfile=pdb_file,
        autoconstraints=None,
        periodic=not implicit,
        rigidwater=False,
        numcores=cores,  # OpenMM CPU threads; ASH defaults to 1
    )


def make_qmmm(qmmm_inputs: Dict[str, Any], mm, fragment, qm_atoms, boundary_excluded,
              cores: int, runmode: Optional[str] = None):
    xtb = xTBTheory(xtbmethod=qmmm_inputs.get("xtb_method", "GFN2"), numcores=cores,
                    runmode=runmode or qmmm_inputs.get("xtb_runmode", "inputfile"))
    return QMMMTheory(
        qm_theory=xtb,
        mm_theory=mm,
        fragment=fragment,
        qmatoms=qm_atoms,
        excludeboundaryatomlist=boundary_excluded,
        embedding="electrostatic",
        printlevel=1,
    )


def resolve_actatoms(actatoms_file, inputs: Dict[str, Any], mm, fragment, qm_atoms, boundary_excluded) -> List[int]:
    """Active region: an active_atoms file (e.g. the reactant's), else all atoms within
    `actradius` of `actregion_origin`, else just the QM atoms."""
    if actatoms_file:
        return read_intlist_from_file(actatoms_file)
    if inputs.get("actradius"):
        originatom = inputs.get("actregion_origin")
        if originatom is None:
            if not boundary_excluded:
                raise ValueError("actradius is set but there is no actregion_origin and no "
                                 "non-standard QM residue to centre the active region on")
            originatom = boundary_excluded[0]
        actregiondefine(mmtheory=mm, fragment=fragment, radius=inputs.get("actradius"), originatom=originatom)
        return read_intlist_from_file("active_atoms")
    return list(qm_atoms)


def distance(coords, i: int, j: int) -> float:
    return float(np.linalg.norm(np.asarray(coords[i]) - np.asarray(coords[j])))


def bond_constraints(pairs, coords) -> List[list]:
    """[[i, j]] or [[i, j, d]] -> [[i, j, d]]; a missing d is the current distance (Å)."""
    out = []
    for p in pairs or []:
        i, j = int(p[0]), int(p[1])
        d = float(p[2]) if len(p) > 2 and p[2] is not None else round(distance(coords, i, j), 4)
        out.append([i, j, d])
    return out


def check_in_active(pairs, actatoms, what: str = "constraint"):
    """geomeTRIC only sees the active region, so constrained atoms must be in it."""
    act = set(int(a) for a in actatoms)
    missing = sorted({int(a) for p in pairs for a in p[:2]} - act)
    if missing:
        raise ValueError(f"{what} atoms {missing} are not in the active region")


def water_constraints(mm, actatoms) -> List[list]:
    """O-H and H-H pairs of active-region TIP3P waters (kept rigid during optimisation)."""
    return [[int(p[0]), int(p[1])] for p in
            getwaterconstraintslist(openmmtheoryobject=mm, atomlist=actatoms, watermodel='tip3p')]


def optimizer_constraints(water_pairs, user_pairs, coords):
    """
    Constraints dict and constrainvalue flag for ASH's Optimizer.

    Without user constraints this is the original behaviour: the waters' bonds are
    frozen ($freeze). With user constraints geomeTRIC's $set block needs a value for
    every bond, so the waters are held at their current lengths, and so are user pairs
    given without a value. A fresh list is returned each time because ASH converts the
    indices to active-region indices in place.
    """
    if not user_pairs:
        return {'bond': copy.deepcopy(water_pairs)}, False
    bonds = bond_constraints(water_pairs, coords) + bond_constraints(user_pairs, coords)
    return {'bond': bonds}, True


def clean_geometric_files():
    """geomeTRIC appends to its trajectory files, so remove any left from earlier runs."""
    for f in ("geometric_OPTtraj_Full.xyz", "geometric_OPTtraj_QMregion.xyz"):
        Path(f).unlink(missing_ok=True)


def write_pdb(mm, coords, filename: str) -> str:
    """Full system with the MM topology (keeps residue names/ids, e.g. for NEB endpoints)."""
    with open(filename, "w") as f:
        PDBFile.writeFile(mm.topology, np.asarray(coords) * unit.angstrom, f, keepIds=True)
    return Path(filename).resolve().as_posix()



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

    forcefields, implicit = mm_settings(inputs, md)

    print("[STEP] Building MM system...")
    mm = make_mm(forcefields, pdb_file, implicit, inputs.get("cores", 1))

    fragment = Fragment(pdbfile=pdb_file)

    print("[STEP] Building QM region...")
    qm_atoms, boundary_excluded = resolve_qm_region(inputs, preprocess, pdb_file)
    result: Dict[str, Any] = {
        "qm_atoms": qm_atoms,
        "boundary_excluded": boundary_excluded,
    }
    save_result(result, filename=args.config, section=section)

    qmmm_theory = make_qmmm(inputs, mm, fragment, qm_atoms, boundary_excluded, inputs.get("cores"))

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

        # An active_atoms file reuses another run's active region (e.g. the reactant's for the product)
        actatoms = resolve_actatoms(getattr(args, "actatoms", None), inputs, mm, fragment,
                                    qm_atoms, boundary_excluded)

        # Optional bond constraints: [[i, j]] holds the current distance, [[i, j, d_Å]] sets it
        user_pairs = inputs.get("opt_constraints") or []
        check_in_active(user_pairs, actatoms, "opt_constraints")
        constraints, constrainvalue = optimizer_constraints(water_constraints(mm, actatoms), user_pairs,
                                                            fragment.coords)
        if user_pairs:
            result["opt_constraints"] = bond_constraints(user_pairs, fragment.coords)
            print("[INFO] Bond constraints (i, j, Å): " + ", ".join(
                f"{i}-{j} {d:.3f}" for i, j, d in result["opt_constraints"]))

        clean_geometric_files()

        print("[STEP] Running QM/MM Optimisation...")
        opt = Optimizer(fragment=fragment, theory=qmmm_theory, ActiveRegion=True, actatoms=actatoms,
                        maxiter=inputs.get("opt_maxiter", 200), constraints=constraints,
                        constrainvalue=constrainvalue, charge=inputs.get("charge"), mult=inputs.get("mult"))

        # Full optimised system with the MM topology, e.g. as an NEB endpoint
        result["opttraj"] = Path("geometric_OPTtraj_Full.xyz").resolve().as_posix()
        result["optimized_pdb"] = write_pdb(mm, fragment.coords, "QM_MM_optimized.pdb")
        if getattr(opt, "energy", None) is not None:
            result["opt_energy"] = float(opt.energy)
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
