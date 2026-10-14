#!/usr/bin/env python3
"""
QM/MM workflow (xTB) with automatic QM region selection.

Pipeline:
1. Generate onechain PDB from trajectory
2. Select QM residues (metal coordination)
3. Build QM/MM system
4. Run short QM/MM MD
"""

from __future__ import annotations

import argparse
import os
from typing import List, Set, Dict, Any, Tuple
from openmm.app import PDBFile
from pathlib import Path
from pyscf import gto, dft


from preprocess import (
    save_result,
    load_section,
)

from ash import *


# ==============================
# Constants
# ==============================

BACKBONE_ATOMS = {
    "N", "H", "H1", "H2", "H3",
    "CA", "HA", "HA2", "HA3",
    "C", "O", "OXT",
}

AMINO_ACIDS = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
    "HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
    "THR","TRP","TYR","VAL",
}


# ==============================
# QM Region Builder
# ==============================
def build_qm_region(
    pdb_file: str,
    qm_residue_dict: Dict[str, Dict[int, str]],
    backbone_atoms=BACKBONE_ATOMS,
    amino_acids=AMINO_ACIDS,
):
    pdb = PDBFile(pdb_file)
    topology = pdb.topology

    qm_residue_dict = {
        chain: {str(k): v for k, v in residues.items()}
        for chain, residues in qm_residue_dict.items()
    }

    qm_atoms = set()
    boundary_excluded = set()

    print("\n===== QM Region Construction =====")

    for res in topology.residues():
        chain = res.chain.id
        resid = res.id.strip()   # ← FIXED

        if chain in qm_residue_dict and resid in qm_residue_dict[chain]:
            print(f"[QM] Including residue {res.name} {chain}:{resid}")

            for atom in res.atoms():
                if atom.name not in backbone_atoms:
                    qm_atoms.add(atom.index)

                if res.name not in amino_acids:
                    qm_atoms.add(atom.index)
                    boundary_excluded.add(atom.index)

    print(f"\n[INFO] QM atoms: {len(qm_atoms)}")
    print(f"[INFO] Boundary excluded atoms: {len(boundary_excluded)}")

    return sorted(qm_atoms), sorted(boundary_excluded)



def run(args) -> Dict[str, Any]:

    preprocess = load_section(args.config, "preprocess")
    md = load_section(args.config, "md")
    qmmm = load_section(args.config, "qmmm")
    inputs = qmmm.get("inputs")

    pdb_file = inputs.get("pdb_file")
    if pdb_file is None:
        pdb_file = md["results"].get("lastframe")

    forcefields = inputs.get("forcefields")
    if forcefields is None:
        forcefields = md["results"].get("forcefields")

    qm_residues = inputs.get("qm_residues")
    if qm_residues is None:
        qm_residues = preprocess["results"].get("qm_residues")

    print(qm_residues)

    print("[STEP] Building MM system...")
    mm = OpenMMTheory(
        xmlfiles=forcefields,
        pdbfile=pdb_file,
        autoconstraints=None,
        periodic=True,
        rigidwater=False,
    )

    fragment = Fragment(pdbfile=pdb_file)

    print("[STEP] Building QM region...")
    qm_atoms, boundary_excluded = build_qm_region(pdb_file, qm_residues)
    save_result(
        {
            "qm_atoms": qm_atoms,
            "boundary_excluded": boundary_excluded,
        },
        filename=args.config,
        section="qmmm"
    )

    xtb = xTBTheory(xtbmethod="GFN2", numcores=inputs.get("cores"))

    qmmm = QMMMTheory(
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
            theory=qmmm,
            timestep=inputs.get("timestep"),
            simulation_time=inputs.get("time"),
            traj_frequency=inputs.get("time"),
            temperature=inputs.get("temp"),
            integrator='LangevinMiddleIntegrator',
            trajfilename="QM_MM",
            coupling_frequency=1,
            charge=inputs.get("charge"),
            mult=inputs.get("mult"),
        )

        save_result(
            {
                "trajfile": Path("QM_MM.dcd").resolve().as_posix(),
                "lastframe": Path("QM_MM_lastframe.pdb").resolve().as_posix(),
            },
            filename=args.config,
            section="qmmm"
        )

    if not inputs.get("skip_optimize"):
        
        fragment = Fragment(pdbfile="QM_MM_lastframe.pdb")

        if inputs.get("actradius") != 0:
            actregiondefine(mmtheory=mm, fragment=fragment, radius=inputs.get("actradius") , originatom=boundary_excluded[0])
            actatoms = read_intlist_from_file("active_atoms")
        else:
            actatoms=qm_atoms

        waterconlist = getwaterconstraintslist(openmmtheoryobject=mm, atomlist=actatoms, watermodel='tip3p')
        waterconstraints = {'bond': waterconlist}

        print("[STEP] Running QM/MM Optimisation...")
        Optimizer(fragment=fragment, theory=qmmm, ActiveRegion=True, actatoms=actatoms, maxiter=200,
            constraints=waterconstraints, charge=inputs.get("charge") , mult=inputs.get("mult") )

        save_result(
            {"opttraj": Path("geometric_OPTtraj_Full.xyz").resolve().as_posix()},
            filename=args.config,
            section="qmmm"
        )

    with open("geometric_OPTtraj_QMregion.xyz") as f:
        lines = f.readlines()

    natoms = int(lines[0].strip())
    block_size = natoms + 2

    last_frame = lines[-block_size:]

    with open("OPTtraj_lastframe.xyz", "w") as f:
        f.writelines(last_frame)

    pyscf = PySCFTheory(scf_type="RKS", functional="wb97x-v", basis="def2-svp", solvation="ddCOSMO", solvation_eps=78)
    qm_region = Fragment(xyzfile="OPTtraj_lastframe.xyz")
    Singlepoint(theory=pyscf, fragment=qm_region, charge=inputs.get("charge"), mult=inputs.get("mult"))

    mol = pyscf.mol
    mf = pyscf.mf
    dm = mf.make_rdm1()
    pyscf.cubegen_density(mol, "density.cube", dm)
    pyscf.cubegen_mep(mol, "mep.cube", dm)

    homo = mf.nelec[0] - 1
    lumo = mf.nelec[0]
    pyscf.cubegen_orbital(mol, "HOMO.cube", mf.mo_coeff[:, homo])

    pyscf.cubegen_orbital(mol, "LUMO.cube", mf.mo_coeff[:, lumo])

    return


