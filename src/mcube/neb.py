from __future__ import annotations

import argparse
import os
from typing import List, Set, Dict, Any, Tuple
from openmm.app import PDBFile

from .preprocess import (
    save_result,
    load_section,
)

from .qmmm import build_qm_region
from ash import *


def run(args) -> Dict[str, Any]:
    result: Dict[str, Any] = {}

    qmmm = load_section(args.input, "qmmm")
    md = load_section(args.input, "md")
    preprocess = load_section(args.input, "preprocess")

    qm_residues = preprocess["results"]["qm_residues"]

    numimages=8

    print("[STEP] Building MM system...")
    mm = OpenMMTheory(
        xmlfiles=qmmm["inputs"]["forcefields"],
        pdbfile=args.reactant,
        autoconstraints=None,
        periodic=True,
        rigidwater=False,
    )

    print("[STEP] Building QM region...")
    qm_atoms, boundary_excluded = build_qm_region(args.reactant, qm_residues)
    result["qm_atoms"] = qm_atoms
    result["boundary_excluded"] = boundary_excluded

    #Theory to use for NEB
    xtb = xTBTheory(xtbmethod="GFN2", numcores=args.cores, runmode='library')

    ################################################
    # Defining reactant and product ASH fragments
    #################################################
    react=Fragment(pdbfile=args.reactant, charge=qmmm["inputs"]["charge"], mult=qmmm["inputs"]["mult"])
    prod=Fragment(pdbfile=args.product, charge=qmmm["inputs"]["charge"], mult=qmmm["inputs"]["mult"])

    qmmm = QMMMTheory(
        qm_theory=xtb,
        mm_theory=mm,
        fragment=react,
        qmatoms=qm_atoms,
        excludeboundaryatomlist=boundary_excluded,
        embedding="electrostatic",
        printlevel=1,
    )


    print("Reactant atoms:", react.numatoms)
    print("Product atoms: ", prod.numatoms)
    print("Reactant attributes:", [a for a in dir(react) if "coord" in a.lower()])
    print("Product attributes: ", [a for a in dir(prod) if "coord" in a.lower()])
    
    #Run NEB to find saddlepoint. Returns an ASH Results object
    NEB_result = NEB(reactant=react, product=prod, theory=qmmm, images=numimages, runmode='parallel', numcores=args.cores)
    print(NEB_result)

    #Optional NumFreq job on saddlepoint to confirm that a saddlepoint was found.
    NumFreq(theory=xtbcalc, fragment=NEB_result.saddlepoint_fragment)
