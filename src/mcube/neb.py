from __future__ import annotations

from typing import Dict, Any

from .preprocess import (
    save_result,
    load_section,
)

from .qmmm import resolve_qm_region
from ash import *


def run(args) -> Dict[str, Any]:
    result: Dict[str, Any] = {}

    qmmm_cfg = load_section(args.config, "qmmm")
    md = load_section(args.config, "md")
    preprocess = load_section(args.config, "preprocess")
    qmmm_inputs = qmmm_cfg.get("inputs") or {}
    neb_inputs = load_section(args.config, "neb").get("inputs") or {}

    charge = qmmm_inputs["charge"]
    mult = qmmm_inputs["mult"]
    numimages = neb_inputs.get("images", 8)

    # ASH runs up to `numcores` images in parallel, each with the theory's own numcores,
    # so split the cores instead of giving both the full count
    neb_workers = max(1, min(numimages, args.cores))
    theory_cores = max(1, args.cores // neb_workers)
    runmode = 'parallel' if neb_workers > 1 else 'serial'

    forcefields = qmmm_inputs.get("forcefields")
    if forcefields is None:
        forcefields = md["results"].get("forcefields")

    implicit = qmmm_inputs.get("implicit", (md.get("inputs") or {}).get("implicit", False))

    print("[STEP] Building MM system...")
    mm = OpenMMTheory(
        xmlfiles=forcefields,
        pdbfile=args.reactant,
        autoconstraints=None,
        periodic=not implicit,
        rigidwater=False,
    )

    print("[STEP] Building QM region...")
    qm_atoms, boundary_excluded = resolve_qm_region(qmmm_inputs, preprocess, args.reactant)
    result["qm_atoms"] = qm_atoms
    result["boundary_excluded"] = boundary_excluded

    #Theory to use for NEB
    xtb = xTBTheory(xtbmethod=qmmm_inputs.get("xtb_method", "GFN2"), numcores=theory_cores, runmode='library')

    ################################################
    # Defining reactant and product ASH fragments
    #################################################
    react=Fragment(pdbfile=args.reactant, charge=charge, mult=mult)
    prod=Fragment(pdbfile=args.product, charge=charge, mult=mult)

    if react.numatoms != prod.numatoms:
        raise ValueError(f"Reactant ({react.numatoms} atoms) and product ({prod.numatoms} atoms) differ")

    qmmm_theory = QMMMTheory(
        qm_theory=xtb,
        mm_theory=mm,
        fragment=react,
        qmatoms=qm_atoms,
        excludeboundaryatomlist=boundary_excluded,
        embedding="electrostatic",
        printlevel=1,
    )

    # Active region: atoms within actradius of the origin atom, else just the QM atoms
    actradius = neb_inputs.get("actradius", qmmm_inputs.get("actradius"))
    if getattr(args, "actatoms", None):
        actatoms = read_intlist_from_file(args.actatoms)
    elif actradius:
        originatom = neb_inputs.get("actregion_origin", qmmm_inputs.get("actregion_origin"))
        if originatom is None:
            if not boundary_excluded:
                raise ValueError("actradius is set but there is no actregion_origin and no "
                                 "non-standard QM residue to centre the active region on")
            originatom = boundary_excluded[0]
        actregiondefine(mmtheory=mm, fragment=react, radius=actradius, originatom=originatom)
        actatoms = read_intlist_from_file("active_atoms")
    else:
        actatoms = qm_atoms
    print(f"[INFO] NEB active region: {len(actatoms)} atoms")

    #Run NEB to find saddlepoint. Returns an ASH Results object
    NEB_result = NEB(reactant=react, product=prod, theory=qmmm_theory, images=numimages,
                     runmode=runmode, numcores=neb_workers,
                     ActiveRegion=True, actatoms=actatoms, charge=charge, mult=mult)
    print(NEB_result)

    saddlepoint = getattr(NEB_result, "saddlepoint_fragment", None)
    if saddlepoint is None:
        print("[WARNING] NEB did not return a saddlepoint (not converged or CI disabled)")
        save_result(result, filename=args.config, section="neb")
        return result

    saddlepoint.write_xyzfile(xyzfilename="saddlepoint.xyz")
    result["saddlepoint_xyz"] = "saddlepoint.xyz"
    result["saddlepoint_energy"] = float(NEB_result.energy)

    #Optional NumFreq job on saddlepoint (QM atoms only) to confirm that a saddlepoint was found.
    if neb_inputs.get("numfreq", True):
        freq_result = NumFreq(theory=qmmm_theory, fragment=saddlepoint, hessatoms=qm_atoms,
                              charge=charge, mult=mult, numcores=neb_workers, runmode=runmode)
        result["frequencies"] = [float(f) for f in freq_result.frequencies]

    save_result(result, filename=args.config, section="neb")
    return result
