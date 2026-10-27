#!/usr/bin/env python3
"""
QM/MM single point (PySCF) on a finished geometry, with density, MEP and
HOMO/LUMO cube files.

A separate stage from `qmmm` because the two want opposite machines. A DFT
single point is one big SCF: it scales across cores and wants a lot of memory
(measured on Bunya: 96% CPU on 16 cores, 56 GB for ~120 QM atoms). The xTB
QM/MM MD and optimisation in front of it saturate around four threads and use
under a gigabyte. Running both in one job means either paying for 64 idle cores
through a 25 ps MD or starving the SCF, so they are submitted separately and
sized separately.

Nothing about the result changes: this is the same `QMMMTheory` single point
that used to run at the end of `qmmm`, on the same geometry.

    mcube sp config.yaml --section sp_product   # sp_product.inputs.pdb_file

Inputs are read from `qmmm.inputs`, then `sp.inputs`, then the `--section`
inputs, so the force fields, charge, multiplicity and QM region do not have to
be repeated and a config written before this stage existed still works. The
level is set by the `sp_*` keys documented in the README:

    sp:
      inputs:
        sp_functional: wb97x-v
        sp_basis: def2-svp
        sp_dispersion: d3bj        # post-SCF D3 on the QM region
        sp_cores: 64
        sp_memory: 192000
        cubes: true
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from .preprocess import load_section, resolve_cores, save_result
from .qmmm import (
    apply_stage_metal_parameters,
    d3_energy,
    make_dft,
    make_mm,
    mm_settings,
    pyscf_memory,
    resolve_qm_region,
    resolve_stage_atoms,
    resolve_start_structure,
)

from ash import *


def sp_inputs(config: str, section: Optional[str] = None) -> Dict[str, Any]:
    """`qmmm.inputs`, overridden by `sp.inputs`, overridden by `section`'s inputs.

    Layered this way so a single point inherits the model it is scoring -- same
    force fields, same charge, same QM region -- and only the level and the
    resources have to be stated.
    """
    inputs: Dict[str, Any] = dict(load_section(config, "qmmm").get("inputs") or {})
    for name in ("sp", section):
        if name and name != "qmmm":
            inputs.update(load_section(config, name).get("inputs") or {})
    return inputs


def write_cubes(pyscf, mult: Optional[int]) -> list:
    """density, MEP and HOMO/LUMO cubes for the QM region as PySCF saw it."""
    mol, mf = pyscf.mol, pyscf.mf
    dm = mf.make_rdm1()
    if mult != 1:
        dm = dm[0] + dm[1]  # total density from alpha + beta
    pyscf.cubegen_density(mol, "density.cube", dm)
    pyscf.cubegen_mep(mol, "mep.cube", dm)

    # Alpha orbitals for UKS
    mo_coeff = mf.mo_coeff if mult == 1 else mf.mo_coeff[0]
    homo = mol.nelec[0] - 1
    pyscf.cubegen_orbital(mol, "HOMO.cube", mo_coeff[:, homo])
    pyscf.cubegen_orbital(mol, "LUMO.cube", mo_coeff[:, homo + 1])
    return [Path(f).resolve().as_posix()
            for f in ("density.cube", "mep.cube", "HOMO.cube", "LUMO.cube")]


def run(args) -> Dict[str, Any]:
    preprocess = load_section(args.config, "preprocess")
    md = load_section(args.config, "md")

    section = getattr(args, "section", None) or "sp"
    inputs = sp_inputs(args.config, section)
    if not inputs:
        raise KeyError(f"No qmmm.inputs or sp.inputs section in {args.config}")

    # A single point has to be told which geometry to score: there is no sensible
    # default, and silently taking the last MD frame would score the wrong thing.
    # So no fallback here, unlike every other stage.
    pdb_file, pdb_source = resolve_start_structure(
        args.config, section, inputs, stage="sp", what="geometry to score")
    print(f"[INFO] Geometry: {pdb_file} (from {pdb_source})")

    forcefields, implicit = mm_settings(inputs, md)
    mult = inputs.get("mult")

    print("[STEP] Building MM system...")
    # The MM side only supplies point charges here, so it needs no threads of its own.
    mm = make_mm(forcefields, pdb_file, implicit, 1)
    inputs = resolve_stage_atoms(inputs, mm)
    fragment = Fragment(pdbfile=pdb_file)

    print("[STEP] Building QM region...")
    qm_atoms, boundary_excluded = resolve_qm_region(inputs, preprocess, pdb_file)
    # The same metal parameters as md, before the QMMMTheory is built
    apply_stage_metal_parameters(inputs, preprocess, md, mm, forcefields, qm_atoms)
    result: Dict[str, Any] = {
        "pdb_file": Path(pdb_file).resolve().as_posix(),
        "qm_atoms": qm_atoms,
        "boundary_excluded": boundary_excluded,
    }
    save_result(result, filename=args.config, section=section)

    # No continuum solvation: ASH's PySCF interface ignores it when point charges
    # are present, so asking for it would quietly do nothing.
    sp_cores = resolve_cores(inputs.get("sp_cores"), "sp_cores")
    print(f"[INFO] PySCF single point: {sp_cores} threads, "
          f"max_memory {pyscf_memory(inputs, 'sp')} MB")
    pyscf = make_dft(inputs, sp_cores, mult, "sp")
    theory = QMMMTheory(
        qm_theory=pyscf,
        mm_theory=mm,
        fragment=fragment,
        qmatoms=qm_atoms,
        excludeboundaryatomlist=boundary_excluded,
        embedding="electrostatic",
        printlevel=1,
    )

    print("[STEP] Running PySCF QM/MM single point...")
    sp = Singlepoint(theory=theory, fragment=fragment, charge=inputs.get("charge"), mult=mult)
    result["sp_energy"] = float(sp.energy)
    result["sp_level"] = f"{inputs.get('sp_functional', 'wb97x-v')}/{inputs.get('sp_basis', 'def2-svp')}"

    if inputs.get("sp_dispersion"):
        # Post-SCF D3 on the QM region (link atoms included) exactly as PySCF saw it
        disp = d3_energy([pyscf.mol.atom_symbol(i) for i in range(pyscf.mol.natm)],
                         pyscf.mol.atom_coords(unit="Bohr"),
                         inputs.get("sp_dispersion_method", inputs.get("sp_functional", "wb97x-v")),
                         inputs["sp_dispersion"])
        result["sp_dispersion"] = inputs["sp_dispersion"]
        result["sp_disp_energy"] = disp
        result["sp_energy_disp"] = result["sp_energy"] + disp
        print(f"[INFO] {inputs['sp_dispersion']} dispersion: {disp:.8f} Eh; "
              f"E(SCF + disp) = {result['sp_energy_disp']:.8f} Eh")
    save_result(result, filename=args.config, section=section)

    if inputs.get("cubes", True):
        result["cubes"] = write_cubes(pyscf, mult)
        save_result(result, filename=args.config, section=section)

    return result
