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

The DFT single point that used to end this stage is now `mcube sp` (sp.py): it
scales across cores and wants tens of gigabytes, while the xTB work here
saturates about four threads, so the two are sized and submitted separately.
"""

from __future__ import annotations

import copy
import os
from typing import Dict, Any, List, Optional
from pathlib import Path

import numpy as np

from openmm import unit
from openmm.app import PDBFile

from . import ions
from .preprocess import (
    resolve_cores,
    save_result,
    save_input,
    load_section,
    build_qm_region,
    qm_region_settings,
    resolve_atom_inputs,
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

    # `extra_qm_residues: {chain: [resid, ...]}` adds residues (e.g. second-shell waters)
    # to the preprocess QM region for this stage only
    extra = inputs.get("extra_qm_residues")
    if extra:
        base = {str(ch): {str(k): v for k, v in res.items()} for ch, res in (results.get("qm_residues") or {}).items()}
        topology = PDBFile(pdb_file).topology
        for ch, resids in extra.items():
            wanted = {str(r) for r in resids}
            found = {r.id.strip(): r.name for r in topology.residues()
                     if r.chain.id == str(ch) and r.id.strip() in wanted}
            missing = wanted - set(found)
            if missing:
                raise ValueError(f"extra_qm_residues: chain {ch} residues {sorted(missing)} not in {pdb_file}")
            base.setdefault(str(ch), {}).update(found)
            print(f"[INFO] Extra QM residues for this stage: chain {ch} {sorted(found.items())}")
        return build_qm_region(pdb_file, base, **settings)

    if results.get("qm_atoms") is not None:
        print(f"[INFO] Using {len(results['qm_atoms'])} QM atoms saved by preprocess")
        return results["qm_atoms"], results.get("boundary_excluded", [])

    return build_qm_region(pdb_file, results["qm_residues"], **settings)


# ==============================
# Helpers shared by the QM/MM stages (qmmm, scan, neb, ts)
# ==============================

# Input keys that hold atoms; each atom may be an index or a "chain:resid:atom" selector
ATOM_KEYS = ("actregion_origin", "opt_constraints", "atoms", "atoms2", "constraints", "hessian_atoms",
             "partial_hessian_atoms", "monitor")


def resolve_stage_atoms(inputs: Dict[str, Any], mm) -> Dict[str, Any]:
    """Atom selectors in a stage's inputs -> indices, using the MM topology."""
    return resolve_atom_inputs(inputs, mm.topology, ATOM_KEYS)


def resolve_start_structure(config: str, section: str, inputs: Dict[str, Any],
                            md_results: Dict[str, Any], pdb_arg: Optional[str]):
    """The structure this stage runs on, and where it came from.

    Precedence: `--pdb`, then this section's `pdb_file`, then the md stage's
    `lastframe`. That last one is a convenience for the usual
    preprocess -> md -> qmmm chain, and it is refused when it reaches outside
    the config's own directory, because it is then almost certainly a frame
    from a different run - and QM/MM started on the wrong structure converges
    and looks entirely normal, so nothing downstream would catch it.
    """
    for value, source in ((pdb_arg, "--pdb"),
                          (inputs.get("pdb_file"), f"{section}.inputs.pdb_file"),
                          (md_results.get("lastframe"), "md.results.lastframe")):
        if value:
            break
    else:
        raise KeyError(
            "No starting structure for the qmmm stage. Set `pdb_file` in "
            f"`{section}.inputs`, pass `--pdb`, or run the `md` stage first so that "
            "`md.results.lastframe` exists.")
    path = Path(value)
    if not path.is_file():
        raise FileNotFoundError(
            f"Starting structure named by {source} does not exist: {value}")
    if source == "md.results.lastframe":
        home = Path(config).resolve().parent
        if home not in path.resolve().parents:
            raise ValueError(
                f"The qmmm stage would fall back to {value} (md.results.lastframe), "
                f"which is outside this config's directory {home}. That is usually a "
                "frame from another run, and QM/MM started on the wrong structure "
                "converges and looks fine, so set `pdb_file` in "
                f"`{section}.inputs` to the frame you mean, or pass --pdb.")
    return path.resolve().as_posix(), source


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


def make_mm(forcefields, pdb_file: str, implicit: bool, cores: int, platform: str = "CPU",
            properties=None):
    # `platform` is honoured only here: ASH ignores OpenMM_MD(platform=) for QM/MM
    return OpenMMTheory(
        xmlfiles=forcefields,
        pdbfile=pdb_file,
        autoconstraints=None,
        periodic=not implicit,
        rigidwater=False,
        numcores=cores,  # OpenMM CPU threads; ASH defaults to 1
        platform=platform or "CPU",
        properties=properties,
    )


def apply_stage_metal_parameters(inputs: Dict[str, Any], preprocess: Dict[str, Any],
                                 md: Dict[str, Any], mm, forcefields, qm_atoms) -> str:
    """
    Give the metals the same nonbonded parameters in QM/MM as in md.

    Before this, md swapped in IOD (or 12-6-4) radii but every QM/MM stage built
    its MM system from the force-field XMLs alone, so a QM metal met the MM
    atoms through Amber's stock Lennard-Jones terms. Called after the QM region
    is known and before the QMMMTheory is built: 12-6-4's C4 term is switched
    off for QM atoms, whose polarisation the embedding already supplies.
    """
    which = ions.metal_parameter_set(inputs, md)
    metals = (preprocess.get("inputs") or {}).get("metals")
    if not metals or which == "none":
        return which
    polarizabilities = inputs.get("polarizabilities", (md.get("inputs") or {}).get("polarizabilities"))
    ions.apply_metal_parameters(mm.system, mm.topology, metals, which, forcefields=forcefields,
                                qm_atoms=qm_atoms, polarizabilities=polarizabilities)
    return which


def level_option(inputs: Dict[str, Any], prefix: str, name: str, default=None):
    """
    A QM level option, read as `<prefix>_<name>` and falling back to `sp_<name>`.

    The fallback means a config that already defines the DFT level for single points
    does not have to repeat it to run the gradients at the same level.
    """
    for key in (f"{prefix}_{name}", f"sp_{name}"):
        if key in inputs:
            return inputs[key]
    return default


def pyscf_memory(inputs: Dict[str, Any], prefix: str = "sp") -> int:
    """
    PySCF max_memory in MB. ASH's default (3100) sends density-fitting integrals to
    disk; by default use 80% of the SLURM job's memory so they stay in RAM.
    """
    memory = level_option(inputs, prefix, "memory")
    if memory is None:
        slurm_mb = os.environ.get("SLURM_MEM_PER_NODE")
        memory = int(0.8 * int(slurm_mb)) if slurm_mb and slurm_mb.isdigit() else 3100
    return int(memory)


def ecp_from_file(path: str) -> Dict[str, Any]:
    """{element: ECP} for every element with an ECP block in an NWChem-format file
    (e.g. written by basis_set_exchange), parsed for PySCF."""
    import re
    from pyscf import gto
    symbols = sorted(set(re.findall(r"^\s*([A-Z][a-z]?)\s+nelec\s+\d+", Path(path).read_text(), re.M)))
    return {s: gto.basis.load_ecp(path, s) for s in symbols}


def dft_ecp(inputs: Dict[str, Any], prefix: str, basis):
    """
    The effective core potentials for a DFT theory.

    `<prefix>_ecp_file` (NWChem format, from basis_set_exchange) wins, then `<prefix>_ecp`
    (a name PySCF knows). Otherwise a def2 basis gets the def2 ECPs of the same name: def2
    basis sets for Rb onwards are built for those ECPs, and PySCF applies them only to the
    elements that have one, so lighter elements are unaffected. Without this, a La or Lu
    QM region would be computed all-electron in a basis meant for a 46- or 28-electron core.
    """
    ecp_file = level_option(inputs, prefix, "ecp_file")
    if ecp_file:
        return ecp_from_file(ecp_file)
    ecp = level_option(inputs, prefix, "ecp")
    if ecp is None and isinstance(basis, str) and basis.lower().startswith("def2"):
        ecp = basis
    return ecp


def make_dft(inputs: Dict[str, Any], cores: int, mult: Optional[int], prefix: str = "sp"):
    """
    PySCF theory for QM/MM energies or gradients.

    `dispersion` is ASH's own keyword, so it enters the gradient and a DFT/MM
    optimisation is dispersion-corrected. The explicit `sp_dispersion` correction
    applied after a single point (`d3_energy`) is separate and stays post-SCF.
    """
    functional = level_option(inputs, prefix, "functional", "wb97x-v")
    basis_file = level_option(inputs, prefix, "basis_file")
    basis = None if basis_file else level_option(inputs, prefix, "basis", "def2-svp")
    return PySCFTheory(scf_type="RKS" if mult == 1 else "UKS",
                       functional=functional,
                       basis=basis,
                       basis_file=basis_file,
                       ecp=dft_ecp(inputs, prefix, basis),
                       gridlevel=level_option(inputs, prefix, "gridlevel", 3),
                       densityfit=level_option(inputs, prefix, "densityfit", True),
                       auxbasis=level_option(inputs, prefix, "auxbasis"),
                       dispersion=level_option(inputs, prefix, "scf_dispersion"),
                       platform=level_option(inputs, prefix, "platform", "CPU"),
                       # On GPU the MM point-charge gradient can also run on the device,
                       # which is most of the QM/MM gradient for a large MM environment
                       GPU_pcgrad=bool(level_option(inputs, prefix, "gpu_pcgrad", False)),
                       memory=pyscf_memory(inputs, prefix),
                       numcores=cores)


def qm_level_label(inputs: Dict[str, Any]) -> str:
    """How the QM energies in a stage's results were produced.

    Recorded by every stage that reports an energy, because a barrier assembled
    from `opt_energy`, a scan's `energy`, `mep_energies` or `saddlepoint_energy`
    is otherwise indistinguishable from one computed at another level. That
    matters here specifically: GFN2-xTB is a path-finding method and its
    energies are not barriers. Published practice is to find the path
    semiempirically and refine the stationary points at DFT - Gervasoni et al.,
    Proteins 2022, 90, 372 refine with B3LYP/MM, and Guven et al., J. Phys.
    Chem. B 2024, 128, 12416 sample with DFTB3 and refine with DFT.
    """
    engine = str(inputs.get("qm_engine", "xtb")).lower()
    if engine in ("dft", "pyscf"):
        return "%s/%s" % (level_option(inputs, "dft", "functional", "wb97x-v"),
                          level_option(inputs, "dft", "basis", "def2-svp"))
    return "%s-xTB" % inputs.get("xtb_method", "GFN2")


def make_qm(inputs: Dict[str, Any], cores: int, runmode: Optional[str] = None):
    """
    The QM theory for the gradient-driven stages (optimisation, scan, NEB, TS, frequencies).

    `qm_engine: xtb` (the default) is GFN2-xTB; `qm_engine: dft` runs those stages at
    DFT/MM instead, at the `dft_*` level (falling back to the `sp_*` keys). DFT gradients
    cost orders of magnitude more than xTB ones, so a DFT engine is for refining
    stationary points that an xTB path has already located, not for exploring.
    """
    engine = str(inputs.get("qm_engine", "xtb")).lower()
    if engine in ("dft", "pyscf"):
        dft_cores = resolve_cores(inputs.get("dft_cores", cores), "dft_cores")
        theory = make_dft(inputs, dft_cores, inputs.get("mult"), "dft")
        print(f"[INFO] QM engine: DFT ({theory.functional}/"
              f"{level_option(inputs, 'dft', 'basis', 'def2-svp')}), {dft_cores} threads, "
              f"platform {level_option(inputs, 'dft', 'platform', 'CPU')}")
        return theory
    if engine != "xtb":
        raise ValueError(f"qm_engine must be 'xtb' or 'dft', not {engine!r}")
    runmode = runmode or inputs.get("xtb_runmode", "inputfile")
    if str(runmode).lower() == "library":
        raise ValueError(
            "xtb_runmode 'library' cannot be used here: mcube only ever runs xTB as the "
            "QM region of a QMMMTheory, and every stage that does so needs gradients. "
            "ASH's library path cannot return a point-charge gradient, because xtb-python "
            "does not wrap xtb_getPCGradient; it calls ashexit() there, which terminates "
            "the process rather than raising. Refused at config time instead, so the key "
            "to change is named. Use 'inputfile'.")
    # xTB's thread placement comes from the job script's environment, which it
    # inherits as a separate process.
    # SCC controls, ASH's defaults unless set. A highly charged metal next to bare MM point
    # charges can close the HOMO-LUMO gap and stop the SCC converging mid-MD (an Ln3+ aqua
    # ion did after ~6 ps); a higher electronic temperature or softer point charges is the
    # usual remedy.
    scc = {key: inputs[name] for key, name in (("electronic_temp", "xtb_etemp"),
                                               ("maxiter", "xtb_maxiter"),
                                               ("hardness_PC", "xtb_pc_hardness"))
           if inputs.get(name) is not None}
    return xTBTheory(xtbmethod=inputs.get("xtb_method", "GFN2"), numcores=cores, runmode=runmode,
                     **scc)


QM_LEVEL_KEYS = ("qm_engine", "xtb_method", "xtb_runmode", "xtb_etemp", "xtb_maxiter",
                 "xtb_pc_hardness")


def qm_level_inputs(qmmm_inputs: Dict[str, Any], stage_inputs_: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """
    `qmmm.inputs` with the QM-level keys overridden by a stage's own section.

    The scan, neb and ts stages take the QM theory from `qmmm.inputs`, so without
    this a `qm_engine: dft` set on one of those sections would be silently ignored
    and the stage would quietly run at xTB. Everything else about the QM/MM setup
    (charge, mult, force fields) still comes from `qmmm.inputs`.
    """
    merged = dict(qmmm_inputs)
    for key, value in (stage_inputs_ or {}).items():
        if key in QM_LEVEL_KEYS or key.startswith(("dft_", "sp_")):
            merged[key] = value
    return merged


def make_qmmm(qmmm_inputs: Dict[str, Any], mm, fragment, qm_atoms, boundary_excluded,
              cores: int, runmode: Optional[str] = None):
    return QMMMTheory(
        qm_theory=make_qm(qmmm_inputs, cores, runmode),
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


def water_constraints(mm, actatoms, qm_atoms=()) -> List[list]:
    """
    O-H and H-H pairs of active-region *MM* TIP3P waters, to keep them rigid while optimising.

    Waters in the QM region are left free: freezing their O-H bonds would stop them
    giving up or taking a proton (nucleophilic water, proton relays).
    """
    qm = set(int(a) for a in qm_atoms)
    return [[int(p[0]), int(p[1])] for p in
            getwaterconstraintslist(openmmtheoryobject=mm, atomlist=actatoms, watermodel='tip3p')
            if int(p[0]) not in qm and int(p[1]) not in qm]


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


def last_traj_frame(xyzfile: str, natoms: int):
    """Coordinates (Å) and energy of the last frame of ASH's geometric_OPTtraj_Full.xyz,
    read from the end of the file (it holds the full system at every step)."""
    path = Path(xyzfile)
    with open(path, "rb") as f:
        f.seek(0, 2)
        size = f.tell()
        f.seek(max(0, size - (natoms + 2) * 120))
        lines = f.read().decode().splitlines()
    block = lines[-(natoms + 2):]
    if len(block) != natoms + 2 or block[0].strip() != str(natoms):
        raise ValueError(f"Could not read the last {natoms}-atom frame of {xyzfile}")
    coords = np.array([[float(x) for x in line.split()[1:4]] for line in block[2:]])
    words = block[1].split()
    energy = float(words[words.index("Energy") + 1]) if "Energy" in words else None
    return coords, energy


def optimize(fragment, theory, actatoms, inputs, constraints, constrainvalue, maxiter,
             coordsystem=None):
    """
    ASH/geomeTRIC optimisation that survives running out of steps: geomeTRIC
    raises at maxiter and ASH then leaves the fragment at its start, so on
    non-convergence the fragment is set to the last step of the trajectory.
    Returns (energy, converged).

    `coordsystem` is geomeTRIC's internal-coordinate system, from
    `opt_coordsystem` (`tric`, `hdlc`, `dlc`, `prim` or `cart`). Unset leaves
    ASH's own default, so existing runs are unaffected. It is worth reaching for
    when an optimisation reports `SVD Error ... ( nan -> nan)`: the coordinate
    system has degenerated, and geomeTRIC then retries in a loop, writing a
    ~23 MB dump each time until the job is killed for memory rather than
    failing. `cart` cannot degenerate that way, at the cost of slower
    convergence.
    """
    try:
        extra = {"coordsystem": coordsystem} if coordsystem else {}
        opt = Optimizer(fragment=fragment, theory=theory, ActiveRegion=True, actatoms=actatoms,
                        maxiter=maxiter, constraints=constraints, constrainvalue=constrainvalue,
                        charge=inputs.get("charge"), mult=inputs.get("mult"), **extra)
    except Exception as e:
        if type(e).__name__ != "GeomOptNotConvergedError":
            raise
        print(f"[WARNING] Optimisation not converged in {maxiter} steps; keeping its last geometry")
        coords, energy = last_traj_frame("geometric_OPTtraj_Full.xyz", fragment.numatoms)
        fragment.coords = coords
        return energy, False
    energy = getattr(opt, "energy", None)
    return (float(energy) if energy is not None else None), True


def d3_energy(symbols, coords_bohr, method: str, damping: str = "d3bj") -> float:
    """
    DFT-D3 dispersion energy (Eh) from the simple-dftd3 library (`dftd3-python`).
    damping: d3bj (rational) or d3zero. `method` is the functional's name in
    simple-dftd3's parameter table (e.g. wb97x, b3lyp).
    """
    from dftd3.interface import DispersionModel, RationalDampingParam, ZeroDampingParam
    from pyscf.data.elements import charge
    damping = damping.lower()
    if damping not in ("d3bj", "d3zero"):
        raise ValueError(f"sp_dispersion must be d3bj or d3zero, not {damping!r}")
    param = (RationalDampingParam if damping == "d3bj" else ZeroDampingParam)(method=method.replace("_", "-"))
    model = DispersionModel(np.array([charge(s) for s in symbols]), np.asarray(coords_bohr, dtype=float))
    return float(model.get_dispersion(param, grad=False)["energy"])


def write_pdb(mm, coords, filename: str) -> str:
    """Full system with the MM topology (keeps residue names/ids, e.g. for NEB endpoints)."""
    with open(filename, "w") as f:
        PDBFile.writeFile(mm.topology, np.asarray(coords) * unit.angstrom, f, keepIds=True)
    return Path(filename).resolve().as_posix()



def run(args) -> Dict[str, Any]:

    preprocess = load_section(args.config, "preprocess")
    md = load_section(args.config, "md")
    qmmm = load_section(args.config, "qmmm")
    if qmmm.get("inputs") is None:
        raise KeyError(f"No qmmm.inputs section in {args.config}")
    # Results go to this section, so several structures (e.g. NEB endpoints) can share one config;
    # its own `inputs` (if any) override qmmm.inputs, e.g. skip_equilibrate or opt_constraints per run
    section = getattr(args, "section", None) or "qmmm"
    inputs = stage_inputs(args.config, "qmmm", section)
    md_results = md.get("results") or {}

    pdb_file, pdb_source = resolve_start_structure(
        args.config, section, inputs, md_results, getattr(args, "pdb", None))
    print(f"[INFO] Starting structure: {pdb_file} (from {pdb_source})")
    if not inputs.get("pdb_file"):
        save_input(args.config, section, "pdb_file", pdb_file)

    forcefields, implicit = mm_settings(inputs, md)

    print("[STEP] Building MM system...")
    cores = resolve_cores(inputs.get("cores"))
    platform, properties = inputs.get("platform", "CPU"), inputs.get("properties")
    mm = make_mm(forcefields, pdb_file, implicit, cores, platform, properties)
    inputs = resolve_stage_atoms(inputs, mm)

    fragment = Fragment(pdbfile=pdb_file)

    print("[STEP] Building QM region...")
    qm_atoms, boundary_excluded = resolve_qm_region(inputs, preprocess, pdb_file)
    result: Dict[str, Any] = {
        "qm_atoms": qm_atoms,
        "boundary_excluded": boundary_excluded,
        "pdb_file": pdb_file,
        "pdb_file_source": pdb_source,
    }
    result["qm_level"] = qm_level_label(inputs)
    result["metal_parameters"] = apply_stage_metal_parameters(inputs, preprocess, md, mm,
                                                              forcefields, qm_atoms)
    save_result(result, filename=args.config, section=section)

    qmmm_theory = make_qmmm(inputs, mm, fragment, qm_atoms, boundary_excluded, cores)

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
            # Fresh MM and QM/MM objects: after OpenMM_MD, ASH's QMMMTheory no longer has the
            # attributes the optimiser logs (AttributeError: 'QMMMTheory' object has no attribute
            # 'MMenergy'), which killed job 29937641 at the first optimisation step
            print("[STEP] Rebuilding the QM/MM system from the QM/MM MD end point for the optimisation...")
            mm = make_mm(forcefields, "QM_MM_lastframe.pdb", implicit, cores, platform, properties)
            apply_stage_metal_parameters(inputs, preprocess, md, mm, forcefields, qm_atoms)
            qmmm_theory = make_qmmm(inputs, mm, fragment, qm_atoms, boundary_excluded, cores)

        # An active_atoms file reuses another run's active region (e.g. the reactant's for the product)
        actatoms = resolve_actatoms(getattr(args, "actatoms", None), inputs, mm, fragment,
                                    qm_atoms, boundary_excluded)

        # Optional bond constraints: [[i, j]] holds the current distance, [[i, j, d_Å]] sets it
        user_pairs = inputs.get("opt_constraints") or []
        check_in_active(user_pairs, actatoms, "opt_constraints")
        constraints, constrainvalue = optimizer_constraints(water_constraints(mm, actatoms, qm_atoms), user_pairs,
                                                            fragment.coords)
        if user_pairs:
            result["opt_constraints"] = bond_constraints(user_pairs, fragment.coords)
            print("[INFO] Bond constraints (i, j, Å): " + ", ".join(
                f"{i}-{j} {d:.3f}" for i, j, d in result["opt_constraints"]))

        clean_geometric_files()

        print("[STEP] Running QM/MM Optimisation...")
        energy, converged = optimize(fragment, qmmm_theory, actatoms, inputs, constraints, constrainvalue,
                                     inputs.get("opt_maxiter", 200),
                                     inputs.get("opt_coordsystem"))

        # Full optimised system with the MM topology, e.g. as an NEB endpoint. If the
        # optimisation ran out of steps, this is its last geometry (converged: false).
        result["opttraj"] = Path("geometric_OPTtraj_Full.xyz").resolve().as_posix()
        result["optimized_pdb"] = write_pdb(mm, fragment.coords, "QM_MM_optimized.pdb")
        result["converged"] = converged
        if energy is not None:
            result["opt_energy"] = energy
        save_result(result, filename=args.config, section=section)

    # The DFT single point is its own stage (`mcube sp`): it scales across cores and
    # wants tens of gigabytes, while the xTB MD and optimisation above saturate about
    # four threads, so one job cannot size both. See sp.py.
    if inputs.get("singlepoint") is not None:
        target = result.get("optimized_pdb") or result.get("lastframe") or pdb_file
        print("[INFO] `singlepoint` is no longer part of the qmmm stage. To run it:\n"
              f"         mcube sp {args.config} --pdb {target}\n"
              "       Submit it separately so the SCF can have the cores and memory it "
              "scales with (see sp.py).")

    return result
