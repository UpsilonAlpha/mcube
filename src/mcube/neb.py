"""
Climbing-image NEB between two QM/MM minima (xTB/MM), via ASH's KNARR interface.

Settings come from neb.inputs, overridden by the inputs of the section named by
--section (results are written to that section, default `neb`), so several
elementary steps can share one config. Any NEB keyword not given keeps ASH's
default. `method: nebts` runs ASH's NEBTS instead: a loose NEB-CI followed by a
geomeTRIC TS optimisation of the climbing image with a partial Hessian.
"""

from __future__ import annotations

from typing import Dict, Any

from .preprocess import (
    atom_names,
    save_result,
    load_section,
)

from .qmmm import (
    apply_stage_metal_parameters,
    resolve_qm_region,
    stage_inputs,
    mm_settings,
    make_mm,
    make_qmmm,
    resolve_actatoms,
    write_pdb,
    resolve_stage_atoms,
    qm_level_inputs,
    qm_level_label,
)
from ash import *

# neb.inputs keys passed straight to ASH's NEB / NEBTS (see ash/interfaces/interface_knarr.py)
NEB_KEYS = (
    "CI", "free_end", "maxiter", "conv_type", "tol_scale", "tol_max_fci", "tol_rms_fci",
    "cichange_freqcheck", "tol_turn_on_ci", "interpolation", "idpp_maxiter", "idpp_springconst",
    "restart_file", "TS_guess_file",
)
NEB_ONLY_KEYS = ("energy_weighted", "springconst", "springconst2", "zoom")
NEBTS_ONLY_KEYS = (
    "OptTS_maxiter", "OptTS_convergence_setting", "OptTS_conv_criteria", "OptTS_coordsystem",
    "hessian_for_TS", "modelhessian", "tsmode_tangent_threshold", "partial_hessian_atoms",
    "partial_npoint_choice",
)


def check_same_system(reactant_pdb: str, product_pdb: str, react, prod) -> None:
    """Refuse endpoints that are not the same atoms in the same order.

    NEB interpolates atom i of the reactant toward atom i of the product, so a
    product whose atom order differs gives a path between atoms that are not
    counterparts - a nonsense reaction coordinate that still converges. Equal
    atom counts do not establish the correspondence: two waters exchanged, or a
    product rebuilt through its own `preprocess` run, leaves the count intact
    and the pairing broken.

    A product built with `scan` passes by construction, because scan writes
    every point through the reactant's own MM topology. One solvated or
    parameterised separately will not, however carefully it was prepared.
    """
    if react.numatoms != prod.numatoms:
        raise ValueError(f"Reactant ({react.numatoms} atoms) and product "
                         f"({prod.numatoms} atoms) differ: NEB needs one system in "
                         "two geometries. Build the product with `scan` from the "
                         "reactant rather than preparing it separately.")
    left, right = atom_names(reactant_pdb), atom_names(product_pdb)
    bad = [i for i, (a, b) in enumerate(zip(left, right)) if a != b]
    if bad:
        detail = "".join("    atom %d: reactant %s, product %s%s" % (i, left[i], right[i], chr(10))
                         for i in bad[:5])
        raise ValueError(
            f"Reactant and product are the same {react.numatoms} atoms in a "
            f"different order ({len(bad)} positions differ):" + chr(10) + detail +
            chr(10) + "NEB would interpolate between atoms that are not "
            "counterparts. Build the product with `scan` from this reactant so "
            "the atom order is identical by construction.")
    print(f"[INFO] Endpoints are the same {react.numatoms} atoms in the same order")


def check_path_shape(energies) -> None:
    """Warn when an interior image lies below both endpoints.

    Either reading is worth knowing. It can be a genuine intermediate, in which
    case the step is not elementary and belongs as two NEBs with that minimum
    optimised as the endpoint of each - a single NEB over two barriers returns
    one saddle and hides the other.

    Or it is the water treatment. `make_mm` builds every QM/MM system with
    `rigidwater=False`, and `qmmm`, `scan` and `ts` restore rigidity for
    active-region MM waters as *optimizer* constraints. ASH's NEB takes no
    constraints, so those waters are flexible along the path and were rigid in
    the endpoint optimisations that produced it. Interior images can then shed
    energy by relaxing water geometry that has nothing to do with the reaction,
    and the dip is the visible symptom of endpoints and path sitting on
    different surfaces.
    """
    if len(energies) < 3:
        return
    interior, ends = energies[1:-1], (energies[0], energies[-1])
    low = min(interior)
    if low < min(ends):
        drop = (min(ends) - low) * 627.509474
        print(f"[WARNING] an interior NEB image is {drop:.1f} kcal/mol below both "
              "endpoints. Either the step is not elementary (split it at that "
              "minimum), or the endpoints were optimised with active-region MM "
              "waters held rigid while the path leaves them free - see "
              "check_path_shape.")


def run(args) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    section = getattr(args, "section", None) or "neb"

    qmmm_cfg = load_section(args.config, "qmmm")
    md = load_section(args.config, "md")
    preprocess = load_section(args.config, "preprocess")
    qmmm_inputs = qmmm_cfg.get("inputs") or {}
    neb_inputs = stage_inputs(args.config, "neb", section)
    qmmm_inputs = qm_level_inputs(qmmm_inputs, neb_inputs)  # qm_engine/dft_* may be set per section

    charge = qmmm_inputs["charge"]
    mult = qmmm_inputs["mult"]
    numimages = neb_inputs.get("images", 12)
    method = str(neb_inputs.get("method", "neb")).lower()
    if method not in ("neb", "nebts"):
        raise ValueError(f"neb method must be 'neb' or 'nebts', not {method!r}")

    # ASH runs up to `numcores` images in parallel, each with the theory's own numcores,
    # so split the cores instead of giving both the full count
    neb_workers = max(1, min(numimages, args.cores))
    theory_cores = max(1, args.cores // neb_workers)
    runmode = 'parallel' if neb_workers > 1 else 'serial'
    if neb_workers > 1:
        # Parallel images each start their own xTB, so this stage must run with
        # no thread pinning: every image would otherwise bind to the same cores.
        print(f"[INFO] {neb_workers} NEB images run in parallel: run this stage "
              "with OMP_PROC_BIND/OMP_PLACES unset (the xtb-parallel recipe)")

    forcefields, implicit = mm_settings(qmmm_inputs, md)

    print("[STEP] Building MM system...")
    mm = make_mm(forcefields, args.reactant, implicit, theory_cores)  # OpenMM threads per image
    neb_inputs = resolve_stage_atoms(neb_inputs, mm)
    qmmm_inputs = resolve_stage_atoms(qmmm_inputs, mm)

    print("[STEP] Building QM region...")
    qm_atoms, boundary_excluded = resolve_qm_region({**qmmm_inputs, **{k: neb_inputs[k] for k in ('qm_residues', 'extra_qm_residues') if k in neb_inputs}}, preprocess, args.reactant)
    # The same metal parameters as md, before the QMMMTheory is built
    apply_stage_metal_parameters(qmmm_inputs, preprocess, md, mm, forcefields, qm_atoms)
    result["qm_level"] = qm_level_label(qmmm_inputs)
    result["qm_atoms"] = qm_atoms
    result["boundary_excluded"] = boundary_excluded

    ################################################
    # Defining reactant and product ASH fragments
    #################################################
    react=Fragment(pdbfile=args.reactant, charge=charge, mult=mult)
    prod=Fragment(pdbfile=args.product, charge=charge, mult=mult)

    check_same_system(args.reactant, args.product, react, prod)

    # Not the xTB library: ASH's library path cannot return a QM/MM point-charge
    # gradient and calls ashexit() when asked for one. Parallel images each get their
    # own working directory from ASH, so inputfile mode does not collide.
    qmmm_theory = make_qmmm(qmmm_inputs, mm, react, qm_atoms, boundary_excluded, theory_cores)

    # Active region: --actatoms file, else actradius around actregion_origin, else the QM atoms
    act_inputs = {**qmmm_inputs, **{k: neb_inputs[k] for k in ("actradius", "actregion_origin") if k in neb_inputs}}
    actatoms = resolve_actatoms(getattr(args, "actatoms", None), act_inputs, mm, react,
                                qm_atoms, boundary_excluded)
    print(f"[INFO] NEB active region: {len(actatoms)} atoms")

    keys = NEB_KEYS + (NEBTS_ONLY_KEYS if method == "nebts" else NEB_ONLY_KEYS)
    options = {k: neb_inputs[k] for k in keys if k in neb_inputs}
    if method == "nebts":
        options.setdefault("hessian_for_TS", "partial")
    result["method"] = method
    result["images"] = numimages
    result["options"] = options
    print(f"[INFO] {method.upper()} with {numimages} images; options: {options or 'ASH defaults'}")

    #Run NEB to find saddlepoint. Returns an ASH Results object
    neb_function = NEBTS if method == "nebts" else NEB
    NEB_result = neb_function(reactant=react, product=prod, theory=qmmm_theory, images=numimages,
                              runmode=runmode, numcores=neb_workers,
                              ActiveRegion=True, actatoms=actatoms, charge=charge, mult=mult, **options)
    print(NEB_result)

    energies = getattr(NEB_result, "MEP_energies_dict", None)
    if energies:
        result["mep_energies"] = {str(k): float(v) for k, v in energies.items()}
        check_path_shape(list(energies.values()))

    saddlepoint = getattr(NEB_result, "saddlepoint_fragment", None) if NEB_result is not None else None
    if saddlepoint is None:
        print("[WARNING] NEB did not return a saddlepoint (not converged or CI disabled)")
        save_result(result, filename=args.config, section=section)
        return result

    saddlepoint.write_xyzfile(xyzfilename="saddlepoint.xyz")
    result["saddlepoint_xyz"] = "saddlepoint.xyz"
    # Full system with the MM topology: the starting point for `mcube ts`
    result["saddlepoint_pdb"] = write_pdb(mm, saddlepoint.coords, "saddlepoint.pdb")
    result["saddlepoint_energy"] = float(NEB_result.energy)

    #Optional NumFreq job on saddlepoint to confirm that a saddlepoint was found.
    #Cost is ~6 gradients per Hessian atom, so `hessian_atoms` can be set to the reacting atoms.
    if neb_inputs.get("numfreq", True):
        hessatoms = neb_inputs.get("hessian_atoms") or qm_atoms
        freq_result = NumFreq(theory=qmmm_theory, fragment=saddlepoint, hessatoms=hessatoms,
                              charge=charge, mult=mult, numcores=neb_workers, runmode=runmode)
        result["frequencies"] = [float(f) for f in freq_result.frequencies]
        result["imaginary_frequencies"] = [f for f in result["frequencies"] if f < 0]

    save_result(result, filename=args.config, section=section)
    return result
