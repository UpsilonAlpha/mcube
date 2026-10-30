#!/usr/bin/env python3
"""
Which ion does the constriction site prefer, with the real protein around it? Ion exchange
energies at the Glu ring (docs/channel/question-and-criteria.md) by QM/MM "vertical swaps":
the same MM frame carries each ion in turn, so the environment (protein, POPE, water, Na/Cl,
the salt bridges of Glu300 and Glu220 that the clusters of site_exchange.py could not
represent) enters through electrostatic embedding and is identical for the ions it compares.

For M1 -> M2 in Ca->La, Ca->Lu, La->Lu

    ddE = < E(M2@site) - E(M1@site) >_site frames  -  < E(M2@water) - E(M1@water) >_water frames

Each difference is taken on ONE frame with only the ion's identity (residue, atom name, element,
charge, MM metal parameters) changed. The water term is the same swap on frames of the ion in a
bulk TIP3P box (docs/channel/bulk), with the same protocol, so the long-range 2+ vs 3+ terms and
the Born-type self energy largely cancel. Negative ddE = the site prefers M2. It is an energy,
not a free energy: no thermal or entropic terms, one local minimum per frame, and the average is
over a handful of frames, so the spread across frames is part of the answer.

Two-sided: frames come from La-bound and Lu-bound MM runs (and La and Lu bulk boxes), and every
ion goes into every frame. A ranking that depends on which ion the frames were sampled with shows
up as a difference between the per-source rows.

Per (frame, ion) job, through mcube's own stages and one config.yaml (exchange.slurm):
  preprocess  skip_solvate on the swapped frame (detection report only; no solvation, no
              neutralisation, so the MM environment of a frame is the same for every ion)
  qmmm        NO QM/MM MD (skip_equilibrate): a local QM/MM optimisation from the shared MM
              frame, atoms within --actradius of the ion free, at r2SCAN/def2-SVP (La/Lu
              lcecp-1-SVP) by default. MD per swapped system would let the environments drift
              apart and destroy the pairing the difference relies on. This is a deliberate
              exception to mcube's usual "QM/MM MD before optimisation" order.
  sp          r2SCAN/def2-TZVP (La/Lu lcecp-1-TZVP) QM/MM single point on the optimised geometry

The QM region is fixed once per frame from the SOURCE ion's coordination, by mcube's rule
(preprocess.coordination_shell: the closest donor of each residue within the cutoff, minimum
image) with --qm-cutoff + --qm-margin, and the same `qm_residues` goes into all three ions'
configs. Coordinating amino acids are cut at CA-CB as mcube does; waters come in whole.

  python exchange.py frames --out $S/channel/qmmm \
      --site LA=$S/channel/system/LA/umbrella_res_035 --site LU=$S/channel/system/LU/umbrella_res_chbe_035 \
      --water LA=$S/channel/bulk/LA_1264_frames --water LU=$S/channel/bulk/LU_1264_frames \
      --basis-svp $S/channel/dft/def2-svp_lcecp1.nw --basis-tzvp $S/channel/dft/def2-tzvp_lcecp1.nw
  python exchange.py analyse --out $S/channel/qmmm
  python exchange.py structure --out $S/channel/qmmm
  python exchange.py clean --out $S/channel/qmmm [--delete]

A --site directory holds an umbrella window's start.pdb, window.dcd (and cv.csv for the frame
spacing); a --water directory holds bulk_ion.py's system.pdb and the frames.dcd that
`bulk_ion.py run --dcd-every` writes. `frames` writes <out>/<kind>/<SRC>_f<frame>/<ION>/ with
frame.pdb and config.yaml, the job lists jobs_site.txt / jobs_water.txt, and frames.json.
`analyse` reads each job's config.yaml (sp.results.sp_energy at the TZVP level, and
qmmm.results.opt_energy at the SVP level as a check) and writes <out>/report.json. `structure`
reads each job's QM-region optimisation trajectory (steps, energy drop, QM RMSD, the ion-O shell
before and after) into <out>/structure.json; `clean` then deletes the bulky files of final jobs
that structure.json covers (dry run unless --delete).
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

HARTREE_KCAL = 627.509474
IONS = {"CA": ("Ca", 2), "LA": ("La", 3), "LU": ("Lu", 3)}       # element, charge
EXCHANGES = (("CA", "LA"), ("CA", "LU"), ("LA", "LU"))
# Formal charge of the part of a residue that goes into QM. Amino acids are cut at CA-CB, so
# this is the side chain's charge; the ff14SB backbone (N H CA HA C O) sums to -0.12 for Glu,
# which rounds the side chain to -1 as the formal charge says.
RES_CHARGE = {"GLU": -1, "ASP": -1, "CYM": -1, "LYS": 1, "ARG": 1, "HIP": 1, "NA": 1, "K": 1,
              "CL": -1}
WATER = {"HOH", "WAT", "TIP3"}
LIPIDS = {"POP", "DLP", "DMP", "DOP", "DPP", "PC", "PE", "PA", "OL"}
AMINO_ACIDS = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
               "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HID", "HIE", "HIP", "HSD",
               "HSE", "HSP", "CYX", "CYM", "ASH", "GLH", "LYN"}            # = mcube's AMINO_ACIDS
BACKBONE = {"N", "H", "H1", "H2", "H3", "CA", "HA", "HA2", "HA3", "C", "O", "OXT"}
SITE_FF = ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml", "amber14/lipid17.xml"]
WATER_FF = ["amber14/tip3p.xml"]


# ==============================
# Frames
# ==============================

def find_ion(topology):
    """The one single-atom CA/LA/LU residue of a frame."""
    ions = [a for a in topology.atoms() if a.residue.name in IONS
            and len(list(a.residue.atoms())) == 1 and a.element is not None
            and a.element.symbol == IONS[a.residue.name][0]]
    if len(ions) != 1:
        raise SystemExit(f"expected one Ca/La/Lu ion residue, found {len(ions)}")
    return ions[0]


def load_topology(top_pdb, bonds_from=None):
    """OpenMM topology of a source run. An mcube/ASH-written frame of a membrane system may
    lack CONECT records, so its lipids come back without bonds and cannot be typed (decision
    log F10); `bonds_from` (the system's processed_pdb, same atoms in the same order) supplies
    them, as umbrella.load_structure does."""
    from openmm.app import PDBFile
    pdb = PDBFile(str(top_pdb))
    top = pdb.topology
    if bonds_from:
        ref = PDBFile(str(bonds_from)).topology
        names = lambda t: [(a.name, a.residue.name) for a in t.atoms()]
        if names(ref) != names(top):
            raise SystemExit(f"{bonds_from} does not match {top_pdb} atom for atom")
        ref.setPeriodicBoxVectors(top.getPeriodicBoxVectors())
        top = ref
    unbonded = [r for r in top.residues() if r.name in LIPIDS and not list(r.bonds())]
    if unbonded:
        raise SystemExit(f"{top_pdb}: {len(unbonded)} lipid residues have no bonds (no CONECT "
                         "records); pass --bonds-from SRC=<the system's processed_pdb>")
    return top


def frame_spacing_ps(src: Path, kind: str, n_frames: int, frame_ps):
    """ps between saved frames: the option, else (site) cv.csv's last time / frames as
    residence.py does, else (water) the `dcd_every_ps` bulk_ion.py records in report.json."""
    if frame_ps:
        return float(frame_ps)
    cv = src / "cv.csv"
    if kind == "site" and cv.exists() and cv.stat().st_size > 100:
        last = cv.read_text().strip().splitlines()[-1].split(",")[0]
        return float(last) / n_frames
    rep = src / "report.json"
    if kind == "water" and rep.exists():
        every = json.loads(rep.read_text()).get("dcd_every_ps")
        if every:
            return float(every)
    raise SystemExit(f"{src}: cannot tell the time between frames; pass "
                     f"--{'frame-ps' if kind == 'site' else 'water-frame-ps'}")


def pick_frames(n_total: int, dt_ps: float, n: int, spacing_ns: float, start: float):
    """`n` frame indices, `spacing_ns` apart, counted back from the last frame, all at or after
    the fraction `start` of the run (frame i is at (i + 1) * dt)."""
    step = max(1, int(round(spacing_ns * 1000.0 / dt_ps)))
    idx = [n_total - 1 - k * step for k in range(n)]
    first_ok = start * n_total - 1 - 1e-6
    if idx[-1] < first_ok:
        fit = int((n_total - 1 - math.ceil(first_ok)) // step) + 1
        raise SystemExit(f"{n} frames {spacing_ns} ns apart do not fit after {start:.0%} of a "
                         f"{n_total * dt_ps / 1000:.2f} ns run (at most {fit}); lower --n, "
                         "--spacing-ns or --start")
    return sorted(idx)


def imaging(topology, anchor: str):
    """mdtraj topology and the anchor molecules for `read_frame`."""
    import mdtraj
    mdtop = mdtraj.Topology.from_openmm(topology)
    if anchor == "protein":
        chains = [c for c in mdtop.chains if any(r.name in AMINO_ACIDS for r in c.residues)]
        return mdtop, [set(c.atoms) for c in chains]
    return mdtop, [{mdtop.atom(find_ion(topology).index)}]


def read_frame(dcd: Path, index: int, mdtop, anchors):
    """Coordinates (A) and box vectors (nm) of one DCD frame, every molecule imaged about the
    anchor: the protein chains for a site frame, the ion for a water frame. ASH applies no
    periodic images to the QM region or its point charges, so the ion, its shell and the
    environment around them must be whole and on one side of the box (as umbrella `energies`
    does before its QM/MM single points)."""
    import mdtraj
    traj = mdtraj.load_frame(str(dcd), index, top=mdtop)
    traj.image_molecules(inplace=True, anchor_molecules=anchors)
    return traj.xyz[0] * 10.0, traj.unitcell_vectors[0]


def read_pdb_frame(pdb: Path, mdtop, anchors):
    """`read_frame` for a single-frame PDB (an umbrella window's lastframe.pdb), on the
    topology of the window's start.pdb (which carries the bonds the imaging needs)."""
    import mdtraj
    from openmm import unit
    from openmm.app import PDBFile
    p = PDBFile(str(pdb))
    if p.topology.getNumAtoms() != mdtop.n_atoms:
        raise SystemExit(f"{pdb}: {p.topology.getNumAtoms()} atoms, start.pdb has {mdtop.n_atoms}")
    xyz = np.array(p.positions.value_in_unit(unit.nanometer), dtype=np.float32)[None]
    box = np.array(p.topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer), dtype=np.float32)[None]
    traj = mdtraj.Trajectory(xyz, mdtop)
    traj.unitcell_vectors = box
    traj.image_molecules(inplace=True, anchor_molecules=anchors)
    return traj.xyz[0] * 10.0, traj.unitcell_vectors[0]


def window_time_ns(window: Path):
    """Length of an umbrella window's run (cv.csv's last time), the time of its lastframe."""
    cv = window / "cv.csv"
    last = cv.read_text().strip().splitlines()[-1].split(",")[0] if cv.exists() else "0"
    return float(last) / 1000.0


def environment(topology, xyz, ion):
    """How bulk a reference frame is: nearest protein heavy atom, lipid P and the box faces
    of the frame as written (the non-periodic embedding ends there), from the ion, in Å."""
    r = np.linalg.norm(xyz - xyz[ion.index], axis=1)
    prot = [a.index for a in topology.atoms() if a.residue.name in AMINO_ACIDS and a.element is not None
            and a.element.symbol != "H"]
    lip_p = [a.index for a in topology.atoms() if a.residue.name in LIPIDS and a.element is not None
             and a.element.symbol == "P"]
    lo, hi = xyz.min(axis=0), xyz.max(axis=0)
    return {"protein_min_A": round(float(r[prot].min()), 1) if prot else None,
            "lipid_P_min_A": round(float(r[lip_p].min()), 1) if lip_p else None,
            "edge_min_A": round(float(min((xyz[ion.index] - lo).min(), (hi - xyz[ion.index]).min())), 1)}


def shell(topology, xyz, box_A, ion, cutoff: float, donors):
    """mcube's coordination rule (preprocess.coordination_shell): the closest donor of each
    residue within `cutoff` of the ion, minimum image. [(residue, donor atom, distance)]"""
    closest = {}
    for a in topology.atoms():
        if a.index == ion.index or a.element is None or a.element.symbol.upper() not in donors:
            continue
        d = xyz[a.index] - xyz[ion.index]
        d -= box_A * np.round(d / box_A)
        r = float(np.linalg.norm(d))
        if r <= cutoff and (a.residue not in closest or r < closest[a.residue][1]):
            closest[a.residue] = (a, r)
    return sorted(((res, a, r) for res, (a, r) in closest.items()), key=lambda t: t[2])


def qm_region(topology, xyz, box_A, ion, args):
    """`qm_residues` ({chain: {resid: name}}), the residue charge it carries besides the ion,
    and the coordination it came from. Lipids are refused (a whole lipid would enter QM);
    monatomic MM ions in reach stay MM by default (no Na/Cl in the lcecp basis files)."""
    rows = shell(topology, xyz, box_A, ion, args.qm_cutoff + args.qm_margin, set(args.donors))
    residues, charge, info, notes = {}, 0, [], []
    residues.setdefault(ion.residue.chain.id, {})[int(ion.residue.id)] = ion.residue.name
    for res, atom, r in rows:
        if res is ion.residue:
            continue
        row = {"residue": f"{res.chain.id}:{int(res.id)}:{res.name}", "donor": atom.name,
               "distance_A": round(r, 2)}
        if res.name in LIPIDS:
            raise SystemExit(f"lipid {row['residue']} {atom.name} is {r:.2f} A from the ion; a whole "
                             "lipid would go into QM. Lower --qm-cutoff or treat this frame by hand")
        if len(list(res.atoms())) == 1 and not args.qm_counterions:
            notes.append(f"{row['residue']} at {r:.2f} A left in MM (--qm-counterions to include)")
            continue
        if res.name in AMINO_ACIDS and atom.name in BACKBONE:
            notes.append(f"{row['residue']} coordinates through backbone {atom.name}, which mcube's "
                         "CA-CB cut leaves in MM")
        residues.setdefault(res.chain.id, {})[int(res.id)] = res.name
        charge += RES_CHARGE.get(res.name, 0)
        info.append(row)
    return residues, charge, info, notes


def qm_elements(topology, residues):
    """Elements of the QM atoms as mcube builds the region (amino acids without backbone), plus
    H for link atoms."""
    out = {"H"}
    for res in topology.residues():
        if int(res.id) in residues.get(res.chain.id, {}):
            for a in res.atoms():
                if res.name not in AMINO_ACIDS or a.name not in BACKBONE:
                    out.add(a.element.symbol)
    return out


def basis_elements(path):
    """Elements with a basis block in an NWChem-format file (None if the file is not readable
    here, e.g. `frames` run off Bunya)."""
    import re
    p = Path(path)
    if not p.is_file():
        return None
    return set(re.findall(r"^\s*([A-Z][a-z]?)\s+[SPDFGHI]+\s*$", p.read_text(), re.M))


def swap(topology, ion, target):
    """Rename the ion residue/atom and set its element: the only change between the three ions'
    frames. Returns a function that undoes it."""
    from openmm.app import Element
    old = (ion.residue.name, ion.name, ion.element)
    ion.residue.name, ion.name = target, target
    ion.element = Element.getBySymbol(IONS[target][0])

    def undo():
        ion.residue.name, ion.name, ion.element = old
    return undo


def write_pdb(topology, xyz_A, box_nm, path):
    """Residues renumbered 1..N per chain, chains relettered in order (OpenMM's writer without
    keepIds, as ASH writes): what preprocess's numbering check requires."""
    from openmm import Vec3, unit
    from openmm.app import PDBFile
    topology.setPeriodicBoxVectors(unit.Quantity([Vec3(*map(float, v)) for v in box_nm], unit.nanometer))
    with open(path, "w") as fh:
        PDBFile.writeFile(topology, np.asarray(xyz_A) / 10.0 * unit.nanometer, fh)


def job_config(args, kind, ion_name, chain_ion, residues, charge, forcefields):
    """The job's config.yaml: preprocess -> qmmm (optimisation only) -> sp, all on frame.pdb."""
    cutoff = round(args.qm_cutoff + args.qm_margin, 3)
    svp, tzvp = (Path(b).resolve().as_posix() for b in (args.basis_svp, args.basis_tzvp))
    opt = {"qm_engine": args.opt_engine}
    if args.opt_engine == "dft":
        opt.update({"dft_functional": args.functional, "dft_basis": "def2-SVP/lcecp-1-SVP",
                    "dft_basis_file": svp, "dft_ecp_file": svp,
                    "dft_cores": args.cores})
    else:
        opt.update({"xtb_method": "GFN2"})
    return {
        "preprocess": {"inputs": {
            # no solvation and no neutralisation: the frame's MM environment is the same for
            # every ion; the detection here only reports (qmmm.inputs.qm_residues decides)
            "protein_pdb": "frame.pdb", "skip_solvate": True, "metals": [ion_name],
            "cutoff": cutoff, "donor_elements": list(args.donors),
            "default_forcefields": list(forcefields)}},
        "qmmm": {"inputs": {
            "pdb_file": "frame.pdb", "forcefields": list(forcefields),
            "qm_residues": residues, "charge": int(charge), "mult": 1,
            "metal_parameters": str(args.metal_parameters), "cores": args.cores,
            # Deliberate exception to mcube's "QM/MM MD before optimisation": MD of each
            # swapped system would let the three environments drift apart and destroy the
            # pairing that makes E(M2) - E(M1) on one frame meaningful. So: local optimisation
            # from the shared MM frame only.
            "skip_equilibrate": True,
            "skip_optimize": False,
            "actradius": float(args.actradius),
            "actregion_origin": f"{chain_ion}:{ion_name}",
            "opt_maxiter": int(args.opt_maxiter),
            **opt}},
        "sp": {"inputs": {
            "pdb_file": "QM_MM_optimized.pdb",
            "sp_functional": args.functional, "sp_basis": "def2-TZVP/lcecp-1-TZVP",
            "sp_basis_file": tzvp, "sp_ecp_file": tzvp,
            "sp_cores": args.cores, "cubes": False}},
    }


HEADER = """\
# exchange.py ({kind} {src} frame {frame}, ion {ion}): QM/MM vertical swap.
# qmmm.skip_equilibrate is true on purpose: no QM/MM MD, only a local optimisation from the MM
# frame every ion shares, so the environments stay paired (docs/channel/qmmm/README.md).
# Comments are lost when a stage first writes its results back.
"""


def parse_sources(items, what):
    """{label: (ion, DIR)} from ION=DIR or ION@LABEL=DIR. The label names the source in frame
    keys and in the analysis, so two runs with the same ion (another C4 set, a Ca seeded at a
    La site) stay apart; it defaults to the ion."""
    out = {}
    for item in items or []:
        key, _, path = item.partition("=")
        ion, _, label = key.partition("@")
        if ion not in IONS or not path or (label and not label.isalnum()):
            raise SystemExit(f"--{what} wants ION=DIR or ION@LABEL=DIR with ION in {sorted(IONS)} "
                             f"and an alphanumeric LABEL, got {item!r}")
        out[label or ion] = (ion, Path(path))
    return out


def frames(args):
    import yaml
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    have = {"svp": basis_elements(args.basis_svp), "tzvp": basis_elements(args.basis_tzvp)}
    bonds = {k: p for k, (_, p) in parse_sources(args.bonds_from, "bonds-from").items()}
    plan = []
    for kind, sources, top_name, dcd_name, n, spacing, start, ff, anchor in (
            ("site", parse_sources(args.site, "site"), "start.pdb", "window.dcd", args.n,
             args.spacing_ns, args.start, args.site_forcefields or SITE_FF, "protein"),
            # a bulk box holding a ligand bound to the ion (bulk_ion.py --ligand-from): read like a
            # water box, analysed as a site (the calibration of docs/channel/qmmm/README.md)
            ("site", parse_sources(args.complex, "complex"), "system.pdb", args.water_dcd,
             args.water_n, args.water_spacing_ns, args.water_start,
             args.site_forcefields or SITE_FF, "ion"),
            ("water", parse_sources(args.water, "water"), "system.pdb", args.water_dcd,
             args.water_n, args.water_spacing_ns, args.water_start,
             args.water_forcefields or WATER_FF, "ion")):
        for src, (src_ion, path) in sources.items():
            plan.append((kind, src, src_ion, path, top_name, dcd_name, n, spacing, start, ff, anchor))
    # The water reference inside the membrane box itself (decision log C22, second round): the
    # last frames of the PMF's outer umbrella windows, ion in the bulk water layer. Same periodic
    # system as the site frames, so the net-charge and Ewald-convention terms that a separate
    # water cube does not share are common to site and reference.
    # A single window with window.dcd (a bulk umbrella window, swap_seed.py --move-z) is read
    # like a site run with the water frame picks; DIR1:DIR2:... takes each window's lastframe.pdb.
    bulk = parse_sources(args.bulk, "bulk")
    for src, (src_ion, path) in bulk.items():
        dcd = "window.dcd" if (path / "window.dcd").exists() else None
        plan.append(("water", src, src_ion, path, "start.pdb", dcd, args.water_n, args.water_spacing_ns,
                     args.water_start, args.site_forcefields or SITE_FF, "ion"))
    if not plan:
        raise SystemExit("nothing to do: give --site, --complex, --water and/or --bulk")

    meta = json.loads((out / "frames.json").read_text()) if (out / "frames.json").exists() else {}
    lists = {"site": [], "water": []}
    for kind, src, src_ion, path, top_name, dcd_name, n, spacing, start, ff, anchor in plan:
        if dcd_name is None:
            # --bulk: DIR is DIR1:DIR2:..., one frame (lastframe.pdb) per umbrella window
            windows = [Path(p) for p in str(path).split(":")]
            path = windows[0].parent
            topology = load_topology(windows[0] / "start.pdb")
            mdtop, anchors = imaging(topology, anchor)
            frames_of = [(int(w.name.rsplit("_", 1)[-1]), lambda w=w: read_pdb_frame(w / "lastframe.pdb", mdtop, anchors),
                          window_time_ns(w)) for w in windows]
            print(f"[{kind} {src}] {len(windows)} umbrella windows: "
                  + ", ".join(w.name for w in windows), flush=True)
        else:
            from mdtraj.formats import DCDTrajectoryFile
            dcd = path / dcd_name
            topology = load_topology(path / top_name, bonds.get(src) if kind == "site" else None)
            with DCDTrajectoryFile(str(dcd)) as fh:
                n_total = len(fh)
            box = top_name == "system.pdb"   # a bulk_ion.py box: water, or a --complex
            dt = frame_spacing_ps(path, "water" if box else "site", n_total,
                                  args.water_frame_ps if box else args.frame_ps)
            mdtop, anchors = imaging(topology, anchor)
            chosen = pick_frames(n_total, dt, n, spacing, start)
            print(f"[{kind} {src}] {n_total} frames, {dt:g} ps apart: frames {chosen} "
                  f"({', '.join(f'{(i + 1) * dt / 1000:.2f}' for i in chosen)} ns)", flush=True)
            frames_of = [(i, lambda i=i: read_frame(dcd, i, mdtop, anchors), (i + 1) * dt / 1000)
                         for i in chosen]
        ion = find_ion(topology)
        if ion.residue.name != src_ion:
            raise SystemExit(f"{path} carries {ion.residue.name}, not {src_ion}")
        for index, read, t_ns in frames_of:
            xyz, box_nm = read()
            box_A = np.diag(box_nm) * 10.0
            residues, res_charge, coord, notes = qm_region(topology, xyz, box_A, ion, args)
            # ASH sees raw positions: every QM atom must sit next to the ion as written
            reach = max(float(np.linalg.norm(xyz[a.index] - xyz[ion.index]))
                        for res in topology.residues() if int(res.id) in residues.get(res.chain.id, {})
                        for a in res.atoms())
            if reach > args.qm_cutoff + args.qm_margin + 6.0:
                raise SystemExit(f"{kind} {src} frame {index}: a QM residue reaches {reach:.1f} A "
                                 "from the ion as written; the imaging failed")
            elems = qm_elements(topology, residues)
            for label, got in have.items():
                if got is not None and not elems <= got:
                    raise SystemExit(f"{kind} {src} frame {index}: QM elements {sorted(elems - got)} "
                                     f"have no basis in the {label} file")
            frame_key = f"{src}_f{index:04d}"
            entry = {"kind": kind, "source": src, "source_ion": src_ion,
                     "source_dir": str(path.resolve()),
                     "frame": index, "time_ns": round(t_ns, 3),
                     "qm_residues": residues, "coordination": coord, "notes": notes,
                     "residue_charge": res_charge, "qm_elements": sorted(elems), "jobs": {}}
            if kind == "water" and src in bulk:
                entry["environment"] = environment(topology, xyz, ion)
                print(f"  frame {index}: " + ", ".join(f"{k} {v}" for k, v in entry["environment"].items()))
                if entry["environment"]["protein_min_A"] < args.bulk_min_A:
                    raise SystemExit(f"{src} window {index}: protein {entry['environment']['protein_min_A']} A "
                                     f"from the ion, under --bulk-min {args.bulk_min_A}; not bulk")
            for target in args.ions:
                job = out / kind / frame_key / target
                job.mkdir(parents=True, exist_ok=True)
                lists[kind].append(job.as_posix())
                entry["jobs"][target] = job.as_posix()
                cfg = job / "config.yaml"
                if cfg.exists() and not args.force:
                    continue                 # a job may already have results in it
                undo = swap(topology, ion, target)
                try:
                    write_pdb(topology, xyz, box_nm, job / "frame.pdb")
                finally:
                    undo()
                # chain and resid as the written file has them (OpenMM's writer reletters/renumbers)
                from openmm.app import PDBFile
                written = PDBFile(str(job / "frame.pdb")).topology
                w_ion = find_ion(written)
                w_res = qm_region_from_written(written, topology, residues)
                config = job_config(args, kind, target, f"{w_ion.residue.chain.id}:{int(w_ion.residue.id)}",
                                    w_res, IONS[target][1] + res_charge, ff)
                cfg.write_text(HEADER.format(kind=kind, src=src, frame=index, ion=target)
                               + yaml.safe_dump(config, sort_keys=False))
            meta[f"{kind}/{frame_key}"] = entry
            print(f"  frame {index}: QM {sum(len(v) for v in residues.values())} residues, "
                  f"charge {res_charge:+d} + ion; " + ", ".join(
                      f"{c['residue']}:{c['donor']} {c['distance_A']}" for c in coord), flush=True)
            for note in notes:
                print(f"  [NOTE] {note}")
    (out / "frames.json").write_text(json.dumps(meta, indent=1))
    for kind, jobs in lists.items():
        if jobs:
            path = out / f"jobs_{kind}.txt"
            old = path.read_text().split() if path.exists() else []
            path.write_text("\n".join(dict.fromkeys(old + jobs)) + "\n")
            print(f"{path}: {len(dict.fromkeys(old + jobs))} jobs")


def qm_region_from_written(written, source, residues):
    """`residues` (keys from the source topology) mapped onto the written file's chain letters
    and residue numbers, by residue order; names as written (the ion's is the target's)."""
    src_res, out_res = list(source.residues()), list(written.residues())
    if len(src_res) != len(out_res):
        raise SystemExit("the written frame has a different residue count")
    out = {}
    for a, b in zip(src_res, out_res):
        if int(a.id) in residues.get(a.chain.id, {}):
            out.setdefault(b.chain.id, {})[int(b.id)] = b.name
    return out


# ==============================
# Analysis
# ==============================

def job_energies(job: Path):
    import yaml
    cfg = job / "config.yaml"
    if not cfg.exists():
        return None
    data = yaml.safe_load(cfg.read_text()) or {}
    q = (data.get("qmmm") or {}).get("results") or {}
    s = (data.get("sp") or {}).get("results") or {}
    p = (data.get("preprocess") or {}).get("results") or {}
    given = ((data.get("qmmm") or {}).get("inputs") or {}).get("qm_residues") or {}
    detected = p.get("qm_residues")
    same = None
    if detected is not None:
        norm = lambda d: {str(c): sorted(int(r) for r in v) for c, v in d.items() if v}
        same = norm(detected) == norm(given)
    return {"sp": s.get("sp_energy"), "sp_disp": s.get("sp_energy_disp"), "opt": q.get("opt_energy"),
            "converged": q.get("converged"), "preprocess_agrees": same}


def stats(values):
    v = np.asarray(values, dtype=float)
    if len(v) == 0:
        return None
    sd = float(v.std(ddof=1)) if len(v) > 1 else None
    return {"n": int(len(v)), "mean": round(float(v.mean()), 2),
            "sd": round(sd, 2) if sd is not None else None,
            "sem": round(sd / math.sqrt(len(v)), 2) if sd is not None else None}


def analyse(args):
    out = Path(args.out).resolve()
    meta = json.loads((out / "frames.json").read_text())
    level_names = {"sp": "r2SCAN/def2-TZVP(lcecp-1) QM/MM single point",
                   "opt": "QM/MM optimisation energy (opt level, def2-SVP by default)"}
    report = {"note": "kcal/mol; dE = E(M2) - E(M1) on one frame; ddE = <dE>_site - <dE>_water; "
                      "negative = the site prefers M2. Gas-phase QM with electrostatic embedding; "
                      "no thermal terms. lra = linear response, the mean of <dE> sampled with M1 "
                      "and with M2 (site and water alike).", "frames": {}, "levels": {}}
    per = {}       # (level, kind, source, pair) -> [dE]
    src_ion = {}   # (kind, source) -> the ion the source was sampled with
    levels = [l for l in ("sp", "sp_disp", "opt")
              if any((job_energies(Path(p)) or {}).get(l) is not None
                     for entry in meta.values() for p in entry["jobs"].values())]
    level_names["sp_disp"] = level_names["sp"] + " + D3(BJ) on the QM region"
    for key, entry in sorted(meta.items()):
        src_ion[(entry["kind"], entry["source"])] = entry.get("source_ion", entry["source"])
        e = {ion: job_energies(Path(p)) for ion, p in entry["jobs"].items()}
        row = {"time_ns": entry["time_ns"], "qm_residues": entry["qm_residues"],
               "converged": {i: (v or {}).get("converged") for i, v in e.items()},
               "preprocess_agrees": {i: (v or {}).get("preprocess_agrees") for i, v in e.items()},
               "dE": {}}
        for level in levels:
            for m1, m2 in EXCHANGES:
                a, b = e.get(m1), e.get(m2)
                ok = a and b and a[level] is not None and b[level] is not None
                if ok and not args.include_unconverged:
                    ok = a["converged"] is not False and b["converged"] is not False
                if not ok:
                    continue
                d = (b[level] - a[level]) * HARTREE_KCAL
                row["dE"].setdefault(level, {})[f"{m1}->{m2}"] = round(d, 2)
                per.setdefault((level, entry["kind"], entry["source"], f"{m1}->{m2}"), []).append(d)
        report["frames"][key] = row
    if args.embed_radius or args.periodic:
        # First order in the far field: the QM region's charge changes by exactly
        # dq = q(M2) - q(M1), so a change of the MM potential at the ion by d(phi) changes that
        # frame's dE by dq * d(phi) (potential.json).
        #   --embed-radius R: drop the MM charges beyond R (d(phi) = -phi_R), the same finite
        #     embedding for site and water;
        #   --periodic: the PME potential in place of the non-periodic direct sum
        #     (d(phi) = phi_pme - phi_direct), the electrostatics the MD ran with.
        pot = json.loads((out / "potential.json").read_text())["frames"]
        R = f"{args.embed_radius:g}" if args.embed_radius else None
        missing = sorted(key for key in meta if key not in pot or (args.periodic and "phi_pme_V" not in pot[key]))
        if missing:
            raise SystemExit(f"potential.json lacks {len(missing)} frames (run `potential` again): {missing[:3]}")
        report["embed_radius_A" if R else "periodic"] = args.embed_radius or True
        for key, entry in meta.items():
            p = pot[key]
            phi = p["phi_V"][R] if R else p["phi_direct_V"] - p["phi_pme_V"]
            for level, pairs in report["frames"][key]["dE"].items():
                for pair, d in pairs.items():
                    m1, m2 = pair.split("->")
                    pairs[pair] = round(d - EV_KCAL * (IONS[m2][1] - IONS[m1][1]) * phi, 2)
        per = {}
        for key, entry in meta.items():
            for level, pairs in report["frames"][key]["dE"].items():
                for pair, d in pairs.items():
                    per.setdefault((level, entry["kind"], entry["source"], pair), []).append(d)
    lra_sets = {}
    for item in args.lra or []:
        name, _, members = item.partition("=")
        lra_sets[name] = members.split("+")
    for level in levels:
        lv = {"what": level_names[level], "site": {}, "water": {}, "ddE": {}, "lra": {}}
        for m1, m2 in EXCHANGES:
            pair = f"{m1}->{m2}"
            pooled = {}
            for kind in ("site", "water"):
                rows = {src: stats(v) for (l, k, src, p), v in per.items()
                        if l == level and k == kind and p == pair}
                allv = [d for (l, k, src, p), v in per.items()
                        if l == level and k == kind and p == pair for d in v]
                pooled[kind] = stats(allv)
                # frames of one source are correlated (contacts that never break), so the spread
                # of the per-source means is the honest error of the pooled mean
                means = [st["mean"] for st in rows.values() if st]
                lv[kind][pair] = {"by_source": rows, "all": pooled[kind],
                                  "sd_between_sources": round(float(np.std(means, ddof=1)), 2)
                                  if len(means) > 1 else None}
            dd = {}
            if pooled["site"] and pooled["water"]:
                dd["all"] = ddE(pooled["site"], pooled["water"])
                for src, st in lv["site"][pair]["by_source"].items():
                    # a site source is referenced to the water sampled with the same ion
                    w = water_by_ion(lv["water"][pair]["by_source"], src_ion, src_ion.get(("site", src)))
                    if st and w:
                        dd[f"{src} frames"] = ddE(st, w)
            lv["ddE"][pair] = dd
            for name, members in lra_sets.items():
                ions = [src_ion.get(("site", s)) for s in members]
                if sorted(ions) != sorted([m1, m2]):
                    continue
                site = [lv["site"][pair]["by_source"].get(s) for s in members]
                water = [water_by_ion(lv["water"][pair]["by_source"], src_ion, i) for i in ions]
                if all(site) and all(water):
                    lv["lra"][f"{pair} {name}"] = lra(site, water)
        report["levels"][level] = lv
    name = (f"report_R{args.embed_radius:g}.json" if args.embed_radius
            else "report_periodic.json" if args.periodic else "report.json")
    (out / name).write_text(json.dumps(report, indent=1))
    print_report(report)
    print(f"\n{out / name}")


def water_by_ion(by_source, src_ion, ion):
    """The water statistics of the sources sampled with `ion`, pooled (weighted by n)."""
    rows = [st for src, st in by_source.items() if st and src_ion.get(("water", src)) == ion]
    if not rows:
        return None
    if len(rows) == 1:
        return rows[0]
    n = sum(r["n"] for r in rows)
    mean = sum(r["mean"] * r["n"] for r in rows) / n
    sem = math.sqrt(sum((r["sem"] or 0) ** 2 * r["n"] ** 2 for r in rows)) / n
    return {"n": n, "mean": round(mean, 2), "sd": None, "sem": round(sem, 2)}


def lra(site, water):
    """Linear-response exchange energy: the mean of the two endpoint ensembles' <dE>, at the
    site and in water, so the environment's reorganisation to each ion enters (Warshel's LRA)."""
    half = lambda a, b: (0.5 * (a["mean"] + b["mean"]),
                         0.5 * math.hypot(a["sem"] or 0, b["sem"] or 0))
    (s, se), (w, we) = half(*site), half(*water)
    return {"ddE": round(s - w, 2), "err": round(math.hypot(se, we), 2),
            "site_endpoints": [x["mean"] for x in site], "water_endpoints": [x["mean"] for x in water],
            "site": round(s, 2), "water": round(w, 2)}


def ddE(site, water):
    err = None
    if site["sem"] is not None and water["sem"] is not None:
        err = round(math.hypot(site["sem"], water["sem"]), 2)
    return {"ddE": round(site["mean"] - water["mean"], 2), "err": err,
            "n_site": site["n"], "n_water": water["n"]}


def read_xyz_frames(path: Path):
    """[(comment, elements, coords Å)] of a multi-frame xyz (geomeTRIC's optimisation trajectory,
    comment "Iteration N Energy E")."""
    lines = path.read_text().splitlines()
    out, i = [], 0
    while i < len(lines) and lines[i].strip():
        n = int(lines[i].split()[0])
        rows = [l.split() for l in lines[i + 2:i + 2 + n]]
        out.append((lines[i + 1], [r[0] for r in rows],
                    np.array([[float(x) for x in r[1:4]] for r in rows])))
        i += 2 + n
    return out


def last_opt_energies(path: Path):
    """QM/MM energies (Eh) of the last run in optimization_energies.log, which ASH appends to on
    every run, each run under its own 'Iteration ...' header."""
    vals = []
    for l in path.read_text().splitlines():
        t = l.split()
        if not t:
            continue
        if t[0] == "Iteration":
            vals = []
        elif t[0].isdigit():
            vals.append(float(t[-1]))
    return vals


def shell_geometry(elements, xyz, cutoff):
    """Ion-O shell of one QM-region geometry: O within `cutoff` of the ion, split into water O
    (two H within 1.2 Å) and other O (carboxylate, carbonyl, hydroxyl)."""
    metal = [k for k, e in enumerate(elements) if e.capitalize() in ("Ca", "La", "Lu")]
    if len(metal) != 1:
        return None
    m = xyz[metal[0]]
    hyd = np.array([k for k, e in enumerate(elements) if e.upper() == "H"], dtype=int)
    water, other = [], []
    for k, e in enumerate(elements):
        if e.upper() != "O":
            continue
        r = float(np.linalg.norm(xyz[k] - m))
        if r > cutoff:
            continue
        nh = int((np.linalg.norm(xyz[hyd] - xyz[k], axis=1) < 1.2).sum()) if len(hyd) else 0
        (water if nh >= 2 else other).append(r)
    allr = water + other
    return {"cn": len(allr), "cn_water": len(water), "cn_other": len(other),
            "mean_MO": round(float(np.mean(allr)), 3) if allr else None,
            "min_MO": round(float(min(allr)), 3) if allr else None}


def structure(args):
    """What the optimisations did to each ion's shell, from the QM-region trajectory (only QM
    atoms move, so it holds the whole path): steps, energy drop, QM RMSD, and the ion-O shell of
    the MM frame (first step) and of the optimised geometry (last step)."""
    out = Path(args.out).resolve()
    meta = json.loads((out / "frames.json").read_text())
    rows, per = {}, {}
    for key, entry in sorted(meta.items()):
        for ion, p in entry["jobs"].items():
            job = Path(p)
            traj = job / "geometric_OPTtraj_QMregion.xyz"
            if not traj.exists():
                rows[f"{key}/{ion}"] = {"missing": traj.name}
                continue
            # The trajectory may also hold steps of an earlier, interrupted run; every run starts
            # from the same MM frame, so the first frame is the start either way.
            fr = read_xyz_frames(traj)
            (_, el0, x0), (c1, el1, x1) = fr[0], fr[-1]
            t = c1.split()
            steps = int(t[1]) if len(t) > 1 and t[0] == "Iteration" and t[1].isdigit() else len(fr) - 1
            e = job_energies(job) or {}
            en = job / "optimization_energies.log"
            vals = last_opt_energies(en) if en.exists() else []
            drop = round((vals[-1] - vals[0]) * HARTREE_KCAL, 2) if len(vals) > 1 else None
            r = {"kind": entry["kind"], "source": entry["source"], "ion": ion,
                 "steps": steps, "converged": e.get("converged"), "opt_drop_kcal": drop,
                 "qm_rmsd_A": round(float(np.sqrt(((x1 - x0) ** 2).sum(1).mean())), 3),
                 "start": shell_geometry(el0, x0, args.cutoff),
                 "end": shell_geometry(el1, x1, args.cutoff)}
            rows[f"{key}/{ion}"] = r
            if r["end"]:
                per.setdefault((entry["kind"], entry["source"], ion), []).append(r)
    summary = {}
    for (kind, src, ion), rs in sorted(per.items()):
        f = lambda sel: round(float(np.mean([sel(x) for x in rs])), 2)
        summary[f"{kind} {src} frames {ion}"] = {
            "n": len(rs), "steps": f(lambda x: x["steps"]), "qm_rmsd_A": f(lambda x: x["qm_rmsd_A"]),
            "cn_start": f(lambda x: x["start"]["cn"]), "cn_end": f(lambda x: x["end"]["cn"]),
            "cn_water_end": f(lambda x: x["end"]["cn_water"]),
            "cn_other_end": f(lambda x: x["end"]["cn_other"]),
            "mean_MO_end": f(lambda x: x["end"]["mean_MO"]),
            "cn_changed": sum(x["start"]["cn"] != x["end"]["cn"] for x in rs)}
    res = {"note": f"ion-O shell within {args.cutoff} Å; water O = O with two H within 1.2 Å",
           "jobs": rows, "summary": summary}
    (out / "structure.json").write_text(json.dumps(res, indent=1))
    for k, s in summary.items():
        print(f"{k:24s} n {s['n']}  steps {s['steps']:6.1f}  rmsd {s['qm_rmsd_A']:.2f} Å  "
              f"CN {s['cn_start']:.2f} -> {s['cn_end']:.2f} (water {s['cn_water_end']:.2f}, "
              f"other {s['cn_other_end']:.2f}; changed in {s['cn_changed']})  "
              f"<M-O> {s['mean_MO_end']:.3f} Å")
    missing = [k for k, r in rows.items() if "missing" in r]
    if missing:
        print("no QM-region trajectory: " + ", ".join(missing))
    print(f"\n{out / 'structure.json'}")


# Bulky files a finished job no longer needs: the full-system optimisation trajectory (only QM
# atoms move, so geometric_OPTtraj_QMregion.xyz holds the same path), the last-step geometry dumps,
# and geomeTRIC's HT.v_<eigenvalue>.txt Hessian-eigenvector dumps (one per negative eigenvalue).
CLEAN = ("geometric_OPTtraj_Full.xyz", "Fragment-currentgeo.xyz", "Fragment-currentgeo.ygg",
         "HT.v_*.txt")


def clean(args):
    """Delete CLEAN files in jobs that are final (a TZVP energy, or an optimisation that did not
    converge) and that `structure` has already read. Dry run unless --delete."""
    out = Path(args.out).resolve()
    meta = json.loads((out / "frames.json").read_text())
    sfile = out / "structure.json"
    read = json.loads(sfile.read_text())["jobs"] if sfile.exists() else {}
    total, kept = 0, []
    for key, entry in sorted(meta.items()):
        for ion, p in entry["jobs"].items():
            job, name = Path(p), f"{key}/{ion}"
            e = job_energies(job) or {}
            final = e.get("sp") is not None or e.get("converged") is False
            if not final or "steps" not in read.get(name, {}):
                kept.append(name)
                continue
            for pat in CLEAN:
                for f in job.glob(pat):
                    total += f.stat().st_size
                    if args.delete:
                        f.unlink()
    print(f"{'deleted' if args.delete else 'would delete'} {total / 1e9:.2f} GB")
    if kept:
        print("left alone (not final, or not in structure.json): " + ", ".join(kept))


def variant(args):
    """Copy finished jobs (optimised geometry and config) to --dest with new `sp` inputs, so a
    different single point (another functional, D3) runs on the same geometries without touching
    the originals. Writes dest/frames.json and dest/jobs_<kind>.txt; exchange.slurm then runs only
    the sp stage (preprocess and qmmm results are in the copied config)."""
    import shutil
    import yaml
    out, dest = Path(args.out).resolve(), Path(args.dest).resolve()
    meta = json.loads((out / "frames.json").read_text())
    sets = {}
    for item in args.set or []:
        k, _, v = item.partition("=")
        sets[k] = yaml.safe_load(v)
    keep = set(args.sources or [])
    taken, new_meta, lists = {}, {}, {"site": [], "water": []}
    for key, entry in sorted(meta.items()):
        if keep and entry["source"] not in keep:
            continue
        group = (entry["kind"], entry["source"])
        if args.frames_per_source and taken.get(group, 0) >= args.frames_per_source:
            continue
        jobs = {}
        for ion, p in entry["jobs"].items():
            src = Path(p)
            e = job_energies(src) or {}
            if e.get("sp") is None:
                continue
            job = dest / entry["kind"] / Path(key).name / ion
            job.mkdir(parents=True, exist_ok=True)
            for name in ("QM_MM_optimized.pdb", "frame.pdb"):
                shutil.copy2(src / name, job / name)
            cfg = yaml.safe_load((src / "config.yaml").read_text())
            cfg["sp"]["inputs"].update(sets)
            cfg["sp"].pop("results", None)
            (job / "config.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
            jobs[ion] = job.as_posix()
            lists[entry["kind"]].append(job.as_posix())
        if len(jobs) == len(entry["jobs"]):
            taken[group] = taken.get(group, 0) + 1
            new_meta[key] = {**entry, "jobs": jobs}
    (dest / "frames.json").write_text(json.dumps(new_meta, indent=1))
    for kind, jobs in lists.items():
        if jobs:
            (dest / f"jobs_{kind}.txt").write_text("\n".join(jobs) + "\n")
            print(f"{dest / f'jobs_{kind}.txt'}: {len(jobs)} jobs")
    print(f"sp inputs set: {sets}; {len(new_meta)} frames")


COULOMB_V_A = 14.399645       # V Å per e: the potential of 1 e at 1 Å
EV_KCAL = 23.060548           # kcal/mol per (e x V)


def potential(args):
    """Far-field MM potential at the ion, as the non-periodic embedding sees it (raw positions,
    every MM point charge, no images): phi(R) = sum over MM atoms beyond R of q / r. The QM atoms
    all sit within ~6 Å, so beyond R >= 8 Å this is exactly what ASH embeds. A 2+ -> 3+ swap
    picks up Δq·phi, so a site-vs-water difference of the far field that is not the protein's
    own (membrane dipole, slab and cube surfaces, net charge of the embedded set) enters ddE at
    23 kcal/mol per volt (Rocklin et al., J. Chem. Phys. 2013, 139, 184103)."""
    import yaml
    from openmm import NonbondedForce, unit
    from openmm.app import ForceField, Modeller, NoCutoff, PDBFile
    out = Path(args.out).resolve()
    meta = json.loads((out / "frames.json").read_text())
    radii = [float(r) for r in args.radii]
    rows, groups = {}, {}
    for key, entry in sorted(meta.items()):
        ion = entry.get("source_ion", entry["source"])
        job = Path(entry["jobs"].get(ion) or next(iter(entry["jobs"].values())))
        cfg = yaml.safe_load((job / "config.yaml").read_text())
        pdb = PDBFile(str(job / "frame.pdb"))
        target = find_ion(pdb.topology)
        r_ion = np.array(pdb.positions[target.index].value_in_unit(unit.angstrom))
        mod = Modeller(pdb.topology, pdb.positions)
        mod.delete([target.residue])
        system = ForceField(*cfg["qmmm"]["inputs"]["forcefields"]).createSystem(
            mod.topology, nonbondedMethod=NoCutoff)
        nb = next(f for f in system.getForces() if isinstance(f, NonbondedForce))
        q = np.array([nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
                      for i in range(nb.getNumParticles())])
        xyz = np.array(mod.positions.value_in_unit(unit.angstrom))
        r = np.linalg.norm(xyz - r_ion, axis=1)
        # Whole residues in or out, by their centre: an atom-based cut splits waters at the
        # boundary into bare partial charges and swings phi by volts from frame to frame
        centre = np.empty(len(q))
        for res in mod.topology.residues():
            idx = [a.index for a in res.atoms()]
            centre[idx] = np.linalg.norm(xyz[idx].mean(axis=0) - r_ion)
        comp = np.array([("water" if a.residue.name in WATER else "lipid" if a.residue.name in LIPIDS
                          else "protein" if a.residue.name in AMINO_ACIDS else "ions")
                         for a in mod.topology.atoms()])
        phi = lambda sel: round(float(COULOMB_V_A * (q[sel] / r[sel]).sum()), 4)
        row = {"kind": entry["kind"], "source": entry["source"], "net_charge": round(float(q.sum()), 3),
               "r_max_A": round(float(r.max()), 1),
               "phi_V": {f"{R:g}": phi(centre > R) for R in radii},
               "by_component_V": {c: {f"{R:g}": phi((centre > R) & (comp == c)) for R in radii}
                                  for c in sorted(set(comp))}}
        if args.pme:
            row["phi_direct_V"] = phi(np.ones(len(q), bool))
            row["phi_pme_V"] = round(pme_potential(mod, cfg["qmmm"]["inputs"]["forcefields"], r_ion), 4)
        rows[key] = row
        groups.setdefault((entry["kind"], entry["source"]), []).append(row)
    summary = {}
    for (kind, src), rs in sorted(groups.items()):
        summary[f"{kind} {src}"] = {
            "n": len(rs), "net_charge": rs[0]["net_charge"],
            "phi_V": {k: round(float(np.mean([x["phi_V"][k] for x in rs])), 4) for k in rs[0]["phi_V"]},
            "phi_sd_V": {k: round(float(np.std([x["phi_V"][k] for x in rs])), 4) for k in rs[0]["phi_V"]}}
        if args.pme:
            d = [x["phi_pme_V"] - x["phi_direct_V"] for x in rs]
            summary[f"{kind} {src}"]["pme_minus_direct_V"] = {"mean": round(float(np.mean(d)), 4),
                                                              "sd": round(float(np.std(d)), 4)}
    res = {"note": "phi(R): potential at the ion from MM charges beyond R (V), raw positions as "
                   "the non-periodic embedding sees them; a 1 e change of the ion's charge picks up "
                   f"{EV_KCAL:.2f} kcal/mol per V", "frames": rows, "summary": summary}
    (out / "potential.json").write_text(json.dumps(res, indent=1))
    print("phi(R) at the ion, V (mean over frames; sd)")
    print(f"{'':16s}" + "".join(f"{'R>' + k:>16s}" for k in summary[next(iter(summary))]["phi_V"]))
    for name, s in summary.items():
        print(f"{name:16s}" + "".join(f"{v:9.3f} ({s['phi_sd_V'][k]:.3f})" for k, v in s["phi_V"].items())
              + f"   net {s['net_charge']:+.2f} e"
              + (f"   PME - direct {s['pme_minus_direct_V']['mean']:+.3f} ({s['pme_minus_direct_V']['sd']:.3f}) V"
                 if "pme_minus_direct_V" in s else ""))
    print(f"\n{out / 'potential.json'}")


def pme_potential(mod, forcefields, r_ion_A, delta=1.0):
    """Potential (V) at r_ion from the MM charges under PME, the periodic electrostatics the MD
    ran with: a probe particle (no LJ) at the ion's position, E(+delta) - E(-delta) over 2 delta.
    The energy is exactly quadratic in the probe charge, so the central difference is exact at
    any delta; 1 e keeps it well above the single-precision noise of a 124k-atom energy. The
    probe's self and image terms are even in delta and cancel."""
    import openmm
    from openmm import unit
    from openmm.app import ForceField, PME
    system = ForceField(*forcefields).createSystem(mod.topology, nonbondedMethod=PME,
                                                   nonbondedCutoff=1.0 * unit.nanometer)
    nb = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    nb.setEwaldErrorTolerance(1e-5)
    probe = system.addParticle(0.0)
    nb.addParticle(0.0, 0.1, 0.0)
    pos = list(mod.positions.value_in_unit(unit.nanometer)) + [openmm.Vec3(*(np.asarray(r_ion_A) / 10.0))]
    ctx = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("CPU"))
    ctx.setPositions(pos)
    energies = []
    # The total energy: only the NonbondedForce sees the probe (no bonds, no other forces on it).
    # The first version read the NonbondedForce's force group alone, and getState(groups=...)
    # returned 0 on this OpenMM, so every frame read 0 V
    for q in (delta, -delta):
        nb.setParticleParameters(probe, q, 0.1, 0.0)
        nb.updateParametersInContext(ctx)
        energies.append(ctx.getState(getEnergy=True).getPotentialEnergy()
                        .value_in_unit(unit.kilojoule_per_mole))
    phi = (energies[0] - energies[1]) / (2 * delta) / 96.4853   # kJ/mol/e -> V
    if phi == 0.0:
        raise RuntimeError("PME probe saw no potential: check the box vectors of the frame")
    return phi


def print_report(report):
    for key, row in report["frames"].items():
        de = row["dE"].get("sp", {})
        bad = [i for i, c in row["converged"].items() if c is False]
        differ = [i for i, c in row["preprocess_agrees"].items() if c is False]
        print(f"{key:24s} {row['time_ns']:6.2f} ns  " + "  ".join(f"{p} {v:8.2f}" for p, v in de.items())
              + (f"  (opt not converged: {', '.join(bad)})" if bad else "")
              + (f"  (preprocess detected another shell: {', '.join(differ)})" if differ else ""))
    for level, lv in report["levels"].items():
        print(f"\n== {lv['what']}")
        for pair, dd in lv["ddE"].items():
            s, w = lv["site"][pair]["all"], lv["water"][pair]["all"]
            fmt = lambda st: f"{st['mean']:8.2f} sd {st['sd']} (n {st['n']})" if st else "      --"
            line = f"{pair:8s} site {fmt(s)}  water {fmt(w)}"
            if dd:
                line += "  ddE " + "; ".join(f"{k} {v['ddE']:+.2f} +- {v['err']}" for k, v in dd.items())
            sb = lv["site"][pair].get("sd_between_sources")
            if sb is not None:
                line += f"  (site sd between sources {sb})"
            print(line)
        for name, v in lv.get("lra", {}).items():
            print(f"  LRA {name}: ddE {v['ddE']:+.2f} +- {v['err']}  (site endpoints "
                  f"{v['site_endpoints']}, water endpoints {v['water_endpoints']})")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="action", required=True)
    f = sub.add_parser("frames")
    f.add_argument("--out", required=True)
    f.add_argument("--site", action="append", help="ION=umbrella window dir (start.pdb, window.dcd)")
    f.add_argument("--water", action="append", help="ION=bulk run dir (system.pdb, --water-dcd)")
    f.add_argument("--complex", action="append",
                   help="ION@LABEL=bulk run dir with a bound ligand (bulk_ion.py --ligand-from), "
                        "read like --water and analysed as a site")
    f.add_argument("--bulk", action="append",
                   help="ION@LABEL=WIN1:WIN2:... umbrella windows of the membrane box with the ion in bulk "
                        "water (one frame each, lastframe.pdb), read as the water reference in the same box")
    f.add_argument("--bulk-min", dest="bulk_min_A", type=float, default=8.0,
                   help="refuse a --bulk frame whose ion is nearer the protein than this (Å)")
    f.add_argument("--bonds-from", action="append",
                   help="ION=processed_pdb whose bonds a site start.pdb lacks (F10)")
    f.add_argument("--ions", nargs="+", default=list(IONS), choices=sorted(IONS))
    f.add_argument("--n", type=int, default=5, help="site frames per source")
    f.add_argument("--spacing-ns", type=float, default=1.5)
    f.add_argument("--start", type=float, default=0.4,
                   help="earliest frame as a fraction of the run (5 x 1.5 ns in 10 ns needs 0.4)")
    f.add_argument("--frame-ps", type=float, help="ps between site frames (default: from cv.csv)")
    f.add_argument("--water-dcd", default="frames.dcd")
    f.add_argument("--water-n", type=int, default=5)
    f.add_argument("--water-spacing-ns", type=float, default=0.25)
    f.add_argument("--water-start", type=float, default=0.5)
    f.add_argument("--water-frame-ps", type=float, default=None,
                   help="ps between water frames (bulk_ion.py --dcd-every)")
    f.add_argument("--qm-cutoff", type=float, default=3.3, help="A, ion-donor distance")
    f.add_argument("--qm-margin", type=float, default=0.0, help="A added to --qm-cutoff")
    f.add_argument("--donors", nargs="+", default=["N", "O", "S"])
    f.add_argument("--qm-counterions", action="store_true",
                   help="put Na+/Cl- in reach into QM (needs Na/Cl in the basis files)")
    f.add_argument("--actradius", type=float, default=8.0, help="A around the ion free to move")
    f.add_argument("--opt-engine", choices=("dft", "xtb"), default="dft")
    f.add_argument("--opt-maxiter", type=int, default=150)
    f.add_argument("--functional", default="r2scan")
    f.add_argument("--basis-svp", required=True)
    f.add_argument("--basis-tzvp", required=True)
    f.add_argument("--metal-parameters", default="1264",
                   help="MM set for the ion's LJ (its C4 is off: it is a QM atom)")
    f.add_argument("--site-forcefields", nargs="+")
    f.add_argument("--water-forcefields", nargs="+")
    f.add_argument("--cores", type=int, default=32)
    f.add_argument("--force", action="store_true", help="rewrite existing job configs")
    a = sub.add_parser("analyse")
    a.add_argument("--out", required=True)
    a.add_argument("--include-unconverged", action="store_true")
    a.add_argument("--lra", action="append",
                   help="NAME=SRC1+SRC2: two site sources sampled with the two ions of an exchange "
                        "(e.g. lasite=CAla+LA); their LRA ddE is reported per exchange")
    a.add_argument("--embed-radius", type=float,
                   help="Å: report ddE as if both site and water embedded only MM residues within this "
                        "radius (first order, from potential.json; writes report_R<radius>.json)")
    a.add_argument("--periodic", action="store_true",
                   help="report ddE with the PME potential at the ion in place of the non-periodic sum "
                        "(first order; potential --pme first; writes report_periodic.json)")
    s = sub.add_parser("structure")
    s.add_argument("--out", required=True)
    s.add_argument("--cutoff", type=float, default=3.1, help="Å, ion-O shell")
    v = sub.add_parser("variant")
    v.add_argument("--out", required=True)
    v.add_argument("--dest", required=True)
    v.add_argument("--set", action="append", help="KEY=VALUE for sp.inputs (YAML value)")
    v.add_argument("--sources", nargs="+", help="only these source labels")
    v.add_argument("--frames-per-source", type=int, help="only the first N frames of each source")
    pt = sub.add_parser("potential")
    pt.add_argument("--out", required=True)
    pt.add_argument("--radii", nargs="+", default=["8", "10", "15", "20", "30"])
    pt.add_argument("--pme", action="store_true",
                    help="also the PME potential at the ion (OpenMM, CPU), for `analyse --periodic`")
    c = sub.add_parser("clean")
    c.add_argument("--out", required=True)
    c.add_argument("--delete", action="store_true", help="delete (default: dry run)")
    bs = sub.add_parser("bulk-survey", help="environment of umbrella windows' last frames, to pick --bulk ones")
    bs.add_argument("windows", nargs="+")
    args = p.parse_args()
    if args.action == "bulk-survey":
        topology = load_topology(Path(args.windows[0]) / "start.pdb")
        mdtop, anchors = imaging(topology, "ion")
        ion = find_ion(topology)
        for w in map(Path, args.windows):
            xyz, _ = read_pdb_frame(w / "lastframe.pdb", mdtop, anchors)
            env = environment(topology, xyz, ion)
            print(f"{w.parent.name}/{w.name} " + " ".join(f"{k} {v}" for k, v in env.items()), flush=True)
    elif args.action == "frames":
        args.donors = [d.upper() for d in args.donors]
        frames(args)
    elif args.action == "structure":
        structure(args)
    elif args.action == "clean":
        clean(args)
    elif args.action == "variant":
        variant(args)
    elif args.action == "potential":
        potential(args)
    else:
        analyse(args)


if __name__ == "__main__":
    main()
