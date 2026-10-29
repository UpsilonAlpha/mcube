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
    out = {}
    for item in items or []:
        key, _, path = item.partition("=")
        if key not in IONS or not path:
            raise SystemExit(f"--{what} wants ION=DIR with ION in {sorted(IONS)}, got {item!r}")
        out[key] = Path(path)
    return out


def frames(args):
    import yaml
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    have = {"svp": basis_elements(args.basis_svp), "tzvp": basis_elements(args.basis_tzvp)}
    bonds = parse_sources(args.bonds_from, "bonds-from")
    plan = []
    for kind, sources, top_name, dcd_name, n, spacing, start, ff, anchor in (
            ("site", parse_sources(args.site, "site"), "start.pdb", "window.dcd", args.n,
             args.spacing_ns, args.start, args.site_forcefields or SITE_FF, "protein"),
            ("water", parse_sources(args.water, "water"), "system.pdb", args.water_dcd,
             args.water_n, args.water_spacing_ns, args.water_start,
             args.water_forcefields or WATER_FF, "ion")):
        for src, path in sources.items():
            plan.append((kind, src, path, top_name, dcd_name, n, spacing, start, ff, anchor))
    if not plan:
        raise SystemExit("nothing to do: give --site and/or --water")

    meta = json.loads((out / "frames.json").read_text()) if (out / "frames.json").exists() else {}
    lists = {"site": [], "water": []}
    for kind, src, path, top_name, dcd_name, n, spacing, start, ff, anchor in plan:
        from mdtraj.formats import DCDTrajectoryFile
        dcd = path / dcd_name
        topology = load_topology(path / top_name, bonds.get(src) if kind == "site" else None)
        ion = find_ion(topology)
        if ion.residue.name != src:
            raise SystemExit(f"{path} carries {ion.residue.name}, not {src}")
        with DCDTrajectoryFile(str(dcd)) as fh:
            n_total = len(fh)
        dt = frame_spacing_ps(path, kind, n_total,
                              args.frame_ps if kind == "site" else args.water_frame_ps)
        mdtop, anchors = imaging(topology, anchor)
        chosen = pick_frames(n_total, dt, n, spacing, start)
        print(f"[{kind} {src}] {n_total} frames, {dt:g} ps apart: frames {chosen} "
              f"({', '.join(f'{(i + 1) * dt / 1000:.2f}' for i in chosen)} ns)", flush=True)
        for index in chosen:
            xyz, box_nm = read_frame(dcd, index, mdtop, anchors)
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
            entry = {"kind": kind, "source": src, "source_dir": str(path.resolve()),
                     "frame": index, "time_ns": round((index + 1) * dt / 1000, 3),
                     "qm_residues": residues, "coordination": coord, "notes": notes,
                     "residue_charge": res_charge, "qm_elements": sorted(elems), "jobs": {}}
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
    return {"sp": s.get("sp_energy"), "opt": q.get("opt_energy"), "converged": q.get("converged"),
            "preprocess_agrees": same}


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
                      "no thermal terms.", "frames": {}, "levels": {}}
    per = {}       # (level, kind, source, pair) -> [dE]
    for key, entry in sorted(meta.items()):
        e = {ion: job_energies(Path(p)) for ion, p in entry["jobs"].items()}
        row = {"time_ns": entry["time_ns"], "qm_residues": entry["qm_residues"],
               "converged": {i: (v or {}).get("converged") for i, v in e.items()},
               "preprocess_agrees": {i: (v or {}).get("preprocess_agrees") for i, v in e.items()},
               "dE": {}}
        for level in ("sp", "opt"):
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
    for level in ("sp", "opt"):
        lv = {"what": level_names[level], "site": {}, "water": {}, "ddE": {}}
        for m1, m2 in EXCHANGES:
            pair = f"{m1}->{m2}"
            pooled = {}
            for kind in ("site", "water"):
                rows = {src: stats(v) for (l, k, src, p), v in per.items()
                        if l == level and k == kind and p == pair}
                allv = [d for (l, k, src, p), v in per.items()
                        if l == level and k == kind and p == pair for d in v]
                pooled[kind] = stats(allv)
                lv[kind][pair] = {"by_source": rows, "all": pooled[kind]}
            dd = {}
            if pooled["site"] and pooled["water"]:
                dd["all"] = ddE(pooled["site"], pooled["water"])
                for src, st in lv["site"][pair]["by_source"].items():
                    wst = lv["water"][pair]["by_source"].get(src)
                    if st and wst:
                        dd[f"{src} frames"] = ddE(st, wst)
            lv["ddE"][pair] = dd
        report["levels"][level] = lv
    (out / "report.json").write_text(json.dumps(report, indent=1))
    print_report(report)
    print(f"\n{out / 'report.json'}")


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
            print(line)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="action", required=True)
    f = sub.add_parser("frames")
    f.add_argument("--out", required=True)
    f.add_argument("--site", action="append", help="ION=umbrella window dir (start.pdb, window.dcd)")
    f.add_argument("--water", action="append", help="ION=bulk run dir (system.pdb, --water-dcd)")
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
    s = sub.add_parser("structure")
    s.add_argument("--out", required=True)
    s.add_argument("--cutoff", type=float, default=3.1, help="Å, ion-O shell")
    c = sub.add_parser("clean")
    c.add_argument("--out", required=True)
    c.add_argument("--delete", action="store_true", help="delete (default: dry run)")
    args = p.parse_args()
    if args.action == "frames":
        args.donors = [d.upper() for d in args.donors]
        frames(args)
    elif args.action == "structure":
        structure(args)
    elif args.action == "clean":
        clean(args)
    else:
        analyse(args)


if __name__ == "__main__":
    main()
