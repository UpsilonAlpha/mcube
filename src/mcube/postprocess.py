import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.analysis.rms import RMSF
from MDAnalysis.analysis import align
from MDAnalysis.coordinates.DCD import DCDWriter


import numpy as np
import networkx as nx
from collections import defaultdict
import polars as pl
import matplotlib.pyplot as plt
import os
from pathlib import Path
from typing import Any, Dict
from ash import *
import mdtraj

from .preprocess import (
    save_result,
    load_section,
    ligand_resname,
    config_relative,
    coordination_rows,
    coordination_from_rows,
    resolve_atoms,
    separation,
)
from .md import md_membrane
from .membrane_analysis import report as membrane_report, options_from_config

# =========================
# Helper functions
# =========================

def residue_id(atom):
    """Return a unique residue identifier tuple"""
    r = atom.residue
    return (r.segid, r.resid, r.resname)

def centroid(atoms):
    """Compute centroid of a group of atoms"""
    return atoms.positions.mean(axis=0)

def ring_normal(atoms):
    """Compute approximate normal vector of aromatic ring"""
    coords = atoms.positions
    v1 = coords[1] - coords[0]
    v2 = coords[2] - coords[0]
    n = np.cross(v1, v2)
    return n / np.linalg.norm(n)


def resolve_overlaps(pos, min_dist=0.08, max_iter=100):
    nodes = list(pos.keys())
    for _ in range(max_iter):
        moved = False
        for i in range(len(nodes)):
            for j in range(i+1, len(nodes)):
                n1, n2 = nodes[i], nodes[j]
                d = pos[n1] - pos[n2]
                dist = np.linalg.norm(d)

                if dist < min_dist:
                    direction = d / (dist + 1e-9)
                    shift = (min_dist - dist) * 0.5 * direction
                    pos[n1] += shift
                    pos[n2] -= shift
                    moved = True
        if not moved:
            break
    return pos
# =========================
# Main pipeline
# =========================

def coordination_occupancy(pdb_file, trajectory, inputs, radius=None,
                           per_residue=False, margin=1.0, stride=1):
    """Per-ligand occupancy and the coordination-number distribution.

    `md` reports the sphere before and after a run, which is two frames, and a
    single frame cannot tell a ligand that is bound from one that happens to be
    close at that instant. The coordination *number* is also the quantity
    semiempirical methods are known to get wrong (D114: AM1/d out by 7.55
    kcal/mol on a 4 -> 5 change), so it is reported as a distribution and not
    as a mean - a site flipping between 4 and 5 and a site sitting at 4.5 are
    different things.

    The rule comes from `coordination_from_rows`, so the cutoff, the donor
    elements and the minimum image are the ones `preprocess` used to build the
    restraints rather than a second set chosen here. `per_residue` defaults to
    False, unlike the restraint path: a carboxylate binding with both oxygens
    is a different sphere from one binding with one, and collapsing to the
    residue hides exactly that (D166 measured it in two runs).

    Donors just past the cutoff are counted too, because a ligand at 2.7 A with
    a 2.6 cutoff is a borderline call about the site and not an absence.
    """
    import numpy as np

    rows, cell = coordination_rows(pdb_file)
    # MDAnalysis notices a topology that does not match the trajectory, but it
    # reports it from inside its reader, which reads as a corrupt file rather
    # than as the wrong pair of files. Said here instead: zipping the wrong
    # topology onto a trajectory would report one atom's distances under
    # another atom's name, which no later check would catch.
    try:
        u = mda.Universe(pdb_file, trajectory)
    except Exception as exc:
        raise ValueError(
            f"{trajectory} cannot be read against {pdb_file} ({len(rows)} "
            f"atoms): {exc}. The topology must be the structure the trajectory "
            "was produced from, so its atoms and the frames agree.") from exc
    if len(u.atoms) != len(rows):
        raise ValueError(
            f"{pdb_file} parses as {len(rows)} atoms here but MDAnalysis reads "
            f"{len(u.atoms)}. The two disagree about the same file, so the "
            "rows and the coordinates cannot be matched up.")

    cutoff = float(radius if radius is not None else (inputs.get("cutoff") or 2.8))
    seen, counts, frames = {}, [], 0
    for _ in u.trajectory[::stride]:
        xyz = u.atoms.positions
        frame_rows = [r[:5] + (tuple(map(float, xyz[i])),)
                      for i, r in enumerate(rows)]
        shell = coordination_from_rows(frame_rows, cell, inputs,
                                       cutoff + margin, per_residue)
        frames += 1
        inside = 0
        for row in shell:
            key = (row["metal"], row["donor"], row["resname"])
            rec = seen.setdefault(key, {"n": 0, "d": []})
            rec["d"].append(row["distance"])
            if row["distance"] <= cutoff:
                rec["n"] += 1
                inside += 1
        counts.append(inside)

    out = []
    for (metal, donor, resname), rec in seen.items():
        d = np.array(rec["d"])
        out.append({"metal": metal, "donor": donor, "resname": resname,
                    "occupancy": round(100.0 * rec["n"] / frames, 1),
                    "min": round(float(d.min()), 2),
                    "mean": round(float(d.mean()), 2),
                    "max": round(float(d.max()), 2)})
    out.sort(key=lambda r: (-r["occupancy"], r["mean"]))
    hist = {n: counts.count(n) for n in sorted(set(counts))}
    return {"frames": frames, "cutoff": cutoff, "margin": margin,
            "per_residue": per_residue, "ligands": out,
            "coordination_number": {str(k): v for k, v in hist.items()}}


def contact_occupancy(pdb_file, trajectory, monitor, cutoff=3.5, stride=1):
    """Named atom pairs over a trajectory: occupancy and distance spread.

    `coordination_occupancy` answers what is on the metal. This answers
    everything else - which residue holds the substrate's carboxylate, whether
    a glutamine reaches the carbonyl - and it answers it for the pairs the
    config names, so a contact is reported because the config asked for it and
    not because someone wrote a one-off script for this protein.

    `monitor` is {name: [selector, selector]}, the selectors being the same
    `chain:resid:atom` that restraints and `zero_charges` take: `resolve_atoms`
    is the one parser for them, so a config means the same thing in every
    stage. Occupancy is the fraction of frames within `cutoff`, 3.5 A by
    default - the usual polar-contact convention for a heavy-atom pair.

    A pair whose atoms are not in this topology is **reported as missing, not
    skipped**: a contact absent from the system is a fact about the model, and
    a selector that silently resolves to nothing is how a monitor goes
    unwatched (D120, where the oxygen that migrated onto the metal was in no
    monitor at all).
    """
    import numpy as np

    rows, cell = coordination_rows(pdb_file)
    try:
        u = mda.Universe(pdb_file, trajectory)
    except Exception as exc:
        raise ValueError(
            f"{trajectory} cannot be read against {pdb_file}: {exc}. The "
            "topology must be the structure the trajectory was produced "
            "from.") from exc

    named, idx, missing = [], [], []
    for name, sel in (monitor or {}).items():
        if not (isinstance(sel, (list, tuple)) and len(sel) == 2):
            missing.append({"pair": name,
                            "why": "needs exactly two selectors, got %r" % (sel,)})
            continue
        try:
            i, j = resolve_atoms(list(sel), pdb_file)
        except ValueError as exc:
            missing.append({"pair": name, "why": str(exc)})
            continue
        named.append(name)
        idx.append((int(i), int(j)))

    dist = [[] for _ in named]
    frames = 0
    for _ in u.trajectory[::stride]:
        xyz = u.atoms.positions
        frames += 1
        for k, (i, j) in enumerate(idx):
            dist[k].append(separation(tuple(map(float, xyz[i])),
                                      tuple(map(float, xyz[j])), cell))

    out = []
    for k, name in enumerate(named):
        d = np.array(dist[k])
        out.append({"pair": name,
                    "occupancy": round(float((d < cutoff).mean()) * 100.0, 1),
                    "min": round(float(d.min()), 2),
                    "mean": round(float(d.mean()), 2),
                    "max": round(float(d.max()), 2)})
    out.sort(key=lambda r: (-r["occupancy"], r["mean"]))
    return {"frames": frames, "cutoff": cutoff, "pairs": out,
            "missing": missing,
            "_distances": {name: dist[k] for k, name in enumerate(named)}}


def contact_profile(rep, want, width=0.05, temperature=300.0):
    """-RT ln P(r) for every monitored pair whose name contains `want`.

    Unbiased sampling already carries the free energy along any distance it
    explored: the histogram of a 100 ns trajectory is a potential of mean force
    over the range it visited, for nothing but the cost of reading the frames.

    It is only as good as the sampling, so the count per bin is reported beside
    the energy - a bin with a handful of frames is noise wearing a number, and
    the top of a profile like this is always the thin end. This is not a
    substitute for umbrella sampling, which is what reaches the range the
    molecule does not visit on its own; it is the part already paid for.
    """
    import numpy as np

    rt = 0.0019872041 * float(temperature)        # kcal/mol
    out = []
    for name, values in (rep.get("_distances") or {}).items():
        if want.lower() not in name.lower():
            continue
        d = np.array(values)
        lo = np.floor(d.min() / width) * width
        hi = np.ceil(d.max() / width) * width
        counts, edges = np.histogram(d, bins=np.arange(lo, hi + width, width))
        centres = 0.5 * (edges[:-1] + edges[1:])
        keep = counts > 0
        # Jacobian: a distance shell grows as r^2, so the radial density is
        # divided by it before the log, or every profile tilts outward.
        p = counts[keep] / (centres[keep] ** 2)
        g = -rt * np.log(p / p.max())
        bins = [{"r": round(float(c), 3), "frames": int(n),
                 "dG": round(float(e), 2)}
                for c, n, e in zip(centres[keep], counts[keep], g)]
        out.append({"pair": name, "temperature": float(temperature),
                    "bin": width, "minimum": round(float(
                        centres[keep][int(np.argmin(g))]), 2),
                    "sampled_frames": int(counts.sum()), "bins": bins})
    return out


def report_contacts(rep, label: str) -> None:
    print("\n===== Monitored contacts over %d frames: %s ====="
          % (rep["frames"], label))
    print("  cutoff %.2f A" % rep["cutoff"])
    print("  %-26s  occ%%    min  mean   max" % "pair")
    for row in rep["pairs"]:
        print("  %-26s %5.1f  %5.2f %5.2f  %5.2f"
              % (row["pair"], row["occupancy"], row["min"], row["mean"],
                 row["max"]))
    for row in rep["missing"]:
        print("  %-26s     -  %s" % (row["pair"], row["why"]))


def report_contact_profiles(profiles) -> None:
    for prof in profiles:
        print("\n  -- %s: -RT ln P at %.0f K, %.2f A bins, %d frames sampled"
              % (prof["pair"], prof["temperature"], prof["bin"],
                 prof["sampled_frames"]))
        print("     %6s %8s %10s" % ("r", "frames", "dG kcal"))
        for b in prof["bins"]:
            print("     %6.2f %8d %10.2f" % (b["r"], b["frames"], b["dG"]))
        print("     minimum at %.2f A" % prof["minimum"])


def report_coordination_occupancy(rep, label: str) -> None:
    print("\n===== Metal coordination over %d frames: %s ====="
          % (rep["frames"], label))
    print("  cutoff %.2f A, +%.1f A margin shown, per_residue=%s"
          % (rep["cutoff"], rep["margin"], rep["per_residue"]))
    for row in rep["ligands"]:
        mark = "" if row["occupancy"] > 0 else "   (never inside the cutoff)"
        print("  %-12s %-16s %-4s %6.1f%%  %4.2f /%4.2f /%4.2f%s"
              % (row["metal"], row["donor"], row["resname"], row["occupancy"],
                 row["min"], row["mean"], row["max"], mark))
    print("  coordination number: "
          + ", ".join("%s in %d frames" % (k, v)
                      for k, v in rep["coordination_number"].items()))

def run(args):

    result: Dict[str, Any] = {}

    preprocess = load_section(args.config, "preprocess")
    lig = (preprocess.get("results") or {}).get("ligand_resname") or ligand_resname(preprocess.get("inputs") or {})

    # Which trajectory was analysed is part of the result, so it comes from
    # `postprocess.inputs` (read against the config's directory) rather than
    # from a flag; with neither key set it is the run the config already
    # records. `--trajectory` used to be accepted here and then ignored.
    section = getattr(args, "section", None) or "postprocess"
    pp_inputs = load_section(args.config, section).get("inputs") or {}
    given_pdb = config_relative(args.config, pp_inputs.get("pdb_file"))
    given_traj = config_relative(args.config, pp_inputs.get("trajectory"))
    for key, value in (("pdb_file", given_pdb), ("trajectory", given_traj)):
        if value and not Path(value).is_file():
            raise FileNotFoundError(
                f"postprocess.inputs.{key} does not exist: {value}. Paths in "
                f"`inputs` are read against the config's directory "
                f"({Path(args.config).resolve().parent}), not the directory the "
                "stage runs in.")

    if args.qmmm:
        qmmm = load_section(args.config, "qmmm")
        md = load_section(args.config, "md")

        pdbfile = given_pdb or qmmm["results"]["lastframe"]
        trajectory = given_traj or qmmm["results"]["trajfile"]
        optimization = qmmm["results"]["opttraj"]

        u = mda.Universe(pdbfile, optimization)

        # 2. Setup the DCD writer
        output_dcd = 'output.dcd'
        with DCDWriter(output_dcd, u.atoms.n_atoms) as W:
            # 3. Iterate through frames and write to DCD
            for frame in u.trajectory:
                W.write(u)

    else:
        # Which md run to analyse. A run submitted as `md md_prod100` saves its
        # results under that name, so defaulting to a section called `md` would
        # find nothing and the paths would have to be repeated by hand.
        md_section = pp_inputs.get("md_section") or "md"
        md = load_section(args.config, md_section)
        if not (md.get("results") or {}) and not (given_pdb and given_traj):
            raise KeyError(
                f"`{md_section}` has no results in {args.config}, so there is "
                "no trajectory to analyse. Name the run with "
                f"`{section}.inputs.md_section` (the section the md job was "
                "submitted as), or give `pdb_file` and `trajectory` outright.")
        pdbfile = given_pdb or md["results"]["lastframe"]
        trajectory = given_traj or md["results"]["trajfile"]
        csv = md["results"]["csv"]

        # The coordination sphere over the whole trajectory, by the same rule
        # that built the restraints. Its settings come from `preprocess.inputs`
        # and what it measured is saved into `postprocess.results`, so the
        # config records the analysis as well as the run - which `--cutoff` and
        # `--metals` on the command line never did (D167).
        coord_spec = pp_inputs.get("coordination")
        if coord_spec is not False:
            spec = coord_spec if isinstance(coord_spec, dict) else {}
            pre_inputs = preprocess.get("inputs") or {}
            if not (pre_inputs.get("metals") or []):
                print("[INFO] No `preprocess.inputs.metals`: no coordination "
                      "analysis. Set it, or `postprocess.inputs.coordination: "
                      "false` to say so deliberately.")
            else:
                rep = coordination_occupancy(
                    pdbfile, trajectory, pre_inputs,
                    radius=spec.get("cutoff"),
                    per_residue=bool(spec.get("per_residue", False)),
                    margin=float(spec.get("margin", 1.0)),
                    stride=int(spec.get("stride", 1)))
                report_coordination_occupancy(rep, Path(trajectory).name)
                result["coordination"] = rep
                result["coordination_trajectory"] = Path(trajectory).resolve().as_posix()
                save_result(result, filename=args.config, section=section)

        # The named pairs, which is everything the metal shell does not cover.
        # `monitor` sits in this section's inputs rather than being discovered,
        # so a contact is reported because the config asks for it.
        monitor = pp_inputs.get("monitor")
        if monitor:
            spec = pp_inputs.get("contacts")
            spec = spec if isinstance(spec, dict) else {}
            crep = contact_occupancy(
                pdbfile, trajectory, monitor,
                cutoff=float(spec.get("cutoff", 3.5)),
                stride=int(spec.get("stride", 1)))
            report_contacts(crep, Path(trajectory).name)
            want = spec.get("profile")
            profiles = (contact_profile(crep, want,
                                        width=float(spec.get("bin", 0.05)),
                                        temperature=float(
                                            spec.get("temperature", 300.0)))
                        if want else [])
            if profiles:
                report_contact_profiles(profiles)
            # The per-frame distances are megabytes and the config is read by
            # hand; the summary is what belongs in it.
            crep.pop("_distances", None)
            result["contacts"] = crep
            if profiles:
                result["contact_profiles"] = profiles
            result["contacts_trajectory"] = Path(trajectory).resolve().as_posix()
            save_result(result, filename=args.config, section=section)

        df = pl.read_csv(csv)

        steps = df['#"Step"']
        time = df['Time (ps)']
        temperature = df['Temperature (K)']

        #Creating data list
        data_list=[temperature]
        #Density and Volume only present for NPT
        try:
            density = df['Density (g/mL)']
            data_list.append(density)
            volume = df['Box Volume (nm^3)']
            data_list.append(volume)
        except:
            pass

        #Looping over data_list and plot
        for pd_col in data_list:
            np_array = pd_col.to_numpy()
            label=pd_col.name
            label_no_unit = label.split('(')[0].replace(' ','')
            print(label_no_unit)
            eplot = ASH_plot(label, num_subplots=1, x_axislabel="Steps", y_axislabel=label)
            eplot.addseries(0, x_list=steps.to_numpy(), y_list=np_array, label=label, color='blue', line=True, scatter=True)
            eplot.savefig(label_no_unit)

        # A membrane system gets the membrane report (box, bilayer, partner pose) in place of the
        # ligand-centred analysis below; its options are `postprocess.inputs.membrane` (optional).
        if md_membrane(md.get("inputs") or {}, preprocess):
            spec = pp_inputs.get("membrane") or {}
            out = spec.get("out", "membrane_report")
            rep = membrane_report(pdbfile, trajectory, out, **options_from_config(spec))
            result["membrane_report"] = [Path(f).resolve().as_posix() for f in rep["files"]]
            save_result(result, filename=args.config, section=section)
            return rep


    #Loading using mdtraj
    system = mdtraj.load(pdbfile)
    traj = mdtraj.load(trajectory, top=system)

    print(f"This trajectory contains {traj.n_frames} frames")

    #Calculating full RMSD (flawed) w.r.t. first frame
    rmsd_all= mdtraj.rmsd(traj, traj[0], 0)

    #Sub-system selection: Defining heavy atoms

    #Selection: All non-H atoms (also flawed because of solvent)
    heavy_atoms = [atom.index for atom in traj.topology.atoms if atom.element.symbol != 'H']

    #Selection: All non-H atoms in protein
    heavy_protein_atoms = traj.topology.select("protein and (element !=  H)")

    #RMSD for heavy protein atoms  w.r.t. first frame
    rmsd_heavy_protein = mdtraj.rmsd(traj, traj[0], 0, atom_indices=heavy_protein_atoms)

    #Plotting using ASH-plot (matplotlib)
    x_label="Frames in trajectory"
    y_label="RMSD (nm)"
    filelabel="RMSD"
    eplot = ASH_plot(filelabel, num_subplots=1, x_axislabel=x_label, y_axislabel=y_label)
    eplot.addseries(0, x_list=traj.time, y_list=rmsd_heavy_protein, label=y_label, color='blue', line=True, scatter=True)
    eplot.savefig(filelabel)



    # Load system (also needed for RMSF when redrawing)
    u = mda.Universe(pdbfile, trajectory)
    protein = u.select_atoms("protein")

    if args.redraw == False:

        # Define atom groups
        acidic = u.select_atoms("resname ASP GLU and name OD1 OD2 OE1 OE2")
        basic  = u.select_atoms("resname LYS ARG and name NZ NH1 NH2")

        aromatic = u.select_atoms("resname PHE TYR TRP HIS")
        hydrophobic = u.select_atoms("resname ALA VAL LEU ILE MET PHE TRP PRO")

        metals = u.select_atoms(f"name {' '.join(args.metals)}")
        donors = u.select_atoms("element N O S")

        # Precompute aromatic atoms per residue
        rings = []
        for res in aromatic.residues:
            ring_atoms = res.atoms.select_atoms("name CG CD1 CD2 CE1 CE2 CZ NE1")
            if len(ring_atoms) >= 5:
                rings.append(ring_atoms)

        # Storage for interaction counts
        edge_counts = defaultdict(lambda: defaultdict(int))
        n_frames = len(u.trajectory)

        # =========================
        # Hydrogen bond analysis
        # =========================

        print("Calculating hydrogen bonds")

        hbond_analysis = HydrogenBondAnalysis(
            u,
            donors_sel=f"(protein or resname {lig}) and (name N* or name O* or name S*) ",
            hydrogens_sel=f"(protein or resname {lig}) and name H*",
            acceptors_sel=f"(protein or resname {lig}) and (name O* or name N* or name S*)",
            d_a_cutoff=3.2,
            d_h_a_angle_cutoff=150,
        )
        hbond_analysis.run()

        # Convert results to Polars DataFrame
        df = pl.DataFrame(
            hbond_analysis.results.hbonds,
            schema=["frame", "donor", "hydrogen", "acceptor", "distance", "angle"]
        ).with_columns([
            pl.col("frame").cast(pl.Int64),
            pl.col("donor").cast(pl.Int64),
            pl.col("acceptor").cast(pl.Int64),
        ])

        # Count H-bond occurrences per frame
        for frame, group in df.group_by("frame", maintain_order=True):

            seen = set()

            donors = group["donor"].to_numpy()
            acceptors = group["acceptor"].to_numpy()

            for d, a in zip(donors, acceptors):
                res1 = residue_id(u.atoms[d])
                res2 = residue_id(u.atoms[a])

                # prevent self-interactions
                if res1 != res2:
                    edge = tuple(sorted([res1, res2]))
                    seen.add(edge)

            for edge in seen:
                edge_counts[edge]["hbond"] += 1            


        # =========================
        # Frame-wise interactions
        # =========================

        print("Calculating other interactions")
        frame_edges = []
        for ts in u.trajectory:

            seen_frame = {
                "salt": set(),
                "pi_pi": set(),
                "cation_pi": set(),
                "coordination": set(),
                "hydrophobic": set()
            }

            coord_residues = set()

            # ---- Coordination bonds ----
            donor_atoms = u.select_atoms("element N O S")
            for m in metals:
                distances = np.linalg.norm(donor_atoms.positions - m.position, axis=1)
                nearby = donor_atoms[distances < args.cutoff]
                for atom in nearby:
                    res = residue_id(atom)
                    coord_residues.add(res)
                    metal_res = ("METAL", m.index, m.name)
                    edge = tuple(sorted([metal_res, res]))
                    seen_frame["coordination"].add(edge)


            # ---- Salt bridges ----
            for a in acidic:
                for b in basic:
                    if residue_id(a) != residue_id(b):
                        if np.linalg.norm(a.position - b.position) < 4.0:
                            edge = tuple(sorted([residue_id(a), residue_id(b)]))
                            seen_frame["salt"].add(edge)
            
            # ---- Pi-Pi stacking ----
            for i, ring1 in enumerate(rings):
                c1 = centroid(ring1)
                n1 = ring_normal(ring1)

                for ring2 in rings[i+1:]:
                    c2 = centroid(ring2)
                    n2 = ring_normal(ring2)

                    dist = np.linalg.norm(c1 - c2)
                    angle = np.degrees(np.arccos(np.clip(np.dot(n1, n2), -1, 1)))

                    if dist < 5.0 and (angle < 30 or 60 < angle < 120):
                        edge = tuple(sorted([residue_id(ring1[0]), residue_id(ring2[0])]))
                        seen_frame["pi_pi"].add(edge)

            # ---- Cation–Pi ----
            for ring in rings:
                c = centroid(ring)
                for b in basic:
                    if residue_id(ring[0]) != residue_id(b):
                        if np.linalg.norm(c - b.position) < 6.0:
                            edge = tuple(sorted([residue_id(ring[0]), residue_id(b)]))
                            seen_frame["cation_pi"].add(edge)
            '''
            # ---- Hydrophobic contacts ----
            for i, a in enumerate(hydrophobic):
                for b in hydrophobic[i+1:]:
                    if residue_id(a) != residue_id(b):
                        if np.linalg.norm(a.position - b.position) < 4.0:
                            edge = tuple(sorted([residue_id(a), residue_id(b)]))
                            seen_frame["hydrophobic"].add(edge)
            '''
            # Accumulate counts
            for itype, edges in seen_frame.items():
                for edge in edges:
                    edge_counts[edge][itype] += 1

        # =========================
        # Build network
        # =========================

        print("Building interaction network")

        # build graph and dataframe
        G = nx.Graph()
        rows = []
        for edge, types in edge_counts.items():
            res1, res2 = edge
            for itype, count in types.items():
                occ = count / n_frames
                if occ > args.occupancy:
                    G.add_edge(
                        res1,
                        res2,
                        interaction=itype,
                        occupancy=occ,
                        color={
                            "hbond": "lightblue",
                            "salt": "blue",
                            "pi_pi": "purple",
                            "cation_pi": "red",
                            "coordination": "green",
                            "hydrophobic": "orange"
                        }[itype]
                    )

                    rows.append({
                        "res1_segid": res1[0],
                        "res1_resid": res1[1],
                        "res1_resname": res1[2],

                        "res2_segid": res2[0],
                        "res2_resid": res2[1],
                        "res2_resname": res2[2],

                        "interaction_type": itype,
                        "occupancy": occ,
                        "pair": f"{res1[2]}-{res2[2]}",
                    })

        # =========================
        # Save network plot and csv
        # =========================
        df_out = pl.DataFrame(rows)
        df_out = df_out.sort(["interaction_type", "occupancy"],descending=[False, True])
        df_out.write_csv("interaction_network.csv")

    else:
        print("Redrawing")
        G = nx.Graph()

        df = pl.read_csv("interaction_network.csv")

        for row in df.iter_rows(named=True):
            if row["occupancy"] > args.occupancy:

                res1 = (row["res1_segid"], row["res1_resid"], row["res1_resname"])
                res2 = (row["res2_segid"], row["res2_resid"], row["res2_resname"])


                G.add_edge(
                    res1,
                    res2,
                    interaction=row["interaction_type"],
                    occupancy=row["occupancy"],
                    color={
                        "hbond": "lightblue",
                        "salt": "blue",
                        "pi_pi": "purple",
                        "cation_pi": "red",
                        "coordination": "green",
                        "hydrophobic": "orange"
                    }[row["interaction_type"]]
                )


    plt.figure(figsize=(10, 10))

    labels = {
        node: f"{node[2]}{node[1]}"  # resname + resid
        for node in G.nodes
    }

    node_type = dict()
    for n in G.nodes:
        resname = n[2]

        if resname in ["ASP", "GLU"]:
            node_type[n] = "#b3c7ff"

        elif resname in ["LYS", "ARG", "HIS"]:
            node_type[n] = "#ffb3b3"

        elif resname in ["PHE", "TYR", "TRP"]:
            node_type[n] = "#ffd699"

        elif resname in ["HOH", "WAT"]:
            node_type[n] = "#b3e6ff"

        elif resname in ["ZN", "MG", "FE", "CU"]:  # adjust for your system
            node_type[n] = "#b3ffb3"

        elif resname in [lig, "LIG", "SUB"]:  # ligand
            node_type[n] = "#d9b3ff"

        else:
            node_type[n] =  "#dddddd"

    
    node_colors = [node_type[n] for n in G.nodes]

    node_sizes = [
        100 + 100 * sum(d["occupancy"] for _,_,d in G.edges(n, data=True))
        for n in G.nodes
    ]
    pos = nx.spring_layout(G, seed=42, k=args.spring_constant/ np.sqrt(len(G.nodes)), iterations=args.iterations)
    pos = resolve_overlaps(pos)

    colors = [d["color"] for _, _, d in G.edges(data=True)]
    edge_widths = [2 * d["occupancy"] for _,_,d in G.edges(data=True)]
    nx.draw(G, pos, node_color=node_colors, node_size=node_sizes, edgecolors="white", edge_color=colors, labels=labels, font_weight="bold", font_size=3, width=edge_widths)
    plt.savefig("interaction_network.png", dpi=500)
    plt.close()

    # =========================
    # RMSF calculation
    # =========================

    print("Calculating RMSF")

    # Align trajectory to first frame using Cα atoms
    align.AlignTraj(
        u,
        u,  # align to itself
        select="protein and name CA",
        in_memory=True
    ).run()

    rmsf = RMSF(protein.select_atoms("name CA")).run()

    # Map RMSF to B-factor
    for res, val in zip(protein.residues, rmsf.results.rmsf):
        for atom in res.atoms:
            atom.tempfactor = val

    # Save structure
    output_pdb = "rmsf_bfactor.pdb"
    protein.write(output_pdb)

    print(f"Saved RMSF structure to {output_pdb}")
    print("Saved interaction network to interaction_network.png")

    result["rmsf_pdb"] = Path("rmsf_bfactor.pdb").resolve().as_posix()
    result["network_graph"] = Path("interaction_network.png").resolve().as_posix()
    result["network_csv"] = Path("interaction_network.csv").resolve().as_posix()
    save_result(result, filename=args.config, section=section)


    return output_pdb



#ffmpeg -framerate 10 -pattern_type glob -i "network_frames/frame_*.png" -c:v libx264 -pix_fmt yuv420p out.mp4