#!/usr/bin/env python3
"""
Are the ring Glu that do not coordinate the ion charged? (decision log C22, next step 1)

The cluster exchange energies flip from Ca-selective to Ln-selective when the uncoordinated ring
Glu (pocket -1/-2 -> -4) are included, so the answer turns on whether those Glu are charged and
unscreened in the protein. For each MD trajectory this reports, per ring Glu:

  - PROPKA 3 pKa on frames spread through the run (protein only: PROPKA ignores the metal ion,
    the membrane and explicit water, so a bound frame gives the pKa the conformation alone
    implies; the ion's own field lowers it further, see coulomb_shift)
  - the environment of its carboxylate O atoms in every --stride-th frame: distance to the ion,
    Arg/Lys/His N within 4 A (salt bridges), water O within 3.5 A, Na+ within 4 A, other Asp/Glu
    carboxylate O within 3 A (a shared proton is possible), lipid P/N within 6 A
  - coulomb_shift: a crude pKa shift from the ion's charge at the carboxyl C,
    -332 q / (eps r) / 1.364, at eps 20 and 80

  python ring_pka.py --top LA/umbrella_res_035/start.pdb --traj LA/umbrella_res_035/window.dcd \
      --label la_def --out pka

Needs mdtraj (ASH env) and propka (pip; `propka3` on PATH).
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path

import numpy as np

ION_NAMES = {"CA": 2, "LA": 3, "LU": 3}
AMBER_TO_PDB = {"HIE": "HIS", "HID": "HIS", "HIP": "HIS", "CYX": "CYS", "CYM": "CYS",
                "ASH": "ASP", "GLH": "GLU", "LYN": "LYS"}
BASIC_N = ("(resname ARG and name NE NH1 NH2) or (resname LYS and name NZ) or "
           "(resname HIS HIE HID HIP and name ND1 NE2)")
ACID_O = "(resname ASP ASH and name OD1 OD2) or (resname GLU GLH and name OE1 OE2)"


def ring_residues(top, ring, offset):
    out = {}
    for num in ring:
        hits = [r for r in top.residues if r.is_protein and r.resSeq == num - offset]
        if len(hits) != 1:
            raise SystemExit(f"residue {num} (mcube {num - offset}): found {len(hits)}")
        out[num] = hits[0]
    return out


def find_ion(top):
    ions = [r for r in top.residues if r.name in ION_NAMES and r.n_atoms == 1]
    return (ions[0].atom(0).index, ION_NAMES[ions[0].name]) if ions else (None, 0)


def environment(traj, rings, ion):
    import mdtraj as md
    top = traj.topology
    sel = {"basic": top.select(BASIC_N), "water": top.select("water and element O"),
           "na": top.select("resname NA Na and element Na"), "acid": top.select(ACID_O),
           "lipid": top.select("not protein and not water and (element P or element N)")}
    cut = {"basic": 0.40, "water": 0.35, "na": 0.40, "acid": 0.30, "lipid": 0.60}
    rows = {}
    for num, res in rings.items():
        oe = [a.index for a in res.atoms if a.name in ("OE1", "OE2")]
        cd = [a.index for a in res.atoms if a.name == "CD"][0]
        rec = {}
        for key, hay in sel.items():
            hay = np.setdiff1d(hay, [a.index for a in res.atoms])
            if len(hay) == 0:
                rec[key] = [0] * traj.n_frames
                continue
            nb = md.compute_neighbors(traj, cut[key], np.array(oe), haystack_indices=hay, periodic=True)
            # distinct residues (water molecules, salt-bridge partners)
            rec[key] = [len({top.atom(int(i)).residue.index for i in f}) for f in nb]
        if ion is not None:
            pairs = np.array([[o, ion] for o in oe] + [[cd, ion]])
            d = md.compute_distances(traj, pairs, periodic=True) * 10.0
            rec["d_ion_O"] = d[:, :2].min(axis=1).tolist()
            rec["d_ion_CD"] = d[:, 2].tolist()
        rows[num] = rec
    return rows


def run_propka(traj, frames, workdir):
    """PROPKA on the protein of each frame; returns {resSeq: [pKa per frame]} for GLU/ASP."""
    workdir.mkdir(parents=True, exist_ok=True)
    prot = traj.topology.select("protein")
    out = {}
    for i in frames:
        pdb = workdir / f"frame_{i:05d}.pdb"
        traj[i].atom_slice(prot).save_pdb(str(pdb))
        text = pdb.read_text()
        for amber, std in AMBER_TO_PDB.items():
            text = re.sub(rf"(?m)^(ATOM  .{{11}}){amber}", rf"\g<1>{std}", text)
        pdb.write_text(text)
        subprocess.run(["propka3", "--quiet", pdb.name], cwd=workdir, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        summary = (workdir / pdb.with_suffix(".pka").name).read_text().split("SUMMARY OF THIS PREDICTION")[1]
        for line in summary.splitlines():
            m = re.match(r"\s+(GLU|ASP)\s+(\d+)\s+\S\s+(-?[\d.]+)\s+(-?[\d.]+)", line)
            if m:
                out.setdefault(int(m.group(2)), []).append(float(m.group(3)))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--top", required=True)
    p.add_argument("--traj", help="DCD; omit to use --top alone")
    p.add_argument("--label", required=True)
    p.add_argument("--out", default="pka")
    p.add_argument("--ring", default="156:220:240:300", help="channel numbering")
    p.add_argument("--offset", type=int, default=130, help="channel number = mcube number + offset")
    p.add_argument("--stride", type=int, default=10, help="every n-th DCD frame for the environment")
    p.add_argument("--propka-frames", type=int, default=20)
    args = p.parse_args()

    import mdtraj as md
    traj = md.load(args.traj, top=args.top, stride=args.stride) if args.traj else md.load(args.top)
    top = traj.topology
    ring = [int(r) for r in args.ring.split(":")]
    rings = ring_residues(top, ring, args.offset)
    ion, q = find_ion(top)
    env = environment(traj, rings, ion)
    frames = np.unique(np.linspace(0, traj.n_frames - 1, min(args.propka_frames, traj.n_frames)).astype(int))
    pka = run_propka(traj, frames, Path(args.out) / args.label / "propka")

    report = {"label": args.label, "top": args.top, "traj": args.traj, "frames": traj.n_frames,
              "stride": args.stride, "ion": top.atom(ion).residue.name if ion is not None else None,
              "propka_frames": frames.tolist(), "ring": {}}
    for num, res in rings.items():
        e = env[num]
        vals = pka.get(res.resSeq, [])
        row = {"resname": res.name,
               "pka_mean": round(float(np.mean(vals)), 2) if vals else None,
               "pka_sd": round(float(np.std(vals)), 2) if vals else None,
               "pka_range": [min(vals), max(vals)] if vals else None,
               "frac_frames_pka_above_7": round(float(np.mean(np.array(vals) > 7.0)), 2) if vals else None}
        for key in ("basic", "water", "na", "acid", "lipid"):
            row[f"{key}_mean"] = round(float(np.mean(e[key])), 2)
            row[f"{key}_frac_any"] = round(float(np.mean(np.array(e[key]) > 0)), 2)
        if ion is not None:
            row["d_ion_O_mean"] = round(float(np.mean(e["d_ion_O"])), 2)
            row["coordinating_frac"] = round(float(np.mean(np.array(e["d_ion_O"]) < 3.0)), 2)
            r = float(np.mean(e["d_ion_CD"]))
            row["coulomb_shift"] = {f"eps{eps}": round(-332.0 * q / (eps * r) / 1.364, 1) for eps in (20, 80)}
        report["ring"][f"GLU{num}"] = row
    out = Path(args.out) / args.label
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
