"""
The exchange stage: QM/MM "vertical swap" ion-exchange energies, analysed from a config.

A vertical swap scores one MM frame with each ion in turn (same environment, only the ion's
identity changed), at a site and in a water reference, so for M1 -> M2

    ddE = < E(M2@site) - E(M1@site) >_site frames  -  < E(M2@water) - E(M1@water) >_water frames

Negative ddE = the site prefers M2. It is an energy, not a free energy: no thermal terms, one
local minimum per frame. Each (frame, ion) job is an ordinary mcube run (preprocess -> qmmm
optimisation only -> sp) in its own directory; the run directory's `frames.json` lists them.
The frames and jobs are still written by docs/channel/qmmm/exchange.py `frames`.

This module is the analysis side, run as a stage so that every choice that changes a number
is in `exchange.inputs` and every number is saved in `exchange.results`:

    mcube exchange config.yaml analyse   [--section S]   # ddE per electrostatics treatment
    mcube exchange config.yaml potential [--section S]   # MM potential at the ion (potential.json)
    mcube exchange config.yaml structure [--section S]   # the optimisations' ion-O shells

README "exchange.inputs" lists the keys.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, Optional

import numpy as np
import yaml

HARTREE_KCAL = 627.509474
IONS = {"CA": ("Ca", 2), "LA": ("La", 3), "LU": ("Lu", 3)}       # element, charge
EXCHANGES = (("CA", "LA"), ("CA", "LU"), ("LA", "LU"))
WATER = {"HOH", "WAT", "TIP3"}
LIPIDS = {"POP", "DLP", "DMP", "DOP", "DPP", "PC", "PE", "PA", "OL"}
AMINO_ACIDS = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
               "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HID", "HIE", "HIP", "HSD",
               "HSE", "HSP", "CYX", "CYM", "ASH", "GLH", "LYN"}            # = mcube's AMINO_ACIDS
def find_ion(topology):
    """The one single-atom CA/LA/LU residue of a frame."""
    ions = [a for a in topology.atoms() if a.residue.name in IONS
            and len(list(a.residue.atoms())) == 1 and a.element is not None
            and a.element.symbol == IONS[a.residue.name][0]]
    if len(ions) != 1:
        raise SystemExit(f"expected one Ca/La/Lu ion residue, found {len(ions)}")
    return ions[0]


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
                # "all" weights each source by its frame count, and sources sampled with one ion
                # can outnumber the other's water frames (the La -> Lu value of C22 was moved by
                # that); "balanced" gives each site source, referenced to its own ion's water, one
                # vote. Its error is the spread between sources, since a source's frames are
                # correlated.
                per_source = [v["ddE"] for k, v in dd.items() if k.endswith(" frames")]
                if per_source:
                    sd = float(np.std(per_source, ddof=1)) if len(per_source) > 1 else None
                    dd["balanced"] = {"ddE": round(float(np.mean(per_source)), 2),
                                      "err": round(sd, 2) if sd is not None else None,
                                      "n_sources": len(per_source)}
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
    return report, (out / name).as_posix()


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
    return res, (out / "structure.json").as_posix()

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
            if args.components:
                # which part of the box makes PME differ from the direct sum (membrane dipole or not)
                masks = {str(c): comp == c for c in sorted(set(comp))}
                parts = pme_potential(mod, cfg["qmmm"]["inputs"]["forcefields"], r_ion, masks=masks)
                row["phi_pme_V"] = round(parts["all"], 4)
                row["phi_pme_by_component_V"] = {c: round(parts[c], 4) for c in masks}
                row["phi_direct_by_component_V"] = {c: phi(sel) for c, sel in masks.items()}
            else:
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
            if "phi_pme_by_component_V" in rs[0]:
                summary[f"{kind} {src}"]["pme_minus_direct_by_component_V"] = {
                    c: {"mean": round(float(np.mean(dc)), 4), "sd": round(float(np.std(dc)), 4)}
                    for c in rs[0]["phi_pme_by_component_V"]
                    for dc in [[x["phi_pme_by_component_V"][c] - x["phi_direct_by_component_V"][c] for x in rs]]}
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
    if any("pme_minus_direct_by_component_V" in s for s in summary.values()):
        print("\nPME - direct by component, V (mean over frames; sd)")
        for name, s in summary.items():
            print(f"{name:16s}" + "".join(f"  {c} {v['mean']:+.3f} ({v['sd']:.3f})"
                                          for c, v in s["pme_minus_direct_by_component_V"].items()))
    print(f"\n{out / 'potential.json'}")
    return res, (out / "potential.json").as_posix()


def pme_potential(mod, forcefields, r_ion_A, delta=1.0, masks=None):
    """Potential (V) at r_ion from the MM charges under PME, the periodic electrostatics the MD
    ran with: a probe particle (no LJ) at the ion's position, E(+delta) - E(-delta) over 2 delta.
    The energy is exactly quadratic in the probe charge, so the central difference is exact at
    any delta; 1 e keeps it well above the single-precision noise of a 124k-atom energy. The
    probe's self and image terms are even in delta and cancel.

    `masks` ({name: bool array over the MM atoms}): also the potential of each subset alone, the
    other charges zeroed, returned as {"all": phi, name: phi, ...}. Only probe-particle terms are
    odd in delta, so exceptions and the subset's own energy cancel, and the parts sum to "all"."""
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

    def probe_phi():
        energies = []
        # The total energy: only the NonbondedForce sees the probe (no bonds, no other forces on it).
        # The first version read the NonbondedForce's force group alone, and getState(groups=...)
        # returned 0 on this OpenMM, so every frame read 0 V
        for q in (delta, -delta):
            nb.setParticleParameters(probe, q, 0.1, 0.0)
            nb.updateParametersInContext(ctx)
            energies.append(ctx.getState(getEnergy=True).getPotentialEnergy()
                            .value_in_unit(unit.kilojoule_per_mole))
        return (energies[0] - energies[1]) / (2 * delta) / 96.4853   # kJ/mol/e -> V

    phi = probe_phi()
    if phi == 0.0:
        raise RuntimeError("PME probe saw no potential: check the box vectors of the frame")
    if masks is None:
        return phi
    params = [nb.getParticleParameters(i) for i in range(probe)]
    out = {"all": phi}
    for name, sel in masks.items():
        if len(sel) != probe:
            raise ValueError(f"mask {name} covers {len(sel)} atoms, the system has {probe}")
        for i, (q, sigma, eps) in enumerate(params):
            nb.setParticleParameters(i, q if sel[i] else 0.0 * q, sigma, eps)
        out[name] = probe_phi()
    return out


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




# ==============================
# The stage
# ==============================

# Every choice that changes a reported number, with the value the analysis used before it was a
# stage (docs/channel/qmmm/exchange.py's flags), except `pooling`: the headline is the
# source-balanced ddE, not the frame-weighted "all" (decision log C22, review corrections).
DEFAULTS: Dict[str, Any] = {
    "run_dir": ".",
    "level": "sp",
    "analyse": {"electrostatics": ["direct"], "embed_radii": [], "lra": {},
                "include_unconverged": False, "pooling": "balanced"},
    "potential": {"radii": [8, 10, 15, 20, 30], "pme": True, "components": False},
    "structure": {"cutoff": 3.1},
}
ELECTROSTATICS = ("direct", "periodic")
POOLING = ("balanced", "all")


def exchange_inputs(config: str, section: Optional[str] = None) -> Dict[str, Any]:
    """DEFAULTS < exchange.inputs < `section`'s inputs; `analyse`, `potential` and `structure`
    merge key by key, so a section can change one setting of one step."""
    with open(config) as fh:
        data = yaml.safe_load(fh) or {}
    layers = [(data.get("exchange") or {}).get("inputs") or {}]
    if section and section != "exchange":
        layers.append((data.get(section) or {}).get("inputs") or {})
    inputs = {k: (dict(v) if isinstance(v, dict) else v) for k, v in DEFAULTS.items()}
    for layer in layers:
        for key, value in layer.items():
            if isinstance(inputs.get(key), dict) and isinstance(value, dict):
                inputs[key].update(value)
            else:
                inputs[key] = value
    unknown = sorted(set(inputs) - set(DEFAULTS))
    if unknown:
        raise ValueError(f"exchange.inputs: unknown keys {unknown} (known: {sorted(DEFAULTS)})")
    for step in ("analyse", "potential", "structure"):
        bad = sorted(set(inputs[step]) - set(DEFAULTS[step]))
        if bad:
            raise ValueError(f"exchange.inputs.{step}: unknown keys {bad} (known: {sorted(DEFAULTS[step])})")
    a = inputs["analyse"]
    if isinstance(a["electrostatics"], str):
        a["electrostatics"] = [a["electrostatics"]]
    bad = [t for t in a["electrostatics"] if t not in ELECTROSTATICS]
    if bad:
        raise ValueError(f"exchange.inputs.analyse.electrostatics: {bad} not in {ELECTROSTATICS} "
                         "(finite embeddings go in embed_radii)")
    if a["pooling"] not in POOLING:
        raise ValueError(f"exchange.inputs.analyse.pooling must be one of {POOLING}")
    for name, members in (a["lra"] or {}).items():
        if not isinstance(members, (list, tuple)) or len(members) != 2:
            raise ValueError(f"exchange.inputs.analyse.lra.{name}: two site sources, e.g. [CAla, LA]")
    return inputs


def run_directory(config: str, inputs) -> Path:
    """`run_dir`, read against the config's directory (README "Where paths are read from")."""
    path = Path(str(inputs["run_dir"]))
    if not path.is_absolute():
        path = Path(config).resolve().parent / path
    if not (path / "frames.json").exists():
        raise FileNotFoundError(f"exchange.inputs.run_dir: no frames.json in {path} "
                                "(docs/channel/qmmm/exchange.py frames writes it)")
    return path.resolve()


def headline(report, level: str, pooling: str) -> Dict[str, Any]:
    """The numbers a reader quotes, per exchange: the pooled ddE (`pooling`), each source's,
    and the linear-response pairs."""
    lv = report["levels"].get(level)
    if lv is None:
        return {}
    out = {}
    for pair, dd in lv["ddE"].items():
        if not dd:
            continue
        out[pair] = {"ddE": dd.get(pooling), "pooling": pooling,
                     "by_source": {k[:-len(" frames")]: v for k, v in dd.items() if k.endswith(" frames")}}
        lras = {name.split(" ", 1)[1]: v for name, v in lv["lra"].items() if name.startswith(pair + " ")}
        if lras:
            out[pair]["lra"] = lras
    return out


def run_analyse(config, section, inputs, out: Path) -> Dict[str, Any]:
    from .umbrella import save_nested
    a, level = inputs["analyse"], inputs["level"]
    lra = [f"{name}={'+'.join(m)}" for name, m in (a["lra"] or {}).items()]
    treatments = [(t, None) for t in a["electrostatics"]] + [(f"R{float(r):g}", float(r))
                                                            for r in a["embed_radii"] or []]
    result = {"level": level, "pooling": a["pooling"], "inputs_used": a, "treatments": {}}
    for name, radius in treatments:
        print(f"\n######## {name}")
        args = SimpleNamespace(out=out, include_unconverged=bool(a["include_unconverged"]), lra=lra,
                               embed_radius=radius, periodic=(name == "periodic"))
        report, path = analyse(args)
        if level not in report["levels"]:
            raise ValueError(f"exchange.inputs.level {level!r}: the jobs have {sorted(report['levels'])}")
        result["level_name"] = report["levels"][level]["what"]
        result["treatments"][name] = {"report": path, "ddE": headline(report, level, a["pooling"])}
    # replace, not merge: a treatment dropped from the inputs must not linger in the results
    save_nested(config, section, ["analyse"], None)
    save_nested(config, section, ["analyse"], result)
    return result


def run_potential(config, section, inputs, out: Path) -> Dict[str, Any]:
    from .umbrella import save_nested
    p = inputs["potential"]
    if p["components"] and not p["pme"]:
        raise ValueError("exchange.inputs.potential.components needs pme: true")
    res, path = potential(SimpleNamespace(out=out, radii=[f"{float(r):g}" for r in p["radii"]],
                                          pme=bool(p["pme"]), components=bool(p["components"])))
    result = {"file": path, "inputs_used": p, "summary": res["summary"]}
    save_nested(config, section, ["potential"], None)
    save_nested(config, section, ["potential"], result)
    return result


def run_structure(config, section, inputs, out: Path) -> Dict[str, Any]:
    from .umbrella import save_nested
    s = inputs["structure"]
    res, path = structure(SimpleNamespace(out=out, cutoff=float(s["cutoff"])))
    result = {"file": path, "inputs_used": s, "summary": res["summary"]}
    save_nested(config, section, ["structure"], None)
    save_nested(config, section, ["structure"], result)
    return result


def run(args) -> Dict[str, Any]:
    config = str(Path(args.config).resolve())
    section = getattr(args, "section", None) or "exchange"
    inputs = exchange_inputs(config, section)
    out = run_directory(config, inputs)
    steps = {"analyse": run_analyse, "potential": run_potential, "structure": run_structure}
    if args.action not in steps:
        raise SystemExit(f"unknown exchange action {args.action!r}")
    return steps[args.action](config, section, inputs, out)
