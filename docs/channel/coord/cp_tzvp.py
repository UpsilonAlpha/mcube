#!/usr/bin/env python3
"""
Larger-basis, counterpoise-corrected single points on the carboxylate-exchange geometries.

carboxylate_exchange.py optimises each species at r2SCAN/def2-SVP (lanthanides on
lcecp-1-SVP) and writes geom/<species>.xyz. Here every species gets an r2SCAN single point
in a larger basis (def2-TZVP, lanthanides on lcecp-1-TZVP), and every complex a
Boys–Bernardi counterpoise correction over all of its fragments (the ion, each water, the
acetate):

  E_cp(complex) = E(complex) + sum_f [ E_f(own basis) - E_f(all basis functions, rest as ghosts) ]

each term at the complex geometry. Reaction energies (one or two waters swapped for acetate)
are then given uncorrected and corrected, with the La−Ca, Lu−Ca, Lu−La differences, the
numbers decision log C17 compares with 12-6-4 (the MM values are read from report.json).

  python cp_tzvp.py --run coord_check --basis-file def2-tzvp_lcecp1.nw --out coord_cp
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np

HARTREE_KCAL = 627.509474
ION_CHARGE = {"Ca": 2, "La": 3, "Lu": 3}


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    rows = [l.split() for l in lines[2:2 + n]]
    return [r[0] for r in rows], np.array([[float(v) for v in r[1:4]] for r in rows])


def fragments(elems):
    """Atom index groups: the ion, each water (O H H), the acetate (C H H H C O O)."""
    out, i = [], 0
    while i < len(elems):
        if elems[i] in ION_CHARGE:
            out.append(("ion", [i])); i += 1
        elif elems[i] == "O" and elems[i + 1:i + 3] == ["H", "H"]:
            out.append(("water", [i, i + 1, i + 2])); i += 3
        else:
            out.append(("acetate", list(range(i, i + 7)))); i += 7
    return out


def frag_charge(kind, elems, idx):
    if kind == "ion":
        return ION_CHARGE[elems[idx[0]]]
    return -1 if kind == "acetate" else 0


def energy(elems, xyz, charge, ghosts, basis_file, cores, functional):
    """r2SCAN energy (kcal/mol) of the real atoms, with basis functions on `ghosts` too."""
    from pyscf import dft, gto, lib
    lib.num_threads(cores)
    atoms, basis, ecp = [], {}, {}
    for k, (e, x) in enumerate(zip(elems, xyz)):
        label = f"GHOST-{e}" if k in ghosts else e
        atoms.append((label, tuple(x)))
        basis[label] = gto.basis.load(basis_file, e)
        if k not in ghosts and label not in ecp:
            try:
                entry = gto.basis.load_ecp(basis_file, e)
            except Exception:
                entry = None
            if entry:
                ecp[label] = entry
    mol = gto.M(atom=atoms, charge=charge, spin=0, basis=basis, ecp=ecp, unit="Angstrom",
                verbose=0, max_memory=int(os.environ.get("PYSCF_MAX_MEMORY", 60000)))

    def scf():
        mf = dft.RKS(mol)
        mf.xc = functional
        mf.grids.level = 3
        mf.conv_tol = 1e-9
        mf.max_cycle = 300
        mf.damp = 0.3
        mf.level_shift = 0.2
        return mf.density_fit()

    mf = scf()
    e = mf.kernel()
    if not mf.converged:
        mf = scf().newton()
        e = mf.kernel()
    return float(e) * HARTREE_KCAL, bool(mf.converged)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", required=True, help="carboxylate_exchange.py output dir (report.json, geom/)")
    p.add_argument("--basis-file", required=True)
    p.add_argument("--functional", default="r2scan")
    p.add_argument("--cores", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    p.add_argument("--out", required=True)
    args = p.parse_args()
    run = Path(args.run).resolve()
    basis_file = str(Path(args.basis_file).resolve())
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    src = json.loads((run / "report.json").read_text())
    report_path = out / "report.json"
    report = json.loads(report_path.read_text()) if report_path.exists() else {"species": {}}
    report["settings"] = {**vars(args), "geometries_from": str(run)}

    def save():
        report_path.write_text(json.dumps(report, indent=1))

    for name in src["species"]:
        if name in report["species"] and report["species"][name].get("converged"):
            continue                                   # resume
        geo = run / "geom" / f"{name}.xyz"
        if not geo.exists():
            report["species"][name] = {"error": f"no geometry {geo}"}; save(); continue
        elems, xyz = read_xyz(geo)
        frags = fragments(elems)
        charge = sum(frag_charge(k, elems, idx) for k, idx in frags)
        t0 = time.time()
        e, conv = energy(elems, xyz, charge, set(), basis_file, args.cores, args.functional)
        row = {"E": round(e, 3), "converged": conv, "cp": 0.0}
        if len(frags) > 1:
            for kind, idx in frags:
                q = frag_charge(kind, elems, idx)
                own_e = [elems[i] for i in idx]
                e_own, c1 = energy(own_e, xyz[idx], q, set(), basis_file, args.cores, args.functional)
                rest = set(range(len(elems))) - set(idx)
                e_ghost, c2 = energy(elems, xyz, q, rest, basis_file, args.cores, args.functional)
                row["cp"] += e_own - e_ghost
                row["converged"] = row["converged"] and c1 and c2
            row["cp"] = round(row["cp"], 3)
        row["E_cp"] = round(row["E"] + row["cp"], 3)
        row["seconds"] = round(time.time() - t0, 1)
        report["species"][name] = row
        print(f"[{name}]", json.dumps(row), flush=True)
        save()

    sp, mm = report["species"], src["species"]
    rx = {}
    for ion in ("CA", "LA", "LU"):
        for k, nw in (("w-1", 1), ("w-2", 2)):
            names = (f"{ion}_{k}", f"{ion}_aqua", "water", "acetate")
            if not all(n in sp and sp[n].get("converged") for n in names):
                continue
            row = {}
            for key in ("E", "E_cp"):
                row[f"dE_{key}"] = round(sp[f"{ion}_{k}"][key] + nw * sp["water"][key]
                                         - sp[f"{ion}_aqua"][key] - sp["acetate"][key], 2)
            row["dE_mm"] = (src.get("reactions", {}).get(f"{ion}_{k}") or {}).get("dE_mm")
            rx[f"{ion}_{k}"] = row
    rel = {}
    for k in ("w-1", "w-2"):
        for a, b in (("LA", "CA"), ("LU", "CA"), ("LU", "LA")):
            if f"{a}_{k}" in rx and f"{b}_{k}" in rx:
                rel[f"{k}: {a}-{b}"] = {m: round(rx[f"{a}_{k}"][m] - rx[f"{b}_{k}"][m], 2)
                                        for m in ("dE_E", "dE_E_cp", "dE_mm")
                                        if rx[f"{a}_{k}"].get(m) is not None and rx[f"{b}_{k}"].get(m) is not None}
    report["reactions"], report["relative"] = rx, rel
    save()
    print(json.dumps({"reactions": rx, "relative": rel}, indent=1))


if __name__ == "__main__":
    main()
