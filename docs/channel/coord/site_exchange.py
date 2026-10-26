#!/usr/bin/env python3
"""
Which ion does the constriction site prefer? Ion exchange free energies at the Glu ring
(docs/channel/question-and-criteria.md), after Angelova, Dudev et al. (PCCP 2023, PMC10094714):
DFT on a cluster of the site, solvation by a continuum, and the reaction

    M2(H2O)n2 + [M1-site]  ->  [M2-site] + M1(H2O)n1 + (n2 - n1) H2O

for Ca->La, Ca->Lu and La->Lu. Negative means the site prefers M2.

The site comes from an MD frame in which an ion is bound (the residence runs, decision log
C21). Each listed residue (channel numbering) is cut at CA, which keeps its HA and gets two H
caps along the old CA-N and CA-C bonds (Glu -> butanoate, Asp -> propanoate, as Dudev). Waters
whose O lies within --shell of the ion come along, and so does any other protein residue with a
side-chain O there (reported). Residues that coordinate the ion keep their side chain free except
CA and its H atoms. Residues that do not coordinate it (e.g. Glu300) are frozen whole, so they keep
the pocket's charge and shape without collapsing onto the ion in the gas phase. The same cluster
carries each ion in turn, and each is optimised (r2SCAN/def2-SVP, La/Lu on lcecp-1-SVP).

Energies: r2SCAN single points in the larger basis (def2-TZVP, lcecp-1-TZVP), gas phase and
IEF-PCM. The site is at eps 4 (buried) and eps 29 (solvent-exposed), the dielectrics Dudev
use; aqua ions and water at eps 78.36. Dudev use SMD. PySCF's pip build lacks the SMD
non-electrostatic (CDS) library, so this is electrostatics only. The site complexes differ only in
the ion, so their CDS terms nearly cancel; the aqua ions' do not quite. Ion cavity radii are set
explicitly (--radii), because PySCF's default for La/Lu is a 2.0 A placeholder. Aqua ions and
water are the r2SCAN/SVP geometries from carboxylate_exchange.py (--aqua: Ca 8, La 9, Lu 8 waters).
No thermal or entropy corrections (Dudev include them). The water counts are balanced, so they
largely cancel, but say so with any number.

  python site_exchange.py build --frame LA/umbrella_res_035/lastframe.pdb --residues 156:220:240:300 --out site_la
  python site_exchange.py run --site site_la --aqua coord_check/geom --basis-svp def2-svp_lcecp1.nw \
      --basis-tzvp def2-tzvp_lcecp1.nw --cores 32
  python site_exchange.py calibrate --aqua coord_check/geom --basis-tzvp def2-tzvp_lcecp1.nw --out calib

`calibrate` runs the same PCM-water protocol on the acetate complexes of carboxylate_exchange.py
(M(H2O)n-k(OAc) + k H2O vs M(H2O)n + OAc-) for comparison with experiment: Ln - Ca about
-1.4 to -1.7 kcal/mol, La - Yb about +0.2 (c4-tuning.md). It is criterion 4 of the doc.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np

HARTREE_KCAL = 627.509474
IONS = {"CA": ("Ca", 2, 8), "LA": ("La", 3, 9), "LU": ("Lu", 3, 8)}   # element, charge, aqua waters
CHARGED = {"GLU": -1, "ASP": -1, "LYS": 1, "ARG": 1}
BACKBONE = {"N", "H", "C", "O", "H1", "H2", "H3", "OXT"}
SOLVENTS = {"eps4": 4.0, "eps29": 29.0}          # site dielectrics (Dudev)
WATER = {"water": 78.3553}
# cavity radii (A) for the ions; Shannon ionic radius + 1.0 A water-contact shell is a common
# choice; they mostly sit inside a water shell, so the result is not sensitive to them
DEFAULT_RADII = {"Ca": 2.0, "La": 2.2, "Lu": 1.98}


# ==============================
# Site cluster
# ==============================

def read_frame(pdb_file):
    from openmm.app import PDBFile
    from mcube.preprocess import pdb_resseq
    pdb = PDBFile(str(pdb_file))
    xyz = np.array(pdb.getPositions()._value) * 10.0
    return pdb.topology, xyz, pdb_resseq


def build(args):
    topology, xyz, pdb_resseq = read_frame(args.frame)
    atoms = list(topology.atoms())
    ions = [a for a in atoms if a.residue.name in IONS and len(list(a.residue.atoms())) == 1]
    if len(ions) != 1:
        raise SystemExit(f"expected one ion in {args.frame}, found {len(ions)}")
    ion = ions[0]
    r_ion = xyz[ion.index]
    # the frame is wrapped atom by atom: bring everything to its image nearest the ion, which
    # keeps each residue of the cluster whole (all of them lie well within half a box of it)
    cell = topology.getUnitCellDimensions()
    if cell is not None:
        box = np.array([cell[i]._value for i in range(3)]) * 10.0
        d = xyz - r_ion
        xyz = r_ion + d - box * np.round(d / box)
    wanted = {int(r) - args.offset for r in args.residues.replace(",", ":").split(":") if r}

    def near(res):
        return any(at.element is not None and at.element.symbol == "O"
                   and at.name not in BACKBONE
                   and np.linalg.norm(xyz[at.index] - r_ion) < args.shell for at in res.atoms())

    residues, waters = [], []
    for res in topology.residues():
        if res.name in ("HOH", "WAT"):
            o = [a for a in res.atoms() if a.element.symbol == "O"][0]
            if np.linalg.norm(xyz[o.index] - r_ion) < args.shell:
                waters.append(res)
        elif res.chain.id == ion.residue.chain.id:
            continue
        elif pdb_resseq(res.id) in wanted or (res.name not in ("POPE",) and near(res)):
            if any(a.name == "CA" for a in res.atoms()):
                residues.append(res)

    elems, coords, frozen, info = [], [], [], []
    elem, charge, _ = IONS[ion.residue.name]
    elems.append(elem); coords.append(r_ion)
    total = 0
    for res in residues:
        coordinating = near(res)
        num = pdb_resseq(res.id) + args.offset
        info.append({"residue": f"{res.name}{num}", "coordinating": coordinating,
                     "listed": pdb_resseq(res.id) in wanted})
        total += CHARGED.get(res.name, 0)
        byname = {a.name: a for a in res.atoms()}
        ca = byname["CA"]
        for a in res.atoms():
            if a.name in BACKBONE:
                continue
            idx = len(elems)
            elems.append(a.element.symbol); coords.append(xyz[a.index])
            if not coordinating or a.name in ("CA", "HA", "HA2", "HA3"):
                frozen.append(idx)
        for partner in ("N", "C"):        # H caps on CA, 1.09 A along the cut bonds
            if partner in byname:
                v = xyz[byname[partner].index] - xyz[ca.index]
                idx = len(elems)
                elems.append("H"); coords.append(xyz[ca.index] + 1.09 * v / np.linalg.norm(v))
                frozen.append(idx)
    for w in waters:
        for a in w.atoms():
            elems.append(a.element.symbol); coords.append(xyz[a.index])
    missing = wanted - {pdb_resseq(r.id) for r in residues}
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    meta = {"frame": str(Path(args.frame).resolve()), "source_ion": ion.residue.name,
            "site_charge": total, "residues": info, "waters": len(waters),
            "missing_residues": sorted(m + args.offset for m in missing),
            "frozen": frozen, "n_atoms": len(elems), "shell_A": args.shell}
    write_xyz(out / "site.xyz", elems, np.array(coords), json.dumps({"ion": ion.residue.name}))
    (out / "site.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta | {"frozen": f"{len(frozen)} atoms"}, indent=1))


# ==============================
# DFT
# ==============================

def make_mol(elems, xyz, charge, basis_file):
    from pyscf import gto
    basis, ecp = {}, {}
    for e in set(elems):
        basis[e] = gto.basis.load(basis_file, e)
        try:
            entry = gto.basis.load_ecp(basis_file, e)
        except Exception:
            entry = None
        if entry:
            ecp[e] = entry
    return gto.M(atom=[(e, tuple(x)) for e, x in zip(elems, xyz)], charge=charge, spin=0,
                 basis=basis, ecp=ecp, unit="Angstrom", verbose=0,
                 max_memory=int(os.environ.get("PYSCF_MAX_MEMORY", 60000)))


def scf(mol, functional, eps=None, radii=None):
    from pyscf import dft
    mf = dft.RKS(mol)
    mf.xc = functional
    mf.grids.level = 3
    mf.conv_tol = 1e-9
    mf.max_cycle = 300
    mf.damp = 0.3
    mf.level_shift = 0.2
    mf = mf.density_fit()
    if eps:
        from pyscf.data import radii as R
        from pyscf.data.elements import charge as Z
        from pyscf.solvent import pcm
        mf = mf.PCM()
        mf.with_solvent.method = "IEF-PCM"
        mf.with_solvent.eps = float(eps)
        if radii:
            # PCM's default table (1.2 x modified Bondi, Bohr, indexed by atomic number; the
            # cavity is built by element, so ECP core charges do not shift the lookup), with
            # the ions' entries replaced by the given cavity radii
            table = 1.2 * pcm.modified_Bondi.copy()
            for e, r in radii.items():
                table[Z(e)] = r / R.BOHR
            mf.with_solvent.radii_table = table
    return mf


def energy(elems, xyz, charge, basis_file, functional, eps=None, radii=None):
    mol = make_mol(elems, xyz, charge, basis_file)
    mf = scf(mol, functional, eps, radii)
    e = mf.kernel()
    if not mf.converged:
        mf = scf(mol, functional, eps, radii).newton()
        e = mf.kernel()
    return float(e) * HARTREE_KCAL, bool(mf.converged)


def optimise(elems, xyz, charge, frozen, basis_file, functional, maxsteps, workdir,
             eps=None, radii=None):
    """Constrained optimisation, in the gas phase or (eps) in IEF-PCM. The ring clusters are
    anionic (four Glu-, two of them uncoordinated), and in the gas phase their SCF flips
    between solutions from step to step (C22); a dielectric binds the extra charge."""
    from pyscf.geomopt.geometric_solver import kernel as geomopt
    mol = make_mol(elems, xyz, charge, basis_file)
    cons = None
    if frozen:
        cons = Path(workdir) / "constraints.txt"
        cons.write_text("$freeze\nxyz " + ",".join(str(i + 1) for i in frozen) + "\n")
    mf = scf(mol, functional, eps, radii)
    try:
        conv, mol = geomopt(mf, maxsteps=maxsteps, constraints=str(cons) if cons else None)
    except RuntimeError as err:
        print(f"[WARN] optimisation stopped: {err}; retrying with second-order SCF", flush=True)
        conv, mol = geomopt(scf(mol, functional, eps, radii).newton(), maxsteps=maxsteps,
                            constraints=str(cons) if cons else None)
    return bool(conv), mol.atom_coords(unit="Angstrom")


# ==============================
# Actions
# ==============================

def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    rows = [l.split() for l in lines[2:2 + n]]
    return [r[0] for r in rows], np.array([[float(v) for v in r[1:4]] for r in rows])


def write_xyz(path, elems, xyz, comment=""):
    with open(path, "w") as fh:
        fh.write(f"{len(elems)}\n{comment}\n")
        for e, (x, y, z) in zip(elems, xyz):
            fh.write(f"{e} {x:.6f} {y:.6f} {z:.6f}\n")


def load_report(path):
    return json.loads(path.read_text()) if path.exists() else {}


def save(path, report):
    path.write_text(json.dumps(report, indent=1))


def solution_energies(report, key, elems, xyz, charge, args, solvents):
    """Gas and PCM single points in the large basis, cached in report[key]."""
    row = report.setdefault(key, {})
    for label, solvent in [("gas", None)] + list(solvents.items()):
        if label in row and row[label]["converged"]:
            continue
        t0 = time.time()
        e, ok = energy(elems, xyz, charge, args.basis_tzvp, args.functional, solvent, args.radii)
        row[label] = {"E": round(e, 3), "converged": ok, "seconds": round(time.time() - t0, 1)}
        print(f"[{key} {label}] {row[label]}", flush=True)
    return row


def aqua_species(args, report, solvents=None):
    aq = Path(args.aqua)
    out = {}
    for ion, (elem, charge, n) in IONS.items():
        e, x = read_xyz(aq / f"{ion}_aqua.xyz")
        out[ion] = solution_energies(report, f"{ion}_aqua", e, x, charge, args, solvents or WATER)
    e, x = read_xyz(aq / "water.xyz")
    out["water"] = solution_energies(report, "water", e, x, 0, args, solvents or WATER)
    return out


def run(args):
    site = Path(args.site)
    meta = json.loads((site / "site.json").read_text())
    elems, xyz = read_xyz(site / "site.xyz")
    rpath = site / "report.json"
    report = load_report(rpath)
    report["meta"] = meta
    if report.setdefault("opt_eps", args.opt_eps) != args.opt_eps:
        raise SystemExit(f"{site} was optimised with opt_eps {report['opt_eps']}; use a fresh site directory")
    sp = report.setdefault("species", {})
    geom = site / "geom"; geom.mkdir(exist_ok=True)
    for ion, (elem, charge, _) in IONS.items():
        key = f"{ion}_site"
        e = [elem] + elems[1:]
        q = charge + meta["site_charge"]
        gfile = geom / f"{key}.xyz"
        if not gfile.exists():
            t0 = time.time()
            conv, x = optimise(e, xyz, q, meta["frozen"], args.basis_svp, args.functional,
                               args.maxsteps, site, args.opt_eps, args.radii)
            write_xyz(gfile, e, x, json.dumps({"converged": conv}))
            sp.setdefault(key, {})["opt"] = {"converged": conv, "seconds": round(time.time() - t0, 1)}
            save(rpath, report)
        e, x = read_xyz(gfile)
        d = sorted(np.linalg.norm(x[1:] - x[0], axis=1))
        sp.setdefault(key, {})["first_shell_A"] = [round(float(v), 2) for v in d[:10] if v < 3.0]
        solution_energies(sp, key, e, x, q, args, SOLVENTS)
        save(rpath, report)
    aq = aqua_species(args, sp)
    save(rpath, report)

    rx = {}
    for m1, m2 in (("CA", "LA"), ("CA", "LU"), ("LA", "LU")):
        n1, n2 = IONS[m1][2], IONS[m2][2]
        row = {}
        for label in ["gas"] + list(SOLVENTS):
            s1, s2 = sp[f"{m1}_site"][label], sp[f"{m2}_site"][label]
            a1, a2, w = (aq[m1]["gas" if label == "gas" else "water"], aq[m2]["gas" if label == "gas" else "water"],
                         aq["water"]["gas" if label == "gas" else "water"])
            if not all(t["converged"] for t in (s1, s2, a1, a2, w)):
                row[label] = None
                continue
            row[label] = round(s2["E"] + a1["E"] + (n2 - n1) * w["E"] - s1["E"] - a2["E"], 2)
        rx[f"{m1}->{m2}"] = row
    report["exchange_kcal"] = rx
    report["note"] = ("negative = site prefers the second ion; eps4/eps29 = site in IEF-PCM eps 4/29, "
                      "aqua ions and water in IEF-PCM water; electrostatics only; no thermal corrections")
    save(rpath, report)
    print(json.dumps({"meta": {k: meta[k] for k in ("source_ion", "site_charge", "residues", "waters")},
                      "exchange_kcal": rx}, indent=1))


def calibrate(args):
    """Acetate binding in PCM water with the same protocol, against experiment (criterion 4)."""
    aq = Path(args.aqua)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    rpath = out / "report.json"
    report = load_report(rpath)
    sp = report.setdefault("species", {})
    water = WATER
    a = aqua_species(args, sp, water)
    e, x = read_xyz(aq / "acetate.xyz")
    ace = solution_energies(sp, "acetate", e, x, -1, args, water)
    save(rpath, report)
    bind = {}
    for ion, (elem, charge, n) in IONS.items():
        for k, nw in (("w-1", 1), ("w-2", 2)):
            f = aq / f"{ion}_{k}.xyz"
            if not f.exists():
                continue
            e, x = read_xyz(f)
            c = solution_energies(sp, f"{ion}_{k}", e, x, charge - 1, args, water)
            save(rpath, report)
            vals = {}
            for label in ("gas", "water"):
                parts = (c[label], a["water"][label], a[ion][label], ace[label])
                vals[label] = (round(c[label]["E"] + nw * a["water"][label]["E"] - a[ion][label]["E"]
                                     - ace[label]["E"], 2) if all(p["converged"] for p in parts) else None)
            bind[f"{ion}_{k}"] = vals
    best = {ion: min((v["water"] for k, v in bind.items() if k.startswith(ion) and v["water"] is not None),
                     default=None) for ion in IONS}
    rel = {}
    for a1, b1 in (("LA", "CA"), ("LU", "CA"), ("LU", "LA")):
        if best[a1] is not None and best[b1] is not None:
            rel[f"{a1}-{b1}"] = round(best[a1] - best[b1], 2)
    report["binding_kcal"], report["best_mode_water"], report["relative_water"] = bind, best, rel
    report["experiment"] = {"LA-CA": "-1.4 to -1.7", "LU-LA (Yb stands in)": "about -0.2"}
    save(rpath, report)
    print(json.dumps({"binding_kcal": bind, "relative_water": rel, "experiment": report["experiment"]}, indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="action", required=True)
    b = sub.add_parser("build")
    b.add_argument("--frame", required=True)
    b.add_argument("--residues", default="156:220:240:300", help="channel numbering, ':'-separated")
    b.add_argument("--offset", type=int, default=130, help="channel number = mcube number + offset")
    b.add_argument("--shell", type=float, default=3.0)
    b.add_argument("--out", required=True)
    for name in ("run", "calibrate"):
        q = sub.add_parser(name)
        q.add_argument("--aqua", required=True, help="carboxylate_exchange.py geom/ directory")
        q.add_argument("--basis-tzvp", required=True)
        q.add_argument("--functional", default="r2scan")
        q.add_argument("--cores", type=int, default=os.cpu_count())
        q.add_argument("--radii", default=",".join(f"{k}={v}" for k, v in DEFAULT_RADII.items()),
                       help="ion cavity radii in A, e.g. Ca=2.0,La=2.2,Lu=1.98")
        if name == "run":
            q.add_argument("--site", required=True)
            q.add_argument("--basis-svp", required=True)
            q.add_argument("--maxsteps", type=int, default=150)
            q.add_argument("--opt-eps", type=float, default=None,
                           help="optimise in IEF-PCM at this dielectric (default: gas phase)")
        else:
            q.add_argument("--out", required=True)
    args = p.parse_args()
    if args.action == "build":
        return build(args)
    from pyscf import lib
    lib.num_threads(args.cores)
    args.radii = {k: float(v) for k, v in (kv.split("=") for kv in args.radii.split(",") if kv)}
    (run if args.action == "run" else calibrate)(args)


if __name__ == "__main__":
    main()
