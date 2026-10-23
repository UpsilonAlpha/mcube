#!/usr/bin/env python3
"""
Does 12-6-4 rank Ca2+, La3+ and Lu3+ the way r2SCAN does when a carboxylate replaces
first-shell water? Glu300, the conserved carboxylate at the channel's constriction, is the
likely site, and an MM PMF is only as good as this exchange.

For each ion, the gas-phase reactions

  w-1:  M(H2O)n + OAc-  ->  M(H2O)n-1(OAc) + H2O
  w-2:  M(H2O)n + OAc-  ->  M(H2O)n-2(OAc) + 2 H2O

(n = 8 Ca, 9 La, 8 Lu). The acetate starts one-oxygen (w-1) or two-oxygen (w-2) bound, but
each method settles it where it likes; in 12-6-4 it closes to two-oxygen binding in both, so
the report gives the denticity each method ends with (acetate O within 2.9 A of the ion).

with every species optimised by each method on its own:

  MM     amber14 TIP3P (rigid) + acetate built from ff14SB's Asp side chain (types CT/HC/C/O2,
         Asp's side-chain charges with the methyl carbon taking the remainder to -1), the ion
         with mcube's metal parameters (12-6-4 by default, --set iod for the 12-6 set)
  DFT    r2SCAN/def2-SVP, density fitting, def2 ECPs (basis and ECP from a BSE file, as in
         dft/bse_basis.py), geomeTRIC, starting from the MM minimum

What matters is not the absolute reaction energy (gas-phase ion-anion Coulomb dominates it,
and def2-SVP has BSSE) but the differences between ions: does MM get La - Ca, Lu - La and
the w-1/w-2 preference (how many waters stay) with the same sign and roughly the same size as r2SCAN?
A proton moving from a water to the acetate in the DFT optimisation is flagged, because
then the two methods are no longer describing the same species.

  python carboxylate_exchange.py --basis-file def2-svp.nw --out coord_check
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np

HARTREE_KCAL = 627.509474
ION = {"CA": ("Ca", 2, 2.46, 8), "LA": ("La", 3, 2.54, 9), "LU": ("Lu", 3, 2.33, 8)}

ACETATE_XML = """<ForceField>
 <Residues>
  <Residue name="ACT">
   <Atom name="C1" type="protein-CT" charge="-0.2332"/>
   <Atom name="H1" type="protein-HC" charge="0.0122"/>
   <Atom name="H2" type="protein-HC" charge="0.0122"/>
   <Atom name="H3" type="protein-HC" charge="0.0122"/>
   <Atom name="C2" type="protein-C" charge="0.7994"/>
   <Atom name="O1" type="protein-O2" charge="-0.8014"/>
   <Atom name="O2" type="protein-O2" charge="-0.8014"/>
   <Bond atomName1="C1" atomName2="H1"/>
   <Bond atomName1="C1" atomName2="H2"/>
   <Bond atomName1="C1" atomName2="H3"/>
   <Bond atomName1="C1" atomName2="C2"/>
   <Bond atomName1="C2" atomName2="O1"/>
   <Bond atomName1="C2" atomName2="O2"/>
  </Residue>
 </Residues>
</ForceField>
"""
FORCEFIELDS = ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]


# ==============================
# Starting structures
# ==============================

def sphere(n):
    """n roughly even directions (Fibonacci sphere), the first at the north pole."""
    i = np.arange(n) + 0.5
    z = 1 - 2 * i / n
    phi = np.pi * (1 + 5 ** 0.5) * i
    pts = np.stack([np.sqrt(1 - z ** 2) * np.cos(phi), np.sqrt(1 - z ** 2) * np.sin(phi), z], 1)
    return pts[np.argsort(-pts[:, 2])]


def frame(u):
    u = u / np.linalg.norm(u)
    a = np.array([1.0, 0, 0]) if abs(u[0]) < 0.9 else np.array([0, 1.0, 0])
    v = np.cross(u, a)
    v /= np.linalg.norm(v)
    return u, v, np.cross(u, v)


def water_at(u, r):
    """Water with O at r*u, both H pointing away from the ion (HOH 104.52, OH 0.9572)."""
    u, v, _ = frame(u)
    o = r * u
    half = np.radians(104.52 / 2)
    h1 = o + 0.9572 * (np.cos(half) * u + np.sin(half) * v)
    h2 = o + 0.9572 * (np.cos(half) * u - np.sin(half) * v)
    return [("O", o), ("H", h1), ("H", h2)]


def acetate_at(u, r, bidentate):
    """Acetate on direction u: one O at distance r (k1), or both O at r (k2)."""
    u, v, w = frame(u)
    co, half = 1.26, np.radians(125.0 / 2)
    if bidentate:
        dc = co * np.cos(half) + np.sqrt(max(r ** 2 - (co * np.sin(half)) ** 2, 0.1))
        c2 = dc * u
        o1 = c2 - co * np.cos(half) * u + co * np.sin(half) * v
        o2 = c2 - co * np.cos(half) * u - co * np.sin(half) * v
        axis = u
    else:
        o1 = r * u
        c2 = o1 + co * u
        # O2 points sideways and outward
        d2 = np.cos(np.radians(180 - 125)) * u + np.sin(np.radians(180 - 125)) * v
        o2 = c2 + co * d2
        axis = -((o1 - c2) / co + (o2 - c2) / co)
        axis /= np.linalg.norm(axis)
    c1 = c2 + 1.52 * axis
    a2, b2, c3 = frame(axis)
    hs = [c1 + 1.09 * (np.cos(np.radians(70.5)) * a2 + np.sin(np.radians(70.5)) *
                       (np.cos(t) * b2 + np.sin(t) * c3)) for t in np.radians([0, 120, 240])]
    return [("C", c1), ("H", hs[0]), ("H", hs[1]), ("H", hs[2]), ("C", c2), ("O", o1), ("O", o2)]


def complex_atoms(ion, kind):
    sym, _, r0, n = ION[ion]
    dirs = sphere(n)
    atoms = [(sym, np.zeros(3))]
    if kind == "aqua":
        waters = dirs
    elif kind == "w-1":
        waters = dirs[1:]
    else:
        # the acetate takes the pole and the site nearest it
        waters = np.delete(dirs, [0, 1], axis=0)
    for d in waters:
        atoms += water_at(d, r0)
    if kind != "aqua":
        u = dirs[0] if kind == "w-1" else dirs[0] + dirs[1]
        atoms += acetate_at(u, r0, kind == "w-2")
    return atoms


# ==============================
# MM
# ==============================

def mm_topology(atoms, ion):
    from openmm.app import Element, Topology
    top = Topology()
    chain = top.addChain()
    i = 0
    xyz = []
    while i < len(atoms):
        sym = atoms[i][0]
        if sym in ("Ca", "La", "Lu"):
            res = top.addResidue(ion, chain)
            top.addAtom(ion, Element.getBySymbol(sym), res)
            xyz.append(atoms[i][1]); i += 1
        elif sym == "O" and i + 2 < len(atoms) and atoms[i + 1][0] == "H" and atoms[i + 2][0] == "H":
            res = top.addResidue("HOH", chain)
            o = top.addAtom("O", Element.getBySymbol("O"), res)
            h1 = top.addAtom("H1", Element.getBySymbol("H"), res)
            h2 = top.addAtom("H2", Element.getBySymbol("H"), res)
            top.addBond(o, h1); top.addBond(o, h2)
            xyz += [atoms[i][1], atoms[i + 1][1], atoms[i + 2][1]]; i += 3
        else:                                           # acetate: C1 H H H C2 O O
            res = top.addResidue("ACT", chain)
            names = ["C1", "H1", "H2", "H3", "C2", "O1", "O2"]
            made = [top.addAtom(nm, Element.getBySymbol(atoms[i + k][0]), res) for k, nm in enumerate(names)]
            for a, b in ((0, 1), (0, 2), (0, 3), (0, 4), (4, 5), (4, 6)):
                top.addBond(made[a], made[b])
            xyz += [atoms[i + k][1] for k in range(7)]; i += 7
    return top, np.array(xyz)


def mm_minimise(atoms, ion, which):
    import openmm
    from openmm import unit
    from openmm.app import ForceField, NoCutoff, Simulation
    from mcube import ions as mions
    xml = Path("acetate_ff14sb.xml")
    xml.write_text(ACETATE_XML)
    ffs = FORCEFIELDS + [str(xml.resolve())]
    top, xyz = mm_topology(atoms, ion)
    system = ForceField(*ffs).createSystem(top, nonbondedMethod=NoCutoff, constraints=None,
                                          rigidWater=True)
    if any(a.element is not None and a.element.symbol in ("Ca", "La", "Lu") for a in top.atoms()):
        mions.apply_metal_parameters(system, top, [ion], which, forcefields=ffs)
    sim = Simulation(top, system, openmm.VerletIntegrator(0.001),
                     openmm.Platform.getPlatformByName("Reference"))
    sim.context.setPositions(xyz * 0.1 * unit.nanometer)
    openmm.LocalEnergyMinimizer.minimize(sim.context, 1e-4, 20000)
    st = sim.context.getState(getEnergy=True, getPositions=True)
    pos = st.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    e = st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    return float(e), [(a[0], np.array(p)) for a, p in zip(atoms, pos)]


# ==============================
# DFT
# ==============================

def dft_optimise(atoms, charge, basis_file, cores, functional, maxsteps):
    from pyscf import dft, gto, lib
    from pyscf.geomopt.geometric_solver import kernel as geomopt
    lib.num_threads(cores)
    elems = sorted({a[0] for a in atoms})
    basis = {e: gto.basis.load(basis_file, e) for e in elems}
    ecp = {}
    for e in elems:
        try:
            entry = gto.basis.load_ecp(basis_file, e)
        except Exception:
            entry = None
        if entry:
            ecp[e] = entry
    mol = gto.M(atom=[(a[0], tuple(a[1])) for a in atoms], charge=charge, spin=0, basis=basis,
                ecp=ecp, unit="Angstrom", verbose=0,
                max_memory=int(os.environ.get("PYSCF_MAX_MEMORY", 60000)))
    mf = dft.RKS(mol)
    mf.xc = functional
    mf.grids.level = 3
    mf.conv_tol = 1e-9
    mf = mf.density_fit()
    if len(atoms) > 1:
        conv, mol = geomopt(mf, maxsteps=maxsteps)
    else:
        conv = True
    mf = dft.RKS(mol)
    mf.xc = functional
    mf.grids.level = 3
    mf = mf.density_fit()
    e = mf.kernel()
    xyz = mol.atom_coords(unit="Angstrom")
    return float(e) * HARTREE_KCAL, bool(conv and mf.converged), [(a[0], np.array(p)) for a, p in zip(atoms, xyz)]


def geometry_report(atoms):
    """Ion-O distances by ligand, and a flag if a proton has left its water."""
    xyz = np.array([a[1] for a in atoms])
    syms = [a[0] for a in atoms]
    out = {"ion_O_water": [], "ion_O_acetate": [], "acetate_denticity": None, "proton_moved": False}
    i = 1
    while i < len(atoms):
        if syms[i] == "O" and i + 2 < len(atoms) and syms[i + 1] == "H" and syms[i + 2] == "H":
            out["ion_O_water"].append(round(float(np.linalg.norm(xyz[i] - xyz[0])), 3))
            oh = max(np.linalg.norm(xyz[i + 1] - xyz[i]), np.linalg.norm(xyz[i + 2] - xyz[i]))
            out["proton_moved"] |= bool(oh > 1.25)
            i += 3
        else:
            out["ion_O_acetate"] = [round(float(np.linalg.norm(xyz[i + k] - xyz[0])), 3) for k in (5, 6)]
            out["acetate_denticity"] = sum(d < 2.9 for d in out["ion_O_acetate"])
            i += 7
    out["ion_O_water"].sort()
    return out


# ==============================
# Driver
# ==============================

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ions", default="CA,LA,LU")
    p.add_argument("--set", default="1264", help="metal parameter set for the MM side")
    p.add_argument("--basis-file", required=True)
    p.add_argument("--functional", default="r2scan")
    p.add_argument("--maxsteps", type=int, default=200)
    p.add_argument("--skip-dft", action="store_true")
    p.add_argument("--cores", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    p.add_argument("--out", required=True)
    args = p.parse_args()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    os.chdir(out)
    basis_file = str(Path(args.basis_file).resolve())
    report = {"settings": vars(args), "species": {}, "reactions": {}}

    def save():
        Path("report.json").write_text(json.dumps(report, indent=1, default=str))

    def run(name, atoms, ion, charge):
        t0 = time.time()
        e_mm, mm_atoms = mm_minimise(atoms, ion, args.set)
        row = {"E_mm": round(e_mm, 3), "mm_geometry": geometry_report(mm_atoms) if ion else None}
        if not args.skip_dft:
            e_dft, conv, dft_atoms = dft_optimise(mm_atoms, charge, basis_file, args.cores,
                                                  args.functional, args.maxsteps)
            row.update({"E_dft": round(e_dft, 3), "dft_converged": conv,
                        "dft_geometry": geometry_report(dft_atoms) if ion else None})
        row["seconds"] = round(time.time() - t0, 1)
        report["species"][name] = row
        print(f"[{name}]", json.dumps(row), flush=True)
        save()
        return row

    water = run("water", water_at(np.array([0, 0, 1.0]), 0.0), None, 0)
    acetate = run("acetate", acetate_at(np.array([0, 0, 1.0]), 0.0, False), None, -1)
    for ion in args.ions.split(","):
        q = ION[ion][1]
        sp = {k: run(f"{ion}_{k}", complex_atoms(ion, k), ion, q - (k != "aqua"))
              for k in ("aqua", "w-1", "w-2")}
        for k, nw in (("w-1", 1), ("w-2", 2)):
            rx = {}
            for m in ("E_mm", "E_dft"):
                if m in sp[k]:
                    rx[m.replace("E_", "dE_")] = round(sp[k][m] + nw * water[m] - sp["aqua"][m] - acetate[m], 2)
            report["reactions"][f"{ion}_{k}"] = rx
        save()

    # the differences that decide a relative comparison
    rel = {}
    for k in ("w-1", "w-2"):
        for a, b in (("LA", "CA"), ("LU", "CA"), ("LU", "LA")):
            ra, rb = report["reactions"].get(f"{a}_{k}"), report["reactions"].get(f"{b}_{k}")
            if ra and rb:
                rel[f"{k}: {a}-{b}"] = {m: round(ra[m] - rb[m], 2) for m in ra if m in rb}
    for ion in args.ions.split(","):
        r = report["reactions"]
        if f"{ion}_w-2" in r:
            rel[f"w-2 minus w-1: {ion}"] = {m: round(r[f"{ion}_w-2"][m] - r[f"{ion}_w-1"][m], 2)
                                            for m in r[f"{ion}_w-2"] if m in r[f"{ion}_w-1"]}
    report["relative"] = rel
    save()
    print(json.dumps({"reactions": report["reactions"], "relative": rel}, indent=1))


if __name__ == "__main__":
    main()
