#!/usr/bin/env python3
"""
Is g-xTB a better QM method than GFN2-xTB for an Ln3+ aqua ion? A test on clusters,
which needs no point-charge embedding (g-xTB has none, decision log F7).

GFN2-xTB/MM of La3+ and Lu3+ in water lost SCC convergence after 6-9 ps, with the
HOMO-LUMO gap falling to ~0.001 Eh, and with the second shell in QM it pulled 11-12
waters into the first shell (F4, F7). This takes frames from those failing xTB/MM
trajectories, cuts the ion out with its first shell (S1) and with its first and second
shells (S2), and compares, against r2SCAN/def2-SVP (the DFT level of C16):

  - the HOMO-LUMO gap: does it close in GFN2 or g-xTB where r2SCAN keeps it open?
  - the gradient: RMS error against r2SCAN, and the mean radial force on the first-shell
    oxygens (positive = pushed away from the ion). A method that pulls the shell in
    too tightly shows a more negative radial force than r2SCAN on the same geometry.
  - for S1 also in the field of the rest of the box as bare point charges (TIP3P), for
    GFN2 and r2SCAN only: is the gap collapse a feature of the embedded problem itself,
    or of GFN2. A g-xTB run with point charges is kept as a check: the g-xTB build of
    xtb (v2.0.1) accepts the charges and silently ignores them (same energy to every digit).
  - S1 optimised in the gas phase with each method: Ln-O distances.

  python cluster_test.py --run LA=/scratch/.../LA_1264 --run LU=/scratch/.../LU_1264 \\
      --gxtb /scratch/user/uquadhik/opt/gxtb/xtb-6.7.1/bin/xtb \\
      --basis-file /scratch/user/uquadhik/channel/dft/def2-svp.nw --out gxtb_test

Each run directory needs system.pdb and qmmm/QM_MM.dcd (the xTB/MM leg of bulk.slurm).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

ION = {"CA": ("Ca", 2), "LA": ("La", 3), "LU": ("Lu", 3)}
TIP3P = {"O": -0.834, "H": 0.417}
BOHR = 0.529177210903          # A
HARTREE_BOHR_TO_KCAL_A = 627.509474 / BOHR
EV = 27.211386


# ==============================
# Clusters
# ==============================

def frame_clusters(xyz, box, ion, waters, cut1, cut2, env):
    """
    Ion-centred clusters from one frame (A). Each water is moved whole by its oxygen's
    minimum image. Returns S1, S2 and, for S1, the point charges of every other water
    whose O lies within `env` A.
    """
    centre = xyz[ion]
    shift = lambda r: r - box * np.round(r / box)
    elems_by_water, coords_by_water, r_o = [], [], []
    for o, h1, h2 in waters:
        d = shift(xyz[o] - centre)
        move = d - (xyz[o] - centre)
        coords_by_water.append(np.array([xyz[o], xyz[h1], xyz[h2]]) - centre + move)
        r_o.append(np.linalg.norm(d))
    r_o = np.array(r_o)
    order = np.argsort(r_o)

    def cluster(cut):
        chosen = [i for i in order if r_o[i] < cut]
        elems = [ION_SYMBOL] + ["O", "H", "H"] * len(chosen)
        coords = np.vstack([np.zeros((1, 3))] + [coords_by_water[i] for i in chosen])
        return {"elems": elems, "xyz": coords, "n_water": len(chosen),
                "first_shell_O": [1 + 3 * k for k, i in enumerate(chosen) if r_o[i] < cut1]}

    s1, s2 = cluster(cut1), cluster(cut2)
    pc = [i for i in order if cut1 <= r_o[i] < env]
    s1["pc_xyz"] = np.vstack([coords_by_water[i] for i in pc])
    s1["pc_q"] = np.array([TIP3P["O"], TIP3P["H"], TIP3P["H"]] * len(pc))
    return s1, s2


def write_xyz(path, elems, xyz):
    with open(path, "w") as fh:
        fh.write(f"{len(elems)}\n\n")
        for e, (x, y, z) in zip(elems, xyz):
            fh.write(f"{e} {x:.6f} {y:.6f} {z:.6f}\n")


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    rows = [l.split() for l in lines[2:2 + n]]
    return [r[0] for r in rows], np.array([[float(v) for v in r[1:4]] for r in rows])


# ==============================
# Methods
# ==============================

def run_xtb(binary, flags, cl, charge, work, cores, pc=False, opt=False):
    """One xtb run (GFN2 or g-xTB by `flags`): energy, gap, gradient (Eh/bohr), convergence."""
    work = Path(work)
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    write_xyz(work / "c.xyz", cl["elems"], cl["xyz"])
    cmd = [binary, "c.xyz", *flags, "--chrg", str(charge), "--uhf", "0",
           "--opt" if opt else "--grad", "--parallel", str(cores)]
    if pc:
        # the way ASH embeds: ORCA-style pcharge file in A, hardness 1000 (bare charges)
        with open(work / "pcharge", "w") as fh:
            fh.write(f"{len(cl['pc_q'])}\n")
            for q, (x, y, z) in zip(cl["pc_q"], cl["pc_xyz"]):
                fh.write(f"{q} {x} {y} {z} 1000\n")
        (work / "xtbinput").write_text("$embedding\ninterface=orca\nend\n")
        cmd += ["--input", "xtbinput"]
    env = dict(os.environ, OMP_NUM_THREADS=str(cores), OMP_STACKSIZE="2G")
    t0 = time.time()
    with open(work / "xtb.out", "w") as fh:
        rc = subprocess.run(cmd, cwd=work, stdout=fh, stderr=subprocess.STDOUT, env=env).returncode
    out = (work / "xtb.out").read_text(errors="replace")
    res = {"rc": rc, "seconds": round(time.time() - t0, 1),
           "converged": rc == 0 and "did not converge" not in out}
    gaps = re.findall(r"HL-Gap\s+(-?[\d.]+)\s+Eh", out)
    res["gap_eV"] = round(float(gaps[-1]) * EV, 4) if gaps else None
    energies = re.findall(r"TOTAL ENERGY\s+(-?[\d.]+)\s+Eh", out)
    res["energy_Eh"] = float(energies[-1]) if energies else None
    # GFN2 prints ":: number of electrons  72 ::", g-xTB "number of electrons  7.19...E+01 e"
    nel = re.findall(r"number of electrons\s+([\d.Ee+-]+)", out, flags=re.I)
    res["electrons"] = int(round(float(nel[0]))) if nel else None
    res["pc_mentioned"] = bool(re.search(r"point charge|pcem|embedding", out, flags=re.I))
    if opt:
        if (work / "xtbopt.xyz").exists():
            res["opt"] = read_xyz(work / "xtbopt.xyz")
    elif (work / "gradient").exists():
        res["grad"] = read_turbomole_gradient(work / "gradient", len(cl["elems"]))
    return res


def read_turbomole_gradient(path, n):
    lines = Path(path).read_text().splitlines()
    start = max(i for i, l in enumerate(lines) if "cycle" in l)
    rows = lines[start + 1 + n:start + 1 + 2 * n]
    return np.array([[float(v.replace("D", "E")) for v in r.split()[:3]] for r in rows])


def pyscf_mol(cl, charge, basis_file):
    from pyscf import gto
    elems = sorted(set(cl["elems"]))
    basis = {e: gto.basis.load(basis_file, e) for e in elems}
    ecp = {}
    for e in elems:
        try:
            ecp[e] = gto.basis.load_ecp(basis_file, e)
        except Exception:
            pass
    atom = [(e, tuple(x)) for e, x in zip(cl["elems"], cl["xyz"])]
    return gto.M(atom=atom, charge=charge, spin=0, basis=basis, ecp=ecp, unit="Angstrom",
                 verbose=0, max_memory=int(os.environ.get("PYSCF_MAX_MEMORY", 60000)))


def run_dft(cl, charge, basis_file, cores, functional, pc=False):
    from pyscf import dft, lib, qmmm
    lib.num_threads(cores)
    mol = pyscf_mol(cl, charge, basis_file)
    mf = dft.RKS(mol)
    mf.xc = functional
    mf.grids.level = 3
    mf.conv_tol = 1e-8
    mf = mf.density_fit()
    if pc:
        mf = qmmm.mm_charge(mf, cl["pc_xyz"], cl["pc_q"])
    t0 = time.time()
    energy = mf.kernel()
    occ, e = mf.mo_occ, mf.mo_energy
    res = {"converged": bool(mf.converged), "energy_Eh": float(energy),
           "gap_eV": round(float(e[occ == 0].min() - e[occ > 0].max()) * EV, 4),
           "electrons": int(mol.nelectron)}
    res["grad"] = np.asarray(mf.nuc_grad_method().kernel())
    res["seconds"] = round(time.time() - t0, 1)
    return res


def opt_dft(cl, charge, basis_file, cores, functional):
    from pyscf import dft, lib
    from pyscf.geomopt.geometric_solver import optimize
    lib.num_threads(cores)
    mf = dft.RKS(pyscf_mol(cl, charge, basis_file))
    mf.xc = functional
    mf.grids.level = 3
    mf = mf.density_fit()
    t0 = time.time()
    mol = optimize(mf, maxsteps=150)
    return {"opt": (cl["elems"], mol.atom_coords(unit="Angstrom")),
            "seconds": round(time.time() - t0, 1)}


# ==============================
# Comparison
# ==============================

def radial_force(cl, grad):
    """Mean force on the first-shell O along ion->O (kcal/mol/A); + = pushed out."""
    f = -np.asarray(grad) * HARTREE_BOHR_TO_KCAL_A
    out = []
    for o in cl["first_shell_O"]:
        u = cl["xyz"][o] - cl["xyz"][0]
        out.append(float(f[o] @ (u / np.linalg.norm(u))))
    return round(float(np.mean(out)), 2)


def summarise(cl, res, ref):
    row = {k: res.get(k) for k in ("converged", "gap_eV", "electrons", "seconds", "rc",
                                   "energy_Eh", "pc_mentioned") if k in res}
    if res.get("grad") is not None:
        row["radial_force_O1"] = radial_force(cl, res["grad"])
        if ref is not None and ref.get("grad") is not None:
            d = (np.asarray(res["grad"]) - ref["grad"]) * HARTREE_BOHR_TO_KCAL_A
            row["grad_rmse_vs_dft"] = round(float(np.sqrt((d ** 2).sum(1).mean())), 2)
            core = [0] + cl["first_shell_O"]
            row["grad_rmse_ion_O1"] = round(float(np.sqrt((d[core] ** 2).sum(1).mean())), 2)
    return row


def ln_o(elems, xyz, cut=3.3):
    r = np.linalg.norm(np.asarray(xyz)[1:] - np.asarray(xyz)[0], axis=1)
    r_o = np.sort([d for e, d in zip(elems[1:], r) if e == "O"])
    inner = r_o[r_o < cut]
    return {"cn": int(len(inner)), "mean_ln_o": round(float(inner.mean()), 3),
            "ln_o": [round(float(v), 3) for v in r_o[:len(inner) + 1]]}


# ==============================
# Driver
# ==============================

def main():
    global ION_SYMBOL
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", action="append", required=True, help="ION=run_dir, e.g. LA=.../LA_1264")
    p.add_argument("--qmdir", default="qmmm")
    p.add_argument("--dt", type=float, default=0.02, help="ps between saved frames (bulk.slurm DT)")
    p.add_argument("--gxtb", required=True, help="the g-xTB build of xtb")
    p.add_argument("--gfn2", default=shutil.which("xtb"), help="xtb for GFN2 (default: on PATH)")
    p.add_argument("--basis-file", required=True)
    p.add_argument("--functional", default="r2scan")
    p.add_argument("--frames", type=int, default=6, help="frames per run, the last included")
    p.add_argument("--cut1", type=float, default=3.3)
    p.add_argument("--cut2", type=float, default=5.1)
    p.add_argument("--env", type=float, default=14.0, help="point-charge radius for S1 (A)")
    p.add_argument("--skip-dft", action="store_true")
    p.add_argument("--skip-s2", action="store_true")
    p.add_argument("--skip-opt", action="store_true")
    p.add_argument("--cores", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    p.add_argument("--out", required=True)
    args = p.parse_args()

    import mdtraj
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    report = {"settings": {k: v for k, v in vars(args).items()}, "runs": {}}
    gxtb, gfn2 = ["--gxtb"], ["--gfn", "2"]

    def save():
        (out / "report.json").write_text(json.dumps(report, indent=1, default=str))

    for spec in args.run:
        name, run_dir = spec.split("=", 1)
        ION_SYMBOL, charge = ION[name]
        run_dir = Path(run_dir)
        traj = mdtraj.load(str(run_dir / args.qmdir / "QM_MM.dcd"), top=str(run_dir / "system.pdb"))
        t = traj.topology
        ion = next(a.index for a in t.atoms if a.residue.name == name)
        waters = []
        for res in t.residues:
            if res.name == "HOH":
                o = [a.index for a in res.atoms if a.element.symbol == "O"]
                h = [a.index for a in res.atoms if a.element.symbol == "H"]
                waters.append((o[0], h[0], h[1]))
        picks = sorted(set(np.linspace(10, traj.n_frames - 1, args.frames).round().astype(int)))
        runs = report["runs"][f"{name}:{run_dir.name}"] = {"n_frames": int(traj.n_frames),
                                                         "frames": {}}
        for k in picks:
            xyz = traj.xyz[k] * 10.0
            box = traj.unitcell_lengths[k] * 10.0
            s1, s2 = frame_clusters(xyz, box, ion, waters, args.cut1, args.cut2, args.env)
            w = out / f"{name}_{run_dir.name}_f{k}"
            row = runs["frames"][int(k)] = {"time_ps": round(k * args.dt, 2),
                                           "S1_waters": s1["n_water"], "S2_waters": s2["n_water"],
                                           "S1_ln_o": ln_o(s1["elems"], s1["xyz"], args.cut1)}
            ref1 = None if args.skip_dft else run_dft(s1, charge, args.basis_file, args.cores, args.functional)
            ref1pc = None if args.skip_dft else run_dft(s1, charge, args.basis_file, args.cores,
                                                        args.functional, pc=True)
            row["S1"] = {"dft": summarise(s1, ref1, None) if ref1 else None,
                         "gfn2": summarise(s1, run_xtb(args.gfn2, gfn2, s1, charge, w / "s1_gfn2", args.cores), ref1),
                         "gxtb": summarise(s1, run_xtb(args.gxtb, gxtb, s1, charge, w / "s1_gxtb", args.cores), ref1)}
            row["S1_pc"] = {"dft": summarise(s1, ref1pc, None) if ref1pc else None,
                            "gfn2": summarise(s1, run_xtb(args.gfn2, gfn2, s1, charge, w / "s1pc_gfn2",
                                                          args.cores, pc=True), ref1pc),
                            "gxtb": summarise(s1, run_xtb(args.gxtb, gxtb, s1, charge, w / "s1pc_gxtb",
                                                          args.cores, pc=True), None)}
            if not args.skip_s2:
                ref2 = None if args.skip_dft else run_dft(s2, charge, args.basis_file, args.cores, args.functional)
                row["S2"] = {"dft": summarise(s2, ref2, None) if ref2 else None,
                             "gfn2": summarise(s2, run_xtb(args.gfn2, gfn2, s2, charge, w / "s2_gfn2", args.cores), ref2),
                             "gxtb": summarise(s2, run_xtb(args.gxtb, gxtb, s2, charge, w / "s2_gxtb", args.cores), ref2)}
            print(f"[{name} frame {k}]", json.dumps({c: {m: (v or {}).get("gap_eV") for m, v in row[c].items()}
                                                     for c in ("S1", "S1_pc", "S2") if c in row}), flush=True)
            save()

        if not args.skip_opt:
            # the gas-phase S1 of the first picked frame, optimised by each method
            xyz = traj.xyz[picks[0]] * 10.0
            s1, _ = frame_clusters(xyz, traj.unitcell_lengths[picks[0]] * 10.0, ion, waters,
                                   args.cut1, args.cut2, args.env)
            w = out / f"{name}_{run_dir.name}_opt"
            opts = {"gfn2": run_xtb(args.gfn2, gfn2, s1, charge, w / "gfn2", args.cores, opt=True),
                    "gxtb": run_xtb(args.gxtb, gxtb, s1, charge, w / "gxtb", args.cores, opt=True)}
            if not args.skip_dft:
                try:
                    opts["dft"] = opt_dft(s1, charge, args.basis_file, args.cores, args.functional)
                except Exception as err:          # e.g. geomeTRIC missing; the rest stands
                    opts["dft"] = {"error": repr(err)}
            runs["opt_S1"] = {"waters": s1["n_water"]}
            for m, r in opts.items():
                runs["opt_S1"][m] = ({**ln_o(*r["opt"], cut=args.cut1), "seconds": r.get("seconds")}
                                     if r.get("opt") else {"failed": True, "rc": r.get("rc")})
            print(f"[{name} opt]", json.dumps(runs["opt_S1"]), flush=True)
            save()
    save()
    print(json.dumps(report, indent=1, default=str))


if __name__ == "__main__":
    main()
