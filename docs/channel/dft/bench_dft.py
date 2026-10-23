#!/usr/bin/env python3
"""
How long does one DFT/MM gradient take for an Ln3+ ion and its first hydration shell?
That decides how many ps a DFT/MM umbrella window over the desolvation step can afford.

Builds the QM region the way an umbrella QM/MM window does (ion + the shell preprocess
would detect), the MM system with mcube's make_mm and the same metal parameters, and the
DFT theory with mcube's make_dft (so the def2 ECP default is what is exercised). Times
`--repeat` DFT/MM gradients (the first includes the SCF from scratch; later ones may
start from the previous orbitals) and one xTB/MM gradient on the same system.

  python bench_dft.py --pdb LA_1264/system.pdb --ion LA --set 1264 --out bench_LA_32
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pdb", required=True)
    p.add_argument("--ion", required=True)
    p.add_argument("--charge", type=int, default=3)
    p.add_argument("--set", default="1264", help="metal parameter set used for the MM side")
    p.add_argument("--cutoff", type=float, default=3.3, help="first-shell cutoff (A)")
    p.add_argument("--functional", default="r2scan")
    p.add_argument("--basis", default="def2-svp")
    p.add_argument("--basis-file", help="NWChem basis file (e.g. for Lu, from basis_set_exchange)")
    p.add_argument("--ecp-file", help="NWChem ECP file")
    p.add_argument("--cores", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    p.add_argument("--repeat", type=int, default=3)
    p.add_argument("--platform", default="CPU", help="CPU, or GPU for gpu4pyscf (an NVIDIA node)")
    p.add_argument("--skip-xtb", action="store_true")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    from ash import Fragment, QMMMTheory, Singlepoint
    from mcube import ions
    from mcube.qmmm import make_dft, make_mm, make_qm
    from mcube.umbrella import qm_charge, umbrella_qm_region

    pdb = str(Path(args.pdb).resolve())
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    os.chdir(out)
    forcefields = ["amber14/tip3p.xml"]
    preprocess = {"inputs": {"metals": [args.ion], "cutoff": args.cutoff}}

    def system():
        mm = make_mm(forcefields, pdb, False, args.cores)
        ion = next(a.index for a in mm.topology.atoms() if a.residue.name == args.ion)
        qm_atoms, boundary, residues = umbrella_qm_region(pdb, {}, preprocess, ion)
        ions.apply_metal_parameters(mm.system, mm.topology, [args.ion], args.set,
                                    forcefields=forcefields, qm_atoms=qm_atoms)
        return mm, qm_atoms, boundary, residues

    mm, qm_atoms, boundary, residues = system()
    charge, raw = qm_charge(mm.system, qm_atoms)
    fragment = Fragment(pdbfile=pdb)
    dft_in = {"dft_functional": args.functional, "dft_basis": args.basis,
              "dft_densityfit": True, "dft_gridlevel": 3, "dft_platform": args.platform}
    if args.basis_file:
        dft_in["dft_basis_file"] = str(Path(args.basis_file).resolve())
    if args.ecp_file:
        dft_in["dft_ecp_file"] = str(Path(args.ecp_file).resolve())
    dft = make_dft(dft_in, args.cores, 1, "dft")
    theory = QMMMTheory(qm_theory=dft, mm_theory=mm, fragment=fragment, qmatoms=qm_atoms,
                        excludeboundaryatomlist=boundary, embedding="electrostatic", printlevel=1)
    report = {"ion": args.ion, "qm_atoms": len(qm_atoms), "qm_residues": residues,
              "charge": charge, "charge_from_mm": raw, "functional": args.functional,
              "basis": args.basis_file or args.basis, "cores": args.cores,
              "platform": args.platform, "dft_s": []}
    for _ in range(args.repeat):
        t0 = time.time()
        res = Singlepoint(theory=theory, fragment=fragment, charge=charge, mult=1, Grad=True)
        report["dft_s"].append(round(time.time() - t0, 1))
        report["dft_energy_Eh"] = float(res.energy)
        print("[BENCH] DFT/MM gradient", report["dft_s"][-1], "s", flush=True)

    # xTB/MM on a fresh MM system (QMMMTheory edits the forces it is given)
    if args.skip_xtb:
        steady = min(report["dft_s"][1:] or report["dft_s"])
        report["ps_per_day_at_1fs"] = round(86400.0 / steady / 1000.0, 2)
        Path("report.json").write_text(json.dumps(report, indent=1))
        print(json.dumps(report, indent=1))
        return
    mm, qm_atoms, boundary, _ = system()
    xtb = make_qm({"xtb_method": "GFN2"}, min(args.cores, 4))
    theory = QMMMTheory(qm_theory=xtb, mm_theory=mm, fragment=fragment, qmatoms=qm_atoms,
                        excludeboundaryatomlist=boundary, embedding="electrostatic", printlevel=1)
    t0 = time.time()
    Singlepoint(theory=theory, fragment=fragment, charge=charge, mult=1, Grad=True)
    report["xtb_s"] = round(time.time() - t0, 1)

    steady = min(report["dft_s"][1:] or report["dft_s"])
    report["ps_per_day_at_1fs"] = round(86400.0 / steady / 1000.0, 2)
    Path("report.json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
