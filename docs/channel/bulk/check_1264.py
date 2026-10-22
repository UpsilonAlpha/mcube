#!/usr/bin/env python3
"""
Does mcube's 12-6-4 give the same energies as Amber's own implementation?

The reference is built entirely with AmberTools: tleap writes a prmtop for the box with
the ion's 12-6-4 12-6 terms (frcmod.ions234lm_1264_tip3p), ParmEd's `add12_6_4` adds
the C4 matrix (LENNARD_JONES_CCOEF) from lj_1264_pol.dat, and OpenMM's AmberPrmtopFile
builds the System from that file. mcube's System is built the way the md and umbrella
stages build it: amber14/tip3p.xml plus `ions.apply_metal_parameters(..., "1264")`.

Both are evaluated on the same frames (no constraints, PME, same cutoff), and the
potential energies, the ion's force, and the energy of moving the ion are compared.
A match means a disagreement with published 12-6-4 hydration numbers is not an
implementation error.

  python check_1264.py --pdb LU_1264/system.pdb --ion LU --out check_LU
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np

TLEAP = """source leaprc.water.tip3p
loadamberparams frcmod.ions234lm_1264_tip3p
m = loadpdb {pdb}
set m box {{ {a:.4f} {b:.4f} {c:.4f} }}
saveamberparm m sys.prmtop sys.inpcrd
quit
"""


def amber_pdb(src: str, dst: str):
    """The box with residue names tleap knows (HOH -> WAT), one TER per residue."""
    out = []
    for line in open(src):
        if line.startswith(("ATOM", "HETATM")):
            line = line.replace(" HOH ", " WAT ")
            if line[17:20].strip() == "CL":
                line = line[:12] + " CL " + line[16:]     # OpenMM writes "Cl", leap wants "CL"
            out.append(line)
            if line[12:16].strip() in ("H2", "LU", "LA", "CA", "CL", "NA") or \
                    line[17:20].strip() in ("LU", "LA", "CA", "CL", "NA"):
                out.append("TER\n")
    Path(dst).write_text("".join(out) + "END\n")


def box_of(pdb: str):
    for line in open(pdb):
        if line.startswith("CRYST1"):
            return [float(line[6:15]), float(line[15:24]), float(line[24:33])]
    raise ValueError("no CRYST1")


def energy_and_ion_force(system, positions, box_nm, ion):
    import openmm
    from openmm import unit
    ctx = openmm.Context(system, openmm.VerletIntegrator(0.001),
                         openmm.Platform.getPlatformByName("Reference"))
    ctx.setPeriodicBoxVectors(*[openmm.Vec3(*v) for v in np.diag(box_nm)])
    ctx.setPositions(positions)
    st = ctx.getState(getEnergy=True, getForces=True)
    e = st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    f = st.getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole / unit.angstrom)[ion]
    return e, f


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pdb", required=True)
    p.add_argument("--ion", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    import openmm
    from openmm import unit
    from openmm.app import AmberPrmtopFile, ForceField, PDBFile, PME
    from mcube import ions

    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    pdb_src = Path(args.pdb).resolve()
    os.chdir(out)
    a, b, c = box_of(str(pdb_src))
    amber_pdb(str(pdb_src), "amber.pdb")
    Path("leap.in").write_text(TLEAP.format(pdb="amber.pdb", a=a, b=b, c=c))
    subprocess.run(["tleap", "-f", "leap.in"], check=True, stdout=open("tleap.log", "w"),
                   stderr=subprocess.STDOUT)
    pol = Path(os.environ["CONDA_PREFIX"]) / "dat" / "leap" / "parm" / "lj_1264_pol.dat"
    Path("parmed.in").write_text(f"add12_6_4 :{args.ion} watermodel TIP3P polfile {pol}\n"
                                 "outparm sys1264.prmtop\nquit\n")
    subprocess.run(["parmed", "-i", "parmed.in", "sys.prmtop"], check=True,
                   stdout=open("parmed.log", "w"), stderr=subprocess.STDOUT)

    settings = dict(nonbondedMethod=PME, nonbondedCutoff=1.0 * unit.nanometer,
                    constraints=None, rigidWater=False)
    box_nm = np.array([a, b, c]) / 10.0

    # Amber reference
    prmtop = AmberPrmtopFile("sys1264.prmtop")
    ref = prmtop.createSystem(**settings)
    # mcube
    pdb = PDBFile(str(pdb_src))
    mine = ForceField("amber14/tip3p.xml").createSystem(pdb.topology, **settings)
    changed = ions.apply_metal_parameters(mine, pdb.topology, [args.ion], "1264",
                                          forcefields=["amber14/tip3p.xml"])
    ion = changed[0]
    # tleap regroups residues, so match atoms by the coordinates it wrote (inpcrd, 1e-7 A)
    from openmm.app import AmberInpcrdFile
    ref_xyz = np.asarray(AmberInpcrdFile("sys.inpcrd").getPositions(asNumpy=True)
                         .value_in_unit(unit.angstrom))
    pos = pdb.getPositions(asNumpy=True)
    my_xyz = np.asarray(pos.value_in_unit(unit.angstrom))
    # tleap also translates the system, so match after removing the centroid shift
    ref_xyz = ref_xyz - ref_xyz.mean(axis=0) + my_xyz.mean(axis=0)
    from scipy.spatial import cKDTree
    dist, order = cKDTree(my_xyz).query(ref_xyz)   # ref atom k is mcube atom order[k]
    order = [int(i) for i in order]
    if sorted(order) != list(range(len(my_xyz))) or dist.max() > 0.01:
        raise SystemExit(f"could not match every prmtop atom to a PDB atom by position "
                         f"(max mismatch {dist.max():.3f} A)")
    ion_ref = order.index(ion)

    rows = []
    # Frame as is, then the ion displaced along x: the energy *change* isolates the ion's terms
    for shift in (0.0, 0.2, 0.4, -0.2):
        p2 = pos.copy()
        p2[ion] = p2[ion] + openmm.Vec3(shift / 10.0, 0, 0) * unit.nanometer
        p_ref = [p2[i] for i in order]
        e_ref, f_ref = energy_and_ion_force(ref, p_ref, box_nm, ion_ref)
        e_mine, f_mine = energy_and_ion_force(mine, p2, box_nm, ion)
        rows.append({"shift_A": shift, "E_amber": e_ref, "E_mcube": e_mine,
                     "F_ion_amber": [round(float(x), 4) for x in f_ref],
                     "F_ion_mcube": [round(float(x), 4) for x in f_mine]})
    e0a, e0m = rows[0]["E_amber"], rows[0]["E_mcube"]
    report = {"ion": args.ion, "atoms": len(order),
              "dE_total_kcal": round(e0m - e0a, 4),
              "dE_ion_moves": [{"shift_A": r["shift_A"],
                                "amber": round(r["E_amber"] - e0a, 4),
                                "mcube": round(r["E_mcube"] - e0m, 4)} for r in rows[1:]],
              "F_ion": [{"shift_A": r["shift_A"], "amber": r["F_ion_amber"],
                         "mcube": r["F_ion_mcube"]} for r in rows]}
    Path("report.json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
