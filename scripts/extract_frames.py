#!/usr/bin/env python3
"""
Write chosen trajectory frames as full-system PDBs for later stages
(to name in a section's `pdb_file`).

Frames are written with OpenMM and the topology PDB, so atom order, names,
residue ids, chains and ligand CONECT records match the input exactly, and
CRYST1 holds the frame's own box (the NPT-equilibrated one), not the
topology file's.

Usage:
    python extract_frames.py topology.pdb traj.dcd --frames 300 400 -1 [--outdir frames]
"""

import argparse
from pathlib import Path

import mdtraj
from openmm import unit
from openmm.app import PDBFile


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("topology")
    p.add_argument("trajectory")
    p.add_argument("--frames", type=int, nargs="+", required=True, help="0-based; negative counts from the end")
    p.add_argument("--outdir", default="frames")
    p.add_argument("--prefix", default="frame")
    args = p.parse_args()

    traj = mdtraj.load(args.trajectory, top=args.topology)
    omm_top = PDBFile(args.topology).topology
    if traj.n_atoms != omm_top.getNumAtoms():
        raise SystemExit(f"Trajectory has {traj.n_atoms} atoms, topology {omm_top.getNumAtoms()}")
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    for i in args.frames:
        k = i % traj.n_frames
        if traj.unitcell_vectors is not None:
            omm_top.setPeriodicBoxVectors(traj.unitcell_vectors[k] * unit.nanometer)
        name = out / f"{args.prefix}{k:05d}.pdb"
        with open(name, "w") as f:
            PDBFile.writeFile(omm_top, traj.xyz[k] * 10.0 * unit.angstrom, f, keepIds=True)
        box = traj.unitcell_lengths[k] * 10.0 if traj.unitcell_lengths is not None else None
        print(f"frame {k} (t = {traj.time[k]:.0f} ps) -> {name}; box {box} A")


if __name__ == "__main__":
    main()
