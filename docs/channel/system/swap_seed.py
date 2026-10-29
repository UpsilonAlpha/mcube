#!/usr/bin/env python3
"""
A run directory for one MM umbrella window of a different ion, seeded from a frame of an
existing window: the bound La (or Lu) frame with the ion swapped, so the new ion starts in the
same site and the window samples its own equilibrium there. Made for the Ca endpoint of the QM/MM
exchange energies (decision log C22): no Ca run ever reached the Glu ring, so the Ca->Ln swaps were
one-sided, all on frames arranged for a 3+ ion.

    python swap_seed.py --window $S/channel/system/LA/umbrella_res_035 --frame 159 --ion CA \
        --base $S/channel/system/CA/config.yaml --out $S/channel/system/CA_la \
        --set k=2.5 time=10000 frames_at=[0.0] traj_frequency=6250

Then, from --out: sbatch -J res_CA_la --export=ALL,FIRST=35,LAST=35,SECTION=umbrella_res
umbrella_windows.slurm. The window restrains the seed's first shell during equilibration (as every
window does) and then releases it, so the ion may keep or leave the carboxylates.

- The frame is read whole about the protein (exchange.py's imaging) from the window's DCD, on
  start.pdb's topology (which carries the bonds).
- A 3+ -> 2+ swap leaves the box at -1 e: the Na+ farthest from the ion goes, so the box stays
  neutral under PME. A 2+ -> 3+ swap would need a Cl- removed the same way.
- config.yaml = --base's preprocess, md and umbrella inputs (the new ion's metal parameters),
  with preprocess.processed_pdb pointing at the seed (it carries every bond), and one section
  (--section, default umbrella_res) whose window --window-index has `seed` set.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "qmmm"))
import exchange  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--window", required=True, help="source window dir (start.pdb, window.dcd)")
    p.add_argument("--frame", type=int, required=True, help="DCD frame index")
    p.add_argument("--ion", required=True, choices=sorted(exchange.IONS))
    p.add_argument("--base", required=True, help="config.yaml of a run of --ion (its inputs)")
    p.add_argument("--out", required=True)
    p.add_argument("--section", default="umbrella_res")
    p.add_argument("--window-index", type=int, default=35, help="window whose seed this is (z0 0)")
    p.add_argument("--set", nargs="*", default=[], help="KEY=VALUE section inputs (YAML values)")
    args = p.parse_args()

    from openmm import unit
    from openmm.app import Modeller, PDBFile

    src = Path(args.window)
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    topology = exchange.load_topology(src / "start.pdb")
    mdtop, anchors = exchange.imaging(topology, "protein")
    xyz, box_nm = exchange.read_frame(src / "window.dcd", args.frame, mdtop, anchors)
    ion = exchange.find_ion(topology)
    old = ion.residue.name
    dq = exchange.IONS[old][1] - exchange.IONS[args.ion][1]
    exchange.swap(topology, ion, args.ion)

    mod = Modeller(topology, xyz / 10.0 * unit.nanometer)
    if dq:
        counter = "NA" if dq > 0 else "CL"
        cands = [a for a in mod.topology.atoms() if a.residue.name == counter]
        d = [float(np.linalg.norm(xyz[a.index] - xyz[ion.index])) for a in cands]
        drop = [cands[i] for i in np.argsort(d)[::-1][:abs(dq)]]
        print(f"[INFO] {old} -> {args.ion}: removing {len(drop)} {counter} at "
              + ", ".join(f"{d[cands.index(a)]:.1f} A" for a in drop) + " from the ion")
        mod.delete([a.residue for a in drop])
    from openmm import Vec3
    mod.topology.setPeriodicBoxVectors([Vec3(*map(float, v)) for v in box_nm] * unit.nanometer)
    seed = out / f"seed_{old}{args.frame:04d}.pdb"
    with open(seed, "w") as fh:
        PDBFile.writeFile(mod.topology, mod.positions, fh, keepIds=True)
    new_ion = exchange.find_ion(PDBFile(str(seed)).topology)
    ion_id = f"{new_ion.residue.chain.id}:{int(new_ion.residue.id)}:{args.ion}"

    base = yaml.safe_load(Path(args.base).read_text())
    pre = base["preprocess"]
    extra = {k: yaml.safe_load(v) for k, v in (s.partition("=")[::2] for s in args.set)}
    config = {
        "preprocess": {"inputs": pre["inputs"],
                       "results": {**(pre.get("results") or {}), "processed_pdb": seed.as_posix()}},
        "md": {"inputs": base["md"]["inputs"]},
        "umbrella": {"inputs": base["umbrella"]["inputs"]},
        args.section: {"inputs": {**extra, "ion": ion_id},
                       "results": {"windows": {f"{args.window_index:03d}": {"seed": seed.as_posix()}}}},
    }
    (out / "config.yaml").write_text(
        f"# swap_seed.py: {old} frame {args.frame} of {src} with the ion swapped to {args.ion}\n"
        + yaml.safe_dump(config, sort_keys=False))
    shell = exchange.shell(PDBFile(str(seed)).topology, np.array(PDBFile(str(seed)).positions.value_in_unit(unit.angstrom)),
                           np.diag(box_nm) * 10.0, new_ion, 3.0, {"O", "N", "S"})
    print(f"[INFO] {seed}: ion {ion_id}; shell within 3.0 A: "
          + ", ".join(f"{r.name}{int(r.id)}:{a.name} {d:.2f}" for r, a, d in shell))
    print(f"[INFO] {out / 'config.yaml'}: section {args.section}, window {args.window_index:03d}")


if __name__ == "__main__":
    main()
