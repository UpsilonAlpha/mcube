#!/usr/bin/env python3
"""
Unified pipeline for protein/ligand preprocessing, MD, and QM/MM.
"""

from __future__ import annotations
import argparse
import importlib
import os


def build_parser():
    parser = argparse.ArgumentParser(prog="mcube")

    subparsers = parser.add_subparsers(dest="command", required=True)

    p1 = subparsers.add_parser("preprocess", help="Prepare system")
    p1.add_argument("config", help="Input yaml file")

    p2 = subparsers.add_parser("md", help="Run molecular dynamics")
    p2.add_argument("config", help="Input yaml file")

    p3 = subparsers.add_parser("qmmm", help="Run QM/MM")
    p3.add_argument("config", help="Input yaml file")

    # =========================================================
    # NEB
    # =========================================================
    p4 = subparsers.add_parser("neb", help="Nudged elastic band simulation")

    p4.add_argument("config", help="Input yaml file")
    p4.add_argument("--product", required=True, help="Product PDB file")
    p4.add_argument("--reactant", required=True, help="Reactant PDB file")
    p4.add_argument("--cores", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))


    # =========================================================
    # POSTPROCESS
    # =========================================================
    p5 = subparsers.add_parser("postprocess", help="Analyse trajectory")

    p5.add_argument("config", help="Input yaml file")
    p5.add_argument("--trajectory", help="DCD file")
    p5.add_argument("--cutoff", type=float, default=2.6)
    p5.add_argument("--cores", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    p5.add_argument("--metals", nargs="+", default=["ZN"])
    p5.add_argument("--spring_constant", type=int, default=5)
    p5.add_argument("--iterations", type=int, default=200)
    p5.add_argument("--occupancy", type=float, default=0.1)
    p5.add_argument("--redraw", action="store_true")
    p5.add_argument("--qmmm", action="store_true")


    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    # =========================================================
    # DISPATCH
    # =========================================================
    # Stage modules are imported lazily so that e.g. `preprocess` does not
    # require pyscf or MDAnalysis to be importable.
    stage = importlib.import_module(f".{args.command}", package=__package__)
    result = stage.run(args)

    # =========================================================
    # OUTPUT
    # =========================================================
    if isinstance(result, dict):
        print("\n[RESULT]")
        for k, v in result.items():
            print(f"{k}: {v}")


if __name__ == "__main__":
    main()