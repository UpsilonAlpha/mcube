#!/usr/bin/env python3
"""
Unified pipeline for protein/ligand preprocessing, MD, and QM/MM.
"""

from __future__ import annotations
import argparse
import importlib
import os


# Thread placement has to be in the environment before ASH, OpenMM or xTB start,
# so a job script must set it. But which policy is right depends on the stage and
# on the config, which the job script would have to parse. `mcube env` works it
# out and prints the exports for the job script to apply, so nobody has to
# remember the rule and no stage can be launched under the wrong one.
#
# Two things decide it. xTB is a separate process, so binding its threads one per
# physical core is safe and was worth 2-9x on SMT nodes. PySCF and OpenMM run
# in-process, and binding there pins the main thread so later BLAS pools inherit
# a single-core mask (measured: ~18% CPU per thread). Stages that run several xTB
# processes at once must not bind either, or every worker lands on the same cores.
def thread_policy(stage: str, inputs: dict, cpus: int) -> dict:
    """{"OMP_NUM_THREADS": n, "bind": bool, "why": str} for a stage."""
    dft_engine = str(inputs.get("qm_engine", "xtb")).lower() == "dft"
    half = max(1, cpus // 2)

    if stage == "preprocess":
        return {"OMP_NUM_THREADS": 1, "bind": False,
                "why": "file handling; nothing to parallelise"}
    if stage == "sp" or dft_engine:
        return {"OMP_NUM_THREADS": 1, "bind": False,
                "why": "PySCF runs in-process and sets its own thread count from "
                       "sp_cores/dft_cores; binding would pin its BLAS threads"}
    if stage == "md" or (stage == "umbrella" and str(inputs.get("engine", "mm")).lower() != "qmmm"):
        return {"OMP_NUM_THREADS": half, "bind": False,
                "why": "OpenMM runs in-process; binding would pin it"}
    if stage == "neb":
        return {"OMP_NUM_THREADS": half, "bind": False,
                "why": "NEB images run as parallel xTB processes; binding would "
                       "put every image on the same cores"}
    if stage == "ts" and int(inputs.get("freq_workers", 1) or 1) > 1:
        return {"OMP_NUM_THREADS": half, "bind": False,
                "why": f"freq_workers={inputs.get('freq_workers')} runs parallel xTB "
                       "processes; binding would put every worker on the same cores"}
    return {"OMP_NUM_THREADS": half, "bind": True,
            "why": "one xTB process at a time, so bind one thread per physical core"}


# What each stage can actually use, from seff on real jobs (docs/running-on-bunya.md).
# Asking for more is silently wasted, so `mcube env` says so at submit time.
USEFUL_CPUS = {"preprocess": 2, "md": 16, "qmmm": 8, "scan": 8, "ts": 8, "sp": 32, "umbrella": 8}


def env_exports(config: str, stage: str, section: str | None = None) -> str:
    import yaml

    with open(config) as fh:
        data = yaml.safe_load(fh) or {}
    base = stage if stage not in ("sp",) else "qmmm"
    inputs = dict((data.get(base) or {}).get("inputs") or {})
    for name in (stage, section):
        if name and name != base:
            inputs.update((data.get(name) or {}).get("inputs") or {})

    cpus = int(os.environ.get("SLURM_CPUS_PER_TASK") or os.cpu_count() or 1)
    policy = thread_policy(stage, inputs, cpus)

    lines = [f"export OMP_NUM_THREADS={policy['OMP_NUM_THREADS']}"]
    if policy["bind"]:
        lines += ["export OMP_PROC_BIND=close", "export OMP_PLACES=cores"]
    else:
        lines.append("unset OMP_PROC_BIND OMP_PLACES")

    useful = USEFUL_CPUS.get(stage)
    notes = [f"[INFO] {stage}: OMP_NUM_THREADS={policy['OMP_NUM_THREADS']}, "
             f"binding {'on' if policy['bind'] else 'off'} - {policy['why']}"]
    if useful and cpus > useful * 1.5:
        notes.append(f"[WARNING] {cpus} CPUs allocated but {stage} uses about "
                     f"{useful}; the rest is wasted. See the resources table in "
                     "docs/running-on-bunya.md.")
    lines += [f"echo {note!r} >&2" for note in notes]
    return "\n".join(lines) + "\n"


def build_parser():
    parser = argparse.ArgumentParser(prog="mcube")

    subparsers = parser.add_subparsers(dest="command", required=True)

    p1 = subparsers.add_parser("preprocess", help="Prepare system")
    p1.add_argument("config", help="Input yaml file")

    p2 = subparsers.add_parser("md", help="Run molecular dynamics")
    p2.add_argument("config", help="Input yaml file")
    p2.add_argument("--section", help="Config section for results; its inputs override md.inputs")

    p3 = subparsers.add_parser("qmmm", help="Run QM/MM")
    p3.add_argument("config", help="Input yaml file")
    p3.add_argument("--pdb", help="Starting structure (overrides qmmm.inputs.pdb_file)")
    p3.add_argument("--section", default="qmmm", help="Config section to write results to")
    p3.add_argument("--actatoms", help="File of active-atom indices (e.g. another run's active_atoms)")

    # The single point is a separate stage from qmmm because a DFT SCF scales
    # across cores while the xTB MD in front of it does not: see sp.py.
    p9 = subparsers.add_parser("sp", help="QM/MM single point (PySCF) on a finished geometry")
    p9.add_argument("config", help="Input yaml file")
    p9.add_argument("--pdb", help="Structure to score (e.g. a stage's QM_MM_optimized.pdb)")
    p9.add_argument("--section", default="sp", help="Config section to write results to")

    # Not a stage: prints the thread settings a job script should apply for a
    # stage, so the job script does not have to know the rule. See thread_policy.
    p10 = subparsers.add_parser("env", help="Shell exports for a stage's thread placement")
    p10.add_argument("config", help="Input yaml file")
    p10.add_argument("stage", help="Stage the job is about to run")
    p10.add_argument("--section", help="Config section, if the stage takes one")

    # =========================================================
    # NEB
    # =========================================================
    p4 = subparsers.add_parser("neb", help="Nudged elastic band simulation")

    p4.add_argument("config", help="Input yaml file")
    p4.add_argument("--product", required=True, help="Product PDB file")
    p4.add_argument("--reactant", required=True, help="Reactant PDB file")
    p4.add_argument("--actatoms", help="File of active-atom indices (e.g. the reactant's active_atoms)")
    p4.add_argument("--section", default="neb",
                    help="Config section for results; its inputs override neb.inputs")
    p4.add_argument("--cores", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))

    # =========================================================
    # RELAXED SCAN / TS REFINEMENT + IRC
    # =========================================================
    p6 = subparsers.add_parser("scan", help="Relaxed QM/MM distance scan")
    p6.add_argument("config", help="Input yaml file")
    p6.add_argument("--pdb", help="Starting structure (e.g. the reactant's QM_MM_optimized.pdb)")
    p6.add_argument("--section", default="scan",
                    help="Config section defining the scan (its inputs override scan.inputs) and receiving results")
    p6.add_argument("--actatoms", help="File of active-atom indices (e.g. the reactant's active_atoms)")

    p7 = subparsers.add_parser("ts", help="TS optimisation, frequency check and IRC")
    p7.add_argument("config", help="Input yaml file")
    p7.add_argument("--pdb", help="TS guess (e.g. the NEB's saddlepoint.pdb)")
    p7.add_argument("--section", default="ts",
                    help="Config section for results; its inputs override ts.inputs")
    p7.add_argument("--actatoms", help="File of active-atom indices (e.g. the reactant's active_atoms)")

    # =========================================================
    # UMBRELLA SAMPLING
    # =========================================================
    p8 = subparsers.add_parser("umbrella", help="Umbrella sampling of an ion along z, and its PMF")
    p8.add_argument("config", help="Input yaml file")
    p8.add_argument("action", choices=["pull", "window", "energies", "wham", "stitch"])
    p8.add_argument("--window", type=int, help="Window index (window, energies)")
    p8.add_argument("--section", default="umbrella",
                    help="Config section; its inputs override umbrella.inputs (e.g. umbrella_qmmm)")
    p8.add_argument("--level", default="mm", choices=["mm", "xtb", "dft"],
                    help="energies: the model to score saved frames with")
    p8.add_argument("--stride", type=int, default=1, help="energies: every Nth saved frame")
    p8.add_argument("--qm-section", default="umbrella_qmmm",
                    help="Section whose QM-region rule energies (xtb/dft) and stitch use")
    p8.add_argument("--reweight", choices=["dft"],
                    help="wham: reweight each frame from xTB/MM to this level")
    p8.add_argument("--bootstrap", type=int, default=50, help="wham: block-bootstrap resamples")

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
    # `env` only reads the config: it must not import a stage, because the job
    # script calls it before deciding how to run one.
    if args.command == "env":
        print(env_exports(args.config, args.stage, getattr(args, "section", None)))
        return

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