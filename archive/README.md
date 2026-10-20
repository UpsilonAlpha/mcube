# Archive

Code that worked but is not part of the pipeline any more. It is kept here
rather than deleted because it is referenced in the project's decision log and
because reviving it is a matter of moving files back, not rewriting them.

Nothing in this directory is imported by `mcube`, and the `tests/` suite does
not collect it.

## `threads.py` (archived 2026-10-20)

Thread placement, moved into the job script. It set `OMP_PROC_BIND=close` and
`OMP_PLACES=cores` for xTB alone, through a wrapper script written to
`.mcube_xtb/`, and left the Python process unbound so PySCF's BLAS pools were
not pinned to one core.

The measurements behind it still stand, and are now in the README's three
thread recipes and in `docs/csr/run.slurm`'s named profiles:

- on SMT nodes `--cpus-per-task=N` gives N/2 physical cores, and unbound xTB
  threads shared cores and ran **2–9× slower** with erratic step times;
- binding an in-process OpenMP runtime pins the main thread and later
  numpy/BLAS pools inherit that single-core mask (**~18% CPU per thread** for a
  DFT single point);
- so several xTB processes at once (NEB images, parallel NumFreq) must run
  **unbound**, or each binds to the same first cores.

What it was working around is that these variables have to be set before ASH,
OpenMM or xTB start, which a library cannot reliably do for a process it is
already inside. A job script can, and the policy is per-stage anyway, so that
is where it belongs. `resolve_cores` — the one piece stages still need — moved
to `preprocess.py` and now only resolves a count: the stage's `cores`, else
`$OMP_NUM_THREADS`, else `$SLURM_CPUS_PER_TASK`, else 1. `test_threads.py` is
here too; its `resolve_cores` cases live on in `tests/test_cores.py`.

## τRAMD (archived 2026-10-20)

Random acceleration MD for residence times (Lüdemann 2000; τRAMD protocol,
Kokh 2018/2020). Implemented, tested and documented, but it answers a binding
kinetics question, and the project it was built in is a mechanism study — the
`preprocess → md → qmmm → sp` path, with `scan`, `neb` and `ts` for the
reaction coordinate. It was carrying documentation and a CLI stage for nobody.

| File | What it is |
| --- | --- |
| `ramd.py` | the stage: replica setup, equilibration, trajectory loop |
| `ramd_engine.py` | the force: random direction, COM displacement check, redraw |
| `ramd_analysis.py` | τ per replica, bootstrap spread, KS against the Poisson CDF |
| `test_ramd.py` | unit tests for the engine and the analysis |
| `ramd.md` | the documentation |
| `protein_protein_ramd.yaml` | a worked config: protein–protein complex, no metals |
| `ramd_gpu.slurm` | GPU submission script |

To revive it: move `ramd*.py` back to `src/mcube/`, `test_ramd.py` to `tests/`,
and restore the `ramd` subparser in `cli.py` (it took `--section`, `--replicas`
and `--analyse`). It builds its MM system with `build_md_theory` from `md.py`,
which is still there and unchanged.
