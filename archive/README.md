# Archive

Code that worked but is not part of the pipeline any more. It is kept here
rather than deleted because it is referenced in the project's decision log and
because reviving it is a matter of moving files back, not rewriting them.

Nothing in this directory is imported by `mcube`, and the `tests/` suite does
not collect it.

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
