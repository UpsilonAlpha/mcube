# Ion-channel permeation PMF: umbrella sampling plan

How we go from a sequence to PMFs for Ca²⁺, La³⁺ and Lu³⁺ crossing the channel.
Decisions are in [channel-decision-log.md](channel-decision-log.md) (C*, F*, Q*).

## Working notes (read first in a new session)

- **State lives in files, not in the conversation.** Decisions → the decision log. Protocol and
  status → this file's Status table. Cluster facts → `running-on-bunya.md` and
  `binder-physics/docs/BUNYA_BROWSER_SHELL.md` (authoritative for the browser bridge).
- **Paths.** Local repo `Library/University/Year1/Structures/mcube`; on Bunya
  `/scratch/user/uquadhik/mcube` (pull, never type long payloads). Run directory on Bunya:
  `/scratch/user/uquadhik/channel/` (not yet created).
- **One agent drives the Bunya pane at a time.** Agents that only research or write code run in
  parallel and never touch the pane.
- **Delegate wide reads.** ASH source, literature and long logs go to subagents that return a
  summary; keep file dumps out of the main context.
- Every cluster script goes in `docs/channel/` and is pulled; checksum anything transferred.

## Pipeline

```
0 fold      AF2 (alphafold_custom), monomer + trimer, mature chain        C1, C2
1 orient    barrel axis on z, aromatic girdles at ±~12 Å                   
2 preprocess membrane: true, POPE, metals: [<ion>], ion in protein_pdb    C3, C8
3 md        membrane NPT equilibration (GPU), coordination restrained     C8
4 umbrella  pull → windows → per-window MM sampling (GPU)                 C5–C7
5 umbrella  QM/MM windows over the constriction (CPU xTB)                 C5, C9
6 sp        DFT single points on sampled QM/MM frames                     C5
7 analysis  WHAM/MBAR, MM→QM/MM and xTB→DFT reweighting                   C5
```

## System build per ion

1. Equilibrate the apo channel in POPE once (`preprocess` membrane + `md`).
2. Per ion, put the ion into a bulk water position of that box and re-neutralise, as in the
   README's "new partner in an equilibrated bilayer" route.
3. `preprocess` on that file with `skip_solvate: true`, `membrane: true`, `metals: [LA]`,
   `cutoff`: it detects the hydration shell and defines the QM region. (Checked in the code:
   with `skip_solvate`, `membrane` only adds the lipid force field, which is what is wanted.)
4. `md`, short: minimise, NPT and ~1 ns NVT with the shell restrained (C8). Its `lastframe` is
   where `pull` starts.

## The `umbrella` stage (implemented on branch `umbrella`)

`src/mcube/umbrella.py`. The full input table and the actions are in the README
(`### umbrella.inputs`); a complete config is `config_examples/umbrella_config.yaml`. In short:

```bash
mcube umbrella config.yaml pull                                   # one MM job (GPU)
mcube umbrella config.yaml window --window N                      # MM window N (GPU), one job each
mcube umbrella config.yaml window --window N --section umbrella_qmmm      # QM/MM window (CPU xTB)
mcube umbrella config.yaml energies --window N --level mm|xtb|dft [--section S]
mcube umbrella config.yaml wham [--section S] [--reweight dft]
mcube umbrella config.yaml stitch
```

What differs from the first design, and why (decision log): the bias is a CustomCVForce of two
centroid CVs with explicit element-mass weights (C10); MM windows run on OpenMM directly (C13);
the QM region has a fixed composition, `qm_waters` + `qm_residues` (C9); `pull_speed` in Å/ns
replaces `pull_time` (C14); per-window results use their own locked merge, not `save_result`,
which replaces a whole section's results. There is no MBAR: WHAM is self-contained numpy.

Order of operations for one ion:

1. `pull`, then every MM window (about 141 for ±35 Å at 0.5 Å).
2. `wham` on the MM leg; read the overlap warnings and add windows where they are flagged.
3. Choose `qm_residues` from the MM windows around the constriction (the residues that ever
   coordinate), then every `umbrella_qmmm` window.
4. `energies --level mm` and `--level xtb` on the MM and QM/MM windows at both ends of the QM/MM
   range; `wham --section umbrella_qmmm`; `stitch`. Read the cycle closure.
5. `energies --level dft --stride K` on the QM/MM windows, then `wham --section umbrella_qmmm
   --reweight dft`.

## Status

| step | state |
|---|---|
| decisions C1–C14, F1–F2 recorded | done 2026-10-22 |
| ASH biasing / PLUMED / GPU check | done from source (C10, C11); no HIP in the ASH env, MD on gpu_cuda |
| literature, lanthanide partner, 12-6-4 | done (C4, C5, C12); La³⁺ vs Lu³⁺ chosen |
| `umbrella.py`, 12-6-4, F1 fix | on `fix/pipeline-review` and on Bunya; 168 tests pass in the ASH env; no QM/MM window run yet (Q6) |
| Bunya probes Q1, Q2, Q5 | done 2026-10-22; +2 12-6-4 rows and the AmberTools polarisability table added |
| fold | queued on A100: 30410203 (monomer), 30410204 (trimer); MI355X cannot run AF2 |
| bulk ion check (C15) | La 1264, Lu IOD, Lu 1264 done (decision log); Ca IOD/1264, La IOD rerunning (30416201–03); xTB/MM SCC fixes on trial (30416204–09) |
