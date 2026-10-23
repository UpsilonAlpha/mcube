# Ion-channel permeation PMF: umbrella sampling plan

How we go from a sequence to PMFs for Ca²⁺, La³⁺ and Lu³⁺ crossing the channel.
Decisions are in [channel-decision-log.md](channel-decision-log.md) (C*, F*, Q*).
What the protein is (identity, gene neighbourhood, AF2 model, conservation): [channel/protein.md](channel/protein.md).

## Working notes (read first in a new session)

- **State lives in files, not in the conversation.** Decisions → the decision log. Protocol and
  status → this file's Status table. Cluster facts → `running-on-bunya.md` and
  `binder-physics/docs/BUNYA_BROWSER_SHELL.md` (authoritative for the browser bridge).
- **Paths.** Local repo `Library/University/Year1/Structures/mcube`; on Bunya
  `/scratch/user/uquadhik/mcube` (pull, never type long payloads). Run directory on Bunya:
  `/scratch/user/uquadhik/channel/` (`af2/`, `bulk/`, `dft/`, `gxtb/`, `coord/`, `system/{apo,CA,LA,LU}/`).
- **Each Claude session has its own browser pane.** After the pane is reopened the old
  code-server terminal may silently run nothing: probe with `date > $S/.cc/s.log` and, if the
  file does not change, open a new one (F1 → "Terminal: Create New Terminal").
- **One agent drives the Bunya pane at a time.** Agents that only research or write code run in
  parallel and never touch the pane.
- **Delegate wide reads.** ASH source, literature and long logs go to subagents that return a
  summary; keep file dumps out of the main context.
- Every cluster script goes in `docs/channel/` and is pulled; checksum anything transferred.

## Pipeline

```
0 fold      AF2 monomer, mature chain (done)                              C1, C2
1 orient    crop 110-130, barrel on z, charge-free belt at 0 (done)       C18
2 preprocess apo in POPE (done); per ion, ion in bulk + skip_solvate      C3, C8, C18
3 md        membrane NPT, then short md per ion with the shell restrained  C8
4 umbrella  pull → MM 12-6-4 windows over the whole path (GPU)            C5–C7, C17
5 umbrella  r2SCAN DFT/MM windows over the constriction only              C16, C17
6 analysis  WHAM; MM↔DFT/MM stitching (BAR) at the constriction ends      C5
```
xTB/MM (and the xTB→DFT reweighting) is no longer part of the plan: F7, F8, C17.

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
4. `energies --level mm` and `--level dft` on the MM and DFT/MM windows at both ends of the QM/MM
   range; `wham --section umbrella_qmmm`; `stitch`. Read the cycle closure.
5. (Superseded) `energies --level xtb` and `--reweight dft` belonged to the xTB/MM leg, which
   is out (F7); with `umbrella_qmmm` run at DFT/MM (`qm_engine: dft`) there is nothing to reweight.

## Status

| step | state |
|---|---|
| `umbrella.py`, 12-6-4, F1 fix | done, on `fix/pipeline-review` and Bunya; 168 tests pass in the ASH env; no QM/MM window run yet (Q6) |
| fold, protein identity | done: AF2 monomer (trimer not supported); BLAST, neighbourhood, conservation in `channel/protein.md` |
| bulk ion check (C15, F4) | Ca 1264, La 1264, La IOD, Lu IOD, Lu 1264 done; all done; Lu CN 7.6 vs paper 9.0 confirmed in sander (2 ns, CN 7.64); OpenMM-on-prmtop (30430921) running |
| QM level (F7, F8, C17) | xTB/MM out; g-xTB set aside; r2SCAN DFT/MM at the constriction; 12-6-4 vs r2SCAN carboxylate check running (30454349); GPU r2SCAN benchmark failed on missing CUDA 12 libraries (F11), to fix |
| system (C18) | apo POPE done (30455128); Ca/La/Lu placed; after F9/F10 preprocess+md+pull on gpu_viz L40s (F10, F11 fixed; md done; pulls 30523774 La on viz, 30523775 Ca and 30523776 Lu on short) |
| MM umbrella test | windows −35…+35 Å, 1 Å, 2 ns: after the pulls (`umbrella_windows.slurm`) |
