# Ion-channel permeation PMF: umbrella sampling plan

How we go from a sequence to PMFs for Ca²⁺, La³⁺ and a heavy lanthanide crossing the channel.
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

## The `umbrella` stage (design, to be finalised)

One config section, `umbrella`, read and written like every other stage. Sketch:

```yaml
umbrella:
  inputs:
    ion: "B:1:LA"                       # selector for the permeating ion
    reference: {chain: A, atoms: CA}    # channel atoms whose COM defines z = 0
    start: -35.0                        # Å, CV range along the membrane normal
    end: 35.0
    spacing: 0.5
    k: 10.0                             # kcal/mol/Å², umbrella force constant
    radial: {r_max: 8.0, k: 10.0}       # flat-bottomed cylinder around the axis
    engine: mm                          # mm | qmmm
    equilibrate: 100                    # ps per window, coordination restrained
    time: 5000                          # ps per window, unrestrained sampling
  results: {}                           # windows, CV files, PMF
```

Open design questions: OpenMM force vs PLUMED (ASH check), how windows map to jobs
(`--window N`), and what the WHAM step writes back.

## Status

| step | state |
|---|---|
| decisions C1–C9 recorded | done 2026-10-22 |
| ASH biasing / PLUMED / GPU check | in progress |
| literature, lanthanide partner, 12-6-4 | in progress |
| Bunya probes Q1, Q2 | waiting for the VS Code session |
| fold | not started |
