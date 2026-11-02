# CSR-1: current state and plan

**This file is the single source of truth for what is true now.** The decision
log (`../csr-decision-log.md`) is append-only history: it records how each
conclusion was reached, including the ones later overturned, and it is 180
entries long. Reading it end to end to find the current value of something is
how the same facts get rediscovered, and how a contradiction can sit unnoticed
for a week (see *Contradictions caught late* below).

**The rule.** When something changes, the log gets a new entry **and this file
is edited in place**. A log entry that changes a settled fact or a plan item
without editing this file is incomplete. Status words here are
`settled` / `measured, open` / `void` / `not started` / `in flight`, and
nothing else.

Last updated: 2026-11-02.

---

## 1. Settled facts

These are not to be re-derived. Each gives the entry that settled it and, where
one exists, the entry that wrongly contradicted it later.

### The enzyme and the site

| Fact | Status | Entry |
| --- | --- | --- |
| CSR-1 is a B3-subclass MBL from *Cronobacter sakazakii*, UniProt A7MFY4, published as characterisation only (Pedroso 2020) | settled | D22 |
| **No published kinetics for CSR-1 with meropenem, and no published mechanistic study.** The user has measured kinetics and is **withholding them** so the computation stands alone | settled | D22, D40, **D178** |
| The Zn1 site is **H-R-H**: Arg90 (BBL R118) where canonical B3 has a third His | settled | D22, user correction ×2 |
| The site is **mono-zinc**; no Zn2 to stabilise the anionic leaving nitrogen | settled | D22, D107 |
| **Zn2 is vacant in the wild type, by positive evidence:** 6NC5's `SEQADV` shows HIS A92←GLN and HIS A225←LYS as engineered mutations restoring a canonical Zn2 ligand set. The only way anyone put a metal in Zn2 was to engineer two of its ligands back | settled | D22/D40, audit §3(a) |
| The modelled zinc is **crystallographic**: 6NC5's ZN302, whose own ligands His87/His163 are *not* mutated. The wild-type crystals are metal-free because of affinity and crystallisation, not an absent site | settled | D40 correcting D22 |
| **6NC5 is a Zn2-site mutant and is NOT a wild-type reference for Asp92. 6DN4 is.** | settled | **D108** — then contradicted by D145, re-settled by D179 |
| BBL mapping: H116→His88, R118→Arg90, D120→Asp92, Q121→Gln93, H196→His164, K254/263→Lys226. 6NC5 is our numbering **−1**; 6DN4/6DQ2 are in **BBL numbering** | settled | D108, D179 |

### The crystallographic zinc sphere

| Measurement | Value | Entry |
| --- | --- | --- |
| 6NC5 ZN302, four ligands within 2.8 Å | HOH485 1.90, His163 NE2 2.09, HOH493 2.14, His87 NE2 2.22 | D108 |
| **No aspartate on ZN302** | Asp91 OD2 4.82 / OD1 4.13 Å from His87 NE2 — *in the mutant* | D108, D179 |
| Asp92's conformer across CSR-1's own three crystals | **6DN4 (WT, 1.99 Å, B 15–16): 3.00/3.07 Å from His116 NE2 — close.** 6DQ2 (mutant S60, B 53–66): 4.28/5.49. 6NC5 (Zn2 mutant): 4.13/4.82 | **D179** |
| **Zn–Asp92 OD2 in the wild type: 2.68 Å** — 6NC5's ZN302 superposed into 6DN4 on the two histidine rings alone, fit RMSD 0.184 Å. Our model has 2.80, and all six site distances agree to within 0.13 Å | settled — **the starting structure is right**; D145's 4.20 Å was the Zn2 mutant's | **D184** |
| Gln93 is remote, from the crystal as well as the trajectory | NE2 **7.10 Å** from the zinc in 6DN4; 7.95–9.88 Å from the ligand over 34 ns | D184 + D173 |
| Therefore: the carboxylate is **mobile**, and the only well-ordered wild-type copy has it close | measured, open | D179 + D151 |
| L1 7UHI/7UHJ, the only measured B3 Michaelis complex | intact carbonyl **2.53/2.63 Å** from Zn1; ring N and C-carboxylate on Zn2 | D156 |
| CSR-1 vs L1 similarity | 39.0% identity, 221/249 core, **1.05 Å CA RMSD** — closer than L1 vs SMB-1 | D158 |

### Pipeline rules that are now enforced, not merely intended

| Rule | Status | Entry |
| --- | --- | --- |
| Coordination restraints are detected **per atom**, not per residue | settled (fixed 2026-11-02) | D175, D176 |
| Detection cutoff is **3.0 Å**; everything it finds is restrained at the distance the structure has | settled (user's rule) | D174 |
| Restraints are stored as **selectors**, never atom indices | settled | D115, D116 |
| A chosen restraint in one run may not become a detected one in the next (**laundering refusal**) | settled | D68, D71 |
| `rigid_water: false` by default for gradient stages | settled | D137–D142 |
| OpenFF **renames ligand atoms**; the map is recorded and warned about | settled | D168, D169 |
| Coordination occupancy is measured **over the whole trajectory**, not from `md`'s two frames | settled (stage added 2026-11-01) | D167, D177 |
| Hydrogens are checked on every built system, and the ligand's formula against its SMILES | settled | D172 |

---

## 2. Void — do not reuse these results

| What | Why | Entry |
| --- | --- | --- |
| Everything D28–D70 | **restraint laundering**: run 2 froze run 1's imposed 1.94 Å Zn–carbonyl holonomically. Six constraints where there were four | D67, D68 |
| `scan_H1_f00497`, `scan_H2_f00395`, D81's pair, D105's endpoint | frozen-QM-water bug and/or `$set` rigid-water convergence; **flagged, never re-measured** | D36, D137–D142 |
| D143's "the crystal's 2.80 Å sits 5.4–7.9 kcal/mol up this curve" | rested on 2.80 Å being crystallographic | D145 |
| **D145's own 4.20 Å "crystal" value** | it is 6NC5's — the Zn2 mutant — which D108 had already excluded as a WT reference | **D179** |
| `md_wt100`, `md_wt100b` and everything descended from them | **zinc six-coordinate in 4995/5000 frames**, Asp92 OD2 bound at 2.07 Å in 100% | **D177** |
| `md_prod100` (34 ns before the NaN) | site rearranged: Asp92 to 2.03 Å, ligand turned over, O2x travelled 4.75 Å onto the metal | D173 |
| `prod100b` | withdrawn: it restarted from a frame of the above | D173 |
| Any mechanistic claim resting on the docked pose | no Vina pose puts an atom within 5 Å of the zinc; pose 2 is **7.23 Å**, and the structure being run is **5.15 Å** from L1's measured pose | D126, D156 |

**No reaction path has ever been computed.** `neb` has never completed, `ts`
has never run; there is no barrier, no transition state, no IRC, no
intermediate. (Audit headline, confirmed against D144.)

---

## 3. Contradictions caught late — the failure mode this file exists to stop

| The contradiction | Days unnoticed | How it was caught |
| --- | --- | --- |
| D145 quoted 6NC5's Asp distance as "the crystal's" when **D108 had already ruled 6NC5 out as a WT reference** | 8 | measurement, D179 |
| D120 prescribed the coordination fix; it never worked | 5 | occupancy over a trajectory, D177 |
| D48 settled the resting state as Zn-OH₂; the whole `CSR_nac` chain ran in the other state | 8 | D112 |
| `coordination_from_rows` said restraints are per atom; `find_coordinators` called it per residue | ~6 weeks | D175 |
| A design decision (Gln93 out of the QM region) inferred from config files instead of read from the log | — | user, twice: D165, D171 |

Every one of these was a *reading* failure, not a computing failure.

---

## 4. The plan (audit `audit-2026-11-02.md`, adopted in D178)

Items 1–3 gate everything below them. **Status is maintained here.**

| # | Item | Status |
| --- | --- | --- |
| 0 | Asp92's placement in the starting structure | **settled (D179)** — 2.80 Å supported; `prod100c` unblocked |
| 1+2 | Experiment-anchored Michaelis complex, both endpoints on one topology | **in flight**: `prod100c` (31423830) queued → `prodqm` → `prodopt` → `recombine`. 7UHI is the *check*, not the source |
| 3 | Equilibration that preserves the site | **RUNNING on both arms** on bun076: `prod100c` (31423830) at 13.6/100 ns, **7 detected restraints + 1 configured**; `wt100c` (31483853) at 7.6/100 ns, **5 + 1**. Both lists verified correct, both fit in 5 h (D185) |
| 4 | Nucleophile identity, origin, protonation state | not started — blocked on 1–3. D48 and D107 disagree **in sign** |
| 5 | PROPKA/H++ and the His164 HID/HIE sensitivity check | **blocked on one pip install**: neither `propka` nor `pdb2pqr` is in the ASH env, both are pure Python, so `pip install propka pdb2pqr` is the whole of it — the user's env, so their call (D185). Open since **D7** |
| 6 | scan → NEB → TS → frequencies → IRC, ≥3 snapshots | not started — the paper's actual content |
| 7 | Cluster-model DFT-vs-GFN2 benchmark, Imhof-style | **not started; depends on nothing.** Open since **D9** |
| 8 | Protonation of the leaving nitrogen | not started |
| 9 | Per-residue roles reaching R, Q and K | **design decision due before item 3's equilibration** |
| 10 | Di-zinc comparison at identical level | not started — scope decision |
| 11 | Asp92 free energy | deferred by design — only if the discrepancy survives 1–3 |
| 12 | `md` resume | **done (D182)** — one implementation now serves both stages |

**Item 9's deadline is the non-obvious one.** Arg90 and Lys226 are inside the
QM region in all eleven systems, so charge deletion can never reach them
(D153); only Gln93 (−0.87) and Arg229 (+44.34) are measurable. Whether the
production system puts R, K, Arg229 and Q **all in MM** so the set can be
decomposed consistently, or whether the roles are done by computational
mutagenesis, changes what gets built — so it must be decided *before* the
equilibration.

### The blind-prediction protocol (D178)

The kinetics exist and are withheld. Therefore: the computed barrier, its level
of theory, its spread over snapshots, and the mechanism's qualitative claims go
into the log **before** the numbers are released, and are **not revised
afterwards**. A barrier that moves once the experiment arrives is a fit.
**Do not ask for the kinetics.**

---

## 4a. Two radii, and the gate before any QM/MM (D186)

**They answer different questions and must not share a number.**

| radius | key | question |
| --- | --- | --- |
| **3.0 Å** | `preprocess.inputs.cutoff` | what is close enough to be **affected** by the metal, and so must be held at its crystal distance. Wider than the sphere on purpose: 2.5–3.0 Å is where a non-bond becomes a bond if nothing holds it |
| **2.5 Å** | `postprocess.coordination.cutoff` | what is **bonded** — the coordination number. `LIGAND_MAX`, the longest Zn–O that is still a bond |

At 2.5 Å the site is **four**: two histidines, the hydroxide-derived oxygen,
one water. Asp92 OD2 at 2.80 and the product ligand's O4x/O5x at 2.92/2.94 are
outside it — restrained second-shell contacts, not ligands.

`postprocess.coordination` with no `cutoff` falls back to
`preprocess.inputs.cutoff`, which is now the restraint radius — so **always set
`coordination.cutoff: 2.5` in an analysis section.** `prodpostc` and `wtpostc`
do. `margin: 1.0` still shows everything out to 3.5 Å, which is how the held
contacts stay visible without being counted.

**The gate:**

1. **Four-coordinate at 2.5 Å**, at high occupancy.
2. **Asp92 OD2 near 2.80, not collapsed toward 2.0** — the artefact of D173/D177,
   now restrained, so this tests whether 50 kcal/mol/Å² holds it.
3. **Asp92 OD1 near 3.81**, not pivoted in (D120, D130).
4. **Product only: O4x/O5x near 2.9, O2x away** — its 4.75 Å migration onto the
   metal was the turnover of D173.
5. Hydrogens clean, `missing` list empty.

**Not a gate: a transient bulk water.** MM cannot be trusted for coordination —
that is why the restraints exist — and water comes and goes. What the MM
occupancy tests is whether the *restrained* contacts held. The coordination
behaviour that counts is QM/MM's.

## 5. Open questions, with what would settle each

| Question | What would settle it | Blocked by |
| --- | --- | --- |
| Is the nucleophile hydroxide or water, and how does it form with no Zn2? | equilibrate each state, keep the one that holds the crystal's sphere (D113/D114), plus DFT at a larger basis on unconstrained geometries | items 1–3 |
| Can GFN2-xTB be trusted for coordination number here? | the cluster benchmark — Imhof's own predicts ~7.5 kcal/mol error on 4→5 coordination, which is our contested observation | nothing |
| Does the re-formed carbonyl land where L1 puts it (2.53–2.63 Å)? | `recombine`, then `refstruct.py transfer` 7UHI as a check | items 1–2 |
| What protonates the leaving nitrogen? | step 3 as its own scan/NEB, proton source as a hypothesis set | item 6 |
| Is the N-terminal truncation question material? | residues 1–14 are modelled from no experimental coordinates, 9–13 Å from the zinc | unresolved, a referee will ask |

---

## 6. What is publishable now

**Nothing mechanistic.** Two structural results and a methods paper.

1. **Zn2 is vacant in wild-type CSR-1**, from the `SEQADV` argument above.
2. **CSR-1's Zn1 site measured against the B3 references**, including the
   Asp92 conformer comparison across its own three crystals (D179) with
   resolutions and B-factors.
3. **A methods/software note**, which should stay separate from the enzyme
   paper: restraint laundering, rigid-water constraints becoming convergence
   criteria, per-residue vs per-atom detection granularity, three upstream
   ASH/xtb bug reports, the silent-success class, and measured cost figures
   (750–760 ns/day for 49k atoms on an H100; 2.57 ps/h for GFN2-xTB QM/MM MD
   on 8 cores; ~26 min per ωB97X-D3(BJ)/def2-SVP/MM single point on 119 QM
   atoms).

---

## 7. Operational facts that cost time when forgotten

| Fact | Detail |
| --- | --- |
| The code-server pane dies | it is a Slurm `sys/dash` job with a walltime; the node/port **moves** (bun151:20685 → bun157:36073 on 2026-11-02). Read `location.href` first, every time |
| The compute node has **no outbound internet** | `curl` to files.rcsb.org returns 0-byte files. Download reference PDBs locally and measure them there |
| A GPU submission needs `-q gpu` | `run.slurm` hardcodes `--qos=normal`, which `gpu_cuda` rejects with `Invalid qos specification` |
| The browser shell **blocks** while a command runs | typing the next command sends it to the running command's stdin. Wait for the output file |
| Scratch | 272/300 GB and **889k/1M inodes** — the inode count is the tighter one |
| Local tests | `PYTHONPATH=src <mcube-test python> -m pytest -q tests`; MDAnalysis-gated tests only run in the ASH env on Bunya |
| Reading the GPU queue | `sinfo -p gpu_cuda` showing `mix-` is **PLANNED**, not draining: Slurm is holding those nodes for a higher-priority pending job, which is why a job can sit at reason `Priority` while the partition looks idle. `scontrol show node <n>` prints the real state. `drain$` is a genuinely unavailable node |
| `md` resumes | `resume: true` **in the section's inputs** - a config option, not a CLI flag; `cli.py` has no such flag. Same as `qmmm`, reading the `OpenMM_MD_state.xml` ASH writes as it runs. Opt-in on purpose (for a trajectory, continuing and replacing differ), but an unset `resume` with a state file present now warns rather than silently starting over (D183). `nvt_time` is the total wanted; minimisation and NPT are skipped on a resume (D182). So a long MD can be a chain of short jobs, which schedules far better on a partition whose nodes are held |
| A long walltime hurts backfill | asking 8 h for a 4.1 h run means it only starts when an 8 h window opens. A pending job's `TimeLimit` can be **lowered** (`scontrol update jobid=... TimeLimit=5:00:00`) but a user generally cannot raise it again, so it is not a free change |
| **CSR is not where the disk goes** | all eleven run directories are 14.5 GB of the 269 GB in use. Cleaning CSR is for clarity, not quota; the pressure is other projects (D180) |
| Deleting run output | use `docs/csr/cleanup.sh` - it lists by default, never deletes without `--delete`, and protects every `*_results/` directory. Leave a `DELETED.md` saying what went, when, how much, and which entry to see |
| **The evidence for a live claim is not redundant** | however void the run that produced it. `md_wt100`'s trajectory is void as an equilibration and is the primary evidence for D177, so it stays until a better measurement supersedes the claim (D180) |
