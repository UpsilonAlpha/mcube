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

| Histidine tautomers | **The two zinc ligands are HID with NE2 facing the metal in every run** — so the HE2-on-everything source file D21 found never reached a production system. 53 HIE, 150 HIE agree everywhere. **His112 is HIE in CSR_ligand/CSR_hydrolysed and HID everywhere since**, unattributed to any source; it is 15.5 Å from the zinc and mechanistically irrelevant, so set `112: HIE` in future builds and rebuild nothing (D189) |
| His164 HID vs HIE is still open | a chemistry question, not a parameterisation error: its NE2 is 3.2 Å from the developing oxyanion, so the tautomer decides whether the oxyanion has an H-bond donor (D7, still plan item 5) |

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

## 2a. The substrate pose — settled, and not to be re-litigated by reverting (D194)

**There are two poses in this project and the project chose the second one.**

| pose | Zn–O1x (side-chain OH) | Zn–O2x (carbonyl) | builds |
| --- | --- | --- | --- |
| hand-built near-attack, **the user's own** (D20) | 13.3–14.0 | **2.34–2.55** | `CSR_ligand`, `CSR_hydrolysed` |
| docked pose 2 **pulled onto the metal** by two hand restraints (D68) | 3.82 | 1.94 | `CSR_crystal_v2` |
| **Vina docked pose 2, as docked — chosen deliberately** (`model: 2`, D68) | 3.13 | **5.43** | `CSR_crystal`, `CSR_nac`, `CSR_unrestrained2`, **`CSR_wt`** |

`CSR_ligand`/`CSR_hydrolysed` are the user's own builds and the pose used for
most of the early work. **It was then deliberately set aside for the docked
pose.** Transferring it back is undoing a decision, not fixing a regression —
proposed on 2026-11-03 and withdrawn the same day.

**The carbonyl is unrestrained on purpose.** D98: *"the substrate gets one
restraint on the attack coordinate at the distance the docked pose already has.
Nothing asserts the carbonyl coordinates the zinc."*

**The one restraint that holds the substrate is on the attack coordinate**,
`[C:2:O01, C:1:C4x, 3.47]`. With it, the substrate held through 100 ns
(D100/D101, Zn–O2x 5.47 Å). Without it, it leaves.

**What is still open:** whether the docked pose is right at all. 7UHI/7UHJ put
L1's carbonyl at **2.53/2.63 Å from Zn1** (D123, which withdrew D40's claim
that inner-sphere coordination was unestablished), and the docked pose never
makes the carboxylate–Arg contact every other MBL substrate complex shows. That
is the blocking scientific question (D89) and it is to be settled on evidence,
**not** by reverting to an older build.

**And a substrate leaving an MM run is a statement about the restraint list,
not about the pose.** Mistaken for the latter twice — D89, then D192.

## 2b. Is there a usable QM/MM starting structure? NOT YET, but the blocker has moved (D195 measured, superseded in part by D197)

**Read 2d first.** `mich_eq5` now has the carbonyl distance, the attack
distance and the coordination number in ONE structure, which nothing in the
table below does. What it fails is the attack ANGLE. The table remains the
record for the 84 pre-rebuild frames.

`site.py check` over all 84 last frames. **One frame in the project is
accepted** — `CSR_crystal_v2/qmmm_results/QM_MM_lastframe.pdb` — and it has a
six-coordinate zinc with Asp92 OD2 at 2.04 Å and an attack angle of 26°.

**The two halves of a Michaelis complex are never in the same structure:**

| what | best frame |
| --- | --- |
| attack geometry (3.0–3.5 Å, 95–120°) | `naclong` 3.33/101°, `wt100b` 3.65/102° |
| carbonyl coordination (2.4–3.0 Å) | **only** `CSR_crystal_v2`, and only after QM/MM (2.56 Å) |
| four- or five-coordinate zinc | every frame is 5 or 6; **5 is now allowed** once the substrate is bound (D196), so what disqualifies these frames is the 6, via Asp92 OD2, a bulk water, or Gln93 |

**Confirmed: a carbonyl held coordinated through MM stays under QM/MM.**
`CSR_crystal_v2` went in at 1.94 Å and xTB *relaxed it to* **2.56 Å** — inside
0.07 Å of 7UHI/7UHJ's 2.53/2.63 Å. Every frame whose carbonyl was uncoordinated
in MM lost the substrate under QM/MM (naclong 5.16→6.20, D111). So GFN2 is not
overbinding or expelling it, and **the carbonyl needs restraining in MM** — at
2.4–3.0 Å, not at 1.94.

**Read the previous paragraph as hysteresis, not as adjudication (D197).** Held
in, it stays in; left out, it goes further out. The one inward move on record is
5.16→4.37, which never approaches the sphere. **No QM/MM run in this project has
crossed the coordination boundary in either direction**, so 1.94→2.56 is a
bond-length relaxation *within* the coordinated basin and is **not** evidence
that GFN2 would reject a wrongly-placed donor. The mechanism is in D96's
analysis: MM MD overfavours octahedral zinc (PMID 34748329), DFTB3 QM/MM MD
rectifies it, ours did not, and the reason is that **ligand exchange at Zn(II)
is slower than tens of ps**. A QM/MM MD therefore cannot settle whether a fifth
donor belongs on the metal; only a thermodynamic comparison can, which is P1.

**`CSR_ligand` and `CSR_hydrolysed` cannot be checked:** `B:1:ZN not in the
structure`. The user's own builds predate the current chain naming, so no
`site.yaml` selector resolves. Remapping their chains makes them checkable
without re-running anything, and they are the only frames whose pose has the
carbonyl near the metal with no frozen MM constraint.

## 2c. The coordination number, and the rebuild (D196, the user's decision)

**Five is allowed once the substrate is bound. The target is still four.** The
resting site is four-coordinate — 2 His + 2 waters (D120). On binding, the
substrate **displaces** one of those waters; the other becomes the nucleophilic
hydroxide, and the carbonyl takes the vacated position: 2 His + hydroxide +
carbonyl = **four**. Five is not rejected, for two separate reasons:

- GOB-18's B3 site expands four- to five-coordinate on substrate binding
  (`reference-structures.md`), so a retained second water is a real
  possibility, not a build error.
- Asp92 OD2 sits at **2.80 Å** in the crystal as a second-shell contact (D120).
  A frame where it drifts inside the 2.5 Å measurement radius counts as five
  here while being nothing of the kind chemically.

Both `michaelis_hypotheses` coordination criteria are now `[4, 5]`. The old
`[4, 4]` was rejecting correct structures for being correct.

**`eq5` is the rebuild that follows.** D195's finding was that each run had ONE
of the three Michaelis properties, so `eq5` restrains all three at once, from
`CSR_crystal_v2`'s input pose (1.94 → 2.60 Å is a small move; 5.43 → 2.60 is
not):

| restraint | why |
| --- | --- |
| `[C:2:O01, C:1:C4x, 3.47]` | the attack coordinate, which `wt100c` **replaced** rather than added to (D194) |
| `[B:1:ZN, C:1:O2x, 2.60, 500]` | the carbonyl on the metal. Within 0.04 Å of what GFN2 chose for itself (2.56) and 0.03 Å of L1's 2.53/2.63. **Not 1.94** — nothing measures that, and holding it there is what made crystal_v2's attack angle 26–30° |
| `[B:1:ZN, A:92:OD1, 3.81, 500]` | OD1 pinned out, as on every arm since D130 |

**500 kcal/mol/Å², not 50.** D191 measured a 2.80 Å restraint sitting at a mean
of 2.44 Å — a 0.36 Å offset is a systematic ~36 kcal/mol/Å force, not thermal
spread, and offset goes as F/k. At 500 the restrained mode's period is ~53 fs, so
`eq5` sets `nvt_timestep: 0.002` — 26 steps per period, where the 4 fs this
run directory inherits would give 13. The MM site
is deliberately rigid: MM cannot be trusted for coordination (D174), and QM/MM
is what decides it (confirmed in D195).

**Measured, not assumed (2026-11-03).** Re-running the criterion over all ten
NVT frames moved four coordination verdicts to PASS (`crystal_v2`, `nacgpu`,
`unrestrained`, `unrestrained2` — all 5) and **changed no overall verdict**:
every one still fails `Zn - carbonyl O2x`. The blocker was never the
coordination number. And `CSR_crystal_v2`'s QM/MM sphere is `D:1:HOH 1.98,
Asp92 OD2 2.04, His88 2.06, His164 2.08, hydroxide 2.21, carbonyl 2.56` —
removing `D:1:HOH` leaves exactly five, and that frame then passes both the
carbonyl distance and the coordination number, failing only the attack
geometry. So `remove_residues: [D:1:HOH]` is the right selector for this
system, measured in the structure it will be applied to.

**A rebuild is a new run directory.** `preprocess` takes no `--section`, so
`remove_residues: [D:1:HOH]` and `coordination_force: 500` go in the new
directory's own `preprocess.inputs`, with `eq5` as the `md` section on top. The
water selector must be checked against the **merged input**, not the final
system — preprocess refuses a name mismatch rather than deleting the wrong
residue.

**`preprocess.inputs.remove_residues` is new in the library** (D196). There was
no config route to displacing a water; the only way was to edit the input PDB
by hand. It takes `CHAIN:RESSEQ:RESNAME`, runs before solvation (deleting a
water afterwards just lets the solvent refill the position), takes ANISOU and
CONECT records with it, and refuses the whole request if any selector misses.

## 2d. `mich_eq5`: every restraint held and the structure is still rejected (D197, measured)

**Status: measured, open.** 31721717 COMPLETED in 35.1 min, 10 ns at 2 fs on an
L40S, no NaN.

**All seven restraints held to better than 0.11 A** - carbonyl O2x 2.645 (set
2.60), Asp92 OD1 3.782 (3.81), O01-C4x 3.458 (3.47), hydroxide 2.035, His164
2.166, His88 2.193, Asp92 OD2 2.743. Nothing collapsed, nothing dissociated.

**`site.py check` under `carbonyl_on_metal` at `--at=reactant`: REJECTED.**
Four pass - carbonyl 2.65 (2.4-3.0), nucleophile 2.04 (1.9-2.2), coordination
**4** (4-5), attack distance 3.46 (3.0-3.5). **This is the first structure in
the project to hold all three of D196's criteria at once**, so D195's blocker is
gone. Three fail: the closest substrate oxygen is O1x 2.06 where O2x is wanted,
the **O01-C4x-O2x angle is 56 degrees against a 95-120 window**, and Asp92
OD1-O01 is 5.81 A against 2.5-3.2.

**One unrestrained atom explains the angle.** `mer`'s **O1x** - the
6-(1-hydroxyethyl) alcohol oxygen, identified from O1x-C2x 1.39 A with C2x
bonded to two carbons, against O2x-C4x 1.24 with C4x-N1x 1.35 for the lactam
carbonyl - coordinated the zinc at **2.06 A**, closer than the restrained
carbonyl and as close as the anionic hydroxide. It sits between them, both at
76.8 degrees from it, squeezing the OH-Zn-O2x angle to 76.8 where **106.7-128
degrees** is needed for a legal attack angle. 107 degrees is merely tetrahedral,
so `carbonyl_on_metal` and Burgi-Dunitz are **not** mutually exclusive - the
fifth donor is what closes the angle.

O1x was at 3.13 A in the October pose and outside the 3.0 A sphere at
preprocess, so nothing restrained it: **`fill_restraints` can only restrain what
is inside the radius when it looks**, and a donor arriving from outside during
the MD is the one case detection-time restraints cannot cover. Same failure as
`prod100b` (D173), and `wtpostc` was already watching `Zn-mer O1x`.

**Whether O1x on the metal is real is NOT settled, and MM cannot settle it.**
Not for the reason given here before (corrected in D199): the MD **did** apply
the Li/Merz IOD parameters, so 2.06 A is not a collapse past a missing
repulsion - it is close to the 1.89 A minimum those fitted parameters define,
i.e. **the fitted Zn-O distance applied to the wrong oxygen**. `ions.py` says
why: a 12-6 term knows only distance and atom type, so it cannot tell an
alcohol from a hydroxide and rewards every oxygen in reach; only the
coordination restraints can, and O1x sat 0.13 A outside the 3.0 A detection
radius so none was placed on it. A carbapenem's 1-hydroxyethyl genuinely
chelating the metal is also not absurd. QM/MM decides it, per the standing
rule.

**Asp92 OD1-O01 needs its own restraint.** 5.81 A here against 4.76 in the
QM/MM frame. Pinning OD1 3.81 A from the **zinc** fixes the distance and not the
direction, and it pointed away; the restraint has to be OD1-O01 itself.

**The `md` stage's before/after coordination report is keyed per residue and so
can mislead.** It announced `LEFT B:1:ZN - C:1:O2x (was 1.94)` when O2x is at
2.645 A, inside the sphere and exactly where it was restrained: O1x became
`mer`'s closest atom and displaced O2x from a table that shows one atom per
residue. So the printed "after MD" sphere **omits a genuine donor**, the implied
count is short by one, and a restrained contact can be reported as broken when
it is not. Not an imaging artefact - the imaged and unimaged last frames agree
to three decimals. `postprocess` has `per_residue: false` for this; the `md`
reporter has no such option.

**Measured over all 10,000 frames (31736226), which settles the contamination
question.** At the 2.5 A radius: hydroxide **100%** (mean 2.03), His164 100%
(2.13), His88 100% (2.21), **mer O1x 97.9%** (2.18), Asn225 OD1 17.9% (2.81),
**mer O2x 0.3%** (2.60), Asp92 OD2 0.0% (2.77). Coordination number 4 in 8328
frames, 5 in 1638.

**O1x was there from the start, so no frame of this run is usable** - there is
no early clean window to take to QM/MM. And **the carbonyl is inside 2.5 A in
0.3% of frames**: its mean sits exactly at the 2.60 it was restrained to, and
2.60 is outside the 2.5 A bond radius *by construction*. By this project's own
standard the lactam carbonyl is not a ligand at any point in the run, while an
unrestrained alcohol is one for 97.9% of it.

**[WITHDRAWN 2026-11-03, D197 corrections, on the user's challenge. I called
this a conflict between two settled numbers; it is not one.** D124 already
settled it: thermolysin and CPA hold the nucleophile and the carbonyl on one
zinc, and the condition is that the carbonyl is **loosely** coordinated - one
oxygen tight (2.01, 2.19), the other loose (2.41, 2.98). A carbonyl at 2.60 Å
that is not counted as a bond at 2.5 Å **is** that loose oxygen, and D186's two
radii exist to express exactly this. The 0.3% occupancy is the design working.
Nothing needs changing, and at 2.60 Å "loosely coordinated" and "polarised" are
not a distinction these measurements can resolve.]

It also explains how `site.py` reports `PASS coordination number 4`: the four
are hydroxide + 2 His + **O1x**, with the carbonyl outside. Correct number,
wrong membership, and the criterion cannot see it because it counts donors
without requiring the carbonyl to be one. Testing properties one at a time can
pass a composition that would be rejected if asked directly.

**Asn225 OD1 is a second uninvited donor** - 17.9% occupancy, min 2.02 A,
unrestrained, and on no earlier arm's monitor list.

**Everything restrained behaved, across the whole run.** O01-C4x mean 3.47
(target 3.47), Asp92 OD1 mean 3.80 (3.81), Asp92 OD2 mean 2.77 against the
crystal's 2.80 at **0.0%** inside 2.5. **No ligand turnover**: the carboxylate
stayed at 5.93/7.05, so unlike `prod100b` this is one intramolecular donor, not
a reorientation. Incidentally, Lys226 NZ-O01 is inside 3.5 A in only 10.3% of
frames (mean 4.93), and Gln93 NE2-Asp92 OD1 averages 7.74 A - the Q of RQK has
no contact to decompose on this arm.

**Operational.** `postprocess` takes the **full config key** as its section
(`run.slurm postprocess postprocess_michpost`, not `michpost`); with the short
name the section is silently not found, the stage uses the plain `postprocess`
inputs, and the error reads `` `md` has no results `` as though the MD were
missing. The run wrote **12.7 GB** because every trajectory is kept twice, raw
and imaged, taking scratch from 269 to 286/300 GB.

## 2e. Frame 4374 is the best starting structure in the project (D198, measured)

**Status: measured, open.** Found while building the cluster benchmark, not
looked for: it is the frame of `md_eq5`'s 10,000 with the longest Zn-O1x
(3.74 A), and `site.py check` passes **five of seven** criteria on it against
four on the last frame.

**It is the first structure in the project to pass `closest substrate O`** -
carbonyl O2x at 2.58 A with O1x out at 3.74 - which is the criterion that
encodes what `carbonyl_on_metal` actually claims. Also passing: carbonyl 2.58,
nucleophile 2.07, coordination 5, attack distance 3.50. Still REJECTED on the
attack angle (68 deg) and Asp92 OD1-O01 (5.67 A).

**The attack angle tracks O1x across three independent frames** - 1.92 A/39
deg, 2.06 A/56 deg, 3.74 A/68 deg - so D197's crowding mechanism is confirmed.
**It is also not sufficient.** In frame 4374 the five donors are hydroxide 2.07,
His164 2.13, **Asn225 OD1 2.19**, His88 2.23, carbonyl 2.58: Asn225 took the
site O1x vacated, the sphere stayed five-coordinate, and the angle opened by
only 12 deg. So what closes the angle is not O1x specifically but that this site
keeps a **fifth donor** in MM whatever is done, and a five-coordinate zinc puts
the nucleophile and the carbonyl too close together. Getting into the 95-120
window needs a genuinely four-coordinate sphere, which no MM run has produced.

**The production QM/MM starts here, and that is a selection, not a sample.**
31738709, `qmmm_michqm`, 50 ps of 300. Because QM/MM is hysteretic in this
project (see 2b), the starting donor set is what the run will report, so the
frame was chosen rather than taken from the end of the MD - and an extremum of a
distribution describes the reactant state we chose, not an equilibrium average.
Say so wherever this run is quoted.

## 3. Contradictions caught late — the failure mode this file exists to stop

| The contradiction | Days unnoticed | How it was caught |
| --- | --- | --- |
| D145 quoted 6NC5's Asp distance as "the crystal's" when **D108 had already ruled 6NC5 out as a WT reference** | 8 | measurement, D179 |
| D120 prescribed the coordination fix; it never worked | 5 | occupancy over a trajectory, D177 |
| D48 settled the resting state as Zn-OH₂; the whole `CSR_nac` chain ran in the other state | 8 | D112 |
| `coordination_from_rows` said restraints are per atom; `find_coordinators` called it per residue | ~6 weeks | D175 |
| A design decision (Gln93 out of the QM region) inferred from config files instead of read from the log | — | user, twice: D165, D171 |
| D89's "dissociation implies a bad pose" was examined and **discarded** as an artefact of my own restraint set; **D192 re-derived it** from the same observation | ~3 weeks | user, D194 |
| `wt100c` **replaced** `wt100b`'s attack-coordinate restraint instead of adding to it; the `site.yaml` comment says "**One** restraint beyond the sweep" and never names what it overwrote | 1 | user, D194 |

Every one of these was a *reading* failure, not a computing failure.

---

## 4. The plan (audit `audit-2026-11-02.md`, adopted in D178)

Items 1–3 gate everything below them. **Status is maintained here.**

| # | Item | Status |
| --- | --- | --- |
| 0 | Asp92's placement in the starting structure | **settled (D179)** — 2.80 Å supported; `prod100c` unblocked |
| 1+2 | Experiment-anchored Michaelis complex, both endpoints on one topology | **open as a scientific question, not blocked on tooling (D194).** The docked pose binds with the carbonyl at 5.43 Å and one attack restraint; 7UHI/7UHJ put L1's at 2.53/2.63 Å (D123). Deciding between them is the blocking question (D89) and needs evidence, not a revert to `CSR_ligand` (§2a) |
| 3 | Equilibration that preserves the site | **FAILED, both arms (D191), cause found (D194): the attack-coordinate restraint was replaced, not added to.** Re-run needs `[C:2:O01, C:1:C4x, 3.47]` **and** `[B:1:ZN, A:92:OD1, 3.81, 50.0]` in the list. Still to fix alongside: 50 kcal/mol/Å² does not hold Asp92 OD2 at 2.80 (mean 2.44/2.37), a bulk water holds the free site 97.0%/89.3%, and `prod100c` NaN'd at 22.68 ns — the product arm's second NaN, so it predates D186 |
| 4 | Nucleophile identity, origin, protonation state | not started — blocked on 1–3. D48 and D107 disagree **in sign**. One datum arrived anyway: **Lys226 NZ–O01 at 2.72 Å mean in 100% of 5000 frames** (D191) |
| 5 | PROPKA/H++ and the His164 HID/HIE sensitivity check | **blocked on one pip install**: neither `propka` nor `pdb2pqr` is in the ASH env, both are pure Python, so `pip install propka pdb2pqr` is the whole of it — the user's env, so their call (D185). Open since **D7** |
| 6 | scan → NEB → TS → frequencies → IRC, ≥3 snapshots | not started — the paper's actual content |
| 7 | Cluster-model DFT-vs-GFN2 benchmark, Imhof-style | **scoped (D187)**: 37 atoms of site (Zn, two 4-methylimidazoles, acetate, hydroxide, water), net 0, plus a ~15-atom substrate fragment giving −1. The two `sp` additions are **built** (D188): `cluster: true` and `sp_engine: xtb|dft`. Still needs the cluster PDB and reading the ligand's connectivity from `mer_crystal.pdb`'s CONECT records rather than guessing it. Open since **D9** |
| 8 | Protonation of the leaving nitrogen | not started |
| 9 | Per-residue roles reaching R, Q and K | **decided and applied (D190)**: Arg90, Lys226, Arg229 and Gln93 **all in MM**, so charge deletion decomposes the set as one. QM charge 0 on both arms, set explicitly on all nine QM/MM-family sections. Deletions must still be run at **both ends** of the coordinate |
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

### The gate was run on 2026-11-03 and BOTH ARMS FAILED (D191)

| criterion | reactant `wt100c` (100 ns) | product `prod100c` (22.68 ns, NaN'd) |
| --- | --- | --- |
| 1. four-coordinate at 2.5 Å | **FAIL** — 6 in 3523/5000, 5 in 1477 | **FAIL** — 6 in 852/1134, 5 in 281 |
| 2. Asp92 OD2 near 2.80 | **FAIL** — mean 2.44, max 2.97 | **FAIL** — mean 2.37, max 2.75 |
| 3. Asp92 OD1 near 3.81 | pass — 3.73, 0.3% | pass — 3.82, 0.3% |
| 4. product O4x/O5x near 2.9 | — | pass — 3.11 / 3.15, never inside 2.5 |
| 5. hydrogens clean | pass | pass |

**And the criterion nobody wrote down, which is the one that matters: the
substrate has to still be there.** It is not. `O01–mer C4x` is inside 3.5 Å in
**0.2% of 5000 frames**, with mean 19.95 Å and max 50.69 Å — meropenem leaves
in roughly the first 0.2 ns and diffuses across the box. The product's ligand
goes the same way (`Zn–O6x` 6.44–8.10 Å). Nothing restrains the substrate,
because its closest zinc contact (`mer O1x`, 3.13 Å) is outside the 3.0 Å
radius and was correctly excluded by this section's own rule.

**Criterion 6, from now on: the substrate's attack distance must stay inside
3.5 Å in ≥90% of frames, or the run is not a Michaelis complex and nothing
downstream of it means anything.**

**The bulk water was not transient either.** E:2920 sits on the reactant's zinc
for **97.0%** of 100 ns and E:3189 on the product's for **89.3%** — resident,
not passing through, and the reason both arms read six-coordinate. The
paragraph above excuses a water that comes and goes; it does not excuse these.

**Why criterion 2 failed is now quantitative, and it is the force constant.**
Asp92 OD2 is restrained harmonically at 2.80 Å with k = 50 kcal/mol/Å², and
over 100 ns its distribution's *maximum* is 2.97 while its mean is 2.44: the
distance hardly ever reaches the restraint's own minimum. At 2.44 Å the
restraint pulls with 2·50·0.36 = 36 kcal/mol/Å and still loses, which is what a
formal +2 point charge against a carboxylate oxygen at 2.4 Å looks like without
a bond. D186 settled *which* atoms are restrained; **50 kcal/mol/Å² is the open
question, and holding a mean within ~0.05 Å of 2.80 needs of order 500.**

`prodqm` and `wtqm` were **not** submitted, although the standing authority was
to submit on a pass: 50 ps of QM/MM MD from a frame with no substrate in the
site is an expensive converged meaningless number (the D152 failure mode).

### Is the dissociation real, or a periodic-image artefact? SETTLED: REAL (D192)

The claim above rests entirely on `wtpostc`'s monitored-contacts table:
`O01-mer C4x` 0.2% inside 3.5 A, min 3.14 / mean 19.95 / max 50.69 A, with
`Zn-mer O1x` and `Zn-mer O2x` saying the same thing. **It was written up as
settled before the obvious check was done, and that was wrong.**

`contact_occupancy` computes every distance as `separation(..., cell)`, and
`separation` applies the minimum image **only when `cell` is not None**.
`cell` comes from `coordination_rows(pdb_file)` -> `_coords_and_cell`, which
returns a cell only if the topology PDB carries a `CRYST1` whose three angles
are all 90 within 1e-3; anything non-orthorhombic or missing yields `None` and
**no wrapping correction at all**. Whether
`md_wt100c_results/NVTtrajectory_firstframe.pdb` satisfies that has not been
checked. If it does not, a ligand that never left the site but sits near a box
face reads as tens of Angstrom, which is exactly what maxima clustered at
50.69 / 51.25 / 51.62 A look like.

The maxima cannot discriminate. For a cubic box of side L the largest possible
minimum-image distance is L*sqrt(3)/2; this system is ~10,000 waters, so
L ~ 70 A and the ceiling is ~61 A. 51.6 A is therefore consistent with genuine
diffusion **and** with an uncorrected wrap.

**What does corroborate dissociation, independently of the ligand's absolute
position:** the reactant's coordination listing contains **twenty-plus bulk
waters at 3.1-3.5 A from the zinc** and `mer O1x`/`O2x` nowhere in it. That
scan examines every atom within cutoff+margin of the metal, so the substrate is
absent from a neighbourhood that bulk solvent has filled. A periodic wrap
displaces the ligand; it does not put water into the ligand's pocket.

**The test that settles it, in one read:** the per-frame series, which
`contact_occupancy` already returns as `_distances` and writes to the CSV.
Dissociation gives a smooth random walk away from 3.14 A; a periodic-image
artefact gives instantaneous jumps between ~3 A and ~50 A with no density in
between. Check the `CRYST1` angles of the topology PDB at the same time.

**Checked, and it resolves against the doubt (D192).**
`md_wt100c_results/NVTtrajectory_firstframe.pdb` carries
`CRYST1 79.028 79.028 79.028 90.00 90.00 90.00` — orthorhombic, so
`_coords_and_cell` returns a cell and `separation` applies the minimum image.
The ceiling for a cubic box is L·√3/2 = **68.4 Å**, and the observed maximum of
50.69 Å sits well inside it. **The dissociation is real and no longer
provisional**, and the per-frame series is not needed to establish it.

**The cause is one line of config, not the pose (D194).** `nac` and `wt100b`
carried `restraints: [['C:2:O01', 'C:1:C4x', 3.47]]` — the attack coordinate —
and with it the substrate held through 100 ns. **`wt100c` and `prod100c`
replaced it** with `[['B:1:ZN', 'A:92:OD1', 3.81, 50.0]]`. `restraints` is a
list and it was treated as a slot; both belong in it. The Asp92 OD1 restraint
was itself sound (D120/D130, the bidentate pivot twice), it simply went in
*instead of* rather than *beside*.

D192's geometry is right — carbonyl `O2x` 5.43 Å, `O1x` 3.13 Å nearest,
confirmed by connectivity — but **its diagnosis is withdrawn**: 5.43 Å is
docked pose 2 as docked, chosen deliberately, with the carbonyl unrestrained on
purpose (§2a).

### `constraint_exclude_residues` is not for the molecule under study (D193)

The reactant's ligand was named in `constraint_exclude_residues: {C: [1]}`,
kept deliberately and documented as costing nothing. **The reasoning was wrong
but the line was inert**, and the first version of D193 wrongly blamed it.
`constraint_exclude` is tested inside the loop over `shell`, which is already
filtered to within the cutoff, so the key can only act on residues inside it —
and `mer`'s nearest atom was 3.13 Å. **Neither run's log contains a single
`No restraint for ... (constraint_exclude_residues)` line.** The reactant's
restraint list lacks the ligand because nothing was *found*, not because
anything was *excluded*; the product's contains `O4x`/`O5x` because 2.92/2.94 Å
are inside the cutoff. mcube restrained everything within 3.0 Å on both arms,
exactly as built. **D192 is the sole cause.**

**The rule:** this key is for atoms near the metal *by accident* — D8's
arginine, restrained as though it were a ligand, is why it exists. A substrate
carbonyl is the opposite case. "We are not sure it coordinates" is an argument
for holding it where the structure puts it and letting QM/MM decide, never for
letting it go; omitting a restraint is not neutrality.

**And a warning for "you excluded something inside the cutoff" is worse than
redundant**: it would fire on every exclusion that does anything, D8's arginine
included. The useful guard is the inverse — warn when an excluded residue has
*no* atom inside the cutoff, i.e. dead config that reads as load-bearing. Not
built; offered.

## 5. Open questions, with what would settle each

| Question | What would settle it | Blocked by |
| --- | --- | --- |
| Is the nucleophile hydroxide or water, and how does it form with no Zn2? | equilibrate each state, keep the one that holds the crystal's sphere (D113/D114), plus DFT at a larger basis on unconstrained geometries | items 1–3 |
| Can GFN2-xTB be trusted for coordination number here? | **partly answered (D198)**: on the O1x contact, dd = GFN2 - DFT = **-5.92 kcal/mol**, i.e. GFN2 *under*-binds the hydroxyethyl - the opposite sign to the worry. Same magnitude as Imhof's 7.55 for AM1/d, so still a real error, and still one contact on one geometry pair rather than a statement about coordination number in general. Single points on MM geometries, not binding energies | the cluster benchmark — Imhof's own predicts ~7.5 kcal/mol error on 4→5 coordination, which is our contested observation | nothing |
| Does the re-formed carbonyl land where L1 puts it (2.53–2.63 Å)? | `recombine`, then `refstruct.py transfer` 7UHI as a check | items 1–2 |
| What protonates the leaving nitrogen? | step 3 as its own scan/NEB, proton source as a hypothesis set | item 6 |
| Does meropenem's 1-hydroxyethyl OH really coordinate the zinc, or is the 2.06 A an MM collapse? | QM/MM from this frame and let GFN2 choose, OR restrain O1x out and re-run `eq5` to see whether OH-Zn-O2x opens past 107 degrees. The two answers differ in what the paper's Michaelis complex IS, so this is not a tuning question (D197) | nothing |
| Is the 56 degree attack angle fixable by geometry alone? | an OD1-O01 restraint plus O1x held out; if the angle still will not open with the carbonyl at 2.6, `carbonyl_on_metal` itself is in doubt rather than the structure | the question above |
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
| The code-server pane dies | it is a Slurm `sys/dash` job with a walltime; the node/port **moves** (bun151:20685 → bun157:36073 on 2026-11-02, bun157:61214 on 2026-11-03). Read `location.href` first, every time |
| The compute node has **no outbound internet** | `curl` to files.rcsb.org returns 0-byte files. Download reference PDBs locally and measure them there |
| A GPU submission needs `-q gpu` | `run.slurm` hardcodes `--qos=normal`, which `gpu_cuda` rejects with `Invalid qos specification` |
| Running a command in the browser shell | do **not** type it: put it in an editor tab and use ctrl+shift+p → "Terminal: Run Selected Text In Active Terminal", which writes the selection into the pty with real Enters and no truncation. Typed newlines are not Enter, and long typed payloads truncate silently. Write the script without quotes, braces or parentheses, because monaco auto-pairs them |
| The browser shell **blocks** while a command runs | typing the next command sends it to the running command's stdin. Wait for the output file |
| Scratch | 269/300 GB and 655k/1M inodes on 2026-11-03, after ten redundant `.dcd` files were deleted (6 GB). Both have headroom now; the inode count had been the tighter one at 889k |
| Local tests | `PYTHONPATH=src <mcube-test python> -m pytest -q tests`; MDAnalysis-gated tests only run in the ASH env on Bunya |
| Reading the GPU queue | `sinfo -p gpu_cuda` showing `mix-` is **PLANNED**, not draining: Slurm is holding those nodes for a higher-priority pending job, which is why a job can sit at reason `Priority` while the partition looks idle. `scontrol show node <n>` prints the real state. `drain$` is a genuinely unavailable node |
| **`resume: true` + a killed run in the same output directory** | resumes the **dead run's system**, because a state file carries its own particle count. 31738904 built the right 47,979-atom system then died on `Called setPositions() on a Context with the wrong number of positions`, loading the cancelled 34,311-atom run's `OpenMM_MD_state.xml`. Clear `protdft/<section>/` before the first batch of anything with `resume` pre-set. Loud, but the message names positions, not the state file's provenance (D199) |
| `md` resumes | `resume: true` **in the section's inputs** - a config option, not a CLI flag; `cli.py` has no such flag. Same as `qmmm`, reading the `OpenMM_MD_state.xml` ASH writes as it runs. Opt-in on purpose (for a trajectory, continuing and replacing differ), but an unset `resume` with a state file present now warns rather than silently starting over (D183). `nvt_time` is the total wanted; minimisation and NPT are skipped on a resume (D182). So a long MD can be a chain of short jobs, which schedules far better on a partition whose nodes are held |
| A long walltime hurts backfill | asking 8 h for a 4.1 h run means it only starts when an 8 h window opens. A pending job's `TimeLimit` can be **lowered** (`scontrol update jobid=... TimeLimit=5:00:00`) but a user generally cannot raise it again, so it is not a free change |
| **CSR is not where the disk goes** | all eleven run directories are 14.5 GB of the 269 GB in use. Cleaning CSR is for clarity, not quota; the pressure is other projects (D180) |
| **The run directories live in `CSR/`** | all thirteen, since 2026-11-03: `CSR/CSR_wt`, `CSR/CSR_crystal_v2`, `CSR/CSR_mich` and so on, with **nothing** left at the top of scratch. `CSR_mich` went in last, after its MD and its postprocess had both finished, because Slurm records a `WorkDir` at submit time and moving a running job's directory breaks it. Paths in the decision log and job inventory are as they were when the work ran, so prefix a `CSR_*` path with `CSR/` when following one. The absolute paths **inside** the moved `config.yaml`/`.slurm` files were rewritten (`.premove` copies kept) because the pipeline reads those; `.out` logs were not |
| Deleting run output | use `docs/csr/cleanup.sh` - it lists by default, never deletes without `--delete`, and protects every `*_results/` directory. Leave a `DELETED.md` saying what went, when, how much, and which entry to see |
| **The evidence for a live claim is not redundant** | however void the run that produced it. `md_wt100`'s trajectory is void as an equilibration and is the primary evidence for D177, so it stays until a better measurement supersedes the claim (D180) |
