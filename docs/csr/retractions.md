# CSR-1: what the decision log has overturned

**Why this file exists.** `csr-decision-log.md` is append-only and 207 entries
long. Its corrections live *inside* later entries, and in roughly forty cases the
corrected entry still **leads with the wrong claim** in its bold headline. A grep
therefore lands on the retracted half. That has caused four re-derivations of
discarded arguments — D130, D179, D194, D206 — and it happened again on
2026-11-05, twice in one session.

**How to use it.** Before quoting any decision-log entry, find it here. If it is
listed, read the superseding entry first. For a *current* value, do not use the
log at all — use `state.md`.

Built 2026-11-05 by auditing all 207 entries in five slices. Keyed by entry id,
which is stable; line numbers drift and are omitted deliberately.

---

## 1. The five open contradictions — ALL RESOLVED 2026-11-05 (D208)

Resolved by measuring the deposited structures with `refstruct.py shell`, which
nobody had done, and by reading the restraint builder. Kept here because the
wrong values are still quoted throughout the log.

**(a) Asp92's distance from the zinc: 4.20 A (OD2) / 4.22 A (OD1).** Measured in
6NC5; 6NC5's ZN302 is four-coordinate (His87 NE2 2.22, His163 NE2 2.09, HOH485
1.90, HOH493 2.14) with no aspartate. **D145 was right.** The **2.80 A** quoted
in a dozen entries is our own model's number and appears nowhere in the file.
**D184's 2.68 A is void**: it superposes 6NC5's zinc into 6DN4, and 6DN4 and
6DQ2 both contain **zero metal atoms**. **D151 was right** that 6NC5's ZN301
holds that carboxylate at 1.95 A as a Zn2 ligand, with HOH493 bridging both
zincs — so the 4.20 A carries the caveat that Zn2 is occupied, and it remains
the only measured value because the wild-type files are apo.

**(b) "No B3 enzyme has Asp120 on Zn1" — now measured.** 1SML 4.08, 5K0W 4.17,
5B15 4.05, 7UHJ 4.02 A from Zn1; 2.06–2.13 A from Zn2, monodentate in all four.
6NC5's 4.20 A is within 0.03–0.18 A of every one. Two by-products: **5K0W has
four zincs, both sites in both chains**, so "GOB-18's single zinc sits at Zn2"
is not what the deposit shows; and GOB-18's Zn1 is the HQH site, Gln98 OE1 2.04.

**(c) Does GFN2 over-coordinate the zinc? Open, and now secondary.** All four
closures are retracted or reclassed. But our runs put this contact at 1.77–2.11
A against 4.02–4.20 A in four experimental structures, so the geometry is wrong
whatever produced it; which component to fix is a separate question.

**(d) What put the fifth ligand on the zinc: two necessary causes.** The 3.0 A
cutoff admitted a second-shell contact, **and** the restraint builder's
per-residue default then held only the nearer carboxylate oxygen, leaving the
other free to pivot in. D120/D130/D174 blamed the cutoff alone; D175 named the
default and was substantially right; D176 over-corrected it. Fixed in
`find_coordinators` (`per_residue=False`).

**(e) D152's citation.** The 0.58 A Lys226/O4x overlap is **D148's** finding.
Corrected in place.

**The consequence: Asp92 should never have been in a restraint list at all.** At
4.20 A it is outside the 3.0 A "affected by the metal" radius and far outside
the 2.5 A coordination radius. Every structure built to hold it at 2.80 A has
that contact 1.4 A too short; the ones that collapsed to 1.77–2.11 A are 2.1–2.4
A too short.

---

## 2. Entries whose headline still asserts something later overturned

The grep traps. Read the right-hand column before quoting the left.

| entry | what its headline/lead claims | overturned by |
|---|---|---|
| **Standing rule 3** | xTB over-coordination "is no longer the open question — D58/D60 settled it"; "Asp92 binds the zinc at every level of theory" | D58 retracted (D59, D206). **This is the read-before-any-run block.** *Now marked.* |
| **D8** | Zn coordination cutoff 3.0 Å — "Revisit: **no**" | standing rule 1: "2.5 is defensible, 3.0 is not". The 3.0 Å laundered the carbonyl restraint (D68) |
| **D21** | a table headed **"Crystal"** with His164 NE2 at 1.73 Å | D23 — rebuilt from 6NC5; the `CSR_ligand` lineage is superseded |
| **D29** | Lys226 protonates the Zn–OH⁻ at once; Asp92 bound in 99.5% of frames | D38/D41/D48 — computed in the degraded six-coordinate state |
| **D35** | the transfer is −13.1/−15.5 kcal/mol at DFT, so H4 is dropped | D48 — −5.23; "roughly two-thirds was an artefact"; the hydroxide is "not comfortably excluded, only disfavoured" |
| **D38** | the site is "over-coordinated", the complex "probably wrong"; there is no hydroxide in the model | D41 — "wrong, as stated": the motif is **H-R-H**, two protein ligands, so two waters is canonical |
| **D40** | our mono-zinc site is "inference, not literature"; the non-truncated form is inactive | D41 — the zinc position **is** crystallographic; the inactivity claim was "relayed from a literature summary and not verified" |
| **D48** | "the Zn-bound oxygens still cannot attack" — 38–40°, so H2/H3/H4 are geometrically dead | its own table still had Zn–O2x **constrained at 2.90 Å**; D73 gets O01–C4x 3.01 Å at **84°** once released |
| **D51** | Asp92-as-base is not accessible, "tested twice from two different starting states" | D67/D68 — both tests descend from the condemned 1.94 Å pose; D64 — both sat in one arbitrary docking basin |
| **D52** | the zinc activates the electrophile; the proton leaves down a water wire | D53 (no MBL precedent), D54 ("withdrawn as implausible") |
| **D57** | Asp92 is "a genuine zinc ligand, not loosely associated", minimum 2.18 Å | both axes of that curve are our own model; §1(a), §1(b) |
| **D58** | "the six-coordinate zinc is **real**: DFT binds Asp92 harder than xTB" | D59, D206. *Marked inline; headline still reads the old way* |
| **D60** | DFT at 96.2% efficiency; a DFT/MM single point costs 4–5 h on 16 cores | D70 — a **mixed** job; 37% on 64 cores, ~26 min per single point, "no speedup factor should be quoted". The 4–5 h killed tier-B DFT and drove the GPU detour |
| **D61** | both removable ligands tested, neither helps; 22.4 kcal/mol for the carbonyl | D64 — "the cost of that hinge, not of reaching any alternative bound state" |
| **D63** | the forced 1.94 Å "never propagated"; D59/D62 overstated the constraint | D64, D67, D68. **This is why the clean re-do (job 30205872) was cancelled** |
| **D66** | "nobody chose this" pose; pose 2 was inherited | D68 — chosen deliberately, `model: 2` in all three configs |
| **D100** | the restraint mechanism "works" — Gln93 stayed out, no water intruded | D102 — 100 ns of the same scheme goes six-coordinate; "the 1 ns test was too short to rule it out" |
| **D102** | the attack angle "rotated into the window. It was a sampling question" | D104 — "the near-attack geometry did not survive, and this reverses D102": 101° → 57° on release |
| **D104** | GFN2 removed the classical over-coordination, six ligands to five, "confirmed rather than assumed" | D120 — both arms five-coordinate in 100% of frames against the crystal's four; "neither the QM/MM MD nor GFN2 created this; both inherited it" |
| **D104/D106/D110** | GFN2 will not hand us a near-attack frame, so `scan_attack` is mandatory | D120, D122 (explained by the saturated shell and a pose with no Zn–carbonyl contact, "without any appeal to GFN2"); D126 — 3.17 Å at 106° reached with **no attack restraint at all** |
| **D105** | the fix "needs no new setting: shrink `actradius`" | D106 — "the lesson is the opposite": 5 Å diverges (−240 → +1139 Eh); the 8 Å run was slow but healthy |
| **D111** | "two things are rock solid" — the hydroxide stays on the zinc 100%, Asp92 OD1 holds it 100% | D131 — "measured with Asp92 restrained against the metal… has to be re-established on a correct shell before it can be cited" |
| **D119** | M1 answered: the protonation state "is not the variable that decides this" | D120 — "moving one proton cannot open a coordination site that a restrained carboxylate is occupying, so both arms were always going to lose the substrate" |
| **D120** | "the restraint whose purpose was to preserve the crystal's coordination sphere is what destroyed it" | D133 — "D102's two newcomers were unrestrained, so a restraint is not the mechanism". Drove a cancellation, a cutoff change and the whole `CSR_wt` rebuild |
| **D121** | acceptance criterion 1: a four-coordinate Zn1 with carbonyl O2x at 2.0–2.6 Å, "acceptance criteria, not preferences" | D131 — Zn1 is H-R-H with two coordinating residues, so the canonical-equivalent sphere is two His + hydroxide + one water. Encoded as a hard gate that has already rejected a valid structure |
| **D122** | the docked pose "has exactly one defect", the carbonyl | D123 — the ligand is bound the other way round, O1x facing the metal; D126 — no pose puts any atom within 5 Å of the zinc |
| **D123** | pose 2 puts the carbonyl **2.55 Å** from the zinc | D126 — FALSE: it is the May run's own restraint output, read from a file D123 itself calls unreusable; pose 2's closest heavy atom is 7.23 Å. *Now marked at the headline* |
| **D124** | pose 2's 2.55 Å is inside the loose-coordination range thermolysin/CPA establish | D126 — pairs a real literature result with a false measurement |
| **D125** | the criteria are "a command, not a habit"; "the closest substrate oxygen must be the β-lactam carbonyl"; regression test "pose 2 at 2.55 Å must be accepted" | D127 — that is a criterion for the hypothesis *under test*, and the gate rejected the best attack geometry in the project; D126 — the accept case is a false distance frozen into the suite |
| **D128** | `md_subst` passed `site.py check`, submitted as job 30923516 | D130/D131 (cancelled; restrains Asp92 harder than the run D120 condemned), D135 — "`md_subst` is DO NOT RUN" |
| **D132** | `md_wt100`'s zinc carries the crystal's four ligands and nothing else — "the first time a run of ours has had the right sphere" | D134 — Asp92 bound in 100% of frames, 1.91 Å at the first written frame, CN mean 5.57, four-coordinate never. The sphere lasted under 20 ps |
| **D133** | the standing xTB over-coordination question is answered | D59/D206, and its own next paragraph. *Now marked* |
| **D134** | removing the restraints "produces a worse-coordinated" zinc — "the question D132 was built to answer, and it is answered"; "rule 3 now has a restraint-free number behind it" | D144 — every MD-derived statement about Asp92 coordination is "inadmissible, not merely caveated"; and a non-bonded point-charge MM zinc is the method this log elsewhere says cannot answer coordination |
| **D136** | `scan_aspout` will decide whether a relaxed scan reports the basin it started in | its start is `naclong`'s frame, which D120/D122/D130 all call the condemned shell and D131 says nothing measured on it can be cited |
| **D137** | "a relaxed scan **cannot converge**" with `$set` rigid waters | D141 — "the sweeping version of D137's claim is wrong": 67 points across four scans converged. D142 narrows it to MD-frame starts |
| **D138** | Asp92 bidentate in every frame, six-coordinate throughout; D120's artefact diagnosis "misattributed" | D144 — inadmissible |
| **D143/D144** | "the crystal's 2.80 Å sits 5.4–7.9 kcal/mol up this curve — the honest figure to quote"; restrain at 2.80 | D145 — the crystal geometry "is off the end of both arms, and the cost of reaching it has not been measured at all". §1(a) |
| **D145** | the crystal has Asp92 4.20 Å from the catalytic zinc | D151 — measured in a mutant that **restores Zn2**. §1(a) |
| **D147** | Asp92's position is "worth about 9 kcal/mol", settling D144's question 2 | D149 — "the trend D147 reported does not survive the last two rows": non-monotonic, peaking at 2.75 Å |
| **D153** | Arg229's QM–MM electrostatics is +44.34 kcal/mol | D156/D157 — computed on a structure whose ligand is 5.15 Å from the only experimental Michaelis complex |
| **D157** | Gln93's contribution is −0.87 kcal/mol, "the last word on the Q of RQK"; the decomposition is "complete" | its own §4 — "the deletions and the PMF are both instruments pointed at a complex whose ligand is in the wrong place" |
| **D170** | Gln93 "is in neither QM region" — presented as a defect to fix | D171 — the omission was the user's deliberate decision; Gln93 was removed again |
| **D173** | Gln93 is 7–9 Å from everything, "which settles it" (neutral, remote, non-functional); and the NaN is explained by the rearranged site | D191 — Gln93 OE1 on the zinc in 3.0% of frames, min 1.98 Å; D195 — 1.79 Å in `CSR_unrestrained`. D191 — `prod100c` NaN'd again *with* the site restrained, so the rearrangement does not account for it |
| **D175** | the per-atom option "existed from the day it was split out and was simply never passed" | D176 — "No: detection was never per atom for restraints. Per residue, always" |
| **D184** | 6DN4 puts the zinc 2.68 Å from Asp120 OD2; all six site distances agree within 0.13 Å, so "the starting structure is right" | 6DN4 has **no metal** (D179), so 2.68 Å is our own superposition; and "right" holds for one side chain — D192/D195 condemn the pose and the sphere |
| **D185** | both runs are live and "both restraint lists are exactly right" | D194 — the attack restraint `[C:2:O01, C:1:C4x, 3.47]` was **replaced** rather than added to |
| **D185** | "the zinc stays four-coordinate" cannot be the test — a gate written three times is wrong | D186 — "the four-coordinate gate was never wrong; measuring it at the restraint radius was" |
| **D186** | at the 2.5 Å radius the site is **four**-coordinate, so the gate stands as it was; and bulk water is not a gate | D191 — CN 6 in 3523/5000 reactant frames; and E:2920/E:3189 are resident 97.0%/89.3%, "not coming and going" |
| **D189** | the HE2-on-every-histidine problem "never reached a production system — in every single run" | its own text admits six run directories were unaudited; "every single run" is five structures read |
| **D191** | Lys226 NZ–O01 at 2.72 Å in **100% of 5000 frames** — the converged RQK result | D197 — inside 3.5 Å in only **10.3%** of frames, mean 4.93 Å. Measured in `wt100c`, the run condemned in the same entry |
| **D194** | the cause is one line of config, the replaced attack restraint — "that is the whole of it" | D197 — with the restraint restored in `eq5`, every restraint held and the structure is **still** rejected; the pose's own O1x took the fourth site |
| **D195** | GFN2 is not overbinding: 1.94 → 2.56 Å, "inside 0.07 Å of 7UHI/7UHJ" | D197 — "a bond-length relaxation within the coordinated basin, not a rejection". n = 1, our own run. *Now marked* |
| **D196** | the first build whose resting sphere **is** the Michaelis sphere | D197 — over 10,000 frames the carbonyl is inside 2.5 Å in **0.3%**; the fourth site is held by the unrestrained hydroxyethyl O1x at **97.9%** |
| **D197** | "PASS — Zn coordination number 4", the gate cleared for the first time | same entry, 120 lines later: "the number is correct and its membership is wrong" — the four are hydroxide + 2 His + O1x, carbonyl absent |
| **D197** | the 2.60 Å carbonyl target and the 2.5 Å radius expose a conflict between two settled numbers | its own correction block — "the carbonyl at 2.60 Å is not in conflict with anything. It is D124's design criterion being met" |
| **D198** | "GFN2 does **not** over-coordinate the zinc", answering the blocker open since D9; frame 4374 is the only qualifying structure | its own caveat limits it to one contact, one geometry pair, single points; D200 — batch 1's last frame also passes |
| **D200** | batch 1's zinc site is the designed sphere; a Bürgi–Dunitz angle needs the ligands essentially trans | D201 §4 ("D200 said it was, and that is wrong") and D201 §2 (1TMN and 7CPA are 57–60° **cis**) — 240 lines later |
| **D204** | the reactant is Zn–OH₂ with neutral Lys226 "in a correctly built Michaelis complex and not only in a degraded one" | its own closing section (Zn–Asp92 at 1.78 Å) and D206, which rejects that geometry on the structural criterion |
| **D225** | "**No 100 ns trajectory survives**, so the 10-vs-100 comparison cannot be made directly" | **D226.** True of the DCD, false as written. `wt100c` finished its full 100 ns and its `NVTtrajectory_lastframe.pdb`, `NVTtrajectory_firstframe.pdb`, NPT frames and 5000-row `nvtsim.csv` are all on disk — in `CSR/CSR_wt/`, not `CSR_mich/`, which is why one directory's listing was reported as a general absence. The comparison was made from them. **Its conclusion — 10 ns is enough — survives and is now measured rather than inferred.** |
| **D225 (box drift)** | "box volume drift **+5.396** nm³ ... eight times larger" marks `asppull3` as still relaxing, and *"the box is still accommodating the 1.3–2.2 Å displacement of the Asp92/Lys226 pair"* | **D231.** The statistic is last-decile-minus-first-decile, and NPT box volume steps once then plateaus, so it measures how far the box had to travel, not convergence. All three runs plateau in decile 1 and reach the same ~480.7 nm³; by final-three-decile spread `asppull3` (0.51) beats `eq6` (1.07). The displacement needed no accommodation. **D225's other reason — RMSD still climbing over 2 ns — stands, and is sufficient.** |
| **D216** (mechanism) | *"deprotonating Lys226 removes what holds Asp92 off the metal, so Asp92 moving onto Zn1 is **downstream of the Lys226 proton transfer**"*, offered as the testable explanation of the collapse | **D232.** The test ran. Over 32 ps of unrestrained GFN2 from the crystal geometry, **all three Lys226 NZ hydrogens stay on the nitrogen** (0.98-1.09 A), the hydroxide keeps exactly one H, the salt bridge holds at 2.87 A mean - and Asp92 chelates the zinc anyway, bidentate at 2.10/2.29 A from 8 ps on. The transfer is not a precondition for the collapse. The measurement that D216 was explaining stands; the explanation does not. |

---

## 3. Overturned, but not in the headline

| entry | claim | overturned by |
|---|---|---|
| **D236** (mechanism/method) | *"the DFT gradient at the crystal-consistent geometry does not pull Asp92 onto the metal ... the first direct evidence that GFN2's chelate is the Hamiltonian's and not the enzyme's"* | **D240, six hours later.** The matched GFN2 optimisation - same frame, same QM region, same active region, same optimiser, only `qm_engine` different - **moves Asp92 0.23 A further out where DFT moved it 0.07**, and tightens the same salt bridge to the same place. So `aspdft`'s three steps say nothing about DFT versus GFN2; what survives is that *neither* gradient points the carboxylate at the metal from this start. The numbers in D236 stand; the inference does not. `gfn2-over-coordinates-zn1-from-a-crystal-consistent-start` is narrowed accordingly, as D238 pre-registered it would have to be. |
| D204 | Lys226 hands its proton to the zinc-bound oxygen "within 2 ps of becoming quantum", stated as a property of the enzyme | **D232.** Same QM region, same charge, 32 ps from a start with the carboxylate on the lysine instead of the metal: **no transfer at all**. Real, and conditional on the starting geometry rather than general - which matters because `nucleophile-is-zn-bound-water-with-neutral-lys226` rests on it. |
| D12, D13 | MBLs bind the C3 carboxylate with a Lys, and Lys226 is it | D107 — it is BBL Lys263, the Zn2 position *(marked)* |
| D32 | Arg90 NH1···O1x at 99% occupancy | D148 — 0% over 5,000 frames; pose-dependent |
| D44 | tier B: every published number re-optimised at DFT/MM | D47 — "tier B (D44) is reversed" |
| D40 | no crystallographic Michaelis complex shows a firm Zn–O(carbonyl) bond | D121/D123 — false for B3; 7UHI/7UHJ at 2.53–2.63 Å. *Marked in D40 — but D53 repeats it unmarked* |
| D49 | no water both attacks and donates to Asp92 | self-caveated in its text; the headline is absolute |
| D74–D92 | the dissociation is a pose problem; the pose never made the canonical carboxylate–Arg contact — "the next scientific question" | later entries call that criterion circular and MM dissociation "uninformative by construction". Three weeks of pose-hunting rest on it |
| D97 | the 1.79 Å Zn–Gln93 contact sits ~5.6 kcal/mol up a repulsive wall | measured on "a broken model — run 1's MM end state, carrying the 1.79 Å artefact", with pre-merge zinc radii |
| D107, D112 | "a 2+ cation 4–5 Å away should **raise** a lysine's pKa" | D206 — sign error; the conjugate base is the neutral amine, so a cation **lowers** it. The bookkeeping is indeterminate, not evidence for `lysH` |
| D114 | a water does not fit where the hydroxide sits, so the geometry is available only to a hydroxide | the frame is `md_naclong`'s condemned restrained shell; D119 found the two arms indistinguishable anyway |
| D131 | the `SEQADV` records **prove** Zn2 is vacant in the wild type | they establish the deposited construct is a double mutant; vacancy in the wild type is an inference from that — and the same records are why the 2.80 Å reference may be a mutant artefact |
| D147 | `substituent_on_metal` at 3.17 Å / 106° "with nothing restrained" is the surviving hypothesis | D148 — that structure has Lys226 NZ **0.58 Å** from substrate O4x; "its structure is refuted" |
| D150, D155, D158 | Arg229+Lys226 are "a fixed clamp", 0.08 Å over 45 points — "the first thing that behaves like a real anchor" | D156 — all 45 points descend from pose 2, 5.15 Å from the experimental mode. The rigidity is a property of one docked pose |
| D151, D153 | Asp92 might be protonated, since the di-zinc pKa argument "does not transfer" | D154 — at pH 7 that carboxylate is deprotonated either way; the Zn1 contact itself lowers the pKa |
| D157 §3 | Asp92 stays on the metal at a third level of theory, confirming D102/D134/D138 | D144 — those three are inadmissible, so they cannot corroborate; 3 ps on the condemned pose |
| D159, D165, D166 | the product arm's start has "crystal geometry", Asp92 OD2 out at 2.80 Å | written after D145 retired that number. §1(a) |
| D171 | 34.14 ns answers the zinc-shell question as well as 100 ns would | D173 — the 34 ns "did not equilibrate the complex, it rearranged it" |
| D174 | the 3.0 Å cutoff caused the five-coordinate zinc | D175/D176 — a per-residue default in `find_coordinators`; the same wrong blame as D120 and D130 |
| D175 | `prod100c` must be resubmitted, and the cost is GPU hours already spent | D176 — the job was still PENDING and used no GPU time |
| D178 | 7UHI becomes "the check, not the source" | D192 — "that makes 7UHI the source for the Michaelis geometry, not merely the check" |
| D179 | two independent lines say the carboxylate is mobile | one of the two is our own PMF; and the entry itself says the decisive Zn–Asp number "is geometry, not measurement" |
| D181 | `constraint_exclude_residues` stays — "the substrate must be free to react", and it "costs nothing either way" | D193 — "a QM/MM concern applied to an MM stage", and "costs nothing either way was circular" |
| D187, D188 | the cluster is net −1 with a substrate fragment, ~15 atoms to the substrate | D198 — both clusters are 51 atoms, C15H25N5O5Zn, charge **0** |
| D196 | removing `D:1:HOH` is "confirmed by measurement" to leave the five-donor sphere | D197 — the vacated site was taken by O1x at 97.9%. The "measurement" is a frame D195 accepts while recording CN 6 and a 26° attack angle |
| D198 | the site holds a fifth donor in MM whatever we do; Asn225 took O1x's place | D200 — Asn225 leaves (2.19 → 4.18 Å closest, 5.86 Å mean) |
| D200 | hysteresis narrowed: an anionic or chelating donor on the metal does not leave in 50 ps | D201 — Asp92 is anionic and crossed **inward** in that same run |
| D201 §4 | "Asp92 on the metal is a real feature of this model" *(bolded)* | D206 rule 3. This sentence justified starting `michrk` from a 2.00 Å Asp92 frame and ruling out CphA a third time |
| D204 | the coordination sphere corroborates the protonation state "independently" | both observables come from the same GFN2 trajectory on the same geometry — not independent evidence classes |

---

## 4. Claims stronger than their evidence (not overturned; do not quote as settled)

- **D197** — "no GFN2 benchmark for zinc coordination number in enzymes appears
  to exist": a literature-absence claim with no recorded search.
- **D199** — the O1x contact being a 12-6 artefact makes D197/D198's conclusions
  "stronger, not weaker": an inference that immunises them against D198's own
  caveat.
- **D201 §2** — "the project has been rejecting its own structures on a criterion
  both reference enzymes also fail": 1TMN/7CPA are measured, but the
  angle-versus-distance curve is our own extrapolation to distances those
  crystals do not occupy.
- **D206** — promotes D48 to "best evidence" while its own rule 3 excludes that
  evidence class. Now distinguished in the entry: D48 was *handed* its geometry
  from the crystal and asked only about protons, a legitimate conditional;
  D58's was loaded by a geometry no structure supports.
- **Central hypothesis block** — Arg90 at 4.41 Å and Lys226 at 8.73 Å are called
  "two measured caveats" in the document's most-read section. They are one
  optimised frame of our own model, from a condemned lineage.
- **Cross-cutting, D170–D197**: every structural conclusion in that range except
  the crystal comparisons rests on our own MD, optimisations or single points.
  The only non-model evidence quoted anywhere in it is 6DN4/6DQ2/6NC5,
  7UHI/7UHJ's 2.53/2.63 Å, GOB-18's four-to-five expansion,
  thermolysin/carboxypeptidase A, and Imhof's 7.55 kcal/mol.

---

## 5. Already marked inline

Not grep traps; listed for completeness. D12, D13, D40 (carbonyl bullet), D44,
D58 (body), D105 (count claim), D123 (body), D124 (body), D137 (option 2), D151
(umbrella cost), D153 (pKa sign), D154 (`CSR_hydrolysed`), D155 (240 days), D159
(items 1 and 3), D165 (withdrawn section, **kept verbatim** — a grep for
"restrained at 2.80 Å" still hits it), D170 (Gln93), D178 (item 0), D181
(`md` resume), D191 (provisional), D192 (diagnosis withdrawn), D193 (causal
claim withdrawn), D196 (timestep), D197 (two corrections), D207 (its own
qualification).
