---
type: claim
status: supported
entries: [D156, D158, D159, D178, D191, D192, D194, D195, D197, D198, D200, D201, D207, D214, D231, D246, D247]
checked: [D158, D159, D178, D191, D192, D194, D195, D197, D198, D200, D201, D246, D247]
---

# No usable Michaelis complex has been built, and that is what blocks the mechanism

Eighty-four candidate last frames have been checked with `site.py check`; the
criteria have never all been met in one structure that also survives scrutiny
of its membership.

The record, compressed:

- **`wt100c`, 100 ns:** the substrate **leaves**. O01-C4x inside 3.5 A in 0.2%
  of 5000 frames, mean 19.95 A; confirmed real rather than a periodic-image
  artefact (D191/D192). Cause: one config line, an attack restraint replaced
  rather than added to (D194).
- **`mich_eq5`, 10 ns:** every restraint held to better than 0.11 A and the
  structure is still rejected - the pose's own 1-hydroxyethyl O1x took the
  fourth coordination site at 97.9% occupancy, unrestrained, while the carbonyl
  was inside 2.5 A in 0.3% of frames (D197). `site.py` reported "coordination
  number 4" correctly with the wrong membership.
- **Frame 4374 + QM/MM batch 1:** O1x is an MM artefact and leaves at the GFN2
  level (10% against 97.9%), Asn225 leaves, and the carbonyl coordinates
  unrestrained - but Asp92 bound the metal at 2.00 A in 96% of frames, so the
  sphere is not the designed one (D200 as corrected by D201).

- **`eq6`, 10 ns, 2026-11-06 (D214):** the first build to pass the
  equilibration gate. Every configured restraint held to within 0.03 A, Zn1 is
  **four-coordinate at 2.50 A** (two histidines, hydroxide, water) with the
  carbonyl just outside at 2.57, and the Lys226-Asp92 salt bridge is held at
  2.79 A for the first time in the project. **What it does not settle:** Asp92
  was restrained at both ends, so the run is a valid starting structure rather
  than evidence about the carboxylate; and the underlying *pose* is still the
  open question below.

- **`relax10`, 10 ns, 2026-11-06 (D231):** the first build whose Asp92 is
  where the crystals put it *and* which holds there. Every restraint within
  0.09 A, and Zn-OD1 at 4.33 against 6DQH's 4.28 with nothing restraining it -
  the carboxylate choosing the crystal's rotamer rather than being told to
  ([[relax10-holds-the-crystal-geometry-for-10-ns]]). **What it does not
  settle:** the pose, which it inherits unchanged from the same docked complex
  as every run above; and the sphere, which is five-coordinate with two bulk
  waters where [[csr1-zn1-site-is-h-r-h]] expects one.

**The claim is now narrower than its heading.** An equilibrated complex exists.
What does not exist is a pose anchored to experiment. The docked complex puts
the carbonyl at 5.43 Å with one attack restraint, while 7UHI/7UHJ — the only
experimental B3 Michaelis complexes — put theirs at 2.53/2.63 Å (D123).

**THAT CHOICE WAS MADE, AND IT WAS NEITHER OF THEM (D178), THOUGH THE SCOPING
HAS SINCE MOVED.** The plan of record
builds the Michaelis complex from **5AXO** — SMB-1 with *hydrolysed meropenem*,
our own substrate in a B3 enzyme, measured — by **re-forming the C4x–N1x bond**,
reversing the hydrolysis:

    prod100c -> prodpost gates -> prodqm -> prodopt -> recombine
    (100 ns)    (shell + H)       (50 ps)   (opt)      (re-form C4x-N1x)

The structure `recombine` ends on **is** the Michaelis complex, anchored to
experiment through the product rather than through a docking score, and on the
product's own topology so it is a valid NEB partner by construction. **7UHI
becomes the check, not the source:** if the re-formed carbonyl lands at L1's
measured 2.53–2.63 Å the pose is independently corroborated; if it does not,
that is a result about CSR-1 rather than a modelling failure. D178 also
rejected the alternative this node used to point at — grafting meropenem's
substituents onto 7UHI's transferred moxalactam core — as "strictly harder and
less direct than re-forming a bond in a molecule that is already meropenem".

**D178's scoping of 7UHI was overturned: D192 made it the SOURCE for the
Michaelis geometry rather than the check**, on the ground that a wrong pose
cannot be repaired by restraining it harder. **D192's diagnosis was then
withdrawn by D194**, which is a consolidation entry whose purpose is to stop
one argument being re-derived: *"dissociation-implies-bad-pose"* was raised and
discarded at D89, again at D192, and D194 records it as discarded a third time.
Docked pose 2 with an **uncoordinated** carbonyl, held by one attack restraint,
is what the project deliberately chose to study (D98: *"Nothing asserts the
carbonyl coordinates the zinc"*), and D194 withdraws "transfer a different
pose" as a proposal to undo a decision rather than to fix a regression. **The
substrate loss in `wt100c`/`prod100c` was one config line** — an Asp92 OD1
restraint written *instead of* the attack restraint rather than beside it — not
the pose. **Anyone reaching for the pose to explain a ligand leaving should
stop and read D194 first.**

**Two things survive that argument's collapse, and they are what keeps this
node open.**

1. **D156's measurement, which owes nothing to dissociation.** Against 7UHI's
   transferred moxalactam, pose 2 is the **second-worst of nine** by centroid
   (3.66 Å), with its carbonyl 4.90 Å against L1's 2.27 and its C3 carboxylate
   8.14 Å against 2.82. D156: *"Had this comparison existed in October, pose 5
   or 8 would have been the candidate and pose 2 would have been rejected."*
2. **D197, which superseded D194.** With the attack restraint restored, every
   restraint held to better than 0.11 Å **and the structure was still
   rejected**: the pose's own 1-hydroxyethyl `O1x` took the fourth coordination
   site at 97.9% occupancy, unrestrained, while the carbonyl was inside 2.5 Å
   in 0.3% of frames. That is a pose defect established without the discarded
   argument.

**MEASURED 2026-11-09 (D247): the equilibration bent one coordinate and left
the body behind.** 7UHI's moxalactam transferred into `relax10`'s own frame,
fit RMSD **0.310 Å** — the best in this project's history of doing this,
against 0.622 into `CSR_wt` and 0.851 into `CSR_ligand`:

| quantity | `relax10` | 7UHI transferred |
| --- | --- | --- |
| ligand centroid to Zn1 | **7.28** | **4.77** |
| **centroid separation** | **6.20** | — |
| β-lactam carbonyl to Zn1 | 2.610 (restrained) | 3.15 |
| C3 carboxylate O to Zn1 | 5.656 / 6.813 | 1.85 / 3.98 |

The nine docked poses spanned 0.98–3.66 Å by centroid and `scan_00` gave 5.15,
which D156 called further than any of them. **`relax10` is 6.20 Å.** The
carbonyl reads 2.610 only because `eq6` restrained that one atom; the molecule's
heavy-atom centroid is 2.51 Å further from the metal than L1's, and the C3
carboxylate L1 anchors is 5.7–6.8 Å out. Frames and atom counts differ, so
tenths mean nothing — 2.5 Å does. **This does not revive
"dissociation-implies-bad-pose"**; nothing in it rests on a ligand leaving.

The transfer also puts L1's carboxylate oxygen **1.85 Å** from our zinc — a
clash — where `CSR_wt` gave 2.54, so this frame's metal sits differently with
respect to its own histidine pair and the transferred pose needs D156 §5.2's
clash check before it can be used.

**The 5AXO route is half built and stalled.** `CSR_prod` transferred 5AXO's LMP
into the CSR-1 frame, mapped it onto our atom names by graph isomorphism, and
established that **nothing of the product is within coordination range of
Zn1** — everything that binds, binds the Zn2 CSR-1 lacks — so the product arm
takes no substrate restraint (D159; that entry's separate claim of "crystal
geometry" with Asp92 OD2 at 2.80 Å is retracted, D145 having retired that
number). Then `prod100c` NaN'd at 22.68 ns and failed the gate, as did the
reactant arm (D191). Everything built since — `eq6`, `relax10`, `aspqm2`,
`aspdft`, `aspxtb`, `aspxtb8` — descends from pose 2, with `eq6` restraining
the carbonyl in to ~2.6 Å, which is the step that made it pass.

So "usable" is an assessment of our own structures by our own criteria, which
is why this node is held on `ours` evidence and flagged by the validator's
check 1. That flag is correct and is not resolvable by more computation of the
same kind. **`relax10` does not change this**: a startable structure is not a
correct one, and 10 ns of MM on an unresolved pose is exactly "more computation
of the same kind." Neither does the QM/MM leg: `aspxtb`/`aspxtb8` hold the
carbonyl at 2.21–2.28 Å unrestrained, but `eq6` put it there with a restraint,
so that is the pose being self-consistent at the QM level, not being chosen.

**Two structural lessons worth keeping.** Testing properties one at a time can
pass a composition that would be rejected if asked directly (D197). And an
extremum of a distribution describes the reactant state we chose, not an
equilibrium average - frame 4374 was selected, and that must be said wherever
the run is quoted (D198).

**This is the thesis's real blocker, and it is prior to the choice of
instrument** (D207): a path run from a wrong reactant "would spend its barrier
moving meropenem across the site, and the number at the top would not be a
chemical barrier at all. That failure would look entirely normal in the
output."

## Supports
- [[no-reaction-coordinate-has-ever-been-computed]]
