# Equilibrating the carboxylate-on-metal arrangement: the plan

Written ¶286 (2026-11-15), **not submitted**. The user asked for this to be
planned and will say when the cluster is open. Nothing here is to be run
without that.

## What the run is for

To find out whether the arrangement every geometric argument since ¶283 rests
on - substrate **carboxylate** on the metal, **carbonyl** off it and held by
Arg90, nucleophile in the freed position - survives contact with solvent,
dynamics and QM/MM. At the moment it is a set of distances carried in from
another structure and never relaxed.

## The problem the plan has to solve first

The arrangement comes from 7UHI's ligand **MX0**, which is an **oxacephem**
(C13 H16 N2 O8, no sulfur, 23 heavy atoms) and **not meropenem** (C17 H25 N3 O5
S, 51 atoms). ¶156's transfer renamed it `mer` and the name stuck. So those
coordinates **cannot be given to our force field** - the molecule is not ours and
is not parameterised.

Three ways to reach the arrangement with meropenem, and the third is the one
this plan takes:

| route | what it costs | why not / why yes |
| --- | --- | --- |
| A. Equilibrate MX0 itself | new ligand parameters; a `preprocess` from scratch | Answers a question about a **cephem**, not our substrate. Interesting for selectivity later, not this run. |
| B. Map meropenem onto MX0's shared atoms, then equilibrate | a build; a 3-4 atom fit like ¶156's | It is "an arbitrarily shifted structure" by another name, and the user has twice ruled that out as a starting point. |
| **C. Restrained pull from the existing equilibrated system** | **nothing new** - no parameters, no `preprocess`, no build | Uses our own parameterised, solvated, QM/MM-equilibrated system and moves the ligand along a **physical path**. Starts from an equilibrated structure, as required. |

## The coordinates involved

Measured in `aspxtb8_nowat.pdb`:

    O2x  on C4x    the beta-lactam carbonyl       Zn-O2x  **2.21**  (ON the metal)
    O3x  on C6x    the C3 carboxylate             Zn-O3x    6.70
    O4x  on C6x    the C3 carboxylate             Zn-O4x  **5.42**  (nearer of the pair)
    O1x  on C2x    the other ring carbonyl        Zn-O1x    4.20
    O5x  on C11x   the side-chain amide carbonyl  Zn-O5x   10.27

So the move is: bring **O4x** from 5.42 Å in to about **2.2 Å**, and let **O2x**
leave on its own.

**The target is 2.2 Å, not 1.847.** 1.847 is the transferred value and is a
transfer artefact to within ~1 Å (a global fit gives 2.75). 2.2 Å is 5B15's
**measured** Zn2-carboxylate distance (OAE 2.23); 7UHI's is 2.65. Using a
measured number rather than our own placed one is the point.

## Leg A - `posepull`: pull the carboxylate in, and leave the carbonyl alone

    stage: scan

A `scan` on the single coordinate `B:1:ZN - C:1:O4x`, from 5.4 Å to 2.2 Å.
**Only the carboxylate is driven.** The carbonyl is left completely free.

This is deliberate and is the most informative part of the plan. If the
arrangement is real, the carbonyl should leave the metal by itself as the
carboxylate arrives, because a two-histidine zinc cannot hold both plus a
nucleophile. **If it does not leave, the arrangement is not accessible to our
system**, and that is a result rather than a failed run. Driving both
coordinates would assume the answer.

Expected on success: Zn-O2x relaxes out toward ~3.2 Å (the template's 3.153)
with no restraint telling it to.

## Leg B - `poseeq`: hold it and let everything relax around it

    stage: md

Minimise, NPT, then NVT, following `relax10`'s pattern exactly - `gentle:
false`, `coordination_redetect: false`, CUDA, mixed precision, 10 ns.

**THE RESTRAINT SET AS WRITTEN IS INVALID AND MUST BE SETTLED BEFORE THIS LEG
RUNS (user, 2026-11-20: "we don't use classical MD for coordination - that's a
hard rule").**

As written, `poseeq` restrains one coordinate bond - `[B:1:ZN, C:1:O4x, 2.2,
500.0]` - and says the carbonyl and the hydroxide are "free". But in Leg A's
endpoint both of those are still coordinating: Zn-O2x at ~2.18 Å and Zn-O01 at
~1.75 Å, each well inside the 3.0 Å rule. Leaving them unrestrained hands them
to MM, which is precisely what MM cannot represent, and it hands MM the very
question Leg C exists to decide. `relax10`, which this leg claims to follow,
did restrain Zn-O2x (2.60) - dropping that restraint was not a documented
choice, it was an oversight.

A second gap, pre-existing and not introduced here: **no MD section in this
project has ever restrained Zn-O01**, the zinc-hydroxide bond. `relax10`
restrains the O01-C4x attack distance but not the bond to the metal. That needs
the user's decision, not a default.

So this leg is suspended pending the restraint list. It cannot be written in any
case until `scan_10.pdb` exists, so nothing is blocked by waiting. What the leg
may legitimately do once settled is **re-equilibrate the solvent around a ligand
that has moved**, with the whole first coordination sphere held at Leg A's
geometry. What it may not do is decide anything about coordination - including
whether the pose "survives", which was how it was described to the user and was
wrong.

This is the "solvate and equilibrate" leg. Note the system **is** already
solvated; what this does is re-equilibrate the water around a ligand that has
moved, which is the part that has never been done.

## Leg C - `poseqm`: QM/MM with no restraints at all - the actual test

    stage: qmmm

Modelled on `aspqm2`, but with **`qm_residues: *qmreg_michrks`** - the region
including **Ser186**, approved by the user ¶287 and now in `site.yaml`.
`charge: 2`,
`skip_optimize: true`, `skip_equilibrate: false`, `rigid_water: false`,
`time: 50`, chunked by raising `time` to 100, 150, 200 one at a time.

    require:
      unset: [restraints, opt_constraints, constraints]

**No restraints, by design and by precedent.** The MM legs restrain the
coordinate bond because MM cannot represent it; QM/MM can, so holding it here
would measure the restraint instead of the chemistry - which is why `aspqm2`
refuses them outright.

**This leg is the experiment.** If the carboxylate stays on the metal through 50
ps of unrestrained QM/MM and the carbonyl stays off it, the arrangement is
viable and every argument built on it stands. If the ligand reverts to the
`aspxtb8` arrangement - carbonyl on, carboxylate off - then the carboxylate-on-
metal reading is an artefact of the transfer and a large part of ¶283-¶285
goes with it. Either way it is decisive, which no amount of further geometry
will be.

Watch, in order: Zn-O4x, Zn-O2x, Zn-OH, the attack distance OH...C4x, the
out-of-plane angle, Arg90 NH1-O2x, Ser186 OG to O3x/O4x, and the metal's
coordination number.

## Decisions needed before this is submitted

1. ~~Ser186 in the QM region?~~ **SETTLED (user, ¶287): yes.** `michaelis_rks`
   is written into `site.yaml` - `michaelis_rk` plus `A: {186: SER}` - and the
   **charge is unchanged at 2**, serine being neutral. Leg C uses it.
   `michaelis_rk` is kept alongside, because `aspqm2`'s 300 ps is the
   comparison baseline and changing its region would break the comparison.
2. **The starting structure.** Preferred is the end of the QM/MM leg
   (`qmmm_aspqm2_results`), to honour "start from a QM/MM equilibrated
   structure"; fallback is `md_relax10_results/NVTtrajectory_lastframe.pdb`.
   **Which of these exists on Bunya needs checking before anything is written** -
   `state.md` is not reliable here, since the `prodpull*` sections it reports as
   run were deleted in ¶262.
3. **Nothing else.** No `preprocess`, no new parameters, no new ligand, no
   builder, no change to any default. Leg A's scan is the only non-default
   choice in the plan and `scan` is the stage that exists for it.

## What this plan does not do

- It does not test **hydroxide against water**. That needs a different
  experiment, and 50 ps of dynamics cannot decide a protonation state it was
  handed - though running Legs B and C twice, once with each, and comparing
  Zn-O bond length and the hydrogen-bond pattern against the crystals would
  discriminate indirectly. Worth doing **after** the pose is known to hold, not
  before.
- It does not address the **carbapenem selectivity**, which needs kinetics.
- It does not verify the **metal's position**, which is the thing that differs
  from every crystal and which only a structure with M1 occupied and M2 empty
  could settle.
