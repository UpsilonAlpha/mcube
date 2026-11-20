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

**One restraint is WRITTEN here; it is not the only one applied.** `md_restraints`
builds "the detected metal coordination, plus any asked for" - so
`preprocess.results.coordination`, every donor inside the 3.0 Å radius, is
restrained automatically on every MM leg, and the `restraints:` list in a section
is the *additions*. Zn-O2x and Zn-O01 are in that detected set. They are held.

So the earlier claim here - that this leg leaves the carbonyl and the hydroxide
free, and that no MD section has ever restrained Zn-O01 - was **wrong, and is
retracted (¶293)**. It read a config list as if it were the whole restraint set.
`relax10`'s hand-written `[B:1:ZN, C:1:O2x, 2.60, 500.0]` was an override of a
detected row, which is why it looks like the only place the carbonyl appears.

**The live question is different, and it is `coordination_redetect`.** This leg
sets it `false`, so the detected restraints come from the ORIGINAL structure and
hold Zn-O2x and Zn-O01 at their preprocess distances - while the leg starts from
Leg A's endpoint, where Zn-O2x has gone to ~2.18 Å. Holding the sphere at
distances the structure no longer has is restraining it to a different pose.
`true` would re-measure on `scan_10.pdb` with the same metals, cutoff and donor
elements, which is the documented route for an edited frame - at the cost that
the run is then not numerically comparable with the ones that used the recorded
values. That choice is the user's.

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
