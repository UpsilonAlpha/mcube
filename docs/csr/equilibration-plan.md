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

## Leg B - DROPPED (D295)

**The user, 2026-11-20: "I see no reason to run B if we're doing the scan from an
already optimised starting point. There's no reason to drop down to a lower
level of theory."** Correct. There are two legs.

The scan starts from `aspxtb8`'s converged QM/MM minimum and runs at QM/MM
throughout, so an MM leg between it and Leg C re-describes a QM-derived geometry
with a force field that has no parameters for the species in question - and 10 ns
is long enough to degrade it. The pipeline order is preprocess -> md -> qmmm for
a reason; going back up it is not a step.

**My argument for keeping it (D294) was wrong in its own terms.** The evidence I
cited - `aspxtb8` (D244), where freeing the second shell from 124 atoms to 479
moved the carboxylate rotamer by up to 1.5 Å - is evidence for a **larger active
region in QM/MM**. I used it to argue for an MM leg, which is not what it
measured.

**And the leg could not have been made safe.** `coordination_redetect: true`,
which the leg needed to be valid at all, would re-measure coordination on
`scan_10.pdb` - whose Zn-O4x distance of 2.2 Å was **imposed** by the scan's
driving restraint, a number chosen from 5B15 rather than measured here.
Redetection would enter it as *detected* coordination: D68's laundering
mechanism exactly. And the D71 guard does not reach it, because `preprocess`
compares detection against `md.results.restraints` only, so a scan's driven
coordinate is invisible to it.

**Residual approximation, stated rather than fixed:** nothing beyond the active
region relaxes around the moved ligand. That is the same frozen-shell
approximation `aspdft` and `aspxtb` already run under, and the instrument for it
is `actradius` in Leg C, not a force field.

**And the redetect "exception" written here on 2026-11-20 is withdrawn (D296).**
Always redetect: detection is the check for intruders, which is this project's
largest invalidation (D177, a six-coordinate zinc in 4995 of 5000 frames), and a
restrained pair is the one most worth checking because a restraint holds a
distance without holding the site. The laundering refusal belongs in the D71
guard, which needs extending to see a `scan`'s driven coordinate.

## What Leg A is for, and what its energies are not (D297)

**The scan forces a structure for the geometric checks. That is its whole job**
(user, 2026-11-20: "the scan is to force a structure for geometric checks, a PMF
doesn't help with our main goal"). D296 proposed an umbrella PMF along Zn-O4x
and that proposal is **withdrawn**.

`scan` optimises the active region at each point with the driven distance fixed:
11 constrained minima at GFN2. Those minima are a means of arriving at an
endpoint, not a measurement. **Nothing claims the enzyme reaches this pose by
having the carboxylate dragged onto the zinc from an already-bound structure** -
the substrate binds from solution. So the driven coordinate is an artefact of how
we build the structure, and a free energy along it answers a question nobody
asked. A PMF there would be rigour spent on the wrong quantity.

**So the +22.1 kcal/mol must not be cited as a result,** and D291 and D294 both
quoted it as "the cost". It is the work a restrained optimiser did along an
artificial coordinate, with the region beyond 8 Å frozen. Report it, if at all,
as evidence that the restraint had to push - never as a barrier, a binding
penalty, or an accessibility figure.

**What Leg A legitimately yields:**

1. `scan_10.pdb`, the endpoint structure, for the reactant gates - **`O01 - Zn -
   O2x` 45-95 deg** and **`O01 ... O2x` 2.0-3.3 Å**, both `required_at:
   [reactant]`, both currently `unsupported`, and both able to **reject** the
   structure (`gates.md`).
2. The monitored distances along the path, which are geometry and need no energy
   to read: the carbonyl staying on the metal at 2.18 Å and the attack distance
   opening from 3.457 to 3.977 Å.

**Then Leg C decides**, by unrestrained QM/MM dynamics, which is the only leg
that decides anything. Two legs, one gate between them.

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
