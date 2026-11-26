---
type: observation
status: supported
evidence: ours
entries: [¶1, ¶17, ¶144, ¶207, ¶308, ¶339, ¶341]
checked: [¶144, ¶341]
---

# No reaction coordinate has ever been computed in this project

**The claim still holds and all three of its original grounds are now
obsolete - updated 2026-11-26.** There is still no barrier, no transition state,
no IRC and no intermediate, but not for the reasons first given: the scan HAS
run (¶308, 16 points, all converged, so `scan_attack/scan_15.pdb` exists), and
`ts` HAS run - twice, 33410233 OOM-killed at 32G and 33419978 at 128G failing on
`tsopt_maxiter: 100`. What is now in the way is different and larger.

**The NEB was abandoned** after ¶335 measured it diverging rather than
converging, and ¶334 found no MBL mechanism paper uses NEB at all. CI-NEB had
been in the plan since ¶1 and survived every revision; it is out on measurement,
and the sections are kept in `site.yaml` so the decision is reversible (¶338).

**And the `ts` route's second failure was a result, not a resource problem.**
`saddle1` held the reaction coordinate still and drained 18.41 kcal/mol of MM
strain out of the geometry the scan handed it
([[our-stage-energies-are-not-on-a-common-surface]], ¶341). So the obstacle is
no longer "no path method has completed" but **"no two stages here produce
numbers that can be subtracted"**, which ¶332 has to settle first.

**The diagnosis worth keeping** (¶207): each blocker was investigated on its
own terms rather than as a question about the endpoint. Asp92's coordination
got three campaigns, the protonation state got three, and the loop could not
terminate because every attempt to settle the reactant used our own minima as
the evidence - which ¶206 rules inadmissible for exactly that. **A question
whose only admissible evidence class cannot settle it will reopen
indefinitely**, and it did.

**And a correction to how the number should be framed.** NEB was our
substitution for Imhof's path network, not the field's method. No MBL mechanism
study in this project's survey uses NEB; CphA's ~14 kcal/mol - the mono-zinc
number we would compare against - is a **PMF** barrier, and ¶32 was a user
decision to replicate that. So a minimum-energy barrier quoted beside it
compares two different quantities (¶207, as qualified by the user the same day).
¶334 later evidenced the first half of this and corrected the second: about a
third of MBL mechanism papers DO report potential-energy barriers, so a
potential-energy barrier is publishable in this field - it is NEB that is
unprecedented, not the quantity. **¶334's further claim that Khrenova & Nemukhin
2018 and Levina 2020 used our exact protocol is unverified and withdrawn**
(¶341): Levina defers its whole protocol to a closed paper and states no
equilibration, restraints, frozen atoms or active radius.

## Supports
- [[thesis-mechanism-differs-from-canonical-b3]]

## Supported by
- [[our-stage-energies-are-not-on-a-common-surface]]

## Refutes
- [[rqk-cations-replace-zn2-electrostatically]]
