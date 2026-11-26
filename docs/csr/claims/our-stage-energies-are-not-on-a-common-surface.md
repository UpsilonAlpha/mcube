---
type: observation
status: supported
evidence: ours
entries: [¶333, ¶339, ¶341, ¶344]
checked: [¶341, ¶344]
---

# Energies from different stages of our pipeline are not on a common potential-energy surface

The MM environment relaxes between stages by far more than the chemistry costs,
so a number from one stage cannot be subtracted from a number from another.

**¶344 bounds this claim, and the bound is what makes a barrier recoverable.**
It holds between a FREELY OPTIMISED structure and anything else. It does **not**
hold within one `scan`, because the relaxation is far-field - 88% of the squared
displacement in the 8-12 Å shell - and the far field does not move with the
reaction, so its strain is the same at both ends of a scan and cancels in the
difference. See [[the-barrier-is-scored-along-the-scan]].

**Measured directly on `saddle1` (¶341).** The `ts` run started from
`scan_attack`'s point 11 and, over 101 frames:

- **held the reaction coordinate still** - O01-C4x 2.0304 → 2.0156 Å, range
  1.94-2.07, 0.015 Å start to finish;
- **dropped 20.65 kcal/mol, 18.41 of it in the MM term**, ending 6.38 kcal/mol
  below `scan_attack`'s point 0;
- did so **smoothly, with no discontinuity**: largest single MM step -6.75
  kcal/mol at the first step, five steps above 1 kcal/mol and all inside the
  first nine rows, each of the remaining ~90 steps under 1 kcal/mol;
- traded **ΔQM +1.49 for ΔMM -6.75 on the first step**.

So there was no basin hop and no roll-back toward the reactant. The run relaxed
strain that **the scan's handoff geometry was carrying**, at fixed chemistry.
The active region (476) and QM region (129) were identical in both runs, so this
is not the region mismatch of ¶316/¶317.

**¶333 measured something ADJACENT, not the same thing - the scope was wrong
until ¶344 and `retractions.md` records it.** On the NEB at 9 Å: of 47,979
atoms only 39-90 move past 0.1 Å and the largest displacement anywhere is 1.38
Å, yet the MM energies span ~20 kcal/mol, because *"'correctly minimised' means
the GRADIENT is small - it does not mean the water network found its global
minimum."* **But ¶333's movers were all first or second shell and `saddle1`'s are
8-12 Å out**, which is the opposite radial shape; the NEB's images have
genuinely different chemistry so their first shell must move, while `saddle1`
held O01-C4x at 2.02 Å and could only move environment. The two must not be
pooled as "the environment relaxes": a first-shell term is chemistry and cannot
be frozen away, a far-field term is Ryde's and cancels.

**This is a named problem in the literature, not a defect peculiar to us.** It
is Ryde's **local-minima problem** (J Chem Theory Comput 2017, 13, 5745): *"it is
very hard to ensure that the surrounding enzyme remain in the same local minima
throughout a reaction sequence"*, with ~20 kJ/mol per far-field hydrogen bond
and a compiled spread of **0.6-97 kJ/mol** across 24 published studies. Our
20.65 kcal/mol is 86 kJ/mol, inside that range - but the field reports it as
snapshot-to-snapshot scatter, and **no published instance of a single saddle
optimisation descending through it was found**. **¶344 then met Ryde's criterion
on its own terms**, which is about distance: *"the addition or removal of a single
hydrogen bond FAR FROM THE ACTIVE SITE affects the energy by ~20 kJ/mol in a way
that most likely is irrelevant for the reaction."* Our term is 9-12 Å out, on the
hydrophobic cluster Trp21/Phe18/Leu42/Asn41, with every reacting atom under
0.031 Å. It is Ryde's term, located.

**The field's four remedies, none of which this project does:** report free
energies; keep the surroundings fixed; run the coordinate forward *and*
backward; average over MD snapshots. Lonsdale, Harvey & Mulholland (Chem Soc Rev
2012, 41, 3025) are explicit that when relaxation to a conformationally distinct
path occurs *"the energies obtained should not be used"*. Freezing is cheap:
Calixto, Ramos & Fernandes (J Chem Theory Comput 2017, 13, 5486) measure the
cost of a frozen outer region at **0.8 kcal/mol on average, 1.6 maximum** above
a 6 Å free region, and **our 8 Å is above their overconstraint threshold**.

**What this does not say.** It does not say any individual run is wrong.
`aspxtb8`, `attack`, `attackmin` and `ringopen` all converged and their
component sums reproduce their printed totals to the last digit (¶317's
identity). It is the *differences between them* that are not meaningful, which
is why ¶332 - what environment the barrier is scored in - had to be settled
before any barrier could be quoted. **¶344 settled it: along the scan, with no
freely optimised point in the comparison.**

## Supports
- [[the-barrier-is-scored-along-the-scan]]
- [[no-reaction-coordinate-has-ever-been-computed]]
- [[convergence-gates-must-not-encode-science-criteria]]
