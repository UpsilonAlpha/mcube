---
type: observation
evidence: experimental
status: supported
entries: [¶217]
checked: [¶217]
---

# Asp92's B-factors are below average with no alternate conformer, while Lys226's ammonium is the mobile end

The check demanded before any claim that the crystal shows no equilibrium, and
it splits the two residues.

    6DQH (1.10 A; protein mean B 13.2, median 11.3, side-chain mean 14.4)
      Asp120   CB 9.6  CG 9.9  OD1 10.0  OD2 11.9       <- below average, no altloc
      Lys260   CB 10.6 CG 12.5 CD 12.9 CE 13.6 NZ 14.6  <- rising along the chain
    6DN4 (1.99 A; mean 16.2, median 13.8, side-chain 17.4)
      Asp120   CB 14.2 CG 15.4 OD1 16.4 OD2 16.3        <- about average
      Lys254   CB 13.9 CG 17.0 CD 19.1 CE 22.1 NZ 24.8  <- 1.5x the mean at the tip
    ours1 (1.27 A; mean 13.2, median 11.4, side-chain 15.3)
      Asp78    CB 10.9 CG 10.7 OD1 9.9  OD2 11.5        <- below average, on Zn2

**For the aspartate the claim holds.** Every side-chain atom sits under the
protein median, occupancy is 1.00, and Asp120 is not among the 12 residues in
6DQH that carry an alternate conformer. A carboxylate exchanging between a metal
at ~2.1 A and an ammonium at ~2.8 A would show two conformers or elevated B; at
1.10 A it shows neither. So **4.28 A from the metal is a single real position**,
not a two-state average.

**For the pair it does not.** Lys226's B climbs monotonically out the side chain
and ends above the protein mean in 6DQH and at **1.5x** the mean in apo 6DN4.
The honest statement is **"the carboxylate is pinned and the ammonium reaches
it"**, not "the lysine clamps the carboxylate" - which makes the thesis's verb
weaker than written, and makes the apo structures much thinner evidence for the
salt bridge than their 2.83-2.99 A distances suggest. It is a second reason to
take numbers from 6DQH rather than from 6DN4.

**What B-factors cannot exclude** is a minor state below the density threshold -
fast exchange that is fully shifted. That is a free-energy question, not a
crystallographic one.

## Supports
- [[zinc-bound-crystal-holds-asp92-off-the-metal]]
- [[question-is-asp92-in-equilibrium-between-lys226-and-zn1]]
