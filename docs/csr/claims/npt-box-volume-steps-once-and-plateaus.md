---
type: observation
status: supported
evidence: ours
entries: [¶225, ¶231]
checked: [¶225, ¶231]
---

# NPT box volume steps once and plateaus, so a first-to-last difference measures travel rather than convergence

**¶225 rejected a candidate frame on a box-volume drift of +5.40 nm3. ¶231
retracted that**, having measured the volume trace: under NPT the box steps once
early and then plateaus, so the difference between the first and the last value
reports how far it travelled on the way to equilibrium, not whether it has
equilibrated. `asppull3`'s box was converged by the criterion that matters while
failing the one that was applied.

The frame was rejected for a second reason that **stands** - CA RMSD still
climbing, three highest deciles last - so the decision survived its bad
argument. That is what makes this worth keeping as an observation: a
mis-specified convergence criterion rejected a frame, and only an unrelated
correct criterion stopped the error from mattering.

¶198 is the same failure in the other direction: judging His88 on a single last
frame (2.337 Å, outside the gate band) would have stopped a chain whose
trajectory read 100% occupancy, min 2.09 / mean 2.23 / max 2.57.

**¶225's box-volume criterion is listed in `retractions.md` and is cited here
only as the retracted thing this node is about.**

## Supports
- [[convergence-gates-must-not-encode-science-criteria]]
