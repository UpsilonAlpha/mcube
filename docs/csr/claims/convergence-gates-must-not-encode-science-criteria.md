---
type: claim
status: unsupported
entries: [¶225, ¶231, ¶278]
checked: [¶225, ¶231, ¶278]
---

# A convergence gate may ask whether a run is trustworthy, not whether its answer is right

**Marked `unsupported` deliberately: this is a methodological norm, not an
empirical claim, so no observation can hold it up.** The two episodes below are
illustrations of the cost of ignoring it, not evidence for it. It is on the map
because it governs how every other node gets its status.

**The ¶231 gate mixes two kinds of test and only one of them belongs there.**

**Legitimate run-quality**, and `convergence.py` does exactly this and nothing
else: CA RMSD flat across the last deciles, temperature drift small, no NaN,
RMSF profiles correlating between windows and between independent halves, box
volume plateaued. These ask whether the trajectory can be read at all. They are
cheap, they are independent of the mechanism, and failing them invalidates a run
whatever it shows.

**Not legitimate:** *"unrestrained second-shell distances still matching the
crystal"*. That is a scientific finding wearing a convergence badge. Stated as a
convergence criterion it means a run that **disagrees** with the crystal is
marked unconverged rather than informative - so the gate can only ever confirm
the starting hypothesis, and a real disagreement is discarded as a technical
failure.

The distinction has already cost a decision in each direction. ¶225 rejected a
frame on box-volume drift, which ¶231 retracted: NPT box volume steps once and
plateaus, so a first-to-last difference measures travel, not convergence.
¶198's rule is the other half - judging His88 on a single last frame would have
stopped the chain wrongly, where the trajectory read 100% occupancy with mean
2.23 Å.

**So: run-quality tests gate a run; agreement with a crystal is a result, and
goes on this map with a status.** Where a project needs both - and it does - they
are reported separately and the second never fails a job.

## Supports
- [[every-acceptance-gate-encodes-a-mechanism]]
