---
type: observation
status: supported
evidence: ours
entries: [¶286]
---

# A two-metal wild-type QM/MM did not keep its metal in the M1 site

Reported by the user (2026-11-15): a QM/MM run on wild-type CSR-1 built with
**two** zincs did not retain the metal in the **M1** site.

**Why it matters.** It is the first dynamical evidence for the premise the whole
mechanism rests on - that wild-type CSR-1 is mononuclear because the M1 site as
constituted (two histidines, with Arg90 in place of the third) cannot hold a
metal against a competent M2 site. Everything else supporting mononuclearity is
static: crystals with the metal at M2
([[the-zn2-aspartate-is-42-angstroms-from-the-m1-metal-with-or-without-the-lysine]]),
and the observation that a free Asp92 recruits its own metal at 2.13 Å.

It also speaks to the one thing ¶285 could not test. No crystal has a metal at
M1 and none at M2, so our single-metal arrangement had no structural check, and
our metal's position - between the two canonical sites - was the unverified
part. A di-zinc run shedding the M1 metal is consistent with that site being the
weak one.

**Not yet verified here, and recorded as the user's result rather than a
measurement of ours.** I have not seen the run directory, the starting
structure, the QM region, the charge, the length, or where the metal went - into
bulk, onto Asp92, or to M2. Which of those it was changes the reading
substantially:

- **to M2 or onto Asp92** - the aspartate out-competes the M1 site, which is
  the mononuclear argument directly.
- **into bulk** - weaker: it may say the M1 site is poor, or only that a
  two-metal model with one under-coordinated site is unstable, which is nearly
  tautological.

**Before this is cited in the manuscript it needs:** the run directory and
section name, the QM region and total charge, the simulation length, and the
metal's final coordination. Until then it is a pointer, not a result.

## Supports
- [[thesis-mechanism-differs-from-canonical-b3]]
