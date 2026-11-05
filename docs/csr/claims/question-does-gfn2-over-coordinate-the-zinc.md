---
type: question
status: open
entries: [D9, D58, D59, D133, D195, D197, D206, D208]
---

# Does GFN2-xTB over-coordinate the zinc?

**Open since D9 (2026-10-15), and every closure of it has been withdrawn.**
The three closures - standing rule 3, D133, and project memory - all rest on
D58, which is retracted (D59 said so at the time, then re-asserted it in its
own closing line; D206 retracts it in place). The fourth, D195's 1.94 -> 2.56 A
relaxation, is reclassed by D197 as a bond-length relaxation *within* the
coordinated basin rather than a rejection of a wrongly-placed donor.

**No admissible test exists.** D206 rule 3 is why: a single point cannot
promote the geometry it was handed, and neither can agreement between two
methods on the same handed geometry.

**It no longer gates anything** (D208 section 3). The discrepancy against the
crystals is too large to be a method artefact alone -
[[our-models-put-asp92-on-the-zinc]] against 3.94-4.20 A in five deposited
structures - so the geometry is wrong whatever produced it, and attributing the
error is a secondary question.

**What partly answered it, in the opposite direction to the worry:** on the
meropenem 1-hydroxyethyl contact, dd = GFN2 - DFT = **-5.92 kcal/mol**, i.e.
GFN2 *under*-binds (D198). Same magnitude as Imhof's 7.55 for AM1/d, so a real
error, and one contact on one geometry pair rather than a statement about
coordination number. The cluster-model benchmark (plan item 7, open since D9) is
what would settle it.

**A literature-absence claim to distrust:** D197's "no GFN2 benchmark for zinc
coordination number in enzymes appears to exist" has no recorded search
(`retractions.md` section 4).

## Supports
- [[no-usable-michaelis-complex-exists-yet]]
