---
type: question
status: open
entries: [D9, D58, D59, D133, D195, D197, D206, D208, D231]
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

## An admissible test now exists, and is running (D231)

**Pre-registered before the result, per D178.** `aspqm2` chunk 1, job
32416288: 50 ps of GFN2-xTB/MM dynamics, **nothing restrained**, starting from
`md_relax10`'s last frame, whose Asp92 sits where four crystals put it
([[relax10-holds-the-crystal-geometry-for-10-ns]]).

This escapes D206 rule 3, which closes the door on single points and on two
methods agreeing about a handed geometry. It is neither: it is *dynamics*, free
to move, begun from a geometry the crystals support rather than one already
collapsed. D211's collapse to a 1.93 A mean is not a counter-example because
that lineage **began** collapsed and with Lys226 deprotonated, so it tested
nothing. Starting distances: Zn-OD2 5.15, Zn-OD1 4.33, Zn-CG 5.19, Zn-NZ 5.28,
NZ-OD2 2.88, and a five-donor sphere at 2.50 A (HOH3095 2.02, hydroxide 2.04,
HOH3171 2.05, His164 2.13, His88 2.18).

**What each outcome would mean, written down first so it cannot be fitted
afterwards:**

- **Asp92 stays beyond ~4 A** -> GFN2 does not over-coordinate this zinc from a
  crystal-consistent start. The question closes in the negative, and
  [[our-models-put-asp92-on-the-zinc]] is explained by the *starting
  geometries* those models inherited rather than by the Hamiltonian.
- **Asp92 collapses inside ~2.5 A** -> GFN2 prefers a coordination five
  deposited structures exclude. The question closes in the affirmative, every
  GFN2 geometry in this project becomes suspect, and the def2-SVP leg stops
  being a refinement and becomes the only admissible instrument.
- **It sits between 2.5 and 4 A, or drifts without settling** -> 50 ps is too
  short and the chunked run continues; no closure either way.

**Second, independent reading from the same chunk:** whether **HOH3095
leaves**, which [[csr1-zn1-site-is-h-r-h]] predicts it should, since an H-R-H
site has room for two histidines, a hydroxide and one water. MM cannot decide
this and was not asked to.

**A literature-absence claim to distrust:** D197's "no GFN2 benchmark for zinc
coordination number in enzymes appears to exist" has no recorded search
(`retractions.md` section 4).

## Supports
- [[no-usable-michaelis-complex-exists-yet]]
