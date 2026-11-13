---
type: question
status: void
entries: [¶9, ¶58, ¶59, ¶133, ¶195, ¶197, ¶206, ¶208, ¶231, ¶232]
---

# Does GFN2-xTB over-coordinate the zinc?

**Open since ¶9 (2026-10-15), and every closure of it has been withdrawn.**
The three closures - standing rule 3, ¶133, and project memory - all rest on
¶58, which is retracted (¶59 said so at the time, then re-asserted it in its
own closing line; ¶206 retracts it in place). The fourth, ¶195's 1.94 -> 2.56 A
relaxation, is reclassed by ¶197 as a bond-length relaxation *within* the
coordinated basin rather than a rejection of a wrongly-placed donor.

**No admissible test exists.** ¶206 rule 3 is why: a single point cannot
promote the geometry it was handed, and neither can agreement between two
methods on the same handed geometry.

**It no longer gates anything** (¶208 section 3). The discrepancy against the
crystals is too large to be a method artefact alone -
[[our-models-put-asp92-on-the-zinc]] against 3.94-4.20 A in five deposited
structures - so the geometry is wrong whatever produced it, and attributing the
error is a secondary question.

**What partly answered it, in the opposite direction to the worry:** on the
meropenem 1-hydroxyethyl contact, dd = GFN2 - DFT = **-5.92 kcal/mol**, i.e.
GFN2 *under*-binds (¶198). Same magnitude as Imhof's 7.55 for AM1/d, so a real
error, and one contact on one geometry pair rather than a statement about
coordination number. The cluster-model benchmark (plan item 7, open since ¶9) is
what would settle it.

## ANSWERED 2026-11-07 (¶232), in the affirmative, as pre-registered

**The test below ran and reported, and the prediction it would satisfy was
written first.** `aspqm2` chunk 1 brought Asp92 onto the zinc as a **bidentate
chelate within 8 ps** - OD2 2.10 A at 78% occupancy, OD1 2.29 A at 53%, Zn-CG
2.52 A - and held it for the remaining 24 ps with **Lys226 still protonated and
the salt bridge still formed at 2.87 A**. That is the second of the three
outcomes listed below, word for word, and it closes this question in the
**affirmative**.

Nothing here is revised to fit. The answer now lives in
[[gfn2-over-coordinates-zn1-from-a-crystal-consistent-start]] as a claim that
can be argued with; this node is kept for the record, with the predictions it
made still readable above the result, which is the whole point of having made
them.

Second reading, also as pre-registered: **HOH3095 does leave** - 0% occupancy
inside 2.5 A over 32 frames, and HOH3171 only 3% - so
[[csr1-zn1-site-is-h-r-h]]'s expectation that the extra water goes is upheld.
The sphere is five-coordinate anyway, because the carboxylate's two oxygens took
the vacated places.

## An admissible test now exists, and is running (¶231)

**Pre-registered before the result, per ¶178.** `aspqm2` chunk 1, job
32416288: 50 ps of GFN2-xTB/MM dynamics, **nothing restrained**, starting from
`md_relax10`'s last frame, whose Asp92 sits where four crystals put it
([[relax10-holds-the-crystal-geometry-for-10-ns]]).

This escapes ¶206 rule 3, which closes the door on single points and on two
methods agreeing about a handed geometry. It is neither: it is *dynamics*, free
to move, begun from a geometry the crystals support rather than one already
collapsed. ¶211's collapse to a 1.93 A mean is not a counter-example because
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

**A literature-absence claim to distrust:** ¶197's "no GFN2 benchmark for zinc
coordination number in enzymes appears to exist" has no recorded search
(`retractions.md` section 4).


## Superseded by
- [[gfn2-over-coordinates-zn1-from-a-crystal-consistent-start]]
