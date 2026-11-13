---
type: claim
status: open
entries: [¶35, ¶41, ¶48, ¶201, ¶204, ¶206]
---

# The reactant is a Zn-bound water with a neutral Lys226 as the proton acceptor

The mechanism ¶201 left standing, and the one that makes the RQK lysine
catalytic rather than merely electrostatic: the nucleophile is Zn-OH2, and
neutral Lys226 - the state ¶35/¶48 establish at DFT - takes the proton. It is
an amine 2.82 A from the nucleophile and the only group besides a transient
water in hydrogen-bond reach of it (¶201). **It has never been computed as a
mechanism**, only observed as a resting-state proton transfer, because Lys226
sat outside the QM region from ¶190 until ¶202 put it back.

**The best energetic evidence is ¶48's -5.23 kcal/mol** for the Lys226 -> O01
transfer at wB97X-D3(BJ)/def2-SVP, in a constructed clean four-coordinate
sphere with identical constraints on both states. ¶206 explicitly keeps this
one: ¶48 was *handed* its geometry and asked only about protons, which is a
legitimate conditional, where ¶58's was loaded by a geometry no structure
supports. ¶35's -13.1/-15.5 kcal/mol is retracted - roughly two-thirds was an
artefact of the six-coordinate start (¶48).

**Status is `open` for three reasons.** The QM/MM observation of the transfer
comes from a trajectory whose Zn-Asp92 OD2 reaches 1.78 A
([[lys226-protonates-the-nucleophile-under-qmmm]]), a geometry the structural
criterion rejects. The pKa bookkeeping is **indeterminate**, not favourable:
¶107/¶112's "a 2+ cation 4-5 A away raises a lysine's pKa" is a sign error
corrected in ¶206, and the BcII aqua pKa of ~5.6 needs revising upward for a
mono-zinc site, which pushes back the other way. And the nucleophile is not
poised in any structure we have - O01-C4x 3.72-4.16 A (¶204).

## Supports
- [[thesis-mechanism-differs-from-canonical-b3]]

## Requires
- [[lys226-occupies-the-zn2-site]]
- [[question-what-protonates-the-leaving-nitrogen]]

## Supported by
- [[lys226-protonates-the-nucleophile-under-qmmm]]
