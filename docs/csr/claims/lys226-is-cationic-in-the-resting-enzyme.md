---
type: claim
status: supported
entries: [¶191, ¶204, ¶209, ¶214, ¶216]
checked: [¶191, ¶204, ¶216]
---

# Lys226 is protonated, and therefore cationic, in the resting wild-type enzyme

The thesis needs this: an uncharged lysine cannot restrain a carboxylate by
electrostatics, and the whole substitution-for-Zn2 argument is electrostatic.

**What supports it.** A 2.8-3.0 A contact between NZ and a carboxylate oxygen in
three independent crystals ([[lys226-salt-bridges-asp92-in-wild-type-crystals]]) is
the geometry of a charged hydrogen bond; a neutral amine at that distance from an
anion would be unusual. The user's reading, 2026-11-05: salt bridges in the
crystal mean the lysine starts protonated. `eq6` was built on that and the
residue carries HZ1/HZ2/HZ3 (¶214).

**¶191 is cited above and is overturned**, so what it is cited for matters. ¶191
is in the retractions index for a different claim - Lys226 NZ-O01 at 2.72 A in
100% of 5000 frames, which ¶197 measured at 10.3% of frames and 4.93 A mean. That
was the lysine's contact with the *substrate*, not its protonation state, and
nothing in this node rests on it. The protonated topology ¶191's run was built
with is not what ¶197 overturned.

**What does not support it, and must not be cited as if it did.** The project has
previously modelled this residue **neutral** - `nucleophile-is-zn-bound-water-with-neutral-lys226`
is a standing node - and the QM/MM runs that transferred its proton away (¶204)
are the lineage where Asp92 then collapsed onto the zinc. So the project has
held both states, and the collapse artefact is *downstream of deprotonation*.
That makes this premise doubly load-bearing and it has **no pKa behind it**.

**Ruled `supported` by the user, 2026-11-06**, on the ground that a neutral
lysine could not form this contact in every crystal, stated as its own premise in
[[a-neutral-lysine-would-not-sit-at-28-angstroms-from-a-carboxylate]]. That is
the right shape of argument - an observed geometry plus general structural
chemistry - and it is now what this node rests on rather than the pKa.

**What the ruling does not cover, and the paper should not pretend it does.** The
resting state is settled; the state *during turnover* is not, and that is where
¶204's deprotonation sits. The pKa question stays open as a question rather than
as a gap in this premise ([[question-what-is-lys226-pka]]). PROPKA and pdb2pqr
are still not installed in the ASH environment (plan item 5, open since ¶7, two
pure-Python installs), so even the cheap estimate has not been made.

## Supports
- [[thesis-lys226-restrains-asp92-in-place-of-zn2]]

## Supported by
- [[lys226-salt-bridges-asp92-in-wild-type-crystals]]
- [[a-neutral-lysine-would-not-sit-at-28-angstroms-from-a-carboxylate]]
