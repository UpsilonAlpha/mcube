---
type: claim
status: open
entries: [D191, D204, D209, D214]
checked: [D204]
---

# Lys226 is protonated, and therefore cationic, in the resting wild-type enzyme

The thesis needs this: an uncharged lysine cannot restrain a carboxylate by
electrostatics, and the whole substitution-for-Zn2 argument is electrostatic.

**What supports it.** A 2.8-3.0 A contact between NZ and a carboxylate oxygen in
two independent crystals ([[lys226-salt-bridges-asp92-in-wild-type-crystals]]) is
the geometry of a charged hydrogen bond; a neutral amine at that distance from an
anion would be unusual. The user's reading, 2026-11-05: salt bridges in the
crystal mean the lysine starts protonated. `eq6` was built on that and the
residue carries HZ1/HZ2/HZ3 (D214).

**What does not support it, and must not be cited as if it did.** The project has
previously modelled this residue **neutral** - `nucleophile-is-zn-bound-water-with-neutral-lys226`
is a standing node - and the QM/MM runs that transferred its proton away (D204)
are the lineage where Asp92 then collapsed onto the zinc. So the project has
held both states, and the collapse artefact is *downstream of deprotonation*.
That makes this premise doubly load-bearing and it has **no pKa behind it**.

**Status is `open` and should stay open until a pKa exists.**
[[question-what-is-lys226-pka]] is the gap. PROPKA and pdb2pqr are still not
installed in the ASH environment (plan item 5, open since D7, two pure-Python
installs), so even the cheap estimate has not been made.

## Requires
- [[question-what-is-lys226-pka]]

## Supports
- [[thesis-lys226-restrains-asp92-in-place-of-zn2]]

## Supported by
- [[lys226-salt-bridges-asp92-in-wild-type-crystals]]
