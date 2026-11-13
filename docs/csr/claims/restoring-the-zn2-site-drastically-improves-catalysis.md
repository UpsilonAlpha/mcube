---
type: observation
evidence: literature
status: supported
entries: [¶216]
checked: [¶216]
---

# Giving Asp92 a Zn2 to bridge makes CSR-1 a much better lactamase, not a worse one

Krco et al. (Front Chem 2023;11:1196073), reviewing this enzyme and citing
Pedroso et al. 2020 Supplementary Table S2, verbatim:

> replacing both the glutamine and lysine residues in the CSR-1 Zn2 site with
> their canonical histidine counterparts enhances the affinities of both Zn(II)
> ions and **drastically improves the catalytic efficiency to a level comparable
> to that of other B3 MBLs**

and for the Zn1 substitution:

> replacement of its arginine in position 118 by the canonical histidine
> enhances metal ion binding affinity at the Zn1 site with a **modest
> improvement in catalytic efficiency**

Zinc affinities from Pedroso et al. (ITC): wild type alpha_Kd ~2 uM / beta_Kd
~150 uM; R118H 0.1 / 180; **Q121H/K263H 0.14 / 0.14**; triple 0.06 / 0.4.

**This is the hardest fact in the map for the thesis, and it points the opposite
way to the mechanism the thesis asserts.** The double mutant is precisely the
construct in which Asp92 becomes a true bridging ligand - 6NC5 measures it at
1.95 A from the engineered Zn2 ([[engineered-zn2-holds-asp92-at-195-angstrom]]).
If a carboxylate on a metal were what suppresses hydrolysis, that construct
should be worse. It is dramatically better.

**What survives.** The thesis's structural half is untouched: the lysine does sit
in the Zn2 site, it does hold the carboxylate, and 6DQH measures it
([[zinc-bound-crystal-holds-asp92-off-the-metal]]). What does not survive is the
causal clause *"allowing for MBL activity"*. On this evidence the lysine is what
**limits** CSR-1's activity rather than what permits it, and the defensible
reading is the weaker one: the lysine is a structural stand-in that lets a
mono-zinc B3 retain *some* activity, not the reason it has activity.

**Two reservations, both real.** The numbers are in a supplementary table neither
I nor the review quotes, so the fold changes are unread - "drastically" is the
review's word. And the double mutant changes two things at once: it adds a
second metal *and* removes the lysine. Nothing here isolates the carboxylate's
own contribution, and **no ¶120 mutant of this enzyme or of SPR-1 exists in
either paper** - the claim that removing the bridging aspartate raises activity
in this subclass is not in the deposited literature.

## Refutes
- [[thesis-lys226-restrains-asp92-in-place-of-zn2]]
- [[asp92-on-zn1-interferes-with-hydrolysis]]
