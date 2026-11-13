---
type: thesis
status: open
entries: [¶216, ¶22, ¶40, ¶107, ¶178, ¶209, ¶214]
checked: [¶40, ¶107, ¶178]
---

# Lys226 replaces Zn2 in the CSR-1 wild-type mechanism by restraining the bridging aspartate, allowing MBL activity

The root of this map, replacing the broader
[[thesis-mechanism-differs-from-canonical-b3]], which this entails and which is
kept as a claim.

**The argument in one line.** Canonical B3 enzymes hold Asp92's carboxylate on a
second zinc at ~2.1 A, which parks it ~4 A from Zn1 and leaves Zn1's sphere free
for the nucleophile and the substrate carbonyl
([[asp-zn1-distances-are-set-by-the-zn2-ion]]). CSR-1 has no Zn2
([[csr1-is-mono-zinc-with-zn2-vacant]]). In its place a lysine ammonium holds the
same carboxylate by a salt bridge
([[lys226-salt-bridges-asp92-in-wild-type-crystals]]), performing the one
structural job the missing metal did. Without that restraint the carboxylate
falls onto Zn1, where it is a fifth donor, and a five-coordinate Zn1 cannot
present the nucleophile at an attacking angle
([[asp92-on-zn1-interferes-with-hydrolysis]]). So the lysine is not a spare
cation near a vacancy - it is what makes the enzyme work.

**Status is `open`, and the specific thing that is open is the verb.** Every
structural half of this is measured in a deposited crystal. The claim that the
restraint *allows activity* is causal, and nothing computed yet separates it
from a correlation: no reaction coordinate has ever been run
([[no-reaction-coordinate-has-ever-been-computed]]), so there is no barrier with
and without the restraint to compare.
[[question-does-restraining-asp92-lower-the-barrier]] is that gap, named.

**Three ways this thesis could be wrong**, kept here so they are not lost:

1. The restraint could be real and incidental - Asp92 might sit off Zn1 for
   reasons that have nothing to do with the lysine, which is testable by the
   Lys226 charge-deletion already scoped in plan item 9.
2. Asp92-on-Zn1 might not actually cost anything. Our evidence that a fifth
   donor closes the attack angle is MM-only and was measured with *other* fifth
   donors ([[a-fifth-zinc-donor-closes-the-attack-angle]]).
3. The literature's reading of the aspartate is different from ours
   ([[asp120-mutants-lose-activity-in-b1-and-b3]]): in IMP-1, Asp120 mutants
   lose 10^2-10^4-fold in kcat *with both zincs still bound*, so the residue
   does something past positioning a metal. If that something is what CSR-1's
   lysine cannot supply, the thesis is incomplete rather than wrong.

**The blind-prediction protocol (¶178) applies.** The kinetics exist, are
withheld, and the computed mechanism goes in the log before they are released.

**Standing refutation of the causal clause, 2026-11-06 (¶216).** The structural
half - the lysine occupies the Zn2 site and holds the carboxylate out of the
metal's reach - is now measured at 1.10 A resolution in a metal-bound crystal
([[zinc-bound-crystal-holds-asp92-off-the-metal]]) and is in good shape. The
clause *"allowing for MBL activity"* is refuted as stated:
[[restoring-the-zn2-site-drastically-improves-catalysis]] shows that giving Asp92
a real Zn2 to bridge makes the enzyme **much better**, not worse. The thesis that
survives the evidence is weaker and should be written that way - **the lysine is
a structural stand-in that lets a mono-zinc B3 retain some activity at all,
holding the carboxylate in the second shell where it hydrogen-bonds the
nucleophilic water (2.64 A in 6DQH) instead of closing the metal's open face.**
It is not what makes the enzyme good at hydrolysis; on the published mutagenesis
it is part of what makes it poor.

## Requires
- [[canonical-b3-mechanism-is-dizinc-hydroxide]]
- [[csr1-is-mono-zinc-with-zn2-vacant]]
- [[csr1-zinc-sits-at-zn1]]
- [[csr1-hydrolyses-meropenem]]
- [[lys226-occupies-the-zn2-site]]
- [[lys226-is-cationic-in-the-resting-enzyme]]
- [[asp92-on-zn1-interferes-with-hydrolysis]]
- [[the-truncated-form-is-the-one-modelled]]
- [[question-does-restraining-asp92-lower-the-barrier]]

## Supports
- [[thesis-mechanism-differs-from-canonical-b3]]

## Supported by
- [[lys226-salt-bridges-asp92-in-wild-type-crystals]]
- [[engineered-zn2-holds-asp92-at-195-angstrom]]
- [[asp-zn1-distances-are-set-by-the-zn2-ion]]
- [[asp92-is-a-zn2-ligand-in-b3-structures]]
- [[asp120-mutants-lose-activity-in-b1-and-b3]]
- [[mutant-and-double-mutant-kinetics]]
