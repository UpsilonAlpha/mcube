---
type: claim
status: unsupported
entries: [¶281, ¶282]
checked: [¶282]
---

# The most likely route from our reactant is a Zn1-bound hydroxide with Ser186 as its proton partner - neither canonical B3 nor CphA

**DEMOTED to `unsupported` (¶282), on two counts.** Its Ser186 premise is void -
the residue is **not in the QM region**, so the 1.696 Å contact is an MM hydroxyl
across the QM/MM boundary at a distance shorter than a real hydrogen bond
([[ser186-hydrogen-bonds-the-nucleophile-in-our-reactant]]). And the ranking it
rests on was computed on `aspxtb8`, a pose **6.11 Å from the experimental mode**.

**On the pose with experimental backing the ranking inverts.** There the
metal-bound hydroxide is **4.702 Å** from the scissile carbon, so a Zn-OH
nucleophile is not available at all, and the metal holds the substrate
carboxylate at 1.847 Å instead
([[the-crystal-derived-pose-puts-the-carboxylate-on-the-metal]]).

**Kept because its rejection of the other two routes was right for the wrong
reason, and one half of it survives:** canonical B3 remains unavailable, and now
for a measured reason that holds in BOTH poses - Lys226 is 5.594 Å from the
leaving nitrogen in `aspxtb8` and **5.830 Å** in the crystal-derived pose. The
CphA dismissal was wrong and is replaced by
[[the-cpha-conditions-are-met-in-the-crystal-derived-pose]].


Ranked against the three candidates on the measured geometry of `aspxtb8`
([[aspxtb8-is-substrate-bound-but-not-a-near-attack-complex]]).

**Native B3-like: not available as written.** It needs a hydroxide bridging two
metals and **Zn2 stabilising the anionic nitrogen**
([[canonical-b3-mechanism-is-dizinc-hydroxide]]). We have one metal, a terminal
hydroxide, and a lysine **5.594 Å from N1x**
([[lys226-is-too-far-to-stabilise-the-oxyanion-or-the-leaving-nitrogen]]). The
only surviving version is "Lys226 substitutes for Zn2", which this geometry does
not support by contact.

**B2 / CphA-like: architecture absent.** CphA is a **Zn2-site** enzyme whose
nucleophile is a non-metal-bound water activated by Asp120, with the metal
binding the substrate carboxylate. Here the metal is at **Zn1**, the nucleophile
**is** metal-bound at 1.776 Å, and the metal is 5.42-6.70 Å from both
carboxylate oxygens. [[cpha-route-is-closed]] is `open` because its stated
closure rested on a refuted premise, but the structural objection inside it
stands.

**B1-monometallic-like: closest, and what this complex is a draft of.** A single
zinc at the histidine site carrying the hydroxide that attacks the lactam
carbonyl, as proposed for mono-Zn BcII
([[b1-monozinc-bcii-carries-its-metal-at-zn1]]). Two differences from that
precedent, the first of which is the novelty: our site is **two histidines plus
Arg90** ([[csr1-zn1-site-is-h-r-h]]), and in BcII the substrate carbonyl is not
also a metal ligand.

**The variant the structure points at adds one thing, and it is the novel
part:** Ser186's hydroxyl donates a **1.696 Å** hydrogen bond to the
nucleophile ([[ser186-hydrogen-bonds-the-nucleophile-in-our-reactant]]). That
gives the hydroxide a proton acceptor without needing Asp92 as a general base or
a second metal-bound water, and it sits **between the two literatures**: in kind
it is what the only published mono-zinc mechanism predicts - activation by a
hydrogen-bond network rather than by the metal
([[the-published-mono-zinc-mechanism-has-no-metal-activated-nucleophile]]) -
while disagreeing with it in detail, because ours is on the metal as well.

**Open, and three things would have to be true.** The pose is **6.11 Å** from
the experimental mode, so the attack geometry is a property of a docked pose.
Nothing is positioned to protonate N1x. And the hydrogen bond is a GFN2 number
needing re-measurement at DFT (¶245: GFN2 runs ~0.13 Å short on Zn-O).

**It also makes S214D the diagnostic mutant**, for a different reason than ¶277
allowed: if Ser186 is the proton partner, replacing it with an aspartate removes
that partner and installs a competing base - and the manuscript's own fork,
"the aspartate residue directly prevents binding or inhibits hydrolysis", is the
prediction
([[the-manuscript-reads-s214d-as-inhibition-by-the-bridging-aspartate]]).

## Superseded by
- [[the-cpha-conditions-are-met-in-the-crystal-derived-pose]]

## Requires
- [[question-is-the-reactant-pose-experiment-anchored]]
