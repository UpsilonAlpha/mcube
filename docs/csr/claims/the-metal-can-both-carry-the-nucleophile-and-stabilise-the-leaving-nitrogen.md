---
type: claim
status: open
entries: [¶283]
---

# One zinc can carry the nucleophile and stabilise the leaving nitrogen at once, so neither cation is needed for the chemistry

In the crystal-derived pose the single metal is **3.254 Å from the leaving
nitrogen** - too long for a bond in the reactant, but positioned to close as the
ring opens. That is the role **Zn2** plays in canonical di-zinc B3, and the role
the single metal plays in CphA and GOB-18 in stabilising the anionic
intermediate.

And the same metal has a free coordination position putting a hydroxide **2.63 Å
over the scissile carbon at 83.6° out of the amide plane**
([[a-free-zinc-site-puts-a-hydroxide-263-angstroms-over-the-scissile-carbon]]).

**So one zinc does both jobs**, and the division of labour comes out as:

| role | canonical di-zinc B3 | CSR-1, this pose |
| --- | --- | --- |
| binds the substrate | Zn1 + Zn2 across the molecule | Zn - C3 carboxylate **1.847** |
| anchors that carboxylate | Ser221 / Thr223 | Ser186 OG **2.72 / 3.00** |
| polarises the carbonyl | Zn1 coordination | **Arg90 NH1 2.98** |
| carries the nucleophile | bridging hydroxide | the free Zn site, **2.63 from C1** |
| stabilises the leaving N | **Zn2** | the same zinc, **3.254** |

**Neither RQK cation is needed for the chemistry on this account.** Arg90 does
the polarisation and nothing else - it cannot reach N1
([[arg90-is-on-the-carbonyl-and-cannot-reach-the-leaving-nitrogen]]) - and
Lys226 has **no direct role in the bond-making and bond-breaking steps**,
which is consistent with it being busy holding
Asp92 ([[thesis-lys226-restrains-asp92-in-place-of-zn2]]).

That **sharpens** [[rqk-cations-replace-zn2-electrostatically]] rather than denting
it. That node already predicts the asymmetry - the cations anchor but cannot take
over the catalytic roles - and this measurement is where the line actually falls:
**Arg90 does the carbonyl polarisation, which is catalytic, and the metal does the
leaving-group stabilisation.** So the cations take over more than anchoring and
still not everything.

**"No direct role" is not "no impact on the chemistry", and the earlier wording
of this node overstated it.** Two indirect roles are real and neither is small:

1. **It keeps the enzyme mononuclear.** `refine_12` shows what a free Asp92
   does - it recruits a metal of its own at the M2 site and coordinates it at
   **2.13 Å**. If the salt bridge is what prevents that, it is the reason
   there is one metal rather than two, which is the entire premise of this
   mechanism ([[thesis-lys226-restrains-asp92-in-place-of-zn2]]).
2. **It is a second-shell electrostatic contribution.** A formal +1 at
   **5.04-5.28 Å** from the metal sits on a site that has to stabilise an
   anionic nucleophile and, on Lisa 2017-s account, an anionic intermediate.
   That shifts the hydroxide-s pKa and the intermediate-s energy without the
   lysine ever touching either.

So the honest statement is **no catalytic-step role, large structural and
electrostatic consequence** - and the second part is untested here, because
nothing in this work computes a pKa or an intermediate energy.

**Open, on three things.** The nucleophile position is empty space, not a
structure. Nothing is yet shown to protonate N1 - in this frame a water sits
2.565 Å away, but it belongs to the docked pose's solvation and its position
here is coincidental. And no crystal-derived pose has ever been solvated or
equilibrated ([[no-crystal-derived-pose-has-ever-been-equilibrated]]), which is
the one experiment that would test all of this at once.

**It is B1-like in the nucleophile, CphA/GOB-18-like in the metal's role and in
binding through the carboxylate, and novel in the arginine.** So it is a
composite rather than any of the three routes offered, and it does not require
the mechanism to be unprecedented in order to be new
([[the-novelty-is-a-monometallic-m1-site-with-two-histidines-and-an-arginine]]).

## Supports
- [[thesis-mechanism-differs-from-canonical-b3]]

## Requires
- [[the-crystal-derived-pose-puts-the-carboxylate-on-the-metal]]
- [[no-crystal-derived-pose-has-ever-been-equilibrated]]
