---
type: observation
evidence: ours
status: supported
entries: [¶259, ¶260]
checked: [¶260]
---

# Our models hold the two Zn1 histidines 13–19° too closed, QM/MM opens the angle by only 4°, and Arg90 is displaced 2–2.5 Å from the metal as a consequence

Superposition-free measurement — all distances and angles are internal to each
structure, so none of this depends on a fit.

| structure | N–Zn–N | Zn–His88 NE2 | Zn–His164 NE2 |
| --- | --- | --- | --- |
| 6DQH, R118H crystal, 1.10 Å | **95.9°** | 2.080 | 2.054 |
| **dmCSR-1 2-metal, our own crystal, 1.47 Å** | **102.2°** | 2.230 | 2.046 |
| `aspxtb8`, converged GFN2 QM/MM | **86.6°** | 2.063 | 1.997 |
| `relax10`, 10 ns MM | **82.9°** | 2.178 | 2.126 |

**The bond lengths are right and the angle is not.** Every structure, ours
included, puts Zn–NE2 at 2.0–2.2 Å. The crystals open the N–Zn–N angle to
**96–102°**; our MM frame closes it to **82.9°** and QM/MM recovers only
**3.7°** of that, ending at 86.6° — still 9–16° short.

**Why QM/MM cannot fix it.** The angle is set by the protein backbone holding
the two histidines, which is MM in every QM/MM run here and frozen outright
where `actradius: 0` (¶235). So the metal's *bonding* relaxes and the
*scaffold* does not. That bounds what "QM/MM gives the correct coordination
behaviour" can deliver: it is true of bond lengths and of which donors are
present, and false of the geometry the protein imposes on them.

**Arg90's displacement follows from the same defect.** Arg90 stands at the
third Zn1 ligand position ([[arg90-stands-at-the-third-zn1-ligand-position]]),
and measured CZ-to-metal:

**CORRECTED ¶276. The 6DQH row measured the wrong atom: 6DQH IS the R118H
mutant, so it has no arginine at that position and no CZ. Its 4.14 / 5.01 /
7.58 reproduce exactly as **His118 NE2** to Zn / His116 NE2 / Asp120 CG.
**6NC5 has a real unmutated Arg89 with a metal in the modelled site** and was
omitted. Against the two real arginines WITH a metal present, ours is
**1.7-2.0 A out, not 2.5 A against "two crystals".**

| structure | Arg CZ–Zn | Arg CZ–His88 NE2 | Arg CZ–Asp92 CG |
| --- | --- | --- | --- |
| 6DN4, apo WT (no metal) | — | 5.98 | 7.23 |
| ~~6DQH, R118H~~ **His118 NE2, NOT an Arg** | (4.14) | (5.01) | (7.58) |
| **6NC5, real Arg89, metal at ZN302** | **4.95** | **5.72** | **8.11** |
| dmCSR-1 2-metal, real Arg | 4.60 | 5.32 | 8.05 |
| `aspxtb8` QM/MM | **6.64** | **6.82** | **8.95** |

Ours is the outlier on every coordinate, and worst against the metal: **6.64 Å
where two independent crystals give 4.14 and 4.60.** The apo wild type, which
has no metal at all, still holds the arginine closer to the histidines (5.98)
than our metal-bound model does (6.82), so this is not a consequence of the
metal being present or absent.

**Consequence.** A metal position taken from our own frames inherits the closed
angle, which is why `aspxtb8`'s zinc sits 1.06 Å from our own crystal's Zn1 and
1.31 Å from 5AXO's (¶257/¶259) — one defect, two apparent offsets. Rebuilding
should take the metal **and the histidine geometry** from a crystal rather than
the metal alone.

**What this does not say.** It does not say which crystal is the right source:
6DQH is the R118H mutant at the third ligand position, and our own 2-metal
structure has Zn1 at 0.53 occupancy with two of its three waters bridging to a
Zn2 the wild type lacks (¶259). It also does not say the closed angle changes
any *conclusion* reached so far — the Asp92 and Lys226 results of ¶244/¶245 are
second-shell distances measured on a metal whose own bonding is correct — only
that the site's geometry is not crystallographic.

## Requires
- [[relax10-holds-the-crystal-geometry-for-10-ns]]

## Supports
- [[no-usable-michaelis-complex-exists-yet]]
