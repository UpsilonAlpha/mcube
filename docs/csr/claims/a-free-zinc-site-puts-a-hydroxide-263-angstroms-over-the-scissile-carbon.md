---
type: observation
status: supported
evidence: ours
entries: [¶283]
checked: [¶283]
---

# With the carbonyl off the metal, a free zinc coordination position puts a hydroxide 2.63 angstroms over the scissile carbon

In the crystal-derived pose the whole beta-lactam ring sits ~3.2-3.3 Å from the
metal (Zn-C1 3.300, Zn-O1 3.153, Zn-N1 3.254), and the metal holds the substrate
**carboxylate** rather than the carbonyl
([[the-crystal-derived-pose-puts-the-carboxylate-on-the-metal]]).

**Sphere search.** 20,000 probe points at 1.95 Å from the zinc, keeping only
those >= 2.60 Å from His88 NE2, His164 NE2 and the bound carboxylate O8, making
>= 80° L-Zn-L angles to all three, and landing 2.50-3.40 Å from C1.
**722 survive - 3.6% of the sphere.** The best-oriented one:

    OH ... C1                 **2.63**
    out of the amide plane    **83.6°**      (90 = straight over the carbon)
    OH-C1-O1 / OH-C1-N1        87.3 / 87.6°
    OH - His88 NE2             3.02    (L-Zn-L  94°)
    OH - His164 NE2            3.47    (L-Zn-L 117°)
    OH - carboxylate O8        3.25    (L-Zn-L 118°)

Out-of-plane angles available across all 722: **60.2 to 83.6°**.

**Against the reactant we actually built:**

| | attack distance | out of plane |
| --- | --- | --- |
| `aspxtb8` | 3.626 Å | 35.7° |
| crystal-derived pose | **2.63 Å** | **83.6°** |

**The position is free precisely because the carbonyl is off the metal.** In
`aspxtb8` the carbonyl occupies that site at 2.211 Å and pins the hydroxide
into the carbonyl plane
([[aspxtb8-is-substrate-bound-but-not-a-near-attack-complex]]); here the
carboxylate holds the metal, Arg90 polarises the carbonyl at 2.98 Å, and the
vacant site points over C1.

**It also resolves the coordination-number question** without a gate: 2 His +
carboxylate + hydroxide = **four**
([[the-m1-site-is-four-coordinate-at-rest-and-not-when-substrate-bound]]).

**What it is not.** A feasibility test on empty space, not a structure. No
hydroxide has been placed, nothing solvated or relaxed, and the pose has never
been equilibrated ([[no-crystal-derived-pose-has-ever-been-equilibrated]]). The
probe idealises Zn-O at 1.95 Å with a 2.60 Å clash cutoff and an 80° angle
floor, and 3.6% is a real but modest window.

**One correction to my own method, recorded because it nearly gave the opposite
answer.** Placing the probe *along* the Zn->C1 vector puts it 1.350 Å from C1 -
a bond, not an approach - from which I would have concluded there is no room.
That is the degenerate collinear case; Zn-Nu 1.95 with Nu-C1 2.8 and Zn-C1 3.30
is a valid triangle at a Zn-Nu-C1 angle of ~86°.

## Supports
- [[the-metal-can-both-carry-the-nucleophile-and-stabilise-the-leaving-nitrogen]]
