---
type: observation
status: supported
evidence: ours
entries: [¶247, ¶281]
checked: [¶281]
---

# aspxtb8 is a well-anchored substrate-bound complex, and its attack geometry does not work

**The anchoring is real.** Two good hydrogen bonds hold the substrate:
Asn41 ND2 - O3x **2.859** on the C3 carboxylate and Arg229 NH1 - O5x **2.928**
on the side-chain amide, with Arg90 NH1 - O1x at 4.028 on the
6-alpha-hydroxyethyl oxygen. The metal sphere is clean and four-coordinate:
His88 2.063, His164 1.997, hydroxide 1.776, lactam carbonyl 2.211.

**The attack geometry does not.**

    O01 ... C4x      attack distance       3.626
    O01-C4x-O2x      approach angle         38.5°
    O01 out of the carbonyl plane           35.7°
    O01 - Zn - O2x   bite angle             87.7°
    O01 ... O2x      separation             2.781
    C4x=O2x / C4x-N1x                       1.212 / 1.371

**3.626 Å is not a near-attack conformation** - the convention is <=3.5 Å and
`site.yaml` asks 3.0-3.5. ¶126's substituent pose reached **3.17 Å
unrestrained**, so the preferred arm is further from the scissile carbon than the
arrangement it was meant to beat
([[question-does-the-substituent-on-metal-route-survive]]).

**The direction is the worse problem.** Addition to a carbonyl wants the
nucleophile over the carbon, near the pi* lobe. At 35.7° out of plane and 38.5°
to the C=O vector, with O01...O2x 2.781 and O01...C4x 3.626, the hydroxide sits
**beside the carbonyl oxygen pointing away from the carbon**, and must travel
~2 Å around it to attack. The cause is that both oxygens are ligands of the
same zinc.

One honest distinction from the ¶276 trap: in 1TMN and 7CPA the two oxygens were
**bonded to the same atom** (bites 57.2 / 60.1°), whereas these are a genuinely
separate nucleophile and carbonyl at 87.7°. So this is a real cis pair, not a
gem-diolate mimic - but a cis pair in the carbonyl's plane rather than above it.

**The substrate is unactivated**, which is correct for a Michaelis complex but
means the 2.211 Å Zn-O2x contact is doing little: C4x=O2x 1.212 Å is a free
carbonyl (7UHI 1.21) and C4x-N1x 1.371 Å an intact amide (7UHI 1.37).

**And the pose has not moved.** aspxtb8's meropenem centroid is **6.11 Å** from
7UHI's ligand transferred into this frame, against `relax10`'s 6.20 Å (¶247).
The 8 Å-active-region optimisation moved it **0.09 Å**, so ¶247's verdict -
right coordinate, wrong body - applies here unchanged.

**What follows.** A barrier from this reactant is dominated by reorganisation:
swinging the hydroxide around the carbonyl, or the carbonyl off the metal. That
is a genuine contribution but it is not the chemistry, and it will inflate any
number computed from here.

## Supports
- [[question-is-the-reactant-pose-experiment-anchored]]
- [[no-usable-michaelis-complex-exists-yet]]
