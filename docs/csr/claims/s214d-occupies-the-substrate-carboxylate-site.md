---
type: observation
status: supported
evidence: experimental
entries: [¶284]
checked: [¶284]
---

# S214D inhibits by occupying the substrate's carboxylate site, not by touching the metal

Our frame fitted onto `dmCSR_mutant_2metal_refine_16.pdb` - the same protein, so
the fit is near-exact: **core 250/250, RMSD 0.82 Å**. Carrying the
crystal-derived ligand pose and our Ser186 across:

    our Ser186 OG lands   **1.47 Å** from Asp172 OD2
    substrate O7 lands    **1.66 Å** from Asp172 OD2
    substrate O8 lands    **2.79 Å** from Asp172 OD2    (O8 is the metal-bound oxygen)
    Asp172 OD1 to O7        3.49        to O8  4.91

(`refine_16` numbers residues ours minus 14, so Asp172 is our 186 mutated to
aspartate, and Asp78 is our Asp92.)

**A carboxylate oxygen of the mutant sits 1.66 Å from where a carboxylate
oxygen of the substrate has to sit.** That is the same point in space, anion
against anion - not a perturbation of a hydrogen bond but the removal of the
site. Given that this serine is the canonical anchor
([[ser186-is-the-canonical-carboxylate-anchor]]), S214D is a clean test of the
anchor.

**And it is not a test of the metal.** Asp172 is **4.13 Å** from M1 and 4.06
from M2 in this crystal, and 4.23 Å from the single metal in `refine_12`
([[neither-aspartate-coordinates-the-metal-in-our-own-crystals]]), so the
"bridging aspartate poisons the zinc" reading has no support in either original.
This reverses
[[the-manuscript-reads-s214d-as-inhibition-by-the-bridging-aspartate]] and
answers the standing worry in [[s214d-cannot-isolate-the-serine-role]]: the
mutation isolates the serine's role better than expected, because its side chain
never reaches the metal.

**What it does not explain is the carbapenem selectivity.** Cephalosporins carry
a C-4 carboxylate too, so a block at this site should cost activity on both. The
only published mechanism that would make the loss selective is Lisa 2017's
anionic intermediate, which accumulates for carbapenems and therefore depends on
how well the leaving nitrogen is stabilised - and losing the anchor shifts the
ring, which costs the metal its 3.254 Å grip on N1. **That is a hypothesis with
no structural support here**, and settling it needs kinetics this work does not
use.

**It picks one of the manuscript's own two branches rather than contradicting it.**
[[the-manuscript-reads-s214d-as-inhibition-by-the-bridging-aspartate]] quotes the
paper offering "either that the serine is critical ... **or that the aspartate
residue directly prevents binding**". This is the second branch, measured: the
aspartate prevents binding by sitting in the carboxylate's place. What it does
refute is the *coordination-environment* mechanism the abstract leans on, since
Asp172 never comes within 4.06 Å of either metal.

**Steric or electrostatic? Both, and the steric part is the one that survives
the error bar.** At **1.66 Å** two carboxylate oxygens are not merely repelling
each other, they are inside one another - the sum of two oxygen van der Waals
radii is about 3.0 Å - so the ligand simply cannot occupy that point. Like
charge on top of that only makes it worse. **But 1.66 Å carries the pose-s own
uncertainty**, and a global-fit transfer moves the ligand by mean 1.03 Å (max
1.44), so the honest form is **overlapping within error** rather than 1.66 Å
exactly.

**The solid half of the clash does not involve the ligand.** Asp172 OD2 lands
**1.47 Å** from our **Ser186 OG**, and Ser186-s position is fixed by two
independent deposits on protein atoms alone (L1 Ser206 1.44 Å, SMB-1 Ser175
1.64 Å). So the mutant aspartate demonstrably occupies the **anchor residue-s**
own hydroxyl position, which needs no ligand pose at all, and the substrate
clash follows from the anchor being gone rather than being evidence in its own
right ([[ser186-takes-over-the-canonical-anchor-role]]).

## Supports
- [[ser186-is-the-canonical-carboxylate-anchor]]
- [[the-manuscript-reads-s214d-as-inhibition-by-the-bridging-aspartate]]
