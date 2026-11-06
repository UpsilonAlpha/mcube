---
type: observation
status: supported
evidence: experimental
entries: [D108, D209]
---

# In the Q93H/K226H double mutant a zinc holds Asp92 at 1.95 A, and the two zincs are a bridged pair

Measured from 6NC5 directly (2026-11-06). The deposit holds **six** zincs; only
two are the site:

    ZN301  Asp91 OD2 1.95 | His92 NE2 2.01 | His225 NE2 2.04 | HOH493 2.08
    ZN302  His87 NE2 2.22 | His163 NE2 2.09 | HOH485 1.90 | HOH493 2.14
    ZN303  HOH485 1.88 | HOH504 2.13 | HOH516 1.96          <- no protein ligand
    ZN304  Glu115 OE2 1.87 | PO4 1.90                       <- surface
    ZN305  His111 ND1 2.36 | Glu141 OE2 1.91 | PO4 1.90     <- surface
    ZN306  Glu123 OE1/OE2 2.36/2.50 | Asp129 OD2 2.03       <- surface

(6NC5 numbering, which is ours minus one: its Asp91/His92/His225 are our
Asp92/Gln93->His/Lys226->His.)

**ZN301 is the engineered Zn2** and its ligand set is exactly the canonical one
the two mutations restore: Asp-His-His plus the bridging water. **ZN302 is Zn1**,
with two histidines rather than three because Arg90 stands at the third position
([[arg90-stands-at-the-third-zn1-ligand-position]]). **Zn1-Zn2 is 3.43 A,
bridged by HOH493** at 2.14/2.08 A - the canonical B3 mu-hydroxide geometry.

**This is the comparison the thesis rests on, in one deposit.** The same
carboxylate that a lysine ammonium holds at 2.8-3.0 A in the wild type
([[lys226-salt-bridges-asp92-in-wild-type-crystals]]) is held at **1.95 A by a
metal** once that metal is engineered back. Two cations, one job.

**Two cautions for anyone preprocessing this structure.** ZN303-ZN306 are
crystallisation artefacts and must be dropped explicitly or a metal detector
will restrain surface sites. And ZN303 shares Zn1's terminal water (HOH485 at
1.88/1.90 A from the two metals), so deleting ZN303 leaves that water's position
unsupported - and that water is where the nucleophile goes.

## Supports
- [[thesis-lys226-restrains-asp92-in-place-of-zn2]]
- [[asp-zn1-distances-are-set-by-the-zn2-ion]]
- [[zn2-ligand-set-restored-by-q93h-k226h]]
