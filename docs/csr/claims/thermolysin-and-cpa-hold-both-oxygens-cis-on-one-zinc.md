---
type: observation
status: supported
evidence: experimental
entries: [¶124, ¶201]
checked: [¶124, ¶201]
---

# Thermolysin and carboxypeptidase A hold the nucleophile and the carbonyl oxygen cis on one zinc

Measured from the transition-state-analogue complexes themselves, not read from
reviews (¶124, re-measured with the angles in ¶201):

| | Zn-O(nuc) | Zn-O(carbonyl) | O...O | **O-Zn-O** | O-X-O at the centre |
| --- | --- | --- | --- | --- | --- |
| thermolysin 1TMN | 2.01 | 2.41 | 2.15 | **57.2** | 119.5 |
| carboxypeptidase A 7CPA | 2.19 | 2.98 | 2.68 | **60.1** | 117.2 |

**¶276: THESE ARE NOT A NUCLEOPHILE AND A CARBONYL. Both pairs are the two
oxygens of ONE functional group**, re-measured from the deposits:

    1TMN  OG1 2.009, OG2 2.413  both bonded to CB2 (1.249, 1.234)
          OG1-CB2-OG2 119.5 deg      O-Zn-O bite 57.2 deg   O...O 2.145
    7CPA  OPV 2.194, OPW 2.976  both bonded to P3  (1.603, 1.531)
          OPV-P3-OPW  117.2 deg      O-Zn-O bite 60.1 deg   O...O 2.676

The ligands are a **carboxylate** (an N-carboxymethyl transition-state analogue)
and a **phosphinate** - tetrahedral-intermediate mimics CHELATING the metal.
So "2.01/2.41" is a **chelate bite**, "57.2/60.1" a **bite angle**, and
"119.5/117.2" the **internal angle of that same group**. There is no scissile
carbonyl carbon in either picture and no separate nucleophile, so this node
cannot license a nucleophile-approach window.
| CSR-1 batch 1, last frame | 1.86 | 3.07 | 3.05 | **71.8** | 48.1 |

Both enzymes run the chemistry being proposed here - a zinc-bound water
attacking a zinc-coordinated carbonyl with a carboxylate as general base - on
**one** metal. Thermolysin's zinc is five-coordinate in this complex (His142
1.95, Glu166 1.95, His146 2.04, plus both inhibitor oxygens).

**Two things this settles.** A mono-zinc site *can* hold both the carbonyl and
a metal-bound nucleophile, so ¶123's geometric objection is false as stated
(¶124). And the two reference enzymes are **more cis than we are**, so the
reactant-stage Burgi-Dunitz gate was rejecting our structures on a criterion
the references also fail ([[burgi-dunitz-gate-at-reactant-was-wrong]]).

**The design criterion that survives:** one oxygen tight (2.01, 2.19), the
other **loose** (2.41, 2.98). A carbonyl at 2.4-3.0 A is where a carbonyl that
must still be attacked belongs.

**How far these complexes can be read** (¶124's own caution): they are
post-attack geometries, so 117-119 degrees at the central atom is a consequence
of bonds already formed, not a measured approach angle. Note too that CPA's
O...O of 2.67 A is identical to our failed geometries' 2.67-2.77, so oxygen
separation alone does not discriminate.

## Supports
- [[carbonyl-coordinates-zn1-in-the-michaelis-complex]]
- [[burgi-dunitz-gate-at-reactant-was-wrong]]
- [[question-can-burgi-dunitz-coexist-with-carbonyl-on-one-zinc]]
