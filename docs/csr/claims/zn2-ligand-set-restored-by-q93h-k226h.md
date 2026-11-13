---
type: observation
status: supported
evidence: experimental
entries: [¶39, ¶40, ¶131, ¶208]
checked: [¶40, ¶131]
---

# 6NC5's Zn2 site exists only because two of its ligands were engineered back

`seqcmp.py` against 6NC5: our model is crystal numbering +1, 257/259 residues
identical, with exactly **two** mismatches - crystal His92 -> model Gln93 and
crystal His225 -> model Lys226. 6NC5 is the **Q93H/K226H double mutant**, and
its `SEQADV` records say so (¶39).

A metal appears in the Zn2 pocket only in that construct: ZN301, with Asp91 OD2
1.95, **His92** NE2 2.01, **His225** NE2 2.04, HOH493 2.08 - two of its three
protein ligands being the mutated residues. In wild-type CSR-1 they are a
glutamine and a lysine and cannot ligate a metal.

**The inference, and its limit.** ¶131's wording - that these records *prove*
Zn2 is vacant in the wild type - is overstated and listed as such; what they
establish is that the deposited construct is a double mutant, and wild-type
vacancy is an inference from that (`retractions.md` section 3). The inference is
a strong one: the only way anyone put a metal in Zn2 was to engineer two of its
ligands back.

**Also measured here** (¶39, confirmed ¶208): HOH493 **bridges** ZN301 and
ZN302 at 2.08/2.14 A, the classic di-zinc bridging-hydroxide position, and that
water is the one our model carries as the nucleophile O01.

## Supports
- [[csr1-is-mono-zinc-with-zn2-vacant]]
