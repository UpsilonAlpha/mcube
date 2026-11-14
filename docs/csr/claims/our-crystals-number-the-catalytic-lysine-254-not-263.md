---
type: observation
status: supported
evidence: ours
entries: [¶261, ¶277]
checked: [¶277]
---

# Our wild-type crystal has lysines at BOTH 254 and 263, and the one the double mutant replaces is 254

By pairwise alignment of chain A CA residues, with no numbering assumed:

    CSR_truncated_2metal.pdb  vs  dmCSR_2metal.pdb
        ratio 0.9865, TWO substitutions:   Q121->H    K254->H
    dmCSR_2metal.pdb          vs  dmCSR_mutant_2metal_refine_16.pdb
        ratio 0.9784, ONE substitution:    S214->D

So the residue the double mutant converts to histidine is **254**, and the serine
the aspartate mutant replaces is **214**, both established by sequence rather
than by offset arithmetic.

**This falsifies the hygiene note the mapping rests on.**
[[bbl-mapping-places-lys226-at-k263]] states that "the Lys254/Lys263 discrepancy
is one residue under two schemes - depositor sequential vs BBL". It is not: our
wild-type file contains **LYS at 254 and LYS at 263 simultaneously**, so they are
two different residues in one numbering, not one residue in two. Only the
lysine identification is affected; that node's H-R-H / D-Q-K site assignment is
untouched.

¶261 asked which of the manuscript's `S214D` and `K263H` is not BBL and left it
open. **The answer is `K263H`:** `S214D` names our crystals' serine exactly,
and `K263H` names a lysine that is present but is not the one mutated. To be put
to the lab with the alignment attached before submission.

**Still unsettled, and adjacent:** the second Zn1 histidine is **His192** in this
file and in 6DQ2, while the papers give **His196**, and residue 196 here is ALA.

## Refutes
- [[bbl-mapping-places-lys226-at-k263]]
