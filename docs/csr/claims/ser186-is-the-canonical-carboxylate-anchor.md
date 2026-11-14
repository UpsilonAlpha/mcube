---
type: observation
status: supported
evidence: experimental
entries: [¶284]
checked: [¶284]
---

# Ser186 is the canonical Ser221 carboxylate anchor, by coordinates rather than by numbering

Chain-A CA superposition onto our frame (Needleman-Wunsch, Kabsch, outliers
dropped at max(3.0, 2*RMSD)):

    L1    (7UHI)  Ser206 OG   **1.44 Å** from our Ser186 OG   core 221/258, RMSD 1.22
    SMB-1 (5B15)  Ser175 OG   **1.64 Å** from our Ser186 OG   core 219/254, RMSD 1.25
    L1            Ser208 OG     4.12 Å
    SMB-1         Thr177 OG1    3.61 Å

**1.44 and 1.64 Å are inside the fits' own RMSD**, so the first anchor of the
canonical Ser221/Thr223 pair and our Ser186 are the same position. In both
deposits that serine hydrogen-bonds the metal-bound carboxylate oxygen (L1
Ser206 OG - O8 3.04; SMB-1 Ser175 OG - OAE 2.76), and in our crystal-derived
pose Ser186 OG does the same at 2.72 / 3.00 Å.

**This settles the S214 / S221 ambiguity structurally.** `reference-structures.md`
could not reconcile the manuscript's `S214D` with a BBL `S221` from sequence
offsets. It does not have to: the residue is the canonical anchor whatever the
convention calls it, and that is the statement the mechanism needs.

**The second anchor has no counterpart in CSR-1.** What occupies L1's Ser208 OG
position in our frame is backbone - **Ala188 N 2.44 Å, Ala188 O 2.78** - so this
enzyme has one of the two canonical carboxylate anchors, not two
([[the-substrate-carboxylate-end-is-unanchored-in-our-models]]).

## Supports
- [[the-metal-can-both-carry-the-nucleophile-and-stabilise-the-leaving-nitrogen]]
- [[s214d-occupies-the-substrate-carboxylate-site]]
