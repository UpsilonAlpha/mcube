---
type: observation
status: supported
evidence: ours
entries: [¶277]
checked: [¶277]
---

# The zincs in our in-house wild-type CSR-1 file are the double mutant's, pasted in

`CSR_truncated_2metal.pdb` (WT CSR-1, P 3(2) 2 1) and `dmCSR_2metal.pdb` (the
plain double mutant, P 2(1)2(1)2(1)) carry **byte-identical zinc records**:

    HETATM 1955 ZN    ZN B   1     -12.832  18.175  -5.009  1.00 14.90
    HETATM 1956 ZN    ZN B   2     -14.099  16.870  -7.683  0.53 11.89

Same serials, same coordinates, same occupancies, same B-factors, in two
crystals of **different space groups**. Independent refinements cannot produce
that. The two proteins were superposed into one frame and one pair of zinc atoms
pasted into both.

**The pair belongs to the double mutant**, by chemistry rather than by
assumption: in the wild-type frame the full-occupancy zinc has no D/H/H triad to
sit on and lands **1.98 Å from Lys254 NZ**, 3.02 from Asp120 OD2 and 3.10 from
Gln121 OE1 - occupying the space where the *mutations* would have placed
ligands. In the double mutant the same atom is a proper site (Asp91 OD2 1.92,
His92 NE2 2.23, His225 NE2 2.34).

**So the in-house directory does not dissolve the wild-type build problem.** It
was briefly read as doing so; the B-factors are what caught it. Every statement
about a metal in a wild-type CSR-1 sequence remains a statement about a
placement.

## Supports
- [[wild-type-csr1-binds-no-zinc-in-any-crystal]]
- [[wild-type-csr1-crystals-contain-no-metal]]
