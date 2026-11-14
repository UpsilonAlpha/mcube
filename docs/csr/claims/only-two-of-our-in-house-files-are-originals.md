---
type: observation
status: supported
evidence: ours
entries: [¶277, ¶279]
checked: [¶279]
---

# Only two files in the honours directory are originals, and both are dmCSR-1 S214D

**Two files carry PHENIX refinement records**, both dated 2025-12-16, both from a
folder named `dmCSR S185D (S221D) Mutant`:

    dmCSR_mutant_1metal_refine_12.pdb   r_work 0.1609  r_free 0.1793  292 waters
    dmCSR_mutant_2metal_refine_16.pdb   r_work 0.1633  r_free 0.1864  312 waters
                                        REMARK: Zn occupancy sum 1.53

**Every other file is a superposed, metal-pasted derivative, and `CRYST1` proves
it.** Four distinct crystal forms share one set of zinc coordinates:

    CSR_nontruncated_*     I 4          84.580  84.580  74.982
    CSR_truncated_*        P 32 2 1     74.406  74.406 107.471
    dmCSR_*                P 21 21 21   53.468  74.946  76.095
    dmCSR_mutant_1metal*   P 21 21 21   69.392  77.543  43.698
    dmCSR_mutant_2metal*   P 21 21 21   43.760  69.361  77.388

Only two zinc models exist in the directory, pasted into every protein file: the
one-metal record (occ 0.70, B 9.03, at -13.098) and the two-metal pair (occ 1.00
/ 0.53, B 14.90 / 11.89). `AIM1_*.pdb` - a different protein in P 61 2 2 - sits
in the same coordinate frame with its zincs 0.3-0.5 Å away.

So the directory is **a comparative analysis set**: superposed into one frame,
waters stripped, metals harmonised, `CRYST1` left at each original value.
**`CRYST1` is reliable as provenance and the coordinates are not.**

**The one-metal record is `refine_12`'s own**, which means the wild-type files
carry the **S214D mutant's zinc**
([[our-wt-crystal-zincs-are-grafted-from-the-double-mutant]]).

**What this costs:** there is no original for wild-type CSR-1 and none for the
plain double mutant here, so for those the deposits are the only originals -
6DN4 and 6DQ2 (wild type, metal-free), 6DQH (R118H, one Zn), 6NC5 (plain double
mutant, di-zinc). What it does **not** cost is the sequence work: a superposition
cannot change a residue identity, so the alignments in
[[our-crystals-number-the-catalytic-lysine-254-not-263]] stand.

## Supports
- [[wild-type-csr1-binds-no-zinc-in-any-crystal]]
