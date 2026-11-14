---
type: observation
status: supported
evidence: ours
entries: [¶277]
checked: [¶277]
---

# Our own wild-type crystal is a better starting frame than 6DQH, and needs no reverts

`CSR_truncated_2metal.pdb` carries **1996 ANISOU records**, residues 46-303, and
reads wild type at every position the mechanism turns on: **HIS116, ARG118,
ASP120, GLN121, HIS192, SER214, LYS254, LYS263**. Only its metals are grafted
([[our-wt-crystal-zincs-are-grafted-from-the-double-mutant]]).

Measured with `refstruct.py`, metal at the two-histidine site:

    ZN - His116 NE2   2.25      ZN - His192 NE2   2.15
    ZN - Asp120 OD1   4.67      OD2  5.14      CG  5.36
    ZN - Arg118 CZ    4.36      ZN - Ser214 OG    4.99
    Lys254 NZ - Asp120 OD2      2.83

Against the three criteria ¶217 used to choose 6DQH:

| criterion | 6DQH | `CSR_truncated_2metal` |
| --- | --- | --- |
| metal at the modelled site | full occupancy | grafted, occ 0.53 |
| wild-type catalytic lysine | yes, 2.79 | yes, 2.83 |
| Asp92 outside the 3 Å radius | 4.28 / 5.14 / 5.15 | 4.67 / 5.14 / 5.36 |
| sequence reverts needed | **one (R118H)** | **zero** |
| crystallisation phosphate | **yes** (¶276) | none |

**The agreement on Asp92 is the load-bearing part.** ¶216 voided a transplanted
metal because the transplant *manufactured* a 2.80 Å Asp92 contact in a site
shaped without a metal. This graft manufactures nothing: it leaves the
carboxylate at 4.67-5.36, within 0.4 Å of what an independent full-occupancy
experimental structure measures, **OD2 identical at 5.14**, and **outside** the
3 Å restraint radius - so the rule has nothing to pin and the Asp92
controversy does not re-enter through the starting frame. The frame passes
¶216's own test rather than evading it.

**What the graft still cannot certify** is the histidine geometry:
[[the-wild-type-site-is-12-degrees-more-open-than-the-double-mutant]].

## Supports
- [[the-plain-double-mutant-is-the-mutant-arm]]

## Requires
- [[our-wt-crystal-zincs-are-grafted-from-the-double-mutant]]
- [[the-truncated-form-is-the-one-modelled]]
