---
type: observation
status: supported
evidence: ours
entries: [¶187, ¶188, ¶198]
checked: [¶187, ¶188, ¶198]
---

# GFN2-xTB overstates the metal-coordination preference by 5.9 kcal/mol against DFT on this site

The cluster benchmark (plan item 7) has already been run and its numbers were
never recorded. From `CSR_mich/config.yaml`, two 51-atom site clusters differing
only in whether the substrate's O1x is on the metal, both at charge 0:

    level                  O1x on        O1x off       off - on
    GFN2-xTB            -80.274541    -80.297916    -14.67 kcal/mol
    wb97x-v/def2-svp  -3014.167570  -3014.181518     -8.75 kcal/mol

**Same sign, factor 1.68.** Both levels agree that the alcohol oxygen off the
metal is favoured, and GFN2 exaggerates it by **5.9 kcal/mol** - on exactly the
kind of metal-coordination question the thesis turns on.

**What it licenses and what it does not.** It licenses using GFN2 to find and
order coordination states: the sign is right. It does not license quoting a GFN2
energy difference as a number, and it means any Asp92 approach profile computed
at GFN2 should be expected to be too steep by a similar factor - relevant to
`scan_znasp`'s +8.46 kcal/mol. It is one pair of structures on one coordinate, so
5.9 kcal/mol is a single calibration point, not an error bar.

**¶198's conclusion from these same clusters is listed as stronger than its
evidence, and this node deliberately draws the opposite one.** ¶198 read them as
"GFN2 does not over-coordinate the zinc"; the index notes its own caveat limits
that to one contact, one geometry pair and single points. The energies above are
the same numbers, and what they show is a 1.68-fold exaggeration of the same
sign - a calibration, not a clearance.

**¶187/¶188's own numbers were corrected and the corrected ones are used here.**
The retractions index lists both for describing the cluster as net -1 with a
~15-atom substrate fragment; ¶198 measured both clusters at **51 atoms,
C15H25N5O5Zn, charge 0**, which is what the config records and what is quoted
above.

**This is `ours` twice over** - our energies on our geometries - so it is evidence
about the energy functions, not about the site.

## Supports
- [[gfn2-over-coordinates-zn1-from-a-crystal-consistent-start]]
- [[question-does-restraining-asp92-lower-the-barrier]]
