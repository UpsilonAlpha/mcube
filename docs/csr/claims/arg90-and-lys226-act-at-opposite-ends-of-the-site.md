---
type: observation
status: supported
evidence: ours
entries: [D107, D111, D153, D197]
checked: [D107, D111, D153, D197]
---

# Arg90 and Lys226 act at opposite ends of the site, so they do not jointly replace one charge

Measured in `QM_MM_opt_last.pdb` and carried in the log's central-hypothesis
block: **Arg90 NH1 is 4.41 A from the zinc** and 5.13 A from the nucleophile,
but **Lys226 NZ is 8.73 A from the zinc** and points at the substrate
carboxylate (4.80 A) instead. A Zn2 metal would sit ~3.5-4 A from Zn1, so the
two cations act at opposite ends of the site rather than jointly replacing one
charge.

And **Arg229 NH2 is 4.38 A from the substrate carboxylate and 6.61 A from the
leaving nitrogen** - closer than Lys226 on both counts - which means the
electrostatic role the hypothesis assigns to the RQK set may belong partly to a
residue outside it.

**This refutes the hypothesis as literally worded** ("together these two supply
the +2"), not the broader idea that cations substitute for Zn2
electrostatically. What it removes is the geometric picture of two cations
standing in for one ion at one position.

**Evidence class `ours`, and weakly held even within that.**
`retractions.md` section 4 lists exactly this: the two distances are called
"two measured caveats" in the document's most-read section, and they are **one
optimised frame of our own model, from a condemned lineage**. Lys226's position
is also pose-dependent - NZ to the nucleophile is inside 3.5 A in 10.3% of
frames in `eq5`, mean 4.93 A (D197), against D191's retracted "100% of 5000
frames".

## Refutes
- [[rqk-cations-replace-zn2-electrostatically]]
