---
type: claim
status: unsupported
entries: [¶200, ¶201, ¶276, ¶278]
---

# The reactant must hold the nucleophile and carbonyl cis on the metal, O-Zn-O 45-95 degrees and O...O 2.0-3.3 angstroms

**The gate `site.yaml` applies as `O01 - Zn - O2x, range [45, 95]` and
`O01 ... O2x separation, range [2.0, 3.3]`, both sourced "measured 1TMN 57.2,
7CPA 60.1" and "measured 1TMN 2.15, 7CPA 2.68" (¶201).** These were written as
the replacement for the Burgi-Dunitz gate ¶201 removed from the reactant - what
the angle "was meant to ask, asked at the reactant instead".

**Marked `unsupported`: the provenance is a misread** and ¶276 established it by
re-measuring the deposits. In both reference structures the two oxygens are
**the two oxygens of a single chelating functional group**, not a nucleophile
and a substrate carbonyl:

    1TMN  OG1 2.009, OG2 2.413   both bonded to CB2 (1.249, 1.234)
    7CPA  OPV 2.194, OPW 2.976   both bonded to P3  (1.603, 1.531)

A carboxylate and a phosphinate - transition-state analogues chelating the
metal. So "57.2 / 60.1 degrees" is a **chelate bite angle** and "2.15 / 2.68 Å"
an **intramolecular O...O distance**, and the windows describe an
**already-formed gem-diolate** rather than a reactant. Any cis pair of oxygens
on the zinc passes them.

**What this does not do.** It does not restore the Burgi-Dunitz gate, which
¶201 removed for an independent and still-valid reason
([[burgi-dunitz-gate-at-reactant-was-wrong]]). It leaves the reactant with **no
justified geometric gate on the nucleophile-carbonyl relationship**, which is an
honest state and not a licence to accept any pose.

**What survives as a measurement, threshold-free:** the hydroxide sits **31-35
degrees out of the lactam carbonyl plane** in all three tracked frames (relax10
30.7, aspdft 32.2, aspxtb8 35.1) - near the nodal plane of C=O pi* and on the
oxygen's side rather than over the carbon. That is a statement about position
relative to the pi system, so it needs no window
([[question-can-burgi-dunitz-coexist-with-carbonyl-on-one-zinc]]).

## Supports
- [[every-acceptance-gate-encodes-a-mechanism]]

## Refuted by
- [[thermolysin-and-cpa-hold-both-oxygens-cis-on-one-zinc]]
