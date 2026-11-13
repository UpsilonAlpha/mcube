---
type: observation
status: supported
evidence: ours
entries: [¶51, ¶102, ¶201, ¶204, ¶208]
checked: [¶51, ¶102, ¶201, ¶204]
---

# Our models put Asp92 on the zinc at 1.77-2.11 A

Reproducible across methods and across five weeks of runs: MM at 100 ns
(¶102), unconstrained GFN2-xTB at 25 ps (¶51), QM/MM batch 1 with no MM
restraints - OD2 2.77 -> **2.00 A**, inside 2.5 A in 96% of 50 frames (¶201) -
and 1.78 A in the `michrk` trajectory (¶204). The range over all of it is
**1.77-2.11 A** (¶208 section 3).

**Evidence class `ours`, and that is the whole point of recording it this way.**
Against 3.94-4.20 A in five independent deposited structures, a 2+ A error is a
wrong geometry whatever produced it (¶208). Which component to fix - the MM
start, the restraint set, or xTB - is a real and secondary question
([[question-does-gfn2-over-coordinate-the-zinc]]).

Two necessary causes are now identified and both are fixed at the library
(¶208 section 4): the 3.0 A detection cutoff admitted a second-shell contact
into the restraint list at all, **and** the restraint builder's per-residue
default then held only the nearer carboxylate oxygen, leaving the other free to
pivot onto the metal. ¶120, ¶130 and ¶174 each blamed the cutoff alone; ¶175
named the default; nobody had both halves.

## Supports
- [[asp92-is-a-zn1-ligand-in-our-models]]
- [[gfn2-over-coordinates-zn1-from-a-crystal-consistent-start]]
