---
type: question
status: open
entries: [¶124, ¶200, ¶201, ¶204]
---

# Can a nucleophilic attack trajectory coexist with the carbonyl on the same single zinc?

**REOPENED ¶276**, the reference that closed it being a misread
([[burgi-dunitz-gate-at-reactant-was-wrong]]). What survives independently of
any angle threshold, measured on the three tracked frames: the hydroxide sits
**31-35 deg out of the lactam carbonyl plane** (relax10 30.7, aspdft 32.2,
aspxtb8 35.1), near the nodal plane of C=O pi* and on the oxygen-s side rather
than over the carbon, where pi* addition wants ~75-85 deg. That is a statement
about position relative to the pi system, so it does not depend on a threshold.
It does not prove the pose wrong - the Zn-O(carbonyl) bond can lengthen along
the path - but a barrier from this reactant will be dominated by
reorganisation.

The central geometric question of the mechanism, and the one that decides
whether CSR-1 can run a Zn1-bound-nucleophile route at all.

**The obstruction, measured.** With the nucleophile and the carbonyl oxygen
both on the metal and cis, 82 degrees apart and 2.86 A from each other, the
attack angle is fixed near 29 degrees (¶200). Across every structure we have:
24-54 degrees over 791 frames (¶34), 28.9 mean over 50 ps (¶200), 16-33 in
`michrk` (¶204).

**The gate that made it look impossible was wrong** - see
[[burgi-dunitz-gate-at-reactant-was-wrong]]. The two mononuclear hydrolases
that run this chemistry are **more cis than we are**
([[thermolysin-and-cpa-hold-both-oxygens-cis-on-one-zinc]]), and 95 degrees
arrives at Nu-C = 1.64 A in thermolysin, 2.06 in CPA and 2.47-2.69 in ours, so
our geometry reaches the window *earlier* along the approach than either
reference.

**So the question is live rather than answered either way.** What has not been
done is to drive the approach and read the angle off the result instead of
gating a reactant on it - which is ¶207's exit: "a reactant whose proton
interconverts over a barrier of order kT does not need its protonation state
settled in advance, because the path settles it."

**Three ways out that no run has yet chosen between** (¶200): the pose is wrong
and a trans arrangement is accessible on a five-coordinate zinc (axial-axial in
a trigonal bipyramid, two His and Asp92 equatorial); the carbonyl is polarised
without being a ligand (the position ¶53 held and ¶121/¶123 withdrew); or the
mono-zinc assumption fails and the second site matters after all (plan item
10, a di-zinc comparison at identical level - **not started, a scope
decision**).

**That third branch bears directly on the thesis**, since it is the one outcome
that would return CSR-1 to the canonical mechanism.


**`aspxtb8`, the converged GFN2 QM/MM minimum (¶272):** O01-Zn-O2x **87.7 deg**,
O01...O2x 2.781, O01...C4x 3.626, **O01-C4x-O2x 38.5 deg**, the hydroxide 2.04 A
(34 deg) above the beta-lactam ring plane. Worse than batch 1-s 48.1, and
**still not a gate** - thermolysin and CPA reach 119.5 and 117.2 deg from
O-Zn-O of 57.2 and 60.1, i.e. more cis at the metal than we are. Logged as a
number. This argument has now been re-derived five times; the record caught it
each time.

## Supports
- [[no-usable-michaelis-complex-exists-yet]]
