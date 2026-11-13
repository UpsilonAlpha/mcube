---
type: question
status: open
entries: [¶124, ¶200, ¶201, ¶204]
---

# Can a nucleophilic attack trajectory coexist with the carbonyl on the same single zinc?

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

## Supports
- [[no-usable-michaelis-complex-exists-yet]]
