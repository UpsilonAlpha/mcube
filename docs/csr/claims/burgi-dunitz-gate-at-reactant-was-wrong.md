---
type: claim
status: open
entries: [¶33, ¶34, ¶200, ¶201, ¶276]
checked: [¶200, ¶201]
---

# The 95-120 degree Burgi-Dunitz gate, applied at the reactant, rejected the zinc-bound-nucleophile mechanism by construction

**DEMOTED to `open`, ¶276.** The narrow point stands: a reactant should not
be gated on a transition-state angle. But the REPLACEMENT windows this node
licensed were taken from 1TMN/7CPA, and those two "Zn-O(nucleophile) /
Zn-O(carbonyl)" pairs are **both oxygens of one chelating functional group**
([[thermolysin-and-cpa-hold-both-oxygens-cis-on-one-zinc]]). So the windows
describe an already-formed gem-diolate, any cis pair of oxygens on the zinc
passes them, and the inversion argument extrapolates along a coordinate in
which the two oxygens are bonded to the same atom. ¶273 used this node to
dismiss a 38.5 deg attack angle; that dismissal is withdrawn.

From ¶33 onwards, candidate structures were required to show an O01-C4x-O2x
angle of 95-120 degrees **at the reactant state**. Nothing ever met it: 24-54
degrees over 791 frames (¶34), 28.9 degrees mean with 0% in window over 50 ps
of QM/MM (¶200), 16-33 degrees in `michrk` (¶204). ¶200 then inverted the
geometry and concluded that a 105 degree angle requires the nucleophile and the
carbonyl oxygen to be **trans** across the zinc.

**¶201 measured the two reference enzymes and the gate failed.** 1TMN and 7CPA
hold the nucleophile and the carbonyl-derived oxygen **more cis than we do** -
O-Zn-O 57.2 and 60.1 degrees against our 71.8
([[thermolysin-and-cpa-hold-both-oxygens-cis-on-one-zinc]]). Their 117-119
degrees at the central atom is the angle of a tetrahedral centre whose bonds
have already formed, not an approach angle.

**The error was one of stage, not of value.** The attack angle is a function of
the approach distance; ¶200 held the reactant approach distance fixed while
demanding a transition-state angle. The criterion is now demoted to a stage
`site.py check` never runs, and replaced at `reactant` by O01-Zn-O2x 45-95
degrees and O01...O2x 2.0-3.3 A, both measured from those two structures. Batch
1 passes both.

**One caveat to carry** (`retractions.md` section 4): "the project has been
rejecting its own structures on a criterion both reference enzymes also fail"
is partly our own extrapolation - the angle-versus-distance curve runs out to
distances those crystals do not occupy.

## Supports
- [[carbonyl-coordinates-zn1-in-the-michaelis-complex]]
- [[question-can-burgi-dunitz-coexist-with-carbonyl-on-one-zinc]]

## Supported by
- [[thermolysin-and-cpa-hold-both-oxygens-cis-on-one-zinc]]
