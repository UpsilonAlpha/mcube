---
type: claim
status: supported
entries: [D33, D34, D200, D201]
checked: [D200, D201]
---

# The 95-120 degree Burgi-Dunitz gate, applied at the reactant, rejected the zinc-bound-nucleophile mechanism by construction

From D33 onwards, candidate structures were required to show an O01-C4x-O2x
angle of 95-120 degrees **at the reactant state**. Nothing ever met it: 24-54
degrees over 791 frames (D34), 28.9 degrees mean with 0% in window over 50 ps
of QM/MM (D200), 16-33 degrees in `michrk` (D204). D200 then inverted the
geometry and concluded that a 105 degree angle requires the nucleophile and the
carbonyl oxygen to be **trans** across the zinc.

**D201 measured the two reference enzymes and the gate failed.** 1TMN and 7CPA
hold the nucleophile and the carbonyl-derived oxygen **more cis than we do** -
O-Zn-O 57.2 and 60.1 degrees against our 71.8
([[thermolysin-and-cpa-hold-both-oxygens-cis-on-one-zinc]]). Their 117-119
degrees at the central atom is the angle of a tetrahedral centre whose bonds
have already formed, not an approach angle.

**The error was one of stage, not of value.** The attack angle is a function of
the approach distance; D200 held the reactant approach distance fixed while
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
