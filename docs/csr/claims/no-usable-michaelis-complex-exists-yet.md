---
type: claim
status: supported
entries: [D191, D195, D197, D198, D200, D201, D207]
checked: [D191, D195, D197, D198, D200, D201]
---

# No usable Michaelis complex has been built, and that is what blocks the mechanism

Eighty-four candidate last frames have been checked with `site.py check`; the
criteria have never all been met in one structure that also survives scrutiny
of its membership.

The record, compressed:

- **`wt100c`, 100 ns:** the substrate **leaves**. O01-C4x inside 3.5 A in 0.2%
  of 5000 frames, mean 19.95 A; confirmed real rather than a periodic-image
  artefact (D191/D192). Cause: one config line, an attack restraint replaced
  rather than added to (D194).
- **`mich_eq5`, 10 ns:** every restraint held to better than 0.11 A and the
  structure is still rejected - the pose's own 1-hydroxyethyl O1x took the
  fourth coordination site at 97.9% occupancy, unrestrained, while the carbonyl
  was inside 2.5 A in 0.3% of frames (D197). `site.py` reported "coordination
  number 4" correctly with the wrong membership.
- **Frame 4374 + QM/MM batch 1:** O1x is an MM artefact and leaves at the GFN2
  level (10% against 97.9%), Asn225 leaves, and the carbonyl coordinates
  unrestrained - but Asp92 bound the metal at 2.00 A in 96% of frames, so the
  sphere is not the designed one (D200 as corrected by D201).

**Two structural lessons worth keeping.** Testing properties one at a time can
pass a composition that would be rejected if asked directly (D197). And an
extremum of a distribution describes the reactant state we chose, not an
equilibrium average - frame 4374 was selected, and that must be said wherever
the run is quoted (D198).

**This is the thesis's real blocker, and it is prior to the choice of
instrument** (D207): a path run from a wrong reactant "would spend its barrier
moving meropenem across the site, and the number at the top would not be a
chemical barrier at all. That failure would look entirely normal in the
output."

## Supports
- [[no-reaction-coordinate-has-ever-been-computed]]
