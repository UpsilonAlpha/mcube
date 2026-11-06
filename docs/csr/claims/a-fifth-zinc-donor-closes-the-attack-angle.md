---
type: observation
status: supported
evidence: ours
entries: [D197, D198, D199]
checked: [D197, D198, D199]
---

# In our MM frames the nucleophilic attack angle tracks the fifth donor's distance: 39, 56 and 68 degrees

From `md_eq5`, three independent frames:

    Zn-O1x 1.92 A -> attack angle 39 deg
    Zn-O1x 2.06 A -> attack angle 56 deg
    Zn-O1x 3.74 A -> attack angle 68 deg

against the 95-120 deg window a nucleophilic attack needs. D197's crowding
mechanism is what this measures: a crowded zinc puts the nucleophile and the
carbonyl too close together. Even the best of the three frames is **27 deg short**
of the window.

**The heading is the measurement. The generalisation it was written with is
retracted.** D198 concluded that the site holds a fifth donor in MM whatever is
done, on the strength of frame 4374 - where **Asn225 OD1 took the site at 2.19 A**
the instant O1x left, the sphere stayed five-coordinate, and the angle opened by
only 12 deg. The retractions index overturns that: **D200 measured Asn225
leaving** over the QM/MM trajectory, to 4.18 A at closest and 5.86 A mean. So
"some ligand always fills the fifth site" is not established, and this node must
not be cited for it.

What survives is narrower and still useful: **in these MM frames the angle is
closed, and it opens as the fifth donor retreats.** Whether anything keeps a
fifth donor there once the model is QM/MM is open, and D211 gives a reason to
doubt it - O1x on the metal is an MM artefact that QM/MM does not reproduce
(mean Zn-O1x 4.24 A over 50 frames).

**Two further caveats from the index, read.** D197's own "PASS - coordination
number 4" is overturned by its own correction block: the four donors were
hydroxide, two histidines and O1x, with the carbonyl absent. D199's claim that
O1x being a 12-6 artefact makes the conclusions "stronger, not weaker" is listed
as an inference that immunises them, and is not relied on here. The index's
cross-cutting caveat for D170-D197 applies in full: all of this is our own MD.

## Supports
- [[asp92-on-zn1-interferes-with-hydrolysis]]
- [[no-usable-michaelis-complex-exists-yet]]
