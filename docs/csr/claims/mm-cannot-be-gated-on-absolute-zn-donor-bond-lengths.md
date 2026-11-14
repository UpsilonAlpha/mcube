---
type: claim
status: supported
entries: [¶48, ¶118, ¶196, ¶278]
checked: [¶48, ¶196, ¶278]
---

# Gating an MM structure on a 1.9-2.2 angstrom Zn-donor bond tests the force field, not the structure

**The gate: `Zn - nucleophile O01, range [1.9, 2.2]`, sourced in `site.yaml` as
"tight; thermolysin 2.01, CPA 2.19 - crystal values".** Those are crystal
values, and the Li/Merz 12-6 zinc does not reproduce crystal bond lengths. In a
grafted frame the whole sphere is systematically off - hydroxide **1.86 Å**
(short), His88 **2.22** and His164 **2.32** (long) - and both ¶48 and ¶118
measured that bond relaxing from 1.83-1.86 Å to **1.99-2.05 Å** as soon as
QM/MM took over.

So the window is correct about crystals and uninformative about an MM frame,
which is why `required_at` was moved to `reactant` only. **The criteria that
remain at `start` are the ones that DISCRIMINATE between the hypotheses** -
which oxygen is on the metal, and the substrate geometry - where the difference
is ~1.8 Å and no force-field artefact accounts for it.

This is the general form of the division of labour the user stated in ¶275: MM
equilibrates everything around the site, restraints preserve the crystal
geometry while it does, and QM/MM captures the site accurately. A bond length is
a QM/MM observable, so it is a QM/MM gate.

**The same reasoning applies to anything measured at a fixed radius on an MM
trajectory**, including the coordination count
([[the-michaelis-complex-is-four-coordinate-at-the-zinc]]).

The experimental half of the comparison is
[[thermolysin-and-cpa-hold-both-oxygens-cis-on-one-zinc]], cited here **for its
Zn-O distances only** - 2.01 and 2.19 Å - and not for the angular window it
was written to license, which ¶276 withdrew.

## Supports
- [[every-acceptance-gate-encodes-a-mechanism]]

## Supported by
- [[thermolysin-and-cpa-hold-both-oxygens-cis-on-one-zinc]]
