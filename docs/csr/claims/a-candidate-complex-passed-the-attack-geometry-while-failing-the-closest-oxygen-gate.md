---
type: observation
status: supported
evidence: ours
entries: [¶125, ¶126]
checked: [¶125, ¶126, ¶278]
---

# A candidate complex satisfied the attack geometry at 3.17 angstroms and 106 degrees while failing the closest-oxygen gate

¶125 wrote the criterion *"the closest substrate oxygen must be the beta-lactam
carbonyl"* as though it were a fact about Michaelis complexes. ¶126 then
measured the grafted pose-2 structure:

    Zn - O1x (6-alpha-hydroxyethyl)   2.66     <- the metal holds THIS
    Zn - O2x (beta-lactam carbonyl)   4.50     <- points away
    O01 - C4x (attack)                3.17     <- unrestrained
    O01 - C4x - O2x                   106°      <- inside the 95-120 window
    Asp92 OD1 - O01                   2.87

**So the gate rejects this structure, and the structure satisfies the attack
geometry the gate exists to find** - at 106 degrees, which the preferred arm
never reached in 791 frames (¶34), with nothing restraining that coordinate.

This is the concrete case that shows an acceptance criterion can encode a
mechanism rather than test a structure: the rejected arrangement is a competing
hypothesis ([[question-does-the-substituent-on-metal-route-survive]]), and the
criterion assumed the other one.

## Supports
- [[every-acceptance-gate-encodes-a-mechanism]]
- [[question-does-the-substituent-on-metal-route-survive]]
