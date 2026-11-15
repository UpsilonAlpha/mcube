---
type: observation
status: supported
evidence: experimental
entries: [¶286]
checked: [¶156, ¶283, ¶284, ¶285]
---

# The crystal-derived pose is MX0, an oxacephem, and not meropenem - so it is a template, not a ligand

From 7UHI's own header:

    HETNAM  MX0  (1R,6R,7R)-7-[(2R)-2-CARBOXYPROPANAMIDO]-7-METHOXY-3-METHYL-
                 8-OXO-5-OXA-1-AZABICYCLO[4.2.0]OCT-2-ENE-2-CARBOXYLIC ACID
    FORMUL  MX0  **C13 H16 N2 O8**

and from `l1_in_relax10.pdb`, the transferred pose: **13 C, 2 N, 8 O, no
sulfur**, 23 heavy atoms named C1-C13, N1-N2, O1-O8.

Our ligand is meropenem: **17 C, 3 N, 5 O, 1 S** and 25 H, 51 atoms, named
`C1x`-`C17x`, `S1x`, `N1x`-`N3x`, `O1x`-`O5x`.

**These are different molecules of different classes.** MX0 is a
**5-oxa-1-azabicyclo[4.2.0]octene** - an oxacephem, the moxalactam core - and
meropenem is a carbapenem. Different ring system, different ring size, no
sulfur in MX0, eight oxygens against five.

**What this invalidates and what it does not.** ¶156's transfer script already
renames the transferred residue to `mer` (`--poseres=mer`), and every entry since
has called the result "the crystal-derived pose" of our substrate. That naming is
wrong. What survives:

- **The pharmacophore positions are still a legitimate template.** A carboxylate,
  a beta-lactam carbonyl and a ring nitrogen are common to both classes, and
  those are the only atoms any of the geometric arguments used (O8, O1, N1, C1).
- **The contacts measured on it are measurements of where a bound substrate's
  carboxylate, carbonyl and leaving nitrogen sit in a real B3 complex**, which is
  what they were used for.

**What it does invalidate: any statement that equilibrating "the crystal-derived
pose" means equilibrating our system.** It does not. The MX0 coordinates cannot
be given to our force field, because the molecule is not parameterised and is
not our substrate. Reaching that arrangement with meropenem requires either a
build (mapping meropenem onto the shared atoms) or a restrained pull, and both
are constructions that need deciding rather than assuming
([[no-crystal-derived-pose-has-ever-been-equilibrated]]).

**And it sharpens the selectivity problem rather than solving it.** MX0 is
cephalosporin-class. So 7UHI is a structure of a **cephem** bound to a B3 MBL
with its carboxylate on Zn2 at 2.65 Å - which means the carboxylate-on-metal
arrangement our pose inherits is demonstrably available to the substrate class
that S214D does **not** block. That does not explain the selectivity, but it
does say the explanation cannot be "cephems bind differently at the
carboxylate".

**Why this does not refute [[the-crystal-derived-pose-puts-the-carboxylate-on-the-metal]].**
That node says a carboxylate sits on the metal in this pose, and that is true:
MX0 has one and it does. The node already records that the 1.847 Å contact is
a transfer artefact. What this observation corrects is the pose's **scope** - 
whose molecule it is - not its content, so the edge here is a `Supports` on
the node that tracks the gap between the pose and our system, and the scope
correction is recorded in both bodies.

## Supports
- [[no-crystal-derived-pose-has-ever-been-equilibrated]]
