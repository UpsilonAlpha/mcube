---
type: claim
status: open
entries: [¶119, ¶278]
---

# Asp92 OD1 hydrogen-bonds the nucleophile at 2.5-3.2 angstroms in the reactant

**The gate: `Asp92 OD1 - O01, range [2.5, 3.2], required_at: [reactant]`,
sourced "H-bond not coordination; ¶119 measured 2.82-2.94 at 100% frames".**

The measurement is real and the distance is a hydrogen bond rather than
coordination, which the gate says explicitly. **What makes it a hypothesis
rather than a check** is what the contact is for: it casts Asp92 as the general
base orienting or deprotonating the nucleophile, which is one of the competing
roles for that carboxylate, and the one the rest of the map argues against.

The map holds the alternatives and they are not compatible with this one being
required:

- [[the-lysine-holds-asp92-where-a-zn2-would]] and
  [[thesis-lys226-restrains-asp92-in-place-of-zn2]] - Asp92 is held **away** from
  the metal by the lysine, 4.3-5.4 Å out across every converged optimisation.
- [[asp92-is-never-a-zn1-ligand-in-b3-structures]] - and if it is 5 Å from the
  zinc it is not 2.9 Å from a zinc-bound hydroxide.
- ¶119's 100% occupancy was measured on a frame family whose Asp92 position was
  **imposed** rather than found (¶276), so the occupancy is a property of the
  restraint schedule.

So the gate requires, at the reactant, a contact that this project's own
converged geometries do not produce. **Open**, and it should not fail a
structure until the Asp92 question is settled
([[question-how-far-is-asp92-from-zn1-in-wild-type-csr1]]).

## Supports
- [[every-acceptance-gate-encodes-a-mechanism]]

## Requires
- [[question-how-far-is-asp92-from-zn1-in-wild-type-csr1]]
