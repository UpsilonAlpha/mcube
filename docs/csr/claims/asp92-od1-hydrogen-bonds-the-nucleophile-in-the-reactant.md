---
type: claim
status: open
entries: [¶119, ¶278, ¶302]
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

**¶301: THE WINDOW NOW HAS EXPERIMENTAL BACKING, where before it had only ¶119's
own frames.** Asp120 carboxylate oxygen to the nucleophile oxygen, measured
from deposited coordinates:

    L1 1SML   (resting di-Zn)   2.83 / 2.97
    SMB-1 5AXO (product)        2.76 / 3.07
    GOB-18 5K0W (di-Zn)         3.61 / 3.86

So **2.8-3.9 Å across three B3 enzymes**, and in each the aspartate is also a
metal ligand (OD2-Zn 2.07-2.19). The gate's 2.5-3.2 Å window sits inside that
range at the tight end. The contact is real and experimental; it is the one gate
on this map whose provenance improved.

**Our structures are 7.2 Å, and ¶301 locates why.** Measured on
`aspxtb8_nowat.pdb`: Zn-OD1 **5.830**, Zn-OD2 **5.092**, OD1-O01 **7.211**. The
aspartate is not ligating anything - in every reference enzyme it is held by
Zn2, which we do not have
([[asp92-is-a-zn2-ligand-in-b3-structures]]). So the 7.2 Å is a consequence of
the vacant subsite and not an independent pose error, and this gate cannot be
satisfied while the aspartate is unanchored.

**¶302 WITHDRAWS THE CONCLUSION DRAWN ABOVE, while keeping the measurements.**
The experimental 2.8-3.9 Å is a **bridging**-hydroxide geometry: in 1SML, 5AXO
and 5K0W the nucleophile sits BETWEEN two metals (Zn1-OH 1.82-2.06, Zn2-OH
2.06-2.07) and Asp120 is close to it because Asp120 **ligates the second metal**.
We have one metal and a **terminal** hydroxide, so there is no bridge to occupy
and no second metal for the aspartate to hold. Asp92 at 5.1-5.8 Å from the metal
and 7.2 Å from the hydroxide is therefore the expected consequence of the vacant
subsite, not an independent pose error - and this gate, as written, asks a
mono-zinc structure to reproduce a di-zinc second shell. **It should not fail a
structure on that basis.**

**A rival diagnosis, tested and rejected.** It was put to this project that 7.2 Å
is what a hydroxide placed on the solvent-exposed face of M1, *anti* to Asp92,
would give - the right ligand on the wrong face. Measured, the angles are
**OD1-Zn-O01 135.8 deg** and **OD2-Zn-O01 114.1 deg**: tilted away from the
aspartate face, but nowhere near anti (180 deg). The hydroxide is not on the
wrong face; the aspartate is simply 5-5.8 Å out.

## Supports
- [[every-acceptance-gate-encodes-a-mechanism]]## Supports
- [[every-acceptance-gate-encodes-a-mechanism]]

## Requires
- [[question-how-far-is-asp92-from-zn1-in-wild-type-csr1]]
