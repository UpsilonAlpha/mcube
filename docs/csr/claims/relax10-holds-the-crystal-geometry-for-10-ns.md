---
type: observation
evidence: ours
status: supported
entries: [¶228, ¶231]
checked: [¶231]
---

# A 10 ns MM trajectory holds the 6DQH Asp92/Lys226 geometry, and the carboxylate's rotamer matches the crystal without being restrained to

`md_relax10` (job 32405173), 10 ns NVT through NPT from `asppull3`'s last
frame, `gentle: false`, six restraints unchanged. Measured on its last frame:

| pair | restrained at | measured | dev |
| --- | --- | --- | --- |
| Zn–OD2 | 5.14 | 5.15 | +0.01 |
| Zn–NZ | 5.27 | 5.28 | +0.01 |
| NZ–OD2 | 2.79 | 2.88 | +0.09 |
| O01–C4x | 3.47 | 3.48 | +0.01 |
| Zn–O2x | 2.60 | 2.61 | +0.01 |
| Zn–O1x | 3.82 | 3.86 | +0.04 |
| **Zn–OD1** | **not restrained** | **4.33** | 6DQH: **4.28** |
| **Zn–CG** | **not restrained** | **5.19** | 6DQH: **5.15** |

Convergence: CA RMSD's last four deciles 1.41 / 1.38 / 1.38 / 1.40, spread
0.03 Å; T 300.292 ± 1.399, drift +0.094 K; site CA RMSF 0.30–0.40 Å. Zn1 is
five-coordinate at 2.50 Å — HOH3095 2.02, hydroxide 2.04, HOH3171 2.05,
His164 2.13, His88 2.18 — with Asp92 far outside it.

**What this is worth, stated carefully, because it is `ours` evidence and
¶206 ranks that last.** The honest claim is about the **rotamer**, not the
distance. Zn–OD2, Zn–NZ and NZ–OD2 are all restrained, which pins OD2 relative
to both the metal and the ammonium; OD1 and CG are then free only to rotate
about the Cβ–Cγ axis, and that rotation can swing Zn–OD1 by more than a
ångström in either direction. So 4.33 Å is **not** an independent measurement
of a crystal distance — it is a measurement that, given where OD2 is held, the
carboxylate chooses the crystal's orientation rather than turning its other
oxygen toward the metal. That is the thing ¶120/¶130 got wrong twice by
asserting an orientation with a hand-written restraint, and it is what
`asppull3` first showed over 2 ns (¶224) and this run sustains over ten.

**It does not support [[the-lysine-holds-asp92-where-a-zn2-would]].** That node
is `experimental`, resting on four crystals of three constructs, and hanging a
restrained MM trajectory under it would weaken it rather than strengthen it.
This node stands beside it as the statement that a *startable* structure
reproducing it exists.

**It does not resolve [[no-usable-michaelis-complex-exists-yet]]** either. That
node's blocker is the **pose** — the docked carbonyl at 5.43 Å against
7UHI/7UHJ's 2.53/2.63 Å — and this run inherits the same pose it was built
from. A stable starting structure is not a correct one.

**The five-coordinate sphere is unresolved and is MM's to leave alone.**
[[csr1-zn1-site-is-h-r-h]] implies four donors — two histidines, hydroxide,
one water — and two bulk waters are present instead. MM cannot decide
coordination, so whether HOH3095 leaves is the QM/MM question recorded in
[[question-does-gfn2-over-coordinate-the-zinc]].

**WHAT THIS FRAME IS NOT, MEASURED 2026-11-09 (¶247).** Its ligand's
heavy-atom centroid is **6.20 Å** from 7UHI's transferred moxalactam — further
than all nine docked poses (0.98–3.66 Å) and further than `scan_00` (5.15 Å).
`eq6` restrained `Zn–O2x` in to ~2.6 Å and this run held it, so the carbonyl
reads correctly while the body sits 2.51 Å further from the metal than L1's
does, with the C3 carboxylate 5.7–6.8 Å out. **Everything measured on this
frame inherits that**, which is why the QM/MM observations now declare it as a
`requires` premise and why
[[question-is-the-reactant-pose-experiment-anchored]] is open. The Asp92 and
Lys226 distances below are measurements of the **resting** site — they match
substrate-free 6DQH to 0.13 Å — and not of a Michaelis complex.

## Supports
- [[csr1-zn1-site-is-h-r-h]]
