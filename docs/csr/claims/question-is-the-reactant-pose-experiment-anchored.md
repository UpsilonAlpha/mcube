---
type: question
status: open
entries: [D89, D123, D156, D178, D191, D192, D194, D197, D247]
---

# Is the reactant every QM/MM result is computed on an experiment-anchored Michaelis complex, or a docked pose with one coordinate restrained?

Open, and **this node exists because the map had no way to ask it.** Until
2026-11-09 nothing in the argument map recorded which structure a QM/MM
observation was measured on, so the project ran a 1,500 core-hour def2-SVP leg
on a reactant whose defect was already measured and written down eight days
earlier, and the map could not have flagged it. Provenance edges now exist
([[relax10-holds-the-crystal-geometry-for-10-ns]] is `required by` the QM/MM
observations), and this question is the standing premise they lead to.

**The measurement, D247.** 7UHI's moxalactam transferred into `relax10`'s own
frame, fit RMSD **0.310 Å** on His88 NE2 / His164 NE2 / Asp92 CG — the best fit
this project has achieved, against 0.622 Å into `CSR_wt` and 0.851 Å into
`CSR_ligand`:

| quantity | `relax10` | 7UHI transferred |
| --- | --- | --- |
| ligand centroid to Zn1 | **7.28** | **4.77** |
| **centroid separation** | **6.20** | — |
| β-lactam carbonyl to Zn1 | 2.610 (restrained) | 3.15 |
| C3 carboxylate O to Zn1 | 5.656 / 6.813 | 1.85 / 3.98 |

The nine docked poses spanned 0.98–3.66 Å by centroid (D156); `scan_00` gave
5.15 Å, which D156 called further than any of them. `relax10` is **6.20 Å**.
`eq6` restrained `Zn–O2x` in to ~2.6 Å and `relax10` held it, so the one
coordinate anybody inspects reads correctly while the body sits 2.51 Å further
from the metal than L1's does. **Right coordinate, wrong body.**

**Two arguments that are NOT available here, and must not be re-derived.**

1. **"The substrate dissociated, so the pose is wrong."** Discarded at D89,
   re-derived and discarded at D192, consolidated as discarded at D194, and
   re-derived a fourth time in a draft of
   [[no-usable-michaelis-complex-exists-yet]] on 2026-11-09 (D246). The
   `wt100c`/`prod100c` substrate loss was **one config line** — an Asp92 OD1
   restraint written *instead of* the attack restraint rather than beside it.
2. **"The uncoordinated carbonyl is a regression."** It is a design decision
   (D98: *"Nothing asserts the carbonyl coordinates the zinc"*), and D194
   withdraws "transfer a different pose" as a proposal to undo a decision.

What survives without either: D156's centroid measurement, D197's rejection of
the equilibrated structure on `O1x` occupancy (97.9% unrestrained in the fourth
coordination site, carbonyl inside 2.5 Å in 0.3% of frames), and D247 above.

**What would answer it.** The route of record, D178, builds the Michaelis
complex from **5AXO** — SMB-1 with hydrolysed meropenem, our own substrate in a
B3 enzyme, measured — by re-forming the `C4x–N1x` bond:

    prod100c -> prodpost gates -> prodqm -> prodopt -> recombine
    (100 ns)    (shell + H)       (50 ps)   (opt)      (re-form C4x-N1x)

so the structure `recombine` ends on **is** the complex, anchored through the
product and on the product's own topology, hence a valid NEB partner by
construction. D178 scoped 7UHI as the check; **D192 made it the source**, and
that part of D192 stands. D178 also rejected the alternative — grafting
meropenem's substituents onto 7UHI's transferred moxalactam core — as "strictly
harder and less direct than re-forming a bond in a molecule that is already
meropenem", so that route is not to be re-proposed (D246).

The route is half built: `CSR_prod` transferred LMP into our frame and mapped it
onto our atom names by graph isomorphism (D159), and `prod100c` then NaN'd at
22.68 ns and failed the gate (D191). **One practical obstacle on top of the
NaN:** D247's transfer lands L1's carboxylate oxygen 1.85 Å from our zinc, a
clash, so the clash check D156 §5.2 required is not optional.

**What does not depend on the answer.** The method-comparison results are
common-mode — both Hamiltonians saw the same structure — so the matched
DFT-versus-GFN2 agreement on Lys226 and disagreement on the carboxylate (D245)
hold regardless. The Asp92 and Lys226 distances are good measurements of the
**resting** site, matching substrate-free 6DQH to within 0.13 Å on four of five
coordinates; they cannot be captioned as Michaelis-complex geometry. **The
barrier is what this question gates**, and no barrier has been computed
([[no-reaction-coordinate-has-ever-been-computed]]).

## Required by
- [[question-does-restraining-asp92-lower-the-barrier]]
- [[thesis-mechanism-differs-from-canonical-b3]]

## Refutes
- [[carbonyl-coordinates-zn1-in-the-michaelis-complex]]
