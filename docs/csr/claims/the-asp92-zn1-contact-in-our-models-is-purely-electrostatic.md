---
type: observation
evidence: ours
status: supported
entries: [¶152, ¶153, ¶154]
checked: [¶153]
---

# Protonating Asp92 makes it leave Zn1 on its own and restores the crystal's four-coordinate sphere, so the Asp92–Zn1 contact in our models is driven by charge rather than by a coordination preference

`qmmm_asph`, job 31186037: free GFN2-xTB/MM optimisation from
`scan_08_A92_ASH.pdb` — `scan_08` with Asp92 protonated (HD2 on OD2, anti, in
plane, by `site.py protonate`) — **converged in 12 min 03 s, charge 3, nothing
restrained on Asp92 or on anything else.**

| ligand | protonated start | optimised |
| --- | --- | --- |
| **Asp92 OD1** | **1.90** | **2.99** |
| Asp92 OD2 | 3.43 | 4.73 |
| Asp92 CG | 2.95 | 3.81 |
| hydroxide O01 | 1.94 | **1.80** |
| His88 NE2 | 2.04 | 2.00 |
| His164 NE2 | 2.19 | 2.02 |
| HOH485 O | 2.19 | 2.08 |
| **ligands within 2.8 Å** | **5** | **4** |

**Remove the charge and the carboxylate leaves unrestrained**, from 1.90 to
2.99 Å, and the three remaining protein/water ligands all *tighten* — the
four-coordinate sphere 6NC5 measures (His88 2.22, His164 2.09, HOH485 1.90,
HOH493 2.14) comes back. So in our models the Asp92–Zn1 contact is
**electrostatic**: it is the anionic carboxylate against the dication, not a
coordination preference the Hamiltonian insists on. That is what explains
[[asp92-is-a-zn1-ligand-in-our-models]] rather than merely restating it.

**What this does NOT settle, and the reason this node is not a claim about the
resting state.** ¶153 is explicit that `asph` is **"a diagnostic, not a
candidate reactant."** Whether Asp92 is protonated in the resting enzyme is
tied to [[question-what-is-lys226-pka]], which is open: the pair is either
(Lys226⁺ / Asp92⁻) or (Lys226⁰ / Asp92H), and 6DQH's 2.79 Å salt bridge is
measured on the first reading. ¶153 also corrected its own earlier argument
here — it had reasoned that a buried carboxylate 1.9 Å from a dication is "the
textbook environment for an upshifted pKa", which is **backwards**, since a
nearby cation stabilises the conjugate base and so pushes the pKa *down*.

**Consequence for choosing a starting structure (¶258).** `qmmm_asph`'s
optimised geometry is the only four-coordinate QM/MM sphere in the project that
is not on a driven coordinate, which makes it tempting as a build source. It
should not be used as one: it is a different protonation state, so a build from
it would commit an arm to the unsettled half of
[[question-what-is-lys226-pka]] without saying so.

## Supports
- [[asp92-is-a-zn1-ligand-in-our-models]]
- [[the-lysine-holds-asp92-where-a-zn2-would]]

## Requires
- [[question-what-is-lys226-pka]]
