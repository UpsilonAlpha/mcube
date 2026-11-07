---
type: observation
evidence: ours
status: supported
entries: [D235, D236]
checked: [D236]
---

# Through three wb97x-v/def2-SVP QM/MM optimisation steps from the crystal-consistent frame, Asp92 moves away from Zn1 rather than toward it, and the Lys226 salt bridge tightens onto the crystal value

`aspdft`, job 32453561: wb97x-v/def2-SVP QM/MM, 124 QM atoms, active region cut
to the QM atoms, **nothing constrained**, started from
`md_relax10_results/NVTtrajectory_lastframe.pdb`
([[relax10-holds-the-crystal-geometry-for-10-ns]]).

| | step 0 | step 1 | step 2 | 6DQH |
| --- | --- | --- | --- | --- |
| Zn–OD2 | 5.153 | 5.106 | **5.087** | 5.14 |
| Zn–OD1 | 4.327 | 4.509 | **4.544** | 4.28 |
| Zn–CG | 5.192 | 5.216 | **5.198** | 5.15 |
| Zn–NZ | 5.276 | 5.287 | **5.269** | 5.27 |
| NZ–OD2 | 2.880 | 2.799 | **2.766** | 2.79 |
| His88 NE2 | 2.178 | 2.216 | 2.222 | |
| His164 NE2 | 2.126 | 2.118 | 2.114 | |

Energy falls monotonically, −5201.1362 → −5201.3065 → −5201.3341 → −5201.3519
Eh; trust radius steady at 1.414e-01; step quality 0.645 then 0.871; gradient
rms 1.199e-02 → 6.038e-03. Zero `SVD Error` lines in the `.out` and zero
`LinAlgError` in the `.err`, counted separately because conflating those two
produced a wrong conclusion here once.

**The atom mapping was checked before the numbers were used.** The QM-region
geometry blocks were mapped through `results.qm_atoms` and the input PDB, and
the first block reproduces `relax10`'s independently recorded start to 0.003 Å.

**What this does not establish, stated because the temptation is large.** Three
steps is a *direction*, not a result: total displacement rms 2.7e-02 Å, gradient
two orders of magnitude off the convergence criterion, no stationary point
reached. And it is **not a matched comparison with GFN2**: that was 50 ps of
300 K dynamics, this is a zero-temperature optimisation with no thermal energy
to cross a barrier it might find. The matched test — a GFN2 optimisation from
this same frame — has not been run.

**What it does establish** is narrow and new: the DFT gradient at a geometry four
crystals support does not point the carboxylate at the metal, where GFN2's
dynamics took it from 4.34 to 2.01 Å within 8 ps
([[aspqm2-chunk1-asp92-chelates-zn1-under-gfn2]]). Until now the case that this
was the Hamiltonian's doing rested on a cluster benchmark and on the crystals;
this is the first measurement of the DFT gradient on the protein itself.

## Supports
- [[gfn2-over-coordinates-zn1-from-a-crystal-consistent-start]]
