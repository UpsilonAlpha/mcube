---
type: observation
evidence: ours
status: supported
entries: [D235, D236, D240, D242, D244, D245]
checked: [D236, D240, D242, D244, D245]
---

# Through twenty-five wb97x-v/def2-SVP QM/MM optimisation steps from the crystal-consistent frame, the gradient holds Asp92 at the crystal distance from Zn1 rather than pulling it on, and moves Lys226 onto the metal instead

`aspdft`, job 32453561: wb97x-v/def2-SVP QM/MM, 124 QM atoms, active region cut
to the QM atoms, **nothing constrained**, started from
`md_relax10_results/NVTtrajectory_lastframe.pdb`
([[relax10-holds-the-crystal-geometry-for-10-ns]]).

| | step 0 | step 4 | step 6 | step 10 | step 13 | 6DQH |
| --- | --- | --- | --- | --- | --- | --- |
| Zn–OD2 | 5.153 | 5.087 | 5.073 | 5.143 | **5.134** | 5.14 |
| Zn–OD1 | 4.327 | 4.544 | 4.536 | 4.527 | **4.546** | 4.28 |
| Zn–CG | 5.192 | 5.198 | 5.206 | 5.224 | **5.237** | 5.15 |
| Zn–NZ | 5.276 | 5.269 | 5.213 | 5.145 | **5.131** | 5.27 |
| NZ–OD2 | 2.880 | 2.766 | 2.699 | 2.740 | **2.713** | 2.79 |
| His88 NE2 | 2.178 | 2.222 | 2.235 | 2.180 | 2.162 | |
| His164 NE2 | 2.126 | 2.114 | 2.097 | 2.085 | 2.073 | |

**Zn–OD2 turns over at step 6 and comes back.** The first three steps read as
an outward drift; by step 13 the carboxylate sits 0.006 Å from where 6DQH puts
it, flat to ±0.005 Å over the last four steps, with Zn–CG and Zn–OD1 flat on
the same span. The gradient neither pushes Asp92 out nor pulls it in.

**Lys226 moved onto the metal and then STOPPED (D244), so the distance is a
number and not a direction.** Zn–NZ fell 5.276 → 5.120 by step 16, came back,
and reads **5.140 at steps 23, 24 and 25** — flat to the third decimal. The
ammonium takes up position in the vacant Zn2 subsite, unrestrained, holding
its salt bridge at 2.73: the motion the thesis asserts, 0.13 Å closer to the
zinc than 6DQH's 5.27. **Letting the job run past step 13 is what turned this
into a measurement** — at step 13 it read 5.131 and was still falling.

| | step 13 | step 16 | step 22 | step 25 | 6DQH |
| --- | --- | --- | --- | --- | --- |
| Zn–OD2 | 5.134 | 5.133 | 5.117 | **5.123** | 5.14 |
| Zn–OD1 | 4.546 | 4.553 | 4.535 | **4.534** | 4.28 |
| Zn–CG | 5.237 | 5.237 | 5.220 | **5.223** | 5.15 |
| Zn–NZ | 5.131 | 5.120 | 5.138 | **5.140** | 5.27 |
| NZ–OD2 | 2.713 | 2.726 | 2.735 | **2.729** | 2.79 |

**Four of those five are within 0.14 Å of the crystal**, the fifth (OD1) 0.25 Å
out, from an optimisation with nothing constrained. Through step 24 the energy
is monotonic at every step (−5201.1362 → −5201.3906 Eh), quality 1.15–1.35,
trust growing again to 1.198e−01, gradient rms **3.820e−04**. Zero `SVD Error`
and zero `LinAlgError`. The step rate did not keep improving — 124/132/116 min
against 91 at step 13 — so the 48 h wall lands near step 27 and writes no
`QM_MM_optimized.pdb`.

Energy falls monotonically at every one of thirteen steps, −5201.1362 →
−5201.3871 Eh; step quality 1.07–1.20 over the last six; trust radius grown
from 1.000e-01 and capped at 3.000e-01; gradient rms 3.656e-02 → 1.339e-03,
the last five steps decaying only ~8% each, so the 1e-04 criterion is ~30
steps away and the 48 h wall (≈29 steps) will not reach it. Zero `SVD Error` lines in the `.out` and zero
`LinAlgError` in the `.err`, counted separately because conflating those two
produced a wrong conclusion here once.

**The atom mapping was checked before the numbers were used.** The QM-region
geometry blocks were mapped through `results.qm_atoms` and the input PDB, and
the first block reproduces `relax10`'s independently recorded start to 0.003 Å.

**What this does not establish, stated because the temptation is large.** No
stationary point has been reached: displacement rms is 1.2e-02 Å and the
gradient is an order of magnitude off the criterion. The distances are settled
where the gradient is not, which is weaker than convergence and is the reason
the first three steps' reading had to be withdrawn — a trend over three steps
reversed at the sixth. And it is **not a matched comparison with GFN2**: that was 50 ps of
300 K dynamics, this is a zero-temperature optimisation with no thermal energy
to cross a barrier it might find. The matched test — a GFN2 optimisation from
this same frame — has not been run.

**The inference first drawn from this is RETRACTED (D240).** D236 read it as the
first direct evidence that GFN2's chelate is the Hamiltonian's rather than the
enzyme's. **GFN2 does the same thing and more of it**: its own converged minimum
from this identical start moves Asp92 0.23 Å further out where DFT moved it 0.07,
and tightens the same salt bridge to the same place
([[gfn2-minimum-puts-asp92-further-out-than-the-crystal]]). So these three steps
say nothing about DFT versus GFN2. The numbers stand; the inference does not.

**What survives** is that neither method's gradient, from a crystal-consistent
start with the environment frozen, points the carboxylate at the metal - which
is about where the chelate is *not*, and is the half of
[[gfn2-over-coordinates-zn1-from-a-crystal-consistent-start]] that locates the
bias in the barrier rather than in the minimum.

**THE CROSS-METHOD SPREAD IS UNEVEN, AND ONLY HALF OF IT IS SAFE (D245).**
This run has `actradius: 0` — 124 active atoms, MM environment frozen — so the
**matched** comparison is with `aspxtb` (also 124), not with `aspxtb8` (479).
There is no DFT run with an 8 Å active region.

| pair | DFT/124 | GFN2/124 | GFN2/479 | 6DQH | matched diff |
| --- | --- | --- | --- | --- | --- |
| Zn–NZ | 5.140 | 5.195 | 5.044 | 5.27 | **0.055** |
| NZ–OD2 | 2.729 | 2.718 | 2.749 | 2.79 | **0.011** |
| Zn–OD2 | 5.123 | 5.383 | 5.092 | 5.14 | 0.260 |
| Zn–OD1 | 4.534 | 4.952 | 5.830 | 4.28 | **0.418** |
| Zn–CG | 5.223 | 5.564 | 5.777 | 5.15 | 0.341 |
| His88 NE2 | 2.162 | 2.083 | 2.063 | | 0.079 |
| His164 NE2 | 2.073 | 2.010 | 1.997 | | 0.063 |

**Method-independent:** where Lys226 sits, and that Asp92 is not on the metal.
The salt bridge spans 0.031 Å across all three setups and falls short of the
crystal in every one. **Not method-independent:** how Asp92 is *turned* —
Zn–OD1 differs by 0.418 Å matched and spans 1.296 Å across all three, and that
rotamer is the coordinate the mechanism runs through. Any statement about
Asp92's orientation must name its run.

**GFN2 over-binds the zinc and still keeps Asp92 out.** It is 0.079 and
0.063 Å shorter than DFT on the two histidine nitrogens — the known direction
of GFN2 error, and the same direction Mehmood & Kulik report for too-small QM
regions. The bias and the result point opposite ways, which is worth more here
than two methods agreeing.

## Supports
- [[gfn2-over-coordinates-zn1-from-a-crystal-consistent-start]]

## Requires
- [[relax10-holds-the-crystal-geometry-for-10-ns]]
- [[question-is-the-reactant-pose-experiment-anchored]]
