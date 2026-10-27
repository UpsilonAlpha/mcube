# The question, what answers it, and when an answer counts

Written 2026-10-25 to keep the project on track. Read this before choosing the next calculation;
the decision log (C17, C20, C21) has the evidence behind each line.

## The question (user, 2026-10-23; narrowed 2026-10-25)

**Which ion does the constriction site prefer: Ca²⁺, La³⁺ or Lu³⁺?** Relative, not
publication-grade. **Permeation is out of scope (user, 2026-10-25):** no exit barriers, no
profiles along the pore, no trap-vs-channel question.

The quantity is the exchange free energy at the site:
M₁(site) + M₂(water) → M₂(site) + M₁(water), for Ca→La, Ca→Lu and La→Lu.
It needs no profile along z.

The site is the carboxylate ring Glu156/Glu240/Glu300 near z ≈ 0, where La binds Glu156 (and
Glu240) and Lu binds Glu240 (and Glu220) in the residence runs (C21).

## Order of work

1. **Cluster exchange ΔG at the Glu ring (CPU, cheap).** This is the main calculation.
   - Cut out the ring (Glu156/240/300, plus Glu220 if it is in reach) with backbone caps and the
     ion's first-shell waters, from bound frames of each ion (residence-run lastframes, C21).
   - r2SCAN with def2-TZVP, lcecp-1-TZVP for La/Lu (C19); implicit solvent and thermal corrections;
     machinery as `coord/carboxylate_exchange.py` and `cp_tzvp.py`.
   - Report Ca→La, Ca→Lu and La→Lu for several binding modes (1/2/3 carboxylates, one or both
     O of each). The site's preference is the lowest-free-energy mode per ion.
   - Precedent (checked 2026-10-27, C22 "Literature comparison"): Dudev, Chang & Lim, *JACS* 2005, 127,
     4091 (La³⁺ for Ca²⁺); Nikolova, Kircheva, Dobrev, Angelova & Dudev, *Int. J. Mol. Sci.* 2023, 24,
     6297 (PMC10094714; Ca²⁺ vs La–Lu, ε 4 / 29, SMD), which our protocol follows. Their rule: the
     pocket's net charge decides Ca vs Ln (−4/−3 favour Ln, −1 favours Ca). Their Lu vs La gap is
     small (|ΔΔG| ≤ ~2 kcal/mol).
   - It is QM, so the 12-6-4 problem does not enter. It ignores protein flexibility and the
     barrel/membrane electrostatics beyond the cluster (partly mitigated by the continuum).
2. **Sampling only where it changes the answer.** If the ranking depends on which binding mode the
   site adopts (how many carboxylates, mono- vs bidentate), sample the modes in the real protein.
   Use metadynamics on the carboxylate-O count n_COO at z0 ≈ 0 (the mcube stage being built), with
   integer-state free energies (0/1/2/3), to choose the modes that go into step 1. This is a supporting
   calculation, not the answer.
3. **Experiment for the final word:** Ca/La/Lu binding assays on the protein, or mutations of the
   Glu ring. The genomic context (lanmodulin, XoxF5, PedA1/B/C importer next door) is
   circumstantial evidence for Ln.

## When a model answer counts

1. **The ranking holds across plausible binding modes**, not only one hand-picked geometry.
2. **The ranking holds with a larger cluster or second solvation shell** (or one more ring residue):
   the result must not be a cluster-boundary artefact.
3. **The gap exceeds the model error:**
   - Ca vs Ln: more than ~4 kcal/mol (r2SCAN vs 12-6-4 TZVP+CP scatter −6…+5, C17; continuum
     solvation errors for 2+ vs 3+ ions are of this order too).
   - La vs Lu: more than ~2 kcal/mol (the experimental acetate La−Yb gap is only ~0.2,
     c4-tuning.md).
4. **Calibrated on a known case:** the same protocol on acetate (or a lanmodulin EF-hand site) should
   reproduce the experimental Ln-over-Ca direction (acetate: Ln − Ca ≈ −1.5 kcal/mol, c4-tuning.md).
   **Done 2026-10-25 (C22):** the direction is right, but there is a systematic bias. PCM gives La − Ca −6.0 (exp. −1.5)
   and Lu − La −2.8 (exp. ≈ −0.2). So the protocol over-favours Ln over Ca by ~4.5 kcal/mol and Lu over
   La by ~2.6 **per carboxylate**. Read site results against that offset: Ln over Ca counts only beyond
   ~4.5 × (number of coordinating carboxylates), Lu over La only beyond ~2.6 × that number.

## Status of the answer (2026-10-27)

- **Ca vs Ln:** depends on the effective pocket charge. The stripped sites favour Ca, and the −4 ring
  favours Ln, which matches the literature rule. Settling it needs the protein environment (QM/MM snapshots).
- **La vs Lu: unresolved, pending an SMD check.** Our PCM sites favour Lu by 5–13 kcal/mol even after
  the calibration offset. That is about 10× Nikolova *et al.*, and against lanmodulin's experimental
  light-Ln preference. It is probably an artefact of IEF-PCM with hand-set ion radii. Do not report "Lu over La"
  until SMD single points on our geometries reproduce it (C22).

## Limits to state with any answer

- Static cluster plus continuum: no protein flexibility, no explicit long-range electrostatics.
- Fixed Glu protonation. Protonation of the ring may change when a 3+ ion binds.
- An AF2 monomer with residues 110–130 cropped.
- The biology. Free Ca²⁺ is µM–mM in the environment, free Ln³⁺ nM or less, so preference at the
  site is not the same as occupancy in vivo. Lanthanide-using methylotrophs often take up Ln as
  chelator complexes; not yet checked for this organism.

## What not to drift into

- Permeation: profiles along z, exit barriers, more umbrella windows of any kind.
- Full 2D umbrella sampling, or metadynamics profiles along the pore.
- Publication-grade convergence. The target is a relative ranking with the limits above stated.
