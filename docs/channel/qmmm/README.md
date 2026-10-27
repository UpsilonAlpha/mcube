# QM/MM vertical swaps: ion exchange at the Glu ring with the real environment

The question and the criteria are in [../question-and-criteria.md](../question-and-criteria.md).
The cluster exchange energies (decision log C22, `../coord/site_exchange.py`) flip between
Ca-selective and Ln-selective depending on whether the uncoordinated ring Glu are in the cluster.
The PROPKA and environment check (`../coord/ring_pka.py`) found all four Glu charged. But Glu300 is
salt-bridged to Arg302/Arg366 about 12 Å from the ion, Glu220 is partly paired with Arg202 and has
Na⁺ next to it, and every carboxylate has 3–7 waters. A cluster in a continuum cannot represent
that screening. Electrostatic embedding of the whole protein, membrane, water and ions can.

## Method

For M1 → M2 in Ca → La, Ca → Lu and La → Lu:

    ΔΔE = ⟨E(M2@site) − E(M1@site)⟩_site frames − ⟨E(M2@water) − E(M1@water)⟩_water frames

- **Vertical swap.** Each difference is taken on one MM frame. Only the ion changes: residue
  name, atom name, element, QM charge, and its MM Lennard-Jones terms (`metal_parameters`). The
  ion stays at the same position.
- **Water reference.** The same swap on frames of the ion in a bulk TIP3P box, with the same
  protocol. The 2+ vs 3+ long-range and Born-type terms then largely cancel. Negative ΔΔE means
  the site prefers M2.
- **Two-sided.** Site frames come from the La-bound and the Lu-bound residence runs. Water frames
  come from La and Lu bulk boxes. Every ion goes into every frame. `analyse` reports each source
  separately and pooled. If the per-source rows disagree, the answer depends on which ion the
  frames were sampled with.
- **QM region.** Fixed once per frame, from the source ion's coordination, by mcube's rule
  (`preprocess.coordination_shell`): the closest N/O/S donor of each residue within
  `--qm-cutoff` (3.3 Å) + `--qm-margin` (0), with the minimum image. The same `qm_residues` goes
  into all three ions' configs.
  - Coordinating amino acids are cut at CA–CB, as mcube does.
  - Waters come in whole.
  - Lipids are refused.
  - A Na⁺/Cl⁻ in reach stays MM unless `--qm-counterions` is given. The lcecp basis files have no
    Na or Cl (`../dft/bse_basis.py` default elements).
  - The uncoordinated ring Glu, Arg302/366, Arg202, His222 and the Na⁺ stay MM. That is the
    point: they act through the embedding.
- **Per job** (one frame, one ion), mcube's own stages on one `config.yaml`:
  1. `preprocess`: `skip_solvate`, so there is no solvation and no neutralisation, and the
     frame's MM environment is the same for every ion. Its detection only reports; it is checked
     against `qm_residues` in `analyse`.
  2. `qmmm`: no QM/MM MD (see the exception below). A local optimisation from the MM frame, with
     every residue within `--actradius` (8 Å) of the ion free and MM waters rigid.
     `opt_maxiter` 150. A run that does not converge keeps its last geometry, and `analyse` leaves
     that frame's pairs out unless `--include-unconverged` is given.
  3. `sp`: the single point on `QM_MM_optimized.pdb`.
  - The charge of the QM region is the ion's charge plus the formal charge of the QM side chains
    (Glu −1). So it is +1/+2 at a site with one Glu and +2/+3 in water. That is the only
    difference between the ions' inputs besides the ion itself.
- **Level.**
  - Optimisation: r2SCAN/def2-SVP, La/Lu on lcecp-1-SVP (`qm_engine: dft`, with
    `dft_basis_file`/`dft_ecp_file` = `$S/channel/dft/def2-svp_lcecp1.nw`).
  - Single point: r2SCAN/def2-TZVP, La/Lu on lcecp-1-TZVP (`sp_basis_file`/`sp_ecp_file` =
    `def2-tzvp_lcecp1.nw`).
  - These are the cluster calculation's levels and files (C19, C22), so the two can be compared.
  - Why DFT and not xTB for the optimisation: GFN2-xTB is what the bulk Ln³⁺ runs lost SCC
    convergence with (F7). The question also turns on the ion–carboxylate geometry, so the
    geometry should be at the level the energy is computed at. mcube's `qm_engine: dft` takes NWChem
    basis + ECP files, so nothing new is needed. `--opt-engine xtb` is there as a cheaper fallback.
  - QM gas phase with electrostatic embedding: ASH's PySCF interface ignores a continuum when
    point charges are present.
  - No dispersion correction, as in the cluster runs. The ion's own dispersion with its shell
    differs between ions, so state this with any number.
- **MM metal parameters.** `metal_parameters: "1264"`, the set the source runs were sampled with.
  It follows the ion: `ions.apply_metal_parameters` picks the 12-6 terms by element and charge
  (decision log F1, F2). The ion is a QM atom, so its C4 term is off. The default and Ch-BE/scaled
  C4 choices of the source runs therefore do not enter the energies, only the frames.

### The deliberate exception: no QM/MM MD before the optimisation

mcube's usual order is QM/MM MD, then optimisation. Here each swapped system gets only a local
optimisation from the shared MM frame (`skip_equilibrate: true`, set in `exchange.py` with a
comment where it is written, and noted at the top of every config). A QM/MM MD per swapped system
would let the three environments drift apart, and that would destroy the pairing that makes
E(M2) − E(M1) on one frame meaningful. Sampling comes from the frames, not from QM/MM dynamics.
The cost is that each energy is a single local minimum near the MM structure. The spread across
frames shows how much that matters.

## Frames

- **Site.** 5 frames per source run, 1.5 ns apart, counted back from the last frame:
  `--n 5 --spacing-ns 1.5 --start 0.4`.
  - Five frames 1.5 ns apart span 6 ns, which does not fit in the second half of a 10 ns run.
  - So the default starts at 40 % of the run, giving frames at 4, 5.5, 7, 8.5 and 10 ns.
  - `--start 0.5 --n 4` keeps to the second half with four frames.
  - The frame spacing in time comes from `cv.csv` (as in `residence.py`).
- **Site sources.**
  - `$S/channel/system/LA/umbrella_res_035`: La, default C4. La sits on Glu156 throughout and on
    Glu240 after 3.6 ns.
  - `$S/channel/system/LU/umbrella_res_chbe_035`: Lu, scaled C4. Lu sits on Glu240 throughout,
    with Glu220 later.
  - Frames are imaged about the protein chains before writing, because ASH applies no periodic
    images to the QM region or its point charges.
- **Water.** The bulk boxes of `../bulk/bulk_ion.py` keep only their final structure. So they are
  rerun with the new `--dcd-every` option, which writes full-system `frames.dcd`:
  - 2 ns at 310 K (the site runs' temperature) with the 1264 set, a frame every 50 ps.
  - 5 frames 0.25 ns apart from the second half. Water around an ion decorrelates in tens of
    ps, so 0.25 ns is ample.
  - The frames are imaged about the ion.

## Jobs and cost

- **Jobs.**
  - Site: 2 sources × 5 frames × 3 ions = 30 jobs on ~64k-atom frames.
  - Water: 2 × 5 × 3 = 30 jobs on ~4k-atom boxes.
  - 32 CPUs and 96 GB each.
- **QM size.** About 35–45 atoms at a site (the ion, 1–2 Glu side chains, 4–7 waters). About
  25–30 in water (the ion and 8–9 waters).
- **Estimated cost.** These are not measured: no DFT/MM gradient with 64k point charges has been
  timed for this system. Measure on the first job.
  - Site job:
    - roughly 2–4 min per r2SCAN/SVP QM/MM gradient on 32 cores (the 64k point charges cost
      about as much as the SCF);
    - 60–150 optimisation steps;
    - plus a TZVP single point of 0.5–1 h.
    - Total ≈ 4–8 h, or ≈ 130–260 core-h.
  - Water job:
    - ≈ 1 min per step, plus 20 min for the single point;
    - ≈ 1.5–3 h, or ≈ 50–100 core-h.
  - Bulk MM frame runs: 2 × ~3 h on 8 CPUs ≈ 50 core-h.
  - **Total ≈ 6,000–11,000 core-h**, mostly the site optimisations. With 10 running at once, that
    is about 1–2 days of wall time.
- **Smoke test first.** Submit one site job and one water job. After an hour, compare the
  number of optimisation steps in `<job>/job.log` against the time. If a site step takes more
  than ~5 min, consider `--actradius 6`, or `--opt-engine xtb` with the TZVP single point
  unchanged.

## Commands on Bunya

```bash
S=/scratch/user/uquadhik
Q=$S/mcube/docs/channel/qmmm
# $S/mcube must have this commit (branch umbrella): exchange.py, and bulk_ion.py's --dcd-every

# 1. Bulk-water frames (MM only, ~3 h each). New run directories, so the finished bulk runs stay as they are.
cd $S/channel/bulk
for ION in LA LU; do
  sbatch -J bulkf_$ION --export=ALL,RUN=$S/channel/bulk/${ION}_1264_frames,NS=2,SKIP_QM=1,MMOPTS="--temp 310 --dcd-every 50" \
      $S/mcube/docs/channel/bulk/bulk.slurm $ION 1264
done

# 2. When both have finished: extract, swap and write the job directories (a few minutes; ASH env python)
module load miniforge/26.1.0-0; source $ROOTMINIFORGE/etc/profile.d/conda.sh; conda activate ASH
export PYTHONPATH=$S/mcube/src
$CONDA_PREFIX/bin/python $Q/exchange.py frames --out $S/channel/qmmm \
    --site LA=$S/channel/system/LA/umbrella_res_035 --site LU=$S/channel/system/LU/umbrella_res_chbe_035 \
    --water LA=$S/channel/bulk/LA_1264_frames --water LU=$S/channel/bulk/LU_1264_frames \
    --basis-svp $S/channel/dft/def2-svp_lcecp1.nw --basis-tzvp $S/channel/dft/def2-tzvp_lcecp1.nw
#    (the site part can run before the bulk runs finish: leave out the --water options, rerun later with them;
#     existing job configs are kept unless --force)

# 3. Smoke test: one site job and one water job, then read <job>/job.log
sbatch -J qx_smoke --export=ALL,JOB=$(sed -n 1p $S/channel/qmmm/jobs_site.txt) $Q/exchange.slurm
sbatch -J qx_smoke -t 8:00:00 --export=ALL,JOB=$(sed -n 1p $S/channel/qmmm/jobs_water.txt) $Q/exchange.slurm

# 4. Everything (finished jobs exit at once, so the smoke-test jobs are not redone)
bash $Q/submit.sh $S/channel/qmmm both 10

# 5. Collect: per-frame dE, mean and sd per source and pooled, and ddE for the three exchanges
$CONDA_PREFIX/bin/python $Q/exchange.py analyse --out $S/channel/qmmm     # writes $S/channel/qmmm/report.json
```

## Not verified off Bunya

`frames` and `analyse` were run locally on a synthetic system: villin with La³⁺ on a Glu, every
atom wrapped into the box to test the imaging, and a Lu³⁺ water box. The written frames build
with OpenMM's amber14 force field for all three ions. The QM region and charge are identical
across the ions of a frame, and the `analyse` arithmetic is right. Not run here, for want of ASH
and PySCF: the `preprocess` → `qmmm` → `sp` chain itself, the lcecp basis files through ASH's
`basis_file` path in a QM/MM optimisation, the 64k-point-charge timings, and the real umbrella
start.pdb (lipid CONECT records; `--bonds-from` covers a file without them).
