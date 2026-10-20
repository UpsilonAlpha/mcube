# Running mcube efficiently on Bunya

Notes from setting up and running the CSR (meropenem) reactant/product systems in October 2026. The numbers come from `seff`/`sacct` on real jobs unless marked as estimates.

## Quick recipe

```bash
#SBATCH --partition=general --qos=normal --account=a_imb_ccbcd
#SBATCH --cpus-per-task=8          # = 4 physical cores on SMT nodes; leave `cores:` unset (mcube uses 4) or set it to 4
#SBATCH --mem=16G                  # see the memory table below
#SBATCH --time=12:00:00            # shorter limits start sooner (backfill)

set -euo pipefail
ulimit -s unlimited                # xTB recommendation, avoids stack segfaults
export OMP_STACKSIZE=1G
# Threads depend on what dominates the step - see "Resources by job type".
# For an xTB stage: half the allocated CPUs (physical cores), bound.
export OMP_NUM_THREADS=$(( SLURM_CPUS_PER_TASK / 2 ))
export OMP_PROC_BIND=close
export OMP_PLACES=cores
# For a DFT stage, use OMP_NUM_THREADS=1 instead, with no binding, and let
# PySCF set its own thread count from dft_cores / sp_cores.

module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH                 # mcube is pip-installed (editable) in this env

(cd preprocess_results && mcube preprocess ../config.yaml)
```

- Run each stage in its own directory. Stages write many files to the current directory and write results back into `config.yaml`.
- Don't run two stages at the same time against one `config.yaml`. Each stage reads, edits and rewrites the whole file, so concurrent runs can silently overwrite each other's results. Run them one after another in a single job, or give each its own config.
- After a job finishes, run `seff <jobid>` and adjust the next request.

## Resources by job type

All from `seff` on real jobs. The short version: **what dominates a step decides
whether more cores help**, and that is not simply "xTB versus DFT".

| Stage | Dominated by | CPUs | CPU efficiency | Memory used | Request |
| --- | --- | --- | --- | --- | --- |
| `preprocess` | file handling | 1-2 | - | 0.7 GB | 4 GB |
| `md` (classical, CPU) | OpenMM | 8 | - | ~1 GB | 8 GB |
| `qmmm` optimisation, xTB | geomeTRIC linear algebra | 8 | **~94%** | 0.9 GB | 4-8 GB |
| `scan`, xTB | same, once per point | 8 | ~94% | ~1 GB | 8-16 GB |
| `qmmm` MD, xTB | xTB + OpenMM only | 8 | **~30%** | 0.5 GB | 8 GB |
| `qmmm` single point or optimisation, **DFT** | PySCF | 16 | **~96%** | 56 GB | 96 GB+ |

**DFT scales, so give it cores.** PySCF saturated 16 and will use more; 64 is
worth requesting for single points and DFT gradients. It is also the only stage
that needs real memory - tens of GB, because density-fitting integrals are held
in RAM and spill to disk if `sp_memory` / `dft_memory` is too small, which is far
slower. mcube defaults those to 80% of `SLURM_MEM_PER_NODE`, so ask for what you
want it to use.

**xTB does not scale, and the reason differs by stage.** An xTB *optimisation* or
*scan* still reaches ~94%, but that is geomeTRIC's linear algebra over hundreds
of active atoms, not xTB. An xTB *MD* run has no optimiser, so only xTB and
OpenMM are left and the job idles at ~30% of 8 cores. Either way xTB itself
saturates around four threads: making xTB work finish sooner means **running
points in parallel** (an array job over scan points or frames), not asking for
more threads per point.

**Memory requests have been 10-100x too large** for everything except DFT. An MD
job asking 32 GB used 0.5 GB. Smaller requests fit on more nodes and start
sooner.

**Threading, which is not the same for the two cases:**

- **xTB stages** - bind threads to physical cores. On SMT nodes the allocated
  CPUs are half physical cores, so `OMP_NUM_THREADS=$(( SLURM_CPUS_PER_TASK / 2 ))`
  with `OMP_PROC_BIND=close` and `OMP_PLACES=cores`. Worth 2-9x; see below.
  Binding is safe because xTB runs as a separate process.
- **DFT stages** - set `OMP_NUM_THREADS=1` and do **not** bind. PySCF runs
  in-process and sets its own thread count from `dft_cores` / `sp_cores`; a
  conflicting environment value makes it warn and can clamp it, and process-wide
  binding pins its BLAS threads to a single core. ASH's OpenMM object sets
  `OMP_NUM_THREADS` to the MM core count as a side effect, so set it explicitly
  after activation rather than relying on the default.

`docs/csr/run.slurm` implements both as the `xtb` and `dft` profiles.

## CPU threads

| Component | What controls the threads | Notes |
| --- | --- | --- |
| OpenMM (MD, and the MM part of QM/MM) | `numcores` passed to `OpenMMTheory`, which sets the platform `Threads` property | ASH defaults to **1**. Before commit 380d193, mcube's `qmmm` and `neb` didn't pass `cores`, so every MM force call during QM/MM ran on one thread. |
| xTB (QM part) | `OMP_NUM_THREADS`, which `xTBTheory(numcores=...)` sets when it is created | Small QM regions (50–100 atoms) are not expected to scale past about 8 threads; this still needs measuring (see below). |
| ASH environment side effects | Every `OpenMMTheory` and `xTBTheory` you create overwrites `OMP_NUM_THREADS` for the whole process | The theory created last decides the thread count for everything that runs afterwards. mcube now passes the same `cores` to both. |
| Preprocess (pdbfixer, OpenFF) | Mostly serial | Measured: 68% of 8 cores over 4.8 min, peak 0.7 GB. 2–4 CPUs would do. It's only a few minutes, so bundling it into a larger job costs little. |
| NEB | `--cores N` (defaults to `$SLURM_CPUS_PER_TASK`) | Split as `min(images, N)` image workers × `N // workers` threads per image, so the CPUs aren't oversubscribed. |

Keep `cores:` in `config.yaml` equal to `--cpus-per-task`. The `md` and `qmmm` stages take their thread count from the config, not from SLURM.

### Bind OpenMP threads, one per physical core (measured 2026-10-16)

- **What SLURM allocates:** on the `general` partition's EPYC nodes (2 hardware threads per core), `--cpus-per-task=8` gives **4 physical cores plus their SMT siblings** (e.g. CPUs 31, 34–36, 287, 290–292).
- **Unbound threads float:** OpenMP threads that aren't bound move between those CPUs and stack onto the same one. On a busy shared node, xTB then becomes slow and erratic. The QM/MM MD of job 29937049 saw 3.6–131 s per step.

Benchmark (job 29937358): the same QM/MM xTB call (GFN2, 124 QM atoms incl. link atoms, 47,893 point charges, 17 SCF iterations), 5 runs per setup:

| Threads | `OMP_PROC_BIND=close`, `OMP_PLACES=cores` | Mean wall time | SD |
| --- | --- | --- | --- |
| **4** | **yes** | **0.80 s** | 0.00 |
| 8 | yes | 1.42 s | 0.01 |
| 4 | no | 1.9–4.7 s | 1.2–1.8 |
| 8 | no | 4.0–5.4 s | 0.5–2.0 |

So, for xTB/QM-MM jobs:

```bash
#SBATCH --cpus-per-task=8
export OMP_NUM_THREADS=4          # one thread per physical core = cpus-per-task / 2
export OMP_PROC_BIND=close
export OMP_PLACES=cores
```

and set `qmmm.inputs.cores: 4`.

**mcube now binds xTB only, by default:**
- It runs xTB through a wrapper (`.mcube_xtb/xtb` in the run directory) that exports `OMP_PROC_BIND=close` and `OMP_PLACES=cores` for xTB alone. When `cores` is unset it uses the physical core count.
- **Don't export `OMP_PROC_BIND`/`OMP_PLACES` in job scripts.** If the Python process is bound, an in-process OpenMP runtime pins the main thread to one core, and threads created later inherit that one-core mask.
  - Measured: a PySCF ωB97M-V single point (job 29992285) ran with its process mask at one core plus its sibling, ~18% CPU per thread. That was mcube's first default, which set the variables process-wide; now fixed.
  - Exporting them yourself, or `MCUBE_THREAD_BINDING=0`, turns mcube's binding off.
- **Measured with binding:** the CSR QM/MM MD (job 29937641, 8 CPUs, `cores: 4`) runs at about **1.4 s per step**: xTB 1.1–1.2 s, OpenMM 0.08 s. That is roughly 20 h for 50 ps; unbound, it was 14.6 s median and up to 131 s.
- **Why fewer threads:** 8 threads land on SMT siblings and are slower than 4.
- **Filesystem:** made no difference. `$TMPDIR` is `/scratch/temp/<jobid>` (not node-local), and once threads are bound, running there or in the job directory gives the same time.

## Memory

| Job | Requested | Measured peak | Suggested request |
| --- | --- | --- | --- |
| `mcube preprocess` (CSR, ~41k atoms after solvation) | 64 GB | 0.7 GB | 4 GB |
| Other CPU jobs on this account (`tri`, `s2b`, `step2`, `infill`) | 24–64 GB | 3–4 GB | ~2× the measured peak |
| `mcube qmmm` optimisation (xTB + OpenMM, ~41k atoms, 8 Å active region) | 16 GB | 0.91 GB (`sstat` MaxRSS), 0.94 GB in `top` | 4 GB |
| `mcube qmmm` MD (xTB + OpenMM, 25 ps) | 32 GB | 0.5 GB (`seff`) | 8 GB |
| `mcube qmmm` **DFT** single point (ωB97X/def2-SVP, 119 QM atoms) | 96 GB | 56 GB (`seff`) | 96 GB+ |

Memory requests have been 10–100× larger than needed **for everything except
DFT**, which is the one stage that genuinely wants tens of GB: density-fitting
integrals are held in RAM and spill to disk when `sp_memory` / `dft_memory` is
too small, which is much slower. Everywhere else, smaller requests fit on more
nodes and start sooner. See "Resources by job type" above.

## Measured: QM/MM geometry optimisation (job 29837804)

The run: CSR reactant, 41,377 atoms, 71 QM atoms (GFN2-xTB, electrostatic embedding), 431 active atoms (8 Å around the Zn), geomeTRIC (ASH default TRIC coordinates), 8 CPUs, `cores: 8`, on node bun133. Measured during steps 3–5 with `srun --overlap top` (8 samples, 2 s apart) and `sstat`.

- **CPU:**
  - The Python (`mcube`) process ran at 720–780% in 7 of 8 samples.
  - A separate `xtb` process (ASH's default `inputfile` runmode) ran at about 540% in 1 sample.
  - So the job keeps roughly 7.5 of 8 cores busy, and the OpenMM thread fix (commit 380d193) works.
  - Most of each step is spent in the Python process (OpenMM and geomeTRIC's linear algebra), not in xTB.
- **Memory:** 0.91 GB peak, about 6% of the 16 GB requested. 4 GB is plenty.
- **Speed:**
  - Roughly 60–70 s per optimisation step.
  - The job reached step 5 after 8 min 16 s, including about 1.1 min of preprocessing and the setup.
  - The 200-step cap allows up to about 4 h per endpoint, about 8 h for reactant plus product, inside the 12 h limit.
  - If the optimisation hits the cap without converging, the endpoint is not a true minimum.

Things to try next time (not yet tested):

- **MM-minimise the solvated box before QM/MM.** The first geomeTRIC step lowered the energy by about 7 Eh, almost all of it MM energy (−84.2 → −91.4 Eh), and the step-0 maximum gradient was 47.8 au. The freshly placed waters are strained, and relaxing them costs expensive QM/MM steps. A quick OpenMM minimisation with the QM atoms frozen would remove that.
- **Try HDLC instead of TRIC coordinates for large active regions.** TRIC on hundreds of active atoms means large linear algebra every step, which fits the Python process dominating. The `Optimizer` in ASH takes a `coordsystem` argument. Benchmark both on the same system before switching.
- **Try `xtb_runmode: library`.** It avoids writing and reading files every step. It probably saves little here, because xTB isn't the bottleneck.
- **Set `PYTHONUNBUFFERED=1` in job scripts.** Python buffers output when writing to a file, so the `.out` showed nothing from the QM/MM stage while it was running. Follow progress in `geometric_OPTtraj.log` and `optimization_energies.log` in the run directory, or set this variable.

## Queueing

- When a job is pending with reason `QOSGrpCpuLimit`, the `normal` QOS as a whole is full. The cap is 25,150 CPUs shared by everyone on the `general` partition (and 1,536 CPUs per user). It isn't a lab or personal limit, and the wait was usually minutes.
- Smaller `--cpus-per-task` and a shorter `--time` help a job start sooner through backfill.
- Check what's going on with:
  - `squeue -u $USER -o '%i %T %r %C %m'`
  - `sacctmgr -nP show qos normal format=Name,GrpTRES,MaxTRESPU`
- `sacct` can report `COMPLETED` even when the science failed, so always read the `.err` and `.out` files as well.

## Checking a running job

```bash
srun --jobid=<id> --overlap --mem=10M -n1 top -b -n 8 -d 2 -u $USER   # threads busy (%CPU ≈ 100 × threads)
sstat -j <id>.batch -o MaxRSS,AveCPU               # memory so far
seff <id>                                          # after it finishes
```

- **Give the `srun` step a small `--mem`.** Without it, the step asks for memory that the job has already allocated and fails with "Memory required by task is not available".
- **Take several samples.** In ASH's default mode xTB runs as a separate, short-lived `xtb` process, so one snapshot can miss it.
- **Look for single-threaded work.** If `top` shows a Python process at about 100% while `cores` > 1, something is running on one thread.

## GPUs

Checked 2026-10-16 with `sinfo` and `sacctmgr`:
- **NVIDIA:** `gpu_cuda` (A100, H100, L40S; also MIG slices) and `gpu_sxm` (H100 SXM).
- **AMD:** `gpu_rocm` (MI210, MI300X, MI355X).
- **Access:** our account has QOS `gpu` (at most 4 GPUs and 256 CPUs per user), `sxm` (4 GPUs), `mig` and `debug` (1 h).

The ASH env has OpenMM 8.4 built against CUDA 13.1. Its platforms are Reference, CPU and **CUDA**. The HIP plugin fails to load (no ROCm libraries), so use the NVIDIA partitions:

```bash
#SBATCH --partition=gpu_cuda
#SBATCH --qos=gpu
#SBATCH --gres=gpu:1          # or gpu:h100:1 / gpu:a100:1
```

- **Classical stages on the GPU:** `md.inputs.platform: CUDA` (with `properties: {Precision: mixed}`) runs the md stage there.
- **QM/MM stays on the CPU:** xTB has no GPU code, so keep QM/MM jobs on `general`.
- **Checking the context:** OpenMM failing to open a CUDA context means the node's driver is older than the installed CUDA needs; test it before queueing a long run.

## Future option: FMO-xTB for large QM regions

FMO-xTB, the fragment molecular orbital method with GFN1-xTB, makes QM regions of thousands of atoms affordable (Einsele & Mitrić, [arXiv:2606.28022](https://arxiv.org/abs/2606.28022), 2026). It supports FMO2/FMO3, analytic gradients, covalent cuts via hybrid orbital projection, and elements up to Rn.

- **Where it lives:** only in the group's Rust program DIALECT ([dialect-rs/DIALECT-rs](https://github.com/dialect-rs/DIALECT-rs), LGPL-3.0). It isn't in ASH or xtb/tblite as of October 2026.
- **How it runs:** command line, `dialect <file.xyz>` with a `dialect.toml` (`method = "xtb"`, `[fmo] enabled = true`). A `jobtype = "grad"` run writes the energy to `dialect.out` and the gradient to `gs_gradient.npy`.
- **No electrostatic embedding:** DIALECT has no external point charges, so an ASH wrapper (`DIALECTTheory`) could do gas-phase/cluster calculations or mechanical-embedding QM/MM only.
- **Other limits:** fragments are closed-shell, and only GFN1 is available (mcube uses GFN2).
- **Not worth it for CSR:** the QM region is ~70 atoms, which plain xTB handles easily. Revisit if a project needs a QM region of hundreds to thousands of atoms. Electrostatic embedding would first need point-charge support added to DIALECT.

## Preparing systems: problems we hit

- **A ligand already in `protein_pdb`.** If the complex PDB already contains the ligand, set `merge_ligand: false`. Otherwise preprocess adds a second copy (seen in the CSR runs as `mer:102`).
- **Interleaved residues.** If one residue's atoms sit inside another residue's block (here, hydroxide atoms placed inside meropenem to match the product's atom order), OpenMM, pdbfixer and RDKit split the residue into fragments. The error looks like "no template found… matches NME". Use `scripts/make_residues_contiguous.py reactant.pdb product.pdb`, which applies the same reordering to both files.
- **Bond orders for PDB ligands.**
  - RDKit's distance-based bonding added a spurious O–O bond in a strained geometry, and its automatic bond-order assignment produced enolate/S⁺ forms for carbapenems.
  - mcube now reads bonds from CONECT records only. Give the chemistry with `ligand_smiles`, and check the `[INFO] Ligand SMILES:` line in the output.
- **`charge` is the charge of the ligand file only.** In CSR_ligand, meropenem is neutral (0) and the separate OH⁻ residue gets its −1 from `hydroxide.xml`. In CSR_hydrolysed, the 53-atom ligand includes that OH, so it is −1.
- **QM charge and multiplicity.** Add up the formal charges of the QM residues. Zn²⁺ is d¹⁰, so multiplicity 1. For CSR: Zn²⁺ + OH⁻ + meropenem + Asp⁻ + His = 0, and the product also comes to 0.
- **Coordination cutoff.** Zn–ligand bonds are about 2.0–2.3 Å. A cutoff of 4 Å picked up Arg90 (3.5 Å) as a "ligand" and constrained it. Use 2.5–3.0 Å, and check the "Coordinating residues" table in the preprocess output.

## NEB endpoints

NEB needs reactant and product with identical atoms in identical order, solvent included.

1. Preprocess (solvate) the **reactant** only.
2. Build the product inside the same box: `scripts/build_product_box.py finalsystem.pdb product_solute.pdb product_solvated.pdb`. It checks atom identity, aligns on the non-reactive solute atoms and reports solvent clashes. For CSR: RMSD 0.000 Å, no water within 2.0 Å.
3. QM/MM-optimise both endpoints with the same active region:
   - `mcube qmmm config.yaml --pdb ... --section qmmm_reactant`
   - then `... --section qmmm_product --actatoms ../reactant_opt/active_atoms`
4. Run NEB on the two `QM_MM_optimized.pdb` files with the same `--actatoms`.

Running MD on the reactant and the product separately gives each its own water configuration, so those trajectories can't be used as NEB endpoints.

## Driving Bunya through code-server

See `binder-physics/docs/BUNYA_BROWSER_SHELL.md`: keep typed commands under about 100 characters, write output to files under `/scratch/user/uquadhik/...` and fetch them, never use heredocs, and use `printf` or `sed` to build scripts.
