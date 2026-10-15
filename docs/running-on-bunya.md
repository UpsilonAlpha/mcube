# Running mcube efficiently on Bunya

Notes from setting up and running the CSR (meropenem) reactant/product systems in October 2026. The numbers come from `seff`/`sacct` on real jobs unless marked as estimates.

## Quick recipe

```bash
#SBATCH --partition=general --qos=normal --account=a_imb_ccbcd
#SBATCH --cpus-per-task=8          # match `cores:` in config.yaml
#SBATCH --mem=16G                  # see the memory table below
#SBATCH --time=12:00:00            # shorter limits start sooner (backfill)

set -euo pipefail
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
ulimit -s unlimited                # xTB recommendation, avoids stack segfaults
export OMP_STACKSIZE=1G

module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH                 # mcube is pip-installed (editable) in this env

(cd preprocess_results && mcube preprocess ../config.yaml)
```

- Run each stage in its own directory. Stages write many files to the current directory and write results back into `config.yaml`.
- Don't run two stages at the same time against one `config.yaml`. Each stage reads, edits and rewrites the whole file, so concurrent runs can silently overwrite each other's results. Run them one after another in a single job, or give each its own config.
- After a job finishes, run `seff <jobid>` and adjust the next request.

## CPU threads

| Component | What controls the threads | Notes |
| --- | --- | --- |
| OpenMM (MD, and the MM part of QM/MM) | `numcores` passed to `OpenMMTheory`, which sets the platform `Threads` property | ASH defaults to **1**. Before commit 380d193, mcube's `qmmm` and `neb` didn't pass `cores`, so every MM force call during QM/MM ran on one thread. |
| xTB (QM part) | `OMP_NUM_THREADS`, which `xTBTheory(numcores=...)` sets when it is created | Small QM regions (50–100 atoms) are not expected to scale past about 8 threads; this still needs measuring (see below). |
| ASH environment side effects | Every `OpenMMTheory` and `xTBTheory` you create overwrites `OMP_NUM_THREADS` for the whole process | The theory created last decides the thread count for everything that runs afterwards. mcube now passes the same `cores` to both. |
| Preprocess (pdbfixer, OpenFF) | Mostly serial | Measured: 68% of 8 cores over 4.8 min, peak 0.7 GB. 2–4 CPUs would do. It's only a few minutes, so bundling it into a larger job costs little. |
| NEB | `--cores N` (defaults to `$SLURM_CPUS_PER_TASK`) | Split as `min(images, N)` image workers × `N // workers` threads per image, so the CPUs aren't oversubscribed. |

Keep `cores:` in `config.yaml` equal to `--cpus-per-task`. The `md` and `qmmm` stages take their thread count from the config, not from SLURM.

## Memory

| Job | Requested | Measured peak | Suggested request |
| --- | --- | --- | --- |
| `mcube preprocess` (CSR, ~41k atoms after solvation) | 64 GB | 0.7 GB | 4 GB |
| Other CPU jobs on this account (`tri`, `s2b`, `step2`, `infill`) | 24–64 GB | 3–4 GB | ~2× the measured peak |
| `mcube qmmm` optimisation (xTB + OpenMM, ~41k atoms, 8 Å active region) | 16 GB | still to be measured (see below) | — |

Memory requests have been 10–100× larger than needed. Smaller requests fit on more nodes and start sooner.

## Queueing

- When a job is pending with reason `QOSGrpCpuLimit`, the `normal` QOS as a whole is full. The cap is 25,150 CPUs shared by everyone on the `general` partition (and 1,536 CPUs per user). It isn't a lab or personal limit, and the wait was usually minutes.
- Smaller `--cpus-per-task` and a shorter `--time` help a job start sooner through backfill.
- Check what's going on with:
  - `squeue -u $USER -o '%i %T %r %C %m'`
  - `sacctmgr -nP show qos normal format=Name,GrpTRES,MaxTRESPU`
- `sacct` can report `COMPLETED` even when the science failed, so always read the `.err` and `.out` files as well.

## Checking a running job

```bash
srun --jobid=<id> --overlap top -b -n 1 -u $USER   # threads actually busy (%CPU ≈ 100 × threads)
sstat -j <id>.batch -o MaxRSS,AveCPU               # memory so far
seff <id>                                          # after it finishes
```

If `top` shows a Python process at about 100% while `cores` > 1, something is running single-threaded.

## GPUs (not used yet)

Bunya's GPU partition is AMD MI355X (`-p gpu_rocm --gres=gpu:mi355x:1 --qos=gpu`), with at most 4 concurrent jobs per user. OpenMM MD would be far faster there, but it needs OpenMM's HIP platform. mcube always uses `platform='CPU'`, and the HIP plugin hasn't been checked in the ASH env. Adding a `platform:` key to `md.inputs` would be the next step for long production MD. QM/MM runs on the CPU either way because of xTB, so keep those on `general`.

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
