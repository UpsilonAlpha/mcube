# mcube

A pipeline for protein–ligand systems, often metalloproteins, built on [ASH](https://github.com/RagnarB83/ash). It runs:

1. **preprocess**: ligand parameterisation (OpenFF), protonation and solvation (explicit or implicit), metal-coordination bond constraints and QM-region selection
2. **md**: OpenMM minimisation, optional gentle warm-up, NPT box equilibration and production MD
3. **qmmm**: xTB/MM MD and geometry optimisation, with an optional PySCF QM/MM single point and cube files
4. **scan**: relaxed QM/MM distance scan (builds intermediates/products in the reactant's box)
5. **neb**: QM/MM nudged elastic band between a reactant and a product (optionally NEB + TS optimisation), with a frequency check on the saddle point
6. **ts**: TS optimisation with a partial Hessian, frequency check and IRC (quasi-IRC along the TS mode, or Sella)
7. **ramd**: τRAMD relative residence times (random acceleration MD on the md stage's system, bootstrapped τ per replica)
8. **postprocess**: MD data plots, RMSD, RMSF (as B-factors) and a residue interaction network

All stages read one YAML config file and write their results back into it, so each stage picks up the outputs of the stages before it.

## Install

```bash
conda env create -f environment.yml   # conda-forge deps + ASH from git + mcube (editable)
conda activate mcube
```

`environment.lock.yml` is the exact linux-64 environment the pipeline was developed in.

If ASH and the conda packages are already installed, install mcube on its own with:

```bash
pip install -e .
```

## Usage

```bash
mcube preprocess config.yaml
mcube md config.yaml
mcube qmmm config.yaml
mcube scan config.yaml --pdb reactant.pdb --section scan_attack [--actatoms active_atoms]
mcube neb config.yaml --reactant reactant.pdb --product product.pdb [--section neb_step1] [--actatoms active_atoms] [--cores N]
mcube ts config.yaml --pdb saddlepoint.pdb [--section ts_step1] [--actatoms active_atoms]
mcube ramd config.yaml [--section ramd_complexA] [--replicas 0 1] [--analyse]
mcube postprocess config.yaml [--qmmm] [--redraw] [--cutoff 2.6] [--metals ZN MG] [--occupancy 0.1]
```

`python -m mcube ...` also works. Every stage writes its files to the current directory, so run each stage from its own working directory. Paths saved in `results` are absolute, so later stages can be run from anywhere.

`--cores` defaults to `$SLURM_CPUS_PER_TASK` (or 1). The `md`, `qmmm`, `scan`, `ts` and `ramd` stages take their core count from `cores` in the config instead.

**Threads.** xTB runs with its OpenMP threads bound one per physical core (`OMP_PROC_BIND=close`, `OMP_PLACES=cores`), through a small wrapper script that mcube writes to `.mcube_xtb/` in the run directory. The Python process itself (OpenMM, PySCF, numpy) is not bound: binding an in-process OpenMP runtime pins later threads to one core. Setting `OMP_PROC_BIND`/`OMP_PLACES` yourself, or `MCUBE_THREAD_BINDING=0`, turns mcube's binding off.
- **Why:** on SMT nodes, `--cpus-per-task=N` gives N/2 physical cores. Unbound xTB threads then share cores and ran 2–9× slower, with erratic step times (see `docs/running-on-bunya.md`).
- **Default `cores`:** when `cores` is not set, a stage uses the number of physical cores it is allowed to run on (e.g. 4 for `--cpus-per-task=8`). A larger `cores` gives a warning.
- **Parallel workers:** NEB with parallel images and `ts` with `freq_workers` > 1 turn binding off, because each worker would bind to the same cores.

## Config file

The config file has one section per stage, and each section has `inputs` (written by you) and `results` (written by mcube). When an input is missing, the stage falls back to the earlier stage's results: for example, `md` uses `preprocess.results.processed_pdb` unless `md.inputs.pdb_file` is set. Start from a template in [`config_examples/`](config_examples):

- [`pdbqt_config.yaml`](config_examples/pdbqt_config.yaml): explicit solvent, docked ligand from a PDBQT
- [`implicit_config.yaml`](config_examples/implicit_config.yaml): implicit solvent (GBn2), ligand from a PDB
- [`protein_protein_ramd.yaml`](config_examples/protein_protein_ramd.yaml): protein–protein complex without metals, md and τRAMD on a GPU (`scripts/ramd_gpu.slurm`)

Results are written back with `yaml.safe_dump`, so comments in the config are lost the first time a stage saves its results. Writes are locked (POSIX), so several jobs can run stages on one config at once as long as each writes its own `--section`.

**Atom selectors.** Wherever an input takes atoms (`md` `constraints`/`restraints`; `actregion_origin`, `opt_constraints`; scan `atoms`/`constraints`; ts `hessian_atoms`/`monitor`; neb `hessian_atoms`/`partial_hessian_atoms`), each atom can be a 0-based index or a `"chain:resid:atom"` string matched against the stage's structure, e.g. `"B:1:ZN"` or `"A:88:NE2"`. Numbers in the same entry (distances, force constants) are left alone, so `restraints: [["B:1:ZN", "C:1:O2x", 2.2, 20.0]]` works. Selectors keep a config valid when the solvent box (and so every index after the solute) changes, and make configs readable for new systems. The resolved indices are printed.

### `preprocess.inputs`

| Key | Default | Description |
| --- | --- | --- |
| `protein_pdb` | required | Protein PDB file |
| `ligand` | none | `.pdbqt` (Meeko output with `REMARK SMILES`/`SMILES IDX`), `.pdb` (explicit hydrogens and CONECT records required) or `.mol` |
| `ligand_smiles` | none | `.pdb` ligands only: SMILES (with charges and stereo) used as a bond-order template. Recommended. Without it RDKit guesses bond orders and can pick unphysical resonance forms (e.g. enolate/S⁺ for carbapenems). Connectivity always comes from CONECT, never from distances. The SMILES actually used is printed. |
| `merge_ligand` | `true` | `false` if `protein_pdb` already contains the ligand: it is parameterised but not merged in again |
| `ligand_resname` | ligand file stem | Ligand residue name. Use at most 3 characters. |
| `model` | | PDBQT model (pose) number, starting at 1 |
| `charge` | | Ligand formal charge |
| `ligand_xmlfile` | | Existing ligand force-field XML, used when `ligand` is not given |
| `default_forcefields` | `amber14/protein.ff14SB.xml`, `amber14/tip3p.xml` | OpenMM force-field XML files. The ligand XML is appended automatically. |
| `residue_variants` | none | Protonation states, e.g. `{A: {80: HIE}}` |
| `skip_solvate` | `false` | Use the PDB as it is (protein, or `merged.pdb` with a ligand) |
| `implicit` | `false` | Add hydrogens only (no water box), for implicit-solvent MD |
| `pH` | 8.0 explicit / 7.0 implicit | Protonation pH |
| `solvent_padding` | 10.0 | Water padding in Å (explicit only) |
| `ionic_strength` | 0.1 | NaCl concentration in M (explicit only) |
| `metals` | none | Element symbols treated as metals. Empty or absent (e.g. protein–protein complexes): no coordination constraints and no QM region. |
| `cutoff` | required with `metals` | Metal–donor distance cutoff (Å) |
| `donor_elements` | `[N, O, S]` | Elements that count as metal donors (matched on element, not atom name) |
| `require_his_variants` | `false` | Fail if any histidine has no tautomer in `residue_variants`. A warning is printed either way. |
| `extra_qm_residues` | none | `{chain: [resid, ...]}`: residues added to the QM region on top of the metal's coordination sphere (e.g. catalytic Lys/Arg) |
| `constraint_exclude_residues` | none | `{chain: [resid, ...]}`: residues whose metal–donor distance is left unconstrained during MD (e.g. the substrate) |
| `exclude_residues` | none | `{chain: [resid, ...]}`: coordinating residues left out of the QM region. Their metal–donor bond constraints are still applied in classical MD. Not the same as `boundary_excluded`, which lists the non-amino-acid QM atoms that ASH should not place link atoms on. |
| `backbone_atoms` | standard backbone | Atoms of amino-acid QM residues that stay in MM |
| `amino_acids` | 20 standard + Amber variants | Residues cut at CA–CB. Any other residue (metals, ligands, waters) is put fully into the QM region and excluded from boundary detection. |

`results`: `processed_pdb` (`finalsystem.pdb`), `ligand_xml`, `ligand_resname`, `qm_residues` (`{chain: {resid: resname}}`), `qm_atoms` and `boundary_excluded` (atom indices into `processed_pdb`), and `constraints` (metal–donor atom-index pairs).

### `md.inputs`

| Key | Default | Description |
| --- | --- | --- |
| `pdb_file`, `constraints`, `ligand_xml`, `default_forcefields` | from preprocess | Overrides |
| `cores` | | OpenMM CPU threads |
| `platform`, `properties` | `CPU` | OpenMM platform, e.g. `CUDA` with `{Precision: mixed}` on a GPU node (see `docs/running-on-bunya.md`) |
| `temp` | | Temperature (K) |
| `implicit` | `false` | GBn2 implicit solvent. Uses HBond constraints, 1.5 amu hydrogen mass and a 2 nm cutoff. NPT and reimaging are skipped. |
| `implicit_eps` | | Solvent dielectric for implicit solvent |
| `gentle` | `false` | Gentle warm-up (1 → 10 → 300 K) after minimisation |
| `skip_npt` | `false` | Skip NPT box equilibration (timestep `npt_timestep`) |
| `skip_nvt` | `false` | Skip production MD |
| `nvt_timestep`, `nvt_time` | | Production timestep and length (ps) |
| `traj_frequency` | 1000 | Production trajectory save interval (steps) |
| `constraints` | from preprocess | Fixed distances as `[[i, j], …]`, taken from the starting geometry. Override to drop or add pairs. (ASH silently ignores `[i, j, d]` entries with an explicit distance.) |
| `restraints` | none | Harmonic distance restraints `[[i, j, d_Å, k_kcal/mol/Å²], …]`, active during minimisation, warm-up, NPT and NVT. Use them to pull a substrate into place, then run again without them. |

Minimisation always runs. Metal–donor `constraints` are applied as fixed-distance OpenMM constraints.

`results`: `minimized`, `gentle`, `forcefields`, `trajfile`, `lastframe`, `csv`. In explicit solvent, `trajfile` and `lastframe` are the reimaged files.

### `qmmm.inputs`

| Key | Default | Description |
| --- | --- | --- |
| `pdb_file` | `md.results.lastframe` | Starting structure. It must have the same atom order as `processed_pdb`. |
| `forcefields` | `md.results.forcefields` | |
| `qm_residues` | from preprocess | Override the QM region (rebuilt on `pdb_file`) |
| `implicit` | `md.inputs.implicit` | Non-periodic MM when true |
| `charge`, `mult` | required | QM-region charge and multiplicity |
| `cores` | | xTB / PySCF cores |
| `xtb_method` | `GFN2` | |
| `xtb_runmode` | `inputfile` | `library` avoids file I/O at every step (needs xtb-python) |
| `temp`, `timestep`, `time` | | QM/MM MD settings (K, ps, ps) |
| `traj_frequency` | 1000 | QM/MM MD save interval (steps) |
| `skip_equilibrate` | `false` | Skip QM/MM MD |
| `skip_optimize` | `false` | Skip QM/MM optimisation |
| `opt_constraints` | none | Bond constraints during the optimisation: `[[i, j], …]` holds the starting distance, `[[i, j, d_Å], …]` sets it. Atoms must be in the active region. Active-region waters are always kept rigid. |
| `opt_maxiter` | 200 | Optimisation steps. If they run out, the last geometry is kept as `QM_MM_optimized.pdb` and `converged: false` is recorded. |
| `actradius` | 0 | If > 0, optimise all atoms within this radius (Å) of `actregion_origin`. Otherwise only the QM atoms move. |
| `actregion_origin` | first non-standard QM atom | Atom index at the centre of the active region |
| `singlepoint` | `true` | PySCF QM/MM single point (RKS/UKS from `mult`), plus `density.cube`, `mep.cube`, `HOMO.cube` and `LUMO.cube` |
| `sp_functional`, `sp_basis` | `wb97x-v`, `def2-svp` | Single-point level of theory (PySCF names, e.g. `wb97m_v`; the energy is saved as `sp_energy`) |
| `sp_gridlevel` | 3 | PySCF DFT grid level (ASH's own default is 5, which made VV10 functionals very slow) |
| `sp_densityfit`, `sp_auxbasis` | `true`, PySCF default | Density fitting (RI-JK for hybrids) |
| `sp_cores` | physical cores | Threads for the single point (independent of the xTB `cores`) |
| `sp_platform` | `CPU` | `GPU` uses gpu4pyscf (must be installed) |
| `sp_memory` | 80% of `$SLURM_MEM_PER_NODE`, else 3100 | PySCF `max_memory` (MB). Density-fitting integrals that don't fit are written to disk and re-read every cycle; a range-separated hybrid on ~120 QM atoms with def2-SVP needs ~60 GB to stay in RAM. |
| `sp_dispersion`, `sp_dispersion_method` | none, `sp_functional` | Post-SCF DFT-D3 on the QM region (link atoms included): `d3bj` or `d3zero`, parameters for the named functional (simple-dftd3 names, e.g. `wb97x`). Needs `dftd3-python`. Saved as `sp_disp_energy` and `sp_energy_disp`. |

`results` (section `qmmm`, or `--section NAME`): `qm_atoms`, `boundary_excluded`, `trajfile`, `lastframe`, `opttraj`, `optimized_pdb` (`QM_MM_optimized.pdb`, the full optimised system, usable as an NEB endpoint), `opt_energy` (Eh), `converged`, `opt_constraints` (with the distances used), `cubes`.

`mcube qmmm config.yaml --pdb other.pdb --section qmmm_product` runs the same inputs (overridden by `qmmm_product.inputs`, if given, e.g. `skip_equilibrate: true` or `opt_constraints`) on another structure (e.g. each NEB endpoint) without overwriting the other run's results. `--actatoms FILE` (qmmm, scan, neb and ts) reuses an active-atom list such as the reactant run's `active_atoms`, so both NEB endpoints and the NEB share one active region.

`scan`, `neb` and `ts` use `qmmm.inputs` for `charge`, `mult`, `xtb_method`, `forcefields`, `cores` (scan, ts), the QM region and (unless overridden) `actradius`/`actregion_origin`. Each reads its own base section's `inputs` (`scan`, `neb`, `ts`), overridden by the inputs of the `--section` it writes its results to, so several scans or steps can live in one config:

```yaml
scan:
  inputs: {points: 16, opt_maxiter: 200}
scan_attack:            # mcube scan config.yaml --pdb R.pdb --section scan_attack --actatoms active_atoms
  inputs: {atoms: [4255, 4208], start: 3.0, end: 1.45}
```

### `scan.inputs`

| Key | Default | Description |
| --- | --- | --- |
| `atoms` | required | `[i, j]`: the distance scanned (0-based indices) |
| `start` | current distance | First value (Å) |
| `end` | required | Last value (Å) |
| `points` | 10 | Number of points, both ends included |
| `constraints` | none | Extra bonds held during the scan, `[[i, j]]` or `[[i, j, d]]` |
| `opt_maxiter` | 200 | Optimisation steps per point. A point that runs out keeps its last geometry, is marked `converged: false`, and the scan continues. |
| `monitor` | none | Distances `[[i, j], …]` reported at every point (CSV and results) |
| `pdb_file` | `qmmm.results.optimized_pdb` | Starting structure (or `--pdb`) |

Each point is optimised from the previous one with the scanned distance fixed. Writes `scan_XX.pdb` (full system) and `scan_energies.csv`; `results` (`points` with target/actual distance, energy, kcal/mol relative to point 0, PDB; `max_point`) are saved after every point, so a stopped scan can be continued from its last PDB with a new `start`.

### `neb.inputs` (optional section)

| Key | Default | Description |
| --- | --- | --- |
| `images` | 12 | Number of NEB images |
| `method` | `neb` | `nebts`: ASH NEBTS, a loose NEB-CI then a geomeTRIC TS optimisation of the climbing image |
| `CI`, `maxiter`, `conv_type`, `tol_scale`, `tol_max_fci`, `tol_rms_fci`, `tol_turn_on_ci`, `cichange_freqcheck`, `free_end`, `interpolation`, `idpp_maxiter`, `idpp_springconst`, `restart_file`, `TS_guess_file` | ASH | Passed to ASH's `NEB`/`NEBTS` (ASH defaults: CI on, 100 iterations, tol_max_fci 0.026 / tol_rms_fci 0.013 Eh/Å for NEB, 0.10 / 0.05 for NEBTS, IDPP interpolation) |
| `springconst`, `springconst2`, `energy_weighted`, `zoom` | ASH | NEB only |
| `hessian_for_TS`, `partial_hessian_atoms`, `partial_npoint_choice`, `modelhessian`, `tsmode_tangent_threshold`, `OptTS_maxiter`, `OptTS_convergence_setting`, `OptTS_conv_criteria`, `OptTS_coordsystem` | `hessian_for_TS: partial`, else ASH | NEBTS only |
| `actradius`, `actregion_origin` | qmmm values | Active region. With no radius, only the QM atoms move. |
| `numfreq` | `true` | Numerical frequencies on the saddle point |
| `hessian_atoms` | QM atoms | Atoms in the frequency calculation (cost ≈ 6 gradients per atom) |

`--cores` is split between images that run in parallel and xTB within each image.

`results` (section `neb` or `--section`): `qm_atoms`, `method`, `images`, `options`, `mep_energies`, `saddlepoint_xyz`, `saddlepoint_pdb` (full system, the input for `mcube ts`), `saddlepoint_energy`, `frequencies`, `imaginary_frequencies`.

### `ts.inputs` (optional section)

| Key | Default | Description |
| --- | --- | --- |
| `pdb_file` | `neb.results.saddlepoint_pdb` | TS guess (or `--pdb`) |
| `hessian_atoms` | QM atoms | Atoms with an exact (numerical) Hessian, for the TS optimisation and the frequency check. The reacting atoms and their neighbours are enough; cost ≈ 3 (npoint 1) or 6 (npoint 2) gradients per atom. |
| `monitor` | none | `[[i, j], …]` distances reported for the guess, the TS and both IRC ends |
| `tsopt` | `true` | geomeTRIC TS optimisation; the rest of the active region's Hessian comes from `modelhessian` (`unit`, `Almloef`, `Lindh`, `Schlegel`, `xtb`) |
| `tsopt_maxiter`, `convergence_setting` | 100, geomeTRIC default | |
| `hessian_npoint`, `freq_npoint` | 1, 2 | One- or two-sided displacements for the TS-optimisation Hessian and the frequency check |
| `freq` | `true` | Frequencies at the TS. A TS has exactly one imaginary frequency; the atoms moving most in that mode are reported. |
| `freq_workers` | 1 | Displacements run in parallel (xTB `library` mode, `cores` split between them) |
| `imag_threshold` | 50 | Imaginary frequencies smaller than this (cm⁻¹) are reported separately as numerical noise |
| `irc` | `true` | Follow the TS downhill both ways |
| `irc_method` | `mode` | `mode`: displace ±`irc_step` Å (largest atom displacement) along the imaginary mode and minimise each side. `sella`: Sella IRC on the active region (`irc_dx`; needs `sella`, expensive for a large active region). |
| `opt_maxiter` | 300 | Steps per IRC minimisation |
| `rigid_water` | `true` | Keep active-region waters rigid, as in the `qmmm` minimisations |

`results`: `ts_pdb` (`ts_optimized.pdb`), `ts_energy`, `guess_imaginary`, `imaginary_frequencies`, `small_imaginary`, `lowest_frequencies`, `ts_mode_atoms`, monitored distances, and `irc.forward/reverse` (`irc_*.pdb`, energy, barrier = E(TS) − E(end) in kcal/mol, distances).

### `ramd.inputs` (τRAMD)

Random acceleration MD (Lüdemann et al. 2000) with the τRAMD protocol (Kokh et al., *JCTC* 2018; GROMACS implementation, Kokh et al. 2020). A constant force acts on the binder's centre of mass along a random direction; every `eval_interval` the binder COM displacement *relative to the receptor COM* is checked, and if it is below `rmin` a new direction is drawn. A trajectory ends when the binder–receptor COM distance reaches `max_distance`. Each replica starts with its own unbiased equilibration (new velocities, optionally its own snapshot), then runs `trajectories` RAMD runs with new velocities. τ per replica is the time at which half of its trajectories have dissociated (mean of the medians of random 80% subsets, as `tauRAMD.py`); τ is the mean over replicas.

The system is the **md stage's** (same force fields, HBond constraints, rigid water, 1.5 amu hydrogens and metal-site constraints), so a metal site stays held as in every MM stage. md restraints are not applied. A constraint or restraint that ties the binder to anything else is refused, because it would prevent unbinding.

| Key | Default | Description |
| --- | --- | --- |
| `binder` | required | `{chain: [resid, …]}`, `{chain: all}` (e.g. a binder protein) or a list of atoms (indices or `chain:resid:atom`) |
| `receptor` | everything else except water and monatomic ions | Same format |
| `force` | 14 | Force magnitude, kcal/mol/Å (= 585.76 kJ/mol/nm). **Keep it the same for every complex you compare.** |
| `eval_interval` | 0.1 | ps between direction checks |
| `rmin` | 0.025 | Å of binder COM motion (relative to the receptor COM) per check below which the direction is redrawn |
| `max_distance` | receptor radius + binder radius + `exit_margin` | Binder–receptor COM distance (Å) that counts as dissociated. Must be below half the box. |
| `exit_margin` | 10 | Å added to the summed radii for the automatic `max_distance` |
| `replicas`, `trajectories` | 4, 15 | Protocol size (≥ 4 replicas × ≥ 10–15 trajectories) |
| `equil_time` | 100 | ps of unbiased MD per replica before its RAMD runs |
| `max_time` | 10000 | ps per trajectory; runs that do not dissociate are kept as censored at this time |
| `pdb_file` | `md.results.lastframe` | Start structure |
| `replica_pdbs` | none | One start structure per replica (e.g. from `scripts/extract_frames.py`) |
| `timestep`, `temp` | `md.inputs.nvt_timestep`, `md.inputs.temp` | ps, K |
| `friction` | 1.0 | Langevin friction, 1/ps (the thermostat strongly affects dissociation times; keep it fixed) |
| `pressure` | 1.0 | bar, Monte Carlo barostat (explicit solvent); `null` for NVT |
| `remove_cm_motion` | `true` | Keep the system's CMMotionRemover |
| `constraints` | md constraints | Override the metal-site constraints |
| `restraints` | none | Harmonic restraints, as in md |
| `platform`, `properties` | fastest available | OpenMM platform (`CUDA`, `OpenCL`, `CPU`) and its properties (e.g. `{Precision: mixed}`) |
| `cores` | `md.inputs.cores` | CPU threads (CPU platform) |
| `seed` | 2026 | Base seed for velocities and directions (per replica and trajectory) |
| `traj_interval` | 0 | ps between DCD frames of each RAMD trajectory (0: none) |
| `exit_pdb` | `true` | Write the last frame of each trajectory |
| `bootstrap_rounds`, `bootstrap_fraction` | 5000, 0.8 | τ bootstrap |
| `min_trajectories` | 10 | Warn below this many trajectories per replica |

Box size is the main cost: the binder must not reach the receptor's periodic image before it counts as dissociated, so the box edge should exceed `max_distance` + both radii + the cutoff (a warning says by how much). Set `preprocess.solvent_padding` for this before `md`.

Work is restartable and can be split: each replica uses `replica_RR/` (`equilibrated.xml`, `traj_TT_end.pdb`), stop times go to `ramd_times.csv` (one row per trajectory, file-locked), finished trajectories are skipped, and `--replicas 2` runs one replica (one SLURM array task each). `--analyse` re-analyses the CSV, which can also come from another engine (columns `replica`, `trajectory`, `time_ps`, `dissociated`; set `times_csv`).

`results`: `tau_ps`, `tau_sd_ps` (spread over replicas), `mean_bootstrap_sd_ps`, per-replica `tau_ps`/`sd_ps`/`median_ps`/`ks` (KS distance to the Poisson CDF 1 − exp(−t/τ); 0.1–0.2 is typical), and `warnings` (too few replicas or trajectories, replica scatter above the bootstrap error, τ < 100 ps, censored majority, large KS).

## Outputs by stage

| Stage | Main files |
| --- | --- |
| preprocess | `<resname>.mol/.pdb`, `openff_<resname>.xml`, `merged.pdb`, `finalsystem.pdb` |
| md | `frag-minimized.pdb`, `warmup_MD_cycle*`, `relaxbox_NPT.dcd`, `nptsim.csv`, `NVTtrajectory(.dcd, _imaged.dcd, _lastframe(_imaged).pdb)`, `nvtsim.csv`, `system.xml` (implicit) |
| qmmm | `QM_MM.dcd`, `QM_MM_lastframe.pdb`, `active_atoms`, `geometric_OPTtraj_Full.xyz`, `geometric_OPTtraj_QMregion.xyz`, `*.cube` |
| scan | `scan_XX.pdb`, `scan_energies.csv`, `active_atoms` |
| neb | ASH NEB files, `saddlepoint.xyz`, `saddlepoint.pdb` |
| ramd | `ramd_times.csv`, `replica_RR/equilibrated.xml`, `replica_RR/equilibration.dcd`, `replica_RR/traj_TT_end.pdb` (and `traj_TT.dcd`) |
| ts | `Hessian_partial_TS`, `ts_optimized.pdb`, `irc_forward.pdb`, `irc_reverse.pdb`, `irc_*_opttraj.xyz` (or Sella's `sella_*.xyz`) |
| postprocess | `Temperature.png`, `Density.png`, `BoxVolume.png`, `RMSD.png`, `interaction_network.csv/.png`, `rmsf_bfactor.pdb` |

## postprocess

By default, postprocess analyses the MD trajectory (`md.results`). With `--qmmm`, it analyses the QM/MM trajectory and also converts the optimisation trajectory to `output.dcd`. It computes:

- **Plots** of temperature, density and box volume from the MD CSV (MD only), and an RMSD plot of the protein heavy atoms
- **Interaction network** with H-bonds (protein and ligand), salt bridges (< 4 Å), π–π (< 5 Å), cation–π (< 6 Å) and metal coordination (`--cutoff`). Edges whose occupancy is above `--occupancy` are kept. `--redraw` redraws the plot from `interaction_network.csv` without recomputing. `--spring_constant` and `--iterations` control the layout.
- **RMSF** of Cα atoms after alignment, written as B-factors to `rmsf_bfactor.pdb`

## Running on Bunya

See [docs/running-on-bunya.md](docs/running-on-bunya.md) for thread and memory settings, queue behaviour, system-preparation problems we hit, and the NEB endpoint workflow.

## Standalone scripts

All are in `scripts/`:

- `plot_md_data.py <data.csv>`: plots temperature (plus density and volume for NPT) from an OpenMM CSV
- `plot_rmsd_via_mdtraj.py <traj.dcd> <topology.pdb>`: plots RMSD of the protein heavy atoms
- `make_residues_contiguous.py <reactant.pdb> [product.pdb ...]`: makes each residue's atoms contiguous in the reference PDB (readers split interleaved residues into fragments) and applies the same atom permutation to the other files, so atom i still matches atom i. Writes `*_contig.pdb`.
- `move_proton.py <in.pdb> <out.pdb> --resname R --proton H --to N [--extract lig.pdb]`: moves a proton to another atom in the same residue without changing the atom order, updates CONECT, and can write the residue alone for re-parameterisation
- `select_nac_frames.py <topology.pdb> <traj.dcd>`: picks near-attack frames by nucleophile–C distance and Nu–C=O angle, writing a CSV plus full-system PDBs (`--intact bond_frames.csv` keeps only frames whose bonds `check_bonds.py` found intact)
- `nucleophile_candidates.py <top.pdb> <traj.dcd> --carbon C:1:C4x --oxygen C:1:O2x [--metal ..] [--base ..] [--extra ..]`: candidate nucleophile oxygens (all waters + extras) near an electrophilic carbon: distance, attack angle, metal binding and H-bonds to base candidates per frame, summarised per candidate
- `contact_report.py <topology.pdb> <traj.dcd> --target SEL --sel SEL [--cutoff 3.5] [--near 10] [--frames list]`: polar contacts between chosen residues and a target group, with occupancy (optionally on a subset of frames)
- `strip_solvent.py <in.pdb> <out.pdb> [--keep D] [--keep C:2]`: removes bulk water and monatomic ions (keeping chosen chains/residues) so an MD end point can be the next `preprocess.inputs.protein_pdb`
- `extract_frames.py <topology.pdb> <traj.dcd> --frames 300 400 -1`: writes chosen frames as full-system PDBs (same atom order and ids, the frame's own box in CRYST1), e.g. as `mcube qmmm --pdb` inputs
- `d3_correction.py <structure.pdb> config.yaml [--section S] [--method wb97x] [--damping d3bj]`: post-SCF DFT-D3 energy of the QM region (QM atoms + ASH-style link atoms), to add to single points run without dispersion
- `check_bonds.py <topology.pdb> <traj.dcd> --config config.yaml [--tol 0.25] [--contact 1.75]`: checks that bond geometry in the QM region (or `--sel`) is preserved: stretched bonds, new heavy-atom contacts and proton transfers per frame, with an `intact` flag per frame in `bond_frames.csv`
- `site_report.py <topology.pdb> <traj.dcd> [--metal SEL] [--pair A B] [--angle A B C] [--rmsd-sel SEL]`: metal coordination shell (bound fraction, distances), chosen distances/angles and RMSD over a trajectory (or a single PDB as the trajectory). Atoms as mdtraj selections, `chain:resid:atom` or indices.
- `build_product_box.py <reactant_solvated.pdb> <product_solute.pdb> <out.pdb>`: builds a solvated NEB product in the reactant's water box. NEB endpoints must share every atom, waters included, so solvating reactant and product separately does not work. The product solute must have the same atoms in the same order as the reactant solute. The script aligns on the non-reactive solute atoms and reports water clashes. QM/MM-optimise both endpoints before running NEB.

## Notes and limitations

- **Atom order.** QM atoms, constraints and `actregion_origin` are atom indices into `processed_pdb`. Any structure passed to later stages must keep that atom order, which OpenMM and mdtraj outputs from this pipeline do.
- **Fixed-distance metal constraints.** Metal–donor constraints hold each bond at its initial length during classical MD. In QM/MM, the QM region is free.
- **NPT with implicit solvent.** NPT is not possible with implicit solvent and is skipped automatically.
- **No checkpointing.** Each stage reruns from the start. Results are saved after each step, so partial results survive a crash.
