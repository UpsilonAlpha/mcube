# mcube

A pipeline for protein–ligand systems, often metalloproteins, built on [ASH](https://github.com/RagnarB83/ash). It runs:

1. **preprocess**: ligand parameterisation (OpenFF), protonation and solvation (explicit water, implicit solvent or an explicit lipid bilayer), metal-coordination restraints and QM-region selection
2. **md**: OpenMM minimisation, optional gentle warm-up, NPT box equilibration and production MD
3. **qmmm**: xTB/MM MD and geometry optimisation
4. **scan**: relaxed QM/MM distance scan (builds intermediates/products in the reactant's box)
5. **neb**: QM/MM nudged elastic band between a reactant and a product (optionally NEB + TS optimisation), with a frequency check on the saddle point
6. **ts**: TS optimisation with a partial Hessian, frequency check and IRC (quasi-IRC along the TS mode, or Sella)
7. **sp**: PySCF QM/MM single point on a finished geometry, with density, MEP and HOMO/LUMO cubes
8. **postprocess**: MD data plots, RMSD, RMSF (as B-factors) and a residue interaction network

The path through them is linear, with `sp` as the scoring step applied to each
geometry you want an energy for:

```
preprocess -> md -> qmmm -> sp                    the complex, and its energy
                      \-> scan | neb -> ts -> sp   the path, its TS, and the barrier
```

`sp` is not a step the next stage consumes: it reads a finished geometry and
writes an energy, so it is run once per stationary point (reactant, products,
TS) rather than once per pipeline. `scan` and `neb` are alternatives — a relaxed
scan to build endpoints along one coordinate, or a band between endpoints you
already have.

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
mcube sp config.yaml --pdb QM_MM_optimized.pdb [--section sp_product]
mcube postprocess config.yaml [--qmmm] [--redraw] [--cutoff 2.6] [--metals ZN MG] [--occupancy 0.1]
```

`python -m mcube ...` also works. Every stage writes its files to the current directory, so run each stage from its own working directory. Paths saved in `results` are absolute, so later stages can be run from anywhere.

`--cores` defaults to `$SLURM_CPUS_PER_TASK` (or 1). The `md`, `qmmm`, `scan` and `ts` stages take their core count from `cores` in the config instead, and `sp` from `sp_cores`.

**Threads.** You choose `--cpus-per-task`, `--mem` and `--time`; nothing else.
`mcube env <config> <stage> [--section S]` prints the thread settings that stage
needs and the job script applies them before anything starts, so a stage cannot
run under the wrong policy. It also warns when the allocation is larger than the
stage can use. `docs/csr/run.slurm` is a ready-made job script that does this.

mcube takes the thread *count* for a stage from its `cores` (or `sp_cores` /
`dft_cores`), else `$OMP_NUM_THREADS`, else `$SLURM_CPUS_PER_TASK`, else 1.

## Config file

The config file has one section per stage, and each section has `inputs` (written by you) and `results` (written by mcube). When an input is missing, the stage falls back to the earlier stage's results: for example, `md` uses `preprocess.results.processed_pdb` unless `md.inputs.pdb_file` is set. Start from a template in [`config_examples/`](config_examples):

- [`pdbqt_config.yaml`](config_examples/pdbqt_config.yaml): explicit solvent, docked ligand from a PDBQT
- [`implicit_config.yaml`](config_examples/implicit_config.yaml): implicit solvent (GBn2), ligand from a PDB
- [`membrane_config.yaml`](config_examples/membrane_config.yaml): membrane protein in an explicit POPC bilayer (see [`docs/membrane-validation.md`](docs/membrane-validation.md) for the area-per-lipid and barostat checks)
- [`qmmm_dft_refine.yaml`](config_examples/qmmm_dft_refine.yaml): enzyme mechanism explored at xTB, then stationary points, frequencies and IRC refined at DFT/MM (see [docs/qm-levels.md](docs/qm-levels.md))
- [`scan_2d.yaml`](config_examples/scan_2d.yaml): 2D relaxed scan, bond formation coupled to a proton transfer

Results are written back with `yaml.safe_dump`, so comments in the config are lost the first time a stage saves its results. Writes are locked (POSIX), so several jobs can run stages on one config at once as long as each writes its own `--section`.

**The metal coordination a force field has no bonded parameters for is held by harmonic distance restraint.** Every metal–donor contact `preprocess` detects within `cutoff` is restrained at the distance the input structure already has, for the whole MD stage, and released for `qmmm`. This is the [ASH workflow](https://ash.readthedocs.io/en/latest/Metalloprotein-I.html) mcube follows, which sanctions the substitution directly — *"we could also have added harmonic bond restraints instead of rigid constraints"* — and it is how Gervasoni et al. ([Proteins 2022, 90, 372](https://doi.org/10.1002/prot.26227)) hold a dizinc metallo-β-lactamase site through 100 ns of NVT.

Restraints rather than rigid constraints, for three reasons: many constraints sharing one metal over-determine the constraint network and minimisation diverges; OpenMM refuses any constraint touching a frozen atom, which rules out freezing a donor that carries a hydrogen; and a restraint needs no element-specific parameters, so the same rule works for any metal without the parameterisation a bonded, dummy-cation or 12-6-4 model needs. Detection is per atom, so both oxygens of a carboxylate are restrained and it cannot rotate one in.

Expect some residual drift even so. Gervasoni report the coordination number still rising by recruiting waters *with* restraints on, and rely on the QM/MM stage to recover the geometry — which is what `qmmm` does here. `md` prints the coordination sphere before and after the run, so a contact that appeared or vanished is visible rather than inherited silently.

**A distance you imposed does not silently become a measurement.** `md` records the restraints *you* asked for in `md.results.restraints`, separately from the detected coordination it restrained (`md.results.coordination_restraints`). If `preprocess` is then re-run and its detection would report one of yours as a found contact, it stops and asks which you meant: `constraint_exclude_residues` to leave the contact free, or `acknowledge_restrained` to say the coordination is real. Holding a contact at the distance it already had asserts nothing, so the automatic restraints are not checked — only a distance someone chose can look like structure to a later run.

**Atom selectors.** Wherever an input takes atoms (`md` `restraints`; `actregion_origin`, `opt_constraints`; scan `atoms`/`constraints`; ts `hessian_atoms`/`monitor`; neb `hessian_atoms`/`partial_hessian_atoms`), each atom can be a 0-based index or a `"chain:resid:atom"` string matched against the stage's structure, e.g. `"B:1:ZN"` or `"A:88:NE2"`. Numbers in the same entry (distances, force constants) are left alone, so `restraints: [["B:1:ZN", "C:1:O2x", 2.2, 20.0]]` works. Selectors keep a config valid when the solvent box (and so every index after the solute) changes, and make configs readable for new systems. The resolved indices are printed.

**Chain IDs after solvation.** ASH writes its PDB files with OpenMM's `PDBFile.writeFile` and no `keepIds`, so `processed_pdb` (`finalsystem.pdb`) carries OpenMM's own chain letters: `A`, `B`, `C`, … in chain order, whatever the input used. An input with chains `R` and `P` comes back as `A` and `B`, and in a membrane system the lipids, water and ions take the letters after those. This is true of every solvation mode — it simply goes unnoticed when a single input chain is already `A`. So:
- `residue_variants` and `require_his_variants` are matched against the **input** PDB, before solvation, so they use the input's chain IDs;
- everything resolved after solvation — `extra_qm_residues`, `exclude_residues`, `constraint_exclude_residues`, and the `chain:resid:atom` selectors in `md`, `qmmm`, `scan`, `ts` and `neb` — is matched against `processed_pdb`, so it uses the relabelled letters.

A stale chain ID fails loudly in some places and silently in others: a `chain:resid:atom` selector that matches nothing raises, and so does `extra_qm_residues`, but `exclude_residues` and `constraint_exclude_residues` are only ever tested against the residues that were found, so an entry under a chain letter that no longer exists just does nothing.

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
| `membrane` | `false` | Build an explicit lipid bilayer and water box around the protein (OpenMM's `addMembrane`). Mutually exclusive with `implicit`. **The bilayer is built in the xy plane and the protein is not reoriented**, so `protein_pdb` must already have the membrane normal along z and its hydrophobic belt centred on `membrane_center_z` (e.g. a structure from [OPM](http://opm.phar.umich.edu)). The input's chain IDs do not survive solvation (see **Chain IDs after solvation** above): a receptor/peptide pair given as chains `R` and `P` comes back as `A` and `B`, with the lipids and solvent after them. |
| `membrane_lipidtype` | `POPC` | Lipid: `POPC`, `POPE`, `DLPC`, `DLPE`, `DMPC`, `DOPC` or `DPPC` (OpenMM's pre-equilibrated patches). **The lipid residues are not named after the lipid**: Lipid17 truncates them to three characters, so POPC and POPE are both `POP`, and the others are `DLP`, `DMP`, `DOP`, `DPP`. Select on those in analysis scripts and VMD/PyMOL selections. |
| `membrane_padding` | 10.0 | Smallest distance in Å from any protein atom to a box edge (`addMembrane`'s `minimumPadding`). It is a lower bound, not the box size: the bilayer is tiled from whole copies of OpenMM's patch (63.7 × 66.9 Å for POPC), so the xy dimensions come out as integer multiples of the patch and are usually larger than protein + 2 × padding. It sizes the water layer in z as well. |
| `membrane_center_z` | 0.0 | z of the bilayer midplane in Å |
| `membrane_lipid_xml` | `amber14/lipid17.xml` | Lipid force field, appended to `default_forcefields` on the membrane path (the protein/water defaults have no lipid parameters). Checked against the installed OpenMM. |
| `pH` | 8.0 explicit / 7.0 implicit | Protonation pH |
| `solvent_padding` | 10.0 | Water padding in Å (explicit only) |
| `ionic_strength` | 0.1 | NaCl concentration in M (explicit only) |
| `metals` | none | Element symbols treated as metals. Empty or absent (e.g. protein–protein complexes): no coordination restraints and no QM region. |
| `cutoff` | required with `metals` | Metal–donor distance cutoff (Å) |
| `donor_elements` | `[N, O, S]` | Elements that count as metal donors (matched on element, not atom name) |
| `require_his_variants` | `false` | Fail if any histidine has no tautomer in `residue_variants`. A warning is printed either way. |
| `extra_qm_residues` | none | `{chain: [resid, ...]}`: residues added to the QM region on top of the metal's coordination sphere (e.g. catalytic Lys/Arg) |
| `constraint_exclude_residues` | none | `{chain: [resid, ...]}`: residues whose metal–donor distance is left unconstrained during MD (e.g. the substrate) |
| `acknowledge_restrained` | none | `[["B:1:ZN", "C:1:O2x"], ...]`: pairs in `md.results.restraints` that this run may still constrain. Without it, preprocess **refuses**, because an imposed distance would become an apparently structural one |
| `exclude_residues` | none | `{chain: [resid, ...]}`: coordinating residues left out of the QM region. Their metal–donor restraints are still applied in classical MD. Not the same as `boundary_excluded`, which lists the non-amino-acid QM atoms that ASH should not place link atoms on. |
| `backbone_atoms` | standard backbone | Atoms of amino-acid QM residues that stay in MM |
| `amino_acids` | 20 standard + Amber variants | Residues cut at CA–CB. Any other residue (metals, ligands, waters) is put fully into the QM region and excluded from boundary detection. |

`results`: `processed_pdb` (`finalsystem.pdb`), `ligand_xml`, `ligand_resname`, `qm_residues` (`{chain: {resid: resname}}`), `qm_atoms` and `boundary_excluded` (atom indices into `processed_pdb`), and `coordination` — the metal–donor contacts it detected, as `[[i, j, d_Å, k_kcal/mol/Å²], …]`, which `md` restrains. Both numbers are written out, so the distance each contact is held at and the force constant are visible and editable per pair before `md` runs.

### `md.inputs`

| Key | Default | Description |
| --- | --- | --- |
| `pdb_file`, `ligand_xml`, `default_forcefields` | from preprocess | Overrides |
| `cores` | | OpenMM CPU threads |
| `platform`, `properties` | `CPU` | OpenMM platform, e.g. `CUDA` with `{Precision: mixed}` on a GPU node (see `docs/running-on-bunya.md`) |
| `temp` | | Temperature (K) |
| `implicit` | `false` | GBn2 implicit solvent. Uses HBond constraints, 1.5 amu hydrogen mass and a 2 nm cutoff. NPT and reimaging are skipped. |
| `implicit_eps` | | Solvent dielectric for implicit solvent |
| `membrane` | `preprocess.inputs.membrane` | Explicit bilayer: the system is built here (PME, 12 Å cutoff, HBond constraints, rigid water, 1.5 amu hydrogens) with an `openmm.MonteCarloMembraneBarostat` (XYIsotropic, ZFree, zero surface tension) and passed to ASH as a serialised system, so the area per lipid and the bilayer thickness relax independently. ASH's own barostat is isotropic and would couple them. Equilibration runs through `OpenMM_MD` instead of `OpenMM_box_equilibration` for the same reason, and production keeps the barostat (so it is NPT, not NVT). `temp` is required. |
| `npt_time` | 1000 | Membrane NPT length in ps (`npt_timestep`). Box size and area per lipid per frame are written to `membrane_box.csv` and recorded under `results.membrane`. |
| `surface_tension` | 0.0 | Membrane barostat surface tension in bar·nm. 0 is a tensionless bilayer, which is what reproduces experimental areas per lipid. |
| `pressure` | 1.0 | Barostat pressure in bar (membrane only) |
| `barostat_frequency` | 25 | Steps between barostat moves (membrane only) |
| `nonbonded_cutoff` | 12.0 | Periodic nonbonded cutoff in Å (membrane only; ASH's own default elsewhere) |
| `metal_parameters` | `iod` | Lennard-Jones parameters for the metals. `iod` is the Li/Merz 12-6 IOD set (TIP3P), fitted to reproduce experimental ion–oxygen distances: Amber's stock Zn(II) puts the Zn–O pair minimum at ~1.69 Å against a real bond of 2.0–2.1 Å, and IOD moves it to ~1.89 Å. Covers 40 ions (+2/+3/+4, including every lanthanide) and needs nothing specified — the element and the charge the force field assigned pick the entry, and an untabulated metal is left alone with a warning. `none` keeps the force field's own values. Which set was used is recorded in `md.results.metal_parameters`. It does **not** fix over-coordination: a 12-6 potential cannot tell a ligand from a non-ligand, the coordination restraints do that |
| `gentle` | `false` | Warm-up after minimisation, for a system that will not start otherwise. Despite the name it is 10 steps at 1 K and 50 at 10 K (0.055 ps in total) followed by **40 ps of ordinary 300 K MD**, so after a successful minimisation it buys nothing NPT does not, and the 300 K leg is long enough to let a force-field artefact develop before any equilibration you are watching |
| `skip_npt` | `false` | Skip NPT box equilibration (timestep `npt_timestep`) |
| `skip_nvt` | `false` | Skip production MD |
| `nvt_timestep`, `nvt_time` | | Production timestep and length (ps) |
| `traj_frequency` | 1000 | Production trajectory save interval (steps) |
| `restraints` | `preprocess.results.coordination` | Harmonic distance restraints `[[i, j, d_Å, k_kcal/mol/Å²], …]`, active during minimisation, warm-up, NPT and NVT. Everything `preprocess` detected is already restrained, so this is only for *extra* restraints — pulling a substrate into place, say. An entry naming a pair `preprocess` detected replaces that row. To change how the metal site itself is held, edit `preprocess.results.coordination`: its rows carry both the distance and the force constant |

Minimisation always runs, with the metal-site restraints already active.

`results`: `minimized`, `gentle`, `forcefields`, `trajfile`, `lastframe`, `csv`, `coordination_restraints` — the detected contacts this run restrained, kept apart from `restraints` (the ones you chose) because only a distance you chose can look like structure to a later `preprocess` — and `coordination_start` / `coordination_end`, the metal's coordination sphere before and after the run, with md printing what changed between them. A donor that appears was outside `cutoff` and so was never restrained; widening `cutoff` holds it. In explicit solvent, `trajfile` and `lastframe` are the reimaged files.

### `qmmm.inputs`

| Key | Default | Description |
| --- | --- | --- |
| `pdb_file` | `md.results.lastframe` | Starting structure. It must have the same atom order as `processed_pdb`. |
| `forcefields` | `md.results.forcefields` | |
| `qm_residues` | from preprocess | Override the QM region (rebuilt on `pdb_file`) |
| `extra_qm_residues` | none | `{chain: [resid, ...]}` added to the preprocess QM region for this run only (e.g. second-shell waters). Also accepted in `scan`, `neb` and `ts` sections. |
| `implicit` | `md.inputs.implicit` | Non-periodic MM when true |
| `charge`, `mult` | required | QM-region charge and multiplicity |
| `cores` | | xTB / PySCF cores |
| `xtb_method` | `GFN2` | |
| `xtb_runmode` | `inputfile` | `library` avoids file I/O at every step (needs xtb-python) |
| `temp`, `timestep`, `time` | | QM/MM MD settings (K, ps, ps) |
| `traj_frequency` | 1000 | QM/MM MD save interval (steps) |
| `skip_equilibrate` | `false` | Skip QM/MM MD |
| `skip_optimize` | `false` | Skip QM/MM optimisation |
| `opt_constraints` | none | Bond constraints during the optimisation: `[[i, j], …]` holds the starting distance, `[[i, j, d_Å], …]` sets it. Atoms must be in the active region. MM waters in the active region are kept rigid; waters **in the QM region are left free**, so a QM water can donate or accept a proton. |
| `opt_maxiter` | 200 | Optimisation steps. If they run out, the last geometry is kept as `QM_MM_optimized.pdb` and `converged: false` is recorded. |
| `actradius` | 0 | If > 0, optimise all atoms within this radius (Å) of `actregion_origin`. Otherwise only the QM atoms move. |
| `actregion_origin` | first non-standard QM atom | Atom index at the centre of the active region |
| `qm_engine` | `xtb` | QM theory for the gradient stages (optimisation, scan, NEB, TS, frequencies): `xtb` for GFN2-xTB, or `dft` to run them at DFT/MM. A DFT gradient costs orders of magnitude more than an xTB one, so use `dft` to refine stationary points an xTB path has already located, not to explore. Settable per `scan`/`neb`/`ts` section too. |
| `dft_functional`, `dft_basis`, `dft_gridlevel`, `dft_densityfit`, `dft_auxbasis`, `dft_platform`, `dft_memory` | the matching `sp_*` value | Level of theory for `qm_engine: dft`. Each falls back to its `sp_*` counterpart, so a config that already sets the single-point level need not repeat it. |
| `dft_scf_dispersion` | none | Dispersion inside the SCF for `qm_engine: dft` (ASH/PySCF's own keyword: `D3`, `D4`, `VV10`, …), so it enters the gradient. Distinct from `sp_dispersion`, which is applied after a single point. |
| `dft_cores` | stage `cores` | Threads for the DFT gradients |
| `dft_gpu_pcgrad` | `false` | With `dft_platform: GPU`, also compute the MM point-charge gradient on the device. For a large MM environment that is most of the QM/MM gradient. |
`results` (section `qmmm`, or `--section NAME`): `qm_atoms`, `boundary_excluded`, `trajfile`, `lastframe`, `opttraj`, `optimized_pdb` (`QM_MM_optimized.pdb`, the full optimised system, usable as an NEB endpoint), `opt_energy` (Eh), `converged`, `opt_constraints` (with the distances used).

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
| `atoms2`, `start2`, `end2`, `points2` | none | A second distance: makes the scan a 2D grid (e.g. nucleophile–C × proton–acceptor), visited row by row over coordinate 2 with coordinate 1 swept back and forth, so each point starts from a neighbour. PDBs are `scan_KK_XX.pdb`; the CSV gains `i1`, `i2`, `target2_A`, `distance2_A`. |
| `resume` | `true` | Continue from an existing `scan_energies.csv` (skip finished points, start from the last structure) |
| `extra_qm_residues` / `qm_residues` | qmmm values | QM region for this scan (e.g. add second-shell waters) |
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
| `irc_dx` | 0.1 | Sella IRC step size, A (`irc_method: sella` only) |
| `opt_maxiter` | 300 | Steps per IRC minimisation |
| `rigid_water` | `true` | Keep MM active-region waters rigid, as in the `qmmm` minimisations. Waters in the QM region are never constrained, so they stay free to transfer protons. |

`results`: `ts_pdb` (`ts_optimized.pdb`), `ts_energy`, `guess_imaginary`, `imaginary_frequencies`, `small_imaginary`, `lowest_frequencies`, `ts_mode_atoms`, monitored distances, and `irc.forward/reverse` (`irc_*.pdb`, energy, barrier = E(TS) − E(end) in kcal/mol, distances).

### `sp.inputs`

The QM/MM single point, run as its own stage (`mcube sp`). It is separate from
`qmmm` because the two want different machines: a DFT SCF scales across cores
and wants tens of gigabytes (measured: 96% CPU on 16 cores, 56 GB for ~120 QM
atoms), while the xTB MD and optimisation in front of it saturate about four
threads and use under a gigabyte. One job cannot size both, so submit them
separately.

Inputs are read from `qmmm.inputs`, then `sp.inputs`, then the `--section`
inputs, so the force fields, charge, multiplicity and QM region are inherited
from the model being scored and only the level and the resources need stating.
A structure is required — `--pdb` or `pdb_file` — because there is no sensible
default geometry for a single point.

| Key | Default | Description |
| --- | --- | --- |
| `sp_functional`, `sp_basis` | `wb97x-v`, `def2-svp` | Single-point level of theory (PySCF names, e.g. `wb97m_v`; the energy is saved as `sp_energy`) |
| `sp_gridlevel` | 3 | PySCF DFT grid level (ASH's own default is 5, which made VV10 functionals very slow) |
| `sp_densityfit`, `sp_auxbasis` | `true`, PySCF default | Density fitting (RI-JK for hybrids) |
| `sp_cores` | physical cores | Threads for the single point (independent of the xTB `cores`) |
| `sp_platform` | `CPU` | `GPU` uses gpu4pyscf (must be installed) |
| `sp_memory` | 80% of `$SLURM_MEM_PER_NODE`, else 3100 | PySCF `max_memory` (MB). Density-fitting integrals that don't fit are written to disk and re-read every cycle; a range-separated hybrid on ~120 QM atoms with def2-SVP needs ~60 GB to stay in RAM. |
| `sp_dispersion`, `sp_dispersion_method` | none, `sp_functional` | Post-SCF DFT-D3 on the QM region (link atoms included): `d3bj` or `d3zero`, parameters for the named functional (simple-dftd3 names, e.g. `wb97x`). Needs `dftd3-python`. Saved as `sp_disp_energy` and `sp_energy_disp`. |
| `cubes` | `true` | Write `density.cube`, `mep.cube`, `HOMO.cube` and `LUMO.cube` |

`results` (section `sp`, or `--section NAME`): `pdb_file`, `qm_atoms`, `boundary_excluded`, `sp_energy` (Eh), `sp_level`, `sp_dispersion`, `sp_disp_energy`, `sp_energy_disp`, `cubes`.

```bash
# the usual pair: xTB QM/MM on a modest allocation, then the SCF on a big one
sbatch --cpus-per-task=8  --mem=24G  ... mcube qmmm config.yaml --section qmmm_reactant
sbatch --cpus-per-task=64 --mem=192G ... mcube sp   config.yaml --section sp_reactant        --pdb protdft/qmmm_reactant/QM_MM_optimized.pdb
```

## Outputs by stage

| Stage | Main files |
| --- | --- |
| preprocess | `<resname>.mol/.pdb`, `openff_<resname>.xml`, `merged.pdb`, `finalsystem.pdb` |
| md | `frag-minimized.pdb`, `warmup_MD_cycle*`, `relaxbox_NPT.dcd`, `nptsim.csv`, `NVTtrajectory(.dcd, _imaged.dcd, _lastframe(_imaged).pdb)`, `nvtsim.csv`, `system.xml` (implicit) |
| qmmm | `QM_MM.dcd`, `QM_MM_lastframe.pdb`, `active_atoms`, `geometric_OPTtraj_Full.xyz`, `geometric_OPTtraj_QMregion.xyz` |
| sp | `*.cube` (`density`, `mep`, `HOMO`, `LUMO`) |
| scan | `scan_XX.pdb`, `scan_energies.csv`, `active_atoms` |
| neb | ASH NEB files, `saddlepoint.xyz`, `saddlepoint.pdb` |
| ts | `Hessian_partial_TS`, `ts_optimized.pdb`, `irc_forward.pdb`, `irc_reverse.pdb`, `irc_*_opttraj.xyz` (or Sella's `sella_*.xyz`) |
| postprocess | `Temperature.png`, `Density.png`, `BoxVolume.png`, `RMSD.png`, `interaction_network.csv/.png`, `rmsf_bfactor.pdb` |

## postprocess

By default, postprocess analyses the MD trajectory (`md.results`). With `--qmmm`, it analyses the QM/MM trajectory and also converts the optimisation trajectory to `output.dcd`. It computes:

- **Plots** of temperature, density and box volume from the MD CSV (MD only), and an RMSD plot of the protein heavy atoms
- **Interaction network** with H-bonds (protein and ligand), salt bridges (< 4 Å), π–π (< 5 Å), cation–π (< 6 Å) and metal coordination (`--cutoff`). Edges whose occupancy is above `--occupancy` are kept. `--redraw` redraws the plot from `interaction_network.csv` without recomputing. `--spring_constant` and `--iterations` control the layout.
- **RMSF** of Cα atoms after alignment, written as B-factors to `rmsf_bfactor.pdb`

## Choosing the QM level

Every gradient-driven stage (`qmmm` optimisation, `scan`, `neb`, `ts` and its frequencies) runs at GFN2-xTB by default, or at DFT/MM with `qm_engine: dft`. See [docs/qm-levels.md](docs/qm-levels.md) for which to use where, how the `dft_*` keys relate to the `sp_*` ones, why frequencies need a partial Hessian at DFT, and how to install and use GPU PySCF without breaking numpy.

`qm_engine` and the `dft_*` / `sp_*` keys may be set on `qmmm.inputs` or on an individual `scan` / `neb` / `ts` section, which overrides it; everything else about the QM/MM setup (charge, multiplicity, force fields) always comes from `qmmm.inputs`.

## Membrane simulations

`membrane: true` runs the unchanged `preprocess` and `md` stages on a protein in an explicit lipid bilayer. The config is the one in [`config_examples/membrane_config.yaml`](config_examples/membrane_config.yaml): `preprocess.inputs.membrane` builds the bilayer, `md.inputs.membrane` (which defaults to the preprocess value) runs it, and every other key keeps its meaning and default. Nothing else in the config selects the membrane mode.

```bash
mcube preprocess config.yaml   # CPU: bilayer, water and ions around the protein
mcube md config.yaml           # GPU (platform: CUDA): minimise, membrane NPT, production
```

- **Input.** `protein_pdb` (and a `ligand`, or a peptide as a second chain) must already be oriented: membrane normal along z, hydrophobic belt centred on `membrane_center_z`. mcube never rotates or shifts the protein. The [OPM database](http://opm.phar.umich.edu) publishes oriented structures.
- **Defaults** are those of the explicit-solvent path (PME, 12 Å cutoff, HBond constraints, rigid water, 1.5 amu hydrogens, pH 8.0 and 0.1 M NaCl unless set) with amber14 + Lipid17 and a POPC bilayer with 10 Å padding. What you must give is `temp` and, for a usable speed, `platform: CUDA`.
- **Ensemble.** The system carries one `MonteCarloMembraneBarostat` (XYIsotropic, ZFree, tensionless, 1 bar), so the area per lipid and the thickness relax independently. NPT equilibration and production both keep it, so production is NPT, not NVT. ASH's isotropic barostat is never added, and `md` checks that exactly one barostat is present before and after every stage. [`docs/membrane-validation.md`](docs/membrane-validation.md) has the area-per-lipid and thickness checks behind this.
- **Outputs** are the usual ones plus `membrane_box.csv` (box and area per lipid per NPT frame) and `membrane_box_production.csv`, summarised in `md.results.membrane` and `md.results.membrane_production`. `trajfile` and `lastframe` are the reimaged files.
- **Segments and restarts.** `md` starts from `md.inputs.pdb_file` (default: `preprocess.results.processed_pdb`). To continue a run, or split it into short jobs that fit a queue, point `pdb_file` at the previous last frame and set `skip_npt: true`. `relaxbox_NPT_lastframe.pdb` is the equilibrated box (equilibrate alone with `skip_nvt: true`) and `NVTtrajectory_lastframe.pdb` the end of production. Each segment minimises again and draws new velocities, so it restarts the structure, not the trajectory.
- **A new partner in an equilibrated bilayer.** To compare many peptides or designed binders in one bilayer without rebuilding it, equilibrate the receptor once, write each partner's coordinates (superposed onto the receptor) into a copy of that last frame, remove the solvent that overlaps it, re-neutralise, and set `md.inputs.pdb_file` to the result. The minimisation then relaxes the remaining contacts. The partner's chain follows the water in that file, so select it by its chain letter.
- **Other modes.** `membrane` and `implicit` are mutually exclusive and raise at once. With `membrane` unset the explicit and implicit paths are unchanged (`tests/test_membrane_pipeline.py` asserts it). `qmmm` and the later stages read the structure `md` wrote and do not use the membrane barostat; a QM/MM region in a membrane system has not been run.
- **Analysis.** `scripts/membrane_report.py` (below) reads the trajectory of any membrane run.

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
- `membrane_report.py <topology.pdb> <traj.dcd ...> --out PREFIX [--receptor A] [--partner B] [--site RES ...] [--core-exclude LO-HI] [--footprint FILE] [--track N=RES] [--distance N=R1,R2] [--angle N=R1,R2,R3]`: membrane health per frame for any membrane run (box, area per lipid, P-P thickness, leaflet imbalance, lipids and water inside the bundle). With a receptor and a partner chain it adds the partner's pose: RMSF in the receptor frame and internal, RMSD, drift, contacts and their retention, and occupancy of chosen receptor residues and of a reference contact set (`--write-footprint` on a reference run writes one). Receptor readouts are named on the command line, so the same call works for any protein.
- `site_report.py <topology.pdb> <traj.dcd> [--metal SEL] [--pair A B] [--angle A B C] [--rmsd-sel SEL]`: metal coordination shell (bound fraction, distances), chosen distances/angles and RMSD over a trajectory (or a single PDB as the trajectory). Atoms as mdtraj selections, `chain:resid:atom` or indices.
- `build_product_box.py <reactant_solvated.pdb> <product_solute.pdb> <out.pdb>`: builds a solvated NEB product in the reactant's water box. NEB endpoints must share every atom, waters included, so solvating reactant and product separately does not work. The product solute must have the same atoms in the same order as the reactant solute. The script aligns on the non-reactive solute atoms and reports water clashes. QM/MM-optimise both endpoints before running NEB.

## Notes and limitations

- **Atom order.** QM atoms, restraints and `actregion_origin` are atom indices into `processed_pdb`. Any structure passed to later stages must keep that atom order, which OpenMM and mdtraj outputs from this pipeline do.
- **Metal-site restraints.** Metal–donor contacts are held near their initial length by harmonic restraint during classical MD. In QM/MM, the QM region is free.
- **NPT with implicit solvent.** NPT is not possible with implicit solvent and is skipped automatically.
- **No checkpointing.** Each stage reruns from the start. Results are saved after each step, so partial results survive a crash.
