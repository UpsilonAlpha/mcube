# mcube

A pipeline for protein–ligand systems, often metalloproteins, built on [ASH](https://github.com/RagnarB83/ash). It runs:

1. **preprocess**: ligand parameterisation (OpenFF), protonation and solvation (explicit or implicit), metal-coordination bond constraints and QM-region selection
2. **md**: OpenMM minimisation, optional gentle warm-up, NPT box equilibration and production MD
3. **qmmm**: xTB/MM MD and geometry optimisation, with an optional PySCF QM/MM single point and cube files
4. **neb**: QM/MM nudged elastic band between a reactant and a product, with a frequency check on the saddle point
5. **postprocess**: MD data plots, RMSD, RMSF (as B-factors) and a residue interaction network

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
mcube neb config.yaml --reactant reactant.pdb --product product.pdb [--cores N]
mcube postprocess config.yaml [--qmmm] [--redraw] [--cutoff 2.6] [--metals ZN MG] [--occupancy 0.1]
```

`python -m mcube ...` also works. Every stage writes its files to the current directory, so run each stage from its own working directory. Paths saved in `results` are absolute, so later stages can be run from anywhere.

`--cores` defaults to `$SLURM_CPUS_PER_TASK` (or 1). The `md` and `qmmm` stages take their core count from `cores` in the config instead.

## Config file

The config file has one section per stage, and each section has `inputs` (written by you) and `results` (written by mcube). When an input is missing, the stage falls back to the earlier stage's results: for example, `md` uses `preprocess.results.processed_pdb` unless `md.inputs.pdb_file` is set. Start from a template in [`config_examples/`](config_examples):

- [`pdbqt_config.yaml`](config_examples/pdbqt_config.yaml): explicit solvent, docked ligand from a PDBQT
- [`implicit_config.yaml`](config_examples/implicit_config.yaml): implicit solvent (GBn2), ligand from a PDB

Results are written back with `yaml.safe_dump`, so comments in the config are lost the first time a stage saves its results.

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
| `metals` | required | Element symbols treated as metals |
| `cutoff` | required | Metal–donor distance cutoff (Å) |
| `donor_elements` | `[N, O, S]` | Elements that count as metal donors (matched on element, not atom name) |
| `exclude_residues` | none | `{chain: [resid, ...]}`: coordinating residues left out of the QM region. Their metal–donor bond constraints are still applied in classical MD. Not the same as `boundary_excluded`, which lists the non-amino-acid QM atoms that ASH should not place link atoms on. |
| `backbone_atoms` | standard backbone | Atoms of amino-acid QM residues that stay in MM |
| `amino_acids` | 20 standard + Amber variants | Residues cut at CA–CB. Any other residue (metals, ligands, waters) is put fully into the QM region and excluded from boundary detection. |

`results`: `processed_pdb` (`finalsystem.pdb`), `ligand_xml`, `ligand_resname`, `qm_residues` (`{chain: {resid: resname}}`), `qm_atoms` and `boundary_excluded` (atom indices into `processed_pdb`), and `constraints` (metal–donor atom-index pairs).

### `md.inputs`

| Key | Default | Description |
| --- | --- | --- |
| `pdb_file`, `constraints`, `ligand_xml`, `default_forcefields` | from preprocess | Overrides |
| `cores` | | OpenMM CPU threads |
| `temp` | | Temperature (K) |
| `implicit` | `false` | GBn2 implicit solvent. Uses HBond constraints, 1.5 amu hydrogen mass and a 2 nm cutoff. NPT and reimaging are skipped. |
| `implicit_eps` | | Solvent dielectric for implicit solvent |
| `gentle` | `false` | Gentle warm-up (1 → 10 → 300 K) after minimisation |
| `skip_npt` | `false` | Skip NPT box equilibration (timestep `npt_timestep`) |
| `skip_nvt` | `false` | Skip production MD |
| `nvt_timestep`, `nvt_time` | | Production timestep and length (ps) |
| `traj_frequency` | 1000 | Production trajectory save interval (steps) |

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
| `actradius` | 0 | If > 0, optimise all atoms within this radius (Å) of `actregion_origin`. Otherwise only the QM atoms move. |
| `actregion_origin` | first non-standard QM atom | Atom index at the centre of the active region |
| `singlepoint` | `true` | PySCF QM/MM single point (RKS/UKS from `mult`), plus `density.cube`, `mep.cube`, `HOMO.cube` and `LUMO.cube` |
| `sp_functional`, `sp_basis` | `wb97x-v`, `def2-svp` | Single-point level of theory |

`results`: `qm_atoms`, `boundary_excluded`, `trajfile`, `lastframe`, `opttraj`, `cubes`.

### `neb.inputs` (optional section)

| Key | Default | Description |
| --- | --- | --- |
| `images` | 8 | Number of NEB images |
| `actradius`, `actregion_origin` | qmmm values | Active region. With no radius, only the QM atoms move. |
| `numfreq` | `true` | Numerical frequencies (QM atoms only) on the saddle point |

`neb` uses `qmmm.inputs` for `charge`, `mult`, `xtb_method`, `forcefields` and the QM region. `--cores` is split between images that run in parallel and xTB within each image.

`results`: `qm_atoms`, `saddlepoint_xyz`, `saddlepoint_energy`, `frequencies`.

## Outputs by stage

| Stage | Main files |
| --- | --- |
| preprocess | `<resname>.mol/.pdb`, `openff_<resname>.xml`, `merged.pdb`, `finalsystem.pdb` |
| md | `frag-minimized.pdb`, `warmup_MD_cycle*`, `relaxbox_NPT.dcd`, `nptsim.csv`, `NVTtrajectory(.dcd, _imaged.dcd, _lastframe(_imaged).pdb)`, `nvtsim.csv`, `system.xml` (implicit) |
| qmmm | `QM_MM.dcd`, `QM_MM_lastframe.pdb`, `active_atoms`, `geometric_OPTtraj_Full.xyz`, `geometric_OPTtraj_QMregion.xyz`, `*.cube` |
| neb | ASH NEB files, `saddlepoint.xyz` |
| postprocess | `Temperature.png`, `Density.png`, `BoxVolume.png`, `RMSD.png`, `interaction_network.csv/.png`, `rmsf_bfactor.pdb` |

## postprocess

By default, postprocess analyses the MD trajectory (`md.results`). With `--qmmm`, it analyses the QM/MM trajectory and also converts the optimisation trajectory to `output.dcd`. It computes:

- **Plots** of temperature, density and box volume from the MD CSV (MD only), and an RMSD plot of the protein heavy atoms
- **Interaction network** with H-bonds (protein and ligand), salt bridges (< 4 Å), π–π (< 5 Å), cation–π (< 6 Å) and metal coordination (`--cutoff`). Edges whose occupancy is above `--occupancy` are kept. `--redraw` redraws the plot from `interaction_network.csv` without recomputing. `--spring_constant` and `--iterations` control the layout.
- **RMSF** of Cα atoms after alignment, written as B-factors to `rmsf_bfactor.pdb`

## Standalone scripts

Both scripts are in `scripts/`:

- `plot_md_data.py <data.csv>`: plots temperature (plus density and volume for NPT) from an OpenMM CSV
- `plot_rmsd_via_mdtraj.py <traj.dcd> <topology.pdb>`: plots RMSD of the protein heavy atoms
- `build_product_box.py <reactant_solvated.pdb> <product_solute.pdb> <out.pdb>`: builds a solvated NEB product in the reactant's water box. NEB endpoints must share every atom, waters included, so solvating reactant and product separately does not work. The product solute must have the same atoms in the same order as the reactant solute. The script aligns on the non-reactive solute atoms and reports water clashes. QM/MM-optimise both endpoints before running NEB.

## Notes and limitations

- **Atom order.** QM atoms, constraints and `actregion_origin` are atom indices into `processed_pdb`. Any structure passed to later stages must keep that atom order, which OpenMM and mdtraj outputs from this pipeline do.
- **Fixed-distance metal constraints.** Metal–donor constraints hold each bond at its initial length during classical MD. In QM/MM, the QM region is free.
- **NPT with implicit solvent.** NPT is not possible with implicit solvent and is skipped automatically.
- **No checkpointing.** Each stage reruns from the start. Results are saved after each step, so partial results survive a crash.
