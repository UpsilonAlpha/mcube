# Ion-channel permeation PMF: decision log

A running record of the modelling decisions for the umbrella-sampling study of Ca²⁺ and two
lanthanide ions passing through the channel in `docs/channel/sequence.fasta`. Each entry gives
the decision, why, the evidence, and whether to revisit it. Same format as
[csr-decision-log.md](csr-decision-log.md). The protocol is in
[channel-us-plan.md](channel-us-plan.md).

**Aim.** Potentials of mean force for a light and a heavy lanthanide and for Ca²⁺ crossing the
channel, built with mcube's existing stages as closely as possible.

---

## System

**C1 – Fold the mature chain, residues 110–373, not all 373.** (2026-10-22; open)
- Evidence: 85–109 (`MRARAQTFWIMAACTFALAPTCVHA|EEG`) reads as a Sec signal peptide: a
  positive n-region, a hydrophobic h-region and an AXA-like cleavage site. 1–84 are upstream of
  a second Met and look like a mis-called start codon (3 of the 5 Cys, Arg-rich). The mature
  264 aa has no Cys, is Gly-rich (12.5 %) and ends in Phe, the C-terminal aromatic typical of
  Gram-negative outer-membrane β-barrels.
- Decision: fold `channel_mature`. Treat the protein as an outer-membrane β-barrel (orientation,
  lipid choice and oligomer state follow from that).
- Inference, not verified: no BLAST/SignalP was run, so the sequence has not been sent to an
  outside service. Confirming the family and cleavage site is the user's call.
- Revisit: yes, once the AF2 model is in (barrel strand count, open lumen).

**C2 – Fold monomer and homotrimer with AF2 2.3.2 via `alphafold_custom`.** (2026-10-22)
- Evidence: general porins are trimers; a smaller OMP can be a monomer. RCC's
  `alphafold_multimer` passes `.../databses/pdb_seqres/...` (typo; binder-physics TOOLS.md §4.5),
  so it fails at template search.
- Decision: `alphafold_custom` with the corrected `pdb_seqres` path, monomer and 3-mer; choose by
  pLDDT/ipTM and whether a continuous lumen forms. Chai-1 (MSA-free, weights staged) as a second
  opinion.

**C3 – Bilayer: POPE, single lipid type, at 310 K.** (2026-10-22; open)
- mcube already accepts `membrane_lipidtype: POPE` (`membrane.py:LIPID_TYPES`), but only POPC
  has been validated (`docs/membrane-validation.md`). A *mixed* POPE/POPC bilayer is not
  supported: `addMembrane` takes one lipid type.
- POPE's main phase transition is near 25 °C, so 300 K sits on it; run at 310 K.
- The outer membrane's LPS outer leaflet cannot be built by mcube; POPE is the usual stand-in.
- Revisit: whether to add a mixed bilayer (user: "may have to add POPE in addition to POPC"),
  and a POPE area-per-lipid check (literature ~56–60 Å² at 310 K) before production.

**C4 – Ions: Ca²⁺, La³⁺ and Lu³⁺.** (user chose La³⁺ vs Lu³⁺, 2026-10-22)
- Gd³⁺ is 4f⁷, multiplicity 8. Dy³⁺ does not fix this: it is 4f⁹ (multiplicity 6) with orbital
  angular momentum, which is harder, not easier. Lu³⁺ (4f¹⁴) is closed shell and is the heaviest
  lanthanide, so La³⁺ vs Lu³⁺ spans the whole lanthanide contraction with two singlets.
- Recommendation: **La³⁺ vs Lu³⁺**, with Y³⁺ (Shannon CN6 0.900 Å, beside Dy³⁺ 0.912) as an
  optional f-free control. The user chose La³⁺ vs Lu³⁺.
- Ca²⁺ cannot use `metal_parameters: 1264` yet (C12: the +2 12-6-4 rows were not available to
  check), so the Ca²⁺ PMF runs on IOD unless that table is added.
- Checked (2026-10-22): Gd³⁺ ⁸S₇/₂, Dy³⁺ ⁶H₁₅/₂; La³⁺/Lu³⁺ ¹S₀ by configuration. Aqua ions:
  nine-coordinate through the early series, ~8.2 at Lu (Persson et al., *Chem. Eur. J.* 2008,
  10.1002/chem.200701281; Sessa et al., *J. Chem. Phys.* 2016, 10.1063/1.4951714). Those
  EXAFS numbers are the validation targets for the xTB/MM hydration shell.
- **GFN2-xTB treats every lanthanide's 4f shell as core** (5d6s6p basis, no f functions; xtb
  `src/xtb/gfn2.f90`), so all Ln³⁺ are closed shell in xTB: `mult: 1` for any choice. The
  spin problem exists only at the DFT reweighting step, and for Lu not at all.
- DFT basis: PySCF's bundled def2 files have La but no Ce–Lu; Lu needs the `basis_set_exchange`
  package (def2-SVP, 28-electron ECP, 4f in valence). Installing it into the ASH env follows the
  `pip install --dry-run` rule.

## Method

**C5 – Hybrid: MM umbrella sampling over the full path, QM/MM umbrella sampling over the
constriction, MM→QM/MM reweighting, then xTB→DFT single-point reweighting.** (user, 2026-10-22)
- Why: QM/MM windows can afford ~10–20 ps at ~1.4–2 s per step (`running-on-bunya.md`), while
  first-shell water exchange for these ions is on the ns scale, so QM/MM alone cannot sample
  dehydration. MM supplies the sampling; QM/MM corrects the chemistry where it matters.
- **Primary protocol: Li, Jia, Pan, Shao, Mei, *JCTC* 2018, 14, 5583 (10.1021/acs.jctc.8b00571)**
  — umbrella sampling at a cheap QM/MM level reweighted to a high level by weighted
  thermodynamic perturbation along a CV, i.e. our xTB/MM → DFT/MM step — with Part 4, Wang et al.
  *JCTC* 2021, 17, 1318 (10.1021/acs.jctc.0c01149), for a QM region whose membership changes.
- MM leg: Allen, Andersen & Roux, *Biophys. J.* 2006, 90, 3447 (10.1529/biophysj.105.077073) for
  the channel PMF protocol (CV, cylinder, spacing, convergence); Im & Roux, *J. Mol. Biol.* 2002
  (10.1016/S0022-2836(02)00778-7) for OmpF, the closest system. WHAM: Kumar et al. 1992
  (10.1002/jcc.540130812); MBAR: Shirts & Chodera 2008 (10.1063/1.2978177); uncertainty:
  Grossfield et al., *LiveCoMS* 2019.
- **Stitching (refines C5):** following Heimdal & Ryde, *PCCP* 2012 (10.1039/C2CP41005B), MM→QM/MM
  perturbation converges poorly over a whole path. So the constriction is sampled directly at
  xTB/MM, and joined to the MM PMF by one offset at each end from MM↔QM/MM BAR in the overlap
  windows, with overlap diagnostics reported. Back-up estimators if overlap is poor:
  non-Boltzmann Bennett (König et al., *JCTC* 2014, 10.1021/ct401118k), force-matched reference
  potentials (Giese & York, *JCTC* 2019, 10.1021/acs.jctc.9b00401).
- No published xTB/MM or DFTB umbrella-sampling PMF of an ion through a channel was found.
- All DOIs above were checked against Crossref on 2026-10-22.

**C6 – One independent simulation per window.** (user, 2026-10-22)
- Each window is its own job and its own config section, seeded from a frame of the steered pull.
  Embarrassingly parallel, restartable per window, and a failed window costs only itself.

**C7 – A new `umbrella` stage in `src/mcube/umbrella.py`, beside `neb`.** (user, 2026-10-22)
- One transparent, minimal config section, like every other stage: inputs written by the user,
  results written back by mcube. Built out from ASH's biased-sampling example.
- ASH's example (Biased-sampling page) is QM-only, a torsion restraint via
  `OpenMM_MD(restraints=[...])`, all windows in one loop, no CV log and no WHAM ("currently
  missing"; it points to FastMBAR). Read at ASH commit 2c7c056b.

**C10 – The bias is a native OpenMM force that mcube adds itself; PLUMED is not used.**
(2026-10-22; recommendation, revisit if a CV outgrows it)
- ASH keeps added forces. For QM/MM, `OpenMM_MDclass` reuses `theory.mm_theory.system` and
  `create_simulation` builds the `Simulation` on that same object; only barostat/Andersen forces
  are ever removed. `QMMMTheory.__init__` zeroes QM–QM bonded terms in forces *already present*
  (`modify_bonded_forces`), so the bias must be added after the QMMMTheory is built. It never
  touches CustomCentroidBondForce or CustomCVForce.
- Why not PLUMED: ASH's only working route is the openmm-plumed plugin (`OpenMM_MD_plumed`), a
  new dependency not known to be in the ASH env, with 1-based indices and nm/kJ units unlike
  every other mcube input; `plumed_object` on `OpenMM_MD` is non-functional. What PLUMED would
  buy (COM z, cylinder, a COLVAR file) is a dozen lines of CustomCentroidBondForce.
- Why not ASH's helpers: they share global parameter names (`k`, `r0`), `add_flatbottom_centerforce`
  never enables PBC, and `add_centerforce` compares r² with r0. A per-bond-parameter force with
  unique names avoids all three.
- Form (as implemented): a CustomCVForce over two CustomCentroidBondForce CVs, `zcv = z2 - z1` and
  `rcv = sqrt((x2-x1)^2+(y2-y1)^2)`, groups {reference, ion}, energy
  `0.5*mcube_us_k*(zcv-mcube_us_z0)^2 + 0.5*mcube_cyl_k*max(0, rcv-mcube_cyl_r)^2`, force group 31.
  The reference centre is weighted by **element masses given explicitly**, not the System's: those
  carry repartitioned 1.5 amu hydrogens, and the analysis computes the centre from the trajectory.
  The force uses raw (unwrapped) positions, so each window first moves the ion to the image
  nearest the channel; OpenMM never wraps integrated positions.
- CV time series: MM windows read `getCollectiveVariableValues` every `cv_frequency` steps (one
  force evaluation, no position copy). QM/MM windows use ASH's `specialatoms` trajectory (ion +
  reference only) and apply the minimum image afterwards. Writing whole-system DCDs at WHAM
  frequency would be ~70 GB per 5 ns MM window.
- One window per job also sidesteps an ASH trap: every `OpenMM_MD` call on one QMMMTheory adds
  another CustomExternalForce and stacks restraints, so a loop over windows in one process
  would leak forces.

**C11 – GPU.** (2026-10-22; partly verified from source)
- `OpenMMTheory(platform=...)` accepts any OpenMM platform (CUDA, HIP, OpenCL) and nothing in
  the QM/MM MD path is CPU-only. `OpenMM_MD(platform=)` is ignored for QM/MM: it must be set on
  the OpenMMTheory. Fixed: `qmmm.inputs.platform`/`properties` now reach `make_mm` (and every
  umbrella section has its own).
- From the GLP-1R session (2026-10-22): the ASH env's OpenMM has **no HIP platform**, so all MD
  runs on `gpu_cuda`. Measured ~185 ns/day for 215,790 atoms (POPC, 4 fs) on H100 and L40S, A100
  ~3x slower. For ~140 windows x 5 ns per ion that is ~4 GPU-days per ion.
- Per QM/MM step ASH copies positions off the device and updates the QM forces with a Python
  loop over every particle (`update_custom_external_force`), so for a ~150k-atom membrane system
  the GPU may not beat CPU for the MM part. To be measured (Q1). MM windows on GPU are
  unproblematic.

**C8 – Metals stay set in `preprocess`.** (user, 2026-10-22)
- The ion is listed under `metals`, so preprocess detects its coordination sphere and defines the
  QM region. Coordination restraints apply only while each window equilibrates; they are released
  for sampling and are never applied in QM/MM. This matches the rule "restrained for the initial
  MD equilibration, not for QM/MM".
- Consequence: the `md` stage cannot be reused unchanged for the MM windows, because it keeps
  the coordination restraints through production (`md.py:run`, restraints passed to
  `OpenMMTheory` for every stage). The umbrella stage holds them for its equilibration leg only.

**C9 – The QM region is rebuilt per window, with a fixed composition.** (2026-10-22; refined)
- The ion moves along the path, so one fixed QM atom list cannot serve every window. Each QM/MM
  window builds its region on its own starting frame: the ion, the `qm_waters` nearest waters
  (default: the detected shell) and the fixed `qm_residues`, cut at CA–CB by `build_qm_region`.
- Why fixed composition: the stitch compares MM and QM/MM energies at both ends of the QM/MM
  range. The QM energy carries a huge composition-dependent constant, which cancels between the
  two ends only if both have the same atoms in the QM region. So `qm_waters` is a number (La³⁺ 9,
  Lu³⁺ 8 from the EXAFS data in C4) and `qm_residues` is the union of the constriction residues,
  not "whatever coordinates now".
- The QM charge is the rounded sum of the MM charges the region replaces (raw sum logged), unless
  `charge` is given. A first-shell water leaving during a window is reported (`shell_changes`),
  not prevented.

**C13 – MM windows run on OpenMM directly; QM/MM windows through ASH.** (2026-10-22)
- ASH's pure-MM `OpenMM_MD` runs `simulation.step(all)` with no per-step hook, so it cannot log
  the CV. mcube already builds raw OpenMM Systems for the membrane path (`membrane.py`), so the
  MM windows and the pull use `membrane.membrane_system` (or the explicit equivalent) plus the
  same `ions.apply_metal_parameters`, LangevinMiddle, the configured platform.
- QM/MM windows use `make_mm` / `make_qmmm` / `OpenMM_MD` exactly as `mcube qmmm` does, in one
  `OpenMM_MD` call per window (equilibration is its discarded start), because a second call would
  add another external force (C10).

**C14 – The pull seeds every window; windows restrain their own frame's shell.** (2026-10-22)
- `pull` moves z0 from the ion's position to the nearer end of the range, then across it, at
  `pull_speed` (default 10 Å/ns), with no coordination restraints (the ion must be able to shed
  water), saving a seed when z0 passes each centre.
- Each MM window detects the shell in its own seed with `coordination_shell` (the same rule as
  preprocess), restrains it at the distances it has for `equilibrate`, then releases it (C8).
  QM/MM windows start from the finished MM window at the same z0.

## Findings

**F1 – `qmmm` did not apply `metal_parameters`.** (2026-10-22; fixed on branch `umbrella`)
- `md.py:apply_metal_parameters` swaps in IOD LJ terms, but `qmmm.py:make_mm` builds the MM
  system from the force-field XMLs alone. A QM metal therefore interacts with MM atoms through
  the stock Amber LJ in QM/MM and through IOD in MM. The umbrella stage must use one set in
  both; whether to fix `qmmm` too (it affects the CSR runs) is the user's call.
- Fixed (user, 2026-10-22): `ions.apply_metal_parameters` is the one implementation; md calls it,
  and `qmmm.apply_stage_metal_parameters` calls it in qmmm, scan, neb, ts and sp after the QM
  region is known and before the QMMMTheory is built. The set is the stage's
  `metal_parameters`, else `md.results.metal_parameters`, else md's input, else `iod`.
- **This changes the CSR runs** from now on: a QM Zn now meets MM atoms through IOD radii, as in md.

**F2 – OpenMM's stock La³⁺/Lu³⁺ are the 12-6-4 model's 12-6 terms without its C4.** (2026-10-22)
- `amber14/tip3p.xml` (OpenMM 8.4, built from AmberTools 17 `frcmod.ions234lm_126_tip3p`) gives
  La³⁺ sigma 0.31324 nm and Lu³⁺ 0.28954 nm, exactly the Rmin/2 of 1.758 and 1.625 Å in
  `frcmod.ions34lsm_1264_tip3p`; IOD would be 1.718 and 1.588 Å. This is the mislabelling of +3/+4
  rows the ASH session found in AmberTools 2016–2019.
- Consequence: `metal_parameters: none` runs half a 12-6-4 model for lanthanides. Use `iod` or
  `1264`, never `none`, for La/Lu.

## Planned checks, not yet done

- Q1: GPU. OpenMM platforms in the ASH env on Bunya (CUDA on `gpu_cuda`, HIP on MI355X); OpenMM
  MM part of QM/MM on GPU vs CPU; whether openmm-plumed is installed.
- Q2: resolved locally for OpenMM 8.4: `amber14/tip3p.xml` has CA, LA and LU templates
  (`tip3p_standard-Ca2+/La3+/Lu3+`). Confirm the Bunya ASH env has the same OpenMM.
- Q3: one QM/MM window in the membrane system: timing, energy drift, box stability (QM/MM in a
  membrane has never been run in mcube).
- Q4: pore radius profile of the AF2 model before planning windows.
- Q5: check the vendored 12-6-4 data against AmberTools on Bunya: `frcmod.ions34lsm_1264_tip3p`
  (+3/+4 rows), `lj_1264_pol.dat` (sha256 765e516f... for ParmEd's copy), and read the +2 rows of
  `frcmod.ions234lm_1264_tip3p` so Ca²⁺ can use 1264.
- Q6: the QM/MM window path (`engine: qmmm`), `energies --level xtb|dft` and `stitch` on real
  data have not run: they need ASH. Only the MM path and the analysis are tested locally.

**C12 – Implement 12-6-4 as `metal_parameters: 1264`.** (user, 2026-10-22; implemented on branch `umbrella`)
- No reason not to: the IOD 12-6 set under-binds 3+ ions, and the user asked.
- How Amber does it (ParmEd `add1264.py`, read): C4(ion, OW) is tabulated; for any other atom j,
  C4 = C4(ion) · α_j / 1.444 Å³, α_j from `lj_1264_pol.dat` by Amber atom type; ion–ion pairs
  sum. TIP3P C4 (kcal/mol/Å⁴): Ca²⁺ 87.3, La³⁺ 152, Gd³⁺ 198, Dy³⁺ 207, Y³⁺ 216, Lu³⁺ 249. The
  12-6 terms that go with it are a separate set (`frcmod.ions234lm_1264_tip3p`), not IOD.
- OpenMM form: a CustomNonbondedForce `-c4/r^4` over ion × all pairs, with NonbondedForce
  exceptions copied as exclusions and the same cutoff — what `AmberPrmtopFile` itself builds
  from `LENNARD_JONES_CCOEF`. mcube's force fields are XML, so it is built in `ions.py`.
- In QM/MM the ion is QM: switch C4 off for QM–MM pairs (xTB plus electrostatic embedding
  already polarises), keep 12-6. The MM end of the stitching BAR must use the same Hamiltonian.
- Caveat: α-ratio C4 is not fitted to protein ligands and tends to over-bind carboxylates;
  check against the Ln–acetate data in Jafari et al., *JPCB* 2024 (10.1021/acs.jpcb.3c06567).
  Jafari et al., *JCTC* 2025 (10.1021/acs.jctc.4c01535) used 12-6-4 for metal transporters.
- As built (`ions.py`): +3/+4 12-6 terms from ParmEd's copy of `frcmod.ions34lsm_1264_tip3p`, C4
  from ParmEd's `DEFAULT_C4_PARAMS['TIP3P']`, polarisabilities from ParmEd's `lj_1264_pol.dat`
  vendored unmodified in `src/mcube/data/`. Atom classes come from OpenMM's own template matching.
  Lipid17 types, absent from the table, take their GAFF analogue's value (cA->c3, oP->o, pA->p5 ...);
  that mapping is mcube's. Unknown types raise, unless given in `polarizabilities`.
- The C4 force copies **all** NonbondedForce exceptions as exclusions: the CPU and GPU platforms
  refuse nonbonded forces whose exclusions differ (the Reference platform does not check, which
  hid this until a CPU run).
- Under 1264 an untabulated metal raises rather than run on stock terms (F2).
