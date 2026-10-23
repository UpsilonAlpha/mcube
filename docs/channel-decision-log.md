# Ion-channel permeation PMF: decision log

A running record of the modelling decisions for the umbrella-sampling study of Ca²⁺ and two
lanthanide ions passing through the channel in `docs/channel/sequence.fasta`. Each entry gives
the decision, why, the evidence, and whether to revisit it. Same format as
[csr-decision-log.md](csr-decision-log.md). The protocol is in
[channel-us-plan.md](channel-us-plan.md).

**Aim.** Potentials of mean force for a light and a heavy lanthanide and for Ca²⁺ crossing the
channel, built with mcube's existing stages as closely as possible.

**The protein, in one page:** [channel/protein.md](channel/protein.md). It is BD122_14570, a
DUF6662 outer-membrane β-barrel of *B. diazoefficiens*, encoded next to lanmodulin in a
lanthanide methylotrophy cluster.

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
- Submitted 2026-10-22 from `/scratch/user/uquadhik/channel/af2` with `docs/channel/af2/af2.slurm`:
  **30410106 `af2_mono`** (`monomer_ptm`) and **30410107 `af2_trimer`** (`multimer`, one prediction
  per model, 5 models). Both on `gpu_cuda` with `--gres=gpu:a100:1`, because RCC's container
  wrapper notes it is incompatible with H100 (CUDA < 11.8); 8 CPUs, 96 GB, 12 h. Database paths are
  the module's own except `pdb_seqres`, which points at the real
  `.../databases/pdb_seqres/pdb_seqres.txt` (226 MB). Scratch was 241/300 GB at submission.
- Moved to AMD at the user's request ("use the MI355X where possible"), then back (2026-10-22).
  `alphafold/2.3.2-rocm` on MI355X (30410147, 30410148) failed in under a minute: JAX found
  "No visible GPU devices". RCC's container is ROCm 6.1.1, which predates MI355X (gfx950)
  support, so AF2 cannot use those nodes as installed. `sbatch --test-only` estimates then gave
  A100 25 Oct, MI300X 28 Oct, MI210 29 Oct (independent of a 4 h or 12 h limit; the per-user
  4-GPU-job QOS is shared with the GLP-1R session's queued MD). Resubmitted on A100:
  **30410203 `af2_mono` (4 h)** and **30410204 `af2_trimer` (8 h)**. The script defaults to
  `gpu_cuda` A100 again; AMD is chosen on the sbatch line (MI300X or MI210, not MI355X).

- **Monomer result (30410203, 22 min on A100, `channel/af2/monomer/mature/`).** All five models
  agree (pLDDT 95.6–96.7); ranked_0 = model_3, mean pLDDT 96.7 (255/264 residues ≥ 90), pTM 0.91.
  `docs/channel/af2/pore.py` on ranked_0:
  - a **closed 12-stranded β-barrel** (strands 133–143, 149–164, 167–181, 189–206, 214–226,
    233–250, 254–269, 275–292, 295–308, 315–329, 334–348, 363–372; first and last strand Cα
    4.35 Å apart), Cα radius 12.2 Å. This supports C1 (an OMP β-barrel).
  - residues 110–132 (the N-terminal 23, coil, pLDDT ~87) **thread into the lumen** and narrow it.
    Pore radius (heavy atoms, vdW radii) is 0.8–1.1 Å over 9 Å of the axis (z −14 to −5 Å), lined
    by Phe119, Met122, Ala286 and Glu300, and 2.3–3.9 Å over the rest of the membrane span.
  - So in this model the channel is **not open to a hydrated ion**, and only marginally to a bare
    one (ionic radii: Ca²⁺ 1.00–1.12 Å, La³⁺ 1.03–1.22 Å, Lu³⁺ 0.86–0.98 Å). A permeation path
    needs the N-terminal segment to move, or the trimer to differ. Pending: the trimer (30410204),
    and whether to fold without 110–132 or pull it out of the lumen in MD.
- **Trimer result (30410204, 38 min, `channel/af2/trimer/mature_trimer/`): AF2 does not support
  a trimer.** 0.8·ipTM + 0.2·pTM is 0.17–0.21 for all five models; ranked_0 (model_3) has ipTM
  0.16, pTM 0.43. Each chain folds the same barrel as the monomer (per-chain pLDDT 89.6–90.0,
  identical DSSP, 12 strands, closed, segment 110–132 again in the lumen), but the chains barely
  touch: 7 Cα pairs within 8 Å for A–B, 1 for B–C, none for A–C. Chain A's narrowest point is
  0.57 Å at z −5 Å (lining Thr242–Ala244, Asn259, Leu260). **Use the monomer** for the system;
  the FapF match (a trimer) alone is not enough to override AF2 here.

**C1/C2 follow-up – relatives and conservation of the lumen segment.** (2026-10-23)
- Sent to NCBI at the user's request: blastp of all 373 aa against ClusteredNR (RID DNGYGYKX016).
  Top hit `APO51500.1`, hypothetical protein BD122_14570 of *Bradyrhizobium diazoefficiens*,
  373 aa, 98.9 % identity over 100 % — this protein. Every other top hit is a 301-aa
  **DUF6662 family protein** of *Bradyrhizobium* (83–99 %), aligned to query 73–373: the RefSeq
  proteins start at **Met73**, not Met85 as C1 guessed. 1–72 are not part of the protein;
  73–109 is still the signal peptide (`MANAVEAVVVFGMRARAQ` n-region, `TFWIMAACTFALA` h-region,
  `APTCVHA|EE`), so the mature chain 110–373 stands.
- AF2's own searches (no outside service): UniRef90 hits are 352 mostly "uncharacterized"
  Alphaproteobacteria proteins (*Bradyrhizobium* 69, *Methylobacterium* 56, *Bosea* 24,
  *Rhodopseudomonas* 19); a few are annotated "MetA-pathway of phenol degradation". PDB70
  templates: **5O65 FapF** (*Pseudomonas* functional-amyloid secretion channel, a trimer of
  β-barrels, header keywords "amyloid, secretion, plug, barrel"; probability 98.2, query 121–350)
  and **4RL8** ("protein involved in meta-pathway of phenol degradation", *P. putida*, OM β-barrel;
  96.3); then generic OMPs below 90.
- `docs/channel/af2/conservation.py` over AF2's MSAs (2173 sequences, 1792 effective; identity to
  the query is ~10 % on average, so this is a diverse family):

  | Region | Coverage | Identity | Conservation (1 − H/ln 20) |
  |---|---|---|---|
  | 110–132 segment | 0.35 | 0.09 | 0.32 |
  | strands | 0.72 | 0.12 | 0.19 |
  | loops | 0.48 | 0.10 | 0.21 |

  Coverage of the segment climbs from 1 % at Glu110 to 58 % at Glu132: most distant relatives
  have no aligned N-terminal segment (a shorter N-terminus, or jackhmmer's local alignment
  dropping it). Where present, **Phe119** (lumen-facing, at the constriction) is the most
  conserved residue of the segment (0.63; identity 0.19), then Gly131 (0.52; 0.39). The acidic
  residues are not: Asp114 0.40 and Glu116 0.42 at ~15–20 % coverage, Asp126 0.16, Glu132 0.17.
  Entropy over few sequences reads high, so the segment numbers are an upper bound.
- Reading: the segment's position is family-specific but not universal; its constriction
  residue Phe119 is conserved where the segment exists. Together with the FapF match (a plugged,
  trimeric secretion barrel), a gate or plug is plausible, not shown. The acidic residues lining
  the segment are not a conserved cation site.

**C1 follow-up – the gene sits in a lanthanide methylotrophy cluster, next to lanmodulin.** (2026-10-23)
- User priority (2026-10-23): decide whether this is a lanthanide channel; the gate is a side
  question, and the N-terminal segment is cropped unless it takes part in ion transport.
- Genome context (NCBI E-utilities, public records only): BD122_14570 is CP013127.1:3185355..3186476
  (+) in *B. diazoefficiens* USDA 122. Neighbours, RefSeq names:
  - BD122_14565 (+, 3184901..3185293, **61 bp upstream on the same strand**, so probably the same
    transcript): **lanmodulin LanM** (WP_028174364.1). 130 aa, signal peptide, four EF-hand loops
    `DTDNDG`, `DPDHDG`, `DPDRDG`, `NPDKDG`: three with the Pro at loop position 2 that marks LanM.
  - BD122_14530: **lanthanide-dependent methanol dehydrogenase XoxF5** (WP_011088953.1), with a
    cytochrome c (14525), the glutathione formaldehyde pathway (14515–14520) and a further PQQ
    alcohol dehydrogenase (14560, methanol/ethanol family).
  - **lanthanide ABC transporter PedA1 / PedB (LutE) / PedC** (BD122_14640–14655, 10 kb
    downstream: the inner-membrane Ln importer), other ABC transporter genes (14550, 14590–14600)
    and PQQ-catabolism genes (14610–14625). Full table: `docs/channel/identity/`.
- Barrel lumen (side chains facing the axis in the AF2 model, segment excluded): 21 Asp/Glu
  against 11 Arg/Lys, so the lumen is net negative. **Glu300 is the most conserved lumen residue
  in the whole barrel** (identity 0.59 against 0.12 for the strands, coverage 0.79), and it
  lines the constriction with Phe119 and Met122 of the segment and Ala286.
- Reading: a periplasmic Ln-binding protein encoded next to it, an XoxF Ln-dependent methanol
  dehydrogenase nearby, and a conserved carboxylate at the constriction all point to a role in
  lanthanide uptake across the outer membrane. This is genomic and structural inference; how the
  ion crosses (free Ln³⁺, or bound to a chelator) is not settled by it.
- For the segment: the constriction is formed by conserved Glu300 *with* the segment's Phe119,
  so cropping may remove half of the site. Check before cropping (below).

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
- All three ions can use `metal_parameters: 1264` (C12): the +2 rows were added on 2026-10-22
  from the AmberTools in the ASH env, so Ca²⁺ is no longer limited to IOD.
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

**C15 – Desolvation is validated in bulk first, and watched in every window.** (user, 2026-10-22:
"especially the key desolvation step with reference to the literature")
- Bulk targets (experiment): Ca²⁺ Ca–O 2.46 Å, CN 8 (Jalilehvand et al., *JACS* 2001, 123, 431);
  La³⁺ 2.52–2.56 Å, CN ≈ 9; Lu³⁺ 2.33–2.34 Å, CN ≈ 8 (Persson et al., *Chem. Eur. J.* 2008,
  10.1002/chem.200701281; compiled in Alanís-Manzano et al., *ACS Omega* 2024,
  10.1021/acsomega.4c04947). Water exchange: heavy Ln³⁺ 4e7–5e8 s⁻¹ (Cossy, Helm & Merbach,
  *Inorg. Chem.* 1988, 27, 1973; Helm & Merbach, *Chem. Rev.* 2005, 10.1021/cr030726o); no direct
  Lu³⁺ value; Ca²⁺ fast and uncertain (ps to tens of ps; Koskamp et al. 2019, 10.1021/acs.jpcc.9b06403).
- What the models are known to do (Li/Merz papers): for Ca²⁺ in TIP3P, IOD gives 2.45 Å, CN 8.0 but a
  hydration free energy ~17 kcal/mol too weak; 12-6-4 gives 2.46 Å, CN 8.0 and the right HFE. For
  M³⁺ the IOD set under-binds by ~80 kcal/mol on average; 12-6-4 fixes HFE and distance but gives
  Ln³⁺ CN ≈ 9–10, so **Lu³⁺ under 12-6-4 is expected to over-hydrate** (target 8), which would
  inflate its desolvation cost. Which set the PMF uses is decided by the bulk results, not assumed.
- Bulk protocol (`docs/channel/bulk/`): per ion and set, MM NPT in TIP3P with the RDF first peak
  and minimum, CN from the minimum, and first-shell exchange counted with the Impey t* = 2 ps
  rule; then xTB/MM MD of the same box through mcube's preprocess + qmmm, analysed the same way.
  Pass: peak within ±0.05 Å of experiment, CN within ±0.5. Exchange within ~10× of experiment;
  the first short runs give structure only, and a long run (≥ 100 ns, larger box) is needed for
  kinetics.
- In the channel: every window records CN_water and CN_protein-O around the ion. A 1D z PMF is
  kept unless constriction windows show bimodal CN histograms, CN that depends on the seeding
  direction, or fewer than ~5 exchanges per window; then the constriction goes to a 2D (z, CN)
  treatment or replica exchange. QM/MM windows (ps) cannot exchange water at all, so they are
  seeded from MM windows and stitched per CN state; they correct the chemistry of a given
  hydration state and do not sample desolvation.
- Bulk runs 2026-10-22, `/scratch/user/uquadhik/channel/bulk`, general partition, 8 CPUs, 4 ns MM +
  10 ps xTB/MM each. Smoke test **30410350** (La³⁺ IOD, 10 ps MM + 0.2 ps xTB/MM) completed RC 0:
  the whole chain ran (MM; mcube preprocess with skip_solvate building ion + shell as the QM
  region; mcube qmmm QM/MM MD through `apply_stage_metal_parameters`, the first real QM/MM run of
  that code) at **0.11 s per xTB/MM step** for a ~28-atom QM region in a 4.1k-atom box, i.e.
  ~18 min per 10 ps. Even 10 ps of MM gave La³⁺ IOD a 2.51 Å peak and CN 8.95 (target 2.52–2.56, ~9).
  It also exposed two analysis faults, both fixed before the real runs: the RDF "minimum" was the
  first zero bin at the peak's tail (2.69 Å for La³⁺), now the middle of the flat gap (e425c0b);
  and an analysis with nothing left after `--skip` divided by zero (51d2965). The first six jobs
  (30410363–68) started 38 s before the fix reached Bunya, so they were cancelled after a minute
  and resubmitted on 51d2965: **30410512** Ca IOD, **30410513** Ca 1264, **30410514** La IOD,
  **30410515** La 1264, **30410516** Lu IOD, **30410517** Lu 1264.
  A 4 ns MM run measures structure (peak, CN); it is too short for La/Lu exchange kinetics.
- **Results, 2026-10-23** (TIP3P box of 1376 waters, 300 K NPT, 4 ns, samples every 0.5 ps; CN to
  the RDF minimum; exchange with t* = 2 ps). Ca IOD, Ca 1264 and La IOD hit the 6 h limit in MM
  (6-h CPU allocations ran ~15–20 ns/day for this 4.1k-atom box), so they have no result yet.

  | ion / set | peak Å (target) | CN (target) | first-shell residence | verdict |
  |---|---|---|---|---|
  | La³⁺ 1264, MM | 2.51 (2.52–2.56) | 9.65: 9 ×36 %, 10 ×64 % (~9.1) | 157 ps, 57 exchanges/ns | distance passes; CN +0.5, at the limit |
  | Lu³⁺ IOD, MM | 2.33 (2.33–2.34) | 8.88 (8) | 145 ps | distance passes; CN +0.9 fails |
  | Lu³⁺ 1264, MM | 2.33 (2.33–2.34) | 7.62: 7 ×38 %, 8 ×62 % (8) | 71 ps | both pass |
  | La³⁺ 1264, xTB/MM 4.3 ps | 2.43 | 8.83 | none in the run | Ln–O 0.1 Å short |
  | Lu³⁺ IOD box, xTB/MM 8 ps | 2.33 | 7.88 | none | both pass |
  | Lu³⁺ 1264 box, xTB/MM 5.1 ps | 2.25 | 7.03 | none | Ln–O 0.08 Å short, CN −1 |

  Three findings.
  1. **Lu³⁺: in our runs 12-6-4 gives CN 7.6 and IOD 8.9 (experiment ~8).** Withdrawn as "12-6-4
     is the better model" (user, 2026-10-23: a claim against the literature needs strong
     justification). What is established: the implementation is Amber's (F5); one 4 ns run per set
     with ~400 exchanges, so the CN means are converged for this protocol. Not established: how
     this compares to Li/Song/Merz's own Lu³⁺ 12-6-4 value (C15 quoted only their Ln³⁺ average,
     CN ≈ 9–10). Until that number is in hand, no model is preferred for Lu³⁺. La³⁺ 12-6-4 gives
     9.65 (~9.1 expected). The choice per ion also waits for La IOD and both Ca runs.
     **Published values found (2026-10-23)**, Li, Song & Merz 2015 SI Table SI.4 (TIP3P, 12-6-4):
     La³⁺ 2.53 Å / CN 9.7 and Lu³⁺ 2.34 Å / CN 9.0 (read from PMC4306492's SI by an agent; rows
     assigned by order and checked against the HFE column). So our La³⁺ 12-6-4 (2.51 / 9.65)
     reproduces the paper and our Lu³⁺ IOD (2.33 / 8.88) agrees with its 12-6 scan (~2.34 / ~9),
     but **our Lu³⁺ 12-6-4 CN (7.62) does not match the published 9.0**, at the same distance.
     Parameters (C4 249, Rmin/2 1.625, ε 0.09235154) are the paper's, and the energies are Amber's
     (F5), so the difference must be in protocol: theirs was 722 waters, 1 fs, Langevin γ = 5 ps⁻¹,
     500 ps NVT + 500 ps NPT, then 2 ns NPT, 10 Å cutoff, Amber 12. Replication with that protocol
     (MM only), La³⁺ as the control: **30416377** Lu³⁺, **30416378** La³⁺ (`<ION>_1264_lsm/`).
     Both died at once ("No module named 'openmm'"): with `--export=ALL` from a shell where ASH
     was already active, `conda activate ASH` did nothing and bare `python` was miniforge's base.
     The job scripts now call `$CONDA_PREFIX/bin/python` and check it imports OpenMM and ASH
     (e6eb2cf). Resubmitted: **30417999** Lu³⁺, **30418000** La³⁺.
     **Result (2 ns each, 1 fs, friction 5 ps⁻¹, 719 waters):** Lu³⁺ 2.33 Å, **CN 7.61** (7: 39 %,
     8: 61 %), residence 80 ps; La³⁺ 2.51 Å, CN 9.61, residence 167 ps. The protocol changes
     nothing: Lu stays at 7.6 (7.62 before), and La reproduces the published 9.7. With the
     implementation also identical to Amber's (F5), the Lu gap to the published CN 9.0 is not
     ours to close by protocol. Either the SI value was misread (an agent read it; check the SI
     table by hand), or the SI used a different Lu C4/parameter set from AmberTools'
     `frcmod.ions34lsm_1264_tip3p` (C4 249). Against experiment (2.33 Å, CN 8) Lu 1264 is
     0.4 under on CN and Lu IOD (CN 8.9) 0.9 over; neither is exact.
     Checked against the paper itself (Li, Song & Merz, J. Phys. Chem. B 2015, 119, 883; PMC4306492):
     its Table 6 gives exactly our Lu³⁺ parameters (Rmin/2 1.625 Å, ε 0.09235154, C4 249), and the
     text says the Ln³⁺ CNs came out "in the range of 9 to 10" across the series. So CN ≈ 9 for Lu
     is what they report, with the same parameters, and we get 7.6. The SI PDF itself would not
     download (PMC interstitial), so the SI.4 row is still unchecked by hand. One untested
     explanation: Lu³⁺ sits near the CN 8/9 switch, so short runs started from a 9-coordinate
     shell can stay there; ours exchange 90 times per ns and are equilibrated.
     **SI read by hand (Europe PMC supplementary archive, jp505875v_si_003.pdf):** Table SI.4,
     TIP3P, Lu³⁺ HFE −839.3, IOD 2.34 Å, **CN 9.0**; La³⁺ 2.53 Å / 9.7. Every heavier Ln³⁺
     (Nd–Lu) is 9.0. Their protocol (main text): ~722 waters, 1 fs, Langevin 5 ps⁻¹, 500 ps NPT,
     2 ns sampling, 10 Å cutoff, PME, SHAKE, AMBER 12, CN by integrating g(r) to its first
     minimum, i.e. what we replicated. Their own grid (Table SI.2a, TIP3P) puts Rmin/2 = 1.6 Å
     with C4 = 500 at CN 8.4 and without C4 at 9.0, so Lu (1.625 Å, C4 249) should fall
     between 8.4 and 9; **our 7.6 is below that bracket, as if our C4 acted larger than 249.**
     The working assumption (user, 2026-10-23) is a fault on our side, not in the paper.
     F5 cannot rule that out: it evaluated mcube *and* the Amber-built prmtop in OpenMM, so a
     fault in how OpenMM evaluates the C4 term would pass it.
     Decisive test, `docs/channel/bulk/amber_crosscheck.py` (no mcube code; tleap + ParmEd
     prmtop of the Lu box; paper protocol): **30430921** OpenMM on the prmtop, **30430922**
     AmberTools sander with `lj1264 = 1` (confirmed on in its output; serial, 5.3 ns/day,
     100 ps chunks). sander ≈ 9 and OpenMM 7.6 → the fault is in OpenMM's/our C4 evaluation;
     both 7.6 → not our code, and the next step is asking the authors, not alleging an error.
  2. **The MM first shells exchange water about 10–100× faster than experiment** (known; it is why
     the hydration shell is restrained while the MM system is prepared for QM/MM, user): residence
     70–160 ps against ~2–25 ns for heavy Ln³⁺ (Cossy, Helm & Merbach 1988; Helm & Merbach 2005).
     So in MM umbrella windows of several ns, hydration is sampled rather than frozen (the C15 worry
     about frozen windows does not apply to MM). But a shell that is this labile is also probably
     easier to strip, so the MM desolvation barrier is suspect, and correcting it is the reason for
     the QM/MM leg. Unsettled: whether the speed is the model's or the t* = 2 ps rule's (Laage &
     Hynes); a stable-states count would decide.
  3. **xTB/MM of Ln³⁺ aqua ions loses SCC convergence after a few ps.** In La 1264 (step 6328) and
     Lu 1264 (step 7142) xTB stopped with "SCC iterator did not converge": the HOMO–LUMO gap had
     fallen to 0.0006 and 0.0009 Eh, with fractional occupations of the frontier orbitals. An
     Ln³⁺–water cluster should have a large gap. The likeliest cause is the known boundary problem
     of electrostatic embedding: a +3 QM region with low-lying Ln 5d orbitals next to bare positive
     MM charges (the H-bonded second-shell water hydrogens). Every QM/MM umbrella window would meet
     the same failure. Remedies in order: raise xTB's electronic temperature, soften the MM point
     charges (ASH `hardness_PC`), and put the second shell in the QM region (Rode-style QMCF
     practice). The xTB/MM distances above come from short runs and are provisional.
- Resubmitted 2026-10-23 (12 h MM limit, 8c9536d): MM **30416201** Ca IOD, **30416202** Ca 1264,
  **30416203** La IOD. QM/MM-only reruns on the finished boxes, to separate the two SCC remedies:
  xTB electronic temperature 1000 K, first shell QM, `qm_t1000/`: **30416204** La 1264,
  **30416205** Lu IOD, **30416206** Lu 1264; second shell in QM (cutoff 5.2 Å La, 5.0 Å Lu),
  default 300 K, `qm_s2/`: **30416207** La 1264, **30416208** Lu IOD, **30416209** Lu 1264.
- **La³⁺ IOD (30416203), 4 ns MM:** La–O 2.49 Å, CN 9.06 (9: 94 %, 10: 6 %), residence 171 ps,
  42 exchanges/ns. Against La 1264 (2.51 Å, CN 9.65) and the target ~2.54 Å / CN 9: IOD gets
  the CN right and Ln–O 0.05 Å short; 1264 is closer on distance and over-coordinates. Its
  default xTB/MM leg (first shell QM, 300 K, bare charges) lost SCC at step 6362 like the others
  (2.45 Å, CN 9.0 over 4.4 ps) — no new information for F7.

**C16 – DFT/MM for the desolvation step; xTB/MM for the rest of the QM/MM leg.** (user, 2026-10-23;
level and details proposed, not yet decided)
- Why: xTB is an approximation that, for Ln³⁺ aqua ions, shortens Ln–O by up to 0.1 Å and loses
  SCC convergence after a few ps (F4). Desolvation is the step the PMF turns on, so it gets a
  proper DFT; the parts of the path where the hydration shell is intact stay at xTB/MM.
- Where: the constriction windows whose MM hydration profile (`wham` → `<section>_hydration.csv`)
  shows the water count falling, i.e. where waters are being exchanged for protein oxygens.
- Proposed level: r2SCAN (a fast meta-GGA, no exact exchange, so density fitting makes it cheap)
  with def2-SVP and def2 ECPs, D3(BJ) dispersion if ASH carries its gradient. Available in the ASH
  env: PySCF 2.13 with libxc r2SCAN, simple-dftd3. Missing: `basis_set_exchange` (PySCF's bundled
  def2 stops at La, so Lu needs it), dftd4, gpu4pyscf.
- Mechanics: `qm_engine: dft` with the `dft_*` keys in the `umbrella_qmmm` section already runs
  QM/MM MD at DFT/MM; the stitch and reweighting are unchanged.
- First: time one DFT/MM gradient on the bulk La box (ion + first shell) at 16 and 32 cores, then
  a few ps of bulk DFT/MM MD for La³⁺ and Lu³⁺ against EXAFS, which also gives xTB a reference.
- Measured 2026-10-23 (`docs/channel/dft/bench_dft.py`, bulk La³⁺ 12-6-4 box, ion + 9 waters = 28
  QM atoms, r2SCAN/def2-SVP, density fitting, def2 ECP, electrostatic embedding, no dispersion):
  16 cores 38.7 / 14.2 / 12.7 s per gradient (first from scratch), 32 cores 14.1 / 9.5 / 10.2 s.
  That is ~7–9 ps/day per window at 1 fs, against 0.2–0.3 s for one xTB/MM gradient (~40×).
  So a DFT/MM window gets a few ps, enough to reweight or check the shell, not to sample exchange.
- `basis_set_exchange` 0.12 installed into the ASH env (dry run: adds it, regex and Unidecode,
  upgrades nothing; numpy 1.26.4, OpenMM 8.4, PySCF 2.13 unchanged). PySCF still does not fall
  back to BSE for Lu, so the basis comes from a file (`docs/channel/dft/bse_basis.py` writes
  def2-SVP for H,C,N,O,S,Ca,La,Lu; pass it as `dft_basis_file` and `dft_ecp_file`).
- **Open: the def2 ECPs are not the same kind for the two ions.** La has a 46-electron core
  (4f outside, empty), Lu a 28-electron core (4f¹⁴ in the valence). A La/Lu comparison at DFT
  compares two different pseudopotential treatments; options are to accept it (it is the def2
  standard), or put Lu on a 4f-in-core ECP (MWB-type, 4f in core for both) for the comparison.

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

**F5 – mcube's 12-6-4 reproduces Amber's exactly.** (2026-10-23, `docs/channel/bulk/check_1264.py`)
- Reference: the Lu³⁺ 12-6-4 bulk box built with tleap (frcmod.ions234lm_1264_tip3p) and ParmEd
  `add12_6_4 :LU watermodel TIP3P` (AmberTools in the ASH env), read by OpenMM's AmberPrmtopFile.
  mcube: amber14/tip3p.xml + `ions.apply_metal_parameters(..., "1264")`. Same frame, PME, 1 nm.
- The energy change for moving the ion by −0.2, +0.2, +0.4 Å (18.3935, 15.0358, 92.7867
  kcal/mol) and the force on the ion in all four positions agree to every printed digit. The
  totals differ by a constant 0.15 kcal/mol out of ~14,000, unchanged when the ion moves, so it
  does not involve the ion (most likely how the two apply the long-range LJ correction to
  water–water pairs).

**F2 – OpenMM's stock La³⁺/Lu³⁺ are the 12-6-4 model's 12-6 terms without its C4.** (2026-10-22)
- `amber14/tip3p.xml` (OpenMM 8.4, built from AmberTools 17 `frcmod.ions234lm_126_tip3p`) gives
  La³⁺ sigma 0.31324 nm and Lu³⁺ 0.28954 nm, exactly the Rmin/2 of 1.758 and 1.625 Å in
  `frcmod.ions34lsm_1264_tip3p`; IOD would be 1.718 and 1.588 Å. This is the mislabelling of +3/+4
  rows the ASH session found in AmberTools 2016–2019.
- Consequence: `metal_parameters: none` runs half a 12-6-4 model for lanthanides. Use `iod` or
  `1264`, never `none`, for La/Lu.

**F6 – `make_dft` gave PySCF no ECP.** (2026-10-23, fixed in f4fd915)
- `qmmm.make_dft` passed `dft_basis` but never an ECP, so a def2 basis on La or Lu would have run
  all-electron with a valence-only basis. It now defaults the ECP to the basis name for def2
  bases (`dft_ecp` overrides it), and reads `dft_basis_file` / `dft_ecp_file` (NWChem format)
  for elements PySCF lacks. Only DFT runs with heavy elements are affected; CSR runs are xTB.

**F7 – Neither SCC remedy makes xTB/MM of Ln³⁺ stable; the second shell in QM makes the shell wrong.** (2026-10-23)
- The 10 ps xTB/MM reruns on the finished MM boxes (jobs 30416204–09; shell from the RDF minimum
  of the QM/MM trajectory, first 2 ps skipped):

  | Box | Remedy | Result | Ln–O (Å) | CN |
  |---|---|---|---|---|
  | La 1264 | etemp 1000 K | SCC failed at step 7092, gap 0.0042 Eh | 2.39 | 9.0 |
  | Lu IOD | etemp 1000 K | SCC failed at step 8984, gap 0.00005 Eh | 2.23 | 8.0 |
  | Lu 1264 | etemp 1000 K | completed 10 ps | 2.25 | 7.09 |
  | La 1264 | 2nd shell QM | SCC failed at step 6404, gap 0.0035 Eh | 2.63 | 11.3 |
  | Lu IOD | 2nd shell QM | completed 10 ps | 2.51 | 11.3 |
  | Lu 1264 | 2nd shell QM | SCC failed at step 6594, gap 0.0007 Eh | 2.43 | 9.70 |

  Targets (C15): La ~2.54 Å / CN 9, Lu 2.33 Å / CN 8.
- The gap still collapses with the second shell in the QM region, where no bare MM charge is
  near the ion. So the boundary is not the main cause; GFN2-xTB's own description of Ln³⁺
  (4f in core) puts empty states near the occupied ones. The failures came 1–3 ps later than
  in F4, not never.
- With more water in the QM region xTB pulls 11–12 waters inside the first-shell minimum, and
  the Ln–O peak moves 0.1–0.2 Å *out*. With the first shell only and 1000 K, Ln–O is 0.08–0.15 Å
  too *short* and CN is right for La and Lu IOD. Neither is a usable reference structure.
- Consequence: xTB/MM is not reliable for any window where an Ln³⁺ is in the QM region, not only
  the desolvation step (C16). Untested remedy: softer MM point charges (`hardness_PC`).
- g-xTB considered and set aside (user, 2026-10-23): ASH's `gxTBTheory` runs the standalone gxtb
  binary and drops the MM point charges, so it cannot do electrostatic embedding, and the ASH
  env's xtb 6.7.1 has no `--gxtb` flag for the xTBTheory route.
- Submitted 2026-10-23, first shell QM, `xtb_pc_hardness: elements` (xTB's per-element hardness,
  i.e. smeared MM charges; 0d4826d), 10 ps: at 300 K, `qm_hel/`: **30420814** La 1264,
  **30420815** Lu IOD, **30420816** Lu 1264; with 1000 K etemp too, `qm_ht/`: **30420817** La
  1264, **30420818** Lu IOD, **30420819** Lu 1264.
- **Result: `hardness_PC: elements` fails too.** All six stopped with "SCC iterator did not
  converge" ("xTB subprocess gave error", RC=1):

  | Box | Remedy | Failed at step | Last gap (Eh) | Ln–O (Å) | CN |
  |---|---|---|---|---|---|
  | La 1264 | soft PC | 6486 | 0.0139 | 2.49 | 9.02 |
  | Lu IOD | soft PC | 8738 | 0.0017 | 2.33 | 7.87 |
  | Lu 1264 | soft PC | 5921 | 0.0235 | 2.31 | 7.07 |
  | La 1264 | soft PC + 1000 K | 9063 | 0.0018 | 2.41 | 9.04 |
  | Lu IOD | soft PC + 1000 K | 9067 | 0.0041 | 2.33 | 8.00 |
  | Lu 1264 | soft PC + 1000 K | 6902 | 0.0074 | 2.27 | 7.00 |

  Softened MM charges give the best shells xTB has produced (La 2.49 Å / CN 9, Lu IOD 2.33 Å /
  CN 7.9, against 2.54/9 and 2.33/8), so the bare charges were what shortened Ln–O. But every
  remedy and combination still loses SCC within 6–9 ps. **GFN2-xTB/MM is out for any window with
  Ln³⁺ in the QM region**; the fallback is DFT/MM there (~7–9 ps/day per window, C16). Ca²⁺ is untested.

**F8 – The g-xTB build of xtb silently ignores point charges.** (2026-10-23, `docs/channel/gxtb/`)
- g-xTB v2.0.1 (grimme-lab release, `xtb-6.7.1-gxtb-140526-linux-x86_64.tar.xz`, sha256 checked)
  is in `/scratch/user/uquadhik/opt/gxtb/xtb-6.7.1/bin/xtb`. Its README lists `$pcem` as
  unsupported, and no ASH branch adds anything (`NEW` and `master` are the same; `dev` is from 2021).
- Given the same ORCA-style `pcharge` + `$embedding` input ASH writes, `xtb --gxtb` runs, prints
  no warning, and returns the gas-phase energy to every digit (La³⁺ S1 −718.735593 Eh both ways).
  **So ASH's `xTBTheory(xtbmethod="gxtb")` in QM/MM would quietly be mechanical embedding.**
  Never use it with electrostatic embedding until g-xTB supports `$pcem`.
- g-xTB counts no valence electrons on the Ln³⁺ (72 for La³⁺·9H₂O, i.e. 8 per water), like GFN2;
  in the smoke frames its gaps are ~9–11 eV against GFN2's 3.8–5.2 eV (a range-separated-hybrid
  target, so gaps are not comparable in absolute terms; what matters is whether they collapse).
- The cluster test (`cluster_test.py`, frames from the failing GFN2 xTB/MM runs LA_1264/qmmm and
  LU_1264/qmmm, last frame included; S1, S1 in the box's TIP3P charges for GFN2 and r2SCAN, S2;
  gas-phase S1 optimisations) against r2SCAN/def2-SVP: **30421227** (`channel/gxtb/gx_test/`).

## Planned checks, not yet done

- Q1: answered 2026-10-22 (`docs/channel/probe_env.py` on bun147): the ASH env has Python 3.11.14,
  **OpenMM 8.4 with Reference, CPU and CUDA** (no HIP, no OpenCL), pyscf 2.13.0, mdtraj 1.11.1;
  pymbar, openmm-plumed and basis_set_exchange are absent. Jobs import mcube from the scratch
  checkout (`PYTHONPATH=$S/mcube/src`), and `mcube.umbrella` and `mcube.ions` import there.
  Still open: GPU vs CPU speed for the MM side of QM/MM (needs a real window).
- Q2: answered: the ASH env's OpenMM 8.4 `amber14/tip3p.xml` has the CA, LA and LU templates, with
  the same stock radii as the local copy (F2 holds on Bunya: La³⁺ Rmin/2 1.758, Lu³⁺ 1.625 Å).
- Q3: one QM/MM window in the membrane system: timing, energy drift, box stability (QM/MM in a
  membrane has never been run in mcube).
- Q4: pore radius profile of the AF2 model before planning windows.
- Q5: answered (`docs/channel/probe_amber.sh`). AmberTools lives inside the ASH conda env
  (`$CONDA_PREFIX/dat/leap/parm`; there is no amber module). Its `frcmod.ions234lm_1264_tip3p`
  holds all +2/+3/+4 rows; every +3/+4 row equals what was vendored, and the +2 rows (Ca²⁺ 1.642 Å,
  0.10185975 kcal/mol) were added. Its `lj_1264_pol.dat` (sha256 9a64b54f...) is **newer than the
  ParmEd test copy** first vendored (765e516f...): it lists the Lipid17 types itself, with exactly
  the values the GAFF-analogue mapping had given, so that mapping was removed; it adds ion
  polarisabilities (La³⁺ 1.141, Lu³⁺ 0.629 Å³) and sets GAFF `ho` to 0.000 (ParmEd's copy: 0.387).
  The AmberTools file is now the vendored one (5ef051c).
- Q7: a QM region holding a titratable residue next to a metal-bound water can move a proton by
  itself (the CSR log's D96: Lys226 protonated a Zn-bound hydroxide during a QM/MM optimisation).
  For Ln³⁺ near Asp/Glu, hydrolysis of a bound water in a QM/MM window would change the model
  between windows and look like noise. Check the H count per QM oxygen at the start and end of
  each QM/MM window before trusting it.
- Q6: the QM/MM window path (`engine: qmmm`), `energies --level xtb|dft` and `stitch` on real
  data have not run: they need ASH. Only the MM path and the analysis are tested locally.
  Partly answered 2026-10-22: the full suite in the ASH env on bun147 (2b90dfc) gave **168 passed,
  0 skipped** in 82 s, including the four ASH-only test files and the end-to-end MM test. No test
  yet drives a QM/MM window; that still needs one real window.
- F3 (fixed, 2b90dfc): the umbrella stage left OpenMM's CPU thread count at its default, one per
  core of the whole node. On a 2-core allocation the end-to-end test ran ~2050 threads and made
  no visible progress for minutes. The CPU platform now gets the stage's `cores` (else
  `$OMP_NUM_THREADS`, else `$SLURM_CPUS_PER_TASK`), as every other mcube stage does.

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
- As built (`ions.py`): +2/+3/+4 12-6 terms from AmberTools' `frcmod.ions234lm_1264_tip3p` (the +3/+4
  rows first came from ParmEd's copy and were checked identical), C4 from ParmEd's
  `DEFAULT_C4_PARAMS['TIP3P']`, polarisabilities from AmberTools' `lj_1264_pol.dat` vendored
  unmodified in `src/mcube/data/` (Q5). Atom classes come from OpenMM's own template matching.
  Unknown types raise, unless given in `polarizabilities`.
- The C4 force copies **all** NonbondedForce exceptions as exclusions: the CPU and GPU platforms
  refuse nonbonded forces whose exclusions differ (the Reference platform does not check, which
  hid this until a CPU run).
- Under 1264 an untabulated metal raises rather than run on stock terms (F2).
