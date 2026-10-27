# Ion-channel permeation PMF: decision log

A running record of the modelling decisions for the umbrella-sampling study of Ca²⁺ and two
lanthanide ions passing through the channel in `docs/channel/sequence.fasta`. Each entry gives
the decision, why, the evidence, and whether to revisit it. Same format as
[csr-decision-log.md](csr-decision-log.md). The protocol is in
[channel-us-plan.md](channel-us-plan.md).

**Aim.** Potentials of mean force for a light and a heavy lanthanide and for Ca²⁺ crossing the
channel, built with mcube's existing stages as closely as possible.

**Current state (2026-10-25; read [channel/question-and-criteria.md](channel/question-and-criteria.md) first).**
- *Question (C22, narrowed 2026-10-25):* which ion the constriction site (Glu156/240/300 ring)
  prefers, Ca²⁺, La³⁺ or Lu³⁺. Permeation is out of scope. The main calculation is an r2SCAN cluster
  exchange ΔG with SMD (C22); the z-only umbrella work below is history (C21: Ln contacts never break).
- *Protein (C1, C2, protein.md):* BD122_14570, a DUF6662 12-strand outer-membrane barrel of
  *B. diazoefficiens*, next to lanmodulin and a lanthanide ABC importer; AF2 monomer (pLDDT
  96.7), no trimer. Glu300, the most conserved lumen residue, lines the constriction.
- *System (C18):* residues 110–130 cropped, oriented, in POPE (equilibrated); Ca²⁺, La³⁺, Lu³⁺
  placed in periplasmic bulk; MM 12-6-4 umbrella runs (−35 to +35 Å, 1 Å, 2 ns) in progress.
- *QM level:* xTB/MM is out for every ion here (F7: GFN2 dynamics tears a first-shell water
  apart; Ca²⁺ fails too). g-xTB has no point-charge embedding (F8) and is set aside for Ca²⁺ vs
  Ln³⁺ (C17). The QM level at the constriction is r2SCAN/def2-SVP DFT/MM (C16 as revised by
  C17), with La and Lu on 4f-in-core lcecp-1-SVP pseudopotentials (C19); GPU DFT (gpu4pyscf,
  installed, CUDA 12 libraries missing: F11) is to be benchmarked.
- *Model checks:* 12-6-4 reproduces Amber exactly (F5); Ca²⁺ and La³⁺ bulk shells are on target;
  Lu³⁺ 12-6-4 gives CN 7.6 in OpenMM and in Amber's own sander (F4) against 9.0 in the paper's
  SI, and OpenMM-on-prmtop agrees too, so the gap is not our code. Carboxylate check (C17):
  12-6-4 ranks La vs Lu like r2SCAN (1–2 kcal/mol). The 11–13 kcal/mol Ln-over-Ca bias seen at SVP
  drops to −6…+5 kcal/mol, with no consistent sign, at TZVP + counterpoise (C17 result, 2026-10-25).
- *Next (C20):* the default 12-6-4 over-binds La³⁺ to carboxylates (experiment-fitted
  Kantakevičius 2022: C4 59.9 → 9.5; Ca²⁺ barely changes). mcube now takes `c4_set` /
  `c4_pairs`. Plan: Ca/La sensitivity run with the literature set, then fit La/Lu (and Ca) to
  acetate stability constants with one protocol, check against r2SCAN, rerun the constriction.
- *Windows (F12, 2026-10-25):* the k = 10 quick-test windows do not overlap, and the Lu³⁺ windows
  keep the pull's carboxylate contacts, so the provisional PMFs are not results. Soft-k (2.5)
  default and Ch-BE sections are queued for Ca/La (30576902–09).
- *Superseded:* C5's xTB/MM QM leg and xTB→DFT reweighting, and C16's "xTB/MM for the rest".

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
- *Superseded in part (2026-10-23):* the MM leg and the stitching stand; the xTB/MM leg does not
  (F7), and the QM level at the constriction is r2SCAN DFT/MM (C17).
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
     **sander final (30430922, 2 ns): Lu³⁺ 2.33 Å, CN 7.64.** First 1 ns: **2.33 Å, CN 7.66** (7: 34 %, 8: 66 %),
     78 exchanges/ns, residence 69 ps. Amber's own engine, on the Amber-built box with
     `lj1264 = 1`, gives what OpenMM/mcube gives (7.61–7.62). **So the difference from SI.4 is
     not in our code or in OpenMM.** It lies between today's AmberTools/sander and the authors'
     2014 AMBER 12 runs. Next step, the user's call: ask the authors how their Lu run was set up,
     framed as a question about reproducing their value, not as an error.
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
- **Ca²⁺ 1264 (30416202), 4 ns MM:** Ca–O 2.45 Å, **CN 8.04** (8: 96 %), against 2.46 / 8: on
  target. Ca IOD (30416201) timed out at 3.2 of 4 ns (6.7 ns/day on bun147, against 15–20
  elsewhere) before writing its report; rerun at 2 ns, MM only (30450561): **Ca–O 2.43 Å,
  CN 8.00**, also on target. Both Ca²⁺ sets reproduce the bulk shell.
- **La³⁺ IOD (30416203), 4 ns MM:** La–O 2.49 Å, CN 9.06 (9: 94 %, 10: 6 %), residence 171 ps,
  42 exchanges/ns. Against La 1264 (2.51 Å, CN 9.65) and the target ~2.54 Å / CN 9: IOD gets
  the CN right and Ln–O 0.05 Å short; 1264 is closer on distance and over-coordinates. Its
  default xTB/MM leg (first shell QM, 300 K, bare charges) lost SCC at step 6362 like the others
  (2.45 Å, CN 9.0 over 4.4 ps) — no new information for F7.

**C16 – DFT/MM for the desolvation step; xTB/MM for the rest of the QM/MM leg.** (user, 2026-10-23;
level and details proposed, not yet decided)
- *Revised (C17, 2026-10-23):* xTB/MM is out everywhere (F7), so the QM/MM leg is the
  constriction only, at r2SCAN DFT/MM; the rest of the path is MM 12-6-4.
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

**C17 – The question is relative selectivity, and the quantity is the free energy across the constriction.** (user, 2026-10-23)
- Goal (user): not a publication-grade PMF, but which the channel suits better, Ca²⁺ or Ln³⁺,
  and if Ln³⁺, light (La) or heavy (Lu). Relative numbers are enough.
- What decides it: the free energy of taking the ion from the hydrated vestibule into the
  constriction, i.e. dehydration *and* compensation by the site (Eisenman field strength; the
  K⁺/Na⁺ and Ca²⁺ channel and metal-site literature say the same, from memory). Dehydration
  alone would rank Ca²⁺ ≫ La³⁺ > Lu³⁺ for every channel (hydration free energies about −1500,
  −3240, −3500 kJ/mol, Marcus). Both the barrier and any well matter; with Ca²⁺ in mM excess
  over Ln, a selective well can matter as much as the barrier.
- Why QM is still needed: 12-6-4 is fitted to water only. Its C4 toward carboxylate O is the
  water value scaled by polarisability (O2 0.569 vs OW 1.444), unchecked for Glu-type
  coordination, and MM has no charge transfer, which differs between Ca²⁺ and Ln³⁺. MM may
  also sample the wrong binding mode, and then reweighting MM frames with DFT is invalid
  without saying so. So at the constriction, QM has to drive the sampling, not only correct it.
- g-xTB is not usable for Ca²⁺ vs Ln³⁺: it has no electrostatic embedding (F8), so the protein's
  field — the thing being compared — would be missing, and the missing polarisation grows with
  the ion's charge (roughly q²), so the bias differs between 2+ and 3+. A full r2SCAN PMF is
  out of reach (7–9 ps/day per window on CPU).
- **Carboxylate check, Ca and La done (30454349; it stopped at the first Lu species on an
  unconverged SCF; rerun **30473610** reuses finished species, with damping, level shift and a
  Newton fallback).** Gas phase, kcal/mol, each species optimised by each method:

  | Reaction | 12-6-4 | r2SCAN | acetate O bound: MM / DFT |
  |---|---|---|---|
  | Ca, 1 water out | −186.8 | −204.0 | 2 / **1** (2nd O 3.52 Å) |
  | Ca, 2 waters out | −171.8 | −178.7 | 2 / 2 |
  | La, 1 water out | −280.9 | −287.9 | 2 / 2 |
  | La, 2 waters out | −257.4 | −262.9 | 2 / 2 |

  La − Ca: 2 waters out, MM −85.6 vs DFT −84.2 (agree to 1.4); 1 water out, MM −94.1 vs DFT −83.9
  (MM 10 kcal/mol too favourable to La). The difference is a binding mode: r2SCAN binds Ca²⁺ by
  one carboxylate O when a water is left in its place, 12-6-4 always by two. r2SCAN aqua
  geometries: Ca–O 2.44 Å (8 waters), La–O 2.60 Å (9). No proton moved. Interim reading: the
  ion ranking agrees where both methods have the same coordination, and MM can be wrong about
  Ca²⁺'s coordination mode; Lu pending.
- **Carboxylate check, final (2026-10-24; La and Lu on lcecp-1-SVP (C19), Ca, water, acetate
  on def2-SVP; every SCF converged; no proton moved).** Gas phase, kcal/mol, each species
  optimised by each method:

  | Reaction | 12-6-4 | r2SCAN | acetate O bound (DFT) |
  |---|---|---|---|
  | Ca, 1 water out | −186.8 | −204.0 | **1** (MM: 2) |
  | Ca, 2 waters out | −171.8 | −178.7 | 2 |
  | La, 1 water out | −280.9 | −287.2 | 2 |
  | La, 2 waters out | −257.4 | −261.9 | 2 |
  | Lu, 1 water out | −297.2 | −301.5 | 2 |
  | Lu, 2 waters out | −269.0 | −274.5 | 2 |

  | Difference | 12-6-4 | r2SCAN |
  |---|---|---|
  | Lu − La, 1 / 2 waters out | −16.3 / −11.6 | −14.3 / −12.6 |
  | La − Ca, 1 / 2 waters out | −94.1 / −85.6 | −83.3 / −83.3 |
  | Lu − Ca, 1 / 2 waters out | −110.4 / −97.2 | −97.5 / −95.8 |

  r2SCAN aqua ions: Ca–O 2.44 Å (8), La–O 2.60 Å (9), Lu–O 2.36 Å (8).
  **Reading:** 12-6-4 ranks **La against Lu** like r2SCAN, within 1–2 kcal/mol, with both
  binding modes. It ranks **Ln against Ca** like r2SCAN only when Ca binds by two oxygens; when
  r2SCAN puts Ca on one oxygen with a water kept, 12-6-4 (always two oxygens) favours the
  lanthanides over Ca by **11–13 kcal/mol too much**. So an MM PMF can be trusted for La vs Lu
  more than for Ca vs Ln, and if it shows Ca and Ln³⁺ coordinating a single Glu300 differently,
  the Ca/Ln gap it gives is suspect in the lanthanides' favour. Gas phase, no BSSE correction,
  Ca and Ln on different basis families: the absolute numbers are rough; the differences are
  what this check is for.
- **Is the 11–13 kcal/mol Ca bias a small-basis artefact? (user, 2026-10-25; prepared, not run:
  Bunya needs a sign-in.)** The SVP run kept only distances, so the chain is: (1)
  `carboxylate_exchange.py --redo ALL` again, now writing `geom/<species>.xyz`; (2)
  `cp_tzvp.py` (`coord_cp.slurm`, dependent job): r2SCAN single points at def2-TZVP with La and
  Lu on lcecp-1-TZVP, plus a Boys–Bernardi counterpoise correction over every fragment (ion,
  each water, acetate), reaction energies uncorrected and corrected next to 12-6-4. Basis file:
  `bse_basis.py --basis def2-TZVP --override La,Lu=lcecp-1-TZVP --out def2-tzvp_lcecp1.nw`.
  Mixing def2 with lcecp-1 is by design: Lukanowski & Weigend, PCCP 2026 made the lcecp-k sets
  error-consistent with def2.
- **Result (2026-10-25, 30576672 → 30576673; every TZVP SCF converged):** r2SCAN/def2-TZVP
  (lcecp-1-TZVP for La, Lu) at the r2SCAN/SVP geometries, kcal/mol, gas phase:

  | Reaction | 12-6-4 | TZVP | TZVP + CP |
  |---|---|---|---|
  | Ca, 1 / 2 waters out | −186.8 / −171.8 | −193.7 / −176.2 | −190.8 / −174.3 |
  | La, 1 / 2 waters out | −280.9 / −257.4 | −281.8 / −263.2 | −280.6 / −262.5 |
  | Lu, 1 / 2 waters out | −297.2 / −269.0 | −296.9 / −276.7 | −295.1 / −276.4 |

  | Difference (1 / 2 waters out) | 12-6-4 | TZVP + CP | 12-6-4 error |
  |---|---|---|---|
  | La − Ca | −94.1 / −85.6 | −89.8 / −88.3 | −4.3 / +2.7 |
  | Lu − Ca | −110.4 / −97.2 | −104.3 / −102.1 | −6.1 / +4.9 |
  | Lu − La | −16.3 / −11.6 | −14.5 / −13.8 | −1.9 / +2.2 |

  **Reading: the 11–13 kcal/mol Ln-over-Ca bias was mostly a small-basis artefact.** At TZVP + CP,
  12-6-4 over-favours Ln³⁺ against Ca²⁺ by 4–6 kcal/mol in the 1-water step (Ca monodentate in
  r2SCAN) and under-favours them by 3–5 in the 2-water step, so there is no consistent sign. The SVP
  Ca 1-water value (−204.0) was the outlier: TZVP lowers it by 10, while every other reaction moves by 0–5.
  Counterpoise shifts are 0.3–2.9. La vs Lu stays within ~2 kcal/mol at every level. Caveats: the
  Ca 1-water geometry is from the redo, whose final SVP SCF did not converge (monodentate, Ca–O
  2.49 / 3.53 Å, as before). These are gas-phase energies only, with no thermal or solvent terms.
  Consequence for C20: the QM check alone does not demand a large Ca/Ln C4 correction. The
  experimental acetate data (c4-tuning.md: Ln − Ca ≈ −1.5 kcal/mol in water) and the Kantakevičius
  La set remain the reasons to tune.
- **Remedy if the bias holds: Merz-style pair-specific C4 tuning** (user likes the approach,
  2026-10-25): `docs/channel/c4-tuning.md` summarises Jafari et al. JPCB 2024 (per-ion acetate C4
  fitted to experimental ΔG by PMF), Fındık et al. JCTC 2024, Panteva et al. 2015, and the Merz/Hu
  ZIP work (metadynamics; umbrella sampling failed there), with the transfer plan for Ca/La/Lu.
- Plan (user, 2026-10-23): (1) MM 12-6-4 umbrella over the whole path as a quick first ranking,
  on the cropped barrel (110–130 removed) in POPE; (2) 12-6-4 vs r2SCAN on water → carboxylate
  exchange around each ion (`docs/channel/coord/`), to see whether MM ranks the ions the way
  DFT does: **30454349**; (3) GPU DFT for r2SCAN sampling at the constriction: gpu4pyscf 1.8.1
  is already in the ASH env (cupy-cuda12x 13.6, numpy 1.26.4; the dry run installs nothing),
  to be benchmarked. Fallback: alchemical Ca→La→Lu at the site and in water, r2SCAN-corrected
  at the end states.

**C18 – The system: cropped, oriented, POPE.** (2026-10-23)
- `docs/channel/af2/orient.py` on the AF2 monomer: residues 110–130 removed, the rest renumbered
  1..243 (mcube residue i = channel residue i + 130; Glu300 = A:170), barrel axis on z with the
  periplasmic end down. The belt centre is the middle of the longest charge-free stretch of the
  lipid-facing wall: −12 to +12 Å (23.9 Å, the usual outer-membrane hydrophobic thickness).
  The lower aromatic girdle (Phe328, Phe167, Tyr203) sits at −11 to −9 Å and the upper one
  (Trp345, Tyr237, Trp266, Tyr307) at +12 to +17 Å. A Kyte–Doolittle slab score had put the
  centre 6 Å lower, pulled down by polar residues near the periplasmic turns, so it was not used.
  Glu300's backbone is at z −1.8 Å, mid-membrane. The protein spans z −24 to +32 Å.
- Apo build and equilibration (`docs/channel/system/apo_config.yaml`: POPE, 0.15 M NaCl,
  310 K, 2 ns NPT + 3 ns production, GPU): **30455128**, run dir `channel/system/apo/`.
- Per ion (`docs/channel/system/`): `insert_ion.py` swaps the water nearest a point on the pore
  axis in periplasmic bulk for the ion (own chain) and removes 2 or 3 Na⁺; `ion_config.yaml`
  (skip_solvate, 12-6-4, short md, umbrella 1 Å spacing, 2 ns windows for the quick test);
  `ion.slurm` (preprocess, md, pull) and `umbrella_windows.slurm` (a block of windows per GPU).
- **Apo result (30455128, 53 min on one GPU):** 482 POPE, 57 Na⁺ / 54 Cl⁻, box 123.0 × 117.3 ×
  82.7 Å, area per lipid 60.0 Å² (production mean; NPT 60.4), P–P thickness 41.9 Å. Protein Cα
  span −22.4 to +32.2 Å from the bilayer midplane (built at −24 to +32): it stayed in place.
  Last frame `channel/system/apo/NVTtrajectory_lastframe_imaged.pdb` (z 0–82.7 Å, midplane at
  37.6, barrel Cα centroid at 42.1).
- **Ions placed (`channel/system/{CA,LA,LU}/`):** the same water, on the pore axis at lab z
  8.45 Å (CV ≈ −33.6 Å, periplasmic bulk 7 Å below the protein), swapped for the ion (chain E,
  `ion: "E:1:<ION>"`), 2 or 3 Na⁺ removed. Windows z −35 to +35 Å at 1 Å (71 per ion; +35 is
  ~7 Å above the extracellular loops, and the two bulks meet across the periodic boundary).
  Preprocess + md (0.5 ns NPT, 0.5 ns NVT, shell restrained) + pull: **30463653** Ca,
  **30463655** La, **30463656** Lu.

**C19 – La and Lu both on 4f-in-core pseudopotentials (lcecp-1-SVP).** (user, 2026-10-24)
- def2-SVP pairs La with a 46-electron core (4f⁰, nothing in the 4f shell) and Lu with a
  28-electron core (4f¹⁴ explicit), so a La/Lu comparison mixed two treatments, and every Lu
  r2SCAN SCF in the carboxylate check failed to converge (F11), the seven filled near-degenerate
  4f orbitals being the likely cause. Lu³⁺'s 4f¹⁴ is chemically inert; 4f-in-core
  pseudopotentials (Dolg–Stoll type, the lanthanide contraction built into the core) are standard
  for closed-shell Ln³⁺ coordination and hydration. Cost: 4f core polarisation and any small 4f
  share in bonding are dropped.
- BSE carries the **lcecp-1** family (`docs/channel/dft/find_ecp.py`): La core 46, Lu core 60,
  def2-style valence bases. With `lcecp-1-SVP` both Ln³⁺ have 8 explicit electrons and 33 AOs.
  H, C, N, O, S and Ca stay def2-SVP.
- File: `/scratch/user/uquadhik/channel/dft/def2-svp_lcecp1.nw` (`bse_basis.py --override
  La,Lu=lcecp-1-SVP`). Every lanthanide DFT/MM from now on passes it as `dft_basis_file` and
  `dft_ecp_file`; the def2-SVP numbers for La (C16 benchmark, the C17 Ca/La table) are on the
  old basis.
- Carboxylate check rerun with it, La and Lu recomputed (Ca, water, acetate reused):
  **30526116** (the def2 Lu redo 30523294 was cancelled).

**C20 – Tune the ion–carboxylate C4 the Merz/Hay way; mcube ships literature sets and takes
explicit overrides.** (user, 2026-10-25)
- **What the literature gives** (`docs/channel/c4-tuning.md`):

  | Pair (TIP3P, Asp/Glu O2) | Default α-ratio C4 | Kantakevičius 2022 (chelates) | Jafari 2024 (acetate) |
  |---|---|---|---|
  | Ca²⁺–O | 34.4 | 29.2 | 36.27 |
  | La³⁺–O | 59.9 | **9.5** | — |
  | Lu³⁺–O | 98.1 | — | — |

  The Ca²⁺ side is close to experiment-calibrated already (Jafari's acetate fit lands 5 % above
  the default; the chelate fit 15 % below). The La³⁺ side is not: the chelate fit cuts La–
  carboxylate C4 six-fold, and with it Kantakevičius reproduce lanmodulin's La > Y > Ca order
  (La − Ca −40 vs −46 kJ/mol experiment) while absolute binding stays too strong for every ion.
  So the C17 finding (12-6-4 favours Ln³⁺ over Ca²⁺ more than r2SCAN) has an independent,
  experiment-based counterpart: default 12-6-4 over-binds La³⁺ to carboxylates. **No published
  tuned value exists for Lu³⁺** (or any heavy lanthanide).
- **mcube (2026-10-25):** `c4_set: <name>` picks a shipped set from `src/mcube/data/c4_pairs.yaml`
  (`kantakevicius2022_chbe_tip3p`, `jafari2024_acetate_tip3p`, each with source, table, water
  model and fitting target); `c4_pairs: {Ca2+: {O2: 29.2}}` gives explicit pairs and wins over the
  set. Both are opt-in: the default stays the plain Li/Song/Merz 12-6-4, because every tuned set
  is specific to one ligand chemistry and one water model. Implemented as a table added to the
  existing C4 force (override − default for listed pairs only); QM atoms stay off; carried to
  every later stage like `metal_parameters`. 11 new tests (`TestPairC4`).
- **Plan:**
  1. Sensitivity first: rerun the constriction windows (z −15…+14 Å) for Ca²⁺ and La³⁺ with
     `c4_set: kantakevicius2022_chbe_tip3p` and compare their PMFs with the default ones. If the
     Ca/La ranking moves, the C4 choice decides the answer and must be fitted properly.
  2. Lu³⁺ needs its own value. Fit La³⁺ and Lu³⁺ (and refit Ca²⁺ as a check) to experimental
     acetate stability constants with Jafari's PMF protocol (M–carboxyl-C distance, 0.05/0.1 Å
     windows, ≥ 6 ns, 3 replicates, Jacobian and standard-state correction, the 3+/1− Coulomb
     tail handled), so all three ions come from one protocol and one target type. Experimental
     log β₁ for La/Lu acetate still to be taken from NIST SRD 46 / Martell–Smith and corrected to
     I = 0. Interim, clearly labelled as an assumption: scale Lu's default by La's ratio.
  3. Check every fitted set against the r2SCAN exchange energies (C17, TZVP + counterpoise) and
     for the binding mode (Ca²⁺ mono- vs bidentate).
  4. Rerun the constriction windows for all three ions with the chosen set; that pair of PMFs is
     the MM answer to the selectivity question.

**C21 – Reseed the Ln³⁺ windows: restrain the backbone, pull from both bulks, test for
hysteresis.** (proposed 2026-10-25, not yet approved; La soft/Ch-BE 30576904/05/08/09 and Lu k = 10
blocks 30576666/67 held; the Ca jobs run)
- **What the pulls did** (pull_cv.csv and seed frames; channel numbering = mcube + 130):

  | Ion, direction | Lumen contacts carried | Residue moved |
  |---|---|---|
  | Ca²⁺, +35 → −35 | none from +8 to −25 Å (7.3 waters throughout); backbone O of Lys145/Asn359 at the extracellular mouth (+9…+19), Glu209 at the periplasmic mouth (−26…−34) | — |
  | La³⁺, −35 → +35 | Glu209 + Asp207 from −31 to −8 Å (2–2.8 O; ion trails z0 by 1.2–1.7 Å), released at −6; Glu156 from −3 to +6 (trails by up to 2.2 Å); Glu265 from +12 to +20 | the Asp207/Glu209 loop (at rest z ≈ −21…−27) dragged ~8 Å up to −13…−16 (its Cα too) |
  | Lu³⁺, +35 → −35 | Gln306 at +9…+13; Glu240 from +5 to −11 Å (ion trails z0 by up to 2.5 Å), released at −11.5 | the Glu240 strand pulled down 3.4 Å (Cα), carboxylate 7 Å |

  Conserved Glu300 (z ≈ −2…−5) is touched by neither Ln in the seeds. The Ln carboxylates that
  matter in the lumen are Glu156 and Glu240, with Glu300 alongside: a ring of carboxylates near z ≈ 0.
- **Why the windows remember the pull:** the 12-6-4 Ln–carboxylate bond is too strong to break in
  1 ns, the pull distorts the protein so the window inherits a strained state, and the window's
  100 ps equilibration restrains the seed's own coordination shell (`coordination_force`), which
  locks the pulled contacts in before sampling starts.
- **Proposal:**
  1. **Backbone restraint** (new mcube input, e.g. `protein_restraint: {selection: "chain A and name
     CA", k: 1.0}` kcal/mol/Å² to the md lastframe) in the pull and in every window. This stops strands and
     loops being dragged (the 8 Å loop move, the 3.4 Å strand move) and leaves side chains free to
     coordinate. It is common in β-barrel PMFs; it removes real backbone motion, which is a stated limitation.
  2. **Two pulls per Ln, one from each bulk** (+35 → −35 and −35 → +35, 5 Å/ns, with the restraint;
     ~14 ns ≈ 5 h each on an L40). Each window is seeded from the pull that reached it from the
     nearer bulk, so contacts are gathered over half the path at most.
  3. **No shell restraint during window equilibration for the umbrella stage** (an input to switch it off), so the
     seed's contacts are not locked in.
  4. **Hysteresis test before the full set:** run the windows z0 −6…+6 (13) from both pulls with
     k = 2.5. Accept the protocol if cn_protein and the PMF from the two seed sets agree within ~0.3 O and
     ~1 kcal/mol. Otherwise those windows need a coordination CV (2D z, CN) or metadynamics, as for ZIP.
  5. Then the full soft-k (default C4) and Ch-BE sets for La, and the soft set for Lu, on the new seeds.
     Ca stays as queued; it carries no lumen contacts. Its mouth contacts are backbone O, which the
     restraint would not change.
- **Cost:** mcube change ~40 lines plus tests; 4 pulls (~20 GPU-h); 13 × 2 × 2 hysteresis windows
  (~20 GPU-h); then the sets already planned.
- **Revised after discussion (user, 2026-10-25):** prefer a *flat-bottom* Cα restraint (free within
  ~2× md RMSF, harmonic beyond; report how often it engages) to a harmonic one. The two-seed test
  is not sufficient alone: seeds can agree because neither binds or unbinds. Strengthened test: three
  seeds per window (from the pull, from the opposite pull, and a hydrated swap-in seed), 3–5 windows across
  the carboxylate ring, 5 ns, accepted only if the cn_protein distributions converge, the two halves agree,
  and contacts change several times. **First a residence-time diagnostic,** to see whether z-only US can
  converge at all.
- **Residence-time test (queued 2026-10-25, `short` on gpu_viz, 6 h each):** one 10 ns window at z0 = 0,
  k = 2.5, seeded from the bound k = 10 window 035, frames every 25 ps. Sections `umbrella_res`
  (default C4) and `umbrella_res_chbe`. Jobs: **30581006** La default, **30581007** La Ch-BE (La–O2
  9.5), **30581008** Lu default, **30581009** Lu with La–O2 9.5 scaled to Lu: 98.1 × 9.5/59.9 =
  **15.6** (an assumption; no fitted Lu value exists, C20). Analysis: `docs/channel/system/residence.py
  <ION>/umbrella_res_035` (per-residue contact fraction, episodes, mean/longest lifetime, cn_protein
  changes from cv.csv). Reading: lifetimes ≪ 1 ns → z-only US with the flat-bottom restraint;
  lifetimes ≳ window length → 2D (z, carboxylate-O count) or metadynamics at the constriction; Ch-BE
  releasing where default does not → the C4 choice decides sampling as well as the answer.
- **La³⁺, default C4 (30581006, 10 ns at z0 = 0):** **no exchange.** Glu156 is bound for all 10 ns
  (one episode, bound at start and end). Glu240 binds at ~3.6 ns and stays for the remaining 6.4 ns.
  cn_protein goes 1.7 → 3.6 O and waters 7.0 → 5.3, with z −0.4…−1.1 and r ≈ 5–6 Å, so La sits at
  the wall between the two carboxylates. Glu300 is not within 3 Å. The cv.csv "changes" count (14724) is the
  rounded switching-function CN flickering around x.5, not contact events; the frame-based episodes are
  the measure. Reading: under the default C4 a La–Glu contact outlives a 10 ns window, so z-only
  windows cannot equilibrate coordination here, and the window gained a second carboxylate rather
  than releasing the first.
- **La³⁺, Ch-BE C4 (La–O2 9.5; 30581007, 10 ns at z0 = 0):** **no exchange either.** Glu156 is bound
  for all 10 ns (one episode), with no second carboxylate. cn_protein stays 1.7–1.8 (both halves) and waters
  7.0. Cutting C4 six-fold stops the second Glu from binding but does not free the first. What holds a
  3+ ion on a carboxylate is mostly Coulomb, not the C4 term. So for La³⁺, z-only windows will not
  equilibrate coordination at the constriction with either C4. The sampling scheme has to treat the
  ion–carboxylate contact explicitly (C21 options: 2D z × carboxylate-O count, or metadynamics on it).
- **Lu³⁺, interim (still running, read from frames alone; 30581008 default 7.3 ns, 30581009 Lu–O2 15.6
  6.1 ns):** Glu240 has been bound for the whole run in both (one episode each). With the weaker C4,
  Glu220 has joined for the last 2 ns, so Lu gains a carboxylate where La did not. **All four runs: no
  Ln–carboxylate contact has broken in 6–10 ns, with either C4.**
- **Recommendation (pending the user's decision):** z-only US is out for Ln³⁺ at the constriction. The
  slow coordinate is the number of carboxylate O on the ion, so restrain it explicitly with **2D umbrella
  sampling (z × n_COO)**. The bias force already computes coordination numbers as collective variables,
  so the carboxylate-O count is a third CN group with a second harmonic bias; 2D WHAM, then the free
  energy projected on z. Grid: z −6…+6 Å (1 Å), n_COO 0, 1, 2, 3; ~52 windows × 2 ns per ion and C4,
  ≈ 20 GPU-h per set. Outside ±6 Å, 1D soft-k windows with the flat-bottom Cα restraint. Metadynamics
  (OPES/WT-MetaD on n_COO inside z windows) is the alternative, but needs PLUMED or a new OpenMM bias
  in mcube and is harder to check. Ca²⁺ needs no 2D treatment (no lumen contacts), but should get the
  same grid at the constriction for a like-for-like comparison.
- **Ca²⁺ Ch-BE section (30576902/03, 30 windows, k = 2.5), WHAM:** overlap still below 0.05 after z0 −10,
  −5, −1, +3, +13 (min 0.0). The soft spring lets Ca slide to local minima (z0 −1 sits at −1.49,
  z0 0 at +0.41; −4 and −3 both at −3.45), leaving 1–2 Å gaps. Provisional profile: 2–8 kcal/mol across
  the barrel (−15…+8 Å), lowest at +10 Å (backbone-O contacts at the extracellular mouth), with a
  peak of 8.8 at +13. Ca touches protein O in a few windows only (−9, −7, −1, +3/+4, +9…+14). Needs
  extra windows at the gaps (or k ≈ 5 there) before it is a result; compare with umbrella_soft (16/30).

**C22 – Site affinity only; cluster exchange ΔG at the Glu ring is the main calculation.**
(user, 2026-10-25; the reference is [channel/question-and-criteria.md](channel/question-and-criteria.md))
- **Scope:** which ion the constriction site prefers. Permeation (exit barriers, profiles along z) is
  out. Cancelled as permeation work: La soft/Ch-BE 30576904/05/08/09, Lu k = 10 blocks
  30576666/67 and Ca soft 30576906/07 (21–24/30 done; the Ca Ch-BE set and all k = 10 windows stay on
  disk). The Lu residence runs 30581008/09 continue; their frames feed the clusters.
- **Method** (`docs/channel/coord/site_exchange.py`, `site.slurm`), after Angelova, Dudev et al.,
  PCCP 2023 (PMC10094714). They use M06-2X/6-311++G(d,p)/SDD, SMD with ε 4 and 29 for the pocket, and
  Ca-binding EF-hand sites. Their 3-carboxylate site prefers La over Ca by 14–20 kcal/mol, the 4-carboxylate
  site by 30–48, and the 1-carboxylate site prefers Ca or neither. La vs heavier Ln moves by only ±2–7, not
  monotonically (Dy best, Lu slightly negative).
  - Ours: the site is cut from a bound MD frame and minimum-imaged around the ion. Glu becomes butanoate
    capped at CA. Coordinating side chains are free except CA/HA/caps; non-coordinating ring Glu are frozen.
    First-shell waters are within 3.3 Å. Each ion is swapped into the same cluster and optimised at
    r2SCAN/def2-SVP (lcecp-1-SVP); single points at def2-TZVP (lcecp-1-TZVP) in the gas phase and in SMD
    diethylether (ε 4.2) and propanonitrile (ε 29). Aqua ions (Ca 8, La 9, Lu 8 waters, C17 geometries)
    and water are in SMD water. Ion cavity radii are set explicitly (Ca 2.0, La 2.2, Lu 1.98 Å), because
    PySCF's SMD has a 2.0 Å placeholder for La/Lu; the cavity is built by element, so ECP core charges do
    not shift the lookup. No thermal corrections.
  - Variants test the pocket-charge sensitivity (criterion 2): `ring` keeps Glu156/220/240/300 (−4; the
    non-coordinating ones frozen); `min` keeps only the coordinating Glu.
  - Sites from La frames: `la_def_ring`/`la_def_min` (Glu156 + Glu240 coordinating, 6 waters) and
    `la_chbe_ring`/`la_chbe_min` (Glu156 only, 8 waters). Lu frames follow when the residence runs end.
  - A bug found on the way: the frames are wrapped atom by atom, so 4 of La's 7 first-shell waters were
    across the box; fixed by minimum imaging in `build`.
- **Calibration (criterion 4):** `site_exchange.py calibrate` applies the same SMD-water protocol to
  the C17 acetate complexes, against experiment (La − Ca ≈ −1.4…−1.7, La − Yb ≈ +0.2 kcal/mol).
- **Queued (general, 32 cores):** 30598663 calibration, 30598664 `la_chbe_min` (smoke test before
  the larger sites).
- **SMD → IEF-PCM (2026-10-25):** the calibration crashed at its first solvated single point. PySCF's
  pip build lacks the SMD non-electrostatic (CDS) library ("compile with -DENABLE_SMD=ON"), and
  gpu4pyscf's SMD would need a GPU. So solvation is IEF-PCM, electrostatics only, at the same
  dielectrics: site ε 4 and 29, water 78.36. The default radii are 1.2 × modified Bondi, with the
  ion radii set explicitly. The site complexes differ only in the ion, so their CDS terms nearly
  cancel; the aqua ions' do not quite, so this adds an uncertainty of a few kcal/mol for Ca vs Ln.
  Checked: Ca(H₂O)₈²⁺ PCM-water solvation −175 kcal/mol, converged in 78 s.
- **Lu³⁺ default residence, final (30581008, 10 ns):** Glu240 is bound for the whole run (one
  episode). cn_protein is 1.5 → 1.9 and waters 6.9. Lu sites: `lu_def_ring` (−4; Glu240 coordinating,
  Glu156/220/300 frozen; 7 waters) and `lu_def_min` (−1; Glu240 only).
- **Running:** 30599579 calibration (PCM), 30599596–600 `la_def_ring`, `la_def_min`, `la_chbe_ring`,
  `lu_def_ring`, `lu_def_min`. 30598664 `la_chbe_min` is still on the SMD code; it saves its first
  geometry, then fails at the single point, and must be resubmitted (it continues from report.json/geom).
  The Lu scaled-C4 residence run 30581009 is still going.
- **Calibration result (30599579; r2SCAN/TZVP, IEF-PCM water, C17 acetate geometries, no thermal
  terms).** Binding of acetate, kcal/mol (1 water out / 2 waters out):

  | Ion | PCM water | Experiment (c4-tuning.md) |
  |---|---|---|
  | Ca²⁺ | −16.4 / −7.2 | −1.3 to −1.6 |
  | La³⁺ | −22.4 / −14.5 | −3.0 |
  | Lu³⁺ | −25.2 / −15.8 | −3.2 (Yb) |

  Relative, best mode each: **La − Ca −6.0 (exp. −1.4…−1.7); Lu − Ca −8.8; Lu − La −2.8 (exp. ≈ −0.2).**
  - The direction is right (Ln over Ca, and slightly heavy over light).
  - The protocol over-favours Ln³⁺ against Ca²⁺ by ~4.5 kcal/mol, and Lu against La by ~2.6, per carboxylate.
  - The absolute values are 10–20 kcal/mol too negative, mostly because the association entropy (TΔS) and CDS
    terms are missing.
  - So criterion 3 applies as a *systematic offset*, not just scatter: a site result counts for Ln over Ca only if it
    exceeds about 4.5 kcal/mol per coordinating carboxylate, and for Lu over La only if it exceeds about 2.6 per carboxylate.
  - Caveat: the Ca 1-water geometry comes from a non-converged SVP optimisation (C17).
- **Lu³⁺ scaled-C4 residence, final (30581009):** Glu240 is bound all 10 ns. Glu220 joins at ~4.1 ns and
  stays (5.9 ns). cn_protein is 1.3 → 2.1. Sites: `lu_chbe_ring` (−4; Glu220/240 coordinating) and
  `lu_chbe_min` (−2), jobs 30600588/89. `la_chbe_min` was resubmitted on PCM as 30600578.
- **First site result, `lu_def_min`** (pocket −1: Glu240 alone, 7 waters; every optimisation and SCF
  converged). Exchange ΔE, kcal/mol (negative = the site prefers the second ion):

  | | gas | ε 4 | ε 29 |
  |---|---|---|---|
  | Ca → La | −66.5 | **+33.5** | **+8.5** |
  | Ca → Lu | −101.6 | +22.5 | −3.6 |
  | La → Lu | −35.1 | −10.9 | −12.2 |

  - **Ca over La at a one-carboxylate site,** in both dielectrics, and more so after the calibration offset,
    which favours Ln. This matches Dudev's 1-carboxylate site (Ca-selective at ε 4, about neutral at ε 29).
  - **Ca vs Lu depends on the dielectric** (+22.5 / −3.6, both within the offset): no answer.
  - **Lu over La by 11–12 kcal/mol,** about 9 after the 2.6 offset.
  - The gas phase favours the 3+ ions by 35–100; solvation reverses that. This is the expected pattern, and
    it is why the dielectric and the pocket charge decide the result.
  - This is the stripped variant. The ring variant (−4) is the test of criterion 2; Dudev's −4 site
    swings to strongly Ln-selective.
- **`lu_def_ring` is slow:** the first optimisation stopped on an SCF that would not converge (net −2
  cluster with Ca), and it restarted with second-order SCF. Step 5 after ~4 h, so it will need resubmitting at 24 h
  (it continues per ion).
- **All four `min` sites done (2026-10-26).** Raw exchange ΔE in kcal/mol at ε 4 / ε 29 (negative =
  the site prefers the second ion); n = coordinating carboxylates; ✗ = that ion's optimisation hit
  the step limit:

  | Site | n | Ca → La | Ca → Lu | La → Lu |
  |---|---|---|---|---|
  | `la_def_min` (Glu156 + Glu240) | 2 | +7.4 / −0.4 | −0.2 / −8.1 | −7.6 / −7.7 |
  | `lu_chbe_min` (Glu220 + Glu240; Ca ✗) | 2 | +5.5 / −3.5 | −7.4 / −16.2 | −13.0 / −12.7 |
  | `la_chbe_min` (Glu156; La ✗) | 1 | +31.0 / +6.6 | +19.2 / −5.4 | −11.8 / −12.0 |
  | `lu_def_min` (Glu240) | 1 | +33.5 / +8.5 | +22.5 / −3.6 | −10.9 / −12.2 |

  Corrected for the calibration bias (add 4.5n to the Ca → Ln columns, 2.6n to La → Lu):
  - **Ca over La** at every site and in both dielectrics: +8.6…+38 kcal/mol, so the site does *not* prefer La.
  - **Ca vs Lu** is not decided: −7.2…+27, with the sign flipping with the dielectric at the 2-Glu sites.
  - **Lu over La** at every site: −2.4…−9.6, beyond the offset at 3 of 4 sites.
  - Sites from La and Lu frames agree in sign (criterion 1).
  - **But these are the stripped pockets (−1/−2).** The ring variants (−4) test criterion 2, and Dudev's −4
    site is strongly Ln-selective.
- **Ring sites:** `la_def_ring` (30599596) is still running at 9 h, past Ca. `la_chbe_ring`, `lu_def_ring` and
  `lu_chbe_ring` ran out of memory at 96 GB, after their Ca step (PySCF max_memory 80 GB plus the dense IEF-PCM matrices
  for 74–77 atoms). Resubmitted with `--mem=180G` as 30609315–17; they continue from the saved Ca
  geometries. `lu_def_ring`'s Ca geometry did not converge in 150 steps, so treat it with care.
- **Gas-phase optimisation fails for the ring clusters (2026-10-26).** In `la_def_ring` the La optimisation,
  on second-order SCF after a first SCF failure, alternates between two energies 0.34 Eh (~210 kcal/mol)
  apart on successive steps, and its trust radius keeps shrinking: the SCF is landing on different solutions.
  `lu_chbe_ring` restarted the same way. The ring clusters are anionic (four Glu⁻, two uncoordinated
  and frozen; net −2 with Ca, −1 with Ln), and in vacuum the excess charge is barely bound. **Fix:**
  `--opt-eps 4` optimises in IEF-PCM ε 4, which binds it. To keep one protocol for the ring vs min
  comparison (criterion 2), **all eight sites rerun** with PCM optimisation in `<site>_p4/`
  (30615271–78; ring sites 180 GB). The gas-optimised min results above stay as a check on how much the
  geometry protocol matters. Ring jobs 30599596 and 30609315–17 were cancelled.
- **PCM gradients segfault for the ring clusters (2026-10-26):** all four ring `_p4` jobs died at their
  first PCM gradient (RC 139), even for Ca, which has no ECP. With ~71 atoms × 302 cavity points the
  grid-point integrals scale as n_ao² × n_grid (~10¹⁰). The min sites (~40 atoms) run.
  `--opt-lebedev 13` (74 points per atom) is now used for the ring optimisations only; the single points
  keep order 29. Resubmitted: 30622024–27 (`la_def_ring`, `la_chbe_ring`, `lu_def_ring`,
  `lu_chbe_ring`, all `_p4`, 180 GB). The PCM `opt_eps` 4.0 is confirmed in the min reports. `la_def_min_p4`
  had one SCF collapse (a step 7.4 Eh higher); geomeTRIC should reject it, so watch it.
- **Order 13 raised KeyError 74:** PySCF's PCM has switching parameters only for some point
  counts. The ring sites were resubmitted at `--opt-lebedev 15` (86 points) as 30623448–51; they run without crashes.
- **First PCM-optimised result, `lu_def_min_p4`** (Glu240 only, n = 1; all converged). ε 4 / ε 29,
  gas-optimised values in brackets:

  | | ε 4 | ε 29 |
  |---|---|---|
  | Ca → La | +24.9 (+33.5) | +0.4 (+8.5) |
  | Ca → Lu | +15.0 (+22.5) | −10.7 (−3.6) |
  | La → Lu | −10.0 (−10.9) | −11.1 (−12.2) |

  - Optimising in PCM moves Ca vs Ln by ~8 kcal/mol toward the lanthanides; La vs Lu barely changes (≤1).
  - After the offset: Ca over La (+29.4 / +4.9); Ca vs Lu still undecided (+19.5 / −6.2); Lu over La (−7.4 / −8.5).
  - The signs are the same as with gas optimisation.
- **`la_chbe_min_p4`** (Glu156 only, n = 1; all converged), ε 4 / ε 29, gas-optimised in brackets:
  - Ca → La +30.0 / +4.9 (+31.0 / +6.6)
  - Ca → Lu +19.6 / −6.1 (+19.2 / −5.4)
  - La → Lu −10.3 / −11.0 (−11.8 / −12.0)

  Here PCM optimisation changes things by ≤ 1.7, against ~8 for `lu_def_min`. The two one-Glu sites, from a La frame
  and a Lu frame, now agree within ~5 kcal/mol (criterion 1).
  - After the offset: Ca over La (+34.5 / +9.4); Ca vs Lu undecided (+24 / −1.6); Lu over La (−7.7 / −8.4).
- **The session moved** to bun147:19973 (the bun155 session ended). All jobs kept running.
- **Two-Glu min sites, PCM-optimised (2026-10-27); n = 2; in both, Ca's optimisation ended unconverged at 150
  steps.** ε 4 / ε 29, gas-optimised in brackets:

  | Site | Ca → La | Ca → Lu | La → Lu |
  |---|---|---|---|
  | `la_def_min_p4` (Glu156 + Glu240) | **−18.1 / −25.9** (+7.4 / −0.4) | −24.4 / −32.4 (−0.2 / −8.1) | −6.3 / −6.5 (−7.6 / −7.7) |
  | `lu_chbe_min_p4` (Glu220 + Glu240) | +2.9 / −6.7 (+5.5 / −3.5) | −7.5 / −16.9 (−7.4 / −16.2) | −10.4 / −10.2 (−13.0 / −12.7) |

  - `lu_chbe_min` hardly moves and keeps the signs: after the offset, Ca over La (+11.9 / +2.4), Ca vs Lu undecided
    (+1.5 / −7.9), Lu over La (−5.2 / −5.0).
  - **`la_def_min` jumps ~25 kcal/mol toward Ln** in the Ca → Ln columns only (La → Lu barely changes). That
    points at the Ca energy: an unconverged Ca geometry has too high an energy, which pushes exactly these
    columns negative. So it is not taken as a result.
  - The Ca step of both sites is rerunning from its unconverged geometry with `--maxsteps 300`, keeping the
    La/Lu geometries and energies: `<site>_p4c`, jobs **30680249** (`la_def_min`) and **30680250** (`lu_chbe_min`).
  - Ring sites: three have finished Ca, and `la_def_ring`/`la_chbe_ring` have La geometries ✗ (150 steps); `lu_chbe_ring`
    has La done; `lu_def_ring` is on La. Session now bun151:34126.
- **`la_def_min_p4c`: the Ca rerun converged (40 more steps), so the swing was the unconverged Ca.** All three ions
  converged; ε 4 / ε 29:
  - Ca → La +8.2 / +0.0
  - Ca → Lu +1.9 / −6.5
  - La → Lu −6.3 / −6.5

  The −26 kcal/mol swing is gone, and the values are back within ~6 of the gas-optimised ones (+7.4 / −0.4; −0.2 / −8.1).
  After the offset (n = 2: +9.0 for Ca → Ln, +5.2 for La → Lu):
  - **Ca over La** (+17.2 / +9.0)
  - **Ca over Lu** at this site (+10.9 / +2.5)
  - **La vs Lu within the offset** (−1.1 / −1.3), so undecided here.

  Lesson: a site counts only once all three ions have converged; the `_p4` value with unconverged Ca was off by ~26.
- **`lu_chbe_min_p4c`: the Ca rerun converged; the values are unchanged** (Ca → La +2.9 / −6.6, Ca → Lu −7.5 / −16.8,
  La → Lu −10.4 / −10.2). After the offset (n = 2): Ca over La (+11.9 / +2.4), Ca vs Lu undecided (+1.5 / −7.8),
  Lu over La (−5.2 / −5.0).
- **First ring site: `lu_chbe_ring_p4`** (Glu220 + Glu240 coordinating, Glu156/300 frozen, pocket −4, n = 2; all three
  ions converged). ε 4 / ε 29:
  - Ca → La −13.4 / −4.1
  - Ca → Lu −31.8 / −19.9
  - La → Lu −18.4 / −15.8

  After the offset:
  - **Ca vs La is now dielectric-dependent** (−4.4 / +4.9).
  - **Lu over Ca** (−22.8 / −10.9), beyond the offset.
  - **Lu over La** (−13.2 / −10.6).

  **Criterion 2 against its min variant (same frame, same coordinating Glu):** the two extra, uncoordinated Glu⁻ move
  every column toward the 3+ ions, by 16 / 7 (Ca → La), 24 / 3 (Ca → Lu) and 8 / 6 (La → Lu), as Dudev's
  pocket-charge result predicts.
  - Ca over La does **not** survive at ε 4.
  - Ca vs Lu turns into Lu over Ca.
  - Lu over La holds and grows.

  So for Ca vs Ln the answer depends on how much of the second-shell Glu charge reaches the ion. That is the
  question the QM/MM snapshot single points (embedding in the real protein) would settle.
- **Ring status:** `la_def_ring` ended with both La and Lu unconverged at 150 steps (both need a rerun);
  `la_chbe_ring` is on La; `lu_def_ring` has La ✗ and is on Lu.
- **Metadynamics stage (support only):** `mcube metadynamics` (commit 6e85c16) is well-tempered MTD on
  n_COO with z held, or 2D, on `openmm.app.Metadynamics` with umbrella's CVs (ASH's wrapper uses the same
  class with a different CN function and no CV logging). `fes` reports integer-state free energies
  and a 50/75/100 % convergence check. 222 tests pass. It is used only if the cluster ranking depends on the binding
  mode.

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
  (30421227 crashed on an empty ECP entry for H/O, fixed in 3a4853b; rerun **30423347**.)
- **La³⁺ result (6 frames from the GFN2 xTB/MM run that failed at step 6328):**

  | Frame (ps) | Gap S1 eV: r2SCAN / GFN2 / g-xTB | Grad RMSE vs r2SCAN, kcal/mol/Å: GFN2 / g-xTB | Outward force on 1st-shell O: r2SCAN / GFN2 / g-xTB |
  |---|---|---|---|
  | 0.2 | 7.18 / 3.91 / 10.14 | 20.5 / 12.1 | 36 / 10 / 21 |
  | 1.4 | 6.94 / 3.99 / 9.92 | 16.7 / 9.5 | 43 / 23 / 34 |
  | 2.6 | 7.14 / 3.89 / 9.79 | 15.3 / 9.5 | 35 / 16 / 25 |
  | 3.9 | 6.68 / 3.74 / 9.82 | 28.2 / 14.8 | 58 / 27 / 42 |
  | 5.1 | 6.86 / 3.61 / 9.77 | 26.3 / 16.2 | 42 / 11 / 21 |
  | **6.3** | **0.15 (SCF unconverged) / 0.42 / 3.70** | 30.8 / 16.0 | 95 / 51 / 75 |

  - **The SCC failure is a broken O–H bond.** In the last frame, ~20 steps before GFN2 failed,
    one first-shell water has an O–H of **2.60 Å** (all others 0.93–1.05 Å), with the lost H
    ≥ 2.3 Å from any other O in the cluster. r2SCAN's gap collapses there too (0.15 eV, as for an
    O–H radical pair), so the gap closing is real for that geometry. The fault is that GFN2/MM
    dynamics tore a water apart in the first place.
  - **Before that, no method's gap is small.** Bare charges (S1_pc) change GFN2's gap by < 0.3 eV.
  - **Forces:** g-xTB's error against r2SCAN is about half GFN2's (9.5–16 vs 15–31 kcal/mol/Å
    RMSE), but still large. Both under-push the first shell outward (r2SCAN 35–58, GFN2 10–27,
    g-xTB 21–42 kcal/mol/Å), and GFN2 most, which fits its short Ln–O.
  - The Lu frames and the gas-phase optimisations are still running (the job may hit its 12 h
    limit before the optimisations).
- **Ca²⁺ fails the same way.** The Ca 1264 xTB/MM leg (30416202) lost SCC at step 9433
  (gap 0.06 eV), with Ca–O 2.17 Å / CN 7.0 against 2.46 / 8. So GFN2-xTB/MM is unreliable here
  for every ion, not only Ln³⁺.
- Reading: g-xTB is better than GFN2 but not good enough to justify writing point-charge
  embedding for an unfinished method. DFT/MM is available now.

## Planned checks, not yet done

- Q1: answered 2026-10-22 (`docs/channel/probe_env.py` on bun147; the one-off probe was retired
  2026-10-23 and is in git history up to 7ab4f1d): the ASH env has Python 3.11.14,
  **OpenMM 8.4 with Reference, CPU and CUDA** (no HIP, no OpenCL), pyscf 2.13.0, mdtraj 1.11.1;
  pymbar, openmm-plumed and basis_set_exchange are absent. Jobs import mcube from the scratch
  checkout (`PYTHONPATH=$S/mcube/src`), and `mcube.umbrella` and `mcube.ions` import there.
  Still open: GPU vs CPU speed for the MM side of QM/MM (needs a real window).
- Q2: answered: the ASH env's OpenMM 8.4 `amber14/tip3p.xml` has the CA, LA and LU templates, with
  the same stock radii as the local copy (F2 holds on Bunya: La³⁺ Rmin/2 1.758, Lu³⁺ 1.625 Å).
- Q3: one QM/MM window in the membrane system: timing, energy drift, box stability (QM/MM in a
  membrane has never been run in mcube).
- Q4: answered 2026-10-23 (`docs/channel/af2/pore.py`, C2): radius 0.8–1.1 Å over z −14…−5 Å
  with the N-terminal segment in the lumen; 2.3–3.9 Å elsewhere. With 110–130 cropped (C18) the
  constriction is Glu300's ring; its radius in the cropped, equilibrated model is still to measure.
- Q5: answered (`docs/channel/probe_amber.sh`, retired 2026-10-23, in git history up to 7ab4f1d). AmberTools lives inside the ASH conda env
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

**F9 – preprocess could not read residue numbers past 9999.** (2026-10-23, fixed in 05bdd7c)
- All three per-ion jobs (30463653/55/56) died in 15 s in `find_coordinators`: "invalid literal
  for int(): 'B66D'". The membrane box has ~20k waters, and OpenMM (so ASH and mcube) writes
  residue numbers past 9999 in *shifted hexadecimal* (`pdbfile._formatIndex`: 10000 → `A000`,
  10501 → `A1F5`), not hybrid-36, and decodes them on reading. `coordination_shell` and
  `atom_names` read columns 23–26 as plain text. New `preprocess.pdb_resseq` inverts OpenMM's
  encoding (unique up to 34576 residues per chain); `tests/test_pdb_resseq.py` checks it against
  OpenMM's own reading of a 12k-water file (passes in the ASH env). Any metal preprocess with
  `skip_solvate` on a large box hit this; fresh solvation never did. Resubmitted:
  **30473607** Ca, **30473608** La, **30473609** Lu (then F10).

**F10 – An ion inserted into md's last frame lost the lipids' bonds.** (2026-10-24, fixed in 70761bf)
- After F9, preprocess passed, then md failed: "No template found for residue 243 (POP) … no
  bonds between its atoms". md's imaged last frame (written by mdtraj) has no CONECT records,
  and POPE is not a standard residue, so its bonds exist only in CONECT. The same frame also
  writes water residue numbers modulo 10000. `insert_ion.py` now takes the topology and bonds
  from the builder's `finalsystem.pdb` (60250 CONECT, the same 124189 atoms in the same order,
  checked) and only positions and box from the last frame.
- The README route "a new partner in an equilibrated bilayer" (write into a copy of the last
  frame) would hit the same thing for any lipid; worth a note there.
- GPU queue (2026-10-24): `gpu_cuda` estimated a 27 h wait (another project's 16 queued GPU jobs
  share our per-user `gpu` QOS of 4). `sxm` and `mig` were rejected; `gpu_viz` L40s (QOS `viz`,
  2 per user, but the VS Code session holds one) started at once and run OpenMM CUDA. Per-ion
  jobs: **30478619** La (running), **30478620** Lu (next `viz` slot), Ca moved to QOS `short` on
  gpu_viz: **30481824** (estimated 03:07). `running-on-bunya.md` corrected.
- **Second half of F10:** the imaged last frame also has molecules split across the periodic
  boundary. With finalsystem's bonds, 1300 lipid bonds were 3–169 Å long, and La's md minimisation
  ran for 56 min dragging atoms across the box before it was cancelled (30478619, 30478620,
  30481824 cancelled). `insert_ion.py` now makes every molecule whole by walking its bonds to the
  nearest periodic image (longest bond 169.3 → 1.85 Å, and it refuses to write anything > 3 Å).
  Resubmitted: **30486649** La (viz), **30486651** Ca and **30486652** Lu (QOS short on gpu_viz).
- Also: the g-xTB cluster test (30423347) hit its 12 h limit before the Lu frames; not rerun,
  since g-xTB is set aside (C17).

**F11 – Three faults in the umbrella stage on a real membrane system.** (2026-10-24, fixed in
04335bc, 3d9c3c0, 4f61dda; 41 umbrella tests, 178 local tests pass)
- After F10, preprocess and md ran (minimisation 2.4 min on an L40, area per lipid 59.7 Å²), and
  every pull failed (30486649/51/52; then 30523099–101):
  1. `reference: "chain A and name CA"` (as in mcube's own example config) is not mdtraj syntax,
     which knows only `chainid N`. `reference_atoms` now translates `chain <letter>` using the
     topology's chain IDs and raises on an unknown letter.
  2. The pull read md's imaged `lastframe` with `PDBFile`: no lipid bonds and split molecules
     (F10 again, now inside mcube). New `load_structure` takes the bonds from preprocess's
     `processed_pdb` (same atoms, same order, checked) and makes molecules whole; the pull and
     both window types use it.
  3. `image_ion` refused a reference group spanning more than half the box, and `displacement`
     unwrapped the group about its first atom. The barrel's Cα span 55 Å of an 82.7 Å box, so the
     check fired, and the unwrap would have folded the far half back and put the centre in the
     wrong place. The group is now made whole by stepping atom to atom (consecutive Cα are
     3.8 Å apart), and "split" means two consecutive atoms more than half a box apart.
- Pulls running (2026-10-24): **30523774** La (viz), **30523775** Ca, **30523776** Lu (short).
  The ions moved during md: CV at start La −23.6 Å, Ca +27.9 Å and Lu +32.4 Å (the last two
  crossed the periodic boundary into the extracellular bulk; it is one water layer). Each pull
  covers about 77 Å at 10 Å/ns.
- **Pulls done (2026-10-24, ~2.5 h each on an L40, ~70 ns/day):** 71 seeds each, −35…+35 Å.
  **Windows started** on the three usable GPUs (gpu_cuda's estimate was 30 Oct): **30559113** La
  (viz, windows 10–64), **30559114** Ca and **30559115** Lu (short, windows 20–49), centre first
  so the constriction comes in first. Window length cut to **1 ns** (+0.1 ns shell-restrained
  equilibration; ~20 min each) for the quick test: 213 windows × 2 ns would be ~145 GPU-hours,
  about two days on three L40s. Windows near the constriction can be extended. (mcube rewrites
  config.yaml without comments, so edits must match the bare `time: 2000` line.)
- **Lu OpenMM-on-prmtop (30430921, 2 ns):** 2.33 Å, CN 7.65, with sander 7.64 and mcube 7.6.
  Three engines agree; F4 is closed on our side.
- **Pitfall:** `sbatch --export=ALL,OPTS="... --redo LA,LU"` splits at the comma, so only `LA`
  was redone and `LU_aqua` was reused from the cancelled def2 run. Rerun with
  `--redo LU_aqua`: **30559102**. Never put commas in an `--export` value.
- **Carboxylate check, Lu:** every Lu species' r2SCAN SCF failed to converge (energies off by
  10⁵ kcal/mol; the geometries are sane: Lu–O 2.32 Å in Lu(H₂O)₈³⁺). The script now never uses an
  unconverged energy, keeps the level shift for the final energy, falls back to Newton, and
  `--redo LU` recomputes only Lu: **30523294**.
- **GPU benchmark (30469103) failed:** gpu4pyscf could not load `libcusolver.so.11` (no CUDA 12
  runtime libraries in the job). To fix: load a CUDA 12 module in the job, or add the
  `nvidia-*-cu12` wheels (dry run first).

**F12 – The quick-test windows do not overlap, and Lu³⁺ windows remember the pull.**
(2026-10-25, provisional WHAM on windows 20–49 of Ca and Lu)
- **Overlap:** k = 10 kcal/mol/Å² (E = ½kΔz²) gives σ = √(kT/k) = 0.25 Å at 310 K, measured
  0.19–0.31 Å, with windows 1 Å = 4σ apart. WHAM flags overlap < 0.05 after 17 (Ca) and 16 (Lu) of
  30 windows (minimum 0.001 and 0.0). The free energy across each gap rests on tails. Fix:
  k = 2.5 (σ 0.5 Å, neighbours 2σ apart) on the same 1 Å grid.
- **Lu³⁺:** the provisional profile climbs 132 kcal/mol from +14 to −16 Å, which is not physical.
  The ion holds one or two protein O (cn_protein 1.0–2.0) from z −11 to +13, and in every window that
  count equals the value it started with: the coordination is the pull's (pulled +8 → +35 → −35 Å
  at 10 Å/ns, dragging a carboxylate along), not sampled. At z0 −11, −10 and −9 the ion sits at
  −8.9, −7.9 and −7.3 Å: the bound side chain holds it 2–3 Å off its restraint, and the
  −12/−11 pair is 2.9 Å apart with no overlap. This is the failure Zhang/Merz/Hu reported for Zn²⁺
  in ZIP (c4-tuning.md). It is made worse by the default C4, which over-binds Ln³⁺ to carboxylates (C20).
- **Ca²⁺:** windows track z0 (±0.1–0.3 Å) and the ion stays fully hydrated (7.3 waters) through
  the lumen; it touches protein O in a few windows only (−7, +4 Å, and +9…+15 Å, the last from the pull
  start). Provisional profile: ~7–12 kcal/mol across the barrel, bounded by poor overlap.
  Neither profile is a result yet.
- **Done (2026-10-25):** sections `umbrella_soft` (default C4, k 2.5; Ca, La, Lu) and
  `umbrella_chbe` (Kantakevičius set, k 2.5; Ca, La) added with `docs/channel/system/add_section.py`;
  they seed from the k = 10 windows' last frames at the same z0, so the default-vs-Ch-BE comparison uses
  one protocol. `umbrella_windows.slurm` takes `SECTION=`. Queued on gpu_cuda, windows 20–49 in
  blocks of 15 (8 h): Ch-BE **30576902–05** (CA 20, CA 35, LA 20, LA 35), soft **30576906–09** (same
  order). Lu soft is not queued yet.
- **La³⁺ (provisional WHAM, 33 windows, −25…+8 Å):** the same failure. The profile climbs ~150 kcal/mol,
  overlap is below 0.05 after 13 windows (min 0.0), the ion holds 1–2 protein O through the lumen with
  cn_protein equal to its start value in nearly every window, and it lags z0 by up to 1 Å. The soft and
  Ch-BE La sections seed from these frames, so they start with the pulled contacts. Ch-BE's weaker
  La–O C4 (9.5) may let them go; default-soft probably will not.
- **Still open for Lu:** soft k alone will not undo the pull memory. Options: seed windows from both
  pull directions and compare, run longer windows near Glu300, or add a coordination-number CV / metadynamics (as Merz/Hu did).

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
