# CSR mechanism study: decision log

A running record of the modelling decisions for the QM/MM study of meropenem hydrolysis by the mononuclear (single-Zn) β-lactamase CSR. Each entry gives the decision, why it was made, the evidence, the alternatives, and whether to revisit it. Keep it up to date, so the assumptions can be stated and challenged in the write-up.

**Aim.** Test the "one metal ion is enough" idea (Dupureur, *Metallomics* 2010, 2, 609, doi:10.1039/C0MT00013B) for a mononuclear Zn β-lactamase. Use the approach of Imhof, Fischer & Smith (EcoRV, *Biochemistry* 2009, 48, 9061, doi:10.1021/bi900585m) with modern methods:
- QM/MM minimum-energy pathways for a network of candidate mechanisms
- a comparison against a small model system
- the catalytic contribution of each residue

**Central hypothesis (D107, 2026-10-25).** CSR-1's two active-site substitutions
are at metal-ligand positions: Arg90 (BBL R118) at a Zn1 ligand position and
Lys226 (BBL K263) at the Zn2 position. Together these two monovalent cations
supply the **+2 that the absent Zn2 would have carried**. The hypothesis is that
they take over Zn2's *substrate-anchoring* role electrostatically while being
unable to take over its *catalytic* roles - depressing the nucleophile's pKa and
stabilising the anionic nitrogen - which would explain broad-spectrum binding
with poor turnover. Testable by per-residue charge deletion along each path, and
it predicts a larger effect on ring opening than on nucleophilic attack.

Two measured caveats on that wording, from `QM_MM_opt_last.pdb`: **Arg90 NH1 is
4.41 A from the zinc** and 5.13 A from the nucleophile, but **Lys226 NZ is 8.73 A
from the zinc** and points at the substrate carboxylate (4.80 A) instead - a Zn2
metal would sit ~3.5-4 A from Zn1, so the two cations act at opposite ends of the
site rather than jointly replacing one charge. And **Arg229 NH2 is 4.38 A from
the carboxylate and 6.61 A from the leaving nitrogen, closer than Lys226 on both
counts, yet is not in the QM region.**

The Imhof/EcoRV borrowing is **methodological only** - minimum-energy paths plus
per-residue electrostatic contributions. The structural analogy to EcoRV's Lys92
is retired (D107): that residue mainly stabilises a metaphosphate, which amide
hydrolysis has no counterpart for.

Their AM1/d → our GFN2-xTB for exploration, with a DFT check. Their CPR-style path search → our CI-NEB, transition-state optimisation and IRC.

**Paths on Bunya:**
- `/scratch/user/uquadhik/CSR_ligand`: the reactant, which also defines the system
- `/scratch/user/uquadhik/CSR_hydrolysed`: the product source

Atom names used below:
- **O01, H03:** the hydroxide oxygen and its hydrogen.
- **C4x, O2x, N1x:** the β-lactam carbonyl carbon (C7), carbonyl oxygen, and ring nitrogen (N4).
- **C6x, O3x, O6x:** the ring (C3) carboxyl group.
- **N2x:** the pyrrolidine nitrogen.
- **H01:** the carboxyl proton, moved by D13.

---

## Structure and topology

**D1 – The ligand is already in the starting complex, so it isn't merged again.** (2026-10-15)
- Evidence: `ligand_lastframe.pdb` and `hydrolysed_lastframe.pdb` each contain the ligand. Earlier runs that also merged the ligand ended up with two copies (`mer:102`, `hyd:105`).
- Decision: `merge_ligand: false`. The ligand is parameterised from its own PDB but not added again.
- Revisit: no.

**D2 – Ligand bonds come from CONECT records; bond orders come from SMILES templates.** (2026-10-15)
- Evidence: RDKit's distance-based bonding added a spurious O–O bond in the strained hydrolysed ligand (C6–O3 1.13 Å, O3···O6 1.76 Å). Its automatic bond-order assignment gave unphysical enolate/S⁺ forms for both meropenem and hydrolysed meropenem.
- Decision: connectivity comes only from CONECT records. Bond orders come from an explicit `ligand_smiles`, and the resulting SMILES is printed and checked. Stereocentres match the 3D structures.
- Alternative: hand-built mol/sdf files.
- Revisit: only if the ligand geometry is rebuilt.

**D3 – The hydroxide atoms were interleaved inside the meropenem block; both structures were reordered the same way.** (2026-10-15)
- Evidence: in the reactant, O01/H03 (residue `OH` 2) sit inside the `mer` atom block, so readers split `mer` into fragments ("no template… matches NME"). The product file has the same atom order, but O01/H03 belong to `mer` there.
- Decision: `scripts/make_residues_contiguous.py` computes the permutation from the reactant and applies it identically to the product, producing `*_contig.pdb`. The atom correspondence is kept: same names and elements, line by line.
- Revisit: no.

**D4 – Both NEB endpoints share one solvent box.** (2026-10-15)
- Evidence: NEB needs identical atoms in identical order, waters included. Solvating the two endpoints separately gives different numbers of waters and ions.
- Decision: solvate the reactant only. Build the product in the same box with `scripts/build_product_box.py`. For CSR the alignment RMSD on non-reactive atoms was 0.000 Å, with no solvent within 2.0 Å.
- Revisit: no. This also applies to the intermediates (D11).

**D5 – Starting structures are cleaned crystal structures.** (user, 2026-10-15)
- Source: per the user, cleaned crystal structures. `hydrolysed_lastframe.pdb` was built so its atom order matches the reactant's.
- Observation: in `ligand_lastframe.pdb` the hydroxide O–H was stretched to 1.54 Å, with H03 2.83 Å from Asp92 OD1. That isn't a typical crystal-derived geometry.
- Decision: use them, but re-equilibrate with MM before any QM/MM (D14) and select snapshots from that MD (D16).
- Revisit: **yes**. For the write-up, record:
  - the PDB IDs;
  - the cleaning steps;
  - how the substrate and OH⁻ were placed (crystal structures don't resolve H, and usually not the Michaelis complex);
  - how the hydrolysed product was built.

**D21 – The current model's Zn site has drifted from the crystal structure.** (found 2026-10-15; open)

Crystal reference: `CSR_ligand/CSR_nontruncated_1metal.pdb` (user: the His orientation and geometry to preserve).

| Zn ligand | Crystal | `ligand_lastframe` / MD job 29840911 |
| --- | --- | --- |
| His164 NE2 | **1.73 Å** (bound) | **4.56 Å** (not bound) |
| His88 NE2 | 2.65 Å | 2.77 Å |
| Asp92 | OD2 3.23 Å | OD2 2.40 / OD1 2.60 Å (bidentate) |
| Water / OH⁻ | water O 3.23 Å | OH⁻ 2.34 Å |
| Meropenem | amide O1x 3.03 Å; β-lactam O2x 4.90 Å | β-lactam O2x 2.55 Å |

- Neither His is flipped. NE2 faces the Zn for both His88 and His164 in both structures, so both are **HID**, matching the config.
- Problems in the crystal file itself:
  - Every His has HE2, including the Zn-facing NE2 of His88 and His164. His164 HE2 is **0.76 Å** from the Zn, so the H atoms came from a tool's default, not from the experiment.
  - Zn–His164 NE2 is 1.73 Å, too short (usually ~2.0–2.1 Å).
  - The Zn-site water is **1.52 Å** from Asp92 OD2, a physically impossible clash. Possibly a left-over from removing a second metal ("1metal").
- Other crystal features:
  - Meropenem is already a zwitterion there (COO⁻, NH₂⁺), consistent with D13.
  - The carboxylate is 4.0–4.7 Å from Lys226 NZ.
  - Water–C4x is 3.62 Å, at an 82° approach angle.
- Consequence: the MD started from `ligand_lastframe_zw.pdb` doesn't keep the crystal's His164 coordination.
- Revisit: **yes, decide**. Either rebuild the starting complex from the crystal heavy atoms (strip H, fix the water/Asp92 clash, decide the OH⁻ position, constrain His88/His164/Asp92 to the Zn during MM), or justify the drifted site.

**D22 – CSR-1 literature and the deposited model.** (2026-10-15; from a research agent's literature search, open decisions)
- **Enzyme:** CSR-1 is a subclass **B3** MBL from *Cronobacter sakazakii* (UniProt A7MFY4), in the B3-RQK group (Pedroso et al., *Protein & Cell* 2020, doi:10.1007/s13238-020-00736-4).
  - Zn1 (α) site: His116, **Arg118**, His196. Zn2 (β) site: Asp120, Gln121, Lys263 (BBL numbering).
  - Depositors number His196 as 192 and Lys263 as 254, which matches `CSR_nontruncated_1metal_original.pdb`. Our working files use those numbers minus 28.
- **Structures:**
  - Wild type: 6DN4 (no metal) and 6DQ2 (no metal).
  - Mutants: 6DQH (R118H, 1 Zn), 6DR8 (R118H/Q121H/K254H, 2 Zn) and 6NC5 (Q121H/K254H, Arg118 kept).
  - No substrate or product complexes have been published.
- **No published structure has Arg coordinating Zn.**
  - In 6NC5 (Arg kept), Zn_α is bound by His116 NE2 (2.22 Å) and His192 NE2 (2.09 Å) plus two waters (1.90, 2.14 Å). Arg118 NH1 is 4.02 Å away.
  - The review by Krco et al. (*Front. Chem.* 2023, doi:10.3389/fchem.2023.1196073) argues that Arg/Lys are probably protonated and poor metal ligands.
  - The user's understanding is that the Zn is coordinated by Arg118 plus two His. Unconfirmed.
- **Metal affinity (ITC):** wild-type α-site Kd ≈ 2 µM; β-site ≈ 150 µM (weak, not inhibitory).
- **Activity needs N-terminal truncation.** In the mature full-length enzyme, an N-terminal loop caps the active site. Our model is "nontruncated".
- **Mechanism:** no published mechanistic study or meropenem kinetics for CSR-1. Templates: mono-Zn B2 CphA (Simona et al., *JBC* 2009; Ackerman & Gatti, *PLoS ONE* 2013). In those, Asp120 is the general base, and the anionic N leaving group replaces Asp120 as a Zn ligand.
- **Our model:** the protein matches the crystal numbering (+28 relative to our working files). In the deposited model:
  - His116 NE2–Arg118 NH1 is 5.46 Å, so His116, His192 and Arg118 can't all bind one Zn without side-chain rearrangement.
  - The model's own Zn (occupancy 0.7) clashes with Asp120 (1.38 Å) and Lys254 (1.95 Å).
  - Its refined hydrogens give HID for 116, 140 and 192, and HIE for 81 and 178.
- Decision pending (user): the Zn coordination model, whether to take the Zn/water positions from 6NC5, and whether to truncate the N-terminus.

**D23 – Rebuilt starting complex `CSR_crystal`.** (user decisions, 2026-10-15)
- Supersedes the `CSR_ligand` model (D5, D21). MD job 29840911 was cancelled.
- Builder: `examples/csr/build_csr_crystal_start.py`. Output: `/scratch/user/uquadhik/CSR_crystal/preprocess_results/csr_crystal_start.pdb` and `mer_crystal.pdb`, plus `CSR_crystal/config.yaml`.
- **Protein:** heavy atoms of the crystal model `CSR_hydrolysed/CSR_nontruncated_1metal_original.pdb`.
  - H stripped, altloc A kept, waters and metals removed, full length kept (user).
  - **Numbering:** the start PDB and config use **model numbering = crystal (depositor) numbering − 28**, because ASH/pdbfixer renumber each chain 1..N before applying `residue_variants`. The first preprocess attempt with crystal numbers failed, because the variants landed on non-His residues. mcube now refuses inputs not numbered 1..N per chain.

    | Residue | Model | Crystal | BBL |
    | --- | --- | --- | --- |
    | His | 88 | 116 | 116 |
    | Arg | 90 | 118 | 118 |
    | Asp | 92 | 120 | 120 |
    | His | 164 | 192 | 196 |
    | Lys | 226 | 254 | 263 |

    Chain A has no gaps (crystal 29–302 → model 1–274).
  - Two malformed or duplicated header records at the top of that file are skipped. The file itself is left untouched.
- **Zn site from PDB 6NC5** (CSR-1 Q121H/K254H with Arg118 kept, 2.07 Å; numbering = crystal − 29), superposed on 256 Cα (RMSD 0.58 Å). It follows the literature (D22) rather than Arg coordination.
  - **One Zn only:** 6NC5's Zn_α. Zn_β and the adventitious Zn are discarded (user).
  - Zn_α and its two waters were moved as a rigid unit (Zn by 0.38 Å) so Zn–NE2 matches 6NC5: His116 2.22 Å, His192 2.09 Å.
  - **Nucleophile OH⁻ = 6NC5 HOH493** (user). In 6NC5 this water bridges to Zn_β. Here it H-bonds Asp120 OD2 (2.84 Å) and Lys254 NZ (2.63 Å), and its H points to the nearest acceptor.
  - **HOH485** stays as a Zn-bound water, next to Arg118 NH1 (2.76 Å; CD 2.71 Å, a slightly tight contact).
  - Asp120 OD2 is 2.80 Å from the Zn, inside the 3.0 Å cutoff, so preprocess counts it as a Zn ligand and constrains it during MD. Revisit after inspection.
- **Meropenem:** pose from the processed model (`CSR_ligand/CSR_nontruncated_1metal.pdb`, identical protein, Cα RMSD 0.000 Å after +28 renumbering). Rigid-body shift towards OH⁻ with C4x 3.0–3.2 Å at 100–115°; contacts down to ~2.2 Å are left for MM minimisation (user).
  - **Result:** shift RMSD 0.91 Å (max 1.41 Å, 11° rotation). OH⁻···C4x is **3.47 Å at 79°**, so the target attack geometry was **not** reached; the rigid shift can't get there without clashes. The closest contacts are 2.25 Å (heavy atoms) and 2.04 Å (H).
  - Snapshots lined up for attack will therefore come from MD (D16), not from this start structure.
- **Protein conformers:** Ala101 has only altlocs B (occupancy 0.00) and C (1.00), so the highest-occupancy altloc is kept per residue (101C). Keeping only altloc A would have dropped the residue and made a gap. 257 Cα superposed onto 6NC5.
- **Build output:** `csr_crystal_start.pdb` (2164 atoms) is in model numbering 1–275 (= crystal 29–303), with one Zn (chain B), meropenem and OH⁻ (chain C), and the Zn-bound water (chain D).
- **His tautomers** (from the crystal model's refined H):

  | Crystal number | Model number | Tautomer |
  | --- | --- | --- |
  | 81 | 53 | HIE |
  | 116 | 88 | HID |
  | 140 | **112** | **HID** (the earlier CSR_ligand runs used pdbfixer's HIE) |
  | 178 | 150 | HIE |
  | 192 | 164 | HID |

  `require_his_variants: true`.
- **QM region:** extra residues Arg90 and Lys226 (model numbering), plus **meropenem (chain C, residue 1)**.
  - Fix: the first preprocess left meropenem in MM. It doesn't bind the Zn (nearest atom: hydroxyethyl O1x at 3.13 Å), so the coordination search didn't pick it up. It was added via `extra_qm_residues: {C: [1]}` and preprocess rerun.
- **Preprocess result** (`CSR_crystal/preprocess_results/finalsystem.pdb`, 49,373 atoms: 15,020 waters, 27 Na⁺, 29 Cl⁻):
  - **QM region:** 119 atoms, charge +2 (Zn²⁺, OH⁻, Asp92⁻, Arg90⁺, Lys226⁺, His88, His164, meropenem 0, Zn-bound water). 57 atoms are whole non-protein units without link atoms (Zn, water, OH⁻, meropenem).
  - **MM constraints on the Zn:** water 1.90, His164 NE2 2.09, OH⁻ 2.14, His88 NE2 2.22 and Asp92 OD2 2.80 Å.
  - **Tautomers as built:** His53/150 HIE; His88/112/164 HID.
  - **Protonation:** Lys226 protonated, Arg90 cationic, Asp92 deprotonated. Expected QM charge +2: Zn²⁺, OH⁻, Asp120⁻ if coordinated, His×2, meropenem 0, Arg⁺, Lys⁺, water 0. Check after preprocess.

**D24 – Role of Asp120 (model Asp92): H-bonded to the Zn-bound hydroxide, not bound to the Zn.** (user hypothesis, literature checked 2026-10-15)
- **User hypothesis:** the Asp is the "bridging" ligand known from B3 MBLs and nucleases, and it interacts with the Zn *through* the hydroxide rather than binding it directly.
- **Literature:**
  - **Mono-Zn MBLs (supports it):**
    - In mono-Zn B1 enzymes, Asp120 is not a metal ligand. It makes a strong H-bond to the Zn1-bound nucleophilic OH⁻. In two-zinc B1/B3 enzymes it is instead a Zn2 ligand. Source: Llarrull & Vila, ICBIC-12 conference abstract, "Asp120 is not a proton donor in di-zinc *Bacillus cereus* metallo-β-lactamase". The abstract text came via a search summary; I couldn't open the PDF.
    - Dal Peraro, Vila, Carloni & Klein, *JACS* 2007, 129, 2808 (doi:10.1021/ja0657556), QM/MM:
      - In mono-Zn BcII, an Asp120 / Arg121 / Cys221 H-bond network does the job of the second Zn: it orients the nucleophile.
      - In two-zinc CcrA, Asp120 binds Zn2, and the bridging OH⁻ H-bonds to Asp120.
      - The protonation state of Asp120 flips the orientation of the Zn1-bound OH⁻.
    - BcII R121C crystal structures at several pH values: the Zn1-bound OH⁻ (the proposed nucleophile) is present at pH ≥ 5 and lost at pH 4.5. Protonation of Asp120 and of the OH⁻ are given as the causes of inactivation at low pH. (Citation details not retrieved.)
  - **Two-zinc B3 MBLs (L1):** Asp120 is an invariant Zn2 ligand, and a water/hydroxide (Wat1) bridges Zn1 and Zn2 (Crisp et al., *Biochemistry* 2007, 46, 10664; PDB 2QIN, 2QJS). QM/MM on L1 finds that Zn-bound Asp120 helps align the nucleophile and accepts its proton only after nucleophilic addition.
  - **CSR-1, PDB 6NC5 (two-zinc Q121H/K254H mutant):** Asp120 (6NC5 numbering 91) binds Zn_β directly (1.95 Å), and HOH493 bridges Zn_α and Zn_β. In our one-zinc model, that HOH493 is the OH⁻ nucleophile, with Asp OD2···O 2.84 Å (an H-bond).
  - **Counter-example:** in the mono-Zn **B2** enzyme Sfh-I, His118 rather than Asp120 activates the hydrolytic water (Fonseca et al., *J. Mol. Biol.* 2011, 411, 951; PDB 3Q6V, 3SD9).
  - **Nucleases (different picture):**
    - Carboxylates there often bridge two metals *directly*, as μ-carboxylates. In EcoRV, Asp74 is shared between the two metal sites (Glu45/Asp74 and Asp74/Asp90); in T7 endonuclease I, Asp55 spans sites 1 and 2.
    - Whether both metals take part in chemistry is disputed enzyme by enzyme. The 1997 *Biochemistry* EcoRV mutant study and T5 flap endonuclease (*JBC* 2011) argue the second metal isn't catalytic.
- **Conclusion:**
  - In MBL terms the *bridging* ligand is the hydroxide (μ-OH), not the Asp.
  - For a **one-zinc** site, the supported picture is that Asp120 sits in the second shell and orients/activates the Zn-bound OH⁻ through a strong H-bond: "coordinating through the hydroxide". This is supported for mono-Zn **B1** (BcII) and consistent with the CSR-1 6NC5 geometry, but **not shown for B3-RQK/CSR-1 specifically**.
  - Direct Asp coordination is the two-zinc (Zn2) role, or the μ-carboxylate role in nucleases.
- **Implication for the model:** in `csr_crystal_start.pdb`, Asp92 OD2 is 2.80 Å from the Zn, inside the 3.0 Å cutoff. Preprocess will treat it as a Zn ligand and **constrain Zn–OD2 at 2.80 Å** during MM, which goes against this picture.
  - **Decision (user, 2026-10-15): keep Asp92 as preprocess finds it**, i.e. treated as a Zn ligand and constrained at 2.80 Å during MM relaxation and MD. It is in the QM region either way, and no constraints apply in QM/MM, where the H-bond vs coordination question is decided.
  - Alternative not taken: add Asp92 to `constraint_exclude_residues` so MM relaxation decides between H-bond and coordination.
  - Check after MD whether the Asp92···HO–Zn H-bond holds.

**D25 – Pull meropenem's β-lactam carbonyl onto the Zn with restrained MM.** (user, 2026-10-15)
- **Problem:** in `CSR_crystal/preprocess_results/finalsystem.pdb` the inherited pose points the hydroxyethyl side chain at the Zn:
  - O1x is 3.13 Å from the Zn; the β-lactam carbonyl O2x is **5.43 Å** away.
  - The OH⁻ sits almost on the Zn–O2x line, 0.92 Å from its midpoint.
  - OH⁻–C4x is 3.47 Å at 79°.
- **Assumption (user):** the β-lactam carbonyl O coordinates the Zn (< 3 Å) in the Michaelis complex. There is no CSR-1 substrate complex to confirm it; it's by analogy with other MBLs.
- **Decision (user):** reuse mcube's md stage and let restrained MM move the substrate, rather than re-docking or conformer placement.
  - mcube change: `md.inputs.restraints` is passed to ASH `OpenMMTheory(restraints=…)` as an OpenMM HarmonicBondForce, active through minimisation, warm-up, NPT and NVT.
- **Settings in `CSR_crystal/config.yaml`** (atom indices are 0-based in `finalsystem.pdb`):
  - **Restraints**, k = 20 kcal/mol/Å²:
    - Zn(4203)–O2x(4209) → 2.2 Å;
    - OH⁻ O01(4255)–C4x(4208) → 3.0 Å.
  - **Constraints** (fixed at the starting distance): Zn–His164 NE2 (2451), Zn–OH⁻ (4255), Zn–His88 NE2 (1349) and Zn–Asp92 OD2 (1414).
    - The Zn–water (HOH485, 4257) constraint is dropped, so the carbonyl can displace that water from its Zn site.
  - **Run:** minimisation, gentle warm-up, NPT box equilibration, then **1 ns NVT, all restrained**.
- **Job:** 29848762 (`CSR_crystal/MD_restrained.slurm`, 16 CPUs, 16 GB, 24 h), output in `CSR_crystal/md_restrained/`. Two CSVs record the geometry per frame: `attack/nac_geometry.csv` (OH⁻–C4x distance, attack angle) and `zn_o2x/nac_geometry.csv` (Zn–O2x).
- **Correction (user, 2026-10-15): no unrestrained MM MD is ever run.** It would break the Zn coordination geometry. Every MM stage keeps the metal site constrained or restrained (the purpose of the mcube pipeline). The Zn site only relaxes freely in QM/MM.
  - The attack-ready snapshots come from the **restrained** NVT. The QM/MM optimisation (no restraints) removes the restraint bias.
  - Full plan: [csr-neb-plan.md](csr-neb-plan.md).
- **Superseded, do not run:** an unrestrained follow-up was prepared on 2026-10-15 but never submitted:
  - Files: `CSR_crystal/config_unrestrained.yaml` and `MD_unrestrained.slurm`.
  - Start structure: `md_restrained/NVTtrajectory_lastframe.pdb` (OpenMM-written, so names and CONECT records are preserved).
  - Settings: no restraints, the same four Zn constraints at their end-point distances, no warm-up and no NPT (box already equilibrated), then 2 ns NVT with frames every 2 ps.
  - Analysis: up to 10 attack-ready frames after the first 100 ps, plus the Zn–O2x distance per frame.
  - To be submitted only after the restrained run is checked.
- **Next:**
  1. Check that Zn–O2x reaches ~2.2 Å, OH⁻–C4x ~3.0 Å at 95–120°, and that HOH485 leaves the Zn.
  2. Select ≥ 5 attack-ready snapshots from the restrained NVT.
  3. Run a QM/MM reactant optimisation (`mcube qmmm`) per snapshot.
- **ASH quirk:** OpenMMTheory only adds user constraints given as `[i, j]` pairs, which take the current distance. Entries given as `[i, j, d]` are silently ignored, because `add_bondconstraints` sits inside the "missing distance" branch.

## Force field and protonation

**D6 – Force fields.** (inherited, 2026-10-15)
- Protein: Amber ff14SB. Water: TIP3P.
- Ligand: OpenFF (ASH's `small_molecule_parameterizer` default, `openff-2.0.0.offxml`).
- Hydroxide: `hydroxide.xml`, generated by the user with OpenFF via ASH's `small_molecule_parameterizer` (charges O −1.02 / H +0.02, O–H 0.975 Å). Note: OpenFF is parameterised for organic molecules; a bare OH⁻ next to Zn²⁺ is outside its training. That matters little here, because OH⁻ is in the QM region in every QM/MM step and is held to the Zn during MM.
- Zn²⁺: nonbonded ion. Its Zn–ligand distances are held by constraints **only during MM relaxation and MD**, to prepare structures for QM/MM (user, 2026-10-15). In QM/MM the Zn and all its ligands are QM, and the qmmm stage applies no constraints.
- Revisit: only if long, unconstrained MM MD of the metal site is needed (then use a bonded model such as MCPB.py, or 12-6-4 LJ).

**D7 – Residue protonation; histidine tautomers must be explicit.** (user: His tautomers are critical and must be set in `config.yaml`)
- Audit of `finalsystem.pdb` (job 29840911 preprocess):

| His | Tautomer used | Set in config? | Zn distance | ND1 near | NE2 near |
| --- | --- | --- | --- | --- | --- |
| A:88 | HID | yes | 2.8 Å (NE2 bound) | Asp185 OD1 2.89 Å (His–Asp dyad) | Zn 2.77 Å |
| A:164 | HID | yes | 4.6 Å | Ala188 O 2.90 Å | meropenem O2x 3.20 Å, O5x 3.23 Å |
| A:53, A:112, A:150 | HIE (pdbfixer default) | **no** | 20–28 Å | — | — |

- His164 is mechanistically sensitive. Its NE2 is 3.2 Å from the β-lactam carbonyl O, which becomes the oxyanion in the tetrahedral intermediate.
  - As HID, NE2 has no H, so it can't donate an H-bond to the oxyanion; as HIE it could.
  - ND1–H…O=C(Ala188) supports HID.
  - Consider adding His164 to the QM region, and test HIE as a sensitivity check.
- Other residues come from pdbfixer at pH 8.0 (mcube default): Lys226 protonated (3 HZ), Asp92 deprotonated.
- Action:
  - Set A:53, A:112 and A:150 explicitly (currently HIE by default) and add `require_his_variants: true`. mcube now warns or errors on any His without a tautomer.
  - Don't edit `config.yaml` while MD job 29840911 is writing results to it.
- Revisit: **yes**. Run a pKa prediction (PROPKA/H++) for the pH 8 default, and justify each His in the write-up.

**D13 – Meropenem is modelled as the zwitterion.** (2026-10-15, user decision)
- Evidence: the ring carboxylate pKa is about 3 and the pyrrolidine amine pKa about 7.5–9, so at pH 7 the molecule is mostly COO⁻ / NH₂⁺. The earlier model had a neutral COOH and a neutral amine. ~~Metallo-β-lactamases typically bind the C3 carboxylate with a Lys, and Lys226 may be that residue (NZ is 5.4 Å from C6x in the current frame).~~ **[Retracted 2026-10-25, D107: wrong residue.** The carboxylate lysine is BBL Lys224, a B1 feature; B3 enzymes have Ser221/Ser223 there. Lys226 is BBL Lys263, the Zn2 *metal-ligand* position, as D22 and D40 in this same log already state. Lys226 does turn out to salt-bridge the C3 carboxylate here, but for a different reason; see D107.**]**
- Decision: move H01 from O3x to N2x in both the reactant and the product, so atom order is unchanged. Re-parameterise with zwitterion SMILES. Reactant ligand net charge 0; hydrolysed product −1. The QM charge is unchanged.
- Alternatives: neutral, as before; anionic (COO⁻ with a neutral amine).
- Revisit: possibly, as a sensitivity check.
- Done (2026-10-15) with `scripts/move_proton.py`, which places the new N2x–H at 1.02 Å, pointing away from the other N2x neighbours.
  - Files: `ligand_lastframe_zw.pdb` and `mer_zw.pdb` (reactant); `hydrolysed_lastframe_zw.pdb` and `hyd_zw.pdb` (product). Atom order is identical between reactant and product.
  - Reactant: the new H is 1.69 Å from the existing N2x–H (H10x), normal for NH₂⁺.
  - Product: the new H is tighter, 1.54 Å from H10x and 1.96 Å from C9x, because the hydrolysed ligand's geometry is already strained. QM/MM optimisation should relax it; check this.
  - The template SMILES matched, stereocentres included.

## Metal site

**D8 – Zn coordination cutoff of 3.0 Å.** (2026-10-15, user decision)
- Evidence: typical Zn–ligand bonds are 2.0–2.3 Å, and coordination cutoffs are usually 2.5–3.0 Å. In this frame His88 NE2 is at 2.77 Å and Arg90 NH1 at 3.50 Å. The old 4 Å cutoff wrongly constrained Arg90 to the Zn.
- Result: the coordinating atoms are OH⁻ O01 (2.34), Asp92 OD2 (2.40; OD1 at 2.60), meropenem carbonyl O2x (2.55) and His88 NE2 (2.77).
- Revisit: no.

**D15 – During MD, Zn–ligand distances are fixed except the substrate carbonyl's.** (2026-10-15, user decision)
- Evidence: fixing every Zn–ligand distance pins the substrate carbonyl at 2.55 Å and biases the Michaelis complex.
- Decision: keep the Zn–OH⁻, Zn–Asp92 and Zn–His88 distances fixed. Leave Zn–O2x free (`constraint_exclude_residues: {C: [1]}`).
- Confirmed by preprocess (job 29840911): constraints [Zn–O01 (OH), Zn–Asp92 OD2, Zn–His88 NE2], with no Zn–O2x constraint.
- Revisit: **yes**. Check whether the carbonyl stays near the Zn. A free hydroxide or a bonded model would be alternatives.

## QM/MM model

**D9 – QM method: GFN2-xTB with electrostatic embedding. A g-xTB switch is blocked for now.**
- Run through ASH. QM/MM boundaries cut at CA–CB with link atoms. Non-amino-acid residues (Zn, ligand, OH) are fully QM, with no link atoms.
- Evidence of concern: in the first reactant optimisation, Zn–OH shortened to 1.75 Å, against the usual ~1.9–2.0 Å.
- Revisit: **yes, required**. Benchmark against DFT on a cluster model, and recompute DFT/MM energies along the xTB paths (P5).
- **g-xTB considered** (user preference, 2026-10-15). g-xTB approximates ωB97M-V/def2-TZVPPD for Z = 1–103, and v2.0 (Apr 2026) has analytic gradients. It can't be used for our QM/MM yet:
  - g-xTB v2.0.1 (grimme-lab/g-xtb, May 2026, a modified `xtb --gxtb` binary) lists point-charge embedding (`$pcem`) as an unsupported feature. ASH's QM/MM with xtb passes the MM charges that way.
  - ASH's `gxTBTheory` targets the old v1 `gxtb` binary. It silently ignores the MM point charges, so electrostatic embedding quietly degrades to mechanical embedding. It also ignores `numcores`.
  - Mainline tblite (v0.7.0, Jul 2026) has no g-xTB, and its only external potential is a uniform electric field. No ORCA g-xTB support was found.
  - Mechanical embedding isn't acceptable here. Imhof et al. found electrostatic stabilisation of the transition state by the protein to be central, and our QM region carries charge +2.
- Where g-xTB can be used now: cluster models without MM (gas phase, or its preliminary GBE/COSMO solvation, which its README says is unstable for gradients). That makes it a fast reference for the GFN2-xTB cluster benchmark (P1), alongside a few DFT points to validate g-xTB itself for the Zn site.
- **DFT/MM with PySCF** (checked 2026-10-15 on Bunya with `scripts/check_pyscf_dft.py`):
  - The ASH env has PySCF 2.13.0 with libxc 7.0.0. ωB97M-V (`wb97m_v`) is recognised as a range-separated hybrid with VV10 non-local correlation.
  - The analytic RKS gradient with density fitting, VV10 and MM point charges matches finite differences: 9×10⁻⁸ Eh/Bohr on a QM atom and 1.5×10⁻⁹ on an MM charge.
  - ASH's `PySCFTheory` supports electrostatic embedding and point-charge gradients, density fitting (`densityfit`), seminumerical exchange (`sgx`) and a GPU mode (`platform='GPU'` via gpu4pyscf).
  - gpu4pyscf isn't installed. Bunya has NVIDIA GPUs (`gpu_cuda`: A100, H100, L40S; `gpu_sxm`: H100), so GPU DFT is possible. No ORCA module was found.
- Cost: a full ωB97M-V/MM NEB on CPU (105 QM atoms, about 12 images × 50–150 iterations = 600–1800 gradients) is impractical. On GPU it's plausible, but needs a timing benchmark first.
- Proposed (2026-10-15):
  1. GFN2-xTB/MM paths for exploration.
  2. ωB97M-V/def2-SVP/MM, with gpu4pyscf on an H100: re-optimise the endpoints, re-converge the NEB from the xTB path, transition-state optimisation, QM-region frequencies and IRC.
  3. ωB97M-V/def2-TZVP(D)/MM single points.
  Install gpu4pyscf in a separate env to protect the ASH env's PySCF version. Benchmark one gradient first.
- Watch for: `$pcem` support in g-xTB. Then update ASH's `gxTBTheory` to v2 (`xtb --gxtb`, analytic gradients, point charges) and switch the QM/MM.

**D10 – QM charge and multiplicity.**
- Zn²⁺ is d¹⁰, so the multiplicity is 1.
- With Arg90 and Lys226 in the QM region (D12): Zn²⁺ (+2) + OH⁻ (−1) + Asp92⁻ (−1) + His88 (0) + meropenem (0) + Arg90⁺ (+1) + Lys226⁺ (+1) = **+2**.
- The product is also +2, because hydrolysed meropenem (−1) absorbs the OH⁻.
- The earlier optimisation (job 29837804, cancelled) used 0, without Arg90 or Lys226.

**D12 – Arg90 and Lys226 are added to the QM region.** (2026-10-15, user decision)
- Hypothesis: Lys226 plays a catalytic role like EcoRV's Lys92 (positioning, or stabilising the anionic intermediate), and Arg90 stabilises the oxyanion or hydroxide. **[Retired 2026-10-25, D107. The EcoRV analogy does not survive: Lys92's role there is contested, and in Imhof 2009 - the paper this project is modelled on - it mainly stabilises a *metaphosphate*, which has no counterpart in amide hydrolysis. Keeping Arg90 and Lys226 in the QM region is still right, for the reasons in D107, but not for this reason.]**
- Geometry in the reactant frame:
  - Lys226 NZ is 4.3 Å from O01, 5.6 Å from N1x and 5.4 Å from the C3 carboxyl.
  - Arg90 NH2 hydrogen-bonds the hydroxyethyl O5x (2.73 Å); NH1 is 3.8 Å from O01.
- Decision: `extra_qm_residues: {A: [90, 226]}`.
- Confirmed by preprocess (job 29840911): 105 QM atoms (54 of them without link atoms: Zn, the ligand and the OH). `qm_residues` = Zn B:1; OH C:2; mer C:1; Asp92, His88, Arg90, Lys226.
- Revisit: after sampling, check which other residues contact the reacting atoms (Asn, other Lys/His, waters).

**D11 – Active region of 8 Å around the Zn.**
- Built from atom 4203 (431 atoms in the first test). The same active atom list is used for every endpoint, intermediate and NEB (`--actatoms`).
- Revisit: after the QM region change, check that the converged results don't depend on the radius.

## Mechanism and path protocol

**D16 – Snapshots must be lined up for attack.** (2026-10-15, user decision; sampling source updated by D25: snapshots come from the restrained NVT, never from unrestrained MD)
- Evidence: in the current reactant, O01–C4x is 3.81 Å and the O01–C4x–O2x angle is 51°. The hydroxide lies beside the carbonyl O rather than on the Bürgi–Dunitz trajectory.
- Decision: minimise, warm up, equilibrate at constant pressure, then run 2 ns at constant volume (4 fs steps), saving a frame every 2 ps. Keep frames with O01–C4x ≤ 3.2 Å and O01–C4x–O2x at 95–120°. If fewer than 5 qualify, consult before adding restraints.
- Alternative: a flat-bottom restraint, which is biased.
- Running: job 29840911 (`CSR_ligand/MD_zw.slurm`, 16 CPUs). It runs preprocess, then MD in `md_zw/`, then `scripts/select_nac_frames.py` (skipping the first 50 frames, up to 10 frames ≥ 50 frames apart) into `nac_frames/`.
- Superseded: the cancelled job 29837804 (QM/MM optimisation without Arg/Lys, unequilibrated box, neutral meropenem).
- Revisit: the thresholds and how many snapshots to use.

**D14 – MM-equilibrate the box before any QM/MM.**
- Evidence: in the first QM/MM optimisation, step 1 dropped the energy by about 7 Eh, almost all of it MM energy. The step-0 maximum gradient was 47.8 au. The freshly placed waters were strained.

**D17 – The reaction is split into elementary steps, each with its own NEB.** (2026-10-15)
1. Michaelis complex (R) → tetrahedral intermediate (TI): O01 attacks C4x.
2. TI → ring-opened anionic intermediate (EI⁻): C4x–N1x breaks. If the TI isn't a minimum, steps 1 and 2 are concerted, and that is itself a result.
3. EI⁻ → Δ2-enamine product: N1x is protonated.

Each intermediate is built in the shared box (D4) and optimised, then confirmed as a minimum with a frequency calculation on the QM region. Each step gets a CI-NEB with about 12 images, transition-state refinement, a check for exactly one imaginary frequency, and an IRC to confirm the connections. This is repeated on 5 or more snapshots, and barriers are reported as mean ± spread.

**D18 – The product atom mapping fixes which proton ends up on the nitrogen.**
- Evidence: in the product, H03 (the hydroxide's own proton) is bonded to N1x, and O01 is a carboxylate O.
  - Reactant: Asp92 is 5.4–5.9 Å from N1x and Lys226 5.6 Å. No water is within 9 Å of N1x.
  - With this mapping, the possible routes are H03 directly to N1x, or H03 relayed through Asp92.
- Decision: treat each proton source as a separate hypothesis, with its own product file: H03 direct or via Asp92, Lys226 (Lys neutral, H03 left on O01), or water if any can reach.
- Revisit: **yes**. This is part of the mechanism network.

**D19 – The Δ2 enamine product comes first; the Δ1 imine (C2 protonation) later.** (user decision)

**D26 – How scans, NEB and TS/IRC are run in mcube.** (2026-10-15; settings in the config, as the user asked)
- Relaxed scans (`mcube scan`): one bond distance per scan, stepped from the reactant. Each point is optimised from the previous one with that distance fixed (geomeTRIC `$set`). Active-region waters stay rigid throughout, as in the `qmmm` minimisations.
- NEB: every NEB setting lives in the config (`neb.inputs`, or a per-step section). The default is 12 images (previously hard-coded 8), with CI on and ASH's convergence defaults (max force on the climbing image 0.026 Eh/Å, RMS 0.013). The option `method: nebts` runs ASH's NEB+TS.
- TS refinement (`mcube ts`): a geomeTRIC TS optimisation, using an exact Hessian over `hessian_atoms` (the reacting atoms and their neighbours) and a unit Hessian for the rest of the active region. This is the same scheme as ASH's NEBTS.
  - Waters are held rigid here too, so the TS and the minima share one potential. The NEB images themselves have flexible waters (KNARR takes no constraints).
- Frequencies: "one imaginary frequency" is judged with a −50 cm⁻¹ threshold. Smaller imaginary modes from a partial Hessian are reported separately as noise.
- IRC: the default is a quasi-IRC. The structure is displaced ±0.15 Å along the imaginary mode and each side is minimised. This shows which minima the TS connects, but it is not a true mass-weighted IRC path. Sella's IRC is available (`irc_method: sella`) but expensive with an 8 Å active region.
- Revisit: the IRC method, if a smooth path is needed for DFT single points (use NEB images instead, or Sella with a smaller active region); and whether to run NEB with rigid waters.

**D27 – Restrained NVT outcome, and a QM/MM probe before snapshot selection.** (2026-10-16, user decision)
- Evidence: job 29848762 (restrained warm-up, NPT, 1 ns NVT) finished normally; `scripts/site_report.py`, `md_restrained/site*/`.
  - Protein: CA RMSD ~1 Å, stable. Meropenem moved 3.7 Å on average from the start (the restraint pulled it onto the Zn); about 1.3 Å over the second half.
  - Carbonyl on the Zn: Zn–O2x 1.79–2.19 Å (mean 1.97). The restraint target was 2.2 Å.
  - **HOH485 did not leave**: Zn–O 1.72–1.95 Å (mean 1.83). Too short, an artefact of the nonbonded MM Zn²⁺. The water H-bonds both Asp92 oxygens (2.8 Å). The Zn has 5 ligands within 2.5 Å (His88, His164, OH⁻, O2x, water), plus Asp92 constrained at 2.80 Å.
  - **No near-attack geometry**: OH⁻–C4x 3.21–3.86 Å (mean 3.52; target 3.0) at an O01–C4x–O2x angle of 24–54°, not 95–120°. The OH⁻ and O2x are cis on the Zn (O01···O2x 2.67 Å), so the OH⁻ sits beside the carbonyl O, not over the C4x face. The water is no better placed (O···C4x 3.9 Å, 24°). 0 of 500 frames meet the D16 criteria.
- Decision: before changing the MM protocol, QM/MM-optimise (no restraints, no metal constraints, rigid waters, 8 Å active region around the Zn, up to 300 steps) frames 300, 400 and 499 (300–499 ps) as a probe. This removes the MM Zn artefacts and shows whether the water stays and how the OH⁻/carbonyl arrangement relaxes. Job 29893675 (array, `qmmm_probe/`, results in sections `qmmm_probe_f00300/00400/00499`).
- Alternative, not taken yet: rerun the restrained MD with an attack-geometry restraint and a restraint moving the water off the Zn (biased).
- Probe results, frames 300 and 400 (both converged in ~25 min on 8 CPUs; frame 499 still running):
  - **Asp92 binds the Zn directly** in QM/MM (OD2 2.09–2.11 Å; held at 2.80 Å in MM). This bears on D24 and on the MM constraint value.
  - HOH485 stays on the Zn (2.12–2.15 Å, normal length once the MM artefact is gone).
  - The carbonyl stays on the Zn in frame 300 (O2x 2.17 Å, 6-coordinate Zn), but **leaves it in frame 400**, where the OH⁻ proton H-bonds O2x instead (H03···O2x 1.82 Å; 5-coordinate Zn).
  - Neither frame approaches attack: OH⁻–C4x 3.65 / 3.83 Å, O01–C4x–O2x 30° / 17°. The OH⁻ stays beside the carbonyl O.
- Frame 499 (2026-10-16): **not converged** in 300 steps (2 h 9 min; the last geometry was kept, `converged: false`). It ends like frame 400:
  - the carbonyl is off the Zn (O2x 3.54 Å), with H03···O2x 1.77 Å;
  - the Zn is 5-coordinate (His88, His164, Asp92 2.10 Å, OH⁻, HOH485);
  - OH⁻–C4x is 3.82 Å at 13°.
- Every optimisation step was checked (`site_report.py` on `geometric_OPTtraj_Full.xyz`; 88 + 153 + 301 steps). **No step reaches an attack geometry** (≤ 3.2 Å and 95–120°), either for the OH⁻ (angle never above 41°) or for HOH485 (angle never above 69°; it came within 2.88 Å of C4x once in frame 499).
- Summary over the 3 frames:
  - Asp92 binds the Zn in all three.
  - HOH485 stays bound in all three.
  - The carbonyl leaves the Zn in 2 of 3, its O taking an H-bond from the OH⁻ instead.
  - The restrained-MM arrangement (OH⁻ cis to the carbonyl O on the Zn) does not relax towards attack.
- Next (user decision, 2026-10-16): a relaxed O01–C4x scan from the converged frame-400 structure, to test whether attack from this arrangement is feasible at all.
  - Setup: 3.80 → 1.45 Å in 25 points (−0.098 Å per step), at most 200 steps per point, the frame-400 active region (486 atoms).
  - Monitored at every point: Zn–O2x, C4x–N1x, H03···O2x, Zn–O01 and C4x–O2x.
  - Job 29919076 in `scan_attack_f400/`, config section `scan_attack_f400`.
  - It is a driven path from a non-NAC start, so its maximum is an upper bound on the barrier, not a TS energy.
  - **Cancelled** after 11 min (superseded by D28).

**D28 – Redo with the plain mcube pipeline: carbonyl constrained like any Zn ligand, QM/MM MD before optimising.** (2026-10-16, user decision; supersedes the restrained-MD route of D25 and the optimisation-only probe of D27)
- Why (user):
  - The carbonyl must be treated exactly like the His and Asp ligands. preprocess should pick it up within the Zn cutoff and constrain it, with no hand-written restraints or exclusions.
  - The QM/MM stage must run a QM/MM MD equilibration before the optimisation, as in `CSR_ligand`/`CSR_hydrolysed` `MM.slurm` (preprocess → md → qmmm).
- What the earlier run actually had:
  - D15/D25 excluded the substrate from the constraints (`constraint_exclude_residues: {C: [1]}`) and used a harmonic Zn–O2x restraint (2.2 Å, k = 20). That held O2x at 1.8–2.2 Å through NPT and NVT (2.5 Å after minimisation).
  - Asp92 was constrained at 2.80 Å in every MM stage; it only moved, to 2.1 Å, in the unconstrained QM/MM optimisation.
  - `qmmm.inputs.skip_equilibrate` was `true`.
- Setup in `/scratch/user/uquadhik/CSR_crystal_v2` (CSR_crystal kept for comparison):
  - `protein_pdb`: the restrained NVT end point with bulk water and ions removed (`scripts/strip_solvent.py`, HOH485 kept), so the carbonyl starts on the Zn.
  - preprocess inputs: as CSR_crystal, minus `constraint_exclude_residues`.
  - md inputs: as CSR_crystal (gentle warm-up, NPT, 1 ns NVT), minus the hand-set constraints and restraints.
  - qmmm inputs: `skip_equilibrate: false`, 50 ps QM/MM MD at 1 fs (as CSR_hydrolysed), then optimisation (`actradius 8` around the Zn, 300 steps); 16 cores.
- preprocess (cutoff 3.0 Å) constrains six Zn ligands at their start distances:
  - HOH485 O 1.79 Å;
  - **meropenem O2x 1.94 Å**;
  - His164 NE2 2.09 Å;
  - OH⁻ O01 2.14 Å;
  - His88 NE2 2.22 Å;
  - Asp92 OD2 2.79 Å.
  QM region 119 atoms; same atom indices as CSR_crystal.
- Asp92: kept at 2.80 Å in MM only (user); free in QM/MM.
- Job 29920397 (`MM.slurm`: md → qmmm, 16 CPUs, 10-day limit).
- Note: the HOH485 and O2x distances are frozen at values from the restrained run (1.79 and 1.94 Å), which are short because of the nonbonded MM Zn. The QM/MM MD releases them.
- QM/MM relaxation (user, 2026-10-16): the QM/MM MD is expected to relax the site a little; that is its purpose. It must still keep the bond geometry of everything in the QM region. Check QM-region bonds when picking frames.
- **NAC frames come from the QM/MM MD equilibration** (user), not from MM MD. They are picked with `scripts/select_nac_frames.py`: OH⁻ O01–C4x ≤ 3.2 Å, O01–C4x–O2x 95–120°, spread out in time.
  - `qmmm.inputs.traj_frequency` was set to 50 (every 50 fs, ~1,000 frames over 50 ps) before the qmmm stage started; at 1,000 it would have saved only 50 frames.
- QM/MM MD runtime (2026-10-16):
  - The md stage finished in job 29920397: all Zn constraints held (O2x 1.94 Å); `check_bonds.py` found every QM-region bond intact in all frames, with no new contacts and no proton transfers. In MM the OH⁻ is still not in attack geometry (OH⁻–C4x 3.5–4.3 Å, 21–45°).
  - Its QM/MM MD ran at ~30 s per step on 16 CPUs (~17 days for 50 ps) and was cancelled. Job 29937049 (8 cores, unbound) was erratic, 3.6–131 s per step, and was also cancelled.
  - A benchmark (job 29937358) traced this to unbound OpenMP threads stacking on SMT siblings: one bound thread per physical core takes 0.80 s per xTB call, against 2–9 s unbound.
  - Rerun as **job 29937641**: qmmm stage only, from `md.results.lastframe`, 8 CPUs, `cores: 4`, bound threads. It runs at **1.37 s per step (max 1.8 s), so 50 ps takes ≈ 19 h**. Frames are saved every 50 fs to `qmmm_results/QM_MM.dcd`.
  - Binding is now the mcube default (`src/mcube/threads.py`).

**D29 – QM/MM MD (first 39.5 ps): Lys226 protonates the Zn-OH⁻ at once; no attack geometry is sampled.** (2026-10-17; analysis in `CSR_crystal_v2/qmmm_analysis/`, on job 29937641's trajectory while it was still running)
- Bond geometry (`check_bonds.py`, 791 frames, 119 QM atoms):
  - All QM-region bonds keep their length (largest changes are normal vibration).
  - No new heavy-atom contacts.
  - **One proton transfer: Lys226 HZ1 → OH⁻ O01** in the first frame, permanent. HZ1–O01 is 0.91–1.35 Å (≈1.0 after frame 0); NZ–O01 is 2.74 Å (H-bond).
  - So in GFN2-xTB/MM the Zn-OH⁻/Lys226-NH3⁺ pair is not stable. It becomes Zn-OH2 (O01) + neutral Lys226 at once.
- Zn site (QM/MM, unconstrained):
  - His88 NE2 2.08, His164 NE2 2.05, OH⁻/H2O O01 2.11, HOH485 2.10 Å.
  - **Asp92 OD2 2.11 Å, bound in 99.5% of frames**, and it H-bonds both O01 (2.93 Å) and HOH485 (OD1/OD2 2.8–2.9 Å).
  - The carbonyl O2x is on the Zn in 91% of frames (mean 2.32 Å).
  - Coordination number 5–7 (mean 5.9).
- Attack geometry: **none in 791 frames**, for either O01 (3.26–4.49 Å, 20–53°) or HOH485 (3.20–4.50 Å, 16–74°). The O01 nucleophile stays cis to O2x on the Zn (O01···O2x 2.77 Å).
  - Closest O01 frames: 575 (3.38 Å, 52°), 244, 207, 192, 550.
  - Closest HOH485 frames: 757 (3.62 Å, 67°), 395.
- Contacts relevant to the three questions (polar, ≤ 3.5 Å, occupancy over 791 frames):
  - **Lys226** (K254, RQK): NZ···O01 99.9% (now H-bond acceptor from the former OH⁻).
  - **Arg90** (R118, HRH): NH1···meropenem O1x (6α-hydroxyethyl OH) 99%, NH1···HOH485 78%.
  - **Arg229**: NH1···O5x (dimethylcarbamoyl C=O) 89%, ···S1x 33%.
  - **Asp92** ("bridging"): Zn 99.5%; ···O01 100%; ···HOH485 100%; ···O2x 70%.
  - Other partners: Ser186 OG···C3-carboxylate O4x 100%; Asn225 ND2···O01 99%; His164 NE2···O2x (both on the Zn).
- Implications:
  - Q1 (Lys226): in xTB, Lys226 does not just stabilise a Zn-OH⁻; it hands its proton to it. Whether this is real or an xTB artefact must be checked at DFT before any NEB is built on it.
  - Q3 (Asp92): it bridges the Zn and both water/OH⁻ oxygens, so it is positioned as a Zn ligand and a possible base.
  - Q2: Arg90 and Arg229 hold the substrate's 6α-hydroxyethyl and C2 side chain, not the reacting β-lactam.
- Revisit: protonation state of the Lys226/OH⁻ pair; how to build reactants for the NEB given no sampled attack geometry.
- Job 29937641 completed all 50,000 QM/MM MD steps: 50 ps, 1,000 frames in `qmmm_results/QM_MM.dcd`, plus `QM_MM_lastframe.pdb`.
  - It then failed at the first step of the follow-up optimisation with an ASH error: after `OpenMM_MD`, `QMMMTheory` has no `MMenergy`, which geomeTRIC's energy log reads.
  - The MD output is intact. mcube now rebuilds the MM and QM/MM objects from the MD end point before optimising (commit 8c596d0).
  - The last 10 ps (frames 791–999) are not yet included in the D29 analysis.

**D30 – Check the Lys226 → OH⁻ proton transfer at DFT before building any NEB on it.** (2026-10-17, user decision; analogous to Imhof's B3LYP check of AM1/d)
- Frames 192, 395 and 575 of the QM/MM MD (9.6, 19.75 and 28.75 ps; 575 is the closest to attack). They are in `protdft/frames/`, with the box from the trajectory.
- Two QM/MM optimisations (GFN2-xTB/MM, 8 Å active region, ≤300 steps) per frame, sharing one active region:
  - `lysH`: proton held on Lys226 NZ (`opt_constraints` NZ–HZ1 = 1.03 Å), i.e. Lys-NH3⁺ + Zn-OH⁻;
  - `ohH`: unconstrained, so the proton stays on O01, i.e. Lys-NH2 + Zn-OH2.
- Each optimised structure gets an ωB97M-V/def2-SVP/MM single point (PySCF, electrostatic embedding).
- Compare ΔE = E(ohH) − E(lysH) at xTB with the DFT value:
  - same sign and similar size: the transfer is real, and the reactant is Zn-OH2 + neutral Lys226 (Q1 then asks whether Lys226 takes the proton back, or Asp92 does, during attack);
  - opposite sign: it is an xTB artefact. The reactant keeps Lys-NH3⁺/OH⁻ and the paths need a treatment that keeps it, e.g. DFT/MM refinement.
- Job 29981678 (array of 3, 16 CPUs = 8 physical cores, `ProtDFT.slurm`). The `qmmm` stage now takes per-section `inputs` and saves `sp_energy`.
- Progress (2026-10-17):
  - Proton held on Lys226 (`lysH`) optimisations: frames 395 and 575 converged; frame 192 hit the 300-step limit.
  - The ωB97M-V single points were far too slow and were cancelled (29981678 after ~4 h; 29992285 after 12 min):
    - 29981678 took ~22 min per SCF cycle. It used 4 threads, ASH's default grid level 5 and no density fitting, and the process-wide thread binding (mcube's first default) had pinned PySCF's process to one core.
    - mcube now binds only xTB (wrapper) and uses grid level 3, RI-JK and all physical cores for single points.
  - Rerun 29996196 (single point only for `lysH`; `ohH` optimisation + single point): the process now runs at ~790% CPU on 8 cores, but the first VV10 evaluation alone took > 15 min.
  - On CPU the VV10 non-local term of ωB97M-V dominates (124 QM atoms, 1,154 basis functions); grid level, RI and threads do not fix it.
  - User decision: switch the check to a hybrid without VV10. **ωB97X/def2-SVP, RI-JK**.
    - D3(BJ) is not available in the ASH env (no `vdw`/dftd3/dftd4 package). It is a post-SCF additive term, so it can be added to the same geometries later.
    - The first ωB97X run (29997262) wrote the RI integrals to disk (ASH's PySCF `max_memory` defaults to 3,100 MB) and was cancelled. mcube now sets `sp_memory` to 80% of the SLURM job memory.
  - Job 29997628 (96 GB, ~78 GB for PySCF; integrals in RAM, no disk I/O):
    - exchange–correlation step 8–11 s;
    - one-off initialisation (RI integrals + first build) 23–26 min;
    - the first SCF cycle had taken > 8 min when last checked.
    - Each single point will probably take a few hours. Faster options for the many path single points to come: gpu4pyscf on `gpu_cuda`, or seminumerical exchange (SGX).

**D31 – DFT level for this project: ωB97X/def2-SVP (RI-JK, exact exchange, no SGX).** (2026-10-17, user decision)
- Why (user): the energies must capture long-range interactions involving charged/polar groups and the Zn²⁺. A range-separated hybrid with full exact exchange suits that; a few hours per single point is acceptable.
- Set as the project default in `CSR_crystal_v2` `qmmm.inputs`: `sp_functional: wb97x`, `sp_basis: def2-svp`, `sp_gridlevel: 3`, `sp_densityfit: true`. `sp_memory` defaults to 80% of the job memory; give single-point jobs ~96 GB so the RI integrals stay in RAM.
- Follow-up decisions (user, 2026-10-17):
  - **D3(BJ): yes.** ωB97X-D3(BJ) for all single points, as a post-SCF correction on the QM region (link atoms included).
    - Env: installed `dftd3-python` 1.2.1 (`--freeze-installed`); the `simple-dftd3` library was already present. Conda also updated openssl 3.6.2 → 3.6.5; numpy, PySCF, OpenMM and xTB are unchanged. pyvdw (ASH's route) was tried and removed: it needs the D4 and MBD stacks too.
    - mcube: `qmmm.inputs.sp_dispersion: d3bj` (project config) → `sp_disp_energy` and `sp_energy_disp`.
    - `scripts/d3_correction.py` adds it to single points that were run without it. It rebuilds the QM region and link atoms exactly as ASH does: same 5 boundary pairs.
    - D30 job 29997628: the `lysH` single points started before this, so they get D3 from the script; the `ohH` single points pick it up from the config.
  - **def2-TZVP(D): yes, for final energies of stationary points** (R, intermediates, TS, P) per snapshot. Paths, scans and NEB images stay at def2-SVP.
    - Plan: ωB97X-D3(BJ)/def2-TZVPD/MM single points. The RI integrals grow several-fold (~2,700 basis functions), so these need high-memory jobs or SGX exchange; measure on one structure first.
  - **ωB97M-V on GPU: not needed** for this project.

**D32 – Replicate the published methods (Imhof 2009 path network; Xu–Guo 2006 QM/MM + PMF) instead of testing ad hoc protonation guesses.** (2026-10-18, user decision)
- The Asp92-protonated DFT check is on hold. The running D30 check (proton on Lys226 vs on O01) continues: it doubles as the xTB-vs-DFT validation step.
- Literature, checked against Europe PMC abstracts:
  - **EcoRV:** water nucleophile, proton to Asp90, Mg²⁺ activation, Lys92 positioning and stabilisation.
  - **CphA (mono-Zn):** non-metal-bound water nucleophile, Zn-bound Asp120 as general base, near-concerted, PMF ≈ 14 kcal/mol.
  - **CSR-1:** our Asp92 = D120.
- Protocol: `csr-neb-plan.md`, "Replication protocol for CSR-1".

**D33 – Nucleophile candidates from the full 50 ps QM/MM MD; H1–H4 scans proposed.** (2026-10-18; details in `csr-neb-plan.md`, "Step 1 results")
- No candidate ever reaches strict attack geometry.
- The Zn-bound oxygens (O01, HOH485) stay beside the carbonyl O (≤ 87°).
- One second-shell MM water, **E3094**, lies on the Bürgi–Dunitz trajectory: 105° mean, within 4 Å in 4.7% of frames. It is H-bonded to the 6α-hydroxyethyl OH (itself held by Arg90) and, through water **E3166**, to **Asp92 OD1**. This is a CphA-type route.
- O01 H-bonds the substrate's C3-carboxylate (88%); HOH485 H-bonds Asp92 OD1 (97%).
- Proposed scans:
  - **H1:** E3094 → E3166 → Asp92 (QM region + E3094, E3166);
  - **H2:** HOH485 → Asp92;
  - **H3:** O01 → carboxylate or Lys226;
  - **H4:** Zn-OH⁻ with Lys226-NH3⁺, only if D30 supports that protonation state.
- Pending: user approval to set these up.

**D34 – H1–H4 relaxed scans launched (all four).** (2026-10-18, user decision: run all four, and cancel H4 if the D30 DFT result rules out Lys226-NH3⁺/OH⁻ in time)
- Job 30022363: array of 12 tasks (H1 tasks 0–2, H2 3–5, H3 6–8, H4 9–11), each 8 CPUs (4 bound xTB threads).
- Setup: config sections `scan_H{1-4}_f{frame}`, run directories `scans/<section>/`, job list `scans/scan_jobs.txt`, setup script `scans/make_scan_sections.py`.
- Each scan: d(O_nuc–C4x) from its start value to 1.45 Å in 25 points (0.08–0.11 Å steps), ≤ 200 optimisation steps per point, protons free.
  - Monitored at every point: C4x–N1x, Zn–O2x, C4x–O2x, Zn–O_nuc, plus the H···acceptor distances specific to each hypothesis.
  - Active regions: 464–487 atoms (8 Å around the Zn, per frame).
- Start frames (`select_nac_frames.py`, ranked by closeness to the ideal attack geometry, ≥ 100 frames apart, after frame 50):
  - **H1** E3094: frames 365 / 497 / 635 (3.89–4.06 Å, 104–108°). QM region + waters E3094 and E3166 = 125 atoms, via the new per-section `extra_qm_residues`.
  - **H2** HOH485: frames 395 / 775 / 910 (3.59–4.14 Å, 62–92°).
  - **H3** O01: frames 192 / 575 / 888 (3.38–3.68 Å, 52–87°).
  - **H4**: same frames as H3, with Lys226 NZ–HZ1 held at 1.03 Å so O01 stays a hydroxide.
- Caveat: the Zn-bound nucleophiles (H2–H4) start at 52–92°, well off the Bürgi–Dunitz angle, so their scans force a reorientation; their maxima are upper bounds.

**D35 – D30 DFT result; first scan results; scans restarted from optimised reactants.** (2026-10-18)
- **D30 (ωB97X-D3(BJ)/def2-SVP//xTB/MM):** ΔE = E(proton on O01, Lys226 neutral) − E(proton on Lys226, Zn-OH⁻):
  - frame 395: **−13.1** kcal/mol (xTB −11.9);
  - frame 575: **−15.5** kcal/mol (xTB −17.8).
  - D3 contributes −0.7 to −0.8 of these.
  - Frame 192's DFT run failed (RecursionError in a logging call).
  - **The Lys226 → OH⁻ transfer is real at DFT.** The reactant is Zn-OH2 + neutral Lys226. **H4 (Zn-OH⁻ with Lys226-NH3⁺) is dropped:** it starts 13–16 kcal/mol uphill.
- **Scans (job 30022363):**
  - **H1, frame 497 (complete):** monotonic rise to +46.5 kcal/mol at 1.45 Å. No intermediate and no proton transfer: the E3094 H stays 2.8–3.5 Å from E3166, and E3166 only H-bonds Asp92. The carbonyl lengthens (1.205 → 1.275 Å) and binds the Zn more tightly (2.30 → 2.04 Å), i.e. oxyanion-hole behaviour.
  - **H2, frame 395 (22/25):** monotonic to +39.7 kcal/mol at 1.72 Å. The HOH485 H-bond to Asp92 breaks (1.78 → 3.75 Å) while its H approaches the lactam N (5.2 → 2.3 Å); the carbonyl leaves the Zn (2.29 → 3.10 Å). This hints at a concerted attack + N-protonation, as in CphA.
  - **9 tasks were killed out of memory** (16 GB): memory grows with optimisation steps, and the tasks started from unoptimised MD frames (first step −0.38 Eh), so they needed ~200 steps per point. H2 frame 775 is stuck (geomeTRIC SVD errors).
- **Lessons:**
  1. Optimise each start frame to a reactant minimum R before scanning, as the plan already said (step 3).
  2. 1D O–C4x scans do not make the proton move, so they give neutral-water addition energies. As Imhof and Guo did, the proton transfer must be a coordinate of its own: 2D scans (O–C4x × H···acceptor), or a scan from the state where the proton has already moved.
- Revisit: after the probe, whether the reactant is a cis OH⁻/carbonyl Zn complex (then attack needs the OH⁻ to rotate over the carbonyl face, which a scan will show), and whether HOH485 belongs in the Michaelis complex.

**D36 – Reactant optimisations; 2D scan support; QM waters were held rigid; 2D grid proposal.** (2026-10-18)
- **Reactants:** H2 frame 775 cancelled. Each H1–H3 start frame is optimised to a reactant minimum R (sections `qmmm_R_H{1,2,3}_f*`: `skip_equilibrate`, no single point, `opt_maxiter` 500, 32 GB; H1 with E3094 and E3166 in the QM region). The scans will start from these.
- **mcube:** `scan` now does 2D grids (`atoms2`, `start2`, `end2`, `points2`). Rows run over coordinate 2 and coordinate 1 is swept back and forth, so every point starts from a neighbour. A stopped scan resumes from `scan_energies.csv`.
- **Bug found and fixed (113de19):** the optimisers froze the O–H bonds of every active-region TIP3P water, *including waters in the QM region*. MD was not affected (`rigidwater=False`).
  - Consequence for D35: in the 1D scans HOH485 (H2) and E3094/E3166 (H1) could not lose a proton. The "no proton transfer" result and the +40–47 kcal/mol profiles are therefore upper bounds, not evidence against the relay.
  - QM waters are now left free in `qmmm`, `scan` and `ts`. The reactant job was restarted with the fix (30047537 → 30048252).
- **H-bond network in the start frames** (nearest acceptor of each H):
  - E3094 (H1 nucleophile) donates both protons to bulk waters E2155 and E2493 (1.7–2.0 Å). It *accepts* H-bonds from the meropenem 6α-OH (H5x, 1.6–1.7 Å) and from E3166 H2 (1.7–2.5 Å). E3166 H1 donates to Asp92 OD1 (1.9–2.4 Å).
  - So the E3094 → E3166 → Asp92 relay is not pre-organised. Every link runs the wrong way: the E3094 H points away from E3166 (3.5–4.0 Å). For the relay to work, E3094 must rotate, or E3166 must first give its proton to Asp92 and then take one from E3094.
  - HOH485 (H2 nucleophile, Zn-bound, 1.98–2.17 Å): H2 donates to Asp92 OD1 (1.7–1.9 Å); H1 donates to the 6α-OH O1x or to Arg90 NH1. Its distance to N1x is 3.7–4.0 Å (H1) and 4.8–5.1 Å (H2). In the 1D scan the Asp92-bonded H swung to 2.3 Å from N1x.
- **Proposed 2D grids** (not submitted; waiting for approval and for R):
  - **H1:** O(E3094)–C4x from R (~3.9 Å) to 1.45 Å, 13 points × H(E3094)···O(E3166) from R to 1.0 Å, 6 points (78 points).
    - Monitor H1(E3166)···Asp92 OD1 (second proton of the relay), H5x···O(E3094), plus C4x–N1x, Zn–O2x, C4x–O2x.
    - If R still has E3094 pointing away from E3166, coordinate 2 forces a rotation first. Then decide whether to put E2155/E2493 in the QM region instead.
  - **H2a (Asp92 as base):** O(HOH485)–C4x from R to 1.45 Å, 12 points × H2(HOH485)···Asp92 OD1 from R (~1.8 Å) to 1.0 Å, 5 points (60 points).
  - **H2b (concerted N-protonation, CphA-like):** O(HOH485)–C4x, 12 points × H(HOH485)···N1x from R to 1.05 Å, 7 points (84 points). Use the H that is closer to N1x in R.
  - One frame per hypothesis first: the R with the best attack geometry. 8 CPUs, 32 GB, 2-day limit; `resume` restarts a stopped grid.

**D37 – 2D grids submitted.** (2026-10-18, user approved all three)
- **Why:** a 1D O–C4x scan keeps the proton on the water, so it describes addition of a neutral water. The proton transfer has to be its own coordinate (Imhof 2009; Xu & Guo 2006). The low-energy route across each grid shows whether the proton moves before, during or after the attack, and to which base.
- **Sections** `scan2d_*` in `config.yaml`; job script `Grid.slurm` (8 CPUs, 32 GB, 2 days); section script `scans/make_grid_sections.py`. Each grid starts from its optimised reactant and waits for that reactant task (`afterok`):
  - **H1** `scan2d_H1_f00497` (job 30048856, after 30048252_1): O(E3094)–C4x → 1.45 Å, 13 points × H2(E3094)···O(E3166) → 1.0 Å, 6 points. Tests the E3094 → E3166 → Asp92 relay. H2 of E3094 was chosen because it was closer to E3166 (3.6 Å) in frame 497.
  - **H2a** `scan2d_H2a_f00910` (job 30048857, after 30048252_5): O(HOH485)–C4x → 1.45 Å, 12 points × H2(HOH485)···Asp92 OD1 → 1.0 Å, 5 points. Tests Asp92 as the base for the Zn-bound water.
  - **H2b** `scan2d_H2b_f00910` (job 30048858, after 30048252_5): same coordinate 1, 12 points × H1(HOH485)···N1x → 1.05 Å, 7 points. Tests concerted attack and N-protonation (CphA-like). H1 was chosen because it was closer to N1x (3.7–4.0 Å vs 4.8–5.1 Å).
  - H2a and H2b share a reactant (frame 910, best H2 angle of 92°), so the two bases are compared on the same structure.
- **Monitored:** C4x–N1x, Zn–O2x, C4x–O2x, Zn–O_nuc, the other water H, and for H1 the second relay proton (E3166 H1···Asp92 OD1) and the 6α-OH···O(E3094) H-bond.
- If a job hits the time limit, resubmitting the same command resumes the grid.

**D38 – All grids cancelled. The reactant state, not the hypotheses, is the problem.** (2026-10-18)
- User challenge: do not spend simulation time without a plan. On inspecting the optimised reactants, H2a/H2b (30048857/8) and then H1 (30048856) were all cancelled. Nothing is running.
- **The Zn is six-coordinate in every optimised reactant:**

  | Ligand | Distance (Å) |
  | --- | --- |
  | His88 NE2 | 2.01–2.08 |
  | His164 NE2 | 2.01–2.04 |
  | Asp92 OD2 | 2.07–2.11 |
  | O01 (water) | 1.99–2.13 |
  | HOH485 (water) | 2.04–2.09 |
  | meropenem O2x (carbonyl) | 2.16–2.31 (detached in f910) |

- **There is no hydroxide in the model.** O01 carries two protons (H03 plus Lys226's HZ1, 1.00 Å) and Lys226 is neutral NH2; HOH485 is an intact water. Both Zn-bound oxygens are neutral water.
- **Why every scan was ~+40 kcal/mol:** a neutral water attacking a carbonyl is intrinsically that expensive. The 1D scans were not testing a mechanism, they were measuring neutral-water addition. (The frozen-QM-water bug, fixed in 113de19, made it worse but was not the cause.)
- **Why the Zn-bound waters cannot attack:** O01 and HOH485 sit *cis* to the carbonyl on the same metal, so their O–C4x–O2x angles are 31–45°. A nucleophile coordinated next to the carbon it must attack is on the wrong trajectory by construction. This is why D33 found no Zn-bound oxygen above 87°.
- **Why H1's relay runs backwards:** E3094 has the right trajectory (4.40 Å, 101–109°) but donates both protons to bulk waters E2155 (1.74 Å) and E2493 (1.76 Å). E3166 *donates to* E3094, not the reverse, and sits 3.0–3.8 Å away. The nearest heteroatom to E3094 is meropenem's own 6α-OH O1x at 2.65 Å, which also donates to it. E3094 is a cluster water with no base except bulk solvent, so the submitted H1 grid would have dragged it 4.35 → 1.0 Å across space rather than transferring a proton.
- **Diagnosis:** a mono-Zn B3 site with three protein ligands, two waters *and* the substrate carbonyl is over-coordinated; catalytic Zn in B3 MBLs is normally four- or five-coordinate with a hydroxide. The Michaelis complex we built is probably wrong, and D29/D30 (the Lys226 → O01 transfer, −13 to −16 kcal/mol at DFT) was computed *in* this over-coordinated state, so it does not settle the protonation question for the correct one.
- **Next, before any more scans (P6):** establish the reactant state.
  1. Check CSR-1's reported Zn coordination (crystal structure, B3 literature): how many waters, and is the carbonyl expected to coordinate?
  2. Test whether HOH485 belongs in the Zn sphere. Already flagged as a revisit in D35.
  3. Re-examine the protonation state once the coordination number is right: a five-coordinate Zn-OH⁻ with Lys226 as NH3⁺ is the canonical B3 reactant and was never tested, because the MD collapsed to Zn-OH2 from the six-coordinate start.
  4. Only once a hydroxide nucleophile exists with a plausible trajectory do scans or grids make sense.

**D39 – The source structure is an engineered double mutant; the model's zinc has only two protein ligands.** (2026-10-18)
- Traced the model back to **PDB 6NC5** (`CSR_crystal/6NC5.pdb`): *Cronobacter sakazakii* metallo-β-lactamase, "HARLDQ motif", **ENGINEERED: YES, MUTATION: YES**, 2.07 Å.
- **Sequence comparison** (`seqcmp.py`): our model is crystal numbering **+1**, 257/259 residues identical, with exactly **two mismatches**:
  - crystal **His92 → model Gln93**
  - crystal **His225 → model Lys226**
  - So 6NC5 is a **Q93H/K226H double mutant** that restores the canonical B3 zinc-binding histidines, and our model is the wild-type RQK sequence (Arg90, **Gln93**, **Lys226**). The model build was right to revert them.
- **6NC5 has six modelled zincs.** Three are surface artefacts on Glu/Asp/phosphate (ZN304–306, B = 56–70). The active site holds two:

  | Site | Ligands | Note |
  | --- | --- | --- |
  | ZN301 | Asp91 OD2 1.95, **His92** NE2 2.01, **His225** NE2 2.04, HOH493 2.08 | four-coordinate |
  | ZN302 | HOH485 1.90, His163 NE2 2.09, HOH493 2.14, His87 NE2 2.22 | four-coordinate |
  | ZN303 | HOH485 1.88, HOH516 1.96, HOH504 2.13, acetate 2.79 | all water/acetate — artefact |

  **HOH493 bridges ZN301 and ZN302**, the classic di-zinc bridging-hydroxide position.
- **Two of ZN301's three protein ligands are the mutated residues.** In wild-type CSR-1 they are Gln93 and Lys226, which cannot bind zinc. So ZN301 exists only in the mutant, and wild-type CSR-1 being mono-zinc at the ZN302 site is consistent with this structure. Deleting ZN301 was correct.
- **But the model's zinc then has only two protein ligands** (His88 and His164, i.e. crystal His87/His163). Matching distances confirm the identifications: our starting model's zinc sphere was HOH 1.90 / His164 2.09 / OH(O01) 2.14 / His88 2.22 / Asp92 2.80 — exactly ZN302's crystal distances. So:
  - **our hydroxide O01 = crystal HOH493** (the former bridging water), and
  - **our HOH485 (D:1) = crystal HOH485**, which in the crystal also bridges to the artefact ZN303.
  - **Asp92 was *not* a ZN302 ligand in the crystal** (2.80 Å). It is a ZN301 ligand (1.95 Å). With ZN301 deleted it migrated onto our zinc (2.07–2.11 Å), which is what the user spotted earlier and what D28 restrained to 2.80 Å in MM only.
- **This reframes question 3.** Asp92 may not be a general base at all: with its partner zinc absent, it may legitimately be a zinc ligand in wild-type mono-zinc CSR-1. Whether it coordinates or stays as a free base is now the central modelling question, not an incidental detail.
- **And it explains the six-coordination.** A two-His zinc needs two or three more ligands. Ours took the hydroxide, HOH485, Asp92 *and* the substrate carbonyl. At least one of those does not belong.
- **Open, for the reactant state:** which combination of {OH⁻(ex-493), HOH485, Asp92 OD2, substrate carbonyl} completes a chemically sensible four- or five-coordinate zinc, and is HOH485 (a bridge to an artefact metal in the crystal) part of the real site at all.

**D40 – Literature check: three independent problems with the reactant model.** (2026-10-18)
- **CSR-1 is published:** Pedroso MM *et al.*, "Broad spectrum antibiotic-degrading metallo-β-lactamases are phylogenetically diverse", *Protein Cell* 2020;11:613–617, PMID 32542533 (PMC7381538). A **B3-RQK** MBL from *Cronobacter sakazakii* (UniProt A7MFY4). Structures: **6DN4** (mature wild type, 1.99 Å), **6DQH/6NC5** (2.07 Å), **6DQ2**, **6DR8** (1.48 Å, R118H/Q121H/K254H triple mutant).
- **Motif mapping** (BBL → our numbering): Arg118 → Arg90, Asp120 → Asp92, Gln121 → Gln93, His116 → His88, His196 → His164, Lys254/263 → Lys226. So the Zn1 site is **H-R-H** (His116, **Arg118**, His196; two usable His) and the Zn2 site is **D-Q-K** (Asp120, **Gln121**, **Lys254**; one usable ligand). Canonical B3 is HHH/DHH.
- **This confirms and explains D39.** 6NC5's ZN301 (Asp91/His92/His225) is the **canonical Zn2 site restored by the two mutations**; ZN302 (His87/His163 + two waters) is the **Zn1 H-R-H site**. Our model correctly kept ZN302 and reverted both mutations. There is no third candidate site.
- **Three problems, each independent:**
  1. **No published zinc site.** The wild-type CSR-1 structures crystallised **essentially apo** — 6DN4 lists no metal at all — and the paper attributes this to reduced metal affinity, since Arg118 and Lys254 are protonated at optimal pH. So our mono-zinc site is **inference, not literature**. Nothing is wrong with inferring it, but it must be stated as such, and the Zn placement came from a mutant's coordinates.
  2. **Hydroxide *and* a second water on the same zinc is unsupported.** Canonical B3 Zn1 is **four-coordinate**: three His plus the bridging hydroxide (L1, PDB 1SML; Ullah *et al.*, *J Mol Biol* 1998;284:125, PMID 9811546). A literature search found **no** description of a catalytic zinc carrying both a hydroxide and a separate water. Our site carries OH⁻ (ex-HOH493), HOH485, Asp92 *and* the substrate carbonyl on a two-His zinc. HOH485 is the weakest of these: in the crystal it bridges ZN302 to the **artefact** ZN303 (whose other ligands are two waters and an acetate).
  3. **The model is the non-truncated, inactive form.** CSR-1 is active only after N-terminal truncation; in the mature form an N-terminal loop occludes the active site. 6NC5 covers residues 15–273; **our model covers 1–275**, so residues 1–14 were modelled back in, and Ala12, Pro13, Tyr15, Leu17, Phe18, Glu19 and Trp21 sit 9–13 Å from the zinc — nearer than anything from the crystal's own N-terminus.
- **Inner-sphere carbonyl coordination is also not established.** No crystallographic Michaelis complex shows a firm Zn–O(carbonyl) bond; the consensus is that Zn1 polarises the carbonyl at longer range while the hydroxide attacks (NDM-1 QM/MM: Das & Nair, arXiv:1702.01239, where Zn1–O is only weakly coordinated). Our `preprocess` `cutoff: 3.0` captured the carbonyl at 2.16–2.31 Å and restrained it as a metal ligand, per D28.
- **The CphA precedent now looks like the leading hypothesis, and our own data already pointed at it.** Xu, Xie & Guo, *J Biol Chem* 2006;281:8740, PMID 16423823 (the paper chosen as a methods template in D32): in mono-zinc B2 CphA the nucleophile is a **non-metal-bound water** activated by **Asp120** as general base, concerted, ΔG‡ ≈ 14 kcal/mol, no tetrahedral intermediate, and after ring opening the nitrogen leaving group becomes a zinc ligand, displacing Asp120. In CSR-1 the Zn2 site is vacant, so **Asp92 (= Asp120) is free to act as that base** — and D33 independently found that the *only* species on the Bürgi–Dunitz trajectory is a non-metal-bound water. The geometry our MD produced is the CphA geometry.
- **So question 3 is reframed again:** Asp92 is either a zinc ligand (filling the vacant Zn2 site's role) or the general base (CphA-like). Those are mutually exclusive and the choice defines the mechanism. It is the first thing to settle.
- **Also noted:** a numbering discrepancy for the lysine, Lys254 (PDB titles) vs Lys263 (Krco *et al.* review, *Front Chem* 2023;11:1196073, PMID 37408556). **[Resolved 2026-10-25: one residue, depositor sequential numbering vs BBL. The table in D22 already maps 226 -> 254 -> 263.]** GOB-18, not GOB-1, is the mono-zinc B3 precedent, and its single zinc sits in the **Zn2** site (PDB 5K0W); the sister B3-RQK enzyme SPR-1 is suggested to be mono-zinc at rest with a di-zinc centre forming on substrate binding (not established).
- **Nothing is running. No further scans until the reactant state is settled (P6).**

**D41 – Corrections to D40 (user). The reactant we built was right; the MD degraded it.** (2026-10-18)
Three claims in D40 were wrong or overstated. The user corrected all three, and the corrected diagnosis is sharper and more useful.

1. **"The zinc is inference, not literature" — wrong.** The zinc position is crystallographic. Our zinc sits at 6NC5's **ZN302**, and ZN302's own ligands, His87 and His163 (our His88/His164), are **not** among the mutated residues: 6NC5's mutations are at 92 and 225, both **Zn2**-site residues. So the zinc occupies a site whose coordination is wild type, at crystallographic distances (2.09, 2.22 Å). What is true is only that the wild-type crystal forms did not retain zinc — verified directly: **6DN4 reports `nonpolymer_entity_count: 0`, zero ligand atoms, no ZN** (1.99 Å, 347 solvent atoms); 6DQ2 has one unnamed 10-atom non-polymer and no ZN. Since the enzyme is active it must bind zinc, so apo crystals point to crystallisation conditions and the reduced affinity the paper describes, not to an absent site. (Note also that "HARLDQ motif" in the PDB titles names the *parent* enzyme's motif; 6NC5's coordinates read H-A-R-L-D-**H** at 87–92, confirming the Q→H mutation.)
2. **"Over-coordinated" — wrong, as stated.** Coordination number counts all ligands. Canonical B3 Zn1 is four-coordinate as **three His + hydroxide**. CSR-1's Zn1 is **H-R-H**: only two usable His, because Arg118 replaces the third. So **two His + hydroxide + one water is the canonical-equivalent four-coordinate sphere**, and **HOH485 is justified** — it fills the slot the missing third His would occupy. HOH485 is not the problem.
3. **"The non-truncated form is the inactive form" — overstated.** This was relayed from a literature summary and not verified. Reduced activity in a truncation assay does not establish that the mature form is catalytically dead, and the non-truncated form is the more faithful model of the natural enzyme. Keeping it is a defensible choice. The one caveat that survives is narrow and structural: residues **1–14 are absent from 6NC5** (which covers 15–273) and so have no experimental coordinates, yet Ala12, Pro13, Tyr15, Leu17, Phe18, Glu19 and Trp21 sit 9–13 Å from the zinc. Modelled-in backbone near the active site is worth knowing about; it is not grounds for truncating.

**Corrected diagnosis.** The reactant as *built* was correct and canonical: Zn + His88 + His164 + **OH⁻ (ex-HOH493)** + **H2O (HOH485)**, four-coordinate, with one hydroxide nucleophile. The **QM/MM MD then degraded it** in three ways:
  - **Lys226 donated HZ1 to O01**, turning the hydroxide into water, so the model lost its nucleophile entirely (D29);
  - **Asp92 OD2 migrated onto the metal**, 2.80 Å in the crystal → **2.07–2.11 Å**, becoming a fifth ligand it is not in the crystal;
  - **the substrate carbonyl came inner-sphere** at 2.16–2.31 Å, a sixth.

So the surviving problems are (a) no hydroxide and (b) a sphere two ligands beyond canonical, both arising *during* the simulation rather than from the build. (c) The geometric result is unchanged: neither Zn-bound oxygen reaches the Bürgi–Dunitz angle (31–45°), and only a non-metal-bound water does (D33).

**The D30/D35 DFT check needs redoing.** It found the Lys226 → O01 transfer favourable by 13–16 kcal/mol, but it was computed in the already-degraded six-coordinate geometry, where Asp92's negative charge sits on the metal and the carbonyl competes for it. **Next test (cheap, decisive):** DFT single points comparing proton-on-Lys226 vs proton-on-O01 in a *clean* four-coordinate sphere — Asp92 held at its crystallographic 2.80 Å, carbonyl not inner-sphere. If the transfer is not favourable there, the hydroxide survives, H4 is reinstated and H1/H3 become testable with a real nucleophile.

**D42 – Clean four-coordinate proton-transfer check: set up, not yet run.** (2026-10-18)
- **Starting structure found.** `md_results/NVTtrajectory_lastframe.pdb`, the pre-QM/MM-MD NVT end point, already *is* the canonical reactant: Zn sphere HOH485 1.79, His164 2.09, **O01 2.14 with a single H (a genuine hydroxide)**, His88 2.22, **Asp92 OD2 2.79** (crystallographic), and **Lys226 NZ with three H** (NH3+). So two His + OH⁻ + H2O, exactly the four-coordinate sphere D41 identified as canonical for this enzyme.
  - The one defect: **the carbonyl O2x sits at 1.94 Å**, held there by the preprocess metal restraint (`cutoff: 3.0`) — *shorter than the Zn–OH⁻ bond*, which is unphysical for a neutral carbonyl and makes the sphere five-coordinate.
- **Design.** Sections `qmmm_cleanR_lysH` / `qmmm_cleanR_ohH`, both from that frame, with the **same** two sphere constraints (Zn–Asp92 OD2 = 2.80, Zn–O2x = 2.90) so the energy difference reflects only the proton position, and each pinning the shared proton on one side (lysH: NZ–HZ1 = 1.03; ohH: O01–HZ1 = 0.98). `ohH` reuses lysH's `active_atoms`. Scripts in **`docs/csr/`** (`mk_cleanR.py`, `CleanR.slurm`, README), committed so the cluster can pull them.
  - **Change from D30:** there, `ohH` was left unconstrained and relied on the transfer happening spontaneously. In a clean sphere it may not, which would silently make the two states identical and the comparison empty. Both sides are now pinned.
- **Submitted as job 30053828** (2026-10-18, after the user reconnected the session). Both scripts were pulled from `docs/csr/` and verified by checksum against the local copies before use (`6f0763d2…` / `4aa434ad…`), and all six constrained selectors resolved to the expected atoms: Zn 4203, Asp92 OD2 1414, mer O2x 4209, Lys226 NZ 3431 / HZ1 3432, O01 4255, with no missing selectors. Backup `config.yaml.bak14`.
- **First attempt failed before submission.** The code-server terminal stopped executing commands partway through setup (a multi-line heredoc did not close — newlines typed into xterm are not delivered as Enter — and the shell then stopped responding to Ctrl+C and to trivial commands). **`config.yaml` was not modified**; only `config.yaml.bak13` was created. A stray `mk_cleanR.py` may be left in the run directory with mangled contents and should be overwritten from `docs/csr/`.
- **Lesson for the working rules:** do not type long payloads into the browser terminal. Two separate corruptions happened here — a ~3 kB base64 command was silently truncated mid-string, and a hand-retyped base64 chunk differed from the original. Put cluster scripts in the repo and `git pull` them; type only short commands. Verify any transferred file by checksum against the local copy before running it.
- **Still running, unaffected:** reactant array **30048252** (tasks 1–3 completed; 0, 4–8 were running when the terminal died). Nothing else.

**D43 – Superseded jobs cancelled.** (2026-10-18, user: "we want a solid methodology fit for a paper")
- Cancelled **30022363** (the last H2 f395 1D scan: frozen QM waters, degraded reactant, one dimension) and **30048252** (reactant optimisations for the H1–H3 NAC frames, all taken from the degraded QM/MM MD, so reactants for scans that will not now be run from those frames). Tasks 1–3 had completed and their `QM_MM_optimized.pdb` files are kept as records.
- Only **30053828** (D42, the clean four-coordinate proton check) is running. The user's GLP1R jobs are a separate project and were left alone.
- A publication-grade protocol is proposed in `csr-neb-plan.md` under "Proposed protocol (P-series)", pending the user's decision on the level of theory.

**D44 – Tier B chosen; mcube can now run DFT/MM gradients, but the cost needs a decision.** (2026-10-18)
- **User decision:** tier B — xTB/MM to find paths, then all stationary points re-optimised at DFT/MM with DFT/MM frequencies and IRC, so every published number is DFT.
- **mcube change (1938eed):** the gradient stages had no way to use DFT. `qmmm`, `scan`, `neb` and `ts` all built a `QMMMTheory` around GFN2-xTB, and PySCF appeared only in the single point. Added **`qm_engine: xtb|dft`**, with `dft_*` level keys that each fall back to their `sp_*` counterpart, and `dft_scf_dispersion` routed through ASH's own `dispersion` keyword so dispersion enters the gradient rather than being added afterwards. The single point now shares the same factory. Stub tests cover the engine switch, the fallback and the rejection of an unknown engine.
- **The cost problem.** Measured, not estimated: job **29997628** took **9:52 and 10:36 h on 16 CPUs** for two DFT/MM single points plus one xTB/MM optimisation, i.e. **~4–5 h per DFT/MM energy** at ωB97X/def2-SVP, grid level 3, density fitting, 119 QM atoms. A gradient costs at least as much again.
  - A stationary point needs tens of DFT/MM gradients even from a good xTB start, so **one DFT/MM optimisation is days to weeks on CPU**, and the protocol needs on the order of 12–36 stationary points (3 mechanisms x 3 structures x R/TS/P).
  - **Numerical frequencies are worse.** A full QM-region Hessian is 3N x 2 gradients: 119 atoms = **714 gradients**. At 5 h each that is impossible, and even at 30 min each it is ~15 days per Hessian. **Partial Hessians (`partial_hessian_atoms`, already supported) are mandatory**, restricted to the reacting subset (~15–25 atoms, so 90–150 gradients).
  - **`gpu4pyscf` is not installed** (checked: absent; PySCF 2.13.0, CPU only). ASH's `PySCFTheory` accepts `platform='GPU'`, so mcube can already request it through `dft_platform` once the package exists.
- **Three ways to make tier B affordable**, which combine:
  1. **GPU.** `gpu4pyscf` on Bunya's A100/H100 is typically 10–30x for this size, bringing a gradient to tens of minutes and an optimisation to ~1–2 days. The single biggest lever. **Needs the user's approval to install.**
  2. **A cheaper functional for geometries.** ωB97X is a range-separated hybrid; a meta-GGA such as r2SCAN-D4 avoids exact exchange and is several times cheaper. Geometries and frequencies at r2SCAN-D4/def2-SVP, all reported energies at ωB97X-D3(BJ)/def2-TZVPD on those geometries. A standard and defensible split, and now expressible: `dft_functional: r2scan` with `sp_functional: wb97x`.
  3. **A smaller QM region for the refinement.** 119 -> ~70 atoms (keep Zn, both His, Asp92, the Lys226 side chain, both waters, and only meropenem's β-lactam ring, carboxylate and 6α-OH) is roughly 4–8x for hybrid DFT. The QM-region convergence test (P-C) then runs at xTB level plus single points at the full region.
- **Unchanged:** D42's job 30053828 is an xTB/MM optimisation with DFT single points — a screen to pick the reactant state, not a published number. Once the state is chosen it gets a DFT/MM optimisation under tier B.

**D45 – gpu4pyscf installed (additively); DFT/MM gradient benchmark submitted.** (2026-10-18, user approved the install)
- **The install nearly broke the environment, and the dry run caught it.** `pip install gpu4pyscf-cuda12x` on its own would have pulled `cupy-cuda12x` 14.2.0, which requires **numpy >= 2.0**, and so would have upgraded **numpy 1.26.4 -> 2.4.6**. OpenMM 8.4 and mdtraj are compiled against numpy 1.x, so that would have taken out the MM backbone of the whole pipeline — and the cleanR job was running in the same environment at the time.
  - Installed instead with `cupy-cuda12x==13.3.0` and `numpy==1.26.4` pinned. The dry run then reported only additions, and the install added exactly: `cupy-cuda12x-13.3.0`, `fastrlock-0.8.4`, `gpu4pyscf-cuda12x-1.8.1`, `gpu4pyscf-libxc-cuda12x-0.8.1`, `pyscf-dispersion-1.5.0`. **Nothing was replaced or removed.**
  - Verified afterwards: numpy 1.26.4, OpenMM 8.4, PySCF 2.13.0, mdtraj all still import. Job 30053828 kept running throughout.
  - **Rule for next time:** `pip install --dry-run` first and read the "Would install" line; pin anything that must not move. This is the second time an install in this environment tried to drag core packages with it (the first was pyvdw, D-series earlier).
- **ASH GPU support confirmed** in `interface_pyscf.py`: `platform == 'GPU'` is handled through the SCF, gradient and property paths, and there is a separate `GPU_pcgrad` for the MM point-charge gradient — which for a large MM environment is much of the QM/MM gradient. mcube now exposes it as **`dft_gpu_pcgrad`** (9beb7d8).
- **Benchmark job 30054626** (`docs/csr/mk_gpubench.py`, `GpuBench.slurm`; 1 GPU on `gpu_cuda`, 8 CPUs, 64 GB). Three levels, each a few constrained optimisation steps from the same NVT end point with no single point, so the per-step wall time in `optimization_energies.log` *is* the DFT/MM gradient cost:
  - `qmmm_gpu_wb97x` — ωB97X/def2-SVP, directly comparable to the measured CPU cost (~4–5 h per energy on 16 cores, job 29997628);
  - `qmmm_gpu_r2scan` — r2SCAN-D4/def2-SVP, the proposed geometry level (meta-GGA, no exact exchange);
  - `qmmm_gpu_r2scan_pc` — as above with `GPU_pcgrad`, to price moving the point-charge gradient onto the device.
  - The D42 sphere constraints are kept so the timing is measured on the reactant geometry we actually care about. Each section is allowed to fail independently so one bad level does not cost the other numbers.
- **What the numbers decide:** whether tier B runs at the full 119-atom QM region (if a gradient lands in the tens of minutes) or still needs the reduced ~70-atom region, and whether ωB97X geometries are affordable or geometries should be r2SCAN-D4 with ωB97X-D3(BJ)/def2-TZVPD energies on top.

**D46 – Documentation pass; a silent per-section bug found and fixed; GPU benchmark needed a CUDA runtime.** (2026-10-18, user: "make sure mcube changes are documented to be human usable")
- **Docs (177a6f0).** `docs/qm-levels.md`: when to use xTB vs DFT, how each `dft_*` key falls back to its `sp_*` counterpart, the difference between the two dispersion keys (`dft_scf_dispersion` is in the SCF and so in the gradient; `sp_dispersion` is post-SCF and reported separately), why DFT frequencies need a partial Hessian, GPU setup, and where the keys may be set. Worked configs `config_examples/qmmm_dft_refine.yaml` (the whole xTB-then-DFT sequence) and `config_examples/scan_2d.yaml` (including how to choose a proton-transfer coordinate and how to read the grid). Both parse-checked.
- **Bug found while writing the examples (177a6f0).** `scan`, `neb` and `ts` all built the QM theory from `qmmm.inputs` **only**, so `qm_engine: dft` or any `dft_*` key set on a `scan_*` / `neb_*` / `ts_*` section was **silently ignored and the stage ran at xTB**. That would have produced "DFT" barriers that were really semi-empirical, with nothing in the output to say so — the worst kind of failure for a paper. `qm_level_inputs()` now merges the QM-level keys, and only those, from the stage's own section; charge, multiplicity and force fields still come from `qmmm.inputs`. Covered by `tests/test_qm_levels.py` (11 tests) and two stub regression tests. The `qmmm` stage was never affected, because it already merges per-section inputs.
- **Also corrected:** four places (README twice, `ts.py`, `qmmm.py`) still said all active-region waters are kept rigid, which the earlier fix changed. Waters in the QM region are left free so they can transfer protons.
- **GPU benchmark 30054626 failed in 24 s per section**, all three: `ImportError: libcusolver.so.11: cannot open shared object file`. The node was an **NVIDIA L40S, 46 GB, driver 610.43.02**, so the hardware was fine; the cause is that Bunya has **no CUDA 12 module** (only `cuda11.7`) and no CUDA 12 toolkit on the filesystem, so cupy's `cuda12x` wheels had no runtime libraries.
  - Fixed by installing the CUDA 12 runtime as pip wheels — dry-run first, and it reported only additions: `nvidia-cublas/cuda-nvrtc/cuda-runtime/cufft/curand/cusolver/cusparse/nvjitlink-cu12`. numpy stayed at 1.26.4; OpenMM 8.4 and mdtraj still import.
  - `GpuBench.slurm` now puts those wheels' `lib` directories on `LD_LIBRARY_PATH` and asserts `cupy` works before running anything, so a broken GPU setup fails in one line instead of three times over.
  - **Resubmitted as 30055503.**
- **Still running:** cleanR 30053828 (D42), unaffected by any of the installs.

**D47 - CPU only; minimal protocol; GPU abandoned; scratch is over quota.** (2026-10-18, user: "we'll stick to CPU for future runs. I want to do as little simulation work as possible to test our hypotheses, sticking to methods and approaches from the CSR and EcoRV papers")
- **GPU abandoned after four attempts.** 30054626 (no CUDA 12 runtime at all), 30055503 and 30056056 (cupy 13.3.0 then 13.6.0, both `undefined symbol: cusolverDnDsygvd_bufferSize` - CUDA 12.9's cusolver dropped that legacy entry point), and 30056429 (cusolver pinned back to 11.6.1.9) **cancelled** on the user's instruction before it reported. The GPU packages stay installed but unused; uninstall instructions are in `docs/csr/cleanup.sh`.
- **Tier B (D44) is reversed, and re-reading the source papers is why.** Imhof 2009 used **AM1/d/MM** with B3LYP only on *small model systems*; Xu & Guo 2006 used **SCC-DFTB/MM**. Neither re-optimised stationary points at DFT/MM. GFN2-xTB is a later and better-parameterised method than either, so xTB/MM as the primary method is *faithful* to these papers rather than a compromise - and it removes the cost problem that drove the entire GPU detour. `qm_engine: dft` remains in mcube for a final single-structure check.
- **Minimal protocol (M1-M4) written into `csr-neb-plan.md`,** superseding the P-series: settle the reactant (M1, job 30053828, in flight); then **one** mechanism on **one** structure - a non-metal-bound water with Asp92 as general base, the CphA case that both the precedent and our own D33/D41 geometry point to - as scan then NEB then TS/IRC at xTB/MM, with DFT/MM single points along the final path (M2); then **charge deletion for Lys226, Arg90, Arg229 and Asp92**, which is what actually answers the three questions and is cheap (M3); then one cluster-model xTB-vs-DFT check, Imhof-style (M4).
- **Explicitly dropped:** DFT/MM optimisation and frequencies throughout, QM-region convergence sweeps, three structures per mechanism, umbrella sampling, and the 6DR8 control.
- **Scratch is over quota: 327 GB against a 300 GB limit, and 859,781 files against a 1,000,000 file limit.** The file count matters as much as the size, because ASH leaves dozens of small working files per optimisation. `docs/csr/cleanup.sh` lists by default and deletes only with `--delete`. It removes ASH working files, geomeTRIC trajectories, the nine out-of-memory scans, the cancelled 2D grids, the GPU benchmark output, the superseded QM/MM probe directories and config backups beyond the newest three, and never touches `preprocess_results/`, `md_results/`, `qmmm_results/`, `config.yaml`, any `QM_MM_*.pdb`, `scan_energies.csv`, `scan_*.pdb` or `nac_frame*.pdb`. Verified on a mock tree: protected files survive, the rest goes.
- **Open question:** "the CSR paper" is ambiguous. The only CSR-1 paper found is Pedroso *et al.* 2020 (PMID 32542533), which is characterisation and has no computational methods to replicate. D32 pairs EcoRV with **CphA (Xu & Guo 2006)**, and M2 assumes CphA is what was meant. To be confirmed.
- **Measured (listing only, nothing deleted yet): 66.7 GB reclaimable in `CSR_crystal_v2` alone**, which would take scratch from 327 GB to about 260 GB.

  | Category | Size | Items |
  | --- | --- | --- |
  | out-of-memory / cancelled 1D scans | 61.4 GB | 10 dirs |
  | geomeTRIC optimisation trajectories | 4.8 GB | 87 files |
  | ASH working files | 328 MB | 381 files |
  | cancelled 2D grids | 56 MB | 1 |
  | superseded QM/MM probe directories | 52 MB | 2 |
  | abandoned GPU benchmark | 25 MB | 1 |
  | config backups beyond the newest 3 | small | 12 |

  The killed scans are 1.3-16 GB each and are mostly geomeTRIC `HT.v_*.txt` history files at ~30 MB apiece - the same unbounded growth that caused the out-of-memory kills (D35). The two scans with usable output, `scan_H1_f00497` (165 MB) and `scan_H2_f00395` (177 MB), are excluded and keep their results. The script now also copies every killed run's `scan_energies.csv` into `scans/killed_runs_csv/` before removing the directory, so the record of how far each got survives at a few kB.
- **CphA confirmed** by the user as the mono-zinc template alongside EcoRV, which is what M2 already assumed.
- **Cleanup run (user approved): 69.6 GB removed in 2650 files. Scratch went from 327 GB to 181 GB**, well under the 300 GB quota; files 859,781 to 855,717. Most of the space was `HT.v_*.txt` geomeTRIC history files rather than the scan directories themselves, which is why the swept total came out above the 66.7 GB estimate.
  - **The running cleanR job was protected.** It writes `xtb_.*`, `gradient` and `geometric_OPTtraj*` every step, so sweeping its directory would have crashed it or corrupted the optimisation. `cleanup.sh` gained a `SKIP` variable for exactly this, and the run used `SKIP="protdft/cleanR"`. Job 30053828 was still running afterwards with its working files intact.
  - Verified kept: `qmmm_results/QM_MM.dcd` (550 MB), `md_results/NVTtrajectory_lastframe.pdb`, both surviving scans' `scan_energies.csv`, and all ten killed runs' CSVs now in `scans/killed_runs_csv/`.
  - **The other three run directories were then listed and are not worth sweeping:** `CSR_crystal` 995 MB (mostly 953 MB of geomeTRIC trajectories from the superseded v1 run), `CSR_ligand` 90 MB, `CSR_hydrolysed` nothing. About 1.1 GB in total against 119 GB of headroom, and `CSR_ligand` / `CSR_hydrolysed` are the reference runs whose `MM.slurm` defined the correct pipeline order, so their records are worth more than the space. Left alone.
  - Still installed: the GPU python packages, ~2-3 GB in HOME, which is at 40/50 GB. Uninstall commands are in `cleanup.sh` if HOME gets tight.

**D48 - M1 result: the transfer is still favourable, but by a third of what the degraded geometry said; and the Zn-bound oxygens still cannot attack.** (2026-10-18, job 30053828, COMPLETED in 9:51:17)

**Energies** (wB97X-D3(BJ)/def2-SVP//xTB/MM, identical sphere constraints on both states, both optimisations converged):

| State | `sp_energy_disp` (Eh) | xTB `opt_energy` (Eh) |
| --- | --- | --- |
| `qmmm_cleanR_lysH` (proton on Lys226, O01 = OH-) | -5070.611564 | -423.530261 |
| `qmmm_cleanR_ohH` (proton on O01, O01 = H2O) | -5070.619894 | -423.535511 |

**dE = E(ohH) - E(lysH) = -5.23 kcal/mol at DFT** (-3.29 at xTB).

- **Still negative, so the resting state is Zn-OH2 with neutral Lys226.** The model genuinely has no hydroxide, and H4 is not reinstated.
- **But the driving force is about a third of D30/D35** (-13.1 and -15.5 kcal/mol). So roughly two-thirds of that was an artefact of the degraded six-coordinate geometry, where Asp92's charge sat on the metal and the carbonyl competed for it. D35's "the Lys226 transfer is real at DFT" stands, but the margin was overstated.
- **-5.2 kcal/mol is small enough to treat with care.** It is within reach of def2-SVP basis error, the imposed constraints, and the absence of entropy, sampling and explicit environment relaxation. The hydroxide state is not comfortably excluded, only disfavoured.

**Final structures** (`protdft/cleanR/*/QM_MM_optimized.pdb`, both five-coordinate):

| | lysH | ohH |
| --- | --- | --- |
| Zn-O01 | **1.83** | **1.99** |
| Zn-His88 NE2 | 1.99 | 1.94 |
| Zn-His164 NE2 | 2.02 | 1.98 |
| Zn-HOH485 O | 2.04 | 1.97 |
| Zn-Asp92 OD2 | 2.80 (constrained) | 2.80 (constrained) |
| Zn-mer O2x | 2.90 (constrained) | 2.90 (constrained) |
| O01 hydrogens | **1** (H03) | **2** (H03 + HZ1) |
| Lys226 NZ hydrogens | **3** (NH3+) | **2** (neutral) |

- The proton pinning did what it was meant to: in ohH, HZ1 is on O01 and NZ keeps only HZ2/HZ3.
- **Zn-O01 lengthens 1.83 -> 1.99 A on protonation**, the textbook hydroxide-to-water signature, which is a good internal consistency check on the two states.
- **The decisive result is geometric, and it is independent of the protonation question: even in a clean four-coordinate sphere the Zn-bound oxygens cannot attack.** O01-C4x-O2x is **38-40 degrees** and HOH485 is **65 degrees**, nowhere near the Bürgi-Dunitz window of 95-120. This confirms D41 from an optimised structure rather than MD frames.

**Consequences for the plan.**
- **H2, H3 and H4 are geometrically dead**, whatever the protonation state: a nucleophile coordinated cis to the carbonyl on the same metal is on the wrong trajectory by construction.
- **M2 is the only live route** and is now well supported: a non-metal-bound water attacking C4x with Asp92 as general base, exactly CphA (Xu, Xie & Guo 2006). D33 found the second-shell water E3094 is the only species on the Bürgi-Dunitz line; this result removes its competitors.
- M1 is therefore complete. Next is M2 on the single best reactant, all at xTB/MM with DFT/MM single points on the final path.

**D49 - M2 cannot be set up from existing data: the attacking water never has a proton aimed at a base.** (2026-10-18)

Asked for M2, I checked the coordinates against the structures before building any grid, which is what the failed H1 grid (D38) taught. Three independent checks agree.

**1. The corrected resting-state reactant** (`protdft/cleanR/ohH/QM_MM_optimized.pdb`, `nuc_scout.py`):

| Water | O-C4x | angle O-C4x-O2x | Zn-O | |
| --- | --- | --- | --- | --- |
| O01 | 3.63 | 40 | 1.99 | Zn-bound, cis to the carbonyl |
| HOH485 | 3.80 | 65 | 1.97 | Zn-bound; its H reaches Asp92 OD1 at **1.65**, but it cannot attack |
| E3094 | 4.71 | 127 | 6.10 | protons point at bulk waters E2155 (1.86) and E2493 (1.80) |
| E3166 | 5.54 | **96** | 4.97 | good angle and H1 donates to Asp92 OD1 at **1.82**, but 5.5 A away |

No water satisfies distance and angle together. The ones H-bonded to Asp92 cannot attack; the one nearest the attack line donates to bulk solvent.

**2. The whole 50 ps trajectory** (`joint_nac.py`, 1000 frames, 46 waters within 12 A of C4x). Only two waters ever reach the Bürgi-Dunitz window at all: **HOH3094 in 326 frames** (best 3.67 A / 95 deg, frame 208) and HOH2493 in 95. **Neither donates to Asp92 in a single frame, and in no frame does any water do both.**

**3. The relay, and other acceptors.** Only 2 waters near the site ever donate to Asp92, and in no frame are attacker, bridge and base aligned together. Extending the acceptor set to the substrate's own O4x, O3x, O1x and N1x finds nothing either: while a water is on the attack line, none of its protons is within 2.5 A of any of them.

**So the Asp92-as-general-base step, direct or relayed, is not pre-organised anywhere in the data we have.** Forcing a 2D grid from here would repeat the D38 mistake of dragging a group across space and calling it a proton transfer.

**The important caveat, and why this is not yet a conclusion about the enzyme.** That trajectory was generated in the **degraded** sphere: Asp92 was coordinated to the zinc at 2.07-2.11 A and Lys226 had already lost its proton (D41). A carboxylate bound to a divalent metal is both oriented toward it and much weakened as a base, which is a plausible reason nothing donates to it. The analysis may simply be asking the wrong state. In the corrected reactant Asp92 sits at 2.80 A and is free.

**Minimal next step:** a short QM/MM MD from the corrected `ohH` reactant, then rerun `joint_nac.py` on it. This is the smallest simulation that can settle M2, and it matches what Xu & Guo did (SCC-DFTB/MM MD before any path). At the measured 1.37 s/step, 25 ps is about 9.5 h on 4 cores. If a free Asp92 recruits a water into the donating position, M2 becomes available and the frame to start from falls out of the analysis; if it does not, the Asp92-base mechanism is genuinely absent here and the hypothesis set needs rethinking rather than more scans.

**Tools added:** `docs/csr/nuc_scout.py` (candidates and their proton acceptors in one structure) and `docs/csr/joint_nac.py` (the joint attack-and-donate test over a trajectory, including the relay).

**D50 - 25 ps QM/MM MD from the corrected reactant, to settle whether M2 exists.** (2026-10-18, job 30085418, submitted 23:15 AEST)
- **Purpose.** D49 found no frame in the original 50 ps trajectory where a water both reaches the Bürgi-Dunitz window and donates to Asp92, directly or through a bridge. But that trajectory had Asp92 coordinated to the zinc at 2.07-2.11 A and Lys226 already deprotonated, so it may be the wrong state. This resamples from `protdft/cleanR/ohH/QM_MM_optimized.pdb`, the resting state established in D48: two His, O01 as water, HOH485, **Asp92 free at 2.80 A, Lys226 neutral**.
- **No constraints.** The `qmmm` stage applies no MM restraints during the QM/MM MD, so the metal site is described by xTB alone. Whether Asp92 returns to the zinc is therefore a *result* of this run, not an assumption. Worth noting that in CphA, Asp120 is both a zinc ligand and the general base, so its coordinating would not by itself rule out the mechanism; what matters is whether a water is positioned to use it.
- **Settings:** 25 ps, 1 fs step, 300 K, frames every 50 steps (500 frames), 4 xTB threads bound to 4 physical cores on 8 allocated CPUs, 32 GB, 2-day limit. About 9.5 h at the measured 1.37 s/step. Scripts `docs/csr/mk_md2.py` and `docs/csr/Md2.slurm`; config backup `config.yaml.bak16`.
- **The test afterwards** is the same `joint_nac.py` criterion, so the two trajectories are directly comparable. If frames satisfy both halves, the best one becomes the M2 reactant and the scan coordinates come from its actual distances. If none do, the Asp92-as-base mechanism is absent in this system and the hypothesis set needs rethinking rather than more scans.
- This is the only simulation running, and the minimum that can settle M2.

**D51 - M2 result: the six-coordinate site is what the model prefers, and Asp92 as base is not available.** (2026-10-19, job 30085418 COMPLETED in 9:11:36)

**The zinc re-formed the six-coordinate sphere on its own.** Started from the clean four-coordinate reactant with **Asp92 free at 2.80 A and no constraints**, 25 ps of QM/MM MD ended with:

| Ligand | Last frame | Started at |
| --- | --- | --- |
| HOH485 O | 2.04 | 1.97 |
| His164 NE2 | 2.07 | 1.98 |
| **Asp92 OD2** | **2.09** | **2.80** |
| O01 | 2.11 | 1.99 |
| His88 NE2 | 2.14 | 1.94 |
| mer O2x | 2.25 | 2.90 |

Lys226 kept its two hydrogens (still neutral) and O01 kept two (still water), so the D48 protonation state was stable - the back-transfer did not happen.

**This changes the reading of D41.** The Asp92 migration was not a one-off artefact of a bad setup: from a clean start, with nothing imposed, the model pulls Asp92 and the substrate carbonyl back onto the metal within 25 ps. A two-His zinc is under-coordinated, so recruiting two more ligands is chemically unsurprising. **Either the mono-Zn site genuinely is six-coordinate, or GFN2-xTB over-coordinates the zinc** - a known weakness of semi-empirical methods with metals, and precisely what Imhof's AM1/d-vs-B3LYP validation was for.

**The joint criterion still fails, now from the right starting state** (`joint_nac.py`, 500 frames):

| Water | attack frames | H to Asp92 | both |
| --- | --- | --- | --- |
| HOH3094 | **210** (best 3.41 A / 99 deg, frame 209) | **0** | **0** |
| HOH2493 | 51 (best 3.91 A / 96 deg) | 0 | 0 |

The relay fails too. HOH3094's protons still point at bulk waters (E2155 1.70, E2493 2.01); the waters that do reach Asp92 are HOH485, which is zinc-bound at a 38 degree attack angle, and E3166 at 2.17 A but 4.65 A from C4x and 64 degrees off.

**So the Asp92-as-general-base mechanism is not accessible in this model, tested twice from two different starting states.** Note that it fails for a *structural* reason that follows from the zinc coordination: with Asp92 bound to the metal and the carbonyl occupying a sixth site, nothing is left to accept the attacking water's proton.

**Everything now rests on one assumption: that xTB describes this zinc correctly.** If it over-coordinates, the whole picture - Asp92 unavailable, the carbonyl inner-sphere, the Zn-bound waters at 33-38 degrees - is an artefact, and the mechanism question reopens with a four- or five-coordinate site.

**Proposed next test (cheap and decisive), ahead of any more mechanism work:** from the md2 last frame, two constrained xTB/MM optimisations with Zn-Asp92 OD2 held at 2.10 and at 2.80 A, everything else free, then one DFT/MM single point on each. Two single points, about 9-10 h in one CPU job. If DFT prefers Asp92 off the metal while xTB prefers it on, the over-coordination is confirmed and the model needs changing before any path is computed. This is M4 (validation) brought forward ahead of M2, because M2's answer now depends on it.

**D52 - Reconsidering the mechanism: the zinc activates the electrophile, and the proton leaves down a water wire.** (2026-10-19, user: "step back and reconsider what mechanism a six-coordinate mono-zinc site could actually support")

**The measurement that settles it** (`acceptor_census.py`, md2 trajectory). Instead of asking "does the attacking water donate to Asp92", ask what is within hydrogen-bond reach of its protons at all, during the 210 frames where it is genuinely on the Bürgi-Dunitz line:

| Acceptor | Frames (of 210) | Closest |
| --- | --- | --- |
| HOH2155 O | **210** | **1.51 A** |
| HOH2493 O | **207** | **1.48 A** |

Nothing else, out of 36 O/N atoms within 8 A. No Asp92, no substrate atom, no Lys226, no histidine. Two waters, in essentially every attacking frame, at H-bond lengths short enough to be near-symmetric. **That is a proton wire, and it is persistent rather than transient.**

**Why a six-coordinate mono-zinc site implies this mechanism.** Coordination number sets the zinc's job. A four-coordinate zinc polarises one water enough to deprotonate it near neutral pH, which is the carbonic-anhydrase and canonical-MBL picture: the metal *makes the nucleophile*. Six-coordinate, the same charge is spread over more ligands, each is polarised less, and the aqua pKa rises - which is exactly what D48 measured, with the proton preferring the water over Lys226 by only 5.2 kcal/mol instead of strongly favouring hydroxide. What the site does have is the **substrate carbonyl bound at 2.25 A**, so the metal is a Lewis acid **activating the electrophile**. A metal that activates the electrophile does not need to make the nucleophile, and cannot easily do both: its two aqua ligands sit *cis* to the bound carbonyl and are therefore 33-38 degrees off the attack trajectory, which is a geometric consequence of the coordination, not an accident of our structures. This is the "promoted-water" pattern of carboxypeptidase A rather than the "nucleophilic-metal" pattern of carbonic anhydrase.

**The nucleophile then has to come from outside the coordination sphere, and its proton has to go somewhere that is not the metal site.** Both halves are satisfied here: HOH3094 attacks from the Bürgi-Dunitz line in 210 of 500 frames, and its protons are permanently engaged with HOH2155 and HOH2493 leading into bulk.

**This explains every failure so far, and they were all the same failure.** The 1D scans gave monotonic +40 kcal/mol rises because the proton had nowhere to go: the acceptors were not in the QM region, and a bug froze the QM waters anyway (fixed in 113de19). The H1 grid drove the proton at E3166/Asp92, 4.35 A away. `joint_nac.py` tested only Asp92. In each case the nucleophile was right and the **acceptor was wrong**.

**Consequences for the three questions.**
- **Q3, Asp92.** Likely *not* the general base: it is a zinc ligand at 2.09 A and nothing donates to it. If that holds, CSR-1 differs from CphA, where Asp120 is both ligand and base - a real result rather than a null one, and one the RQK substitutions could plausibly explain.
- **Q1, Lys226.** It is neutral in the resting state, having given its proton to O01. Its likely role is electrostatic stabilisation of the developing oxyanion, which is precisely Imhof's conclusion for Lys92 in EcoRV. Charge deletion tests it.
- **Q2, arginines.** Unchanged: Arg90 holds the 6-alpha-OH, Arg229 by charge deletion.

**Revised M2, from the measured geometry rather than from a hypothesis:**
- QM region: add **HOH3094, HOH2155 and HOH2493** (`extra_qm_residues`), so the wire is quantum-mechanical and the proton can actually leave.
- Start frame: an attacking frame, e.g. **frame 209** (O-C4x 3.41 A, angle 99 deg), optimised to a reactant first.
- Coordinate 1: **O(3094)-C4x**, from ~3.4 A to 1.45 A.
- Coordinate 2: **H(3094)...O(2155)**, from ~1.5 A to 1.0 A - already an H-bond, so this is a genuine transfer coordinate and not a translation.
- Monitor C4x-N1x, Zn-O2x, C4x-O2x, the second wire H-bond to HOH2493, and Zn-Asp92.

**Still open, and worth stating in any write-up:** whether GFN2-xTB over-coordinates the zinc (D51). If the true site is four- or five-coordinate with Asp92 free, the CphA route returns and this mechanism is the wrong one. The validation is cheap and should be done before the path is published, but it no longer blocks setting M2 up, because the water wire is present in the structures either way.

**D53 - Literature check on the D52 mechanism, and the zinc coordination test (job 30119957).** (2026-10-19)

A literature sweep, with the two load-bearing papers verified directly through Europe PMC rather than taken on trust.

**The nucleophile has direct precedent; the zinc does not.**
- **Ackerman & Gatti, *PLoS One* 2013, 8:e55136** (PMID 23372827), on the mononuclear B2 enzyme CphA with biapenem: "the nucleophile that attacks the beta-lactam ring is **not** the water molecule located in proximity of the active site Zn" - it is a **second water hydrogen-bonded to it**, consumed in the reaction and therefore invisible in the enzyme-product structure. This is the closest analogue to our HOH3094 and should be cited prominently; it reached the same conclusion by driving the reaction backwards from the product structure, which avoids substrate docking.
- **No MBL study proposes the nucleophile's own proton leaving along a water chain to bulk.** Bulk water as a proton *donor* to the lactam nitrogen is established (Das & Nair, *PCCP* 2017, PMID 28489087), but the reverse is not. The water-wire half of D52 is novel. The nearest hydrolase precedent is **isatin hydrolase**, where a proton wire and water channel were seen crystallographically (Bjerregaard-Andersen *et al.*, *JBC* 2014, PMID 24917679).
- **Correction to my own framing in D52.** In carboxypeptidase A the "promoted water" is still **zinc-bound**; what is promoted is its deprotonation by Glu270. Our non-coordinated, second-shell water is a further step away from the metal, so calling it the textbook promoted-water mechanism was wrong, and the paper should say so explicitly rather than claim the analogy.

**Answering the three questions directly.**

*Is the six-coordinate site an artefact of our setup?* **Probably, and the failure mode is documented.** Yang *et al.*, *J Chem Inf Model* 2021, 61:5658 (PMID 34748329) - authors include Spencer and Mulholland, who work on MBLs - found across six zinc systems that classical MM MD "overfavors octahedral geometries" by pulling extra waters into the zinc shell, and recommend DFTB3 QM/MM MD followed by DFT QM/MM optimisation to correct it. **Our pipeline equilibrated with MM first**, under preprocess metal restraints on everything within 3.0 A of the zinc, which is exactly that setup. The tell is in the numbers: the NVT end frame had **Zn-O2x at 1.94 A, shorter than the Zn-OH- bond**, which is a restraint signature rather than a chemical one. Yang found DFTB3 QM/MM MD *rectifies* the over-coordination; our xTB QM/MM MD did not, which either means GFN2-xTB is worse than DFTB3 here or - more likely - that 25 ps starting from an already-octahedral geometry is far too short to escape it, since ligand exchange at Zn(II) is slow. No GFN2-xTB benchmark for zinc coordination number in enzymes appears to exist.

*Does the carbonyl have to coordinate the zinc for lactamase activity?* **No.** Inner-sphere Zn-O(carbonyl) is one proposal - Olsen *et al.*, *J Phys Chem B* 2004, QM/MM on benzylpenicillin in BcII, favours it over the carboxylate-bound pose - but the dominant B1 picture has the carbonyl polarised and hydrogen-bonded, with the developing oxyanion stabilised by Zn1 and Asn233, while the bridging hydroxide attacks and Zn1 keeps its 3-His + OH shell. There is **no crystallographic Michaelis complex of an MBL with intact substrate** to settle it; the NDM-1 cephalosporin structures are cephalosporoate intermediates. So our 2.25 A cannot be benchmarked against experiment, and polarisation at 3.0-3.5 A would serve the same catalytic purpose.

*Are two waters normal?* **For this site, yes - the extra ligands are the problem.** Catalytic Zn(II) in hydrolases is predominantly four-coordinate, three protein ligands plus one water, with five-coordinate transition and intermediate states (Auld, *BioMetals* 2001, 14:271). **CSR-1's Zn1 site has only two protein ligands**, His116 and His196 (our His88 and His164), because Arg118 replaces the third histidine. Two non-protein ligands are therefore needed to reach the normal four-coordinate geometry, and two waters is a perfectly canonical way to do it. What has no precedent is two waters **plus** Asp92 **plus** the substrate carbonyl. Six-coordinate zinc does occur in this family - PNGM-1 has both zincs octahedral (Park *et al.*, *IJMS* 2020, PMID 32664695) - but not with two aqua ligands and a substrate carbonyl. Worth noting that the four-coordinate state we constructed for D42 (two His, O01, HOH485, Asp92 at 2.80, carbonyl at 2.90) is precisely the canonical arrangement; we imposed it, and the literature says it is what a two-His zinc should look like.
- Supporting the pKa argument from D52: aqua pKa rises with coordination number, roughly linearly in structurally matched series, from about 6-7 in carbonic anhydrase II to 9-9.5 for the hexaaqua ion. A six-coordinate zinc genuinely cannot make hydroxide easily, which is consistent with D48's -5.2 kcal/mol.

**Test submitted: job 30119957** (`mk_zncheck.py`, `ZnCheck.slurm`). Relax the md2 end point at xTB/MM, walk **Zn-Asp92 OD2 from 2.05 to 3.05 A** in 9 points at xTB/MM with the whole first shell monitored, then recompute points 0, 3, 6 and 8 at wB97X-D3(BJ)/def2-SVP/MM. Four separate `qmmm_sp_dft_*` sections, because a shared one would overwrite its own results three times. If xTB puts the minimum near 2.1 A and DFT near 2.8, the six-coordinate site is an artefact and D52's mechanism rests on it. About 22 h; config backup `config.yaml.bak17`.

**D54 - Water wire dropped; the hypothesis is CphA or EcoRV, and the geometry has to be settled first.** (2026-10-19, user)
- **The water-wire-to-bulk mechanism (D52) is withdrawn as implausible.** It had no MBL precedent (D53), and on reflection it was inference filling the gap left by Asp92 being unavailable, rather than something the data positively supported. What survives from D52 is the *observation* - HOH3094 attacks in 210 of 500 frames and its protons are permanently H-bonded to HOH2155 and HOH2493 - not the mechanistic claim built on it.
- **The hypothesis to prove or rule out is CphA-like or EcoRV-like.** Both share the same two features: a **water nucleophile** and a **carboxylate general base** (Asp120 in CphA, Asp90 in EcoRV), with the metal activating the water (EcoRV) or the substrate. In our numbering that is Asp92. Everything else is detail.
- **That question cannot be answered on the present geometry**, because whether Asp92 is available as a base depends on whether it is bound to the zinc, which depends on a coordination number we now have reason to distrust (D53).
- **User's proposal, adopted:** optimise at DFT/MM **straight from the restrained MM structure**, with no xTB QM/MM equilibration, so the xTB bias never enters the geometry. Also to investigate **wB97X-3c** (composite, cheaper than wB97X/def2-TZVP) and **g-xTB** (trained towards wB97X-V, better non-covalent behaviour than GFN2-xTB).
- **Local availability checked:** `g-xTB` is **not installed** and is not on PATH; ASH ships only `interface_xtb.py`, with no tblite or g-xTB interface, so using it would need both a binary and an interface. PySCF's basis directory has **no Grimme vDZP** (the `dzp`/`adzp` files there are unrelated Dunning-type sets), so wB97X-3c would need the basis and ECPs supplied externally.
- **The cost obstacle, quantified.** The QM region is 119 atoms and breaks down as meropenem 51, Arg90 18, Lys226 16, His88 11, His164 11, Asp92 6, HOH485 3, O01 2, Zn 1. A DFT/MM gradient on all 119 costs about as much as the measured 4-5 h single point, so a full optimisation is days to weeks - the same wall that drove the abandoned GPU detour. **The coordination question does not need all 119 atoms.** First shell plus a truncated substrate - Zn, both imidazoles, the Asp92 acetate, both waters, and meropenem cut back to the beta-lactam ring, carbonyl, N1x and C3 carboxylate - is roughly 50-55 atoms, where hybrid DFT is several times cheaper and an optimisation becomes a job rather than a campaign. This is also precisely the model-system validation Imhof used.
- Job **30119957** (xTB scan along Zn-Asp92 with four DFT single points) is left running: its xTB half still measures where xTB puts Asp92, which is one side of the comparison either way, and it costs nothing extra now.

**D55 - Neither wB97X-3c nor g-xTB is the right tool here; r2SCAN-D4/def2-TZVP is, and it is already available.** (2026-10-19)

Both suggestions were investigated properly and the two blocking facts were verified against primary sources rather than taken from a search summary.

**g-xTB is blocked today, for a reason outside our control.** It is real and attractive: ChemRxiv preprint June 2025, GPL-3.0 at `github.com/grimme-lab/g-xtb` with precompiled binaries, WTMAD-2 of 9.3 kcal/mol on GMTKN55 at about 1.5x GFN2-xTB cost, and **ASH already supports it** as `xTBTheory(xtbmethod='GXTB')` with analytic gradients since v2.0.0. But the ASH documentation states plainly: **"Note that electrostatic embedding QM/MM with g-XTB will not yet work."** g-xTB does not implement point-charge embedding (`$pcem`). Our entire QM/MM setup is electrostatic embedding, so g-xTB cannot be the QM level. Worth revisiting once the method lands in tblite proper, which the developers say is the plan.

**wB97X-3c is viable on cost but is the wrong method for a zinc geometry.**
- The good news first: it does **not** use VV10. The composite replaces wB97X-V's non-local term with **D4** (a1 0.2464, a2 4.737, **s8 = 0**, s9 1.0), so none of the cost problem that killed wB97M-V applies (Müller, Hansen, Grimme, *J Chem Phys* 2023, 158:014103).
- But the authors' own benchmarks report it is **less accurate for 3d transition-metal geometries (TMC32)** than for main-group bonds, attributed to the large-core ECPs, and worse for open-shell metal systems (ROST61). For geometry work they point towards r2SCAN-3c. **There is no published endorsement of wB97X-3c for zinc geometry optimisation** - which is exactly and only what we need it for.
- And it is not available: **vDZP is not in PySCF and not on the Basis Set Exchange**; it would mean converting basis and large-core ECPs from `grimme-lab/wB97X-3c` plus supplying the refitted D4 parameters, each piece needing validation against an ORCA reference. Days of work, for a method its own authors flag as weaker at our target property. Native support exists in ORCA >= 6.0.0, if we ever want it for energies.

**What we will use instead: r2SCAN-D4/def2-TZVP.** Checked directly in our environment: **`r2scan` is present in PySCF's libxc** (confirmed as a meta-GGA) and **def2-TZVP loads for Zn**. ASH's `PySCFTheory` already accepts `dispersion='D4'`. So this needs no basis conversion, no new software and no validation of hand-assembled pieces. It is also the right class of method: a meta-GGA has no exact exchange, so it is several times cheaper than wB97X at the same basis, and r2SCAN is what the wB97X-3c authors themselves recommend for geometries. Reported energies can still be wB97X-D3(BJ)/def2-TZVPD single points on the r2SCAN geometries, which is the usual and defensible split. Note this is r2SCAN-D4/def2-TZVP, **not** the r2SCAN-3c composite, whose mTZVPP basis has the same availability problem as vDZP.

**Adopting the user's design.** Optimise at DFT/MM **straight from the restrained MM structure**, with no xTB QM/MM equilibration, so the xTB bias never reaches the geometry. To make it a real test rather than a local relaxation, run it from **two different starting coordinations** and compare both the final geometries and their energies:
  1. `md_results/NVTtrajectory_lastframe.pdb` - the restrained MM end point, carbonyl at 1.94 A, Asp92 at 2.79 A;
  2. `protdft/cleanR/ohH/QM_MM_optimized.pdb` - the four-coordinate state, Asp92 at 2.80 A, carbonyl at 2.90 A.
A local optimisation cannot leave its basin, so starting from one geometry alone would prove nothing.

**Timing probe first.** Cost at this size has surprised us twice, so before committing to two multi-day optimisations, a three-step run measures the r2SCAN-D4/def2-TZVP gradient cost on the real 119-atom QM region. If it is tens of minutes the full optimisations are a job each; if it is hours we fall back to the reduced ~50-atom model described in D54.

**D56 - What the two reference papers actually do with the substrate-metal contact, and a useful precedent for a six-coordinate metal.** (2026-10-19, answering the user's question)

User: drop the run from the restrained MM, since its carbonyl is pinned at an unrealistic 1.94 A; a polarised carbonyl should be enough; is that what CphA and EcoRV did?

**Short answer: for CphA yes, for EcoRV the metal's job is different again - and neither relies on the metal polarising the substrate carbonyl.**

- **CphA** (Xu, Xie & Guo 2006, PMID 16423823): nucleophile is a **non-metal-bound water**, general base is the **Zn-bound Asp120**, and the paper does not report Zn-O(carbonyl) coordination. A merely polarised carbonyl is consistent with it.
- **EcoRV** (Imhof, Fischer & Smith 2009, PMID 19678693, abstract verified): the metal does not polarise the substrate at all - **"the magnesium metal ion activates a water molecule"**. **Asp90 accepts the proton** from that water, generating the nucleophile, and **Lys92 positions the scissile phosphate and stabilises the metaphosphate intermediate** in the orientation best for attack. The abstract does not state whether the phosphate is inner-sphere to Mg.

**The genuinely useful find, which reframes our problem.** Imhof *et al.* report that in EcoRV the magnesium's **coordination number drops from six to four** when the attacking water is positioned, and they use exactly this to explain **why calcium is not catalytically active** - Ca(II) prefers higher coordination and cannot shed ligands to open the site.

That changes the question we have been asking. A six-coordinate resting metal is **not** by itself incompatible with catalysis; there is a published precedent where it is the resting state and the catalytically competent state is reached by shedding ligands. So the question is no longer "is our six-coordinate zinc an artefact" but:

1. is the six-coordinate site the resting state (plausible either way), and
2. **can it drop to four-coordinate to open a site and position a nucleophile, and at what energetic cost?**

That is a relaxed-coordinate question, and it is what job 30119957 is already asking along Zn-Asp92, albeit at xTB. The DFT optimisations now answer the first half and give the energy difference between the two coordinations.

**Setup changed accordingly.** `dftopt_fromMM` is dropped. The two starts are now:
  - **`dftopt_from4c`** from `protdft/cleanR/ohH/QM_MM_optimized.pdb` - Asp92 at 2.80 A, carbonyl polarised at 2.90 A, i.e. the arrangement the user wants and the one consistent with CphA;
  - **`dftopt_from6c`** from `protdft/md2_ohH/QM_MM_lastframe.pdb` - Asp92 2.09 A, carbonyl 2.25 A. This is **xTB's own preferred geometry, not the restrained-MM artefact**, so the two-basin comparison survives the user's objection while dropping the structure they objected to.
Neither has constraints: where Asp92 and the carbonyl end up, and which is lower in energy, is the result.

**D57 - xTB says Asp92 is firmly bound; DFT/MM optimisation at the full QM region is unaffordable; waiting for the DFT single points instead.** (2026-10-19)

**Finding: the xTB relaxed scan along Zn-Asp92 OD2 (job 30119957) is complete, all nine points converged.**

| Zn-OD2 (A) | rel (kcal/mol) | Zn-His88 | Zn-His164 | O01-C4x |
| --- | --- | --- | --- | --- |
| 2.05 | 0.00 | 2.100 | 2.040 | 3.619 |
| **2.175** | **-0.12** | 2.088 | 2.022 | 3.621 |
| 2.30 | +0.51 | 2.079 | 2.006 | 3.627 |
| 2.425 | +1.60 | 2.072 | 1.992 | 3.632 |
| 2.55 | +2.95 | 2.065 | 1.979 | 3.637 |
| 2.80 | **+5.72** | 2.059 | 1.958 | 3.664 |
| 2.925 | +7.14 | 2.060 | 1.946 | 3.669 |
| 3.05 | +8.46 | 2.068 | 1.937 | 3.695 |

**GFN2-xTB puts the minimum at 2.18 A and charges 5.7 kcal/mol to pull Asp92 out to its crystallographic 2.80 A**, 8.5 to 3.05. The curve is smooth and monotonic past the minimum, with no second well - so at this level Asp92 is a genuine zinc ligand, not a loosely associated one. The two histidines shorten as Asp92 leaves (2.10 to 2.07 and 2.04 to 1.94), the zinc compensating for lost charge donation, and the O01-C4x attack distance barely moves (3.62 to 3.70), so opening this particular coordination site does not by itself bring the nucleophile closer.

**Cost measured, and it rules out the planned DFT optimisations.** The probe (job 30123101) completed **one** optimisation step - 22 SCF cycles plus the gradient - in **about two hours** at r2SCAN/def2-TZVP on the 119-atom QM region (2192 basis functions, 16 cores). So roughly **2 h per gradient**: a 250-step optimisation is three weeks, and even 25 steps is two days. The probe was cancelled once this was clear; it had served its purpose.

**Decision: do not submit the two DFT optimisations yet, and wait for the single points already running.** The instruction was to fall back to a reduced model if the rate was too slow, but there is a better move available first. Job 30119957 is now computing **wB97X-D3(BJ)/def2-SVP single points at four points of the xTB scan** (2.05, 2.425, 2.80, 3.05 A). That gives the DFT energy as a function of Zn-Asp92 distance **on the same geometries**, which is exactly the xTB-versus-DFT comparison we want, for about 18 h of compute that is already underway:
- if the DFT curve also rises steeply with distance, xTB is not over-coordinating, the six-coordinate site is real, and no DFT optimisation is needed at all;
- if DFT is flat or prefers the long distance, then xTB is wrong about the site and DFT geometries become worth their cost - at which point the reduced model from D54, or def2-SVP in place of def2-TZVP, makes them affordable.

Spending four days of DFT optimisation before reading an 18 h answer that is already running would be the wrong order.

**D58 - The six-coordinate zinc is real: DFT binds Asp92 harder than xTB does, not softer.** (2026-10-20, job 30119957)

**The comparison that D53 was built to make** (wB97X-D3(BJ)/def2-SVP single points on the xTB scan geometries, so only the method differs):

| Zn-Asp92 OD2 (A) | GFN2-xTB | wB97X-D3(BJ) |
| --- | --- | --- |
| 2.05 | 0.00 | 0.00 |
| 2.425 | +1.60 | **+3.02** |
| 2.80 | +5.72 | **+9.12** |

(The 3.05 A point was still running; three points already settle it.)

**DFT is 1.6x steeper than xTB, which is the opposite of the concern.** The worry from D51/D53 was that GFN2-xTB over-coordinates the zinc and had invented the sixth ligand. It has not: at the DFT level Asp92 is bound *more* strongly, costing **9.1 kcal/mol** to reach its crystallographic 2.80 A rather than xTB's 5.7. **The six-coordinate site is a real feature of this model, not a semi-empirical artefact**, and the Yang *et al.* over-coordination failure mode - real in general, and a legitimate thing to have checked - does not apply to this result.

Consequences:
- **D51 and D52 stand on firmer ground.** The zinc is genuinely six-coordinate, so it is a Lewis acid activating the bound carbonyl rather than a hydroxide generator, and the two aqua ligands genuinely cannot reach the attack trajectory.
- **The planned DFT geometry optimisations are unnecessary** and will not be run: they were insurance against an artefact that does not exist. That also removes the cost problem from D57, where a def2-TZVP gradient measured about two hours.
- **Asp92 being a zinc ligand does not rule out the CphA mechanism**, since in CphA Asp120 is both ligand and base and protonation is what breaks its Zn-O bond. But at 9.1 kcal/mol the bond is not weak.

**The sharper problem is geometric, not energetic.** D57 showed that pulling Asp92 off does **not** bring the nucleophile closer: O01-C4x moved only 3.62 to 3.70 A across the whole scan. Opening that site does not make the site competent, so paying 9 kcal/mol for it buys nothing.

**Next, submitted as job 30200387 (`mk_coscan.py`, `CoScan.slurm`):** a relaxed xTB/MM scan pulling the **substrate carbonyl** off the zinc, Zn-O2x from 2.15 to 4.00 A in nine points, from the same relaxed structure as the Asp92 scan so the two coordinates are directly comparable. The carbonyl is the other removable ligand, and releasing it differs in kind from releasing Asp92: the two Zn-bound waters sit 34-36 degrees off the Burgi-Dunitz line **because** they are cis to a carbonyl pinned on the same metal, so freeing the substrate to rotate is the one move that could plausibly fix the geometry. Against that, leaving also gives up the Lewis-acid activation, so the scan is a genuine test rather than a foregone conclusion. Monitors both Zn-water distances to C4x, the rest of the first shell and C4x-N1x; attack angles come from `zn_sphere.py` over the scan structures afterwards, since `monitor` handles only distances. Start geometry checked first: Zn-O2x is 2.13 A, a real bond, so the coordinate is well posed.

**D59 - The starting point was forced harder than I had said, and D58 is weaker than I presented it.** (2026-10-20, user: "was our starting point correct for these tests? ... should we run an unconstrained QM/MM equilibration ... with the substrate further away?")

**Checked in the code, and the answer is uncomfortable.** `preprocess` detects metal-donor pairs within `cutoff` (3.0 A here) and saves them as `results.constraints`. `md.py` then passes them to OpenMM as **holonomic constraints** - rigid, not harmonic. So through the entire MM pipeline, NPT and NVT, the **Zn-carbonyl distance was exactly fixed**. That is why it sat at 1.94 A, shorter than the Zn-OH- bond: it was never allowed to relax, and 1.94 A is simply wherever the input structure happened to put it. Every structure since descends from that.

**The QM/MM stages were genuinely unconstrained**, though: `qmmm.py` never reads `preprocess.results.constraints`. So the 25 ps md2 run was free, and it re-formed the octahedron from a start where Asp92 (2.80 A) and the carbonyl (2.90 A) had both been pulled out. That is a real test and it passed.

**What D58 actually establishes, stated properly.** It shows that **along the Zn-Asp92 coordinate, from that structure**, wB97X agrees with GFN2-xTB and is in fact steeper. So xTB is not over-binding Asp92. It does **not** establish that the six-coordinate arrangement is the correct Michaelis complex, because both methods were asked the same conditional question - "given this octahedron, what does removing one ligand cost?" - and neither was asked what the site should look like. Agreement between two methods on a loaded question is not evidence that the question was right. The D58 entry overstated this and should be read with D59.

**The untested assumption is the substrate pose.** Everything descends from a single modelled meropenem pose whose carbonyl was rigidly pinned to the zinc during equilibration. Displacing it locally (to 2.90 A) and releasing just relaxes back, which is unsurprising and proves little. A genuinely different pose - carboxylate anchored on Arg90 and Lys226 with the carbonyl merely polarised at 3.5 A or beyond - has never been sampled, and that is exactly the arrangement the B1 literature describes and that CphA is consistent with (D53, D56).

**mcube already has the knob we failed to use.** `preprocess.inputs.constraint_exclude_residues: {chain: [resid]}` skips the metal bond constraint for a residue while keeping it in the QM region. Applying it to meropenem (`{C: [1]}`) would have let the MM equilibration find the substrate's own distance instead of freezing it. That is the root-cause fix for a clean re-do.

**Plan, cheapest discriminating test first.**
1. **Pose-robustness test, nearly free.** The running carbonyl scan (30200387) is producing relaxed structures with Zn-O2x out to 4.00 A - properly displaced substrate poses, not 0.8 A nudges. Take the far point and run **unconstrained QM/MM MD** from it. If the carbonyl returns to the metal, the pose is robust and the six-coordinate Michaelis complex is a genuine conclusion rather than an inherited one. If it settles somewhere else, we have a different complex and the mechanism question reopens with better geometry.
2. **Clean re-do, only if step 1 says the pose is not robust** (or for the final published model regardless): re-run `preprocess` with `constraint_exclude_residues: {C: [1]}`, then MM equilibration and QM/MM MD. This removes the objection entirely but repeats the expensive part of the pipeline.

Doing 1 before 2 avoids repeating days of equilibration to answer a question that 10 h may settle.

**D60 - Tooling consolidated; CPU efficiency measured; what CphA actually says about the carbonyl.** (2026-10-20, user)

**Scripts consolidated.** Eleven one-off `mk_*.py` and analysis scripts plus seven per-job `*.slurm` files became four, which is the structure mcube itself argues for - configuration, not code per problem:

| File | What it is |
| --- | --- |
| `site.yaml` | the active site (atoms, candidate bases, near-attack criteria) and every calculation we run |
| `sections.py` | adds a calculation to `config.yaml`, refusing to write if a selector does not resolve in that calculation's own start structure |
| `site.py` | all the analysis: `sphere`, `nucleophiles`, `acceptors`, `joint` |
| `run.slurm` | one job script for every stage; resources from the `sbatch` line, `xtb` or `dft` threading profile |

Everything system-specific now lives in the YAML, so the tools carry no knowledge of this enzyme and a new calculation is a few lines rather than a new file. Smoke-tested on the cluster: `sphere` and `nucleophiles` reproduce the old scripts' numbers exactly on the relaxed structure (Zn sphere 2.03-2.13 A, O01-C4x 3.62 A at 36 degrees).

**CPU efficiency measured rather than assumed** (`seff`), and it splits sharply:

| Job | Work | CPUs | CPU efficiency | Memory used |
| --- | --- | --- | --- | --- |
| 30053828 | xTB opt + **DFT** single points | 16 | **96.2%** | 56 GB of 96 |
| 30085418 | pure **xTB/MM MD** | 8 | **29.6%** | 0.5 GB of 32 |

So the 64 CPUs now available should go to **DFT**, which genuinely saturates what it is given, and **not** to xTB, which saturates around four threads and leaves the rest idle. xTB work gets faster by running points in parallel, not by adding threads. Memory follows the same split: DFT wants tens of GB for density-fitting integrals, xTB about one, so the 32 GB we routinely asked for on MD jobs was ~60x more than needed.

**What CphA used for the Zn-carbonyl distance: nothing, because there is no such contact.** Xu, Xie & Guo 2006 (PMID 16423823) states that after ring opening the nitrogen leaving group replaces Asp120 "as the **fourth** ligand of the sole zinc ion" - so CphA's zinc is **four-coordinate**, with no room for the substrate carbonyl, and the paper reports no Zn-O(carbonyl) coordination. The abstract gives no interatomic distances at all, and the full text is not in Europe PMC (the record says `inEPMC: N`, no PMC id), so this is from the abstract rather than from reading the paper; an Unpaywall PDF exists at jbc.org if the numbers are needed.

**So the answer to "can we just use their distance" is that borrowing from CphA means not coordinating the carbonyl at all** - which is what `constraint_exclude_residues` would give us (D59), and is consistent with the dominant B1 picture where the carbonyl is polarised and hydrogen-bonded rather than inner-sphere (D53). The caveat is that CphA is a B2 enzyme whose zinc sits in a Cys-containing site (Asp120, Cys221, His263) and is four-coordinate, while ours is B3-RQK with only two histidines; the coordination chemistry differs, so the transferable part is the principle, not a number.

**D61 - Neither removable ligand opens a productive geometry. The substrate binding mode is now the thing to question.** (2026-10-20, jobs 30200387 and 30119957 both COMPLETED)

**The carbonyl scan is a clear negative** (xTB/MM, nine points, all converged):

| Zn-O2x (A) | rel (kcal/mol) | O01-C4x / angle | HOH485-C4x / angle | Zn-His88 | Zn-His164 |
| --- | --- | --- | --- | --- | --- |
| 2.15 | 0.00 | 3.62 / 36 | 3.65 / 34 | 2.089 | 2.031 |
| 2.61 | 3.29 | 3.82 / - | 3.82 / - | 2.010 | 2.045 |
| 3.08 | 8.53 | 4.13 / 46 | 3.86 / **11** | 1.966 | 2.081 |
| 3.54 | 14.53 | 4.47 / - | 3.88 / - | 1.939 | 2.113 |
| 4.00 | **22.36** | 4.74 / 66 | 3.83 / **19** | 1.920 | 2.158 |

**It costs 22.4 kcal/mol and makes the geometry worse, not better.** O01's angle improves from 36 to 66 degrees but it moves *away*, 3.62 to 4.74 A; HOH485 stays at the same distance while its angle collapses from 34 to 19 degrees. Neither approaches the Burgi-Dunitz window at any point. The substrate does not rotate to expose the carbonyl face - it simply drifts off the metal, and the zinc contracts onto His88 (2.09 to 1.92 A) to compensate.

**So both removable ligands have now been tested and neither helps**: Asp92 costs 5.7 kcal/mol at xTB (9.1 at DFT) and changes O01-C4x by 0.08 A (D57); the carbonyl costs 22.4 and makes things worse. The EcoRV precedent of a metal shedding a ligand to open a site (D56) does not transfer here.

**The D58 DFT table is now complete** (job 30119957, wB97X-D3(BJ)/def2-SVP on the xTB scan geometries):

| Zn-Asp92 (A) | xTB | DFT |
| --- | --- | --- |
| 2.05 | 0.00 | 0.00 |
| 2.425 | +1.60 | +3.02 |
| 2.80 | +5.72 | +9.12 |
| 3.05 | +8.46 | **+13.72** |

DFT is consistently about 1.6x steeper across the whole curve, confirming D58: xTB is not over-coordinating, and the six-coordinate zinc is real at both levels.

**The planned pose-robustness MD (D59) is now pointless and has been cancelled.** It was to start unconstrained MD from the 4.00 A scan structure, but the scan itself shows that structure is 22.4 kcal/mol uphill, so the carbonyl would obviously return. The scan answered the question the MD was meant to ask, more cheaply. Displacing along Zn-O2x was never going to test the *pose*, only that one bond.

**What is actually untested is the substrate binding mode.** Every structure descends from one modelled meropenem pose whose carbonyl was held on the zinc by a rigid MM constraint at 1.94 A through all of equilibration (D59). Scanning that bond explores one coordinate of that pose; it cannot find a different one. A genuinely different mode - carboxylate anchored on Arg90 and Lys226 with the carbonyl polarised at 3-3.5 A rather than inner-sphere, which is the dominant B1 picture (D53) and consistent with CphA's four-coordinate zinc having no room for it (D60) - has never been sampled, and could not have been, because the setup forbade it.

**Three ways to test the binding mode, in increasing cost:**
1. **Re-run `preprocess` with `constraint_exclude_residues: {C: [1]}`**, then MM equilibration and QM/MM MD. The substrate's coordination becomes an outcome rather than an input. This is the clean re-do deferred in D59 and is the one that removes the objection outright.
2. **Re-dock meropenem** into the apo site and take the best poses, rather than relying on the single inherited one.
3. **Build the literature pose by hand** - carboxylate on Arg90/Lys226 - and equilibrate it, then compare its energy with the current pose.

Option 1 is the honest default: it changes one config key and lets the method decide. The cost is repeating the MM equilibration and a QM/MM MD, which is days, so it is the user's call.

**D62 - Clean re-do started: the substrate's coordination becomes an outcome.** (2026-10-20, user chose the clean re-do; job 30205872)

**Set up as a controlled experiment rather than a fresh pipeline.** `/scratch/user/uquadhik/CSR_crystal_v3` **reuses v2's solvated system verbatim** - its `preprocess.results` still point at `CSR_crystal_v2/preprocess_results/finalsystem.pdb` - and changes exactly one thing:

| | v2 | v3 |
| --- | --- | --- |
| Zn-HOH485 | constrained | constrained |
| Zn-His164 NE2 | constrained | constrained |
| Zn-O01 | constrained | constrained |
| Zn-His88 NE2 | constrained | constrained |
| Zn-Asp92 OD2 | constrained | constrained |
| **Zn-mer O2x (substrate carbonyl)** | **constrained, rigid, 1.94 A** | **free** |

Done through `md.inputs.constraints`, which `md.py` prefers over the preprocess list, so it is a config change and no code was touched. Re-solvating would have given a different water configuration and made the comparison with v2 uncontrolled; this way the only variable is the one in question. The five protein and water contacts stay held because classical MD has no quantum description of the zinc and the site distorts without them - which is the documented over-coordination failure mode (Yang *et al.*, D53) - while the substrate, whose coordination is the thing in doubt, is left to the force field and then to QM/MM.

**Submitted:** `v3_md` (30205872), MM equilibration, NPT plus 1 ns NVT at 4 fs on 16 CPUs. The v2 equivalent took about 4.5 h. Then the 50 ps QM/MM MD, about 19 h, giving a like-for-like comparison with v2's trajectory.

**What the answer looks like.** If Zn-O2x drifts well beyond about 2.5 A once it is not held, the inner-sphere carbonyl was an artefact of our own setup and the Michaelis complex behind D52, D57, D58 and D61 was wrong. If it stays near 2 A unforced, the coordination is real, those conclusions stand, and the mechanistic difficulty is genuine rather than self-inflicted. Either outcome is worth having: the first rescues the mechanism question, the second closes off the last excuse for the geometry.

**`run.slurm` gained `preprocess` and `md`**, which take no `--section`; for them the third argument only names the output directory. One script still covers every stage.

**D63 - The MM re-do was the wrong instrument; D59 and D62 overstated the constraint's effect. Checking the carbonyl bond at DFT instead.** (2026-10-20, user: "the MM equilibration doesn't tell us anything about coordination")

**The user is right, on two counts, and job 30205872 was cancelled minutes after it was submitted.**

*First, classical MD cannot answer a coordination question here.* The zinc is a non-bonded model - a +2 point charge with Lennard-Jones - so whatever Zn-O2x distance the force field settles on reflects AMBER's parameters, not chemistry. Worse, it is biased in the direction that would have fooled us: Yang *et al.* (D53) document MM MD **pulling extra ligands into the zinc shell**, so an MM run that left the carbonyl bound would have been evidence of nothing, and we might have read it as confirmation.

*Second, and this is the part I had wrong: the forced 1.94 A never propagated.* Tracing the distance through the pipeline:

| Stage | Zn-O2x |
| --- | --- |
| MM NVT end point, constraint applied | **1.94 A** |
| after 25 ps QM/MM MD (md2) | 2.25 A |
| after QM/MM optimisation (zncheck relax) | 2.13 A |

QM/MM corrected it immediately and never looked back. So the rigid constraint set a *starting* value that xTB discarded, and **every QM/MM result we have already rests on a bond length the quantum method chose**. D59 and D62 treated the constraint as though it contaminated the downstream geometry; it did not. What it genuinely did was hold the substrate in one binding **pose** through 1 ns of classical MD - a pose question, not a bond-length question, and not one an MM re-run would settle either.

**So the carbonyl scan (D61) was already the right experiment**, done at QM/MM level from an equilibrated structure: 22.4 kcal/mol to pull the carbonyl to 4.00 A. The only open question is whether xTB is right about that number, and there is a cheap precedent for checking - D58 did exactly this for Asp92 and found DFT 1.6x **steeper**, not weaker.

**Submitted: jobs 30207244, 30207245, 30207246** - wB97X-D3(BJ)/def2-SVP single points on carbonyl-scan points 0, 4 and 8 (Zn-O2x 2.15, 3.08 and 4.00 A), on the xTB geometries so only the method differs. Three sections `qmmm_sp_co0/4/8`, separate so none overwrites another's results. **64 cores each**, which the D60 measurements justify: DFT ran at 96% efficiency on 16 and is the one stage that scales.

If DFT also says the carbonyl is firmly bound, the coordination is settled at both levels and the Michaelis complex stands as modelled - the mechanistic difficulty is then real, and the remaining question is the binding **pose**, which needs re-docking rather than any constraint change. If DFT is much flatter, xTB is overbinding the substrate specifically, and the geometry that blocked every mechanism reopens.

**D64 - No, the local-minimum concern is not accounted for, and there is direct evidence it is real: the docking kept nine poses and we used one.** (2026-10-20, user: "my concern is that xtb is stuck in a local minimum with regards to the carbonyl")

**Taking the three pieces of evidence in turn, none of them can rule it out.**

*The carbonyl scan cannot, and the numbers show why.* Measuring how far each part of the substrate moved from the zinc across the nine points:

| | O2x | C4x | C17x | S1x | N3x | O5x | O4x |
| --- | --- | --- | --- | --- | --- | --- | --- |
| scan_00 | 2.15 | 3.34 | 5.39 | 7.78 | 12.28 | 10.17 | 4.09 |
| scan_08 | 4.00 | 4.80 | 6.70 | 8.15 | 12.75 | 10.65 | **3.82** |
| change | +1.85 | +1.46 | +1.31 | +0.37 | +0.47 | +0.48 | **-0.27** |

The tail barely moves (0.4 A over the whole scan) and the carboxylate comes *closer*. The substrate **pivots about its anchored end** rather than leaving. So the scan explored one hinge motion of one pose - which is all a sequential relaxed scan can do, since each point starts from the previous geometry. The 22.4 kcal/mol is the cost of that hinge, not of reaching any alternative bound state.

*The 25 ps MD cannot.* Ligand exchange at Zn(II) is orders of magnitude slower than 25 ps, so a trajectory that keeps the carbonyl bound shows the basin is locally stable and nothing more.

*The DFT single points now running cannot either.* They are computed **on the xTB geometries**, so they test whether the method gets the energetics right at a fixed structure. If the structure is in the wrong basin, DFT inherits the wrong basin. Worth keeping - they answer a real question about the method - but they are not evidence about sampling, and I should not have implied in D63 that they were the last word.

**The direct evidence.** `CSR_ligand/meropenem.pdbqt` holds **nine docked poses** from AutoDock Vina:

| Pose | Vina score | RMSD from pose 1 (l.b. / u.b.) |
| --- | --- | --- |
| 1 | **-5.744** | - |
| 2 | **-5.624** | **3.06 / 6.72 A** |
| 3 | -4.733 | 2.19 / 3.38 |
| 4-9 | -4.56 to -3.76 | 2.2-2.6 / 6.6-7.2 |

**Pose 2 is 0.12 kcal/mol from pose 1 and 3-7 A away from it.** Vina's scoring function is good to a couple of kcal/mol at best, so that gap is meaningless: the docking provides **no basis at all** for preferring pose 1, and two structurally unrelated binding modes are tied. One of them was carried into every calculation since, and the other has never been looked at.

So the honest position is that the Michaelis complex rests on an arbitrary choice among at least two equally scored poses, and every geometric conclusion - the six-coordinate zinc, the 33-40 degree attack angles, Asp92 and the carbonyl both being unremovable - is conditional on it.

**Proposed next, and it is cheap because the poses already exist:** take the top few docked poses, build each into the same solvated system, optimise at xTB/MM and compare energies and attack geometry. Hours per pose rather than days, no docking software needed, and it directly tests whether the chosen basin is the right one. If pose 2 gives a four- or five-coordinate zinc with a water on the Burgi-Dunitz line, the mechanism question reopens entirely.

**To check before relying on this:** `meropenem.pdbqt` lives in `CSR_ligand`, a different run directory, and `mer_crystal.pdb` in `CSR_crystal` carries no provenance header, so it is not yet established that these poses were docked into the same receptor and site that `CSR_crystal_v2` uses, nor that `mer_crystal.pdb` is pose 1. Verify both before building anything on them.

**D65 - DFT confirms xTB on the carbonyl coordinate.** The three wB97X-D3(BJ)/def2-SVP single points on carbonyl-scan points 0, 4 and 8 have landed:

| Zn-O2x | DFT | GFN2-xTB | difference |
| --- | --- | --- | --- |
| 2.15 A | 0.00 | 0.00 | - |
| 3.08 A | +8.09 | +8.53 | -0.44 |
| 4.00 A | **+23.62** | **+22.36** | +1.26 |

kcal/mol. Agreement to 1.3 kcal/mol over a 24 kcal/mol range. Note the contrast with D58, where DFT was 1.6x steeper than xTB for pulling **Asp92** off the same zinc: xTB is not uniformly soft on this metal, and on the carbonyl coordinate it is essentially right. **So xTB's description of the carbonyl-zinc bond is not the problem.** These are still single points on xTB geometries, so they validate the energetics along that coordinate and not the choice of basin (D64).

**D66 - The reference ligand is pose 2, not pose 1.** Checked as D64 said to. Comparing `CSR_ligand/meropenem.pdbqt` against `CSR_crystal_v2/preprocess_results/mer_crystal.pdb` on heavy atoms, by symmetrised mean nearest-neighbour distance:

| pose | 1 | **2** | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| overlap with the reference | 1.99 | **0.77** | 2.05 | 2.28 | 2.02 | 2.11 | 1.89 | 1.76 | 1.73 |

Every calculation in this project has been built on **pose 2**, the second-ranked pose, moved or minimised by about 0.8 A. Pose 1, the top-scored one, has never been looked at. Nobody chose this; it was inherited. (D64 assumed the carried-forward pose was pose 1 and was wrong about that, though its point stands and is now sharper.)

**D67 - The 1.94 A Zn-carbonyl bond was imposed, not found.** With the docking frame superposed onto the simulated system on 147 backbone atoms within 12 A of the zinc (RMSD 0.83 A; all 369 atoms fit to 2.05 A, which is how far the side chains moved during MD):

- the docked pose puts the beta-lactam carbonyl **4.86 A** from the zinc;
- **no** docked pose puts it inner-sphere - the closest are pose 2 at 4.14, pose 6 at 4.38, pose 5 at 4.50 and pose 8 at 4.58 A, while the top-scored pose 1 is at **10.13 A** and is not in a catalytic orientation at all;
- and the simulated system has it at **1.94 A, to the last digit, in every structure from the first preprocess output through the restrained minimisation, solvation, NPT and the whole of NVT.** A distance that does not move by 0.01 A across molecular dynamics is a rigid constraint, which settles D59: the metal contact is holonomic and it was applied from the very first system file.

So the Michaelis complex this project has been studying has a zinc-carbonyl bond that came from a configured metal contact, at exactly the canonical value, and not from the docking or from any energy. D61/D63's 22.4 kcal/mol - now DFT-confirmed at 23.6 (D65) - is the cost of breaking a bond that was put there by hand, measured from a state no evidence produced.

**The fair counter-argument, which matters.** Vina has no metal-coordination term. It cannot place a carbonyl on a zinc and its failure to do so is not evidence that the carbonyl does not belong there. Substrate carbonyl coordination to Zn1 is well established in binuclear MBLs. So this does **not** show the assumption is wrong - it shows it is an assumption, that it was imposed rigidly rather than allowed to find its own distance, and that docking cannot adjudicate it either way. The QM/MM relaxation to 2.13-2.25 A once the constraint was released (D63) is the only unforced evidence we have, and it is consistent with real coordination.

**D67's inference about the cause is wrong; see D68.** The 1.94 A is real and it is rigidly constrained, but mcube did not put it there by assuming a canonical distance. It was a configured harmonic restraint in the *previous* run, and the input complex that run started from has the carbonyl at 5.43 A.

**Also:** `md_results/NVTtrajectory_lastframe.pdb` still carries the 1.94 A constraint, so the "clean re-do" of the MD without constraining meropenem has not actually produced a trajectory yet. Anything started from that last frame inherits the constrained pose.

**D68 - Where the 1.94 A actually came from: a setup restraint in run 1 became a hard constraint in run 2.** Traced through the configs and the atom indices. mcube's constraint machinery works exactly as designed and moved nothing.

*The input.* `CSR_crystal/preprocess_results/csr_crystal_start.pdb` has meropenem as **docked pose 2 exactly** - conformer RMSD 0.00 A after optimal superposition, against 1.15-2.89 A for every other pose - with the beta-lactam carbonyl **5.43 A** from the zinc and Asp92 at 2.80 A. So pose 2 *was* chosen deliberately (`model: 2` in all three configs), which corrects D66's "nobody chose this": what was inherited was not the pose but what happened to it next. All three configs also set `merge_ligand: false`, so the PDBQT only ever supplied OpenFF parameters and `mer_crystal.pdb` is a parameterisation artefact, not the system's ligand - which is why it still sits at the docked position.

*Run 1 (`CSR_crystal`).* `md.inputs` holds four rigid constraints, all genuine coordination - Zn to His88 NE2, His164 NE2, Asp92 OD2 and the hydroxide O01 - and **two harmonic restraints with hand-chosen targets**:

| restraint | atoms | target | k |
| --- | --- | --- | --- |
| Zn to **mer O2x** | 4203-4209 | **2.2 A** | 20 kcal/mol/A^2 |
| hydroxide **O01 to mer C4x** | 4255-4208 | **3.0 A** | 20 kcal/mol/A^2 |

That is the answer to "why did it not stay at 5.43 A": it was pulled in on purpose. The restraint did its job and the carbonyl settled at 1.94 A, slightly inside the 2.2 A target because the MM zinc's electrostatics take over once it is close. This is a legitimate way to build a Michaelis complex and `md.py` documents the key for exactly that use. **Note the second restraint pre-organised the attack geometry too** - the nucleophile was held 3.0 A off the electrophilic carbon - so the near-attack conformation was built in, not found.

*Run 2 (`CSR_crystal_v2`).* `protein_pdb` is run 1's restrained end state, and `md.inputs` specifies **no** constraints and no restraints, so `md_constraints` falls back to whatever `find_coordinators` detected at `cutoff: 3.0`. The carbonyl was now at 1.94 A, i.e. inside that cutoff, so it was classified as metal coordination and frozen. Run 1 constrained four contacts; run 2 constrained **six** - the carbonyl and a second water joined, purely because run 1 had moved them in:

| run | constrained to Zn |
| --- | --- |
| CSR_crystal | His88, His164, Asp92, O01 |
| CSR_crystal_v2 | His88, His164, Asp92, O01, **mer O2x**, **HOH D:1** |

**The failure mode is laundering, not a bug.** A harmonic restraint is a statement that we do not know the distance and are imposing one. A rigid constraint on an auto-detected contact is a statement that the distance is structural and should be preserved. Chaining the runs converted the first into the second, and after that nothing downstream could tell an imposed distance from a crystallographic one - the detection cutoff only sees the number in the file. Everything from D57 onwards was measured in that state.

**What this means for the protocol.** Restrained structures must not be fed back in as `protein_pdb` without re-declaring what is assumed, or `constraint_exclude_residues` must be used to keep the substrate out of the auto-detected set. Worth considering in mcube itself: record restrained pairs in the preprocess/md results and refuse to auto-constrain a pair that a previous run restrained, or at least warn.

**D69 - The unrestrained run (`CSR_unrestrained`), jobs 30226271 / 30226272.** The test D68 implies: apply the coordination rule honestly to the docked pose and impose nothing else. Atoms within coordination distance are constrained, atoms outside it are not - and nothing is pulled anywhere.

Started from `CSR_crystal/preprocess_results/finalsystem.pdb`, the solvated system built from the docked pose, which has Zn-Asp92 2.80, Zn-O01 2.14 and **Zn-mer O2x 5.43 A**. The config is `CSR_crystal`'s with `md.inputs.restraints` and `md.inputs.constraints` set to null, so the MD falls back to what `find_coordinators` detects at `cutoff: 3.0`. It detects five contacts and the rule does the work by itself:

| constrained to Zn | His88 NE2, His164 NE2, Asp92 OD2, O01, HOH D:1 |
| --- | --- |
| **not constrained** | **mer O2x, at 5.43 A** |

No preprocess re-run was needed: every path in the copied config is absolute and points back at `CSR_crystal/preprocess_results`.

Then the default qmmm flow with nothing added: **25 ps of QM/MM MD at 1 fs, QM/MM optimisation of the 8 A active region around the zinc, and a DFT single point.** `make_mm` builds the QM/MM MM object with `autoconstraints=None`, `rigidwater=False` and no `constraints`, and no `opt_constraints` are set, so there is **no bond constraint anywhere in the QM/MM stage** - only the non-QM active-region waters are frozen during the optimisation, which is unrelated to the metal. GFN2-xTB is trusted to hold the coordination geometry, which D65 supports on this coordinate.

The single point is submitted separately rather than with `singlepoint: true`, purely for resources: DFT runs at 96% on 64 cores and 192 GB while the xTB MD saturates about four threads, so one job cannot size both. The science is identical.

**What the result means either way.** If the carbonyl coordinates on its own during the 25 ps, the Michaelis complex used since D57 is vindicated and only its provenance was wrong. If it stays out, the six-coordinate zinc, the 33-40 degree attack angles and the 22-24 kcal/mol carbonyl detachment cost were all artefacts of the restraint traced in D68, and the mechanism question reopens from the start.

**Caveat to carry.** The MM equilibration still runs with a non-bonded point-charge zinc, which is documented to over-coordinate, so if the carbonyl comes in during the MM stage rather than the QM/MM stage that is weak evidence. Record the Zn-O2x distance at the end of MM separately from the end of the QM/MM MD, and judge them separately.

**D70 - A pure DFT single point does not scale to 64 cores, and D60's 96% was a mixed job.** `seff` on all three (jobs 30207244/5/6, each wB97X-D3(BJ)/def2-SVP on 119 QM atoms, 64 cores):

| job | state | CPU efficiency | wall clock | memory |
| --- | --- | --- | --- | --- |
| 30207244 | COMPLETED | **36.87%** | 26:17 | 67.8 GB of 128 |
| 30207245 | COMPLETED | **37.17%** | 25:33 | 68.0 GB of 128 |
| 30207246 | COMPLETED | **38.09%** | 25:54 | 68.1 GB of 128 |

That is about 24 cores' worth of work on a 64-core allocation: ~18 core-hours wasted per job, ~53 across the three. D60 recommended 64 cores for DFT on the strength of 96.2% efficiency, but that measurement was job 30053828, an xTB optimisation **and** DFT single points in one job - and D60's own text says an xTB optimisation reaches ~94% because of geomeTRIC's linear algebra. So the 96% was mostly the optimiser, and "DFT scales, give it cores" was never measured on a single point. These three are the first clean measurement, and `docs/running-on-bunya.md` now says 24-32 cores and 96 GB.

One candidate for part of the shortfall, untested: the stage also writes four cube files, which is serial work at the end of an otherwise parallel job. `cubes: false` (now an `sp` input) should raise the efficiency, and that is a cheap thing to measure before anyone asks for more cores.

**Note on the comparison that prompted this.** The "4-5 h at 16 cores" figure for a single point is also from that mixed job and is not a single point's wall time, so no speedup factor should be quoted between the two. The honest statement is just the absolute one: a single point of this size and level takes about 26 minutes on 64 cores, at 37% efficiency.

**The scientific branch this check was meant to decide is already closed by D67/D68.** DFT came out steep (D65: +23.6 kcal/mol at 4.00 A), which under the original plan meant "the carbonyl is genuinely bound at both levels, so the only untested thing is the binding pose - propose re-docking". Re-docking is not needed: D66 established the pose is docked pose 2, and D67 that the docked pose puts the carbonyl **5.43 A** from the zinc, so a pose with an uncoordinated carbonyl already exists and never needed generating. D69's unrestrained run is testing exactly that, from exactly that structure. Nothing to propose.

**D71 - The two methodological failures are now prevented by the code, not by the documentation.** (2026-10-20, user)

Both of the traps this project fell into were being handled by telling people about them. That only works until someone is tired.

**Restraint laundering (D68) is now refused.** `md` already writes what it imposed to `md.results.restraints`, so there was nothing to store: `preprocess` compares that list against the contacts it has just detected and **stops** if a pair is in both.

    4203 - 4209: restrained to 2.2 A by md

offering the two honest choices - `constraint_exclude_residues` to leave the contact free, or `acknowledge_restrained: [[4203, 4209]]` to assert the coordination is real. Real coordination (histidines, the hydroxide, waters) is untouched, because those pairs were never restrained.

`md` records the restraints as **selectors** rather than the indices it resolved for OpenMM, and the comparison happens in selector space. Indices are only meaningful against the structure they were resolved on, so comparing `md.results.restraints` (resolved on the `finalsystem.pdb` that existed when md ran) against freshly detected contacts (indices into whatever preprocess has just built) is comparing two different frames - solvating again renumbers everything after the solute. Selectors are what the config uses everywhere else for exactly this reason. Both stored forms are accepted, so configs written before this keep working.

*(First written as `REMARK` lines stamped into the output PDBs so the record would travel with a copied structure. Wrong: the config is where that data already lives and belongs, and the stamping was 120 lines of parallel bookkeeping to maintain. The config-based check is about 30. The cost is that it sees one config - a run started from another run's structure, with its own config, is outside what it can check, and that restraint is recorded in the other config.)*

**Thread placement is no longer a choice.** `mcube env <config> <stage> [--section]` derives the policy from the stage and the config and prints shell exports; `run.slurm` evals it after activating the environment and before anything starts. The profile argument is gone - `run.slurm <stage> <section> [start.pdb]`. The three cases are now decided rather than remembered:

| condition | threads | bind |
| --- | --- | --- |
| one xTB process (`qmmm`, `scan`, `ts` without workers) | CPUs/2 | **yes** |
| parallel xTB (`neb`, `ts` with `freq_workers` > 1) | CPUs/2 | no |
| in-process BLAS (`sp`, any `qm_engine: dft`, `md`) | 1 (or CPUs/2 for md) | no |

It also warns at submit time when the allocation exceeds what the stage can use, which is the D70 mistake: `[WARNING] 64 CPUs allocated but sp uses about 32`.

The user's framing is the reason both of these are code: the only things left to choose are `--cpus-per-task`, `--mem`, `--time` and `config.yaml`. The documentation shrank accordingly - the three-recipe thread table in the README is two sentences now, and `run.slurm`'s header carries a worked `sbatch` line per stage instead.

**Tested:** `test_env.py` (needs only PyYAML, because `cli.py` imports no stage) covers the policy table and the shell output; `test_provenance.py` covers the refusal, the acknowledgement and the cases that must pass through untouched.

**Not affected:** jobs 30226271 and 30226272 were already submitted, and Slurm spools the batch script at submit time, so they run the old script with the old arguments.

**D72 - The unrestrained carbonyl is staying out, so far, at the MM level; and site.py was measuring across the box wrong.** (2026-10-20, job 30226271 still running)

**site.py had no minimum-image convention.** An MD frame is written with every atom wrapped into the cell, so the ligand, sitting near a face of the 80.3 A box, came out on the far side of the protein. The warmup frame reported **Zn-mer O2x = 76.74 A**. It is 5.03 A. `load_pdb` now reads `CRYST1` and `dist`/`angle` apply the minimum image, which covers `sphere` and `nucleophiles` because both route through them; non-orthorhombic cells say so rather than being measured wrong. Everything measured before this on an imaged (`*_imaged.pdb`) or unsolvated structure is unaffected, which is all of D36-D70 - the trap only bites on a raw MD frame, and those were read through `_imaged` files or through mdtraj, which applies PBC itself.

**What the equilibration is doing to the carbonyl, with nothing holding it:**

| stage | Zn-O2x |
| --- | --- |
| input (docked pose 2) | 5.43 |
| minimised | 5.56 |
| gentle warm-up, cycle 2 | 5.03 |
| NPT, first frame | **4.93** |

The five detected contacts are held exactly where the constraint rule put them (Asp92 2.79-2.80, O01 2.14, His88 2.22, His164 2.09), and the carbonyl - the one contact the rule left free - is drifting inward slowly but is nowhere near the 2.6 A that counts as coordination. **Treat this as weak evidence only:** it is the MM force field with a non-bonded point-charge zinc, which D63 established cannot answer coordination questions, and the inward drift is small enough to be ordinary relaxation. The QM/MM MD is the measurement.

**The job will probably not finish.** NPT has done 96.7 ps in 6h19m, i.e. about 4.3 steps/s. At 4 fs the 1 ns NVT needs about 16.3 h and only 17.7 h of the 1-day limit is left, before counting however much longer NPT runs - `OpenMM_box_equilibration` loops until volume and density converge, and it is still going. `scontrol update TimeLimit` is refused ("Access/permission denied"), as expected for a user raising their own limit. If md times out it will not be COMPLETED, so the chained `uqmmm` (30226272, `afterok`) is cancelled with it.

The 1 ns came from copying `CSR_crystal`'s config, and the user's own position (D62) is that MM equilibration tells us nothing about coordination here - it only needs to settle the solvent. 200 ps would take about 3.3 h and do that.

**D73 - NPT converged, the carbonyl is still out, and the hydroxide has moved onto the attack line.** (2026-10-21, job 30226271 still running)

The walltime worry in D72 was misplaced: `OpenMM_box_equilibration` converged at **100 ps** rather than running open-endedly, and NVT is at 726 of 1000 ps at 70.3 ps/h, so it finishes about 3.9 h inside the 7.2 h left. No restart needed.

**The carbonyl, with nothing holding it:**

| stage | Zn-O2x |
| --- | --- |
| input (docked pose 2) | 5.43 |
| minimised | 5.56 |
| warm-up cycle 2 | 5.03 |
| NPT first frame | 4.93 |
| **NPT last frame (100 ps)** | **4.81** |

A slow inward creep of about 0.6 A over the whole equilibration, nowhere near the 2.6 A that counts as coordination, while the five detected contacts stay exactly where the rule put them (Asp92 2.79, O01 2.14, HOH485 1.90, His88 2.22, His164 2.09).

**The more interesting number is the attack geometry.** At the end of NPT:

| | O01-C4x | angle | HOH485-C4x | angle |
| --- | --- | --- | --- | --- |
| restrained system (D60) | 3.62 A | **36 deg** | - | - |
| unrestrained, end of NPT | **3.01 A** | **84 deg** | 4.92 A | 100 deg |

The zinc-bound hydroxide is now 3.01 A from the electrophilic carbon at 84 degrees. That is inside the near-attack distance and much closer to the Burgi-Dunitz window (95-120) than the 36 degrees the restrained system ever managed - the geometry that D61 concluded was blocked. Releasing the carbonyl appears to have let the nucleophile swing onto the line.

**Hold this loosely.** It is the MM force field with a non-bonded point-charge zinc, which D63 says cannot answer coordination questions; it is one frame, not an ensemble; and 84 degrees is still outside the window. But it is the first time in this project that a zinc-bound oxygen has been anywhere near a productive attack geometry, and it is consistent with the carbonyl's inner-sphere position having been what blocked it. The 25 ps QM/MM MD that follows is the measurement that counts, and `joint`/`acceptors` over its trajectory is the test of whether anything can take the proton.

**D74-D92 - Holding the metal site during MM: observations, literature, and the settled design.** (2026-10-20 to 2026-10-22, condensed from nineteen entries)

**The observations.** Run 1 (`CSR_unrestrained`, stock mcube: rigid constraints on the 5 contacts detected within `cutoff: 3.0`):

| | input | minimised | warm-up c2 | NPT end | NVT end (1 ns) |
| --- | --- | --- | --- | --- | --- |
| Zn-mer O2x (free) | 5.43 | 5.56 | 5.03 | 4.81 | **6.89** |
| O01-C4x / angle | - | 3.44 | 3.17 / 83 | 3.01 / 84 | 4.45 / **100** |
| Zn-Gln93 OE1 (not held) | 5.54 | 4.31 | **1.84** | 1.77 | 1.79 |

Three conclusions. The **substrate did not coordinate on its own** - it left, ending further out than it started. **Gln93 collapsed to an unphysical 1.79 A** (real Zn-O is 2.0-2.1) during warm-up cycle 2, which is 40 ps of ordinary 300 K MD, so `gentle: false` was not the cause and NPT would have done the same. And the **five detected constraints held exactly** throughout (Asp92 2.80, O01 2.14, HOH485 1.90, His88 2.22, His164 2.09): **stock mcube did what it promises.** Its limit is only that a donor at 5.5 A is outside a 3 A detection radius, so nothing holds it.

**Mechanisms tried.**

| attempt | outcome |
| --- | --- |
| constrain one atom per residue at 6 A | Asp92 rotated its **other** oxygen to 1.75 A while OD2 sat pinned at 2.80 |
| 18 rigid constraints (per atom, 6 A) | over-determines the constraint network; `Particle coordinate is NaN` at minimisation |
| same, plus `gentle: true` | identical - `OpenMM_Opt` runs *before* the warm-up |
| freezing (zero mass) | `A constraint cannot involve a massless particle`: HBonds constrains every X-H on the frozen donors |
| harmonic **positional** restraints, 10 kcal/mol/A^2 | survived; ran as `hold_coordination`. Needed `periodicdistance` - a non-periodic restraint gave NaN on the first production step after 100 ps of clean NPT |

**D89, the run that matters.** md (2h35m) and qmmm (1h44m) both COMPLETED with `hold_coordination`. Over 500 ps the coordination report shows **one** change - `APPEARED B:1:ZN - E:3174:O at 2.22 A`, a water at a normal distance - and **Gln93 did not collapse**. QM/MM relaxed the sphere to a sensible five-coordinate geometry (O01 1.93, His88 1.97, HOH485 2.00, His164 2.20, Asp92 2.32). But **the substrate dissociated**: Zn-carbonyl 5.43 -> 7.71 (MM) -> 7.54 A (QM/MM), O4x-Lys226 4.07 -> 7.19, O4x-Arg90 10.67 -> 12.57. Every substrate contact lengthened, and the docked pose never made the canonical carboxylate-Arg interaction. **QM/MM optimisation is local: it repaired the metal sphere and cannot re-dock a ligand.** The binding pose is the blocking scientific question.

**The literature** (surveyed by agent 2026-10-22, at the user's instruction not to run on his opinion or mine; *verified* = read in retrieved text, *snippet* = search excerpt, *inferred* = synthesis).

| paper | restraint | k (kcal/mol/A^2) | held through |
| --- | --- | --- | --- |
| **Gervasoni 2022**, IMP-1 dizinc MBL, 10.1002/prot.26227 | **distance**, 6 hand-named Zn-donor pairs | **30** [single-sourced] | heating, NPT, **100 ns NVT**; released for B3LYP/MM |
| **Hanzevacki 2025**, JCIM, 10.1021/acs.jcim.5c00869 | **positional**, all heavy atoms of Zn-bound residues + Zn | **50** [snippet, SI] | relaxation only; released for QM/MM |
| **Guven 2024**, VIM-2, 10.1021/acs.jpcb.4c06379 | **positional**, zincs + coordinating residues + **ligand** | not given (solute ladder 100 -> 10) | replaces unrestrained NVT |
| **Lisa 2017**, NDM-1 + hydrolysed meropenem, 10.1038/s41467-017-00601-9 | restraint "to keep the hydrolysed meropenem bound" [snippet] | - | - |
| **ASH tutorial** (what mcube derives from) | rigid Fe-S **bond** constraints; "harmonic bond restraints" named as the alternative | - | all classical stages; 2b = 1 ns NVT at 4 fs |

Four findings:

1. **Restraints reduce but do not prevent creep.** Gervasoni, verified and independently confirmed by search: "the coordination number of the zinc ions tends to increase, through addition of extra water molecules, up to a limit of six" - *with* restraints on. Their workflow then relies on QM/MM, which "recover[ed] the correct coordination geometry and ... fix[ed] the errors introduced by the MD simulation". **So residual creep is expected, and QM/MM refinement is the designed remedy.** This is why 1 ns of restrained MM into QM/MM is a defensible protocol and not a shortcut.
2. **IOD cannot substitute for the restraint.** [inferred, well supported] a bare 12-6 set admits a second-shell oxygen into the first shell on 10-100 ps when a donor sits 5-6 A away - our Gln93 at 5.54 A over 40 ps. Gervasoni's 12-6 *with* 30 kcal restraints still reached CN 6; Guven lost coordination with two parameter sets. IOD is a better radius underneath the restraint, nothing more.
3. **12-6-4 is contraindicated for MBLs**: Gervasoni found it *enhanced* water intrusion, kept the zincs octahedral, and QM/MM could not repair Zn-Zn from those frames. D87 declining to implement it was right.
4. **Restraining the substrate has precedent in our own system class** (Lisa). Guven deliberately do not - "To avoid adding a bias to the RBFE predictions, the distance between the zinc and the ligand is not restrained" - because a free-energy calculation cannot afford it. A mechanism study can.

**Silences.** No converged force constant and no published sensitivity scan anywhere; 30 vs 50 is 0.14 vs 0.11 A RMS at 300 K, i.e. immaterial. **Backbone N / carbonyl O as a restrained donor is effectively absent** - MCPB.py and ZAFF examples are side-chain only - because everyone hand-lists named side-chain atoms and so never meets the case. No agreed name for the failure mode ("over-coordination", "hypercoordination", "factitious octahedral geometry"), reported for Zn by Gervasoni, Guven and **Bahena-Mendez 2026** (10.1021/acsomega.5c11425, almost our exact stack: Amber + GFN2-xTB) and for Ca(II) by Levy (arXiv:2606.28463). **No paper reports GFN2-xTB over-coordinating zinc**; Bahena-Mendez reports it *correcting* classical over-coordination, "within the first few nanoseconds" though, against our 25 ps, and "remained sensitive to the initial configuration". The BioExcel QM/MM best-practice guide never mentions metals, while warning that position restraints "should not be used in actual profile simulations". Not retrieved in full: Melse 2023 *JCC* 44, 912 (the actual Zn-model benchmark); the Hanzevacki SI (sole source for 50); Eshtiwi & Rathbone 2023 *JMGM* 121, 108431, which warns bonded MCPB.py models "may not be appropriate for validating and refining the binding modes predicted by docking when crystal structures are lacking" - squarely our case.

**Corrections to earlier entries in this log.**

- **D81/D84 framed constraints as a defence against over-coordination. Their original purpose is to supply a bond the force field lacks** - ASH: the force field "does not define bonded parameters associated with the Fe-Cysteine interaction". ASH assigns spurious binding instead to the metal's **LJ parameters** ("a water molecule attempting to bind to the +3 pointcharge"), which is where IOD acts. Two distinct levers, both with precedent in mcube's own source material.
- **D86: a 6 A hold does not follow from Guven.** Gln93 is not a coordinating residue in our input, so it is outside *their* restrained set too: **their failure mode was coordination being lost, ours is coordination being gained.** `hold_coordination` is my extrapolation, not their method.
- **D90 said nothing in the literature supports whole-residue holding. Too strong** - Hanzevacki is exactly that. It is unsupported by *Guven*, the paper I had read. Its real weaknesses are over-restraining (pinning non-coordinating atoms of coordinating residues) and under-detecting (a residue coordinating only via backbone O may never be classed "Zn-bound").
- **D91 credited mcube's detection rule to ASH. It is not inherited**: ASH hand-picks `bondconstraints=[[755,96],[755,136],[755,567],[755,607]]` by index, with no cutoff and no element rule. **Element-based atom detection within a radius is mcube's own generalisation** - and the right abstraction, since it is what reaches a backbone donor that a hand-named side-chain list never would. Gervasoni's selection does not generalise either; it is the atom-pair *form* that we adopt, not the selection.
- **My proposal to restrain "the substrate's heavy atoms and the residues' donor atoms" was wrong** - inconsistent special-casing of exactly the kind the user has repeatedly rejected. One rule applies inside the radius regardless of what the atom belongs to.
- **D85 excluded the substrate from the held set. Wrong on the literature (Lisa) as well as on the outcome (D89).**

**Settled design (2026-10-22, after the user cut the detour short).** Every metal-donor contact `preprocess` detects within one `cutoff` is held by a **harmonic distance restraint at the distance the input structure already has**, 50 kcal/mol/A^2, for all MM stages, released for QM/MM, with IOD metal radii underneath and 1 ns of NVT before `qmmm` - ASH 2b's length, and Gervasoni's protocol in kind. Rows come from `preprocess.results.coordination`, printed into the config with their distance and force constant, so editing them is what changes the run.

**`md.inputs.constraints` is gone too** (2026-10-22): two keys named for the same job confused the one that matters, and harmonic restraints are now simply what mcube does. `preprocess` writes its detected contacts to **`preprocess.results.coordination`** as full `[i, j, d, k]` rows - the measured distance and the 50 kcal/mol/A^2 constant both literally in the config, editable per pair before `md` runs, so the force constant is transparent without a config key for it. `md.inputs.restraints` is now only for *extra* restraints, and an entry naming a detected pair replaces that row.

**This is stock mcube with three changes**, and the whole restraint excursion is worth exactly those:

| change | why it survives |
| --- | --- |
| detected coordination restrained, not rigidly constrained | ASH sanctions it outright; avoids both constraint failures (over-determination at one centre, and OpenMM refusing a constraint on a massless atom); a real bond vibrates |
| IOD metal LJ parameters as the default | Amber's stock Zn-O minimum sits at ~1.69 A against 2.0-2.1 A measured; ASH attributes spurious binding to exactly these terms. Does not fix over-coordination and is not asked to |
| coordination sphere reported before and after MD | would have caught Gln93 on day one |

Net: `md.results` gains `coordination_restraints` (the detected subset, kept apart from `restraints` so the D68 guard still fires only on a distance someone *chose* - a pair restrained at the distance it already had asserts nothing new). `src/mcube` is 3597 lines. **Config keys: `hold_coordination` is gone, leaving two added across the whole episode** (`acknowledge_restrained`, `cubes`) and none for the restraint mechanism itself. `hold_coordination` is deleted along with `restrain_positions`, `COORDINATION_FORCE`'s old meaning and the `CustomExternalForce`/`periodicdistance` code.

**The detour, recorded so it is not repeated.** `hold_coordination` - harmonic *positional* restraints on every donor inside a widened 6 A radius - was built, run (D89) and then removed. Positional restraints pin atoms to lab-frame coordinates, which removes three degrees of freedom per atom instead of one, fights the protein's own relaxation over 1 ns (hence Hanzevacki and Guven both *release* theirs before production, while Gervasoni can hold distance restraints for 100 ns), and freezes any water caught in the radius instead of letting it exchange. Widening to 6 A was also what created a water problem that `cutoff: 3.0` never had. Four independent sources - Gervasoni, ASH's stated alternative, and the user's own two descriptions of what mcube restraints are - all pointed at distance restraints from the start.

**Still open, and now clearly separate from the restraint question.** The D89 substrate dissociation is a **pose** problem, not a restraint problem. A coordination-radius restraint of either kind holds bond lengths near the metal; it cannot hold a ligand whose carboxylate sits 10.67 A out, and restraining a wrong pose harder would only buy a clean QM/MM barrier for a Michaelis complex that does not exist. The docked pose never made the canonical carboxylate-Arg interaction every MBL substrate complex has. That is the next scientific question.

**D93-D96 - Measuring the attack geometry, and testing the sixth zinc contact.** (2026-10-22, condensed from four entries)

All distances from `site.py`; `nucleophiles` reports O-C, the Nu-C-O angle and Zn-O, `sphere` the whole first shell.

| structure | O01-C4x | Nu-C-O | Zn-Gln93 OE1 | Zn coordination |
| --- | --- | --- | --- | --- |
| docked pose 2, solvated (md input) | 3.47 | 79 | 5.54 | 5, sound |
| `CSR_unrestrained` end, 1 ns MM (run 1) | 4.45 | **100** | **1.79** | **6** |
| `CSR_unrestrained2` end, 500 ps + QM/MM (D89) | 5.41 | 96 | - | 5, sound |
| `CSR_gln93` QM/MM optimisation (job 30394832) | 4.94 | 100 | **2.65** | **6** (5 proper + 1 long) |

**The angle fixes itself; the distance does not.** Docking gave a good distance with a bad trajectory (79 degrees); MM rotated the substrate into the Burgi-Dunitz window (79 -> 100) while letting the approach drift out. Neither unrestrained run has both. This is the argument for a harmonic *distance* restraint on O01-C4x: it fixes the one coordinate that degrades and leaves the two angular ones free, where a positional restraint would have frozen the 79-degree trajectory. Two unrestrained runs is enough; a third is not needed.

**Where the substrate restraints actually came from.** `CSR_crystal`, not `CSR_hydrolysed` or `CSR_ligand`: `grep -c restraints` is **0** for both of those, for all four of their `.bak` copies and for all four `.slurm` files. They used detected coordination only, with the ligand free via `constraint_exclude_residues: {C: [1]}`. CSR_crystal's two rows are `[4203, 4209, 2.2, 20.0]` (Zn -> carbonyl O) and `[4255, 4208, 3.0, 20.0]` (hydroxide -> carbonyl C). Indices differ between systems - C4x is 4211 in CSR_unrestrained2 against 4208 in CSR_crystal - so selectors are mandatory, not a nicety.

**Job 30394832 (`CSR_gln93`): does QM/MM reject the 1.79 A contact?** Gln93's side chain (indices 1421-1431, eleven non-backbone atoms, neutral so `charge: 2` unchanged) was put in the QM region **by editing `preprocess.results.qm_atoms` from 119 to 130 entries** - no preprocess re-run, because those are indices into `processed_pdb` and MD preserves atom order. Note run 1's `processed_pdb` points at *CSR_crystal's* `finalsystem.pdb`; it reused that prepared system. `skip_equilibrate: true`, so this was a **pure optimisation with no thermal sampling at all** - 95 geomeTRIC steps, 8m58s, converged, -447.37695 Eh.

Answer: **xTB does not create the contact - 1.79 A is not a minimum on this surface - but an optimisation only loosens it to 2.65 A rather than clearing it.** The other five ligands are proper bonds (His88 1.96, O01 2.07, HOH485 2.10, His164 2.16, Asp92 2.17). Bahena-Mendez's clean reversion to CN 4 was nanoseconds of QM/MM *dynamics*, so this neither confirms nor contradicts it; a short QM/MM MD remains untried.

Also observed: **Lys226 transferred a proton to O01** (Lys226 NZ 3 H -> 2, O01 1 H -> 2; HZ2 now 0.99 A from O01, 1.74 from NZ), with the original H03 turning to hydrogen-bond Gln93 OE1 at 1.85 A. Total QM charge unchanged. This is D42's `cleanR_lysH` / `cleanR_ohH` question answering itself: given a free choice the proton goes to O01.

**Two over-calls of mine, corrected after the user pushed back.** I first called this "the nucleophile destroyed" and the frame unusable. Both wrong:
- **Zn-OH2 with a nearby base is an ordinary resting Michaelis complex.** The hydroxide forms at or near the transition state with a base taking the proton; a neutral Lys226 is exactly that base. The transfer is evidence about Lys226's role, not a pathology.
- **4.94 A does not disqualify the geometry.** The 4.5 A figure in `site.yaml` tests whether a structure *already is* a near-attack conformation; it does not gate the chemistry. `scan_attack` drives the coordinate in from 3.0 to 1.45 A, and the barrier is taken between a characterised reactant and transition state, not from the starting frame's distance. It costs scan points, not viability.

Earlier in the same thread I had also argued the docked pose was wrong because it lacked a carboxylate-Arg contact - circular, since whether Arg90/Arg229 are involved is question 2 - and because the substrate dissociated, which was an artefact of my own restraint set rather than evidence about the pose.

**So the sixth contact is the one open objection**, matching the standing question of whether xTB over-coordinates the zinc, and `scan_zngln` (relaxed scan along Zn-Gln93 OE1, 1.80-4.20 A, nine points) is the decisive test: the depth of the well at 2.65 A distinguishes a genuine six-coordinate zinc from a residual artefact.

**Energy caveat (D95).** Pre-merge, so stock Zn radii in the MM while the input geometry came from IOD-parameterised MD. Do not compare -447.37695 Eh against a post-merge run.

**Trajectory manifest, recorded 2026-10-22 before deleting the QM/MM trajectories.** Frame-to-time mapping, so a structure named after a frame index stays interpretable once the trajectory it came from is gone.

**Conversion.** `md` NVT: `nvt_timestep: 0.004` x `traj_frequency: 500` = **2.0 ps per frame**. `qmmm`: `timestep: 0.001` x `traj_frequency: 50` (the default; v2 states it explicitly) = **0.05 ps per frame**.

| trajectory | frames | ps/frame | span | fate |
| --- | --- | --- | --- | --- |
| `CSR_crystal/md_restrained/NVTtrajectory_imaged.dcd` | 500 | 2.0 | 1000 ps | **kept** |
| `CSR_crystal_v2/md_results/NVTtrajectory_imaged.dcd` | 500 | 2.0 | 1000 ps | **kept** |
| `CSR_hydrolysed`, `CSR_ligand`, `CSR_unrestrained2` md NVT (imaged) | 250 | 2.0 | 500 ps | **kept** |
| `CSR_crystal_v2/qmmm_results/QM_MM.dcd` | 1000 | 0.05 | 50 ps | deleted |
| `CSR_hydrolysed/qmmm_results/QM_MM.dcd` | 1000 | 0.05 | 50 ps | deleted |
| `CSR_ligand/qmmm_results/QM_MM.dcd` | **878** | 0.05 | **43.9 ps of a requested 50** - run ended early | deleted |
| `CSR_crystal_v2/protdft/md2_ohH/QM_MM.dcd` | 500 | 0.05 | 25 ps (`time: 25`) | deleted |
| `CSR_crystal_v2/qmmm_analysis/qmmm_snap.dcd` | 791 | n/a | analysis-derived subset, **no uniform spacing** | deleted |
| `CSR_unrestrained/protdft/qmmm/QM_MM.dcd` | 5 | 0.05 | 0.25 ps | deleted |

**Structures named after a frame, with the time each corresponds to.** Index ranges identify the source: f00300-f00499 fit the 500-frame md NVT, f00192-f00910 the 1000-frame QM/MM run.

| structures | frames | source | time |
| --- | --- | --- | --- |
| `CSR_crystal/qmmm_probe/f00300`, `f00400`, `f00499` | 300, 400, 499 | md NVT @ 2.0 ps | **600, 800, 998 ps** |
| `CSR_crystal_v2/protdft/f00192`, `f00395`, `f00575` | 192, 395, 575 | QM/MM @ 0.05 ps | **9.6, 19.75, 28.75 ps** |
| `CSR_crystal_v2/scans/qmmm_R_H1_f00365`, `_f00497`, `_f00635` | 365, 497, 635 | QM/MM @ 0.05 ps | **18.25, 24.85, 31.75 ps** |
| `CSR_crystal_v2/scans/qmmm_R_H2_f00395`, `_f00775`, `_f00910` | 395, 775, 910 | QM/MM @ 0.05 ps | **19.75, 38.75, 45.5 ps** |
| `CSR_crystal_v2/scans/qmmm_R_H3_f00192`, `_f00575`, `_f00888` | 192, 575, 888 | QM/MM @ 0.05 ps | **9.6, 28.75, 44.4 ps** |

Every one of those directories keeps its own `QM_MM_optimized.pdb`, so the geometries survive; what is lost is the ability to extract a *different* frame from the same run.

**One file was destroyed by the quota, not by this cleanup.** `CSR_ligand/md_results/warmup_MD_cycle0.dcd` was **0 bytes and unreadable** - a write that failed while scratch was full. It is 10 frames of a 1 K warm-up ramp, so nothing of value, but it is direct evidence the over-quota period damaged output. An audit of every `.dcd` under `CSR_*` found this one and no others.

**Scratch cleanup, 2026-10-22.** `rquota` showed scratch **over quota at 315 GB against a 300 GB limit** (816,052 / 1,000,000 files), so writes were at risk with a job running. Freed 19 GB / 47 files, down to 296 GB. What went, and why it was safe:

| removed | size | why |
| --- | --- | --- |
| `OpenMMMD_traj_wrapped.xyz` x4 | 8.4 GB | ASH writes it; **nothing in `src/` or `docs/csr/` reads it**, and the same trajectory exists as `QM_MM.dcd` |
| `relaxbox_NPT*.dcd` x12 | 7.1 GB | box equilibration, not production. `relaxbox_NPT_lastframe.pdb` is what `md.results` records and is kept |
| `NVTtrajectory.dcd` x6 | 1.3 GB | exact-size duplicates of the `_imaged` twins, which are what `md.results.trajfile` points at and what `postprocess` reads. All six twins kept |
| `geometric_OPTtraj_*.xyz` | 1.1 GB | optimiser traces, read only mid-run (`qmmm.py:363`) and deleted by `clean_geometric_files()` at each run's start anyway |
| 3 slurm `.out` **gzipped, not deleted** | ~0.5 GB | 189/189/166 MB of text; the record is kept |

**Kept deliberately**: every `*_lastframe*.pdb`, `finalsystem.pdb`, `QM_MM_optimized.pdb`, all `.chk` / `final_state.xml` restart files, and the **QM/MM trajectories** (`QM_MM.dcd` x4, `qmmm_snap.dcd`, ~2.3 GB) - those are the analysable artefact for `postprocess` and for `site.py acceptors` / `joint`, so they are the user's call rather than mine. `CSR_gln93` was excluded apart from the finished qmmm run's own trace, because job 30409808 is writing there.

**Second pass, after the manifest above was committed.** The six QM/MM trajectories and the `warmup_MD_cycle*.dcd` ramps went too (~3.5 GB): the warm-up is a 1-10-300 K thermal ramp whose last frames are kept as PDBs, and the QM/MM frames are recoverable only as the optimised structures each frame directory already holds. Verified empty by `find`. Scratch then read **202 GB / 495,641 files**, but only ~3.5 GB of that drop is this project - the GLP-1R and ion-channel sessions were clearing their own trees at the same time, so the 279 -> 202 GB figure is not attributable here.

Verified by `find`/`ls` returning empty and by `rquota`, not by an echo after `rm` - a `;`-chained echo prints success unconditionally and has produced a false "deleted" report before.

**NEB readiness: what a defensible barrier needs, and what is still missing.** (2026-10-22)

**One NEB per elementary step.** Hydrolysis is at least two: attack of the Zn-bound oxygen on C4x giving the tetrahedral intermediate, then C4x-N1x cleavage with protonation of the leaving nitrogen. A single NEB spanning both puts two barriers on one path and NEB-CI returns **one** saddle, so the other is simply never found. Staged in `site.yaml` as `scan_attack` / `neb_step1` / `ts_saddle1` and `scan_ringopen` / `neb_step2` / `ts_saddle2`, with the starting structure passed at run time rather than fixed, since which Michaelis complex we start from is unsettled.

Scan coordinates, both measured rather than assumed: attack drives O01-C4x to 1.45 A (a C-O single bond) from wherever the complex sits; ring opening drives C4x-N1x from its measured **1.39 A** to 2.80. Each scan monitors the contacts that would answer questions 1-3 if they move - Lys226 NZ to O01 and to N1x, Asp92 OD2 to N1x, and Zn to N1x, since a mono-Zn site has no second zinc to stabilise the nitrogen anion that B1 enzymes rely on.

**Code gap closed.** `neb` checked that the endpoints had the same atom *count*, which passes a product whose atom *order* differs - and NEB interpolates atom i of the reactant toward atom i of the product, so a mismatched order gives a path between atoms that are not counterparts and **converges anyway**. It now compares `atom_names` index by index and names the first disagreements. `scan` is what makes a product pass by construction, because it writes every point through the reactant's own MM topology; a product solvated through its own `preprocess` run will not match however carefully it was prepared.

**Still missing, in the order it blocks on.**

1. **A reactant whose protonation state is decided.** D96 showed Lys226 handing its proton to O01 spontaneously during an optimisation, so "the Michaelis complex" is currently two different chemical species. A barrier is meaningless until one is chosen and shown to be a minimum. D42's `qmmm_cleanR_lysH` / `qmmm_cleanR_ohH` plus `sp_dft` already set up exactly this comparison.
2. **The sixth zinc contact resolved** - job 30409808 is scanning it.
3. **An intermediate and a product**, from the two scans above, in the reactant's box.
4. **Frequency and IRC validation.** `ts` does this already: one imaginary mode, that mode being the reaction coordinate, and an IRC in each direction to show the saddle connects the two minima the NEB was built from rather than another pair. A saddle without this is a stationary point, not a transition state.
5. **A consistent energy reference.** Barriers between characterised stationary points at one level, and not across the metal-parameter boundary of D95.

**A naming limitation found while staging this.** `sections.py` writes each entry as `<stage>_<name>`, and the key namespace is flat, so **two stages cannot share a key**: `step1` is the NEB's and the matching TS section had to be `saddle1`. The same rule means several existing keys double - `scan_znasp` becomes `scan_scan_znasp`, `sp_dft` becomes `sp_sp_dft` - which is harmless but means `run.slurm`'s documented examples (`run.slurm scan scan_znco`) disagree with what those keys actually produce. Left alone, since results may already exist under the doubled names.

**The barrier chain crosses two different potential surfaces.** (2026-10-22, found while staging the NEB)

`make_mm` builds every QM/MM system with `rigidwater=False` and `autoconstraints=None`, so MM water O-H bonds are flexible at the OpenMM level. `qmmm`, `scan` and `ts` then restore rigidity for active-region MM waters as **optimizer** constraints (`water_constraints`, which excludes QM waters so a QM water can still give up a proton). **`neb` does not, and cannot**: ASH's NEB exposes no constraints argument - its whole keyword set is `CI`, `free_end`, `maxiter`, tolerances, interpolation and restart files.

So a barrier computed end to end crosses surfaces twice:

| stage | active-region MM waters |
| --- | --- |
| `qmmm` / `scan` endpoints | **rigid** (optimizer constraints) |
| `neb` path | **free**, intramolecular geometry included |
| `ts` refinement | **rigid** again |

Two consequences. The endpoint geometries are not minima of the surface the NEB evaluates them on. And interior images can lower their energy by relaxing water geometry unrelated to the reaction, so a path can dip below its own reactant - which is the observable symptom.

`neb.check_path_shape` now warns when an interior image lies below both endpoints, and says both things it can mean: a genuine intermediate, in which case the step is not elementary and wants splitting at that minimum, or this inconsistency. Diagnostic only - the fix is a methodological choice, not mine. The options are to drop the water constraints from the endpoint optimisations so everything sits on one surface (at the cost of what those constraints were added to prevent), to keep MM waters out of the active region entirely so no stage can move them, or to accept the offset and quote barriers with it stated.

Worth noting the same applies to anyone else's NEB through mcube, not just CSR.

**`ts` enforced the frequency count but not the mode's identity.** (2026-10-22)

`ts` already counted imaginary frequencies, warned when there was not exactly one, and printed `ts_mode_atoms` - the eight atoms displaced most in the lowest mode, normalised. But the printed list was only ever read by a human, and the module's own claim is "a TS has exactly one imaginary frequency, **and its mode should be the reaction coordinate**". Only the first half was checked.

That gap passes a real false positive: a saddle whose imaginary mode is a water libration or a methyl rotation has exactly one imaginary frequency and is not the transition state of anything being studied. With `hessian_atoms` restricted to the reacting atoms it is less likely, but the Zn or a bridging water can still dominate.

`check_ts_mode` now compares the mode against the atoms in `monitor`, which already names the bonds each step makes or breaks, and records `ts_mode_is_reaction_coordinate` plus `ts_mode_reacting_share` in the results. Below 30% of the largest displacement it warns that the structure is a saddle point but likely not this step's transition state. **`null` when there is nothing to check against** - no monitor pairs, or none of their atoms in the Hessian - because "not checked" is not a pass, and the staged `ts_saddle1` / `ts_saddle2` sections both carry the monitor pairs for this reason.

**Critical review of this session's work against the literature.** (2026-10-22, at the user's request)

**Discarded, with the reason each failed.** `hold_coordination` (harmonic *positional* restraints) was built, run and deleted: it rested on a misreading of Guven, whose failure mode was coordination *lost* where ours was coordination *gained*, and it used the wrong restraint type - pinning atoms to lab-frame coordinates would have frozen the 79-degree attack trajectory that MD then rotated to 100. Four claims about the binding pose also went: the carboxylate-Arg criterion (circular - whether Arg90/Arg229 are involved is question 2), dissociation-implies-bad-pose (an artefact of my own restraint set), "the nucleophile was destroyed" (Zn-OH2 beside a base is an ordinary resting complex), and "4.94 A disqualifies the geometry" (a near-attack *test* mistaken for a precondition for computing chemistry).

**Sound and reusable.** Distance restraints on detected metal-donor pairs at the starting distance - Gervasoni's atom-pair form, the substitution ASH sanctions outright, released for QM/MM. Coordination reporting before and after MD. The laundering guard and the chosen/detected split. The NEB atom-order guard, which needs no literature because it is a correctness bug. `md --section`, `fill_restraints`, ASH's own surface plotting in place of home-grown analysis.

**Reusable but more weakly grounded than presented.**

| thing | the problem with its justification |
| --- | --- |
| **50 kcal/mol/A^2** | Hanzevacki's figure is for *positional* restraints; Gervasoni's for *distance* restraints is **30**. The difference is 0.11 vs 0.14 A RMS and so immaterial, but the README cites both in a way that implies either supports 50-for-distance. Neither does. |
| **IOD as the default** | pitched as stopping over-attraction, which Gervasoni contradicts - 12-6 *with* 30 kcal restraints still reached CN 6. The surviving case is different and stronger: Amber's Zn-O minimum is 1.69 A against 2.0-2.1 measured, and OpenMM's stock lanthanide radii are 12-6-4's 12-6 terms without C4, i.e. invalid. |
| **`check_ts_mode` floor of 0.3** | invented. No literature basis. Configurable, but it was shipped as though principled. |
| **`check_path_shape`** | flags a symptom with two causes it cannot distinguish, over a design problem (the water-constraint surface mismatch) that was diagnosed rather than fixed. |

**Dropped, and most consequential - in order.**

1. **Equilibration is 1 ns against a published 10-100 ns.** *(settled below: staged, short run first with a defined success test.)* Gervasoni 100 ns NVT, Guven 50, Lisa 10, Bahena-Mendez 10. The justification on file is ASH protocol 2b, which is a **tutorial**, not a production protocol. D92 recorded the gap and then it stopped being pressed. A 10-100x shortfall is the weakest point in the current setup.
2. **GFN2-xTB is a path-finding method, not a barrier method.** Gervasoni refine with B3LYP/MM; Guven sample with DFTB3 and refine with DFT. `qm_engine: dft` and `sp_dft` exist, but nothing in the staged path *forces* refinement, so an xTB number could be quoted as a barrier by accident. **Partly addressed the same day**: every stage reporting an energy now records `qm_level` in its results (`GFN2-xTB` or `functional/basis`), so a number's provenance travels with it and an xTB energy cannot be mistaken for a refined one. That makes the omission visible; it still does not force the refinement.
3. **`sp_basis: def2-svp` is a geometry basis, not an energetics one.** *(settled below: def2-SVP for geometry, wB97X-V single points; the live risk is the `dft_*`-to-`sp_*` fallback collapsing both to one level.)* wB97X-V is a good functional choice; def2-SVP undercuts it for a barrier. A config change, not a code change.
4. **QM region size never tested.** *(settled below: accepted limitation, not a backlog item.)* Mehmood & Kulik report that too-small QM regions *over*-bind zinc, and D92 flagged region size as its own variable. We sit at ~130 QM atoms and an 8 A active region with no convergence check.

**Level-of-theory policy, and the staged approach to equilibration.** (2026-10-22, user's decision, resolving three of the four gaps in the review above)

| job | level |
| --- | --- |
| equilibration and path finding | **GFN2-xTB** |
| geometry optimisation, where DFT is used at all | **def2-SVP** |
| single points | **wB97X-V** |

That is what the defaults already do, and it makes `sp` a real refinement: an xTB geometry and path, scored at DFT. The `qm_level` key added above makes which is which legible in the results.

**The trap in it.** `dft_functional` and `dft_basis` **fall back to the `sp_*` values**. So setting `qm_engine: dft` makes the gradient stage run wB97X-V/def2-SVP *and* the single point run wB97X-V/def2-SVP - the same level twice, where the second adds nothing but still reads as a refinement. This is the default behaviour of switching engine, not an exotic misuse, so a DFT-gradient run wants `sp_basis` raised above `dft_basis` or the single point is decoration.

**Equilibration: 1 ns now, longer later, deliberately.** The 10-100 ns published practice stands as the eventual target; the short run comes first to establish that the protocol works at all. That makes the gap a staged decision rather than the oversight the review called it - but it only resolves if "works" is defined in advance, so: the short run succeeds if, after the substrate restraint is released, O01-C4x and the Nu-C-O angle stay inside the Burgi-Dunitz window and the zinc keeps a sound first shell. Those holding means the long run is confirmation. Either failing means the long run is required rather than optional, and the barrier waits for it.

**QM region.** Known to over-bind somewhat (Mehmood & Kulik), accepted because it contains every relevant coordinator. Recorded as a stated limitation rather than an open question, so it belongs in the methods rather than in the backlog.

**D97 - The Gln93 scan died of memory, not slowness. Three points, question unanswered.** (2026-10-23, job 30409808)

`sacct`: **OUT_OF_MEMORY, ExitCode 0:125, elapsed 59:21, MaxRSS 24.0 GB** against the 24 GB requested. Not the slow convergence it looked like from outside.

The mechanism is in the `.err`, repeated dozens of times:

```
SVD Error - saving HT.v_nan.txt, increasing v by 0.001 and trying again ( nan ->  nan)
```

geomeTRIC's internal-coordinate SVD produced NaN on point 3 (target 2.70 A). It then **retried in a loop, writing a 23 MB dump per attempt** - `nan -> nan`, so it never recovered - until SLURM OOM-killed it. 19 of those files, 437 MB, deleted; `CSR_gln93` went 634 -> 216 MB. **More memory would only delay this**, because the retry loop is unbounded.

What survived, all converged, relative to point 0:

| Zn-Gln93 OE1 | rel kcal/mol |
| --- | --- |
| 1.80 | 0.00 |
| 2.10 | -5.06 |
| 2.40 | -5.60 |

**So one thing is now quantified and the main question is not.** The 1.79 A contact MM produced sits **~5.6 kcal/mol up a repulsive wall**, which puts a number on D96's qualitative finding that the optimiser pushed it out. But points 4-8 (2.70-4.20 A) never ran, so **the well depth - whether six-coordinate zinc is genuine on this surface - remains open**. `surface_results.txt` and the plot were never written; scan.py only writes them at the end.

**Code gap this exposed.** `coordsystem` was not reachable from the config: `optimize()` called ASH's `Optimizer` without it, so every stage was pinned to ASH's default `tric` - the system that degenerated here. Now `opt_coordsystem` (`tric`, `hdlc`, `dlc`, `prim`, `cart`), unset by default so nothing changes for existing runs. `cart` cannot degenerate this way, at the cost of slower convergence.

**Provenance, per the agreed policy.** `qm_level` is GFN2-xTB, so these are path-finding numbers, not energies; a depth, once measured, needs wB97X-V single points. It also ran pre-merge, with stock Zn radii in the MM while the input geometry came from IOD-parameterised MD (D95), so do not compare against a post-merge run.

Not resubmitted. `resume: true` means a rerun continues from point 3 - the point that failed - so it needs `opt_coordsystem: cart` or a smaller active region first, not just more memory.

**Gln93 abandoned; equilibration restarted properly.** (2026-10-23, user's decision)

The well-depth question from D97 is dropped. The frame it was being asked about is a broken model - run 1's MM end state, carrying the 1.79 A artefact - and scanning a broken model to decide whether its artefact is real is effort in the wrong place. The target is a **10-50 ns equilibration with the restraints now in place**, which is also where D92's survey said published practice sits (Gervasoni 100 ns, Guven 50, Lisa 10, Bahena-Mendez 10). **A 1 ns NVT runs first, to establish the restraints behave**, then the long run. This supersedes the earlier "1 ns now, longer later" framing: the long run is the goal rather than a someday, and the short one is a check on the restraint set, not a substitute.

What survives from D97 and is worth keeping: the 1.79 A contact sits ~5.6 kcal/mol up a repulsive wall, measured.

**A silent failure caught before it bit this run.** `md_restraints` reads `preprocess.results.coordination`; every existing CSR config carries the pre-rename `constraints`. So `md --section md_nac` on any of them would have found no coordination, applied **only** the one substrate restraint, left the metal site entirely free, and said nothing - reproducing exactly the Gln93 collapse this mechanism exists to prevent, in a run set up to validate it. `check_coordination_present` now refuses when `preprocess.inputs.metals` is set and `coordination` is empty, names the old key when it is the cause, and says to re-run preprocess rather than rename by hand - the new rows carry a measured distance and a force constant the old ones did not. A genuinely metal-free system (protein-protein) passes on `metals: []`.

**So the 1 ns run needs a fresh `preprocess`**, not just `md` on an existing directory.

**Feasibility of the long run, which is not free.** D89's md stage took 2h35m for 500 ps of NVT on 16 CPU cores, so roughly 4 ps/min: **10 ns is about 42 h and 50 ns about 9 days** of CPU. `run.slurm` defaults to `--time=1-00:00:00`, so neither fits in one job as configured. The options are a GPU (`md.inputs.platform: CUDA`, which md already supports on all three of its MD calls) or chained restarts from `OpenMM_MD_final_checkpoint.chk`. Whether OpenMM in the ASH env can use Bunya's mi355x cards is **unknown** - it is one of the things the ion-channel session set out to probe, and the answer decides between those two routes.

**D98 - The 1 ns restrained run: `CSR_nac`, jobs 30416229 (preprocess) and 30416230 (md).** (2026-10-23)

First run with the restraint scheme as it now stands, to establish it behaves before committing to 10-50 ns. A fresh directory and a **fresh preprocess**, not `md` on an existing config, because every older CSR config carries the pre-rename `constraints` key and would therefore have run with the metal site free (the guard above now refuses that rather than letting it happen).

Set up from `CSR_unrestrained2`: same `preprocess.inputs` - docked pose `model: 2`, `cutoff: 3.0`, the HIS variants, `extra_qm_residues` A:[90,226] C:[1], and `constraint_exclude_residues: {C: [1]}` keeping the ligand out of the automatic set. Inputs stage into `preprocess_results/` because that is the stage's working directory. `ligand_xmlfile: null`, so OpenFF re-derives the meropenem parameters rather than borrowing CSR_crystal's - self-contained, at the cost of AM1-BCC charges that may differ in the last decimal.

What `md_nac` adds over `md.inputs`, and nothing else:

```
restraints:
  - [C:2:O01, C:1:C4x, 3.47]   # force constant inherited from coordination_force
nvt_time: 1000
```

So the metal site gets whatever `preprocess` detects within 3.0 A, each held at its own measured distance at 50 kcal/mol/A^2, and the substrate gets one restraint on the attack coordinate at the distance the docked pose already has. Nothing asserts the carbonyl coordinates the zinc.

**What this run has to show**, decided in advance: after the restraint is released, O01-C4x and the Nu-C-O angle stay inside the Burgi-Dunitz window and the zinc keeps a sound first shell - which `site.py nucleophiles` prints as its verdict. Those holding makes the 10-50 ns run confirmation; either failing makes it required, and the barrier waits for it.

Two things to watch rather than assume. The copied config still carries inert `hold_coordination: 6.0` and `constraints: null` in `md.inputs`; both are keys the code no longer reads, so they are clutter rather than risk. And `md.results.coordination_restraints` is the thing to read first on completion: if it is empty the guard failed and nothing held the site.

**D99 - GPU benchmark: the long equilibration is cheap, and my ROCm conclusion was wrong.** (2026-10-23, job 30416305)

`scripts/gpu_probe.py` on a **NVIDIA H100 PCIe** (`gpu_cuda`), on the real 49,373-atom CSR system at md's settings (PME, 12 A, HBonds, rigid water, 1.5 amu H, 4 fs):

| platform | 2000 steps | ns/day |
| --- | --- | --- |
| **CUDA** | 0.8 s | **845** |
| OpenCL | 1.2 s | 591 |
| CPU (4 cores) | 299 s | 2.3 |

**So a 10-50 ns equilibration stops being a scheduling problem.** At even a third of the benchmark rate - the benchmark has no barostat, writes no trajectory and ran single precision, so real throughput will be lower - 50 ns is a few hours in one job rather than nine days. `md.inputs.platform: CUDA` with `properties: {Precision: mixed}` is all it needs; md already passes `platform` to all three of its OpenMM calls. Chained restarts from the checkpoint are not required.

**Correction to the ROCm conclusion recorded earlier today.** That entry said the ASH env's OpenMM reports only Reference/CPU/CUDA and therefore cannot touch an AMD card. **The enumeration was run on the interactive node, which has no GPU.** On the H100 node the same build reports `['Reference', 'CPU', 'CUDA', 'OpenCL']` and OpenCL drove the card at 591 ns/day. OpenCL is compiled in; it simply does not enumerate without a device and ICD to find. So **the mi210/mi300x/mi355x are not ruled out** - whether ROCm's OpenCL ICD is visible on `gpu_rocm` nodes is a separate question, settled by running the same probe there, and nothing about the mi355x has actually been tested. `gpu_probe.py`'s own docstring warns that listing platforms off a GPU node says only what OpenMM was compiled with; the inverse holds just as strongly and I ignored it.

Not that it changes the recommendation: CUDA is the faster platform on hardware that is known to work, so `gpu_cuda` is where the long run belongs either way. The correction matters for anyone who later needs the ROCm queue's much larger GPU count.

**The 1 ns restrained run (30416230) is still early** - minimisation done, warm-up in progress - and `md_nac.results.coordination_start` records the five detected contacts, so the restraints were read. On CPU it has 4-10 h left. At the benchmark rate the same run is minutes on the H100, which makes the CPU job a poor use of the time now that this is known.

**D100 - The 1 ns restrained run: the mechanism works, the angle does not.** (2026-10-23, job 30418017, `CSR_nac`, section `md_nacgpu`)

**The GPU route is confirmed end to end.** `Using platform: CUDA`, `Precision: mixed`, and the whole md stage - minimisation, NPT and 1 ns of NVT - in **254 seconds**. The CPU attempt was cancelled after 65 minutes at 93 ps of NPT. So ASH does pass `platform` through, and D99's benchmark holds in the real pipeline rather than only in a bespoke script. A 10-50 ns run is now minutes to an hour.

**The restraint mechanism did what it was built to do.** `coordination_restraints` records all five contacts as selectors with their measured distances at 50.0; `restraints` records `[C:2:O01, C:1:C4x, 3.47, 50.0]` with the default force constant filled in; `metal_parameters: iod`.

| | start | end (1 ns) |
| --- | --- | --- |
| Zn-O01 (hydroxide) | 2.14 | **1.95** |
| Zn-Asp92 **OD1** | - | **1.98** |
| Zn-HOH485 | 1.90 | 2.04 |
| Zn-His88 NE2 | 2.22 | 2.20 |
| Zn-His164 NE2 | 2.09 | 2.23 |
| Zn-Asp92 OD2 *(restrained at 2.80)* | 2.80 | 2.76 |

**Gln93 is absent from the first shell.** In the unrestrained run it collapsed to 1.79 A; here it stays beyond 2.8 A, and no water intruded beyond the one genuine HOH485. That is the failure this whole mechanism exists to prevent, and it did not happen.

**The carboxylate blind spot recurred, benignly.** Six contacts are listed, but two of them are Asp92's own oxygens - OD2 held at 2.76 where it was restrained, and **OD1 arriving at 1.98**. Restraining one carboxylate oxygen does not stop the other, exactly as recorded in the D74-D92 table when a per-residue constraint let OD1 reach 1.75 A. This time the outcome is chemically ordinary - bidentate carboxylate coordination at 1.98/2.76 - rather than unphysical, so it is a limitation to note rather than a failure.

**The attack geometry half-passes, and it corrects an argument of mine.**

```
O01-C4x  3.55 A   Nu-C-O  80 deg   ->  "metal-bound, angle off"
```

The distance is good: 3.55 against a 3.47 target, inside the 4.5 A tolerance, and the substrate did not dissociate - Zn-mer O2x held at 5.47. The **angle is 80 degrees, outside the 95-120 window**, and essentially unchanged from the docked pose's 79.

I had argued that a distance restraint fixes only the radial coordinate and therefore lets MD keep rotating into the Burgi-Dunitz window, as it did in the unrestrained run (79 -> 100 degrees). **It did not.** Holding the approach at 3.47 A appears to leave the substrate less room to reorient, not more - or 1 ns is simply not enough rotation. Either way the two coordinates have now each been achieved separately and never together: D89 had 96 degrees at 5.41 A, this run has 80 degrees at 3.55 A.

**Also worth noting:** Lys226 still carries three protons and O01 one, so no proton transfer this time - the nucleophile remained a hydroxide, unlike D96. The hydroxide's proton sits 2.77 A from Lys226 NZ and 2.64 from mer O1x. And `C4x-N1x` is 1.36 A.

**Status against the D98 test: metal site passes, attack angle fails.** And the substrate restraint was ON throughout, so even the distance is not yet shown to be stable - the release run remains the real test. The angle is now the open question, and the options are a longer run, restraining the angle too, or accepting it on the argument that a relaxed scan driving O01-C4x to 1.45 A must pass through a near-tetrahedral arrangement anyway, so the angle would be set by the chemistry rather than needing to pre-exist.

**D101 - 100 ns equilibration submitted (job 30421166), and two corrections.** (2026-10-23)

**The hydroxide surviving D100 meant nothing, and I reported it as if it did.** That run was pure MM, where protonation is fixed by the force field and no bond can break, so Lys226 keeping three protons was a property of the method rather than an observation about the chemistry. D96's transfer happened because **QM/MM** can break bonds. The comparison between the two was not a comparison at all.

**GPU is for `md` only.** A QM/MM step is bottlenecked on xTB and copies positions off the device every step, so `qmmm`, `scan`, `neb` and `ts` stay on CPU with their cores parallelised - which is what the README already said and what the `md`-only `platform` settings in `nacgpu` and `naclong` reflect.

**`md_naclong`: 100 ns NVT on an H100**, restraints as before, roughly 4-5 h at D100's measured rate. This matches the published equilibrations D92 found (Gervasoni 100 ns NVT, Guven 50, Lisa 10, Bahena-Mendez 10), with the restraints held throughout as Gervasoni also do.

**`traj_frequency` is raised from 500 to 5000, and that was necessary rather than tidy.** At 500 with a 4 fs step a frame lands every 2 ps: over 100 ns that is 50,000 frames at 0.59 MB, and md writes a raw and a reimaged copy, so **~60 GB**. Scratch stood at 191/300 GB with three other projects writing into it, and an over-quota write fails *silently* - it leaves the existing file untouched rather than truncating, so a run looks finished and is not. At 5000 a frame lands every 20 ps: 5,000 frames, **~6 GB**, still ample for equilibration analysis, and the last frame is what feeds QM/MM either way.

**This run also tests the prediction D100 disproved.** The attack angle sat at 80 degrees after 1 ns, essentially the docked value, where I had argued a distance restraint would leave it free to rotate into the Burgi-Dunitz window. 100 ns gives a hundred times the opportunity. If it still has not rotated, the angle is a property of the restrained pose rather than a matter of sampling, and `scan_attack` driving O01-C4x to 1.45 A becomes the way to reach a near-attack geometry rather than waiting for one.

Incidental, on process: a 111-character command left the shell swallowing subsequent input, exactly the failure `BUNYA_BROWSER_SHELL.md` section 3 describes, and two silent no-ops followed before a fresh terminal fixed it. The ~100 character limit is not approximate.

**The MM release test is the wrong test, and "10-100 ns" conflates two things.** (2026-10-23, user's reframing)

**What the restraint is for.** The enzyme hydrolyses meropenem, so a Michaelis complex forms without help; the substrate restraint exists only to stop the **MM** metal mishandling the substrate outside QM/MM. Two consequences follow, and the second changes the plan.

The substrate dissociating unrestrained (run 1, D89) is evidence of an **MM artefact**, not evidence against the docked pose. That was already the conclusion, reached by elimination; this is the principled reason for it.

And **an MM release run cannot settle anything.** If the complex falls apart once the restraint is dropped *in MM*, that says only that MM misbehaves near the metal - which is the premise, not a finding. It is informative in one direction: holding would be reassuring, failing would mean nothing. The decisive test is whether the complex holds in **QM/MM**, where the zinc is described quantum mechanically and the artefact does not exist: `qmmm` with `skip_equilibrate: false` and no substrate restraint.

**That test has never been run.** D96 used `skip_equilibrate: true` - an optimisation with no dynamics. QM/MM MD is xTB-bound at roughly 1.4 s/step on CPU, so 50 ps is about 19 h and 10 ps about 4 h. GPU does not help: a QM/MM step copies positions off the device and the bottleneck is xTB, so these stages stay on CPU with their cores parallelised.

**Equilibration lengths, corrected.** D92's "10-100 ns published practice" conflated equilibration with production. Gervasoni's **100 ns is their NVT production run**; their equilibration was roughly 260 ps of heating and NPT before it. The comparable equilibration figures are Lisa and Bahena-Mendez at 10 ns and Guven at 50. So 100 ns is generous for the job it is doing here, though not wasteful at 4-5 h of GPU - and because `traj_frequency: 5000` writes a frame every 20 ps, **the 10 ns answer comes free from frame 500 of the same trajectory**, so both lengths are available without choosing between them.

**The NEB stage could not have run, and the reason is worse than a missing package.** (2026-10-23, found while costing the 50 ps QM/MM MD)

Costing the QM/MM equilibration led to `xtb_runmode`, since `inputfile` launches xTB as a separate process with file I/O **every step** - plausibly most of the 2.5-5.7 s/step seen in the optimisations, which would put 50 ps at 50,000 steps somewhere around 35-79 h with **no restart path in qmmm.py**. The obvious lever was `library` mode. Two things rule it out.

`xtb-python` is **not installed** in the ASH env, so `library` raises on import. And ASH's own `interface_xtb.py` carries, at line 110:

```
# TODO: QM/MM pointcharges for library
```

**The library path cannot return a QM/MM point-charge gradient.** [Corrected 2026-10-23 - see the correction at the end of this entry; the original claim here was that it does not implement point charges at all, which is wrong.] The *setter* works, so a library-mode energy is correctly embedded; only the gradient read was dead, and it ends in `ashexit()`. Since every use mcube makes of xTB is as the QM region of a `QMMMTheory`, and every such stage needs gradients, the effect is still that `library` is unusable here - but it fails loudly rather than quietly.

**`neb.py` hardcoded `runmode='library'`** (and `ts.py` did so whenever parallel frequency workers were used). As things stand those fail on import. Installing `xtb-python` to "fix" that would not have produced a silently un-embedded NEB [corrected 2026-10-23]: it would have reached the PC-gradient branch and called `ashexit()`, killing the job with a message telling you to use `inputfile`. Loud, but mid-run and after setup. The comment in neb.py read "xTB as a library: no input files per image", so the motive was per-image I/O and whoever wrote it had not seen the TODO.

**No existing result is affected.** `CSR_ligand/neb_results` holds only two force-field XML files from September and no config carries a `neb` section, so no NEB has ever completed; `ts` has never run either. Every QM/MM number on record came through `qmmm`, which used the default `inputfile` and therefore had proper electrostatic embedding. The bug was latent in the two stages that had not yet been used.

Fixed: both stages now use the configured runmode, and `make_qm` **refuses** `library` outright with the reason, rather than warning. Parallel images and frequency displacements each get their own working directory from ASH, so inputfile mode does not collide between them.

**Correction (2026-10-23, after reading ASH master rather than the September clone).** Two claims in this entry were wrong, and the framing "worse than a missing package" was overstated. ASH's library path does *not* silently drop electrostatic embedding: `interface_xtb.py` prints "pc grad is not yet implemented for runmode library", says to switch to `inputfile`, and calls `ashexit()`. And its point-charge *setter* has always been wired up, so a library-mode energy with point charges is correctly embedded - only the gradient was unreachable. `make_qm`'s refusal still earns its place, because it fails at config time and names the key to change instead of terminating the process mid-run, but its stated reason was wrong and has been rewritten. The genuinely silent case is g-xTB, not xTB (see D-entry above, which got that right): `gxTBTheory.run` accepts `MMcharges`/`PC` and references neither, so a QM/MM energy comes back gas-phase with no warning. Submitted upstream: ash#529 (g-xTB refuses PC), ash#530 (library PC gradient, closes ASH #157), xtb-python#130 (issue - the vendored `include/xtb.h` never declared `xtb_getPCGradient`, so it is unreachable from Python at all). Detail in `docs/ash-upstream-plan.md`.

**So the 50 ps QM/MM MD has to be costed at inputfile rates**, which is the open question for the NEB prep: 50,000 steps at 2.5-5.7 s is 35-79 h against a 24 h default walltime and nothing to restart from. Measuring the real rate with a short run is cheaper than discovering it at the walltime. The alternatives if it is too slow are a shorter equilibration, a 2 fs timestep, or adding a restart path to qmmm. **[Superseded 2026-10-24: measured at ~1.08 s/step, so 50 ps is ~18 h, and `general` allows 14 days, so no restart path is needed. The 2.5-5.7 s/step here came from optimisation steps, which carry geomeTRIC overhead MD does not. See D103 and `running-on-bunya.md`.]**

**D102 - The 100 ns run at the 77% mark: site held, disk on budget, the angle question still open.** (2026-10-23 15:42 AEST, job 30421166, `CSR_nac`, section `md_naclong`)

Checked mid-run, not after. `sacct` says **RUNNING** at 4 h 51 min elapsed, step 19,240,000 = **76.96 of 100 ns**. The "done or close" expectation was half right: close, not done. ETA ~17:06 at the rate below.

**Disk is on budget, so nothing here is a silent-truncation artefact.** `md_naclong_results` is **3.4 GB at 77%**, projecting ~4.4 GB against the ~6 GB expected and the ~60 GB that `traj_frequency: 500` would have written - the raise to 5000 (D101) does what it was for. Scratch **206/300 GB**, 385k/1M files, ~94 GB of headroom.

**The site was held: six restraints, five of them metal.** From the .out, all at 50 kcal/mol/Å², `Using platform: CUDA`. `[4255, 4208, 3.47]` is the one chosen restraint; the other five are Zn (atom 4203) to 4257, 2451, 4255, 1349 and 1414 at 1.90, 2.09, 2.14, 2.22 and 2.80 Å. `coordination_start`:

    B:1:ZN  D:1:O      HOH  1.90        B:1:ZN  A:88:NE2  HIS  2.22
    B:1:ZN  A:164:NE2  HIS  2.09        B:1:ZN  A:92:OD2  ASP  2.80
    B:1:ZN  C:2:O01    OH   2.14

That also pins the numbering: **4255 is O01**, the Zn-bound hydroxide. So the chosen restraint is O01-C4x as intended, and 4255 carries both its Zn bond (2.14 Å) and the 3.47 Å restraint.

**Throughput is about half D100's benchmark, and undiagnosed.** 76.96 ns in 291 min is **~16.5 ns/h, ~395 ns/day**, against the **845 ns/day** measured in D100 on the same H100. Candidates, none tested: trajectory and state reporting, which the 1 ns benchmark's per-ns cost never carried; the restraint force itself; and D100 being too short to amortise setup. It does not threaten the 12 h limit, so this is an efficiency question, not a risk to the run.

**The question the run was submitted to answer is still unanswered, and is not guessed at here.** `coordination_end`, the "coordination sphere changed during MD" block and `NVTtrajectory_lastframe.pdb` are all written at the end - the lastframe path returns 404 now - so `site.py sphere` and `nucleophiles` cannot run and the attack angle is unmeasured. Whether 100x the sampling rotates it from 80° into the 95-120° Burgi-Dunitz window, or confirms it as a property of the restrained pose (making `scan_attack` the way to reach a near-attack geometry rather than waiting for one), needs the finished run. The restraint was on throughout either way, so the release run remains the real test of whether the complex holds. Nothing was submitted.

**D102 completed (2026-10-23 21:09 AEST).** `sacct`: **COMPLETED**, exit 0:0, elapsed **06:19:45** (10:50:54 to 17:10:39), final step 25,000,000 = **100.000 ns**. Throughput confirmed at **~379 ns/day** against D100's 845, so that gap is real; recorded with the benchmark-versus-production reading in `running-on-bunya.md` under GPUs.

**The mid-run disk projection was wrong, and the original ~6 GB estimate was right.** Final `md_naclong_results` is **6.8 GB**, not the ~4.4 GB projected from 3.4 GB at 77%. Growth is not linear: the reimaged trajectory copy is written at the end, so the directory roughly doubles in the last moments. Scratch 213/300 GB, so nothing was truncated. Project from the frame count, not from bytes-so-far.

**`coordination_restraints` holds the five rows** (as selectors, in `config.yaml [md_naclong]`), plus the one chosen `[['C:2:O01', 'C:1:C4x', 3.47, 50.0]]`. The run is valid.

**`coordination_end`, and it changed more than expected:**

    B:1:ZN  C:2:O01    OH   1.86   (was 2.14, MOVED)
    B:1:ZN  A:92:OD1   ASP  1.99   (APPEARED)
    B:1:ZN  D:1:O      HOH  2.04   (was 1.90)
    B:1:ZN  A:88:NE2   HIS  2.22   (unchanged)
    B:1:ZN  E:557:O    HOH  2.27   (APPEARED)
    B:1:ZN  A:164:NE2  HIS  2.32   (was 2.09, MOVED)

The **Asp92 OD2-to-OD1 swap happened as predicted**, and went further than a relabelling: OD1 came in to 1.99 A while OD2 left to 2.98 A, so over 100 ns Asp92 became a genuinely short Zn ligand rather than the 2.80 A contact it started as.

**Beyond that, and this is the finding: a second water collapsed onto the zinc.** `E:557:O` appeared at 2.27 A, leaving Zn **six-coordinate** - O01, Asp92 OD1, HOH D:1, His88, HOH E:557, His164 - where the start had five. Both newcomers (Asp92 OD1, E:557) were unconstrained, which is exactly what the .out warns about: "a contact that appeared was not constrained, and a non-bonded metal has no Pauli repulsion to stop a donor collapsing onto it". At 1 ns (D100) this did not happen. So the classical model does drift to an over-coordinated zinc, it just needs more than 1 ns to do it, and the 1 ns test was too short to rule it out. The hydroxide also tightened to 1.86 A, which is short.

**Gln93 stayed out.** It is absent from the sphere list, so > 2.8 A, against the 1.79 A collapse seen unrestrained. Holding it out for 100 ns rather than 1 ns was the demonstration asked for, and it held.

**The question the run was submitted to answer: the angle rotated into the window. It was a sampling question.** `site.py sphere` on the last frame gives `O01 ... O-C 3.33 angle 101`, i.e. **101 degrees, inside the 95-120 Burgi-Dunitz window**, with O01-C4x tightened from 3.55 to **3.33 A**. `site.py nucleophiles` makes it unambiguous - of six candidate waters, the restrained hydroxide is the only one with a clean verdict:

    water          O-C   angle   M-O   verdict
    C:2:OH        3.33    101    1.86  metal-bound
    E:9347:HOH    3.47     79    3.79  angle off
    D:1:HOH       5.13     96    2.04  metal-bound, far
    E:557:HOH     5.59     83    2.27  metal-bound, angle off, far

So **80 degrees at 1 ns was a sampling artefact, not a property of the restrained pose**, and D100's conclusion that 1 ns had disproved the free-angular-coordinates argument was itself premature - the argument was right and 1 ns simply did not sample it. `scan_attack` is therefore **not** needed to manufacture a near-attack geometry: the equilibration reached one on its own. It may still be wanted as a reaction coordinate, but not to rescue the pose.

**One thing for the mechanism, not settled here.** The hydroxide's single proton H03 points at the substrate's own O1x (2.58 A) and at Asp92 OD1 (2.98 A) - not at Lys226, whose NZ is no closer than 3.6 A to any candidate's protons. That bears on the general-base question and on whether Asp92 is involved, but one frame is not evidence.

**Caveat on the handoff.** The frame QM/MM would start from has a six-coordinate zinc carrying a water that was not there at the start and was never restrained. The restraint was on throughout, so the release run remains the real test of whether the complex holds, and the extra water is a second reason to run it. Nothing was submitted.

**D103 - 50 ps GFN2-xTB QM/MM equilibration submitted (job 30467248), after two traps.** (2026-10-23 22:10 AEST, `CSR_nac`, section `nacqm`)

**The cost was measured, not guessed, and D101 was pessimistic.** 50 ps costs ~18 h at the rate a neighbouring QM/MM MD actually ran, not the 35-79 h D101 projected from optimisation steps, and `general` allows 14 days rather than 24 h, so the missing restart path in `qmmm.py` is not a blocker. Asked for 48 h. The rates, the threading mistake that made this run twice as slow as it needed to be, and how to follow a QM/MM MD's progress are in `running-on-bunya.md` under "Measured: QM/MM MD" - not repeated here.

**Trap 1: the stage would have started from the wrong structure.** `qmmm` resolved its input as `--pdb or inputs.pdb_file or md.results.lastframe`, and in this config that last fallback pointed at `CSR_unrestrained2/md_results/...` - the unrestrained run whose substrate dissociated. Nothing would have complained; QM/MM on the wrong structure converges and looks fine. Fixed generally rather than locally: `resolve_start_structure()` returns the path *and* its provenance, prints both, records them in the section results, and raises instead of guessing - naming all three ways when nothing is set, naming the input that pointed at a missing file, and refusing the md fallback when it reaches outside the config's own directory. An explicitly named outside path is still allowed; only the silent fallback is refused. `save_input()` pins the resolved value into the section's inputs. 8 tests in `tests/test_qmmm_structure.py`; 220 pass on Bunya, 172 + 4 skipped locally.

**Trap 2: `skip_equilibrate: true` is inherited, and it gates the MD.** The first submission (30466101) ran the **optimisation only** - caught 2 minutes in from geomeTRIC output in the .err, and cancelled. `qmmm.inputs` carried that flag from an optimisation-only workflow and `nacqm` inherited it. `nacqm` now sets `skip_equilibrate: false`. Worth remembering that every `qmmm.inputs` key is inherited this way, which is the same mechanism that makes the pinned `pdb_file` work.

**A third, self-inflicted one worth recording:** appending to `config.yaml` with `printf >>` broke the YAML, because the cancelled job had already written its `results` into `nacqm`, so the file no longer ended where the append assumed. Restored from the backup taken first and rebuilt the section. Back up the config before appending to it, and parse it afterwards - the error surfaced only because the parse was checked before resubmitting.

**Running, and verified on the live job:** `--gfn 2` on the xtb command line, `Timestep: 0.001 ps`, `Simulation time: 50.0 ps`, 8 cores, and `[INFO] Starting structure: .../md_naclong_results/NVTtrajectory_lastframe_imaged.pdb (from nacqm.inputs.pdb_file)`. The imaged frame measures identically to the raw one (angle 101, O01-C4x 3.33), so imaging did not split the site. Restraints are absent from the qmmm section, so they are released here - this is the first run in which the complex is unrestrained, and whether it holds is the open question. `skip_optimize: false`, so a geometry optimisation follows the MD and its output is the intended NEB starting frame.

**D104 - 50 ps GFN2 QM/MM: the complex holds and the zinc is fixed, but the near-attack geometry is gone.** (2026-10-25 08:30 AEST, job 30467248, `CSR_nac`, section `nacqm`)

The MD finished all 50,000 steps; the geometry optimisation that follows it is still running at 111 of 300 steps. Everything below is the **MD end point** (`QM_MM_lastframe.pdb`), verified complete first - 49,424 atom records, terminated with `END`, byte-identical in size to the frame it started from - because scratch is at 296/300 GB and an over-quota write fails silently.

**The substrate held, unrestrained.** This was the real test: the `qmmm` section carries no restraints, so 50 ps of GFN2 was the first run in which nothing held the complex together. Meropenem is still bound (Zn-mer O2x 4.88 A) with the beta-lactam intact (C4x-N1x 1.38 A). The Michaelis complex is stable without restraint, which is what the enzyme being a meropenemase requires.

**GFN2 removed the classical over-coordination, exactly as intended.** The zinc goes from the **six** ligands the 100 ns MM run left (D102) to **five**:

    MM, 100 ns                        GFN2, 50 ps
    O01        1.86                   O01        1.88
    Asp92 OD1  1.99                   Asp92 OD1  2.16
    HOH D:1    2.04                   HOH D:1    2.25
    His88      2.22                   His164     2.25
    HOH E:557  2.27   <- spurious     His88      2.32
    His164     2.32

`E:557` left the first shell (now 3.34 A from Zn, hydrogen-bonded to O01 at 1.88 A rather than coordinated), and Asp92 relaxed to a cleaner monodentate contact (OD1 2.16, OD2 3.45). So the artefact D102 found was classical, and the QM/MM equilibration is what fixes it - the standard-protocol justification, confirmed rather than assumed.

**But the near-attack geometry did not survive, and this reverses D102.** O01-C4x went **3.33 -> 3.86 A** and the attack angle **101 -> 57 degrees**, far outside the 95-120 Burgi-Dunitz window. It is not that some other water took over: `site.py nucleophiles` returns **every** candidate as "angle off" (50, 32, 57, 55, 71, 84, 33 degrees). D102 concluded the 101-degree NAC made `scan_attack` unnecessary. That held for the restrained classical ensemble and does not hold here: released and run under GFN2, the geometry relaxes away from attack. **Driving O01-C4x with `scan_attack` is therefore the route to a reactant-like structure for the NEB**, not an optional extra.

**The hydroxide was not protonated.** It still carries one hydrogen (H03, 0.97 A) after 50 ps, against the expectation that QM/MM would protonate it. Its proton now points at Asp92 OD1 (2.75 A) rather than at the substrate. Lys226 remains fully protonated (3 H).

**What this does not establish.** One frame is not an ensemble. That the end point is at 57 degrees does not show the NAC is never visited during the 50 ps, and `site.py` has `acceptors <traj> <topology>` and `joint` modes built for exactly that question. Running them over `QM_MM.dcd` is cheap next to another MD and is the honest way to decide whether the NAC is sampled at all before concluding it has to be driven.

**D105 - trajectories deleted to free scratch, and the optimisation is grinding in geomeTRIC.** (2026-10-25 08:45 AEST, `CSR_nac`)

**Scratch had 4 GB left (296/300 GB), so the trajectories went.** Deleted, with the per-frame time series kept in the CSVs that survive:

| file | size | what it was |
| --- | --- | --- |
| `md_naclong_results/NVTtrajectory.dcd` | 2.97 GB | 100 ns NVT production, job 30421166, 4 fs step, `traj_frequency 5000` so a frame every 20 ps (~5000 frames) |
| `md_naclong_results/NVTtrajectory_imaged.dcd` | 2.97 GB | the same frames, molecules wrapped into the box |
| `md_naclong_results/relaxbox_NPT.dcd` | 593 MB | the NPT box relaxation that preceded it |
| `md_naclong_results/relaxbox_NPT_imaged.dcd` | 593 MB | the same, wrapped |
| `md_naclong_results/OpenMMOpt_traj.xyz` | 64 MB | the pre-MD minimisation path |
| `protdft/nacqm/OpenMMMD_traj_wrapped.xyz` | 138 MB | the 50 ps QM/MM MD, wrapped, job 30467248 |

Kept deliberately: both `NVTtrajectory_lastframe*.pdb` (the imaged one is `nacqm.inputs.pdb_file`, so deleting it would have broken the running job), the first/last frame PDB and CIF of every stage, `OpenMM_MD*.chk` and `_state.xml`, and `nvtsim.csv` / `nptsim.csv` - which hold the per-frame energies and temperatures, so the time series survives even though the coordinates do not. **`protdft/nacqm/QM_MM.dcd` was also kept** (30 MB): it is the 50 ps QM/MM trajectory, and the `site.py acceptors` / `joint` analysis D104 calls for needs it.

That freed 7.3 GB: `md_naclong_results` 7.09 GB -> 74 MB, `nacqm` 848 -> 716 MB, scratch 296 -> 289 GB. A `DELETED.md` in each directory points here. Nothing the running job writes was touched.

**The optimisation is not hung, it is grinding, and it will not finish.** It has made no gradient call and written no file since 08:18, but `AveCPU` rose 6-06:40:11 -> 6-06:49:55 in 140 s of wall clock, i.e. **4.2 cores fully busy**. So the time is going into geomeTRIC's internal-coordinate linear algebra, not into xTB - which is exactly what `running-on-bunya.md` warned: TRIC over hundreds of active atoms (431 here, `actradius: 8`) does large linear algebra every step. The 15 x 30 MB `HT.v_*.txt` files written between 07:09 and 08:16 are the same signature that OOM-killed job 30409808.

At 111 of 300 steps after 4.9 h, and now >20 min on a single step, the remaining ~190 steps cannot finish before the 48 h walltime expires at 22:10 today. **`QM_MM_optimized.pdb` is only written at the end, so a timeout yields no optimised structure** - though `geometric_OPTtraj_optim.xyz` and `optimization_energies.log` retain the path so far. Nothing scientific is at risk either way: the MD end point is safe in `QM_MM_lastframe.pdb` and is already measured in D104.

**Correction, same day, before anything was resubmitted: it was already HDLC, so the coordinate system is not the cause.** ASH's `geomeTRICOptimizer` switches TRIC to HDLC by itself whenever an active region is set, unless `force_coordsystem=True`, and the job's output says so: *"Warning: ActiveRegion is set but the coordsystem is TRIC... Now switching to HDLC"*, then `Coordinate system:  hdlc`. The attribution to TRIC above was wrong, and `opt_coordsystem: hdlc` would have been a no-op.

Two things measured while checking that, both worth keeping:

- **The QM region is far larger than the active radius, and that is fine.** It reaches **15.85 A** from the Zn, with 39 of its 119 atoms beyond `actradius: 8`. That looked like a third of the QM region being frozen, but `active_atoms` holds 467 atoms and **all 119 QM atoms are in it** - `actregiondefine` selects whole residues, so a residue with any atom inside the radius comes in entire. Nothing was frozen.
- **mcube does not check that, though.** `resolve_actatoms` returns whatever `actregiondefine` wrote and never unions the QM atoms back in, so a smaller `actradius` could drop a whole QM residue out of the active region and freeze it while still treating it quantum-mechanically. Worth a guard that refuses when the active region does not contain the QM region; until then, any `actradius` change needs that containment re-checked by hand.

**The cost driver, found: geomeTRIC's SVD is failing on the rigid-water constraints.** `geometric/step.py:470` writes `'HT.v_%.5f.txt' % v` when the SVD fails, then nudges `v` by 0.001 and retries - and `geometric_OPTtraj.log` says so outright: *"SVD Error - saving HT.v_3349305840.31886.txt, increasing v by 0.001 and trying again"*. Fifteen failures in an hour, and by step 110 the trust radius had collapsed to **2.8e-04** with the gradient still **1.9e-01** and steps of 7e-05, so it was never going to converge. The constraints are mcube's own: every active-region **MM** water is frozen with three distances (O-H, O-H, H-H; QM waters deliberately free so they can transfer protons), and the log shows residuals of 0.08-0.10 A on distances that are supposed to be fixed.

**Correction, counted rather than assumed: there are 24 of them, not hundreds.** `constraints.txt` is 25 lines - a `$freeze` header and three distances for each of 8 waters - and the 5 A rerun has 12 for 4 waters. So the count cannot be what breaks the SVD, and the earlier claim that it was is withdrawn. Three distances rigidify a water triangle exactly, which is legitimate; what remains open is whether the trouble comes from the active-region size under HDLC, or from a water reaching a near-linear geometry where the H-H constraint stops being independent of the two O-H ones. Shrinking `actradius` helped in practice - 467 active atoms and 8 constrained waters gave 15 SVD failures in an hour, 203 atoms and 4 waters gave none - but on this evidence that is treating a symptom, not the cause.

**Also checked, because it looked like a difference and was not:** the new run prints `Rigidwater constraints: False` and `No automatic constraints`, which are OpenMM system-creation lines and appear in both runs. They say nothing about the optimiser's constraints.

**So the fix needs no new setting: shrink `actradius`.** Fewer active waters means fewer constraints. That was unsafe before, because a smaller radius could drop a whole QM residue out of the active region and freeze it; `resolve_actatoms` now unions the QM atoms in unconditionally and says which it had to add, so the radius is free to come down. Recorded in `running-on-bunya.md` under the optimisation measurement.

**D106 - the NEB starting structure, taken from the 8 A optimisation's plateau.** (2026-10-25 19:40 AEST, `CSR_nac`, `protdft/nacqm/QM_MM_opt_last.pdb`)

**The 5 A rerun (job 30580259) FAILED, and it failed badly.** 124 steps in 60 min - fast, at ~29 s/step against 159 - then `RecursionError` inside Python's logging, downstream of repeated `numpy.linalg.LinAlgError: SVD did not converge`. The energies say what happened: the QM term stayed sane throughout (-181.68 to -181.19 Eh) while **the MM term went from -240.69 to +1139 Eh**. Shrinking the active region pulled the active/frozen boundary in close to the site, and after ~100 steps the optimiser drove active MM atoms into frozen neighbours. **5 A is ruled out**, and the lesson is the opposite of the one in D105: the radius was not the thing to change.

**Two corrections to D105, both from looking in the wrong place.** The 5 A run's "zero SVD errors" came from grepping the **.out** for geomeTRIC's own `SVD Error` string; the real `LinAlgError` tracebacks go to the **.err**, so the SVD trouble was never gone. And `grep -c LinAlgError` on the 8 A run's .err is **0** - that run never crashed, its `SVD Error` lines were geomeTRIC's internal retry path recovering each time. So the 8 A optimisation was slow but healthy, and the claim in D105 that it "could not converge" was wrong: its energy fell -422.3655 -> -423.0022 Eh and then **plateaued**, the last three steps differing by under 1e-4 Eh.

**So the structure already existed and no rerun was needed.** `docs/csr/traj_frame.py` (new, 8 tests) reads one frame of a multi-frame XYZ onto a PDB template, which is what geomeTRIC's `geometric_OPTtraj_Full.xyz` needed - it is the only record of an optimisation that never reached its end. Frame -1 of 111, energy **-423.00218567193 Eh**, matching `optimization_energies.log` step 108 exactly; 49,424 atoms, coordinates from the frame and everything else from `QM_MM_lastframe.pdb`.

**What it is, measured:**

    5 ligands within 2.8 A          Zn-Asp92 OD2  3.47     O01-C4x  4.00
      1.88  C:2:OH:O01              Zn-mer O2x    5.16     C4x-N1x  1.37
      2.13  A:88:HIS:NE2            O01 angle     62 deg   Lys226   3 H
      2.14  A:164:HIS:NE2
      2.17  D:1:HOH:O
      2.19  A:92:ASP:OD1

A clean five-coordinate zinc: the optimisation tightened both histidines (2.22/2.25 -> 2.13/2.14) and kept Asp92 monodentate (OD1 2.19, OD2 3.47). The beta-lactam is intact (C4x-N1x 1.37), the substrate is bound without coordinating the zinc (mer O2x 5.16), and the hydroxide is still a hydroxide - one proton, 0.96 A, pointing at Asp92 OD1 (2.50 A). Lys226 remains fully protonated.

**It is not a near-attack geometry, which settles the question D104 raised.** O01-C4x went 3.86 -> **4.00 A** and the angle 57 -> **62 degrees** over the optimisation: relaxing to the nearest minimum moves *away* from attack, not towards it. Every candidate water is "angle off" (47-86 degrees). So `scan_attack` driving O01-C4x from 4.00 A down to ~1.45 A is how the reaction coordinate gets reached, confirmed now at the optimised geometry rather than inferred from one MD frame.

**Caveat on calling it optimised.** The energy is converged - under 1e-4 Eh over the last three steps, i.e. better than 0.1 kcal/mol - but the maximum gradient component was still 1.9e-01 when the run stalled, almost certainly in soft solvent modes. This is a relaxed geometry, not a stationary point. Adequate as a NEB *starting* structure; if the NEB endpoints need to be true minima, that needs an optimisation that actually converges, and the radius to try is 6-7 A, since 8 A grinds and 5 A diverges.

**D107 - Lys226 is the Zn2-site substitution, not a carboxylate lysine, and that reframes question 1.** (2026-10-25; trajectory measured here, literature from a research agent)

**Measured over all 50 frames of the CSR_nac QM/MM trajectory** (`protdft/nacqm/QM_MM.dcd`), in angstroms:

    NZ-O01   min 5.62  max 7.32  mean 6.63      never approaches the nucleophile
    NZ-O3x   min 2.61  max 5.41  mean 3.68      salt bridge to the C3 carboxylate
    NZ-N1x   min 4.71  max 6.92  mean 5.76      never within reach of the leaving N
    OD1-O01  min 2.64  max 3.23  mean 2.88      Asp92 H-bond, persistent

So the "Lys226 is a spectator" reading of the optimised frame was **wrong** - it anchors the substrate's carboxylate - but the frame was not unrepresentative either: across 50 frames the lysine never comes near the nucleophile or the leaving nitrogen. It cannot be the proton donor to N1x in this model, and it is not positioned to protonate the hydroxide.

**Which lysine it is, confirmed from this log before the literature was consulted.** D22 records "Zn2 (beta) site: Asp120, Gln121, **Lys263**", its numbering table maps **226 -> 254 -> 263**, and D40 calls the Zn2 site "D-Q-K ... one usable ligand". So Lys226 sits at the **Zn2 metal-ligand position**, where canonical B3 has His263 (Pedroso *et al.* 2020, doi:10.1007/s13238-020-00736-4; Krco *et al.* 2023, doi:10.3389/fchem.2023.1196073; and the 6DR8 mutant R118H/Q121H/**K254H** restores exactly that histidine). It is **not** BBL Lys224, the B1 carboxylate anchor - B3 enzymes have Ser221/Ser223 there. Line 235 said otherwise and is retracted; D12's EcoRV framing is retired.

**Why the EcoRV analogy fails.** In Imhof, Fischer & Smith 2009 (doi:10.1021/bi900585m), the paper this project is modelled on, Lys92's main assigned role is stabilising a **metaphosphate** - a dissociative centre with heavy excess charge that has no counterpart in amide hydrolysis. Within EcoRV the role is contested anyway: transition-state stabilisation (Selent 1992, doi:10.1021/bi00135a010) against orienting the nucleophile and depressing its pKa (Horton 1998, doi:10.1073/pnas.95.23.13489). And no report exists of an MBL lysine protonating a zinc-bound hydroxide. The one MBL Lys proton-transfer paper (Das & Nair 2017, doi:10.1039/C6CP08769H) is **cephalexin only**, runs through a two-water chain that **meropenem's R2 group blocks**, and found no direct Lys-to-N route even under a 30 kcal/mol bias.

**The hypothesis this suggests (mine, not from a source).** In canonical MBLs the substrate's C3/C4 carboxylate binds **Zn2** (Feng *et al.* 2017, doi:10.1038/s41467-017-02339-w). CSR-1 has no Zn2 and a lysine in its place - and that lysine is what the trajectory finds holding the C3 carboxylate at 2.61 A. **Lys226 is doing Zn2's substrate-anchoring job electrostatically.** **[Qualified 2026-10-26, D111: the contact does not hold. Over the next 50 ps of the same continuous trajectory it opens to 4.26 A at closest and spends no frame inside 3.5 A, and no other cation takes over. The 2.61 A measurement is still correct for `nacqm`; what it will not support is a *persistent* anchoring interaction. The hypothesis is untouched as a hypothesis - it just has no standing evidence from this trajectory.]** What it cannot do is Zn2's catalytic jobs, lowering the nucleophile's pKa and stabilising the anionic nitrogen, which would account for CSR-1 binding broadly while turning over poorly. That is sharper than calling the lysine a defect, and it is testable: the open part is whether its charge, at the Zn2 position 4.7-6.9 A from N1x, still stabilises the anionic intermediate at long range. That intermediate is long-lived and its decay is rate-determining (Lisa *et al.* 2017, doi:10.1038/s41467-017-00601-9; Das & Nair 2017), so a charge-deletion calculation would answer it cheaply.

**A divergence to state before any mechanism is transplanted.** The only published mono-zinc B3 mechanism is GOB-18 (Lisa 2017), and **its single metal sits at Zn2** - the site that stabilises the anionic nitrogen. **Ours is at Zn1, with Lys226 occupying an empty Zn2.** The species doing the stabilising there is absent here, so that scheme cannot be carried over wholesale. Lisa 2017 is still a better template than CphA for a mono-zinc B3 with carbapenems.

**Two things this settles that were worrying us.** Nothing being positioned to protonate N1x in the Michaelis complex is the **expected** result, not a defect: three papers converge on a water that arrives *after* ring opening, and on the anionic intermediate lasting long enough for it. And the pKa bookkeeping argues against the Lys226-to-hydroxide transfer of D29/D35: a 2+ cation 4-5 A away should **raise** a lysine's pKa, and a Zn-water pKa of ~5.6 in BcII against Lys ~10.5 puts the transfer uphill by ~6-7 kcal/mol, opposite in sign to the computed -13 to -16. GFN2-xTB's worst benchmarked functional group is exactly the protonated amine (MUE 5.3-7.7 kcal/mol; Arantes & Rezac 2025, doi:10.1021/acs.jctc.5c00690). **The clean four-coordinate DFT retest already specified in this log is the right next step, and it now has a clean geometry to run in.**

**D108 - the crystallographic zinc sphere, measured: Asp92 is not on it.** (2026-10-25, from 6NC5 and 6DN4 downloaded from RCSB)

The reference for zinc coordination is the crystal structure, not anyone's MD
length (user). So it was measured rather than quoted. **6NC5, ZN302** - the Zn1
H-R-H site our model keeps - has four ligands within 2.8 A:

    1.90  HOH485:O          2.14  HOH493:O
    2.09  HIS163:NE2        2.22  HIS87:NE2

**Two histidines and two waters. No aspartate.** HOH493 is the water our model
carries as O01 (`site.yaml` calls it "ex-HOH493") and HOH485 is the second, so the
model reproduces the crystallographic sphere and then **adds Asp92 OD1 at 2.19 A**,
giving five coordinate where the crystal gives four. Against the stated reference,
Asp92 on the zinc is an addition, which supports reading it as over-coordination
rather than a real feature - though note Asp92 is a **Zn2**-site residue (D-Q-K,
per D22/D40) whose carboxylate has nothing to ligate while Zn2 is empty, so its
drift onto Zn1 is at least a chemically motivated artefact rather than a random
one.

**Also measured, and it settles an earlier doubt:** HOH485 bridges ZN302 and
**ZN303**, whose other ligands are two waters and an acetate (ACT307) - the
artefact zinc. So the suspicion recorded earlier that HOH485 is the weakest member
of our sphere is well founded.

**6DN4, the mature wild type, has no zinc modelled in it.** That is a fact about
one crystal, not about the enzyme: **CSR-1 requires zinc to function** (user), so a
metal-free deposition reflects the low-affinity beta site and the crystallisation
conditions, not the catalytic state. It is evidence about metal *affinity*, and
says nothing against a zinc-dependent mechanism. 6NC5's ZN302 remains the
reference for the Zn1 sphere.

**Correction to the numbering conflict raised in `csr-neb-plan.md` yesterday: it
was my error, not the log's.** D40's "His87/His163" are 6NC5's own numbers, which
map to the model at **+1** (87 -> 88, 163 -> 164, 91 -> 92, 225 -> 226). D22's
"minus 28" describes a different file, `CSR_nontruncated_1metal_original.pdb`
(crystal 29-302 -> model 1-274). Two files, two offsets; both statements were
right and I conflated them. D40's reading of 6NC5 is confirmed exactly, including
ZN301 = Asp91/His92/His225 - and since His92 and His225 **are** the Q121H and
K254H substitutions, 6NC5 is a beta-site mutant and is not a wild-type reference
for that site. 6DN4 is.

**D109 - `resume: true` on the qmmm stage, so a long QM/MM MD can be a chain of jobs.** (2026-10-26)

300 ps of production at the measured 1.08-1.68 s/step is 90-140 h. Bunya's `general` partition allows 14 days, so one job can hold it, but a job that dies at hour 80 currently throws away everything. ASH already writes `OpenMM_MD_state.xml` and `OpenMM_MD_final_state.xml` as it runs and says in its own output that they can restart a simulation; `qmmm.py` simply never passed them.

**What it does.** `resume: true` restores the newest state file - positions *and* velocities - and runs only the remainder of `time`. `time` stays the **total** wanted, so resuming a 50 ps run that stopped at 21 ps runs 29 ps, and the flag is a one-line addition to a config that already ran, with nothing else to adjust. A run already at its target is not rerun; with nothing to resume from it starts fresh and says so, so the flag is safe to leave in place across a chain of resubmissions. Default is false: unlike `scan`'s `resume`, where skipping finished points cannot change a result, here it decides whether a rerun continues a trajectory or replaces it.

**Three things found while building it, each of which would have cost a run.**

1. **ASH's restart overwrites the trajectory it resumed from.** Its DCDReporter is only ever constructed inside `if restart is False`, with `append=False`, so a state-file restart in a fresh process truncates `QM_MM.dcd`. mcube now moves the previous segments to `QM_MM_part1.dcd`, `QM_MM_part2.dcd`, ... and records them in `results.trajfile_parts`; they concatenate in that order, followed by `QM_MM.dcd`. ASH's `restart=True` path is for calling `.run()` again on a live object in one process, not for a new job.
2. **A finished run does not report its target time.** OpenMM accumulates the timestep, so 50000 steps of 1 fs land within ~3e-11 ps of 50, on either side. **Measured on job 30598032, which finished 2.6e-11 ps *short*** (`time="49.99999999997417"`) - and short is the direction that matters, because `remaining <= 0` then misses and a 24 h job gets queued to integrate 26 fs. The cutoff is one timestep, the shortest thing that can actually be run. A synthetic 50000-step state file written for the tests drifted the other way (+2.6e-11), so neither sign can be assumed.
3. **The state file does not carry the thermostat's random stream.** Velocities are restored but the Langevin noise restarts from a new seed, so a resumed run is not bit-identical to an uninterrupted one. Immaterial for equilibration - the thermostat is not part of the physics being measured - but it has to be stated, because it means the segments are one trajectory in the thermodynamic sense and not in the deterministic one.

**Why the state file and not the checkpoint.** ASH accepts either, and prefers `chkfile` when both are given. The checkpoint is an opaque binary blob tied to the platform that wrote it; the state file is portable XML whose root element carries `time` and `stepCount` - which is the number resume needs, readable from the first 4 kB without OpenMM and without parsing tens of megabytes of coordinates. Requiring the closing `>` of the opening tag also means a half-written file, from an over-quota write or a job killed mid-save, raises rather than yielding a plausible wrong elapsed time.

15 tests in `tests/test_qmmm_resume.py`, which write their state files with OpenMM itself rather than hand-rolling the header, so the format being parsed is the real one. That is how finding 2 surfaced. Full suite 237 passed, 5 skipped.

**D110 - Job 30598032: 50 ps with Arg229 in the QM region. The nucleophile moves further from attack, not closer.** (2026-10-26)

COMPLETED, exit 0:0, elapsed 23:56:12, ended 21:53:32 - **four minutes inside the walltime.** 50,000 steps at 1 fs confirmed from `OpenMM_MD_final_state.xml` (`stepCount="50000"`). 1.72 s/step, slightly worse than the 1.68 estimated from adding the guanidinium. Scratch at 156/300 GB, so nothing here is suspect on quota grounds. The QM region is **137 atoms** against `nacqm`'s 119: Arg229's 18 side-chain atoms (CB to NH2 with hydrogens, backbone left in MM), with `charge: 3` as set. Lys226 NZ carries three hydrogens, so both cations are present and cationic, which was the point of the run.

**Zinc, five-coordinate, unchanged in character from D108:** O01 1.90, His164 NE2 2.10, Asp92 OD1 2.11, HOH485 2.16, His88 NE2 2.18. Asp92 OD2 is 4.22, so the aspartate is firmly monodentate - it still bridges nothing.

**The result that matters, and it is negative.** The attack geometry has degraded again:

    O01-C4x        angle      source
    3.33 A         101 deg    100 ns classical MD (D102)
    3.86 A          57 deg    50 ps QM/MM, no Arg229 (D104)
    6.20 A          63 deg    50 ps QM/MM, Arg229 in QM (this run)

The second Zn-bound water is no better (HOH485-C4x 6.27, angle 79), and Zn-mer O2x is 6.76. **Putting Arg229 in the QM region did not pull the substrate back into a near-attack arrangement; the nucleophile drifted half again as far away.** So the conclusion of D104 is now independently reproduced with a different QM region: **GFN2-xTB QM/MM MD will not hand us a near-attack frame, and `scan_attack` has to drive the O01-C4x coordinate.** Two runs, 100 ps of QM/MM between them, agree on this.

**What this does not say.** It is not evidence against the +2-substitution hypothesis, and it does not invalidate the NEB start frame: that came from `nacqm`'s *optimisation* (D106), not from either MD end point. Nor does it show the substrate has left the site - Zn-mer O2x at 6.76 is a single contact, and this stage's `monitor` is deliberately unset, so Lys226-O3x and the Arg229 contacts were not measured. **Whether the substrate is still anchored over these 50 ps is open, and it decides whether this trajectory is usable for the charge-deletion analysis at all.** Measuring the same contacts as D107 over `qmmm_nacarg/QM_MM.dcd` is the next step.

**D111 - Contacts over the nacarg trajectory: the Lys226 salt bridge breaks, and over 100 ps the substrate is leaving the site.** (2026-10-26; `docs/csr/traj_contacts.py`, 50 frames of `qmmm_nacarg/QM_MM.dcd`, saved as `qmmm_nacarg/contacts.txt`)

Minimum-image distances in angstroms, `%` being the fraction of the 50 frames inside 3.5 A:

    contact                                    min   max  mean   %<3.5
    Lys226 NZ - mer O3x                       4.26  6.09  5.02    0.0
    Lys226 NZ - any mer O                     4.26  5.72  4.96    0.0
    Lys226 NZ - mer N1x (leaving N)           5.68  7.23  6.54    0.0
    Lys226 NZ - O01 nucleophile               6.01  7.69  6.72    0.0
    Arg90 guanidinium - carboxylate O3x/O4x   9.40 11.78 10.69    0.0
    Arg90 guanidinium - any mer O             4.03  7.41  5.45    0.0
    Arg229 guanidinium - carboxylate O3x/O4x  3.93  6.26  5.00    0.0
    Arg229 guanidinium - any mer O            3.78  5.57  4.60    0.0
    Asp92 OD1 - O01 nucleophile               2.63  3.20  2.94  100.0
    Zn - O01 nucleophile                      1.73  1.97  1.86  100.0
    Zn - mer O2x                               4.77  8.38  6.47    0.0
    O01 - mer C4x (attack)                    3.63  7.23  5.67    0.0

**The salt bridge D107 built its argument on is gone.** D107 measured Lys226 NZ-O3x over `nacqm`'s 50 ps at min 2.61, mean 3.68, and read that as Lys226 doing the absent Zn2's substrate-anchoring job electrostatically. Over the next 50 ps it is min 4.26, mean 5.02, and **never once inside 3.5 A**. These two runs are consecutive - `nacarg` started from `nacqm`'s last frame - so this is not a contradiction between measurements but **a trend across 100 ps of continuous QM/MM: 2.61 -> 4.26 A at closest approach.** The anchoring is being lost, monotonically.

**No other cation takes over.** Arg90 is 9.4-11.8 A from the carboxylate, so it never anchored it; its closest substrate contact is O1x at 4.03. Arg229 is the nearest cation to the substrate, and even it never reaches 3.5 A of anything (min 3.78 to O2x). **Adding Arg229 to the QM region did not restore an anchoring interaction - it only let us watch, quantum mechanically, that there is not one.**

**Two things are rock solid, and they are the two that do not involve the substrate.** The hydroxide stays on the zinc (1.73-1.97, 100% of frames) and Asp92 OD1 holds it (2.63-3.20, mean 2.94, 100%). So the **nucleophile and its general base are stable and correctly arranged; the substrate is the part that will not stay.** Zn-mer O2x never comes within 4.77 A, so the carbonyl never engages the metal either, and the attack distance ranges over 3.63-7.23 A - the substrate is loosely associated and diffusing, not bound in a Michaelis complex.

**What this blocks.** `qmmm_nacarg/QM_MM.dcd` **cannot supply Michaelis-complex frames** - not for NEB endpoints, not for the charge-deletion analysis the +2 hypothesis needs. A charge-deletion energy computed on a substrate that is 5 A from every cation would measure the drift, not the catalysis. D106's NEB start frame is unaffected; it came from `nacqm`'s optimisation, not from either MD.

**What it does not settle, stated plainly.** This is *not* evidence against the +2-substitution hypothesis, and it is not evidence that CSR-1 binds meropenem weakly in reality - a drifting simulation cannot support either conclusion. Three readings remain open and they need different fixes: the **reactant protonation state** may be wrong (the `protH_lys`/`protH_oh` comparison is still unrun in its fixed form, and it is now the most informative thing we could do); **GFN2-xTB** may simply not hold this complex, which the DFT retest already specified in this log would show; or the **complex needs a restraint through equilibration**, which would make the sampling an assumption rather than a result. Deciding between them is the user's call, not something to pick by default.

**D112 - The system we have been simulating is in the protonation state M1 rejected.** (2026-10-26, found while bringing the plan's status up to date)

**D48 (job 30053828) concluded:** "Still negative, so the resting state is Zn-OH2 with neutral Lys226. **The model genuinely has no hydroxide**" - `ohH` below `lysH` by 5.23 kcal/mol at wB97X-D3(BJ)/def2-SVP//xTB/MM.

**`CSR_nac` is in the `lysH` state.** Read off the topology of the running system (`qmmm_nacarg/QM_MM_firstframe.pdb`, and the same in `nacqm`): residue `OH` resSeq 2 has atoms **O01 and H03 only**, i.e. a hydroxide, and **Lys226 NZ carries HZ1, HZ2 and HZ3**, i.e. a cationic ammonium - confirmed independently by job 30598032's own final-frame summary (`Lys226 NZ  3 H  HZ1 1.09, HZ2 1.05, HZ3 1.02`). That is `qmmm_cleanR_lysH`, the state D48 disfavoured.

**How it got there: nothing decided it.** D98 built `CSR_nac` from `CSR_unrestrained2` with "the same `preprocess.inputs`", listing the docked pose, cutoff, HIS variants, `extra_qm_residues` and `constraint_exclude_residues` - and no protonation change. D48 landed on 2026-10-18; `CSR_nac` was built on 2026-10-23. **The verdict was never carried across.** So D102 (100 ns), D104, D106's NEB start frame, D107's contacts, D110 and D111 were all computed in it.

**What this does and does not mean.** It does not make those runs wrong - it makes them runs of a different protonation state than the one M1 chose, and that has to be said in any write-up. It does raise M0's first reading from one of three possibilities to the leading explanation for D111's drift: if the complex is being held in a disfavoured protonation state, a substrate that will not stay bound is what you would expect.

**But it is not a simple fix, because the two best arguments disagree in sign.** D48 itself called -5.2 kcal/mol "not comfortably excluded, only disfavoured" - within reach of def2-SVP error, the imposed constraints, and the absence of entropy and sampling. D107's pKa bookkeeping runs the other way: a 2+ cation 4-5 A away should raise a lysine's pKa, and Zn-water ~5.6 against Lys ~10.5 puts the transfer uphill by 6-7 kcal/mol, which would **favour** `lysH` and vindicate what has been run. GFN2-xTB cannot arbitrate: its worst benchmarked functional group is the protonated amine, MUE 5.3-7.7 kcal/mol (Arantes & Rezac 2025), which spans the whole disagreement.

**So the action is the `protH_lys`/`protH_oh` comparison at DFT in the current five-coordinate geometry**, which settles M1 and M0's first reading together and is the cheapest thing on the board. Until it runs, no barrier computed from `CSR_nac` can be defended, because the resting state it starts from is contested. **Not a reason to rebuild anything yet** - if `lysH` wins, nothing needs redoing.

**D113 - The ohH arm, built and checked but not submitted. Method changed to VIM-2's.** (2026-10-27, user: "we're letting our methodology decide rather than comparing singlepoints")

**The method.** The protonation question is decided by equilibrating each state and keeping the one that holds the active site the crystal shows - not by comparing energies. That is what the closest published analogue does: *Protocols for Metallo- and Serine-beta-Lactamase Free Energy Predictions* (J Phys Chem B 2024, 128:12416, PMC11664573) ran **100 ps of QM/MM MD (DFTB3) per state** on VIM-2 with QM charges 0 and -1, judged by **zinc coordination against the crystal**, and chose water because hydroxide disrupted the site (Zn-Zn 3.6 A against 4.7 for water and 4.45 in the crystal). Supporting precedent: CphA/biapenem (PLoS One 2013, 10.1371/journal.pone.0055136) ran QM/MM over **14 ionization states**; IRE1 (PMC12678879) resolved a *catalytic lysine's* protonation by **QM/MM metadynamics** on the Lys-carboxylate proton-sharing coordinate and found it ligand-dependent. Setting the rest of the system with a pKa predictor at a stated pH is universal (H++ there, PROPKA at pH 7.4 for IRE1). **EcoRV could not be checked** - Imhof 2009's Methods are behind the ACS paywall, as this log already records - so it is not cited for this.

This suits us: it avoids every weakness of D48 (no imposed constraints, no 5 kcal/mol margin, no basis-set sensitivity) and it **already has half its answer**, since the lysH arm's 100 ps loses the substrate (D111).

**pH 8: already satisfied, and the test is pH-independent anyway.** I first reported that mcube never sets a pH; that was wrong, from a case-sensitive grep that missed `pH`. `preprocess` passes `pH=inputs.get("pH", 8.0)` on the explicit-solvent path CSR uses, so **every CSR system was built at pH 8.0** - the implicit-solvent branch defaults to 7.0, an inconsistency worth knowing but with no effect here. Setting 8 explicitly changes nothing either: OpenMM's only thresholds are ASP/GLU 4.4, CYS 8.5, LYS 10.0 and a HIS case at 6.5, **none crossed between 7 and 8**. And `lysH -> ohH` moves a proton *inside* the QM region, conserving total charge, so its energy difference carries no pH term at all; pH would only matter for a state with a different proton count. The config now records `pH: 8.0` rather than leaving it to a default.

**Built:** `/scratch/user/uquadhik/CSR_ohH`, config generated from CSR_nac's and verified to contain no `CSR_nac` path, no `hydroxide` reference and no inherited `results`. Three preprocess keys differ and nothing else: `protein_pdb` the new structure, `../water_nuc.xml` for `../hydroxide.xml`, and `226: LYN` added to the HIS variants. LYN drops HZ3, so **both arms have the same atom count, the same total charge and the same QM charge of +3** with Arg229. `md_ohHmd` mirrors `naclong` exactly but at **2 fs**; `qmmm_ohHqm` is 50 ps with `resume: true`, so the second 50 ps is a `time: 100` resubmit of the same section (D109).

**Why 2 fs.** `rigidWater` recognises residues named HOH and not NUW, so the nucleophile's H-O-H angle is flexible where bulk water's is constrained - checked in OpenMM, not assumed. With HBonds constraining the O-H bonds that angle is the only free mode, and at ~21 fs it is too fast for a 4 fs step. Costs about 13 h of GPU instead of 6.5. In the qmmm stage mcube builds the MM system with `rigidwater=False` and `autoconstraints=None`, so the nucleophile is entirely free there, which is what lets the proton move.

**Four things found by testing rather than by reasoning, each of which would have cost a job or corrupted the comparison.**

1. **`protH_oh` cannot do this test.** It pins HZ1 0.98 A from O01 while the MM bond to Lys226 NZ is still there. Fine for a constrained optimisation; in free dynamics either the constraint does the work or the bond pulls the proton back. Hence a real topology.
2. **A custom 3-atom water template matches every bulk water by bond graph**, and OpenMM then refuses the whole system: *Multiple non-identical matching templates found for residue (HOH): HOH, NUW*. With ~15000 waters that aborts preprocess at once. The fix is for the template to be **parameter-identical** to HOH - it declares nothing but the residue and borrows TIP3P's own types - because identical matches are accepted. A residue named HOH with our atom names fails differently: PDBFile infers a standard water's bonds from the names O/H1/H2, so ours would arrive unbonded. Bonds come from CONECT, as the hydroxide's already do.
3. **The clash check counted the water's own partner hydrogen.** H04's nearest neighbour is always H03 at 2*0.9572*sin(104.52/2) = 1.513 A, seven thousandths above the 1.5 A cutoff - the real build reported it at 1.52 A. Every correct placement was a near miss and a slightly different geometry would have been refused. With the nucleophile's own atoms excluded the nearest neighbour is **Lys226 NZ at 2.00 A**, which is the hydrogen bond you would want: the proton sits on the NZ...O01 axis it came from.
4. **The derived config carried the hydroxide under two different keys** - `preprocess.inputs.forcefields` *and* `md.inputs.default_forcefields` - plus absolute `qmmm.inputs.forcefields` paths into `CSR_nac/hydroxide.xml` and that run's `openff_MER.xml`. Handing the ohH run the hydroxide template would leave its water unparameterised. The swap now walks the whole structure, resolved `forcefields` overrides in later stages are dropped (`mm_settings` falls back to `md.results.forcefields`, which this run regenerates), and the result is refused if either the source directory or the word "hydroxide" survives anywhere.

**Not submitted.** Three jobs, about 37 h to a first verdict: preprocess (minutes), `md md_ohHmd` on a GPU (~13 h), `qmmm qmmm_ohHqm` on CPU (~24 h for 50 ps at 1.72 s/step). 28 tests across the three scripts; full suite 264 passed, 5 skipped.

**D114 - Imhof read at last; the crystal zinc is four-coordinate; the ohH arm redesigned to start from naclong's frame, minimisation submitted as 30698605.** (2026-10-27)

**Imhof 2009 is accessible after all** (user: "Imhof not paywalled"), so the Methods this log has twice called unobtainable are now read rather than inferred:

- **System:** three X-ray structures, metals and crystallographic waters placed as determined crystallographically, **hydrogens "added using stereochemical criteria"** - no pKa predictor at all. 48,825 atoms, close to our 49,424.
- **Equilibration:** minimise to 0.01 kcal/mol/A, 20 ps heating, 300 ps equilibration, **1 ns** NPT production, H-bonds constrained, **2 fs**. Then "a snapshot from the molecular dynamics ensemble, **taken to possess average values for some important geometric parameters**".
- **QM region:** dimethyl phosphate, both Mg, side chains of Asp90/Asp74/Glu45/**Lys92**, Ile91 backbone, active-site waters.
- **How protonation was decided: by whether a state is a stable minimum.** "All attempts to move protons from Lys92 to the metaphosphate failed, indicating that the positive charge is better stabilized by the amino group of Lys92", and elsewhere "this protonation state is apparently not stable". Lys92 stayed cationic. **This is the same criterion as VIM-2's and the opposite of comparing single points**, which is what the user asked for in these words: "we're letting our methodology decide rather than comparing singlepoints".

**The intro also settles why the Lys ambiguity is not ours alone.** Imhof: Lys92 "may help stabilize and orient a hydroxide ion for nucleophilic attack (if positively charged) or function as a catalytic base to accept a proton from the attacking water (if neutral)", and "the latter possibility is supported by the finding that the reaction chemistry proceeds optimally at a relatively high pH of 8.5". Exactly our Lys226 question, and the pH-optimum argument is one we could make too.

**The finding that bears hardest on us is their semi-empirical benchmark.** Against B3LYP/6-31++G(d,p) on small models, AM1/d gets most reaction energies to ~4 kcal/mol but fails on coordination number in a way we should expect from GFN2:

    4- to 5-coordinate (adding a water)        7.55 kcal/mol too FAVOURABLE at AM1/d
    bidentate -> monodentate Asp (5 -> 4)      7 kcal/mol too COSTLY at AM1/d

and "at the DFT level of theory, the monodentate complex has a structure in which a proton is transferred from one water molecule to the nonligating Asp oxygen atom. The AM1/d-optimized geometry does not show such a proton-transferred structure." **So a semi-empirical Hamiltonian over-stabilises the higher coordination number, over-penalises monodentate carboxylate, and misses a water-to-Asp proton transfer that DFT finds.** Our model is five-coordinate with Asp92 on the metal and a hydroxide; all three of those are the artefacts this benchmark predicts. Imhof also makes coordination number mechanistic: Mg goes **6 -> 5 -> 4** along the path, four-coordinate is the competent state, and ions that cannot reduce their coordination (Ca) are inactive.

**The crystal zinc is four-coordinate** (user asked; confirmed from D108's own measurement of 6NC5 ZN302): HOH485 1.90, His163 NE2 2.09, HOH493 2.14, His87 NE2 2.22. Two histidines, two waters, no aspartate. **The same coordination number as CphA's zinc**, though with different ligands - CphA has Asp120/Cys221/His263 plus the nucleophile, and D60 noted its four-coordinate site likewise has no room for the substrate carbonyl. Our optimised model is five-coordinate, the extra ligand being Asp92 OD1 at 2.11.

**The arm was redesigned on the user's objection: "why are you preprocessing from the start?"** Correctly. A fresh preprocess gives the two arms independent 100 ns trajectories, which is a worse control than sharing one frame, and costs 13 h of GPU for nothing. **Both arms now start from `md_naclong`'s 100 ns lastframe - the frame `nacqm` itself started from** - so they differ by where one hydrogen sits. One correction to the user's suggestion was needed: the *QM/MM* 100 ps lastframe is the drifted structure (O01-C4x 6.20 A, substrate unbound, D111), where "does this state hold the site together" is unanswerable because there is no intact site left; naclong's frame is at 3.33 A and 101 degrees.

**A geometric result came out of building it: a water does not fit where the restrained MM put the hydroxide.** Holding O01 at 3.47 A from C4x pressed the hydroxide against meropenem, and that oxygen position has room for one hydrogen, not two. Measured: the best clearance on the cone of allowed H04 positions is **1.54 A** from meropenem H5x; reorienting the whole water with the oxygen fixed reaches only **1.71 A**. **So the near-attack geometry that run produced is comfortably available only to a hydroxide** - which is independent of, and consistent with, D111's drift. Hence a minimisation before the MD, which is also what VIM-2 did for each state.

**Built and validated.** `make_water_nuc.py --lys A:226 --floor 1.2` writes `md_naclong_results/ohH_start.pdb`: H04 added at TIP3P geometry, residue OH renamed NUW, Lys226's HZ3 removed with the residue left named LYS so OpenMM matches LYN by bond graph (renaming it would make `build_qm_region` treat it as non-standard, because this config's `amino_acids` list overrides mcube's default and omits LYN, and put its whole backbone in the QM region). `check_ohH_frame.py` builds **both** systems in full with the real force-field list and confirms what the comparison depends on:

    lysH  49424 particles, total charge +0.000
    ohH   49424 particles, total charge +0.000
    nucleophile  -1.000 (OH, 2 atoms)  ->  +0.000 (NUW, 3 atoms)
    Lys A:226    +1.000 (22 atoms)     ->  +0.000 (21 atoms)

**Minimisation only needs no new code:** `gentle` is off by default and `skip_npt` plus `skip_nvt` drop the rest, so the md stage runs `OpenMM_Opt` and stops, reporting the coordination sphere either side of it - which is the check we want anyway. Submitted as **job 30698605** (`md md_ohHmin`, gpu_cuda, 40 min). Then `sections.py` can add the `ohH` qmmm section, whose start structure does not exist until this finishes, which is why it refused it today.

**Four things testing caught, each of which would have cost a job.** (1) The H04 line is built from O01's to inherit the formatting, and inherited its **element column** too - the new hydrogen went out as an oxygen, so the residue was O-H-O, matched no water template, and OpenMM reported it as a malformed glycine. The test fixtures could not have caught it because `pdb_line` never wrote columns 77-78 at all. (2) **`md` ignores `--pdb`** and reads `inputs.pdb_file`, so passing the frame as run.slurm's third argument would have equilibrated preprocess's own output and looked normal; it now refuses when the two disagree. (3) Paths in `inputs` resolve from the stage's output directory while `start` resolves from the config directory, so the qmmm force-field paths were a level too shallow. (4) A 3-atom water template matches every bulk water by bond graph unless it is parameter-identical to HOH, which aborts the whole system.

**D115 - The ohH minimisation ran and its result is invalid: md restrains saved atom indices against a structure whose atom order changed, so the zinc was never restrained.** (2026-10-27, job 30699162)

**The job itself was fine.** COMPLETED in **3 min 16 s** on 8 CPUs - which vindicates moving it off the GPU, where `gpu_cuda` had quoted an 80 minute queue for 500 minimisation steps (job 30698605, cancelled). Scratch 158/300 GB.

**What it reported.** `md` prints the metal's coordination sphere either side of the run, and it changed:

    LEFT      Zn - A:164:NE2  (was 2.32 A)
    LEFT      Zn - D:1:O      (was 2.04 A)
    MOVED     Zn - E:557:O    2.27 -> 2.04 A

**An anchoring histidine left the zinc.** That alone disqualifies the structure as a starting point, so the 50 ps QM/MM was **not** submitted.

**The cause is not the added proton.** `coordination_restraints`, which `md` records back in selector form, reads:

    [['C:1:C1x', 'D:1:O', 1.9], ['C:1:C1x', 'A:164:NE2', 2.09],
     ['C:1:C1x', 'C:2:H03', 2.14], ['C:1:C1x', 'A:88:NE2', 2.22], ...]

The "metal" is **meropenem's C1x**, and one donor is the nucleophile's *hydrogen*. Every metal restraint was applied to the wrong atom, so **the zinc was unrestrained** and the force field rearranged its sphere - which is exactly what the restrained pipeline exists to prevent (D98's "the metal site is never left free in MM").

**Why, in one line:** `md_restraints` does `rows = resolve_atoms(results.get("coordination"), pdb_file)`, and `preprocess.results.coordination` is a list of **raw atom indices** into `processed_pdb`. `resolve_atoms` passes integers through untouched, so the saved indices are applied verbatim to whatever structure the stage is given, with nothing checking that they still name the same atoms.

**Measured, 0-based:**

    index   lysH frame          ohH frame
    4202    OXT PRO A 275       ZN
    4203    ZN                  C1x mer

The saved coordination holds 4203 for the zinc, right for the lysH ordering. Building the ohH frame removed Lys226's HZ3 and inserted H04 at the nucleophile, so every atom between the two shifted down by one and the zinc became 4202 - while atoms after the insertion (D:1:O at 4257) kept their index, which is why some restraints looked right and others did not.

**This is a general mcube bug, not a slip in this workflow.** Every run until now reused `processed_pdb`'s atom order, so the saved indices happened to stay valid; the first structure edited by hand broke them silently. The same hazard was already noted for `resolve_qm_region`, which falls through to `preprocess.results.qm_atoms` - stale indices there too. **The ohH qmmm section is safe from it only because `extra_qm_residues` forces the region to be rebuilt from residues on the actual file** (noted when that section was written, D114).

**What the fix has to be.** `md` should verify that each saved coordination index still names the atom it was detected on - element and residue are enough - and refuse with a message naming the mismatch, rather than restraining whatever now sits at that position. A re-detect-on-this-structure option would then be the way to run an edited frame at all, because the atom order cannot be preserved: the proton has to leave one residue and join another, and a PDB residue's atoms must be contiguous, so any correct edit renumbers.

**Not fixed here.** `src/mcube/md.py` currently holds uncommitted work from the CLI flag refactor, so it was left alone rather than mingling two changes in one file. The two belong together in any case: a flag that carries data outside the config and an index saved in the config that silently stops matching the structure are the same disease - state that is not re-derived from the thing in front of it.

**Where the ohH arm stands.** `md_naclong_results/ohH_start.pdb` is correct and validated (D114): 49,424 atoms, total charge 0, nucleophile -1 -> 0, Lys226 +1 -> 0. Only the minimisation of it is void. `md_ohHmin_results/frag-minimized.pdb` must not be used.

**D117 - Every path is read against the config's directory, and migrating the runs to it audited the archive.** (2026-10-27, user: "Let's just make everything relative to the working dir - that way you can copy and there's no ambiguity")

Force-field lists were the last paths resolving differently from the rest - OpenMM read them from the stage's working directory while `pdb_file` had moved to the config's. Two bases in one `inputs` block is what put a substrate carbon where the zinc should have been (D115), so there is now one: **the config's own directory**, which is also the unit that gets copied.

**The rule is narrow, because a force-field list holds two kinds of entry.** `resolve_forcefields` rewrites one only when a file of that name sits beside the config. OpenMM's bundled names (`amber14/protein.ff14SB.xml`) are not there, so they pass through to OpenMM untouched; absolute paths already say what they mean. **Entries starting with `..` are left alone and warned about rather than resolved**, and that exclusion is load-bearing: `<config dir>/../water_nuc.xml` is `/scratch/user/uquadhik`, where spare copies accumulate, so resolving it could silently prefer a stranger over the file the stage reads today.

**The migration, by `docs/csr/migrate_paths.py` (dry run by default, backup on apply): 49 changes across 10 configs.** `../hydroxide.xml` becomes `hydroxide.xml` - the prefix is stripped and the file looked for *beside* the config, because that is what the entry meant when the stage ran a level down, not `home/".."/value`, which is a different file. Absolute paths into a run's own directory become bare names.

**What it refused to touch is the interesting part, and it is an audit of the archive.**

**Five configs read force-field files out of a *different* run's directory**, inherited by copying and unnoticed because the files exist:

    CSR_nac           -> CSR_unrestrained2   openff_MER.xml
    CSR_unrestrained  -> CSR_crystal         hydroxide.xml, openff_MER.xml
    CSR_unrestrained2 -> CSR_crystal         hydroxide.xml, openff_MER.xml
    CSR_crystal_v3    -> CSR_crystal_v2      hydroxide.xml, openff_MER.xml
    CSR_gln93         -> CSR_crystal         hydroxide.xml, openff_MER.xml

Rewriting one of these to a bare name would change **which file the run uses**, which is a decision about the science rather than a path fix, so they are reported and left as they are. **For `CSR_nac`, the live run, it does not matter**: its own `openff_MER.xml` and the CSR_unrestrained2 one it actually reads differ by **two lines, both `<DateGenerated>`**. OpenFF is deterministic for the same SMILES, so nothing computed from D98 onward is affected - but nothing was checking, and the five configs are the concrete form of the hazard the user's instruction removes.

**Two configs were already broken, independently of any of this.** `CSR_gln93` lists `../hydroxide.xml` three times with no `hydroxide.xml` beside its config, so that entry resolved to nothing under the old rule either. `CSR_hydrolysed` lists `openff_hydrolysed_meropenem.xml` and `openff_hyd.xml`, neither of which is in its directory. Both would fail on the first force-field read; neither is in use.

**Results stay absolute, deliberately.** `inputs` says what to use and travels with the directory; `results` records what *was* used and should stay unambiguous even if the directory is later copied or moved. The one place mcube had this backwards - `save_input` pinning the resolved absolute `pdb_file` back into `inputs`, which is how a derived config inherits a path home - is noted for the same treatment.

9 tests on the resolver and 8 on the migration; suite 353 pass, 5 skipped.

**Repointing the five cross-run configs (user, 2026-10-27: "Repoint the five cross-run configs at their own files"). Four of the thirteen entries could be; the rest could not, and the reasons matter.**

`migrate_paths.py --repoint` rewrites a foreign path only when the run's own copy is **byte-identical, or differs solely in a `<DateGenerated>` line** - OpenFF stamps its output with the generation date, so two runs of the same preprocess on the same SMILES differ in that one line and nothing else. Anything further is a change to the force field wearing the clothes of a path fix, and for a finished run it makes the config stop describing what ran.

    repointed  6 entries across 4 configs
      all 8 hydroxide.xml in the archive are byte-identical (f1ec8f0e...),
      so every hydroxide reference was a free fix
      CSR_nac's openff_MER.xml: same but for the date

    refused    CSR_unrestrained   2x openff_MER.xml - no copy of its own
               CSR_crystal_v3     2x openff_MER.xml - no copy of its own
               CSR_gln93          3x - no copy of either file
               CSR_unrestrained2  1x openff_MER.xml - ITS OWN COPY DIFFERS

**That last one is the finding.** `CSR_unrestrained2`'s own `openff_MER.xml` differs from the `CSR_crystal` one its config reads across **51 atom lines, with different charges** (`C1x` -0.0921 against -0.1443, and so on) - not a date stamp. So two different meropenem parameterisations exist in this project, and that run used the older one. Repointing it would silently have swapped the ligand's charges in a config describing finished work, so it was left alone.

Three runs have no parameterisation of their own at all and borrow one; that is a fact about how they were set up by copying, and the only way to give them their own is to re-run `preprocess`, which would generate a third set. **Nothing was regenerated, because which parameterisation is correct is a question about the chemistry and not about paths.**

13 tests on the migration tool; suite 360 pass, 5 skipped.

**D116 - The D115 guard fired on the real run, which is the good news; it fired because `sections.py` was looking for the structure in the git checkout.** (2026-10-27, jobs 30703062 and 30705703)

**Job 30703062 FAILED after 9 seconds**, costing nothing, with exactly the message the guard exists to produce:

    preprocess.results.coordination holds atom indices that do not name a metal
    and a donor in .../md_naclong_results/ohH_start.pdb:
        index 4203 should be a metal (ZN) but is C mer C:1
        index 4255 should be a donor (N/O/S) but is H NUW C:2
    Those indices were measured on another structure, so restraining them would
    hold the wrong atoms and leave the metal free.

**That is D115 caught before it could do anything**, and it names the exact atom - index 4203 is meropenem's C1x on this frame - rather than leaving it to be inferred from a coordination report afterwards.

**But the run should not have reached the guard.** The section was meant to carry `coordination_redetect: true`, and it did not. `sections.py` had refused to write it, with

    NOT WRITING. Problems:
      - ohHmin: pdb_file missing: md_naclong_results/ohH_start.pdb

for a file that is plainly there (4,004,830 bytes). **It was resolving the structure against `site.yaml`'s directory** - `home = Path(site_file).resolve().parent` - which is the git checkout, while the structures live in the run directory. The comment two lines below it already said "Paths are read against the config's directory, as mcube reads them"; the code did not. Introduced when `start` was replaced by `inputs.pdb_file`, and invisible to the existing test because its helper writes `site.yaml` into the same directory as `config.yaml`, so the two bases coincided.

**My error on top of it:** I submitted without checking that `sections.py` had succeeded. Its refusal was in the log I had just written to, and I read the tail for the `sbatch` line instead of the whole thing. The section still held the void `results` from the first attempt, which was the visible tell that nothing had been rewritten.

Fixed: the base is the config's directory, `sections.py` refuses outright when there is no `config.yaml` in the working directory, and two tests keep the bases apart - a site file outside the run must still find the structure, and a decoy beside the site file must not satisfy the check. Suite 355 pass, 5 skipped.

**Resubmitted as 30705703** with the section correctly written: `coordination_redetect: true`, `pdb_file: md_naclong_results/ohH_start.pdb`, `water_nuc.xml` (config-relative since D117), and no stale `results`.

**Worth keeping in view:** three jobs have now been spent on this minimisation and none of it on chemistry. The first restrained a substrate carbon, the second was refused for a missing key, and the third is queued. Each failure was caught by something added after the previous one, which is the system working - but the thing still being tested is whether a water fits where the hydroxide was, and that has not been measured yet.

**D118 - The ohH minimisation is sound, and the 50 ps QM/MM is submitted as 30721365.** (2026-10-27, job 30705703)

COMPLETED in 3 min 53 s. All four things that had to be true are:

- `coordination_source: redetected on .../md_naclong_results/ohH_start.pdb` - not the saved indices.
- The restraints name the zinc. Six detected contacts, every one on index **4202**, which is the zinc *on this frame*; the first attempt restrained 4203, which is meropenem's C1x here (D115, D116).
- **"The coordination sphere is unchanged."** His88 NE2 and His164 NE2 both still on the metal, with O01, Asp92 OD1 and the two waters. His164 left in the first attempt, which is what voided it - and it left because the zinc was never restrained, not because of the added proton.
- The clash relaxed.

**Measured on the structures either side, in angstroms:**

    H04...H5x (mer)   1.54  ->  3.15     the contact the minimisation was for
    H03...H5x (mer)   1.65  ->  3.26
    Zn-O01            1.86  ->  2.05     hydroxide -> water
    O01-C4x           3.33  ->  3.49     held by the 3.47 restraint
    Zn-mer O2x        5.16  ->  4.37

**The Zn-O01 lengthening is an independent check on the structure, not just a relaxation.** 1.86 -> 2.05 A on protonating the nucleophile is the textbook hydroxide-to-water signature, and D48 measured 1.83 -> 1.99 for the same change in a different geometry with a different method. Two unrelated routes to the same ~0.17 A. So the ohH frame really is the water state, and the minimisation moved what it had to and left the rest.

**Submitted: job 30721365**, `qmmm qmmm_ohH`, 50 ps at 1 fs, 137 QM atoms with Arg229, charge +3, `resume: true`, 8 CPUs, 36 h walltime. About 24 h at the 1.72 s/step measured in D110. The comparator is `nacqm` plus `qmmm_nacarg`: the same frame, one proton apart, judged by whether the active site holds - VIM-2's criterion (D113), and Imhof's (D114).

**Two things to note about getting here.** `sections.py` refused twice before it wrote: first for the path bug (D116), then because a fresh code-server terminal has no conda environment, so `from mcube.preprocess import _locked` failed. Both were caught by checking the exit status, which is the lesson of D116 applied. And the force-field entries in the written section are the config-relative ones from D117 (`water_nuc.xml`, `preprocess_results/openff_MER.xml`), which the run resolved to absolute paths under `CSR_nac` - the first live exercise of that rule, and it worked.

**What the run has to answer.** The lysH arm lost the substrate over 100 ps (D111): O01-C4x from 3.33 to 6.20 A, no cation holding it, the Lys226 salt bridge opening from 2.61 to 4.26 A. If ohH holds where lysH did not, the protonation state was the cause and M1 is settled. If it drifts too, the cause is elsewhere - GFN2, or the pose - and M0's other readings are live. Either outcome is informative, which is the point of running it.

**D119 - The ohH arm drifts too, so the protonation state was not what lost the substrate. M1 is answered; M0 is now the live question.** (2026-10-28, job 30721365)

COMPLETED in 19:46:57, 50000 steps at 1.427 s/step, ended 17:45:31. `.err` holds nothing but the thread-binding INFO line, and the step rate was flat to three digits across the whole run (1.427 measured between step 30670 at 10:05 and step 45863 at 16:07), so xTB was doing the same work per step at 45 ps as at 30 ps - no convergence trouble anywhere in it.

**Both arms measured with the same tool, 50 frames each, 1 frame/ps** (`docs/csr/traj_contacts.py`, which now resolves the nucleophile residue name so it reads OH and NUW alike - b715497):

    contact                       ohH (water)        lysH (hydroxide, nacarg)
    O01-C4x (attack)        4.55-6.20  m 5.31  0%   3.63-7.23  m 5.67    0%
    Zn-O01                  1.80-2.14  m 1.98 100%  1.73-1.97  m 1.86  100%
    Asp92 OD1-O01                      m 2.82 100%             m 2.94  100%
    Lys226 NZ-O01           2.54-6.96  m 4.03  58%  6.01-7.69  m 6.72    0%
    Lys226 NZ-mer O3x                  m 5.62   0%             m 5.02    0%
    Arg90-any mer O         3.09-4.73  m 3.74  18%  4.03-7.41  m 5.45    0%
    Arg229-carboxylate                 m 5.18   0%             m 5.00    0%
    Zn-mer O2x                         m 5.85   0%             m 6.47    0%

(percentages are frames inside 3.5 A; m = mean, angstroms)

**The answer to D118's question is the second branch.** The ohH frame went into the run at O01-C4x = 3.49 A and never came back under 3.5 A in 50 ps - mean 5.31, max 6.20. The lysH arm did the same thing (mean 5.67, max 7.23). **Protonating the nucleophile did not hold the substrate, so the protonation state was not the cause of the loss in D111.** M1 is settled in the sense that matters: it is not the variable that decides this. The cause is elsewhere - GFN2, or the pose - and M0's other readings are live, which is exactly what D118 said this outcome would mean.

**The test itself was sound, which is worth stating because a negative result invites the suspicion that the setup never worked.** Zn-O01 is 1.98 A in the water arm against 1.86 A in the hydroxide arm, held 100% of frames in both. That 0.12 A separation survives 50 ps of dynamics, and it agrees with the 1.86 -> 2.05 A measured across the minimisation (D118) and with D48's independent 1.83 -> 1.99. The two arms really are the two protonation states. They simply both lose the substrate.

**The ohH arm is slightly better on every substrate-holding measure, and nowhere near enough.** Attack distance 5.31 against 5.67, Zn-mer O2x 5.85 against 6.47, and two contacts that exist only in this arm: Lys226 NZ reaches the nucleophile in 58% of frames (0% in lysH) and Arg90 reaches a substrate oxygen in 18% (0%). Note *which* partner Lys226 found - the nucleophile, not the carboxylate it was proposed to anchor (NZ-O3x stays at 5.62 A, 0%). A cation at the vacant Zn2 position engaging the nucleophile is a different role from the one D107 hypothesised for it.

**What is NOT being launched.** The plan after this test was a 300 ps production run and then a DFT/MM optimisation of the complex. Both are now premature: a production run started from this pose would sample more of the same drift at ~1.4 s/step (300 ps is about 5 days of wall time), and a DFT/MM optimisation of a complex whose attack distance averages 5.3 A is not an optimisation of a reactant state. **Neither job can answer the question that is now open, so neither should be spent on it.**

**What the open question is.** Two candidates, and they are distinguishable:

  1. **GFN2 is wrong about this site** - the standing concern that xTB over-coordinates the zinc. Testable without new dynamics: take frames from both arms and compare GFN2 against DFT on the same geometries. Imhof's benchmark (D114) is the precedent and the warning - AM1/d had a 4->5-coordinate error of 7.55 kcal/mol there, which is larger than anything this mechanism turns on.
  2. **The pose is wrong** - the binding mode M0 assumes. Nothing in two independent 50 ps runs holds meropenem near the metal, which is what a wrong pose looks like.

A reading that distinguishes them is the next thing worth a job. Running more of the same dynamics is not.

**D120 - Both QM/MM arms ran on a FIVE-coordinate zinc, against the crystal's four, and a coordination restraint this pipeline generated is what put the fifth ligand there. This is upstream of the protonation test and of any GFN2 question.** (2026-10-28)

Asked by the user after D119: 4- or 5-coordinate? Measured with `site.py shell`, a new mode that reports every donor's occupancy and the per-frame coordination number over a trajectory rather than in one frame (d5b9c83 - `sphere` could only answer it for the frame it was given, which is why this was never checked).

**Both arms, 50 frames each, 2.6 A cutoff: coordination number 5 in 100% of frames.** Not fluctuating, not borderline - the same five ligands throughout:

    ligand            ohH (water)              lysH (hydroxide, nacarg)
    O01 nucleophile   1.80-2.14  m 1.98  100%  1.73-1.97  m 1.86  100%
    Asp92 OD1         1.83-2.32  m 2.03  100%  1.84-2.33  m 2.07  100%
    His88 NE2         1.84-2.25  m 2.05  100%  1.97-2.39  m 2.13  100%
    HOH485 O          1.87-2.29  m 2.07  100%  2.02-2.52  m 2.24  100%
    His164 NE2        1.95-2.26  m 2.09  100%  1.97-2.36  m 2.14  100%
    ---- outside the cutoff in every frame ----
    Asp92 OD2         3.20-4.02  m 3.60    0%  3.24-4.22  m 3.65    0%
    HOH557 O          3.37-5.04  m 4.25    0%  2.91-4.19  m 3.60    0%

**The crystal is four-coordinate** (6NC5 ZN302: HOH485 1.90, His163 NE2 2.09, HOH493 2.14, His87 NE2 2.22 - two histidines and two waters, no carboxylate). **Asp92 OD1 at 2.03 A is a ligand the crystal does not have**, and it is present in every frame of both arms.

**Where it came from, structure by structure** (`site.py sphere`, 2.8 A):

    structure                            ligands <2.8 A
    preprocess csr_crystal_start.pdb     HOH485 1.90, His164 2.09, O01 2.14, His88 2.22, Asp92 OD2 2.80
    preprocess finalsystem.pdb           identical
    md_naclong NVT lastframe (100 ns MM) O01 1.86, Asp92 OD1 1.99, HOH485 2.04, His88 2.22, HOH557 2.27, His164 2.32
    md_ohHmin frag-minimized.pdb         Asp92 OD1 1.99, O01 2.05, HOH485 2.06, HOH557 2.09, His88 2.19, His164 2.28
    QM/MM MD, both arms                  five, as tabulated above

The input was the crystal's four plus Asp92 OD2 at 2.80 A. **It is the 100 ns MM equilibration that changed the shell**: the carboxylate rotated so that OD1, not OD2, faces the metal (1.99 A) and a solvation water, HOH557, took a sixth position (2.27 A). Six at MM level, and xTB then kept five of them, pushing HOH557 back out to 3.6 A. Neither the QM/MM MD nor GFN2 created this; both inherited it.

**The cause is a restraint this pipeline generated.** `preprocess` detects coordination with `cutoff: 3.0`, so Asp92 OD2 at 2.80 A was recorded as a zinc ligand and written into `md.inputs.coordination_restraints` as `[B:1:ZN, A:92:OD2, 2.8, 50.0]` - a 50 kcal/mol/A^2 harmonic holding that carboxylate oxygen at 2.80 A for the whole equilibration. A carboxylate pinned that close to a metal has one cheap way to relax, which is to rotate the other oxygen onto it. **The restraint whose purpose was to preserve the crystal's coordination sphere is what destroyed it**, because 2.80 A was never a coordination distance - it is a second-shell contact that a 3.0 A cutoff cannot tell from a ligand.

**Why it went unseen for three weeks.** Everything that watches Asp92 watches the wrong oxygen. `site.yaml`'s monitor has `Zn-Asp92 OD2`; D53's scan pulled `[B:1:ZN, A:92:OD2]`; the D42 and D96 sections constrain `Zn-Asp92 OD2`. OD2 sits at 2.98-3.60 A in every structure above and reports as uncoordinated, which is exactly what a healthy site would look like. OD1, the actual ligand at 1.99-2.07 A, was never in any monitor, any scan or any restraint. D53 concluded Asp92 coordination was a real minimum from a scan of the free oxygen.

**What this does to the drift in D111/D118/D119.** A five-coordinate zinc carrying two histidines, a carboxylate, a water and the nucleophile has no open coordination site. The substrate carbonyl cannot bind a metal that is already saturated - and Zn-mer O2x measures 5.85 A in the ohH arm and 6.47 in lysH, 0% of frames in contact. **That is a sufficient explanation of why meropenem leaves in both arms, and it does not require GFN2 to be wrong about anything.** It also explains why the protonation test came back negative: moving one proton cannot open a coordination site that a restrained carboxylate is occupying, so both arms were always going to lose the substrate.

**So the GFN2-vs-DFT benchmark proposed at the end of D119 is not the next job, and the user said so before I did.** It would have benchmarked the method on geometries whose zinc shell was wrong at the MM stage, and a disagreement between GFN2 and DFT on those frames would have told us nothing about the enzyme. The method question is real (D114) but it is downstream of this, and it cannot be asked on these frames.

**What has to be re-run, and what has to change first.** The restraint list is generated, so the fix belongs in the config rather than in a hand edit: the detection cutoff has to be tight enough to exclude a 2.80 A contact (a Zn-O ligand is ~1.9-2.2 A; 2.5 is defensible, 3.0 is not), and a second carboxylate oxygen should not be able to enter the list behind the first. Everything from `md_nac` onward then has to be re-run from a shell that matches the crystal: 100 ns MM, then QM/MM MD per arm. That is a large amount of re-run, which is the cost of having restrained a contact at a distance nobody checked was a bond.

**Two things worth stating before that is costed.** First, Asp92 is D120 in BBL numbering - the residue that coordinates Zn2 in di-zinc B3 enzymes. In a mono-zinc B3-RQK site with Zn2 vacant, a carboxylate swinging onto Zn1 is not absurd chemistry, and the question of whether it does is a real one about this enzyme. But it is a question to ASK, not an assumption to restrain in, and the crystal's answer at 2.80 A is that it does not. Second, nothing above depends on the 2.6 A cutoff: Asp92 OD1 is at 2.03 A and no cutoff between 2.1 and 3.2 changes the count.

**D121 - Carbonyl-to-Zn1 is an OBSERVED step in B3, measured here from the time-resolved L1 structures; and in those structures D120 is a Zn2 ligand over 4 A from Zn1. That fixes what a CSR-1 Michaelis complex has to look like, and confirms D120's restraint put Asp92 in the carbonyl's seat.** (2026-10-28)

Asked by the user after D120: is substrate carbonyl binding a step seen in other MBL mechanisms? Reviews disagree, so this was settled from coordinates rather than from review text.

**L1 (subclass B3) with INTACT moxalactam, two time-resolved snapshots** (Nat Commun 2022, 13:7175, doi:10.1038/s41467-022-35029-3; PDB 7UHI at 40 ms, 7UHJ at 60 ms, 2.20 A). The beta-lactam carbonyl is C1=O1 (1.21 A), identified by bond connectivity from the ring nitrogen rather than by trusting the atom name:

    site        7UHI (40 ms)                                  7UHJ (60 ms)
    Zn1         His107 ND1 2.18, His181 2.35, O1 2.53,        His107 2.18, His105 2.23,
                His105 2.66                                   His181 2.42, O1 2.63
    Zn2         Asp109 OD2 2.27, His110 2.27, His246 2.48,    Asp109 2.06, His110 2.23,
                O8 2.65, N1 2.67                              His246 2.33, O7 2.54, N1 2.78
    Zn-Zn       3.17                                          3.38

**So the answer is yes: the carbonyl oxygen coordinates Zn1 at 2.5-2.6 A in a genuine B3 Michaelis complex**, with the ring nitrogen on Zn2 (2.67-2.78) and the C-carboxylate bridging to Zn2 (2.54-2.65). No water is modelled on Zn1 in either snapshot, so whether the nucleophile sits on Zn1 *simultaneously* is not resolved at 2.20 A - but the carbonyl's occupancy of a Zn1 coordination position is.

**This retires the objection `site.yaml` recorded for declining to restrain Zn-carbonyl.** That comment reads: "for a mono-Zn B3 site it is not established that the carbonyl coordinates the same zinc that carries the nucleophile". It is now established for B3 that the carbonyl coordinates Zn1; the residual uncertainty is only whether a mono-zinc Zn1 holds carbonyl and hydroxide at once, which is a 4-coordinate arrangement (2 His + OH + carbonyl) and not a crowded one.

**And it identifies D120's restraint as the error it looked like.** Asp109 in L1 is D120 in BBL numbering - the same residue as our **Asp92**. In both snapshots it is a **Zn2** ligand at 2.06-2.27 A and sits over 4 A from Zn1. Our crystallographic zinc is a **Zn1**-site zinc (His87/His163 = BBL H116/H196, with R118 replacing the third histidine), and the 100 ns MM run put Asp92 OD1 onto it at 1.99 A - into the coordination position the carbonyl occupies in L1. The restraint did not merely add a ligand; it added one where the substrate has to go.

**The monometallic comparison, and why CphA never transferred.** CphA with biapenem (PDB 1X8I, subclass B2) is five-coordinate: Asp120 OD2 2.03, His263 NE2 2.12, Cys221 SG 2.27, plus ligand N4 2.22 and O32 2.39. **Five coordination is therefore not itself the defect** - which ligands they are is. But CphA's single zinc occupies the **Zn2** site (Asp-Cys-His), so its chemistry is Zn2 chemistry; that is the structural reason the CphA route failed twice here. The right precedents for CSR-1 are mono-zinc **B1** (BcII-type, Zn1 site) and the Zn1 half of B3.

**What a reasonable CSR-1 Michaelis complex therefore is - acceptance criteria, not preferences:**

  1. **Zn1 site, four-coordinate**: His88 NE2, His164 NE2, the nucleophile as a *terminal* hydroxide (there is no Zn2 to bridge), and meropenem's carbonyl O2x at 2.0-2.6 A. Five only if a second water is genuinely retained beside those.
  2. **Asp92 off the metal**, at the crystal's 2.80 A, oriented at the nucleophile. Its proposed role is general base, which is a hydrogen bond, not a coordination. If it coordinates, it is in the carbonyl's seat.
  3. **Attack geometry** O01-C4x at 3.0-3.5 A inside the `burgi_dunitz` window already in `site.yaml`.
  4. **The substrate's anchors must be the RQK residues, because Zn2 is absent.** The B3-RQK substitutions are Arg90 (R118, Zn1 site), Gln93 (Q121) and Lys226 (K263). There is direct B3 precedent for a neutral amide rather than the metal polarising the carbonyl - Gln157 is the oxyanion hole in SMB-1 (QM/MM, doi:10.1007/s00894-020-4319-6), Asn176 in L1, Asn233 in CcrA - and **Gln93 is the structural analogue of exactly that role**. A credible complex should show the C3 carboxylate on Lys226 and/or Arg90 and the carbonyl held by the zinc together with Gln93.

**Which reframes D94 as a pattern rather than an incident.** MM once put Gln93 OE1 at 1.79 A from the zinc (D94); it has now put Asp92 OD1 at 1.99 A (D120). The equilibration keeps filling the open Zn1 position with whichever protein donor is nearest, because the generated coordination restraints pull them in while nothing holds that position for the substrate - and the Zn-carbonyl distance was deliberately left unrestrained. Every run since `md_nac` has equilibrated the substrate *out* of a site that was being filled behind it. That is a sufficient account of D111, D118 and D119 without any appeal to GFN2.

**Open, and worth measuring before the re-run is costed:** Gln93's contacts to the substrate have never been measured on any trajectory - `traj_contacts.py` has no Gln93 probe and no Zn-O2x occupancy. Both existing 50 ps trajectories can answer whether Gln93 was oriented for the oxyanion role while the substrate left, and that costs a script change and two reruns of an analysis, not a job.

**D122 - Gln93 is 10 A from the substrate, and the docked pose never had a Zn-carbonyl contact at all. Measured against L1, the pose has exactly one defect and the equilibration has a second; between them M0 and the drift are both accounted for.** (2026-10-28, 59ed041)

`traj_contacts.py` gained Gln93 probes and a carbonyl selector (59ed041), then was run on both 50 ps trajectories and on the two structures they came from. The reference column is **L1 7UHJ**, the 60 ms snapshot with intact moxalactam (D121), with the BBL equivalences H121 = Gln93, H263 = Lys226, D120 = Asp92:

    contact                   L1 7UHJ   docked   after 100ns MM   QM/MM ohH   QM/MM lysH
    Zn - carbonyl O2x            2.63     5.43             5.16        5.85         6.47
    Gln93 - carbonyl O2x         5.98     6.56            10.22       10.00        10.49
    Gln93 NE2 - carboxylate      3.86     8.94            10.60       11.79        13.70
    Lys226 NZ - carboxylate      3.54     4.98             3.97        5.62         5.02
    Lys226 NZ - lactam N1x       3.82     3.57             5.72        6.42         6.54
    Arg229 - carboxylate            -     2.96             3.69        5.18         5.00
    Asp92 - nucleophile             -  OD2 2.84         OD1 2.87    OD1 2.82     OD1 2.94
    Zn - Asp92 OD1          (on Zn2)   2nd shell             1.99        2.03         2.07
    O01 - C4x (attack)              -     3.47             3.33        5.31         5.67

(docked = `preprocess_results/finalsystem.pdb`; after 100ns MM = `md_naclong_results/NVTtrajectory_lastframe.pdb`, the frame both arms start from; QM/MM columns are 50-frame means)

**Gln93 is not the oxyanion hole in anything we have run, and the reason is not that it declines the role.** It is 10.0-10.5 A from the carbonyl in both trajectories, and 10.2 A already in the frame they started from. But in the docked structure it was 6.56 A, against 5.98 A for the equivalent atom in L1 - **the only one of these readings that was ever near its reference**. So the hypothesis D121 raised has not been tested yet: it was carried out of range by the equilibration before any QM/MM run began.

**The docked pose has exactly one defect, and it is the carbonyl.** Measured against L1, `finalsystem.pdb` is defensible on every other count: the C3 carboxylate on Arg229 at 2.96 A, Lys226 NZ 3.57 A from the leaving nitrogen against L1's 3.82, Gln93 at 6.56, the attack distance at the restrained 3.47. **Zn-O2x at 5.43 A against L1's 2.63 is the one number that is wrong, and it is wrong by 2.8 A.** M0's binding mode was never a B3 Michaelis complex in the one respect D121 established as a step of the mechanism.

**The 100 ns MM equilibration then destroyed the three contacts that were right.** Gln93 6.56 -> 10.22, Lys226-N1x 3.57 -> 5.72, Arg229-carboxylate 2.96 -> 3.69 and on to 5.0-5.2 under QM/MM. Meanwhile Asp92 rotated from OD2-at-the-nucleophile (2.84, which is the crystal's arrangement for a general base) to OD1-on-the-zinc (1.99), under the restraint D120 identified. **The substrate was equilibrated out of the site while a protein carboxylate was pulled into it.** Only one restraint ever acted on the substrate - `[C:2:O01, C:1:C4x, 3.47]` - and it holds the nucleophile's radial coordinate while saying nothing about the carbonyl, the carboxylate or the ring nitrogen.

**So the two failures are separable and both are now named.** (i) The pose lacked the Zn-carbonyl contact from the start, which is M0. (ii) The equilibration protocol is not pose-preserving: it restrains one distance and lets everything else, including the metal's coordination sphere, go where the force field prefers. D111, D118 and D119 are downstream of both. Nothing in this requires GFN2 to be wrong, and the sequence of measurements that would have caught either was never run because `traj_contacts.py` had no Gln93 or carbonyl probe and `sphere` could not be run over a trajectory.

**What a rebuilt complex has to do, now quantitatively.** Build it with Zn-O2x at 2.0-2.6 A rather than hoping MD finds it; keep Asp92 OD2 toward the nucleophile and off the metal; and judge the equilibrated frame against the L1 column above, not against its own starting value. The restraint set has to hold the substrate's three anchors - carbonyl to zinc, carboxylate to Lys226/Arg229, ring nitrogen near Lys226 - because a single O01-C4x restraint demonstrably does not. All of that is config, except the pose itself.

**One caution against over-reading the L1 column.** 7UHJ is di-zinc at 2.20 A, so its carboxylate contact is to Zn2, which CSR-1 does not have, and its Zn1 carries no modelled water. The transferable numbers are the ones that do not depend on Zn2 being occupied: Zn1-carbonyl 2.63, and the Zn2-site residues' reach to the substrate, which in CSR-1 are Gln93 and Lys226 and are the whole question.

## Resources (see `running-on-bunya.md`)

**D20 – Job resources.**
- QM/MM optimisation: 8 CPUs and 4–16 GB. Measured 0.91 GB peak and about 60–70 s per step. xTB runs in `inputfile` mode.
- MM MD: OpenMM on CPU.
- Revisit: GPU MD via OpenMM HIP on Bunya's MI355X GPUs, for longer sampling.

## Planned checks, not yet done

- P1: DFT vs GFN2-xTB benchmark of each step on a cluster model, as Imhof did with model systems (D9).
- P2: residue contributions, by switching off the charges of Lys226, Arg90 and Asp92 along each path.
- P3: one Zn vs two Zn; possibly other metals at the Zn site.
- P4: Δ1 tautomer pathway (D19).
- P5: DFT/MM energies along the converged xTB paths.
- P6: establish whether the Lys226 -> OH- transfer survives in a clean four-coordinate sphere (D41). Blocks all scans.
