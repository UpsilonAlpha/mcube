# CSR mechanism study: decision log

A running record of the modelling decisions for the QM/MM study of meropenem hydrolysis by the mononuclear (single-Zn) β-lactamase CSR. Each entry gives the decision, why it was made, the evidence, the alternatives, and whether to revisit it. Keep it up to date, so the assumptions can be stated and challenged in the write-up.

**Aim.** Test the "one metal ion is enough" idea (Dupureur, *Metallomics* 2010, 2, 609, doi:10.1039/C0MT00013B) for a mononuclear Zn β-lactamase. Use the approach of Imhof, Fischer & Smith (EcoRV, *Biochemistry* 2009, 48, 9061, doi:10.1021/bi900585m) with modern methods:
- QM/MM minimum-energy pathways for a network of candidate mechanisms
- a comparison against a small model system
- the catalytic contribution of each residue

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
- Evidence: the ring carboxylate pKa is about 3 and the pyrrolidine amine pKa about 7.5–9, so at pH 7 the molecule is mostly COO⁻ / NH₂⁺. The earlier model had a neutral COOH and a neutral amine. Metallo-β-lactamases typically bind the C3 carboxylate with a Lys, and Lys226 may be that residue (NZ is 5.4 Å from C6x in the current frame).
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
- Hypothesis: Lys226 plays a catalytic role like EcoRV's Lys92 (positioning, or stabilising the anionic intermediate), and Arg90 stabilises the oxyanion or hydroxide.
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
- **Also noted:** a numbering discrepancy for the lysine, Lys254 (PDB titles) vs Lys263 (Krco *et al.* review, *Front Chem* 2023;11:1196073, PMID 37408556), unresolved. GOB-18, not GOB-1, is the mono-zinc B3 precedent, and its single zinc sits in the **Zn2** site (PDB 5K0W); the sister B3-RQK enzyme SPR-1 is suggested to be mono-zinc at rest with a di-zinc centre forming on substrate binding (not established).
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
