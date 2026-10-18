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
