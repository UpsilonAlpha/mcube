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

**D16 – Snapshots must be lined up for attack.** (2026-10-15, user decision)
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
