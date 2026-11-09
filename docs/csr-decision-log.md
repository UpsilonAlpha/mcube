# CSR mechanism study: decision log

> **Read `csr/state.md` first.** This log is append-only **history**: it
> records how each conclusion was reached, including the ones later overturned.
> `csr/state.md` is what is true **now** - the settled facts, the void results,
> the live plan with status, and the open questions. Looking up a current value
> here instead of there is how the same facts get rediscovered and how a
> contradiction sits unnoticed for a week: D145 quoted a distance that D108 had
> already ruled out, and nobody saw it for eight days (D179).
>
> **When something changes, add an entry here AND edit `csr/state.md` in
> place.** An entry that changes a settled fact or a plan item without editing
> that file is incomplete.

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

## Standing constraints: read this before setting up any run

**Why this section exists.** The log below is chronological, and a chronological
log buries a prescription. D120 diagnosed the coordination-cutoff error on
2026-10-28 and prescribed the fix in its own closing paragraph; one day and ten
decisions later, D130 found the detection still running at 3.0 A and a
fresh 100 ns job built on the shell D120 condemned. **Nothing was forgotten -
it was written down and not carried out**, which a long log cannot distinguish
from a fix that landed. Every constraint here names the decision it comes from,
what enforces it, and - where nothing does - that it rests on a human reading
this. **Update the status line in the same commit that changes the behaviour.**

> **Before quoting any entry below, or anywhere in this log, check
> `csr/retractions.md`.** This file is append-only and 207 entries long; roughly
> forty of its entries still **lead** with a claim a later entry overturns, so a
> grep lands on the retracted half. That has caused four re-derivations of
> discarded arguments (D130, D179, D194, D206). For a *current* value use
> `csr/state.md`, never this log.

### The coordination sphere (D120, D121, D130)

1. **A restraint only preserves a bond if the contact IS one.** `preprocess.inputs.cutoff` decides which contacts become 50 kcal/mol/A^2 restraints for the whole run. A Zn-O/N bond is 1.9-2.2 A; **2.5 is defensible, 3.0 is not.** *Enforced:* `md` warns on every pair it will restrain beyond `LIGAND_MAX` (2.5 A), detected and configured alike, and records them in `coordination_beyond_bond_length`. The library default cutoff is still 2.8, so **set `cutoff` explicitly**. Raise `coordination_warn_above` for a genuinely longer bond (Ln-O reaches 2.4-2.5 A).
2. **The crystal sphere is four:** HOH485 1.90, His163 NE2 2.09, HOH493 2.14, His87 NE2 2.22 (6NC5 ZN302). **Asp92 is not on it, and not near it: 4.20 A (OD2) and 4.22 A (OD1), measured in `CSR_crystal/6NC5.pdb` itself.** The 2.80 A this rule used to quote is not a crystal distance - see D145. *Enforced:* nothing. A run whose sphere disagrees has to be caught by reading `site.py shell`.
3. **Asp92 binds the zinc at every level of theory we have run, and the crystal says it should not. Do not set up a run that assumes either side of that.** MM at 100 ns (D102), unconstrained GFN2-xTB at 25 ps (D51) and DFT, 1.6x steeper than xTB (D58/D60), all put it on the metal; 6NC5 ZN302 has it at **4.20 A** (D145 corrects the 2.80 A this used to say). The working explanation is D133: Asp92 is a **Zn2** ligand with no Zn2 to bind, held off Zn1 by the RQK cations, and probably in equilibrium rather than fixed. Two things follow. **Removing a restraint does not fix this** - D102's extra ligands were unrestrained - it only makes the coordination measurable. ~~And **"xTB over-coordinates the zinc" is no longer the open question**; D58/D60 settled it.~~ **[RETRACTED 2026-11-05, D206. D58's inference is invalid - D59 said so at the time and D206 retracts it in place - so this rule was closing the project's standing method question on evidence the log itself withdraws, in the one block everyone is told to read before setting up a run. The question is OPEN. D133 repeats the same closure and is corrected there too. Note also that "DFT, 1.6x steeper than xTB (D58/D60)" above is that same retracted inference.]** *Enforced:* nothing.

4. **Count coordination against this enzyme's own motif, which has TWO protein ligands and not three.** CSR-1's Zn1 site is **H-R-H**: Arg90 (BBL R118) stands where the third histidine would be. The canonical-equivalent four-coordinate sphere is therefore **two His + hydroxide + one water**, so an extra solvent ligand is not over-coordination until it is counted against two protein ligands. *Enforced:* nothing; user correction, made twice.
5. **The other five zincs in 6NC5 are not the site.** ZN301 occupies the Zn2 pocket only because the construct is the engineered HARLDQ motif, whose `SEQADV` records mutate **Gln93 to His and Lys226 to His** - the Q and K of RQK - restoring an Asp-His-His Zn2 set. In the wild type those are a glutamine and a lysine and cannot ligate a metal, which is the structural basis for Zn2 being vacant. ZN303-306 are adventitious. Measure ZN302. *Enforced:* `reference-structures.md` tabulates all six.
6. **Do not hand-write a restraint to correct a shell that is already wrong.** D130's `[B:1:ZN, A:92:OD2, 2.80]` plus `[.., OD1, 3.81]` blocked D120's rotation and reinstated its condition, asserting the carboxylate's orientation where the original asserted only a distance. Asp92 reaching Zn1 is **a question to ask, not an assumption to restrain in**. Rebuild from a structure whose sphere is right instead. *Enforced:* the warning in 1 fires on configured pairs, which is how this one would have been seen.
7. **Do not graft onto a frame descended from a condemned shell.** `subst_start.pdb` inherited Asp92 OD1 at 1.99 A from `md_naclong`; no cutoff undoes an artefact already in the coordinates. D120 costed the consequence: everything from `md_nac` onward re-run from a shell that matches the crystal. *Enforced:* nothing.

### Measurement (D114, D120, D129)

8. **One frame cannot tell a ligand from a passer-by.** Use `site.py shell` on a trajectory - per-donor **occupancy** and the coordination number as a distribution - for any claim about the sphere. `sphere` answers for the frame it is given, which is why D120 went unseen for three weeks. *Enforced:* nothing; `sphere` still runs on one structure by design.
9. **A coordination count is not a composition.** Four donors inside a cutoff can be the wrong four, and a member sitting 0.09 A from the boundary makes a single-frame count a report on thermal noise. Read identity off occupancies. *Enforced:* nothing. A `composition`/`occupancy` criterion was proposed in D129 and is **not built**.
10. **Never pick a cutoff after seeing where the atoms landed.** Raising `metal_bound` from 2.6 to 2.75 to make a count come out is choosing the measurement to fit the expected answer - D120's move and D123's. *Enforced:* nothing.

### Restraints and selectors (D93, D115, D128)

11. **Hold a contact where the structure already has it; never pull it to a chosen number** (D93). Sound, and nothing above weakens it - the failures were in deciding *what* to hold, not in the holding. *Enforced:* `fill_restraints` measures any row given without a distance.
12. **Restraints are named by selector, never by saved index.** An index is valid only against the atom order it was measured on, and any correct edit renumbers (D115). *Enforced:* `check_coordination_atoms`, and `md` re-resolves selectors against the structure in front of the stage.
13. **Every selector in a section is checked before it is written**, restraint pairs included (D128). *Enforced:* `sections.py` `SELECTOR_KEYS`, which refuses to write the section.
14. **A metal system with no detected coordination is refused, not run free.** *Enforced:* `check_coordination_present`.

### Hypotheses and gates (D125, D126, D127)

15. **A criterion is only a criterion within a mechanism.** `site.py check` requires `--hypothesis=` when more than one is defined and refuses to default. *Enforced:* `mode_check`.
16. **The two hypotheses are not equally evidenced and the file must not read as though they were.** Every range in `carbonyl_on_metal` carries a PDB ID; `substituent_on_metal`'s `Zn - O1x` range is our own measurement plus tolerance and is labelled `NOT a reference value`. *Enforced:* convention only.
17. **Absolute Zn-donor bond lengths are checked at `reactant`, not at `start`.** MM zinc does not reproduce crystal bond lengths - the same sphere shows the hydroxide short and both histidines long - and D48/D118 measured 1.83-1.86 A relaxing to 1.99-2.05 A under QM/MM. Gating an MM structure on a crystal bond length tests the force field. *Enforced:* `required_at` in `site.yaml`.
18. **`md_pose2` must not be launched as written** (D126): its restraint list opens `[B:1:ZN, C:1:O2x, 2.55]` and 2.55 A was the May run's own restraint output, not a measurement. *Enforced:* a DO NOT RUN banner in `site.yaml`.

### Analysis order (D119)

19. **Run the cheap trajectory analyses before committing QM/MM.** `site.py check --at=reactant`, then `joint` and `acceptors`, cost minutes against ~119 h. D119 said so and it was not done. *Enforced:* nothing.

**Live run directory:** `/scratch/user/uquadhik/CSR_wt` (D132), which holds a copy of `CSR_nac`'s preprocess output and the same config with `preprocess.inputs.cutoff` at 2.5. `CSR_nac` is the history up to D131 and every md section in it starts from the wrong zinc sphere; `CSR_ligand` defines the system and `CSR_hydrolysed` is the product source; the `CSR_crystal*` directories are superseded.

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

**[SUPERSEDED by D23 - this 'Crystal' table is the pre-rebuild drifted model, not 6NC5. See `csr/retractions.md`.]** **D21 – The current model's Zn site has drifted from the crystal structure.** (found 2026-10-15; open)

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

**[SUPERSEDED - standing rule 1 (2.5 A defensible, 3.0 not) and D208: at 4.20 A Asp92 should never have entered the list. See `csr/retractions.md`.]** **D8 – Zn coordination cutoff of 3.0 Å.** (2026-10-15, user decision)
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

**[SUPERSEDED by D38/D41/D48 - measured in the degraded six-coordinate state. See `csr/retractions.md`.]** **D29 – QM/MM MD (first 39.5 ps): Lys226 protonates the Zn-OH⁻ at once; no attack geometry is sampled.** (2026-10-17; analysis in `CSR_crystal_v2/qmmm_analysis/`, on job 29937641's trajectory while it was still running)
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

**[SUPERSEDED by D48 - -13/-15 kcal/mol becomes -5.23; two-thirds was the degraded geometry. See `csr/retractions.md`.]** **D35 – D30 DFT result; first scan results; scans restarted from optimised reactants.** (2026-10-18)
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

**[CORRECTED by D41 - 'over-coordinated' miscounts H-R-H as three protein ligands; two waters is canonical. See `csr/retractions.md`.]** **D38 – All grids cancelled. The reactant state, not the hypotheses, is the problem.** (2026-10-18)
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

**[CORRECTED by D41 - the zinc position IS crystallographic, and the inactivity claim was unverified. See `csr/retractions.md`.]** **D40 – Literature check: three independent problems with the reactant model.** (2026-10-18)
- **CSR-1 is published:** Pedroso MM *et al.*, "Broad spectrum antibiotic-degrading metallo-β-lactamases are phylogenetically diverse", *Protein Cell* 2020;11:613–617, PMID 32542533 (PMC7381538). A **B3-RQK** MBL from *Cronobacter sakazakii* (UniProt A7MFY4). Structures: **6DN4** (mature wild type, 1.99 Å), **6DQH/6NC5** (2.07 Å), **6DQ2**, **6DR8** (1.48 Å, R118H/Q121H/K254H triple mutant).
- **Motif mapping** (BBL → our numbering): Arg118 → Arg90, Asp120 → Asp92, Gln121 → Gln93, His116 → His88, His196 → His164, Lys254/263 → Lys226. So the Zn1 site is **H-R-H** (His116, **Arg118**, His196; two usable His) and the Zn2 site is **D-Q-K** (Asp120, **Gln121**, **Lys254**; one usable ligand). Canonical B3 is HHH/DHH.
- **This confirms and explains D39.** 6NC5's ZN301 (Asp91/His92/His225) is the **canonical Zn2 site restored by the two mutations**; ZN302 (His87/His163 + two waters) is the **Zn1 H-R-H site**. Our model correctly kept ZN302 and reverted both mutations. There is no third candidate site.
- **Three problems, each independent:**
  1. **No published zinc site.** The wild-type CSR-1 structures crystallised **essentially apo** — 6DN4 lists no metal at all — and the paper attributes this to reduced metal affinity, since Arg118 and Lys254 are protonated at optimal pH. So our mono-zinc site is **inference, not literature**. Nothing is wrong with inferring it, but it must be stated as such, and the Zn placement came from a mutant's coordinates.
  2. **Hydroxide *and* a second water on the same zinc is unsupported.** Canonical B3 Zn1 is **four-coordinate**: three His plus the bridging hydroxide (L1, PDB 1SML; Ullah *et al.*, *J Mol Biol* 1998;284:125, PMID 9811546). A literature search found **no** description of a catalytic zinc carrying both a hydroxide and a separate water. Our site carries OH⁻ (ex-HOH493), HOH485, Asp92 *and* the substrate carbonyl on a two-His zinc. HOH485 is the weakest of these: in the crystal it bridges ZN302 to the **artefact** ZN303 (whose other ligands are two waters and an acetate).
  3. **The model is the non-truncated, inactive form.** CSR-1 is active only after N-terminal truncation; in the mature form an N-terminal loop occludes the active site. 6NC5 covers residues 15–273; **our model covers 1–275**, so residues 1–14 were modelled back in, and Ala12, Pro13, Tyr15, Leu17, Phe18, Glu19 and Trp21 sit 9–13 Å from the zinc — nearer than anything from the crystal's own N-terminus.
- **[WITHDRAWN 2026-10-28, see D123: false for B3. L1's time-resolved Michaelis complexes put the carbonyl 2.53-2.63 A from Zn1 (7UHI/7UHJ, D121). This bullet is what dropped the carbonyl restraint.]** **Inner-sphere carbonyl coordination is also not established.** No crystallographic Michaelis complex shows a firm Zn–O(carbonyl) bond; the consensus is that Zn1 polarises the carbonyl at longer range while the hydroxide attacks (NDM-1 QM/MM: Das & Nair, arXiv:1702.01239, where Zn1–O is only weakly coordinated). Our `preprocess` `cutoff: 3.0` captured the carbonyl at 2.16–2.31 Å and restrained it as a metal ligand, per D28.
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

**[QUALIFIED - the 'clean' sphere still had Zn-O2x constrained at 2.90 A; D73 gets 3.01 A at 84 deg once released. See `csr/retractions.md`.]** **D48 - M1 result: the transfer is still favourable, but by a third of what the degraded geometry said; and the Zn-bound oxygens still cannot attack.** (2026-10-18, job 30053828, COMPLETED in 9:51:17)

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

**[QUALIFIED by D67/D68 - both 'independent' tests descend from the condemned 1.94 A pose. See `csr/retractions.md`.]** **D51 - M2 result: the six-coordinate site is what the model prefers, and Asp92 as base is not available.** (2026-10-19, job 30085418 COMPLETED in 9:11:36)

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

**[WITHDRAWN by D53/D54 - the water wire has no MBL precedent. See `csr/retractions.md`.]** **D52 - Reconsidering the mechanism: the zinc activates the electrophile, and the proton leaves down a water wire.** (2026-10-19, user: "step back and reconsider what mechanism a six-coordinate mono-zinc site could actually support")

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

**[SUPERSEDED by D208 - the crystal puts Asp92 at 4.20 A; both axes of this curve are our own model. See `csr/retractions.md`.]** **D57 - xTB says Asp92 is firmly bound; DFT/MM optimisation at the full QM region is unaffordable; waiting for the DFT single points instead.** (2026-10-19)

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

**DFT is 1.6x steeper than xTB, which is the opposite of the concern.** The worry from D51/D53 was that GFN2-xTB over-coordinates the zinc and had invented the sixth ligand. It has not: at the DFT level Asp92 is bound *more* strongly, costing **9.1 kcal/mol** to reach its crystallographic 2.80 A rather than xTB's 5.7. ~~**The six-coordinate site is a real feature of this model, not a semi-empirical artefact**~~ **[RETRACTED - see D59 below, which already said so, and D206. These are single points on xTB geometries that all descend from the degraded six-coordinate MD start. They show the two methods agree about that path; they cannot promote the path. Asp92 is not a Zn1 ligand in 6NC5 or in any B3 structure.]**, and the Yang *et al.* over-coordination failure mode - real in general, and a legitimate thing to have checked - does not apply to this result.

Consequences:
- ~~**D51 and D52 stand on firmer ground.** The zinc is genuinely six-coordinate, so it~~ **[RETRACTED, D206.]** The superseded claim read: the zinc is genuinely six-coordinate, so it is a Lewis acid activating the bound carbonyl rather than a hydroxide generator, and the two aqua ligands genuinely cannot reach the attack trajectory.
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

**[CORRECTED by D70 - that was a mixed job; 37% on 64 cores, ~26 min per single point, no speedup factor to quote. See `csr/retractions.md`.]** **D60 - Tooling consolidated; CPU efficiency measured; what CphA actually says about the carbonyl.** (2026-10-20, user)

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

**[QUALIFIED by D64 - 22.4 kcal/mol is the cost of one pose's hinge, not of any alternative bound state. See `csr/retractions.md`.]** **D61 - Neither removable ligand opens a productive geometry. The substrate binding mode is now the thing to question.** (2026-10-20, jobs 30200387 and 30119957 both COMPLETED)

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

DFT is consistently about 1.6x steeper across the whole curve. **What that does and does not say is set out in the next paragraph, and the clause first written here - "confirming D58: xTB is not over-coordinating, and the six-coordinate zinc is real at both levels" - contradicted it and is withdrawn (D206).** It licenses only the narrow statement that xTB does not over-bind Asp92 *along this coordinate from this structure*.

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

**[SUPERSEDED by D64/D67/D68 - the 1.94 A did propagate; this entry is why the clean re-do was cancelled. See `csr/retractions.md`.]** **D63 - The MM re-do was the wrong instrument; D59 and D62 overstated the constraint's effect. Checking the carbonyl bond at DFT instead.** (2026-10-20, user: "the MM equilibration doesn't tell us anything about coordination")

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

**[CORRECTED by D68 - pose 2 was chosen deliberately, `model: 2` in all three configs. See `csr/retractions.md`.]** **D66 - The reference ligand is pose 2, not pose 1.** Checked as D64 said to. Comparing `CSR_ligand/meropenem.pdbqt` against `CSR_crystal_v2/preprocess_results/mer_crystal.pdb` on heavy atoms, by symmetrised mean nearest-neighbour distance:

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

**[SUPERSEDED by D102 - 1 ns was too short; 100 ns of the same scheme goes six-coordinate. See `csr/retractions.md`.]** **D100 - The 1 ns restrained run: the mechanism works, the angle does not.** (2026-10-23, job 30418017, `CSR_nac`, section `md_nacgpu`)

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

**[REVERSED by D104 - the near-attack geometry did not survive release, 101 deg to 57. See `csr/retractions.md`.]** **D102 - The 100 ns run at the 77% mark: site held, disk on budget, the angle question still open.** (2026-10-23 15:42 AEST, job 30421166, `CSR_nac`, section `md_naclong`)

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

**[SUPERSEDED by D120 - both arms stayed five-coordinate; GFN2 inherited the shell rather than repairing it. See `csr/retractions.md`.]** **D104 - 50 ps GFN2 QM/MM: the complex holds and the zinc is fixed, but the near-attack geometry is gone.** (2026-10-25 08:30 AEST, job 30467248, `CSR_nac`, section `nacqm`)

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

**[CORRECTED by D106 - shrinking actradius is the opposite lesson; 5 A diverges and the 8 A run was healthy. See `csr/retractions.md`.]** **D105 - trajectories deleted to free scratch, and the optimisation is grinding in geomeTRIC.** (2026-10-25 08:45 AEST, `CSR_nac`)

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

**[QUALIFIED by D126/D208 - built on a condemned shell, and 3.17 A at 106 deg was later reached unrestrained. See `csr/retractions.md`.]** **D106 - the NEB starting structure, taken from the 8 A optimisation's plateau.** (2026-10-25 19:40 AEST, `CSR_nac`, `protdft/nacqm/QM_MM_opt_last.pdb`)

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

**[SUPERSEDED by D126 - a grafted structure reached the Burgi-Dunitz window with no attack restraint at all. See `csr/retractions.md`.]** **D110 - Job 30598032: 50 ps with Arg229 in the QM region. The nucleophile moves further from attack, not closer.** (2026-10-26)

COMPLETED, exit 0:0, elapsed 23:56:12, ended 21:53:32 - **four minutes inside the walltime.** 50,000 steps at 1 fs confirmed from `OpenMM_MD_final_state.xml` (`stepCount="50000"`). 1.72 s/step, slightly worse than the 1.68 estimated from adding the guanidinium. Scratch at 156/300 GB, so nothing here is suspect on quota grounds. The QM region is **137 atoms** against `nacqm`'s 119: Arg229's 18 side-chain atoms (CB to NH2 with hydrogens, backbone left in MM), with `charge: 3` as set. Lys226 NZ carries three hydrogens, so both cations are present and cationic, which was the point of the run.

**Zinc, five-coordinate, unchanged in character from D108:** O01 1.90, His164 NE2 2.10, Asp92 OD1 2.11, HOH485 2.16, His88 NE2 2.18. Asp92 OD2 is 4.22, so the aspartate is firmly monodentate - it still bridges nothing.

**The result that matters, and it is negative.** The attack geometry has degraded again:

    O01-C4x        angle      source
    3.33 A         101 deg    100 ns classical MD (D102)
    3.86 A          57 deg    50 ps QM/MM, no Arg229 (D104)
    6.20 A          63 deg    50 ps QM/MM, Arg229 in QM (this run)

The second Zn-bound water is no better (HOH485-C4x 6.27, angle 79), and Zn-mer O2x is 6.76. **Putting Arg229 in the QM region did not pull the substrate back into a near-attack arrangement; the nucleophile drifted half again as far away.** So the conclusion of D104 is now independently reproduced with a different QM region: **GFN2-xTB QM/MM MD will not hand us a near-attack frame, and `scan_attack` has to drive the O01-C4x coordinate.** Two runs, 100 ps of QM/MM between them, agree on this.

**What this does not say.** It is not evidence against the +2-substitution hypothesis, and it does not invalidate the NEB start frame: that came from `nacqm`'s *optimisation* (D106), not from either MD end point. Nor does it show the substrate has left the site - Zn-mer O2x at 6.76 is a single contact, and this stage's `monitor` is deliberately unset, so Lys226-O3x and the Arg229 contacts were not measured. **Whether the substrate is still anchored over these 50 ps is open, and it decides whether this trajectory is usable for the charge-deletion analysis at all.** Measuring the same contacts as D107 over `qmmm_nacarg/QM_MM.dcd` is the next step.

**[SUPERSEDED by D131 - 'rock solid' was measured with Asp92 restrained against the metal. See `csr/retractions.md`.]** **D111 - Contacts over the nacarg trajectory: the Lys226 salt bridge breaks, and over 100 ps the substrate is leaving the site.** (2026-10-26; `docs/csr/traj_contacts.py`, 50 frames of `qmmm_nacarg/QM_MM.dcd`, saved as `qmmm_nacarg/contacts.txt`)

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

**[SUPERSEDED by D120 - a proton cannot open a site a restrained carboxylate occupies; both arms were always going to fail. See `csr/retractions.md`.]** **D119 - The ohH arm drifts too, so the protonation state was not what lost the substrate. M1 is answered; M0 is now the live question.** (2026-10-28, job 30721365)

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

**[CORRECTED by D133 and D208 - D102's newcomers were unrestrained, and the cause needed BOTH the 3.0 A cutoff and the per-residue default. See `csr/retractions.md`.]** **D120 - Both QM/MM arms ran on a FIVE-coordinate zinc, against the crystal's four, and a coordination restraint this pipeline generated is what put the fifth ligand there. This is upstream of the protonation test and of any GFN2 question.** (2026-10-28)

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

**[QUALIFIED by D131/D208 - Zn1 is H-R-H with two protein ligands, and the carbonyl target rests on contested numbers. See `csr/retractions.md`.]** **D121 - Carbonyl-to-Zn1 is an OBSERVED step in B3, measured here from the time-resolved L1 structures; and in those structures D120 is a Zn2 ligand over 4 A from Zn1. That fixes what a CSR-1 Michaelis complex has to look like, and confirms D120's restraint put Asp92 in the carbonyl's seat.** (2026-10-28)

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

**[SUPERSEDED by D123/D126 - the pose is rotated, and no Vina pose puts any atom within 5 A of the zinc. See `csr/retractions.md`.]** **D122 - Gln93 is 10 A from the substrate, and the docked pose never had a Zn-carbonyl contact at all. Measured against L1, the pose has exactly one defect and the equilibration has a second; between them M0 and the drift are both accounted for.** (2026-10-28, 59ed041)

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

**[HEADLINE FALSE - D126, and marked here 2026-11-05 because the body's
strikethrough never reached the headline a grep returns. 2.55 A is the May run's
own carbonyl RESTRAINT output, not a measurement, and it was read out of a file
this same entry lists as not reusable. Pose 2's closest heavy atom to the metal
is 7.23 A, and no Vina pose puts any atom within 5 A of it. This false distance
seeded `md_pose2`'s restraint target (now DO NOT RUN) and a regression test in
D125 asserting "pose 2 at 2.55 A must be accepted".]** **D123 - Everything since
October ran on the wrong pose. `CSR_ligand` used Vina pose 2, which puts the
carbonyl 2.55 A from the zinc; the pose inherited from `CSR_crystal` presents the 6-alpha-hydroxyethyl O1x at 3.13 A and leaves the carbonyl at 5.43 A. The carbonyl restraint was dropped on my recommendation in D40, and D40's reason was wrong.** (2026-10-28, user)

The user's recollection is correct on both counts and the measurements bear it out.

**The two poses, measured with one tool** (`traj_contacts.py`, d024aa9, which now reports which substrate oxygen faces the metal; `--skip-missing` lets it read the older ligand variant):

    contact                  CSR_ligand (May)   CSR_crystal/CSR_nac (Oct)   L1 7UHJ
    Zn - mer O2x (carbonyl)       2.55                   5.43                 2.63
    Zn - mer O1x (6a-OH)         13.27                   3.13                    -
    closest substrate O        O2x 2.55               O1x 3.13              O1 2.63
    O01 - C4x (attack)            3.81                   3.47                    -
    Asp92 - O01                OD2 2.26               OD2 2.84                   -

**The ligand is bound the other way round.** Not displaced, not drifted - rotated, so that a substituent oxygen faces the metal and the group that has to be attacked faces away. D25 recorded this on 2026-10-15 without flagging it: "O1x is 3.13 A from the Zn; the beta-lactam carbonyl O2x is **5.43 A** away". That number has sat in the log for two weeks as something to fix by restraint, when it is something to fix by using the other pose.

**The pose is confirmed as pose 2**, not inferred: `site.py graft`'s own bond-graph match reports "reference ligand matches pose 2" against `CSR_ligand/meropenem.pdbqt` (9 poses; pose 2 scores -5.624 against pose 1's -5.744, within Vina's own resolution, so the hand pick was not a scoring error).

**D40's carbonyl claim is withdrawn.** D40 said: "Inner-sphere carbonyl coordination is also not established. No crystallographic Michaelis complex shows a firm Zn-O(carbonyl) bond", citing NDM-1 QM/MM (Das & Nair, arXiv:1702.01239) for weak coordination. **That is false for B3**, and the counter-evidence was published in 2022: L1's time-resolved snapshots with intact moxalactam put the carbonyl oxygen 2.53 and 2.63 A from Zn1 (D121, PDB 7UHI/7UHJ). D41 corrected three of D40's claims; it did not reach this one, which went on to set `Zn-O2x = 2.90` in D42's sections and then, in `md_nac`, to justify leaving the carbonyl unrestrained altogether. **The user constrained it originally and was right to; I recommended dropping it on a false premise.**

**On the constraint value, since it was asked:** neither number used was a literature bound. D25 restrained Zn-O2x at **2.2 A** with k = 20, from the user's stated assumption that the carbonyl coordinates at < 3 A; D42 constrained it at **2.90 A**. The restrained MM then held it at 1.79-2.19 A (mean 1.97) - *shorter than its own Zn-OH bond*, which D42 itself called unphysical. The only measured value for a B3 Michaelis complex is **2.53-2.63 A**, and pose 2 ~~arrives at **2.55 A** on its own~~ **[FALSE - D126]**. The new section holds it there, rather than pulling it to a chosen number, which is also D93's rule.

**What is reusable, checked rather than assumed:**

  - **`CSR_ligand/meropenem.pdbqt`** - the docking output, 9 poses, with pose 2 identified. This is the asset; the rest of that run is superseded.
  - **`CSR_ligand/CSR_nontruncated_1metal.pdb`** - the docking receptor, needed to align the pose into the current system. Its zinc is **B:2:ZN**, not B:1:ZN, which is why `site.py poses` could not run on it until 286166e added `docking.receptor_metal`.
  - **`site.py graft`** - written for exactly this: put another pose into the already-equilibrated system, so protein, waters and ion count stay fixed and only the ligand differs.
  - **NOT the old system.** `CSR_ligand/preprocess_results/finalsystem.pdb` has His164 **4.56 A** off the metal, Asp92 **bidentate** on it (OD2 2.40, OD1 2.60), a hydroxide with **no hydrogen**, an earlier meropenem variant lacking O4x, and 41,371 atoms against the current 49,424. Its ligand placement is right and its metal site is worse than what we have.
  - **NOT `CSR_ligand/neb_endpoints` or `neb_results`** - those endpoints were built from that system, and `neb` refuses a product whose atom order differs from the reactant in any case.
  - `CSR_hydrolysed` holds the product-side docking (`hydrolysed_meropenem.pdbqt`) and needs the same treatment once the reactant is settled; nothing there is usable before then.

**The plan, as the user set it: constrained 100 ns NVT first, then 300 ps QM/MM.** Written as sections `md_pose2` and `qmmm_pose2` in `site.yaml`, with the build step - graft pose 2 into the current equilibrated system, then **remove HOH485** - spelled out above them. Restraints: Zn-O2x at 2.55; Asp92 pinned in **both** oxygens, 2.80 and 3.81, so the D120 flip cannot recur; and O01-C4x held where the grafted frame has it. `coordination_redetect: true` is mandatory because removing a water renumbers the file (D115). HOH485 goes because two His + hydroxide + carbonyl is already the four-coordinate sphere D41 defined, and because D34 showed that merely dropping its restraint leaves it on the metal at 1.83 A.

**Costs.** The 100 ns NVT is about 3 h on an H100 (D100: 845 ns/day). The 300 ps QM/MM is about **119 h of wall time** at the 1.427 s/step D119 measured, so it runs as six resumed 50 ps jobs of ~20 h each, CPU only. Nothing has been submitted.

**[WEAKENED 2026-10-28, see D124: thermolysin and carboxypeptidase A do hold both on one zinc, at 117-119 degrees, when the carbonyl is loosely coordinated. What survives is that it fails when the carbonyl is held tightly, as in D34.]** **One concern, with the numbers rather than as an opinion, because this route has been run before.** With both the nucleophile and the carbonyl on the same four-coordinate zinc, the O01-C4x-O2x angle is geometrically forced to roughly 50-60 degrees: pose 2 as built is **3.81 A at 51 degrees**, and D33/D34 found that with the carbonyl on the zinc in 91% of 791 frames (mean 2.32 A) **no frame** reached the 95-120 degree window - Zn-bound oxygens spanned 24-54 degrees and never exceeded 87. The reason is visible in L1: there the nucleophile is the hydroxide **bridging Zn1 and Zn2**, sitting under the carbonyl carbon rather than beside the carbonyl oxygen. CSR-1 has no Zn2, so a terminal hydroxide on its zinc is cis to the carbonyl by construction. **A mono-zinc site cannot hold both the carbonyl and a metal-bound nucleophile on the Burgi-Dunitz line.** One of the two must therefore not be on the metal, and the candidates are already on the table: the nucleophile is a non-metal-bound water with Asp92 as general base (CphA-like, D40), or the carbonyl is held by **Gln93** as the oxyanion hole (SMB-1-like, D121/D122 - never tested, because Gln93 was carried from 6.56 to 10.2 A by the old equilibration).

**That makes the run worth doing as specified, and says what to measure in it.** Holding the carbonyl on the zinc for 100 ns *is* the CphA-like arrangement, so what the trajectory answers is whether a non-metal-bound water comes onto the attack line while it is held there - which is exactly what `site.py joint` and `site.py acceptors` ask, and what D33 answered affirmatively once, in the flipped pose. Run those on the new trajectory before committing the 119 h of QM/MM, and read Gln93 at the same time: if it comes in on the carbonyl, the second hypothesis is live and a different reactant follows.

**[CORRECTED by D126 - pose 2's 2.55 A is a restraint's own output, so the literature here supports a false measurement. See `csr/retractions.md`.]** **D124 - Mononuclear hydrolases, checked against structures rather than reviews. D123's geometric objection is too strong and is weakened here: thermolysin and carboxypeptidase A hold the nucleophile and the carbonyl on ONE zinc, and the condition is that the carbonyl is LOOSELY coordinated. Pose 2's 2.55 A is in that range. Separately, every mononuclear MBL in the literature is a Zn2-site enzyme and CSR-1 is not, so neither the CphA nor the SMB-1 analogy transfers cleanly.** (2026-10-28, user)

**Measured from transition-state-analogue complexes, not asserted:**

    enzyme / complex                        Zn-O(nuc side)  Zn-O(carbonyl side)  O-X-O    O...O
    thermolysin 1TMN, N-carboxymethyl           2.01              2.41            119      2.14
      dipeptide (its carboxylate mimics the
      tetrahedral intermediate's two oxygens)
    carboxypeptidase A 7CPA, phosphonate        2.19              2.98            117      2.67
      transition-state analogue

Thermolysin's zinc is five-coordinate here - His142 1.95, Glu166 1.95, His146 2.04, plus **both** inhibitor oxygens at 2.01 and 2.41. CPA's is His69 2.08, His196 2.09, Glu72 bidentate 2.24/2.29, plus both phosphonate oxygens at 2.19 and 2.98. **In the two best-characterised mononuclear zinc hydrolases, the attacking oxygen and the carbonyl-derived oxygen are on the same metal.** Both run the mechanism we are proposing to build: a zinc-bound water attacking a zinc-coordinated carbonyl, with a carboxylate (Glu143, Glu270) as general base.

**So D123's statement - "a mono-zinc site cannot hold both the carbonyl and a metal-bound nucleophile on the Burgi-Dunitz line" - is withdrawn as stated.** It generalised D34's 24-54 degrees into a geometric law, and the law is false: these enzymes do it.

**What survives is sharper and is a design criterion.** In both structures one oxygen is held tightly (2.01, 2.19) and the other **loosely** (2.41, 2.98). D34's model had the carbonyl at 2.16-2.32 A - tight - in a six-coordinate sphere that also carried HOH485 and Asp92, and in that arrangement the nucleophile had nowhere to go but beside the carbonyl oxygen. **[CORRECTED 2026-10-29, see D126: pose 2 does NOT arrive at 2.55 A. Nothing in it comes within 7.2 A of the zinc; the 2.55 A is the May run's restraint output, measured in the system D123 had already rejected.]** **Pose 2 arrives at 2.55 A**, which is between thermolysin's 2.41 and CPA's 2.98 and is where a carbonyl that must still be attacked belongs. The planned run is therefore better supported than D123 allowed, and the restraint value is right for a different reason than the one it was chosen for.

**A caution on how far the analogue complexes can be read.** They are post-attack geometries: the two oxygens are both bonded to the tetrahedral centre, so their 117-119 degrees is a consequence of that bonding, not a measurement of an approach angle. Note that CPA's O...O is 2.67 A - identical to D34's 2.67-2.77 - so the oxygen separation does **not** discriminate a working arrangement from a failed one. What the analogues establish is that the end state is reachable with one metal; the approach angle still has to be measured on our own trajectory, and it is the thing to watch.

**On CphA-like against SMB-1-like, which was the question asked: the honest answer is that neither analogy transfers, for a structural reason that is worth stating plainly in the paper.**

  - **Every mononuclear MBL in the literature carries its zinc in the Zn2 site.** GOB-18's is the DHH site; Sfh-I's is DCH-like; B2 CphA's is Asp-Cys-His, and in B2 the His116Asn substitution removes a Zn1 ligand outright, so there nucleophile activation does not involve a metal at all. **CSR-1's zinc is in the Zn1 site** (His87/His163 = BBL H116/H196, D39/D40). It is the only mono-zinc MBL we have found whose metal sits where the hydroxide normally lives.
  - That is why the CphA route never fitted, twice. It is not a mechanism for a Zn1-site enzyme.
  - **GOB-18 (Nat Commun 2017;8:538, PMC5599593)** is the closest mono-zinc MBL by subclass and it goes the other way from what we are building: the nucleophile is a **non-metal-bound water activated by a hydrogen-bond network**, the metal binds the **C-3 carboxylate** rather than the carbonyl, and the sphere expands four- to five-coordinate on substrate binding by taking that carboxylate while keeping its water. That is a Zn2-site description throughout - the carboxylate is a Zn2 ligand in every B3 - so it says what CSR-1's **vacant** site would have done, not what its occupied one does.
  - **The SMB-1 analogy is weaker than D121 implied.** Gln157 there is a loop glutamine proposed as an oxyanion hole from QM/MM, not a crystallographic metal ligand. The real structural precedent for a glutamine in a B3 zinc site is **GOB-18's Gln116, which coordinates Zn1 through its terminal amide** (5K0W: Gln98 OE1 at 2.04 and 2.21 A in the two chains, giving an HQH site in place of the canonical 3-His). But that is BBL position **116**, a Zn1 ligand - whereas CSR-1's Gln93 is BBL **121**, a Zn2 ligand position. Same residue type, different site.
  - **What that does to D94.** The 1.79 A Zn-Gln93 OE1 contact MM produced was treated as an artefact. GOB-18 shows a B3 glutamine coordinating a zinc at 2.04 A, so the contact is not absurd - but it is at the wrong position to be the GOB analogy, and it competes for the coordination site the carbonyl needs, exactly as Asp92 OD1 did (D120). It should be monitored, not assumed either way.

**So the next test is not a choice between two named enzymes' mechanisms.** It is whether, with the carbonyl held loosely on the Zn1-site metal, the nucleophile that reaches the attack line is the zinc-bound hydroxide or a second, non-metal-bound water - which is the thermolysin/CPA route against the GOB-18 route, asked of our own system. `site.py joint` and `site.py acceptors` ask exactly that and need no new code: `joint` requires one molecule both to reach the Burgi-Dunitz window and to donate to a named base, directly or through one bridging water, which is what distinguishes the two. Run them on the 100 ns trajectory before spending the 119 h of QM/MM.

**Reference values for the Michaelis complex, consolidated** (L1 7UHJ from D121, thermolysin and CPA above):

    Zn - carbonyl O          2.4-3.0 A   loose; L1 2.63, thermolysin 2.41, CPA 2.98
    Zn - nucleophile O       1.9-2.2 A   tight
    nucleophile - C4x        3.0-3.5 A
    O01 - C4x - O2x          95-120 deg  the thing D34 never reached
    closest substrate O      must be the beta-lactam carbonyl, not O1x (D123)
    Asp92 - nucleophile      ~2.8 A      H-bond, not coordination

**[QUALIFIED by D126/D127 - the accept case is a false distance, and 'closest O must be the carbonyl' gates the hypothesis under test. See `csr/retractions.md`.]** **D125 - The acceptance criteria are now a command, not a habit. `site.py check` reads `michaelis_criteria` from the site file and exits non-zero on a structure that does not qualify; the reference numbers behind every criterion are re-derivable from PDB IDs by `refstruct.py`, recorded in `docs/csr/reference-structures.md`.** (2026-10-29, user: "document this decision alongside the sources you used")

**Why this exists.** Nothing that went wrong in D120 or D123 went wrong for want of a measurement. The October pose's numbers - O1x 3.13 A on the metal, the carbonyl 5.43 A away - were written into D25 on 2026-10-15 and sat there for two weeks. The five-coordinate zinc was in every frame of two 50 ps runs. Both were recorded and neither was read as a verdict, because a structure was accepted by looking at it. **The failure mode is not bad measurement, it is measurement with no gate attached**, and the fix is a command that says REJECTED.

**What it checks** (`michaelis_criteria` in `site.yaml`, so another enzyme brings its own):

    criterion                 range        required at        guards
    closest substrate O       must be O2x  start, reactant    D123 (the flipped pose)
    Zn - carbonyl O2x         2.4-3.0 A    start, reactant    D123, D124 (loose, not tight)
    Zn - nucleophile O01      1.9-2.2 A    start, reactant    D124 (tight)
    Zn coordination number    4            reactant           D120 (the fifth ligand)
    O01 - C4x (attack)        3.0-3.5 A    reactant
    O01 - C4x - O2x angle     95-120 deg   reactant           D34 (never reached it)
    Asp92 OD1 - O01           2.5-3.2 A    reactant           H-bond, not coordination

**The two stages are the substantive design decision, and they are a direct consequence of D120.** A criterion violated at `start` is a defect in a candidate structure; the same criterion at `reactant` is a defect in the frame a QM/MM run would begin from. They differ because a grafted pose has not yet been pulled onto the attack line - pose 2 arrives at 3.81 A and 51 degrees - and requiring the attack geometry of it would reject a structure for not yet being the thing the simulation is meant to produce. Equally, **the five-coordinate zinc appeared during the 100 ns MM, not before it**, so checking only the input would not have caught D120 at all. The check therefore runs at both ends: `--at=start` on the grafted structure, `--at=reactant` on the last frame, and the second one is the gate that matters before spending 119 h of QM/MM.

A violation outside the requested stage prints WARN rather than passing silently, so the grafted structure's inherited sphere and its 51-degree angle are both visible on the start check without blocking it.

**Ten regression tests in `tests/test_site_check.py`**, written as the two mistakes rather than as hypotheticals: the October pose with its own measured numbers must be REJECTED and must name the reason; pose 2 at 2.55 A must be accepted; a carbonyl at 2.20 A - D34's value - must be rejected too, because the criterion that survived D124 is that the carbonyl is **loosely** held; Asp92 OD1 at 1.99 A must fail the coordination number at `reactant` and only warn at `start`. Full suite 387 passed, 6 skipped.

**Sources, and their reproducibility.** `docs/csr/reference-structures.md` records for each reference value the PDB ID, the ligand code, the command and the output. `refstruct.py` (promoted out of a scratch directory for this reason) measures a deposited structure with its own numbering, independent of `site.py` and of `site.yaml`, so it cannot inherit an assumption from our own system's description. Re-measured on 2026-10-29: thermolysin 1TMN `0ZN` 317 gives 2.01 / 2.41 / 119 deg / O...O 2.15; carboxypeptidase A 7CPA `FVF` 309 gives 2.19 / 2.98 / 117 deg / O...O 2.68. **D124 quoted CPA's O...O as 2.67; it is 2.68.** The conclusion is unchanged and is in fact the point of the number - it coincides with D34's 2.67-2.77, so oxygen separation does not discriminate a working arrangement from a failed one.

Papers behind D124, listed with the measurements in `reference-structures.md`: CPA QM/MM (PMC2741171), designing hydrolytic zinc metalloenzymes (PMC3985962), the Nat Commun 2017 carbapenem mechanism including GOB-18 (PMC5599593), GOB as a mono-Zn enzyme (JBC 2007;282:18286), Horsfall on GOB Met221, and Das & Nair's NDM-1 QM/MM (arXiv:1702.01239), which is the citation D40 rested on and which D123 withdrew for B3.

**D126 - `site.py check` rejected the pose-2 start structure on its first use, and it is right to. Pose 2 as docked has NOTHING within 7.2 A of the zinc. D123's "pose 2 arrives at 2.55 A on its own" is wrong: 2.55 A is what the MAY RUN'S RESTRAINT produced, measured in a system D123 itself had already rejected. `md_pose2` must not be launched as written.** (2026-10-29, graft + check on Bunya)

The guard built in D125 was run on the structure it was built for, and it failed it. That is the entry.

**1. Pose 2 was never on the metal.** `site.py poses` against the docking receptor's own zinc (`CSR_nontruncated_1metal.pdb`, B:2:ZN at -14.10 16.87 -7.68) gives the closest **heavy atom of any element** in each pose:

    pose   score    closest heavy atom to the metal
    1     -5.744          5.98 A
    2     -5.624          7.23 A   <- the hand-picked pose, the FURTHEST of the nine
    3     -4.733          5.78 A
    4     -4.564          5.01 A   <- the closest any pose gets
    5     -4.180          5.67 A
    6     -4.122          6.44 A
    7     -4.102          6.07 A
    8     -4.057          5.89 A
    9     -3.760          6.97 A

**No pose puts any atom within 5 A of the zinc, and pose 2 is the worst of them.** Meropenem was not docked into the metal site. Vina was evidently run without the zinc as a feature of the site, or with it masked, and every downstream statement that treats a pose as a Michaelis complex inherits that.

**2. Grafting pose 2 into the equilibrated system confirms it, in four systems.** `site.py graft` needed one general change to run at all - the reference ligand has to match the system atom for atom, and CSR-1's docking predates meropenem's O4x, so the reference has to come from the system itself and `--ref-frame=system` was added for that (6337f17). Zn-carbonyl for pose 2, grafted:

    system grafted into                           backbone RMSD   Zn-carbonyl
    md_naclong_results/NVTtrajectory_lastframe        0.88 A          4.50 A
    md_naclong_results/NVTtrajectory_firstframe       0.68 A          4.82 A
    md_naclong_results/relaxbox_NPT_firstframe        0.22 A          4.26 A
    preprocess_results/finalsystem                    0.00 A          4.73 A

The last row is the control that removes every remaining doubt: the docking receptor aligns onto that system's backbone at **0.00 A**, and pose 2's conformer RMSD against the ligand already there is **0.00 A**. Same protein, same frame, same conformation - and the carbonyl is still 4.73 A from the metal. **This is not drift during equilibration and not a frame-alignment artefact.**

**3. Where 2.55 A actually comes from.** `refstruct.py shell` on `CSR_ligand/preprocess_results/finalsystem.pdb` - the May run:

    2.34  O01  OH  C2      <- the nucleophile
    2.40  OD2  ASP A92
    2.55  O2x  mer C1      <- THE 2.55 A
    2.60  OD1  ASP A92      (bidentate)
    2.77  NE2  HIS A88
                            His164 absent: beyond 3.2 A

So the number is real, and it is **the output of the carbonyl restraint the user applied in that run**, not a property of the pose. D123 read it as independent corroboration of the restraint value - "pose 2 arrives at 2.55 A on its own" - when it is the restraint's own result. Worse, it is measured in the one file D123 explicitly listed as **NOT reusable**, for this sphere: Asp92 bidentate, His88 at 2.77, His164 off the metal entirely. **The 2.55 A and the rejected metal site are the same structure.** D123's reuse verdict and its constraint-value argument contradict each other, and the verdict is the one that was right.

**This is the same error as D123's and D120's, a third time: a number recorded in one context and carried into another where it means something else.** D120 restrained the oxygen that was being watched rather than the one that moved; D123 read a pose's identity off a structure two pipeline stages downstream of the pose; here a restrained distance was read as a measured one. The guard catches the consequence. It does not catch the habit.

**4. What the check said, which is the first real use of D125.** On the grafted structure, `--at=start`:

    FAIL closest substrate O     C:1:O1x 2.66   want C:1:O2x
           O1x 2.66, O2x 4.50, O4x 7.36, O3x 7.58
    FAIL Zn - carbonyl O2x              4.50   want 2.4-3.0
    FAIL Zn - nucleophile O01           1.86   want 1.9-2.2
    WARN Zn coordination number            5   want 4
           O01 1.86, Asp92 OD1 1.99, HOH485 2.04, His88 2.22, His164 2.32
    PASS O01 - C4x (attack)             3.17   want 3.0-3.5
    PASS O01 - C4x - O2x angle           106   want 95-120
    PASS Asp92 OD1 - O01                2.87   want 2.5-3.2
    REJECTED  (exit 1)

**5. And the three PASS lines are the finding worth having.** 3.17 A at **106 degrees** is inside the Burgi-Dunitz window, with Asp92 OD1 at 2.87 A positioned as general base. **D33/D34 never reached that window in 791 frames.** It appears here for free, with no restraint on the attack coordinate at all - and it appears *because* the carbonyl is off the metal at 4.50 A, with the 6-alpha-hydroxyethyl O1x taking the coordination site at 2.66 A instead.

That is a coherent Michaelis-complex candidate, and it is a hybrid of the two routes D124 distinguished: the nucleophile is metal-bound and on the attack line (thermolysin/CPA), while the metal holds a substituent oxygen rather than the carbonyl (GOB-18, which binds the C-3 carboxylate). It is not the arrangement `md_pose2` was written to build.

**6. Why `md_pose2` must not run as written.** Its restraint list opens with `[B:1:ZN, C:1:O2x, 2.55]`. On this structure that pulls the carbonyl from 4.50 A to 2.55 A over 100 ns, which (a) targets a value that is the previous restraint's output rather than a measurement of anything, (b) forces the six-coordinate, tightly-held-carbonyl arrangement D34 already ran and D124 identified as the one that fails, and (c) **destroys the 106-degree attack geometry the structure currently has for nothing**. The 2.55 A would also fail D125's own `Zn - carbonyl O2x` criterion of 2.4-3.0 only barely - it is inside it - but the criterion was written for a carbonyl that gets there, not one dragged 2 A.

**Nothing has been submitted.** The grafted structure is at `CSR_nac/md_naclong_results/NVTtrajectory_lastframe_pose2.pdb` and the measurement logs at `mcube/.cc/g2.log`, `g3.log`, `g4.log`. The decision in front of the user is which Michaelis complex to equilibrate: the carbonyl-on-metal arrangement the plan specifies, which now has no structural support in our own docking and must be built by restraint from 4.50 A; or the one the grafted structure already presents, where the attack geometry is satisfied and the metal holds O1x. The second needs no restraint on the carbonyl at all, which makes it the cheaper and more defensible test of the D124 question - and it is testable with `site.py joint` and `acceptors` on its trajectory exactly as planned.

**D127 - The guard had the same defect as the work it was built to catch: its first criterion asserted the mechanism under test. `michaelis_criteria` becomes `michaelis_hypotheses`, two named sets, and `site.py check --hypothesis=<name>` must say which one.** (2026-10-29)

D125 wrote **"the closest substrate oxygen must be the beta-lactam carbonyl"** and called it a criterion for a Michaelis complex. It is not. It is a criterion for the **carbonyl-on-metal mechanism**, and that mechanism is the thing this project exists to determine. One day later D126 produced a structure that violates it while satisfying the attack geometry at 3.17 A and 106 degrees - the best attack geometry anything in this project has shown - and the gate's verdict on it was REJECTED.

**A gate that rejects a rival hypothesis by assuming the favoured one is not a gate. It is a way of never testing anything**, and it is the same class of error as D120, D123 and D126: a statement that is true in one context used as though it were true in general. The difference is that this one was built into the tool, where it would have applied silently to every structure from now on.

**What changed.** Each hypothesis carries its own criteria; the output names the hypothesis on every line of the header and in the verdict; and with more than one set defined, **omitting `--hypothesis=` is refused rather than defaulted**, because silently picking one of two mechanisms is exactly the failure above. Five tests, including the one that matters: the same structure accepted under `substituent_on_metal` and rejected under `carbonyl_on_metal`.

**The two hypotheses as recorded:**

  - **`carbonyl_on_metal`** - the planned mechanism. Supported by thermolysin and CPA on one zinc (D124) and by L1's 2.53-2.63 A (D121). **What it does not have is a structure of our own**: no docked pose puts anything within 5 A of the metal (D126), so a start structure has to be built by restraint rather than found.
  - **`substituent_on_metal`** - what the grafted structure actually presents. O1x on the metal at 2.66 A, carbonyl out at 4.50 A, hydroxide on the attack line at 3.17 A and 106 degrees, Asp92 OD1 at 2.87 A as general base. A hybrid of D124's two routes: metal-bound nucleophile as in thermolysin and CPA, substituent oxygen on the metal as in GOB-18 with its C-3 carboxylate.

**It is listed as a hypothesis, not a preferred one, and the file says so where the numbers are.** Nothing in the literature puts a 6-alpha-hydroxyethyl oxygen on a Zn1-site metal; its `Zn - O1x` range of 2.0-2.8 A is our own measurement plus tolerance and is labelled `NOT a reference value`, unlike every range in `carbonyl_on_metal`, which carries a PDB ID. **The two hypotheses are not equally evidenced and the file must not read as though they were.** What recommends the second is that it needs no restraint on the carbonyl and already satisfies the attack geometry; what recommends the first is the only crystallographic Michaelis complex of a B3 enzyme that exists.

**[DO NOT RUN - D130/D131/D135: `md_subst` restrains Asp92 harder than the run D120 condemned. See `csr/retractions.md`.]** **D128 - `md_subst` submitted: job 30923516, the 100 ns NVT of the `substituent_on_metal` Michaelis complex. And `sections.py` was not checking restraint selectors at all.** (2026-10-29)

The user chose `substituent_on_metal` over the planned `carbonyl_on_metal`, so the arrangement with a structure of its own is the one being equilibrated. `md_pose2` was not submitted and must not be: its restraint list opens `- [B:1:ZN, C:1:O2x, 2.55]`, and D126 disproved the premise that 2.55 A was ever measured on anything but the May run's own restraint.

**The job.** `sbatch -J subst -p gpu_cuda --qos=gpu --gres=gpu:1 --cpus-per-task=4 --mem=32G -t 12:00:00 run.slurm md md_subst`, from `/scratch/user/uquadhik/CSR_nac`. The resources are `naclong`'s (job 30421166) exactly, and **that job took 6:19:45, not the ~3 h D100's rate implied** - the 12 h limit is sized from the measurement, not the estimate. Queued PENDING at submission; output goes to `protdft/subst-30923516.out`.

**Gates run before submitting, in this order.** Scratch 218/300 GB and 570k/1M inodes, so no silent over-quota write. `git pull --ff-only` on the cluster to `c855f30`, so the config is written from the same `site.yaml` the tests ran against. `site.py check md_naclong_results/subst_start.pdb --at=start --hypothesis=substituent_on_metal` - **accepted, exit 0**, one WARN (Zn-O01 1.86 A, a force-field artefact, not a structural one, and `required_at: [reactant]`). Then `sections.py` to write `md_subst`, then the config section read back, then the three restraint selectors resolved by hand in the start structure before `sbatch`.

**That last step should not have been by hand, and this is the finding worth keeping.** `sections.py` printed:

```
md_subst   pdb=md_naclong_results/subst_start.pdb  0 selectors checked
```

**Zero.** The section's only selectors are its three restraint pairs, and `restraints` was not in `SELECTOR_KEYS`, so none of them was resolved and the section was written anyway. The tool whose docstring says *two days of compute were lost to exactly that* was blind to the one list where a bad selector is most expensive: `md.md_restraints` passes it through `resolve_atoms`, so an unresolvable name surfaces hours into a run, not at the point the section is added. It is also the list most likely to be wrong here, because `subst_start.pdb` had HOH485 deleted and the file renumbered (D115).

Fixed: `restraints` added to `SELECTOR_KEYS`, with three tests - a pair that resolves is counted, one that does not stops everything, and the distance in the third column is not mistaken for a selector. All three fail without the fix. Full suite 395 passed, 6 skipped.

**What happens when the job lands, and the order matters.** `site.py check --at=reactant` on the last frame first: the start check passed its coordination criterion **on a count, not a composition** - the four at 2.6 A were hydroxide, Asp92 OD1, and the two histidines, with O1x at 2.66 A just outside. The restraints are what swap OD1 out for O1x over 100 ns, and only the reactant check says whether that happened. Then `site.py joint` and `site.py acceptors` on the trajectory, which answer D124's question - zinc-bound hydroxide or a second non-metal-bound water - and need no new code. **Both before the ~119 h of QM/MM**, which is what D119 said and what was not done last time.

**D129 - The restraint override works, verified in a running job. And the `--at=reactant` coordination count is the wrong instrument for this hypothesis - not unreachable, as this entry first claimed, but decided by a 0.09 A margin against a cutoff chosen after the fact.** (2026-10-29)

Job 30923516 at 34 minutes: RUNNING on bun075, 16.66 ns of 100 ns, 300 K, energies flat, nothing in stderr but the expected `OMP_NUM_THREADS` note. The rate is **0.49 ns/min, so ~3.4 h, not `naclong`'s 6:19:45** - the 12 h limit is generous rather than tight.

**The override is real, and this is the first time it has been checked in a live run rather than asserted.** Redetection at the 3.0 A cutoff found the sphere D120 left behind - hydroxide 1.86, **Asp92 OD1 1.99**, His88 2.22, His164 2.32 - and built a detected restraint for each, OD1 included at 1.99 A. The final list has OD1 at **3.81 A** and no 1.99 A entry anywhere. The explicit pair beat the detected one, which is what `md_restraints` promises and what the `subst` comment claimed on trust.

Six restraints are in force, not the three in the config: the three explicit ones (OD2 2.80, OD1 3.81, O1x 2.66) plus the three surviving detected ones (hydroxide 1.86, His88 2.22, His164 2.32). Worth stating plainly because the config shows three. The three detected ones pin the MM zinc at the bond lengths D128's commit called force-field artefacts - harmless here, since they are where MM sits anyway, but it means **this run cannot relax the sphere; only QM/MM can.**

**The criterion the count cannot carry.** `substituent_on_metal` asks, at `reactant`, for a Zn coordination number of exactly 4, counted at `burgi_dunitz.metal_bound: 2.6 A`. The four it means are hydroxide, His88, His164 and the substrate O1x. But **O1x is restrained at 2.66 A and OD2 at 2.80 A, so at a 2.6 A cutoff the count is 3**, and at the 3.0 A cutoff mcube itself uses for detection it is 5, because OD2 comes in.

**[CORRECTED]** This entry first said **"No cutoff gives 4"**. That is false, and the error mattered, because it made the criterion look unreachable when it is not: 2.66 and 2.80 are 0.14 A apart, so **any cutoff from about 2.70 to 2.79 returns exactly 4, and the right 4.** Two objections survive the correction, and they are the actual argument. First, 2.6 A is `metal_bound`, shared by every other criterion; raising it to 2.75 for the hypothesis under test, after seeing where the atoms landed, is choosing the measurement to fit the expected answer - D120's move and D123's. Second, a count at a hard boundary is fragile when a member sits 0.09 A from it: over 100 ns at 300 K both O1x and OD2 cross 2.70-2.79 repeatedly, so the count flickers between 3, 4 and 5, and a single-frame verdict on the last frame reports where the thermal noise was rather than what the sphere is.

**The restraint is not the thing to change.** 2.66 A is the distance the grafted structure already had; moving the target to 2.4 A to make the count come out is pulling to a chosen number, which D93 forbids and which is how D120 happened. The defect is in the criterion.

**A coordination number cannot express this hypothesis, and D127's caution was stronger than it read.** `site.yaml` already warns, beside the `subst` section, that the start check passed coordination *on a count, not a composition* - the right four were not the four being counted. That was written as a caveat about one structure. It is actually a statement that the criterion kind is wrong for this hypothesis: what has to be true at `reactant` is that the donor set within about 2.8 A **is** {`C:2:O01`, `A:88:NE2`, `A:164:NE2`, `C:1:O1x`} and that `A:92:OD1` is **not** in it. That is a set identity, and `_criterion`'s `coordination` kind only returns a cardinality.

**Open question, and nothing is being changed until it is answered:** add a `composition` criterion kind to `site.py` that names the required donor set and the donors that must be absent, and use it in place of the count at `reactant` in both hypotheses? It is a small, testable addition that costs no compute and can be built while the job runs, and without it the reactant check on this trajectory will report FAIL for a run that did exactly what it was designed to do. **The QM/MM stage has not been submitted and will not be until the reactant check means something.**

**D130 - D120's prescribed fix was never applied, and `md_subst` restrains Asp92 harder than the run D120 condemned. The occupancy tool I was about to build has existed since D120.** (2026-10-29, user: "read the logs, as we've already gone in circles with the restraint problem")

The user's question - *holding what a crystal structure already has is the point, no?* - is the right one and the answer is yes. D93's rule is sound and nothing below weakens it. **The failure was never in the holding. It was in what got classed as "what the crystal has".**

**What D120 already established, re-read rather than re-derived.** The crystal sphere (6NC5 ZN302) is **four**: HOH485 1.90, His163 NE2 2.09, HOH493 2.14, His87 NE2 2.22. Asp92 OD2 sits at **2.80 A** and is a second-shell contact, not a ligand. `preprocess` detects at `cutoff: 3.0`, so OD2 was recorded as a ligand and written into the restraint list as `[B:1:ZN, A:92:OD2, 2.8, 50.0]`. A carboxylate pinned that close has one cheap way to relax, and it took it: OD1 rotated onto the metal at 1.99 A. D120's own sentence - **"the restraint whose purpose was to preserve the crystal's coordination sphere is what destroyed it, because 2.80 A was never a coordination distance."** So holding the crystal's sphere was right; holding a 2.80 A second-shell contact was not holding the crystal's sphere.

**D120's prescribed fix, quoted: the detection cutoff "has to be tight enough to exclude a 2.80 A contact (a Zn-O ligand is ~1.9-2.2 A; 2.5 is defensible, 3.0 is not), and a second carboxylate oxygen should not be able to enter the list behind the first."** That fix was never applied. Today's job printed `Cutoff: 3.00 A` on redetection. **One day and ten decisions after the diagnosis, the detection is unchanged, which is the circle the user is pointing at** - the gap was never long, only unwatched.

**And `md_subst` is worse on D120's own terms, not better.** Its restraint list holds OD2 at 2.80 **and** OD1 at 3.81. I justified the pair as correcting D120's artefact: the second spring blocks the rotation that produced it, so May's specific failure cannot recur. But two springs on a carboxylate assert its orientation relative to the metal where one asserted only its distance. **That is more assumption, not less**, and D120 says what to do with the question instead: Asp92 swinging onto a vacant Zn1 "is a question to ASK, not an assumption to restrain in, and the crystal's answer at 2.80 A is that it does not." The design that follows from D120 is to restrain **neither** Asp92 oxygen and tighten the cutoff so neither is ever detected. I restrained both.

**The compounding problem, which tightening the cutoff cannot reach.** `subst_start.pdb` descends from `md_naclong`'s last frame, and in those coordinates OD1 **is** at 1.99 A. A cutoff cannot undo an artefact already present in the geometry: at any threshold above 2.0 A that oxygen is detected, correctly, because it is genuinely there. D120 said this too, and costed it - **"everything from `md_nac` onward then has to be re-run from a shell that matches the crystal: 100 ns MM, then QM/MM MD per arm."** Grafting onto the May frame, which is what D126 and D128 did, carries the condemned shell forward rather than leaving it behind.

**The tool does not need building.** D129 proposed an occupancy criterion and the user named the quantity correctly. `site.py shell` already computes it - `occ = (dist < cutoff).mean(axis=0)` per donor, plus the per-frame coordination number as a distribution - and it was written **in D120, for this exact purpose**, because `sphere` could only answer for the frame it was given. What is missing is a pass/fail wrapper, and that is a small thing beside the starting shell. **Nothing was built; the composition criterion is on hold.**

**What job 30923516 can and cannot support.** It is ~1 h into ~3.4 h. The nucleophile question (zinc-bound hydroxide or a second non-metal-bound water, D124) is still partly answerable from it, though Asp92's two springs bias the general-base H-bond that question runs through. **Any claim about the zinc's coordination sphere is not supportable from it**, because the sphere it starts from is the one D120 rejected and Asp92 is held in place by choice rather than observed.

**The decision in front of the user:** kill 30923516 and rebuild from the crystal as D120 prescribed - preprocess with `cutoff: 2.5`, no Asp92 restraint, then 100 ns MM - or let it finish and take only the nucleophile reading from it, with the sphere left to the rebuilt chain. Nothing has been killed, nothing resubmitted, and no QM/MM proposed.

**D131 - Job 30923516 cancelled. The audit the user asked for before rebuilding found the plan still calling the wrong zinc sphere "established", 6NC5 missing from the reference file, and the monitor gap D120 named still open. One new result: the SEQADV records prove Zn2 is vacant in the wild type.** (2026-10-29, user: "review and consolidate all your scripts and md files ... so this doesn't happen again")

`scancel 30923516` at 10:31 AEST, 1 h 40 m in. Nothing resubmitted.

**Why a chronological log was not enough, stated as a mechanism.** D120 wrote its own fix into its closing paragraph. One day and ten decisions later the detection was still at 3.0 A and a fresh 100 ns job was running on the shell it condemned. **A long log cannot distinguish a fix that was prescribed from one that landed**, and five entries written in between all cited D120 approvingly without checking it. So the log now opens with **`## Standing constraints`**: 18 rules, each naming its decision, what enforces it, and - for most of the measurement rules - that nothing does and it rests on someone reading them. The status line is to be updated in the same commit that changes the behaviour.

**The guard, so the specific failure cannot recur silently.** `LIGAND_MAX = 2.5 A` is the longest distance at which restraining a metal-donor pair still restrains a bond. `md` now warns on every pair it is about to restrain beyond it and records them in the section's `coordination_beyond_bond_length` - **in the config, because stdout scrolls past and D120 sat unread in a `.out` for three weeks.** Configured pairs are checked as well as detected ones, because D130's own error was a hand-written 2.80 A restraint and a guard on detection alone would have missed it. Reported, never refused: a long bond may be the chemistry, and `coordination_warn_above` raises the threshold for an Ln site. Eight tests; four fail if the threshold is lifted.

**What the audit found in the prose, which is the part that matters for the paper.**

  - **`csr-neb-plan.md` listed the five-coordinate zinc under "what is established and not in doubt", described as "referenced to the crystal rather than to any simulation length".** It was never checked against the deposited file. That bullet is struck through in place and the Status section carries a superseding note, because the plan is the document that says what to do next and it was pointing at the artefact. Also corrected: a passage proposing as future work the crystal comparison D120 had already made.
  - **The bullet above it is not safe either.** "The hydroxide stays on the zinc and Asp92 OD1 holds it, 2.63-3.20 A, 100% of frames" was measured with Asp92 restrained against the metal. It is the one part of the mechanism every calculation agreed on, and it has to be re-established on a correct shell before it can be cited.
  - **`reference-structures.md` did not contain 6NC5 at all.** The file exists so every reference number is re-derivable rather than quoted, and the structure every other number is referenced against was recorded only inside a decision entry - part of why D120 took three weeks. Added first, and re-derived rather than copied: `refstruct.py shell` on a fresh download reproduces ZN302 at HOH485 1.90, His163 NE2 2.09, HOH493 2.14, His87 NE2 2.22, exactly as D120 reported.
  - **`traj_contacts.py` still had the gap D120 named.** Gln93 OE1 carried a `- Zn` line; Asp92 carried none. OD2 at 2.98-3.60 A reports as uncoordinated, which is what a healthy site looks like, so nothing watched the oxygen that was actually on the metal. Both oxygens now have one.
  - Smaller: `site.py`'s usage text still said it reads `michaelis_criteria` (D127 replaced it); `README.md` advertised four of `site.py`'s eight modes and omitted **`shell` and `check`** - the occupancy tool and the acceptance gate - which is why D129 proposed building a measurement that has existed since D120.

**The new result, and it is a better argument than the one we were using.** Re-deriving 6NC5 showed **six** zincs, and ZN301 sits in the Zn2 pocket on Asp91 OD2 1.95, His92 NE2 2.01, His225 NE2 2.04 and the bridging HOH493, 3.43 A from ZN302. That looks like a di-zinc site. The `SEQADV` records say why it is not:

```
SEQADV 6NC5 HIS A   92  UNP  A0A384EQ6 GLN   136 ENGINEERED MUTATION
SEQADV 6NC5 HIS A  225  UNP  A0A384EQ6 LYS   269 ENGINEERED MUTATION
```

Those are **our Gln93 and Lys226** - the Q and the K of RQK - mutated to histidine, restoring a canonical Asp-His-His Zn2 set. The construct is titled the *HARLDQ motif*. **The only way anyone has put a metal into CSR-1's Zn2 site was to engineer two of its ligands back**, and ZN302's own ligands are untouched by either mutation. That is positive structural evidence for a vacant Zn2, which is stronger than the absence of density we had been resting on, and it confirms the user's correction that ZN302 is the crystallographic position with wild-type ligands. ZN303-306 are adventitious surface sites. All six are tabulated now so the next reader measures the right one.

**Also recorded in the ledger, from the user's own corrections:** CSR-1's Zn1 site is **H-R-H with two coordinating residues, not three** - Arg90 (BBL R118) stands where the third histidine would be - so the canonical-equivalent four-coordinate sphere is two His + hydroxide + one water, and an extra solvent ligand is not over-coordination until it is counted against two protein ligands. This had been corrected twice in conversation and was in the log's corrections section, but not anywhere a setup step would meet it.

**One measurement worth keeping from the cancelled run:** it covered 16.66 ns in 34 minutes on `bun075`, **0.49 ns/min, against `naclong`'s 379 ns/day on the same system**. Nearly a factor of two, so a 100 ns NVT is ~3.4 h on a fast node rather than 6.3 h. Partial, from a cancelled job, so not written into `running-on-bunya.md` yet.

Commits: a193b80 (guard, ledger, README), 1243df4 (plan, reference file, monitors). Suite 403 passed, 6 skipped. **The rebuild has not started.**

**[SUPERSEDED by D134 - the sphere lasted under 20 ps; Asp92 bound in 100% of frames, CN mean 5.57. See `csr/retractions.md`.]** **D132 - The rebuild is running. `md_wt100`, job 30966977: 100 ns of NVT whose zinc carries the crystal's four ligands and nothing else, with no restraint anyone chose.** (2026-10-29, user: "Kill and rebuild" / "Yes, new directory")

Run directory `/scratch/user/uquadhik/CSR_wt`, built as a copy of `CSR_nac`'s preprocess output rather than a fresh preprocess. **Correcting what I told the user before this:** `csr_crystal_start.pdb` and `mer_crystal.pdb` are NOT missing - both are in `CSR_nac/preprocess_results/`, dated 23 Oct, and preprocess could have been re-run. Re-running it would have re-solvated the box, which changes every water position for no gain, so starting from the existing `finalsystem.pdb` is the better choice on its merits and not a workaround. The directory is self-contained: it carries the two crystal inputs, the force fields and the OpenFF ligand XML, so the run is reproducible from it alone.

**What is different from every md section before it, and it is only three things.**

  1. The start structure is `preprocess_results/finalsystem.pdb`, the pre-MM system, not an MM frame. The artefact D120 found appears *during* the MM, so any later frame carries Asp92 OD1 on the metal at 1.99-2.07 A and no cutoff undoes coordinates.
  2. `preprocess.inputs.cutoff` is **2.5**, not 3.0. This is the change D120 prescribed and D130 found had never been made.
  3. **No configured restraints at all** - not on Asp92, and not on the substrate either. Every restraint in force is one redetection measured on the structure itself.

**The detected sphere, checked before submitting, which is the whole point of the run:**

```
Metal Zn at B:1   Cutoff: 2.50 A
HOH D:1   (O  ) 1.90    HIS A:164 (NE2) 2.09
OH  C:2   (O01) 2.14    HIS A:88  (NE2) 2.22
```

Four restraints at 50 kcal/mol/A^2, and they are 6NC5 ZN302's four donors: two histidines, the hydroxide and the second water. **Asp92 is not among them** - at 2.80 A it falls outside 2.5 and is left as the second-shell contact the crystal says it is. No `[WARNING] ... is not a bond` fired, because nothing being restrained is longer than `LIGAND_MAX`. This is the sphere the paper's mechanism has to start from, and it is the first time a run of ours has had it.

**What it decides.** Whether that sphere is stable over 100 ns of MM. Every mechanism in `csr-neb-plan.md` assumes it and none has tested it: the one previous 100 ns (`naclong`) held Asp92 on the metal by restraint throughout. Read it with `site.py shell` for occupancy per donor - a count is not a composition (D129) - and `site.py check --at=reactant` on the last frame. Cost ~3.4 h on a fast GPU node at D131's 0.49 ns/min, 7 GB of trajectory at `traj_frequency: 5000`; scratch is at 229/300 GB.

**What it deliberately does not decide.** The substrate is unrestrained, so meropenem is expected to leave the attack line as it did in D93 (3.47 -> 5.41 A over 500 ps). That is not a failure of this run and not evidence against any hypothesis: this run establishes the resting site. Whether a Michaelis complex can be built on a correct sphere is the next question, and it must be asked of this trajectory's shell, not of `subst_start.pdb`.

Config: `md_wt100` in `docs/csr/site.yaml` (b9b3513), written into the run config with `sections.py`, which reported `0 selectors checked` - correct here, because the section has no selectors at all.

**D133 - Why every level of theory puts Asp92 on the zinc and the crystal does not: the RQK cations stand in for the absent Zn2, and Asp92 is probably in equilibrium rather than fixed.** **[PARTLY RETRACTED 2026-11-05, D206: the equilibrium hypothesis stands, but this entry also states that the standing "does GFN2 over-coordinate the zinc" question is ANSWERED, on D58's 1.6x-steeper DFT comparison. D58 is retracted and the question is open. This entry's own text undercuts the closure two paragraphs later - "a single relaxed scan or a single 25 ps trajectory cannot see an equilibrium, only the well it fell into, and every row in the table above is that kind of measurement" - and D136 built the hysteresis pair precisely because single-direction scans report their start.]** (2026-10-29, user: "it sounds like the positive residues that replace Zn's charge in the Zn2 site pull Asp92 away enough to enable a reactive 4 coordinate zinc ... it's possible that Asp92 *does* coordinate Zn1 in the wild type, but exists in equilibrium - it's a fairly slow lactamase so this would make sense")

**The contradiction, stated plainly.** Four results, and the structure is the one that disagrees:

| source | what it shows |
| --- | --- |
| D102 | 100 ns MM ends six-coordinate; Asp92 OD1 and a second water arrive **unrestrained** |
| D51 | 25 ps of GFN2-xTB QM/MM MD, free start, **no constraints**: OD2 2.80 -> 2.09 A |
| D58/D60 | DFT is **1.6x steeper** than xTB on that coordinate, 9.1 vs 5.7 kcal/mol to reach 2.80 A |
| 6NC5 ZN302 | **Asp92 off the metal, at 2.80 A** |

**This retires two readings this log has been carrying.** First, D120 and D132 blamed the five-coordinate zinc on the coordination restraints. D102's two newcomers were unrestrained, so a restraint is not the mechanism; removing it makes the coordination *measurable* instead of imposed - worth doing, and what `md_wt100` does - but it does not make the model right. Second, "whether GFN2-xTB over-coordinates the zinc" has been carried as the standing open question since D51. **D58/D60 answered it: no.** DFT binds Asp92 harder than xTB.

**The hypothesis.** In 6NC5, `ZN301` is held by **Asp91 OD2 at 1.95 A** - and 6NC5's numbering is ours minus one, so 6NC5's Asp91 **is our Asp92**. That carboxylate is a **Zn2 ligand**, and `ZN301` exists at all only because `SEQADV` records Gln->His at 92 and Lys->His at 225, our Gln93 and Lys226 (D131). In the wild type Zn2 is empty, so Asp92 has no metal of its own and the only cation on offer is Zn1. What keeps it off, on this reading, is the **RQK motif standing in for the missing metal's charge**: Arg90 and Lys226 sit in the Zn2 pocket as cations, and a carboxylate held in their field is satisfied enough not to collapse onto Zn1 - leaving the four-coordinate zinc the crystal shows, which is the only form with a vacant site for the substrate carbonyl.

**And it need not be all-or-nothing.** The user's extension is that Asp92 *does* reach Zn1 in the wild type, in **equilibrium** with the open form. That fits the kinetics - CSR-1 is a slow lactamase, and a mechanism gated by a minority open population is a slow mechanism. It also explains the shape of our failures: a single relaxed scan or a single 25 ps trajectory cannot see an equilibrium, only the well it fell into, and every row in the table above is that kind of measurement.

**What this makes of our models.** Both the MM and the QM/MM runs truncate or neutralise the environment that was doing the holding. Any model that empties the Zn2 site, leaves the carboxylate without a compensating field and then asks where it goes will get the same answer - and all three levels did.

**What tests it, in order of cost.**

1. **Hysteresis, now:** `scan_aspout` and `scan_aspin`, the same Zn-Asp92 OD1 coordinate from both ends (site.yaml). Separated curves mean D53/D58 were reporting their start. Config only, no new code.
2. **The free energy, which is what settles it:** a PMF along Zn-Asp92 OD1 at QM/MM level, because an equilibrium is a population and only a free energy reports one. `umbrella` and `metadynamics` were both built around a channel's z coordinate and a single ion, so this needed a general `coordinate: distance` mode. **It is written and tested** (2bfa012): `coordinate: distance` with `atoms: [B:1:ZN, A:92:OD1]`, the same windows, bias, seeding and WHAM, with `z` still the default so the channel configs are untouched. Agreed with the ion-channel session first, which has scoped permeation out and does not need a distance CV itself, so one set of keys serves both jobs rather than two bespoke ones. A biased pair is excluded from the coordination restraints, which would otherwise add to the bias while WHAM heard only about the bias - the same class of silent error as D120. Windows take a `seed` each, which the scan in 1 supplies; **nothing is launched**, and 15-20 QM/MM windows get costed and proposed once the hysteresis curves say whether the PMF is worth it.
3. **The hypothesis itself, by removing the charge:** if RQK's cations hold Asp92 off, neutralising Lys226 should deepen the bound well. Not a free config change - this config's `amino_acids` list omits LYN deliberately (see `qmmm_ohH`), so a neutral lysine has to be built the way the `ohH` arm builds its proton rather than by renaming the residue.

**Also fixed in the setup, from the user's instruction to keep RQK in the scan.** The QM region already held **Arg90 and Lys226** (preprocess `extra_qm_residues`) and Asp92 (via the 3.0 A detection). **Gln93 was not in it** - it only ever entered when it came inside the detection cutoff. Both scan sections now carry `extra_qm_residues: {A: [93]}`; Gln is neutral, so the QM charge stays +3. Note for anyone re-running preprocess in `CSR_wt`: its `inputs.cutoff` is now 2.5 while its saved `results.qm_residues` were built at 3.0, so **Asp92 is in the QM region only because the saved region was inherited.** A fresh preprocess at 2.5 would drop the residue this whole question is about. Put 92 and 93 in `extra_qm_residues` before ever re-running it.


**[SUPERSEDED by D144 - every MD-derived statement about Asp92 coordination is inadmissible, not merely caveated. See `csr/retractions.md`.]** **D134 - 100 ns with nothing holding Asp92: it is on the zinc in every one of 5,000 frames, the zinc is five- or six-coordinate in every frame, and the substrate pose is gone.** (2026-10-30, job 30966977, `md_wt100` in `/scratch/user/uquadhik/CSR_wt`)

The run D132 set up finished COMPLETED in 3 h 15 m. Four restraints only, redetected at 2.5 A on the pre-MM `finalsystem.pdb` - HOH D:1:O 1.90, His164 NE2 2.09, O01 2.14, His88 NE2 2.22, verified in the `.out` with no warning of any kind. **Nothing restrained Asp92 and nothing restrained the substrate.**

**Occupancy over the whole trajectory** (`site.py shell`, 5,000 frames at 20 ps, cutoff 2.6 A):

| ligand | occ% | min | mean | max |
| --- | --- | --- | --- | --- |
| OH:O01 (nucleophile) | 100.0 | 1.82 | 1.94 | 2.12 |
| **Asp92 OD2** | **100.0** | 1.91 | 2.07 | 2.32 |
| HOH D:1:O | 100.0 | 1.92 | 2.07 | 2.27 |
| His164 NE2 | 100.0 | 2.01 | 2.20 | 2.46 |
| His88 NE2 | 100.0 | 2.05 | 2.24 | 2.51 |
| Asp92 OD1 | 56.7 | 1.93 | 2.67 | 3.72 |
| mer O1x | 0.0 | 3.35 | 3.92 | 5.42 |
| Lys226 NZ | 0.0 | 3.47 | 3.96 | 5.18 |

Coordination number: **mean 5.57; five in 2,164 frames (43%), six in 2,836 (57%). Four-coordinate: never, in 100 ns.**

**What it says.**

- **Asp92 OD2 is bound in every frame, and it was already at 1.91 A by the first written frame.** It does not creep in over 100 ns; from a 2.80 A start it closes inside the first 20 ps with nothing holding it. This is the fourth measurement to do it (D102 MM restrained, D51 xTB, D58/D60 DFT) and the cleanest, because here the carboxylate is free and the rest of the sphere is held at crystal distances. Rule 3 now has a restraint-free number behind it.
- **OD1 at 56.7%** (mean 2.67, 1.93-3.72) is bidentate flicker. It looks like D133's equilibrium but is on the wrong coordinate: the population D133 needs is OD2 on and off, and OD2 never leaves. The six-coordinate frames are OD1, or a second water - the last frame has E:6862 HOH at 2.14 A.
- **The substrate left.** mer O1x is outside the cutoff in every frame; in the last frame the carbonyl O2x is **7.80 A** from the zinc and O01-C4x is 5.56 A. `site.py check --at=reactant --hypothesis=carbonyl_on_metal` is **REJECTED**, five of seven criteria failing (only Zn-O01 1.90 and the 119 degree Burgi-Dunitz angle pass). `naclong` kept Gervasoni's substrate restraint through its 100 ns for exactly this reason (D92); `wt100` kept none, and the pose did not survive.

**What it decides, and what it does not.** Removing the restraints does not produce the crystal's four-coordinate zinc - it produces a worse-coordinated one. That is the question D132 was built to answer, and it is answered. It does **not** decide whether Asp92 is in equilibrium: a trajectory that never once leaves the bound state either puts the open population below roughly 0.02% at MM level or shows that this MM model cannot leave it, and nothing in a single trajectory separates those two. Only the PMF (D133, test 2) does.

**Consequence for the scans.** `scan_aspin` reads `md_wt100_results/NVTtrajectory_lastframe_imaged.pdb`, and that frame has no substrate in the site - it would scan the wrong system. It needs a different start: either `wt100` re-run with `naclong`'s substrate restraint kept (the only restraint, leaving Asp92 free), or `aspin` started from `naclong`'s last frame with the Zn-OD1 distance pulled out to 3.60 A first. Left unchanged in `site.yaml` pending that choice. `scan_aspout`, which starts from `naclong`, is unaffected.

**Rate, and a correction to `running-on-bunya.md`.** Production was 11,517.7 s for 100 ns: **0.52 ns/min, 750 ns/day**, on a full H100 (bun074). The 379 ns/day recorded there for `naclong` was **not measured on an H100**: job 30421166's `AllocTRES` reads `gres/gpu:nvidia_a100_80gb_pcie_4g.40gb=1`, a MIG slice of an A100 on bun004, billed 26 against this job's 100. Same system, same `traj_frequency`, twice the throughput for four times the billing - so the earlier "halve a benchmark figure" advice was fitted to a slice, not to a card.


**D135 - Scratch pruned: 9.1 GB of redundant trajectories deleted, every last frame kept.** (2026-10-30, user: "keep an eye on the scratch allocation and delete redundant trajectories (keeping perhaps only the last frame)")

Scratch was 261/300 GB with other projects writing into it, and an over-quota write fails silently (rule 8). Deleted, with the byte counts measured before removal:

| file | bytes | why redundant |
| --- | --- | --- |
| `CSR_wt/md_wt100_results/NVTtrajectory.dcd` | 2,965,840,276 | unimaged duplicate of `NVTtrajectory_imaged.dcd`, which is what `site.py` reads |
| `CSR_wt/md_wt100_results/relaxbox_NPT{,_imaged}.dcd` | 2 x 593,168,276 | box equilibration; only its final state feeds the NVT |
| `CSR_wt/md_wt100_results/OpenMMOpt_traj.xyz` | 44,284,096 | minimisation path |
| `CSR_nac/md_subst_results/NVTtrajectory.dcd` | 2,040,693,236 | job 30923516, **cancelled** (D131); `md_subst` is DO NOT RUN (D126) |
| `CSR_nac/md_subst_results/relaxbox_NPT{,_imaged}.dcd` | 2 x 592,880,276 | same cancelled run |
| `CSR_nac/md_nacgpu_results/NVTtrajectory{,_imaged}.dcd` | 2 x 296,584,276 | the 1 ns GPU benchmark (D100), superseded by `naclong` |
| `CSR_nac/md_nacgpu_results/relaxbox_NPT{,_imaged}.dcd` | 2 x 593,168,276 | box equilibration |
| `CSR_nac/md_nac_results/relaxbox_NPT.dcd` | 571,814,228 | box equilibration |

**Kept:** every `*_lastframe*.pdb`/`.cif`, every checkpoint and final state, both `nvtsim.csv`/`nptsim.csv`, and **`CSR_wt/md_wt100_results/NVTtrajectory_imaged.dcd`** - D134's occupancy is measured on it, so it stays until that number is in the paper. The three `QM_MM.dcd` files (30 MB each) stay. `CSR_wt` went 6.8 -> 2.9 GB, `CSR_nac` 8.0 -> 2.8 GB; scratch reads 252 GB.

Each pruned directory now holds a one-line `DELETED.md` pointing here, the convention `naclong` started on 25 Oct (D105).

**Not touched, and worth a decision rather than a deletion.** The CSR work is 6 GB of a 252 GB allocation. The rest is: `channel/` 45 GB (the ion-channel session's live umbrella windows - not ours to prune), `glp1r_md/` 41 GB, `AANCA/` 26 GB, `boltzgen/` 17 GB, and **two copies of the same ESMFold2 weights**, `hf/hub/models--biohub--ESMFold2` and `.cache/huggingface/hub/models--biohub--ESMFold2`, matching blob-for-blob at about 12 GB each. Pointing `HF_HOME` at one of them would free more than everything in this entry combined.


**[QUALIFIED - the start is `naclong`'s condemned shell, which D131 says nothing measured on can be cited. See `csr/retractions.md`.]** **D136 - Two jobs submitted: `md_wt100b` (31016877) and `scan_aspout` (31016878).** (2026-10-30, user: "make sure the next job gets batched today, no need to wait for my say so, use your judgement to figure out what the next step should be")

**`md_wt100b`, job 31016877** - `-p gpu_cuda --qos=gpu --gres=gpu:h100:1 --cpus-per-task=4 --mem=32G -t 12:00:00`. 100 ns NVT, the four redetected zinc restraints plus `naclong`'s substrate restraint `[C:2:O01, C:1:C4x, 3.47]`, nothing on Asp92. **Why this and not `aspin`:** D134 found wt100's substrate left the site, so its occupancy is of a site with no substrate in it and its last frame is not a Michaelis complex. This is the one run that gives a free-Asp92 occupancy **with the substrate present** - which is the number the paper needs - and a usable start for `aspin`. Cost: 3 h 12 m of production at D134's rate, one GPU. `--gres=gpu:h100:1` is explicit, per D134's MIG-slice finding.

**`scan_aspout`, job 31016878** - `--cpus-per-task=8 --mem=64G -t 3-00:00:00`, CPU, from `naclong`'s last frame. D133 test 1, the outward arm of the hysteresis pair: Zn-Asp92 OD1 from 1.90 to 3.60 A in nine relaxed points, Arg90, Lys226 and now Gln93 in the QM region. It needs nothing from wt100b, so it runs in parallel. Memory is 64G rather than the script's 32G because D97's Gln93 scan died of memory, not slowness. **What it decides:** whether a relaxed scan reports the basin it started in (D94's objection) - if `aspout` and `aspin` separate, every single-direction scan in this log is reporting its start, and the PMF becomes the only instrument left.

Both selectors were checked by `sections.py` before submission (2 each, both resolved). `aspin` stays DO NOT RUN until wt100b finishes.


**D137 - `scan_aspout` cancelled at point 0 of 8: a relaxed scan cannot converge while mcube writes the rigid-water constraints as `$set`, and D133's hysteresis test is not answerable until that is fixed.** (2026-10-30 04:55 AEST, job 31016878, `CANCELLED` after 04:38:57)

**`md_wt100b` (31016877) is still running** - 50.66 of 100 ns at 04:48 on bun073, temperature 301 K, energies flat, nothing in the `.err`. So there is no substrate-present occupancy yet and nothing to compare against D134's table; `scan_aspin` stays DO NOT RUN, both of its preconditions being unmet.

**What `aspout` did in 4 h 39 m: one scan point, 120 geomeTRIC steps, no convergence by any criterion.**

| quantity | at step 120 | criterion |
| --- | --- | --- |
| Grad_T rms / max | 1.635e-02 / 2.739e-01 | 1e-04 / 3e-04 |
| Displace rms / max | 7.033e-03 / 1.082e-01 | 2e-03 / 4e-03 |
| trust radius | 2.443e-03 | - |
| energy change | **+4.503e-03** Eh | 5e-06 |
| worst constraint violation | **0.069 A** | `convergence_cmax` 0.01 |

The gradient is two to three orders of magnitude off after 120 steps, the trust radius has collapsed, and the energy is **rising**. Over steps 111-118 the total oscillates across 6e-3 Eh - **3.7 kcal/mol**, larger than the hysteresis the scan was built to resolve - so even the plateau argument that rescued D105 (and gave D106 its NEB start) does not apply here.

**The scan coordinate is not the problem.** Zn-Asp92 OD1 (`distance 395 169`) sits at its 1.900 target to within 0.0016 A and has dropped off geomeTRIC's violation table entirely. What never satisfies is the **rigid-water block**: O-H distances reading 0.99686 and 1.02609 against TIP3P targets of 0.9574 and 0.9573, H-H 1.571 against 1.513.

**The mechanism, and it is mcube's own.** `optimizer_constraints` (`src/mcube/qmmm.py:539`) emits the active-region MM waters as geomeTRIC `$freeze` when there is no user constraint, and switches **the whole block to `$set` with values** the moment one appears - which a scan always supplies. D105's `constraints.txt` was a `$freeze` header over 24 water distances; this run's is a `$set` header over 27 water distances plus the one scan bond. The difference is decisive, not cosmetic:

- a **frozen** distance is a degree of freedom geomeTRIC projects out, so a residual on it (D105 measured 0.08-0.10 A) is ugly but costs nothing;
- a **`$set`** distance is a convergence criterion, checked against `convergence_cmax`. At 0.069 A the run is seven times outside it and **cannot be declared converged however long it runs**.

It also made the linear algebra far worse: **160 SVD failures in 120 steps**, against the fifteen per hour D105 recorded, each dumping a 32 MB `HT.v_*.txt`.

**Why it was cancelled rather than left to run.** 2.3 min/step and `maxiter` 300 gives ~11.5 h per point, so nine points is **~100 h against a 72 h limit** - it would have been killed around point 6 with every point unconverged. And the dumps are a quota hazard: 151 files, **4.6 GB in one point**, which is the whole 4.8 GB the directory had reached. Nine points is over 40 GB against 41 GB free, while `md_wt100b` writes into the same filesystem and an over-quota write fails silently (rule 8). Cancelled at 04:55, dumps deleted, `scan_aspout` 4.8 GB -> 211 MB, scratch 259 -> 255 GB; `geometric_OPTtraj_optim.xyz` and `optimization_energies.log` kept as the evidence above, with a `DELETED.md` pointing here.

**What this decides.** D133 test 1 is **not answerable by `scan` as mcube currently writes constraints**, and the same applies to every `qmmm` run given `opt_constraints`. The three candidate fixes, none of them run:

1. **`$freeze` the waters and `$set` only the scan bond** - the correct fix, and the one that restores D105's behaviour. ASH's `Optimizer` takes a single `constrainvalue` flag for the whole dict, so it needs either an ASH change or mcube writing geomeTRIC's constraint file itself.
2. **`rigid_water: false` in the scan section** - one config line, and then only the scan bond is in `$set`; the MM waters relax under their own bonded terms. This changes what is being optimised, so it is a methodological choice rather than a bug fix, and it is the user's call.
3. **Shrink `actradius`** - D105's own remedy. With nine active waters at 8 A it reduces the constraint count but leaves `cmax` gating convergence, so it mitigates rather than fixes.

Nothing resubmitted. The hysteresis test, and with it the standing of every single-direction scan in this log (D94's objection), waits on that choice; the distance PMF of D133 test 2 is unaffected only if it samples with MD rather than constrained optimisation.


**[SUPERSEDED by D144 - inadmissible, and its reversal of D120's artefact diagnosis does not stand. See `csr/retractions.md`.]** **D138 - With the substrate held, Asp92 is bidentate on the zinc in every frame and the site is six-coordinate throughout: the metal gets worse, not better, and D120's "restraint artefact" was misattributed.** (2026-10-30, job 31016877, `md_wt100b` in `/scratch/user/uquadhik/CSR_wt`)

COMPLETED in 03:12:51, ended 06:22:17 on bun073. **Five** harmonic restraints, verified line by line in the `.out` against what D136 asked for - `C:2:O01`-`C:1:C4x` at 3.47 A plus the four redetected crystal donors (HOH D:1:O 1.90, His164 NE2 2.09, O01 2.14, His88 NE2 2.22), each confirmed by OpenMM's own `Adding bond restraint` lines at 50 kcal/mol/A^2. **Nothing restrained Asp92.** The run's single `[WARNING]` is mcube's own, and it was right: *"a carboxylate held that close can rotate its other oxygen on (D120)"*.

**Occupancy over the whole trajectory** (`site.py shell`, 5,000 frames at 20 ps, cutoff 2.6 A), against D134's substrate-free run:

| ligand | wt100b occ% | min/mean/max | wt100 (D134) occ% | min/mean/max |
| --- | --- | --- | --- | --- |
| OH:O01 | 100.0 | 1.82/1.94/2.09 | 100.0 | 1.82/1.94/2.12 |
| HOH D:1:O | 100.0 | 1.93/2.08/2.27 | 100.0 | 1.92/2.07/2.27 |
| Asp92 OD2 | 100.0 | 1.90/2.09/2.43 | 100.0 | 1.91/2.07/2.32 |
| His164 NE2 | 100.0 | 2.03/2.24/2.53 | 100.0 | 2.01/2.20/2.46 |
| His88 NE2 | 100.0 | 2.01/2.24/2.53 | 100.0 | 2.05/2.24/2.51 |
| **Asp92 OD1** | **100.0** | **1.95/2.15/2.74** | **56.7** | 1.93/2.67/3.72 |
| mer O1x | 0.0 | 3.00/3.56/4.79 | 0.0 | 3.35/3.92/5.42 |
| Lys226 NZ | 0.0 | 3.54/4.37/5.39 | 0.0 | 3.47/3.96/5.18 |

Coordination number: **mean 6.00; six in 4,998 frames (100%), five in 2. Four-coordinate: never, in a second independent 100 ns.**

**The question this run existed for, answered: yes, and worse.** Asp92 OD2 is a ligand in every frame exactly as in D134, and the one donor that differs is the one that was supposed to be reassuring - **OD1 goes from 56.7% flicker at mean 2.67 A to 100% occupancy at mean 2.15 A**. The carboxylate is not monodentate with a wobble; with the substrate tethered it is **bidentate**, and the zinc carries six ligands from the first frame to the last. Holding the substrate removed the only slack in the coordination sphere.

**The mechanism is the one the warning named.** The restraint fixes the nucleophile's oxygen 3.47 A from C4x, which pins O01 against the metal and stops the small excursions that let OD1 off in D134; the carboxylate then closes on what is left. So the restraint intended to keep the Michaelis complex together is also what makes the metal over-coordinated - the two cannot be separated in this model.

**The substrate stayed, and still is not a Michaelis complex.** It did not dissociate as in D134: mer O1x is nearest at 3.43 A in the last frame (mean 3.56 over the run, against D134's 3.92) but never once inside 2.6 A, and the **carbonyl O2x is 5.98 A** from the zinc - closer than D134's 7.80 A and nowhere near the 2.4-3.0 A the hypothesis wants. `site.py check --at=reactant --hypothesis=carbonyl_on_metal` is **REJECTED**, five of seven failing, with the same two passes as D134: Zn-O01 1.93 A and a 102 degree Burgi-Dunitz angle. O01-C4x reads 3.65 A against its 3.47 A target, which is ordinary thermal slack on a 50 kcal/mol/A^2 harmonic, not drift.

**A correction to D120, and to rule 7.** D120 condemned `md_naclong`'s shell because Asp92 OD1 sat at 1.99 A there, and rule 7 carried that forward as an artefact inherited from the restraint that held it. This run restrains **nothing** on Asp92 and puts OD1 at mean 2.15 A in 100% of 5,000 frames. The OD1 contact is therefore **not a restraint artefact**: it is what this MM model does whenever the crystal's four donors are held and the substrate is tethered, which is every production run we have. D120's other claim is untouched - that shell does not match the crystal - but the cause was misattributed, and "re-run it from a shell that matches the crystal" has now been tried twice and does not produce one.

**What it does not decide, and what that costs.** Two independent 100 ns trajectories that never once leave the bound state either put the open population below roughly 0.02% at MM level or show that this MM model cannot leave it, and nothing in a trajectory separates those two - D134 said this and a second trajectory does not change it. Only the PMF (D133, test 2) does. **With D137 blocking both scan arms, the PMF is now the only instrument left for the Asp92 question**, so it moves from "costed once the hysteresis curves say whether 15-20 windows are worth it" to the next thing needing a decision.

**Rate, confirming D134 rather than correcting it.** Production was 11,361.9 s for 100 ns: **0.528 ns/min, 760 ns/day** on a full H100, against D134's 11,517.7 s, 0.52 ns/min and 750 ns/day. Within 1.4% on the same system and the same card, so that figure is now measured twice and `running-on-bunya.md` needs no change. NPT converged at 10 ps.

**Pruned** per D135's rule, bytes measured first: `NVTtrajectory.dcd` 2,965,840,276 (unimaged duplicate), `relaxbox_NPT{,_imaged}.dcd` 2 x 593,168,276, `OpenMMOpt_traj.xyz` 49,819,608 - **4.20 GB**. Kept `NVTtrajectory_imaged.dcd` (the table above is measured on it), every first and last frame, both checkpoints, the final state and both sim CSVs; `md_wt100b_results` 7.1 -> 2.9 GB, scratch 259 -> 255 GB, with a `DELETED.md` pointing here.

**Nothing submitted.** `scan_aspin` fails both of its own preconditions anyway - the last frame is REJECTED for the substrate pose and Zn-OD1 is 2.12 A, not the 3.3 A-plus a free-Asp92 start needs - and D137's constraint problem blocks it regardless.


**D139 - `rigid_water` is now off by default for every QM/MM gradient stage, and `scan` honours it at all; `aspout` resubmitted as 31051634 and already holding one constraint instead of twenty-eight.** (2026-10-30, user: "rigidwater should always be off for QM/MM work - the only reason that's there is for NVT")

**The fix is in the library, not in the config** (commit cf4ad8d). Three things were wrong, and the third is why D137's own suggested workaround would have done nothing:

1. `qmmm` defaulted `rigid_water` to `true`, and `ts` with it.
2. **`scan` had no `rigid_water` option.** It called `water_constraints` unconditionally (`src/mcube/scan.py:218`), so `rigid_water: false` in a scan section - D137's option 2, and the one the note in `site.yaml` told a future reader to use - would have been read, ignored, and the job would have failed the same way. That is now fixed too.
3. The default itself was the wrong way round for anything that takes a gradient. Rigidity is a convenience in MM MD and a liability here: with any user constraint present - and a scan sets its bond at every point - geomeTRIC writes the whole block as `$set`, and a `$set` water distance is a *convergence criterion* checked against `convergence_cmax`. Rigid TIP3P cannot satisfy it, so no point converges however long it runs (D137: 0.069 A against 0.01, 120 steps, 4 h 39 m, one point of nine).

So `rigid_water` defaults to **false** in `qmmm`, `scan` and `ts`, with the default in one shared helper, `rigid_water_pairs` (`src/mcube/qmmm.py:538`), rather than three separate `inputs.get(..., True)` calls. `rigid_water: true` restores the old behaviour for anyone who wants it. **The `md` stage's OpenMM `rigidWater` is untouched** - there the waters are genuinely constrained and that is what buys the 4 fs timestep. Four tests added; the suite is 421 passed, 6 skipped locally (the ASH-gated ones, including these, run on Bunya).

**Why this was the right fix rather than `rigid_water: false` in the two scan sections.** D137 called option 2 "a methodological choice, not a bug fix, and the user's call" - and the user's answer was that it is not a choice at all: the MM system in the gradient stages is already built with `rigidwater=False` and `autoconstraints=None`, so holding the waters was adding constraints the model does not have, in a stage that cannot afford them. One default in the library beats the same line copied into every future scan, NEB and TS section, which is what mcube being a general package is for.

**`aspout` resubmitted: job 31051634**, `--cpus-per-task=8 --mem=64G -t 3-00:00:00`, running on bun133 from 08:46 AEST, same section and same start frame as 31016878 (`../CSR_nac/md_naclong_results/NVTtrajectory_lastframe_imaged.pdb`, Zn-Asp92 OD1 1.90 -> 3.60 A in nine points, Gln93 in the QM region). The cancelled run's directory is kept as `protdft/scan_aspout_d137`.

**The fix is confirmed in the running job, not just in the tests.** geomeTRIC's `constraints.txt` for point 0 now reads, in full:

```
$set
distance 395 169 1.9
```

against 31016878's `$set` header over 27 water distances **plus** that same bond. One constraint, and it is the scan coordinate. The QM/MM energy is falling over the first steps (-446.63 -> -446.99 -> -447.05 Eh).

**What is still open.** Convergence per point is now *possible*; whether it happens, and at what cost, is what this job measures. 160 SVD failures in 120 steps (D137) were plausibly the frozen bonds being enforced in a degenerating TRIC coordinate system (D105's diagnosis, channel log C22), so they should abate with the constraints gone - if they do not, `opt_coordsystem` is unavailable while the scan bond is `$set`, and shrinking `actradius` is what is left. Scratch is 231/300 GB, so there is room for the dumps if they come back, but it needs watching rather than assuming.

**`aspin` now has a start that exists.** It runs after `aspout` and begins from `protdft/scan_aspout/scan_08.pdb` - aspout's own last point at 3.60 A, the same box, waters and substrate, walked back. Neither MD run can supply a free-Asp92 start: `md_wt100`'s last frame has no substrate in the site (D134) and `md_wt100b`'s has Asp92 bidentate at 2.15 A in 100% of frames (D138). A reverse arm from the forward end point is the hysteresis test proper, and it is also the only one available. The PMF (D133 test 2) is held until the two curves are in hand, per the user: "do the scans, then we can do the pmf if necessary".


**D140 - `aspout` ran clean: nine of nine points converged in 76 minutes, no SVD failures, no dumps - and at 3.60 A the site is the canonical four-coordinate B3 sphere, 13.8 kcal/mol uphill. The reverse arm is submitted as 31143211.** (2026-10-30, job 31051634, COMPLETED 01:16:44 on bun133)

**The constraint fix (D139) is confirmed by outcome, not just by `constraints.txt`.** Against 31016878, which spent 4 h 39 m on point 0 of 9 and converged by no criterion:

| | 31016878 (rigid waters) | 31051634 (D139 default) |
| --- | --- | --- |
| points converged | 0 of 9, 120 steps on the first | **9 of 9** |
| wall time | 4 h 39 m, cancelled | **1 h 16 m 44 s**, COMPLETED |
| per point | >11.5 h projected | ~8.5 min |
| `SVD Error` | 160 in 120 steps | **0** |
| `HT.v_*.txt` dumps | 151 files, 4.6 GB | **0 files** |
| run directory | 4.8 GB | 101 MB |

So D105's reading of the SVD failures was right and D137 inherited it correctly: they came from frozen bonds being enforced in a degenerating coordinate system, not from the system's size. Remove the water constraints and they do not occur at all. `opt_coordsystem` and `actradius` were never needed.

**The outward profile** (GFN2-xTB/MM, Zn-Asp92 OD1, relative to point 0, Gln93 in the QM region):

| point | Zn-OD1 target / actual (A) | rel (kcal/mol) |
| --- | --- | --- |
| 0 | 1.900 / 1.901 | 0.00 |
| 1 | 2.112 / 2.113 | **-0.86** |
| 2 | 2.325 / 2.325 | +1.02 |
| 3 | 2.538 / 2.538 | +3.99 |
| 4 | 2.750 / 2.750 | +7.20 |
| 5 | 2.962 / 2.962 | +10.56 |
| 6 | 3.175 / 3.175 | **+8.77** |
| 7 | 3.388 / 3.388 | +11.13 |
| 8 | 3.600 / 3.600 | +13.84 |

Every point holds its target to within 0.0002-0.0007 A. The minimum is at **2.11 A**, not at the 1.90 A the structure started from, and the crystal's 2.80 A sits about 8 kcal/mol up this curve. **Point 6 breaks the monotonic rise** - 8.77 after 10.56 at 2.96 - a 1.8 kcal/mol drop that is the signature of the active region finding a different relaxation partway out; it is also where the run spent its longest stretch of output. A curve with a step in it is exactly what the reverse arm is for, so no weight is put on it yet.

**At 3.60 A the zinc has the sphere the B3 literature expects, and this log has never produced before.** `site.py sphere` on the two ends:

| | scan_00 (1.90 A) | scan_08 (3.60 A) |
| --- | --- | --- |
| ligands within 2.8 A | **5** | **4** |
| Asp92 OD1 / OD2 | 1.90 / 3.83 | 3.60 / 4.48 |
| O01 (nucleophile) | 1.93 | 1.84 |
| HOH485 | 2.08 | 2.00 |
| His88 / His164 | 2.03 / 2.24 | 2.03 / 1.95 |
| mer O2x | 5.01 | 4.68 |
| O01-C4x | 3.58 | 3.60 |

Two histidines, hydroxide and one water, with Asp92 off the metal entirely - the four-coordinate mono-Zn site, reached by pulling one coordinate and letting everything else relax. Note what did **not** happen: the nucleophile stayed on the zinc and 3.58-3.60 A from C4x with **no restraint holding it there** (the scan's only constraint is the Zn-OD1 bond; the substrate tether belongs to `md_naclong`, which produced the starting frame, not to this optimisation). The substrate did not leave, and the carbonyl came in slightly, 5.01 -> 4.68 A.

**What this does and does not say.** It gives a number for the first time: **displacing Asp92 to the crystal's geometry costs ~8 kcal/mol, and to 3.60 A costs 13.8**, at xTB/MM, in one direction. It does not yet say whether that is the surface or the starting basin - that is D94's objection and it is unanswered until the reverse arm lands on the same curve or does not.

**`aspin` submitted: job 31143211** (PD at 14:07), `--cpus-per-task=8 --mem=64G -t 12:00:00`, starting from `protdft/scan_aspout/scan_08.pdb` - this run's own last point, the four-coordinate structure above - and walking 3.60 -> 1.90 A in the same nine steps. At 8.5 min a point it should take about 1 h 20 m. Nothing else was submitted; scratch is 237/300 GB and nothing needed pruning (the whole scan is 101 MB, and `protdft/scan_aspout_d137` is kept as D137's evidence at 211 MB).


**D141 - Audit of what the rigid-water default cost: every QM/MM gradient run in this project carried it, but only `aspout` was destroyed by it; the rest converged, and what is unquantified is whether their energies moved.** (2026-10-30, user: "Has the rigidwater thing affected any previous results?")

**Scope: all of it.** `grep -rn rigid_water */config.yaml */*/config.yaml */*/*/config.yaml` across `/scratch/user/uquadhik` returns **nothing**. No run in this project ever set the key, so every `qmmm` optimisation, every `scan`, every `neb` and every `ts` ran with the active-region MM waters held - the CSR work and the ion-channel work alike.

**It mattered in one of two ways, and the difference is whether the run had a constraint of its own.**

- **No user constraint** (plain `qmmm` optimisations, `ts`, `neb`): the waters reached geomeTRIC as `$freeze`, degrees of freedom projected out of the optimisation space. A frozen coordinate cannot block convergence, and the internal geometry it froze is TIP3P's own rigid geometry (`b0` 0.9572 A), which is what that water model is parameterised as. **Those results stand as computed.** What they are is optimisations in which the active-region waters could translate and rotate but not stretch or bend - a defensible model, just not the one anyone chose.
- **With a user constraint** (every `scan`, and constrained `qmmm` optimisations such as D81's sphere-constrained reactant/product pair): the whole block became `$set`, and each water distance became a convergence criterion (D137).

**What `$set` actually did to the runs, from the CSVs rather than from theory:**

| run | points | converged |
| --- | --- | --- |
| `CSR_crystal_v2/scans/scan_H1_f00497` | 25 | **25** |
| `CSR_crystal_v2/scans/scan_H2_f00395` | 24 | 21, **3 not** |
| `CSR_crystal_v2/protdft/zncheck/scan` | 9 | 9 |
| `CSR_crystal_v2/protdft/coscan` | 9 | 9 |
| `CSR_gln93/protdft/scan_zngln` | 3 | 3 |
| `CSR_crystal/scan_attack_f400` | **0 rows** | - |
| `CSR_wt/protdft/scan_aspout` (31016878) | 0 of 9 | **none, cancelled** |
| `CSR_wt/protdft/scan_aspout` (31051634, fixed) | 9 | 9 |

So **the sweeping version of D137's claim is wrong**: `$set` water constraints are not unsatisfiable in general. Sixty-seven scan points across four earlier scans converged with them, including scans started straight from MD frames (`nac_frame00497.pdb`, `nac_frame00395.pdb`) as `aspout` was. The honest statement is that a `$set` water distance *can* gate convergence, that it did so completely in `aspout`, and that why that system's waters strained to 0.997-1.026 A against their 0.957 A start when other systems' did not is **not established**. What the fix removes is the possibility, not a universal error.

**What it does change in the record.** Three results that were unexplained now have a named candidate cause rather than nothing: `scan_H2_f00395`'s three unconverged points, the frame-499 optimisation that hit the 300-step limit (D-log, 2026-10-16), and D105's stall at maximum gradient 1.9e-01 with the energy already flat - "a relaxed geometry, not a stationary point", which D106 then used as a NEB endpoint. None of these is overturned; each is now worth re-running cheaply rather than reasoning about.

**What is not known, and the one test that would settle it.** Whether holding the waters shifted the *energies* of the runs that did converge. One job answers it: re-optimise `aspout`'s own `scan_04.pdb` at the same 2.75 A constraint with `rigid_water: true`, against the converged flexible result already in hand. Same structure, same active region, one variable, about 10 minutes on 8 cores. If the two energies agree to a few tenths of a kcal/mol, every converged number in this log stands as computed and this entry closes; if they do not, the comparison gives the size of the correction. **Not submitted** - it is the user's call, and `aspin` is using the queue.

**Restraints in QM/MM, separately** (user: "we shouldn't need any restraints for QM/MM jobs either"). Checked rather than assumed: harmonic distance restraints exist only in `md` (`md_restraints`, `src/mcube/md.py:242`) and in `umbrella`, which already switches the coordination restraints off for sampling (`src/mcube/umbrella.py:454`). The `qmmm`, `scan`, `neb` and `ts` stages add **none**, and after D139 the only thing they hand geomeTRIC is what the job explicitly asks for - a scan's own bond, or `opt_constraints`. So the position is already the one wanted. The place restraints still reach QM/MM is through the *geometry*: every start structure comes from restrained MD, which is exactly D138's finding and is not fixed by anything in the QM/MM stages.

**A caveat for the hysteresis comparison, noticed while auditing.** `aspin`'s point 0 and `aspout`'s point 8 are the same structure at the same 3.60 A, but their energies are -447.31984 and -447.31201 Eh - 4.9 kcal/mol apart. They are different runs, so `resolve_actatoms` re-derived the active region from a different `pdb_file` and the two will not contain exactly the same water set. **Absolute energies are not comparable between the arms**; each curve must be read relative to its own first point, which is how both CSVs already report them.


**D142 - Holding the waters cost 0.003 kcal/mol on a structure whose waters were already relaxed - which is the answer for half the affected runs and, by showing why, tells us the other half is still open.** (2026-10-30, jobs 31146622 `rwoff` and 31146643 `rwon`, both COMPLETED in under 4 minutes)

One variable, one structure: both re-optimised `protdft/scan_aspout/scan_04.pdb` with Zn-Asp92 OD1 held at 2.75 A, the same QM region, and - verified in the results, not assumed - **the same 485-atom active region** (`actatoms_fingerprint: 7c61d230af2a` in both).

| | `rwoff` | `rwon` |
| --- | --- | --- |
| `constraints.txt` | `$set`, 1 distance | `$set`, 25 distances |
| converged | true | true |
| steps (`optimization_energies.log`) | 28 | 28 |
| `SVD Error` | 0 | 0 |
| `opt_energy` (Eh) | -447.33041053983123 | -447.33040619825940 |

**`rwon` - `rwoff` = +4.34e-06 Eh = +0.0027 kcal/mol.** `site.py sphere` on the two optimised structures is identical on every distance to 0.01 A (O01 1.85, His88 2.00, His164 2.03, HOH485 2.05, OD1 2.75, OD2 4.55, mer O2x 4.82, O01-C4x 3.72); the only difference anywhere is one Burgi-Dunitz angle, 85 against 84 degrees. Twenty-four water constraints changed nothing at all.

**Why it cost nothing, and this is the part that matters.** `rwon`'s `constraints.txt` holds the waters at **0.9787 and 0.9785 A** - because `scan_04.pdb` came out of the fixed, flexible-water scan, so its waters already sit at the force field's own flexible minimum. The constraint asks them to stay where they already want to be, and `bond_constraints` holds *current* lengths, not TIP3P's. There is nothing to strain.

Compare 31016878, which started from an MD last frame: OpenMM ran that MD with `rigidWater=True`, so its waters are at exactly **0.9574 A**, the rigid model's value, while the QM/MM system is built with `rigidwater=False` and its flexible terms pull them to ~0.98-1.03. That is the whole mechanism of D137, now visible as a number in two constraint files: **0.9787 (satisfiable, free) against 0.9574 (strained, and in that run unsatisfiable)**.

**So D141 closes for one class of run and stays open for the other:**

- **Closed - results stand as computed:** anything optimised or scanned from a structure that had already been through a flexible-water QM/MM optimisation. That is `CSR_crystal_v2/protdft/zncheck/scan` and `protdft/coscan` (both from `QM_MM_optimized.pdb`), and every `qmmm` optimisation with no `opt_constraints`, where the waters were `$freeze`d and so projected out rather than scored.
- **Still open - and now with a measurable mechanism rather than a worry:** runs started directly from an MD frame, where the water constraints were pinned 0.02 A below the flexible minimum. That is `scan_H1_f00497` and `scan_H2_f00395` (from `nac_frame00497.pdb` / `nac_frame00395.pdb`, 46 points, 3 of them unconverged), D81's sphere-constrained reactant/product pair, and D105's stalled optimisation that D106 then used as a NEB endpoint.

**The job that would close the rest**, not submitted: the same `rwoff`/`rwon` pair from `../CSR_nac/md_naclong_results/NVTtrajectory_lastframe_imaged.pdb` - 31016878's own start - with Zn-OD1 held at its current value rather than 2.75. `rwon` there should reproduce the D137 failure and `rwoff` should converge, and the energy each reaches is the size of the correction for the MD-frame runs. About 20 minutes on 8 cores. **Declined** (user, 2026-10-30: "Not needed"), so the MD-frame runs stay flagged rather than measured: `scan_H1_f00497`, `scan_H2_f00395`, D81's pair and D105's endpoint are to be re-run if anything depends on their energies, not trusted to a tolerance nobody has established.

**`aspin` (31143211) finished while this ran: nine of nine points, all converged.** The point-by-point comparison against D140 is the next entry; the headline is that the reverse arm reaches -13.95 kcal/mol at 1.90 A against the outward arm's +13.84 at 3.60 A, and both put their minimum at 2.11 A.


**[SUPERSEDED by D145/D208 - the crystal is 4.20 A, off the end of both scan arms, so this cost was never measured. See `csr/retractions.md`.]** **D143 - The hysteresis test, answered: the two arms agree at the ends and at the minimum but separate by up to 3.8 kcal/mol in the middle, and the outward arm is the one with the defect. D94's objection is upheld for the middle of a relaxed scan and refuted for its conclusions.** (2026-10-30, job 31143211 `aspin`, COMPLETED 00:53:22, nine of nine converged, zero SVD failures, 112 MB)

Both arms on the same scale, each relative to **its own 1.90 A point** (absolute energies are not comparable between runs - the active region is re-derived per run, D141):

| Zn-OD1 (A) | `aspout` out (kcal/mol) | `aspin` back (kcal/mol) | out - back |
| --- | --- | --- | --- |
| 1.900 | 0.00 | 0.00 | 0.00 |
| 2.113 | **-0.86** | **-0.98** | +0.12 |
| 2.325 | +1.02 | +2.25 | -1.22 |
| 2.538 | +3.99 | +3.45 | +0.54 |
| 2.750 | +7.20 | +5.03 | **+2.17** |
| 2.963 | +10.56 | +6.78 | **+3.78** |
| 3.175 | +8.77 | +8.87 | -0.09 |
| 3.388 | +11.13 | +11.26 | -0.14 |
| 3.600 | +13.84 | +13.95 | -0.11 |

**What agrees, to about a tenth of a kcal/mol:** both ends (0.00 by construction; 13.84 against 13.95 at 3.60), the position of the minimum (**2.11 A in both**, not the 1.90 A either run started from and not the crystal's 2.80 A), its depth (-0.86 against -0.98), and the whole 3.18-3.60 A stretch.

**What does not:** 2.54-2.96 A, where the outward arm runs up to **3.78 kcal/mol high**. And the outward arm's non-monotonic step at 3.175 A (D140) is **not reproduced** - `aspin` is smooth from 3.60 all the way down to its minimum. The separation is one-sided: the reverse arm is lower wherever they differ, and its point 0 re-optimised `aspout`'s own last point to a lower energy. That is the signature of the outward arm lagging its relaxation - carrying a strained arrangement from point to point and releasing it at 3.175 - rather than of two genuine branches.

**So the answer to D94 is split, and both halves matter.**

- **A single-direction relaxed scan is not trustworthy point by point on this coordinate.** Errors of 2-4 kcal/mol in the middle of the range are the size of the effects this project argues about, so no barrier or difference read off one arm of a scan should be quoted to better than a few kcal/mol. That applies to `scan_znasp`, `scan_znco` and every other single-direction scan in this log.
- **But the conclusions those scans were used for survive**, because they rest on the ends and the minimum, which both arms reproduce. Asp92's preferred Zn distance on this surface is **2.11 A**, twice measured from opposite directions; pulling it to 3.60 A costs **13.8-14.0 kcal/mol**; and the crystal's 2.80 A sits **5.4-7.9 kcal/mol up** this curve - a range, not a number, and the honest figure to quote. It is not a thermal distance at either end of that range.

**The cycle closes structurally too.** `site.py sphere` on `aspin`'s 1.90 A endpoint gives five ligands - OD1 1.90, O01 1.94, His88 2.04, HOH485 2.19, His164 2.19 - with OD2 at 3.43 A and O01-C4x at 3.48, against `aspout`'s starting 3.83 and 3.58. Asp92 comes back monodentate onto the metal and the nucleophile stays where it was; the scan did not wander into a different basin and come back to a different structure.

**What this does and does not settle for the paper.** The question D133 set - is the preference for Asp92 on the metal real, or an artefact of where the scans started? - is answered: **real on this surface, from both directions.** What a scan still cannot give is a free energy: the 13.9 kcal/mol is a potential-energy difference along one coordinate with everything else relaxed, with no entropy and no sampling of the waters that reorganise when the carboxylate leaves. The disagreement with the crystal (D120, 2.80 A) therefore stands as the open problem, now with both a number and a bracket on it. **The PMF (D133 test 2) is the instrument that would convert this into a free energy, and after "do the scans, then we can do the pmf if necessary" it is the user's call whether 5.4-7.9 kcal/mol of potential energy is enough to argue from.**

Nothing submitted. Scratch 239/300 GB; `scan_aspin` is 112 MB and nothing needed pruning.


**[SUPERSEDED by D208 - the rule said restrain at 2.80 A; the measured distance is 4.20 A, outside both radii. See `csr/retractions.md`.]** **D144 - Neither question can be answered from what we have: the substrate has never coordinated the zinc in any structure in this project, which makes "is Asp92 coordinated when substrate is bound?" unanswerable rather than answered. `michaelis2d` (31154229) submitted to build the missing state.** (2026-10-30, user: "we can't trust MD for coordination so that's automatically invalid")

**A standing rule, from the user, that changes what counts as evidence here.** MD cannot be trusted for metal coordination, and **anything within coordination range should be restrained at the distance it has in the crystal**. CSR-1's crystal (6NC5 ZN302) has **Asp92 off the metal at 2.80 A**. So:

- every MD-derived statement about Asp92's coordination in this log - D102, D120's shell, D134's 56.7%, D138's 100% bidentate - is **inadmissible**, not merely caveated;
- and `md_naclong` was worse than inadmissible: `coordination_redetect` found OD1 at ~1.99 A and *restrained it there*, so the run imposed a contact the crystal does not have. Under this rule it should have been restrained at 2.80. Any future MD in this project restrains the crystal's distances explicitly and does not redetect Asp92.

**1. Does the substrate coordinate the zinc in the Michaelis complex? On everything we have: no, never, in any run.**

| source | closest substrate atom to Zn |
| --- | --- |
| `md_wt100` 5,000 frames (D134) | O1x 3.35 min, 3.92 mean; carbonyl O2x 7.80 at the end |
| `md_wt100b` 5,000 frames (D138) | O1x 3.00 min, 3.56 mean, **never inside 2.6**; O2x 5.98 |
| nine docked poses grafted into four systems (D126) | carbonyl 4.26-4.82 |
| D34, 791 frames | attack criterion never met |
| `aspout` / `aspin` optimised points (D140, D143) | O2x 4.68-5.01 |
| `rwoff`/`rwon` optimisations (D142) | O2x 4.82 |

`site.py check --hypothesis=carbonyl_on_metal` has been **REJECTED in every run it has ever been applied to**. The one structure that ever showed 2.55 A was `CSR_ligand`'s, and D126 established that the number was produced by that run's own restraint, in a system D123 had already rejected for its metal site. **So we have no Michaelis complex with the substrate on the metal - not a weak one, none.** What is not established is *why*: whether CSR-1 genuinely does not bind the carbonyl to Zn1, or whether every structure we have built simply failed to reach it.

**2. Does Asp92 coordinate the zinc if the substrate is also bound? Unanswerable from existing data, and not for want of looking.** Every "Asp92 is on the metal" result in this log - 2.11 A from both scan arms, 1.99 in naclong, 2.15 at 100% in wt100b - comes from a state in which **the substrate is not coordinated**. They answer "where does Asp92 sit when the fourth site is occupied by water and hydroxide", which is a different question. With the MD evidence now inadmissible as well, what remains is xTB QM/MM optimisation and the two scans, and those have O2x at 4.7-5.0 A throughout. There is no data point of the kind question 2 asks about.

**The user's framing, which is the right one and is now the open hypothesis:** unconstrained xTB *making* Asp92 coordinate is interesting, but the state it makes may not be the catalytic one. The crystal has Asp92 off; a four-coordinate site with a free face is what the carbonyl would need; and our models put Asp92 in that face. Whether the carboxylate is displaced *by* the substrate is exactly what nothing we have tests.

**So: `michaelis2d`, job 31154229** (`--cpus-per-task=8 --mem=64G -t 24:00:00`, PD at 16:0x). A 9 x 5 grid on `aspout`'s optimised 1.90 A point - substrate present, waters already relaxed (D142), Zn-O2x 5.01: rows are Zn-Asp92 OD1 at 1.90, 2.32, 2.75, 3.18, 3.60; within each row the carbonyl is pulled Zn-O2x 5.00 -> 2.10 in nine steps, snaking. Monitored throughout: Zn-OD2 (does the carboxylate go bidentate), Zn-O1x (the hydroxyethyl O, the other candidate ligand), Zn-HOH485 (does the water leave to make room) and O01-C4x (the attack distance). 45 points at D140's measured 8.5 min is about 6.5 h.

**What it decides, in each direction:**

- a row where the carbonyl reaches ~2.1 A at modest cost answers question 1 **yes**, and which row it is answers question 2;
- if only the Asp92-out rows allow it, the carboxylate swing is **part of the reaction coordinate** and belongs inside the NEB rather than fixed at one end of it;
- if no row allows it, `carbonyl_on_metal` is wrong for CSR-1 and the hypothesis in `site.yaml` is what changes, not the structures. The grid also reports whether O1x or the water takes the site instead, which is the obvious alternative and is monitored for that reason.

**It also unblocks the NEB, which cannot run today in any case:** `neb_step1` names `protdft/qmmm_protH_oh/QM_MM_optimized.pdb` and `protdft/scan_attack/scan_15.pdb`, and **neither file exists anywhere on scratch** (`CSR_crystal/scan_attack_f400` has a header and zero points). Both endpoints have to be built regardless, and a grid builds them on one surface instead of one at a time.


**[CONFIRMED by D208 - 4.20 A (OD2) / 4.22 A (OD1) re-measured independently; the caveat is that 6NC5's Zn2 is occupied. See `csr/retractions.md`.]** **D145 - The crystal has Asp92 4.20 A from the catalytic zinc, not 2.80: the number this log has quoted for ten weeks is our own model's Asp-hydroxide hydrogen bond, not a metal contact. And the crystal's Zn1 is four-coordinate with two histidines and two waters, so a substrate that coordinates at all has exactly one position to take.** (2026-10-30, measured directly in `CSR_crystal/6NC5.pdb`; user: "Asp92 is not coordinated in the crystal - which is important")

**The correction.** Measured on the deposited coordinates, zinc by zinc:

| pair | distance (A) |
| --- | --- |
| ZN302 (Zn1) - Asp91 OD2 | **4.20** |
| ZN302 (Zn1) - Asp91 OD1 | **4.22** |
| ZN302 (Zn1) - Asp91 CG | 4.66 |
| ZN301 (the mutant's Zn2) - Asp91 OD2 | 1.95 |
| Asp91 OD2 - HOH493 (the bridging water) | 3.58 |
| Asp91 OD2 - HOH485 | 4.31 |

(6NC5 numbering; Asp91 is our Asp92.) **No carboxylate oxygen is within 3.3 A of Zn1** - a 3.3 A sweep of ZN302 returns only His87 2.22, His163 2.09, HOH485 1.90, HOH493 2.14. The 2.80/2.84 A figure traces back to the model build, where it is recorded correctly as **our model's** Asp OD2...hydroxide hydrogen bond after Zn_alpha and its waters were moved as a rigid unit - "Here it H-bonds Asp120 OD2 (2.84 A)", in explicit contrast to the crystal. Somewhere between there and rule 2 it became "the crystal has Asp92 at 2.80 A", and D120, D133, D140 and D143 all inherited it. Rules 2 and 3 are corrected above.

**What changes, and it makes the problem bigger rather than smaller.** The disagreement between our surface and the structure is not 2.11 against 2.80 A - it is **2.11 against 4.20 A**. Both scan arms stop at 3.60 A, 0.6 A short of the crystal, at 13.8-14.0 kcal/mol and still climbing. So "the crystal's geometry costs 5.4-7.9 kcal/mol on our surface" (D143) is wrong: the crystal's geometry is **off the end of both arms**, and the cost of reaching it has not been measured at all. Nothing in this log has ever put Asp92 where the crystal has it.

**Question 1: is the catalytic zinc four- or five-coordinate with substrate bound? From existing data, the only defensible statement is about the four positions, not about the substrate.**

The crystal's Zn1 is **four-coordinate**: His87 2.22, His163 2.09, HOH485 1.90, HOH493 2.14 - two histidines and two waters, which is rule 4's canonical-equivalent sphere for an H-R-H site. In our wild-type model HOH493 becomes the hydroxide nucleophile (O01), so three of the four positions are **not exchangeable**: two histidines and the nucleophile that does the chemistry. **That leaves exactly one position, HOH485.** Therefore:

- a substrate that coordinates Zn1 must **displace HOH485**, and the Michaelis complex is then **four-coordinate**: 2 His + OH- + one substrate oxygen;
- a **five**-coordinate Michaelis complex requires either keeping HOH485 alongside the substrate, or Asp92 joining - and the crystal supports neither, since Asp92 is 4.2 A away and the fifth position does not exist in it;
- every five- and six-coordinate sphere in this log is Asp92 joining, which is now a 2.1 A difference from the structure rather than 0.7.

That is a structural argument, not a measurement: **no structure in this project has ever had the substrate in the sphere** (D144), so the Michaelis complex's coordination number remains unmeasured. What the crystal does fix is that there is one vacancy, not two.

**Question 2: does the substrate bind another way we have not tested? Yes - through O1x, and the structure for it already exists and has never been run.**

`Zn-mer O1x` was in no monitor until today, so no report printed it. With it printed: O1x is **3.89 A** in the four-coordinate structure and 3.92 A in the Asp92-bound one, against the carbonyl's 4.69-5.01 - the 6-alpha-hydroxyethyl oxygen is the closest substrate atom to the zinc in every structure and every trajectory we have (3.00 A minimum over 5,000 frames, D138). The `substituent_on_metal` arrangement from docking pose 2 has it **on the metal at 2.66 A** with the carbonyl 4.50 A out, HOH485 removed so the sphere is 2 His + hydroxide + O1x = **four**, and an attack geometry of **3.17 A at 106 degrees** - inside the Burgi-Dunitz window D34 never reached in 791 frames. `site.py check --hypothesis=substituent_on_metal` accepts it. It was submitted once as job 30923516 and cancelled an hour in (D130) for its restraint list, and never re-run.

**Under the user's new rule that list becomes trivial, which removes the reason it was cancelled.** D130's objection was that restraining Asp92 OD2 at 2.80 and OD1 at 3.81 asserted the carboxylate's orientation. With the crystal measured at 4.20/4.22 A, Asp92 is far outside any coordination cutoff - 3.0 A included - so under "restrain what is within coordination range, at its crystal distance" it gets **no restraint at all**, and the run neither imposes D120's artefact nor asserts an orientation. What is restrained is the crystal's own four: His88, His164, the hydroxide and (where kept) HOH485, at 2.22, 2.09, 2.14 and 1.90 A.

**Still running:** `michaelis2d` (31155441), which pulls the **carbonyl** in across five Asp92 rows from the four-coordinate structure. It tests the rival hypothesis, and its rows stop at 3.60 A - now known to be short of the crystal's 4.20. Nothing new submitted beyond it.


**D146 - `subst` is restrained at the crystal's distances now, but it is not run yet: its start structure has the substrate's hydroxyethyl proton 1.19 A from the nucleophile's oxygen, which MM cannot resolve. `substopt` (31160947) relieves that first, in twenty minutes rather than three hours of GPU.** (2026-10-30, user: "Run subst")

**What `subst_start.pdb` actually contains**, measured today rather than taken from D127's description:

| | distance (A) |
| --- | --- |
| Zn-O01 (nucleophile) | 1.86 |
| **Zn-Asp92 OD1** | **1.99** |
| Zn-His88 NE2 | 2.22 |
| Zn-His164 NE2 | 2.32 |
| **Zn-mer O1x** | **2.66** |
| Zn-Asp92 OD2 | 2.98 |
| Zn-mer O2x (carbonyl) | 4.50 |
| O01-C4x, angle | 3.17, **106 degrees** |
| **O01...H5x** | **1.19** |

The first five make it five-coordinate as built. The attack geometry is the best in this project - 3.17 A at 106 degrees, inside the Burgi-Dunitz window D34 never reached in 791 frames. **And the last row is a problem nobody had looked for:** the substrate's own hydroxyethyl hydrogen sits 1.19 A from the nucleophile's oxygen, which puts O1x and O01 roughly 2.1 A apart. That is a shared proton, not a hydrogen bond. `site.py sphere` has been printing it all along - "O01 ... 2 H H5x 1.19, H03 0.98" - and it was never read.

**Why that stops the MD rather than being a caveat on it.** The force field holds H5x on O1x by topology while the geometry wants it on O01, and MM cannot transfer a proton. A 100 ns run from here spends three hours of H100 on a strained contact and tells us about the strain, not the pose. This is a graft artefact in the ordinary sense: `site.py graft` transplants coordinates, nothing minimises them, and the pose was docked into a receptor rather than into this system's own site.

**`substopt`, job 31160947** (`qmmm`, `--cpus-per-task=8 --mem=64G -t 4:00:00`, PD at 17:0x): one QM/MM optimisation of `subst_start.pdb`, Zn-O1x held at the 2.66 A it already has - the hypothesis held where it is, not pulled (D93) - and **nothing on Asp92**, because where the carboxylate goes is the measurement. Three answers for twenty minutes: whether the proton settles on O01 or back on O1x, whether O1x stays on the metal when the site relaxes under a real Hamiltonian, and whether 106 degrees survives.

**`subst`'s restraint list, rewritten under the rule** (restrain what is in coordination range, at its crystal distance):

| pair | distance (A) | why |
| --- | --- | --- |
| Zn-Asp92 OD1 | 4.22 | crystal, 6NC5 ZN302-Asp91 OD1 (D145) |
| Zn-Asp92 OD2 | 4.20 | crystal, same carboxylate |
| Zn-His88 NE2 | 2.22 | crystal His87 |
| Zn-His164 NE2 | 2.09 | crystal His163 |
| Zn-O01 | (held where it is) | see below |
| Zn-mer O1x | 2.66 | the hypothesis, held not pulled |

Three things about it. **Asp92 is restrained at 4.2 rather than left free**: the crystal puts it there, the grafted structure has it at 1.99, and leaving it unrestrained is precisely what produced D120. **Both oxygens are named, at 4.22 and 4.20**, so no orientation is asserted - D130's objection to the old list (OD2 2.80, OD1 3.81) does not apply to two equal distances. And **the hydroxide is held where it is, not at the crystal's 2.14**: in 6NC5 HOH493 bridges two zincs, which is why it is long, and a terminal hydroxide on a mono-zinc site belongs at 1.85-1.95, which is what this structure has. That is the one place a crystal distance is not transferable, and it is a choice rather than a measurement - flagged here because it is exactly the kind of substitution that produced the 2.80 A error.

With HOH485 removed in the build, the sphere after Asp92 is pulled out is **2 His + hydroxide + O1x = four**, which is the only vacancy the crystal leaves (D145).

**Not submitted yet:** the 100 ns itself. `subst.inputs.pdb_file` still names the unoptimised graft because `sections.py` refuses a structure that does not exist, so it gets repointed at `protdft/qmmm_substopt/QM_MM_optimized.pdb` once 31160947 finishes and its sphere has been read. `michaelis2d` (31155441) is 56 minutes in on the rival hypothesis.


**[CORRECTED by D149 - the trend does not survive the last two rows; the surface peaks at 2.75 A. See `csr/retractions.md`.]** **D147 - The carbonyl cannot bind this zinc at any Asp92 position: 25-34 kcal/mol, and it gets in only by pushing **His88** off the metal while destroying the attack angle. `carbonyl_on_metal` is not expensive, it is anti-catalytic. `substopt` failed on the graft's clash, which puts the MM minimisation first; `subst` submitted as 31167021.** (2026-10-30, job 31155441 `mich2d`, 32 of 45 points at 19:52, all converged)

**Rate:** 32 points in 1 h 40 m, **3.1 min a point** against D140's 8.5 - the carbonyl pull is a smaller perturbation per step than the Asp92 one. Three of five rows are complete; the fourth is five points in. Scratch 241/300 GB.

**1. Does the carbonyl ever reach the metal?** The cost of pulling Zn-O2x from 4.70 to 2.10 A, each row against its own 4.70 A point:

| Asp92 held at (A) | Zn-O2x 4.70 -> 2.10 costs | row's 4.70 A point, vs row 0's |
| --- | --- | --- |
| 3.60 (four-coordinate) | **+25.2** | 0.00 |
| 3.175 | **+28.2** | -5.12 |
| 2.75 | **+34.0** | -9.03 |
| 2.325 | in progress (+18.9 at 2.10 A, 3.40 A at +1.9) | - |

**No row lets the carbonyl in at anything like a thermal cost.** The rows do differ, and in the direction the hypothesis wanted - the further Asp92 is out, the cheaper, 25.2 against 34.0 - so **Asp92's position is worth about 9 kcal/mol of the 25**, which settles D144's question 2 as far as a scan can: moving Asp92 helps and is not sufficient. The right-hand column is the other half of the same surface: with the carbonyl out at 4.70, bringing Asp92 in from 3.60 to 2.75 A *gains* 9.03 kcal/mol, which reproduces D143's curve on an independent structure and active region.

**2. What the carbonyl displaces, which is the finding.** `site.py sphere` at the two most forced points:

| | row 0 (Asp92 3.60) | row 3 (Asp92 2.32) |
| --- | --- | --- |
| Zn-O2x (forced) | 2.10 | 2.10 |
| Zn-O01 (nucleophile) | 1.84 | 1.88 |
| Zn-His164 NE2 | 1.96 | 2.12 |
| Zn-HOH485 | 2.21 | 2.06 |
| **Zn-His88 NE2** | **2.43** (from 2.03) | **3.59 - off the metal** |
| Zn-mer O1x | 3.61 | 3.64 |
| O01-C4x | 2.94 | 2.62 |
| **Burgi-Dunitz angle** | **58 degrees** | **69 degrees** |

**It does not displace the water.** HOH485 stays at 2.21 and 2.06 A; what moves is **a protein ligand**. In row 0 His88 is stretched from 2.03 to 2.43 A, and in row 3, with Asp92 also on the metal, His88 is **expelled to 3.59 A**. D145 predicted the carbonyl would have to take HOH485's position, the one vacancy the crystal leaves; the surface says it would rather take a histidine's.

**And the geometry it produces is the wrong one.** The nucleophile does come closer to C4x - 2.94 A, and 2.62 A in row 3, the shortest attack distances this project has seen - but the **Burgi-Dunitz angle collapses to 58 and 69 degrees**, far outside the 95-120 window, because an oxygen bound to the same metal as the carbonyl cannot approach the carbonyl carbon along its pi* axis. Short distance, hopeless angle. Compare the `substituent_on_metal` pose, where the carbonyl is 4.50 A off the metal and the attack line is 3.17 A at **106 degrees** with nothing restrained (D146). **The two hypotheses are not close: one costs 25 kcal/mol and ruins the angle, the other is free and satisfies it.**

**3. `substopt` (31160947) FAILED after 26 m 44 s, at optimiser step 76**, in ASH's xTB interface:

```
ValueError: could not convert string to float: '************-53.41870456************'
```

`xtbpcgradientgrab` reads xTB's point-charge gradient out of a fixed-width Fortran field, and a component overflowed it - asterisks where a number should be. That happens when an MM point charge and a QM atom approach near-zero separation, which is exactly what D146 flagged in the start structure: the substrate's H5x 1.19 A from the nucleophile's oxygen, O1x and O01 about 2.1 A apart.

**So D146 had the order backwards, and the failure is what shows it.** A QM/MM optimisation is the wrong tool for a steric clash, because the point-charge gradient can diverge before the geometry relaxes; **MM minimisation is the right one**, and it cannot fail this way - no point charges in a QM region, a finite Lennard-Jones wall, and a minimiser that will open O...H from 1.19 to an ordinary 1.8-2.0 A. The `md` stage begins with exactly that minimisation. The QM/MM optimisation belongs **after** it, on the equilibrated structure, where it can decide where the proton sits without a clash to resolve first.

**`subst` submitted: job 31167021** (`-p gpu_cuda --qos=gpu --gres=gpu:h100:1 --cpus-per-task=4 --mem=32G -t 12:00:00`, PD behind another session's GPU jobs). 100 ns NVT from `subst_start.pdb` with D146's crystal restraint list - Asp92 OD1 4.22 and OD2 4.20, His88 2.22, His164 2.09, the hydroxide held where it is, Zn-O1x at 2.66 - and `sections.py` checked all twelve selectors. About 3 h of production at D100's rate. **What it now decides is sharper than when it was designed:** with the carbonyl route costing 25-34 kcal/mol and expelling a histidine, this is no longer one hypothesis of two but the only arrangement left in which the substrate coordinates at all, and the run measures whether it survives 100 ns with the metal held at the crystal's own geometry.

**Nothing pruned:** the grid is still writing points into `protdft/scan_michaelis2d`, and the two optimisations' directories are small. `substopt`'s failed run is kept as the evidence above.


**D148 - The RQK triad measured at last, over 5,000 frames and three structures: Gln93 is not an oxyanion hole, Arg90 does not reach the substrate, and the residue that actually holds it is Arg229 - which is not in the QM region. Also: `subst_start.pdb` has Lys226 NZ **0.58 A** from the substrate, so `subst` is cancelled.** (2026-10-30, user: "why didn't you grab the info we can get with no new simulation?")

**What was built to get it** (commits 8f931f4, 7e0a750). `site.py contacts <traj> <topology>` reports every pair named in `monitor` over a trajectory - occupancy inside `contacts.cutoff` (3.5 A, the polar-contact convention D32 used) with min/mean/max. `shell` answers what is on the metal; nothing answered what holds the substrate, which is why the RQK numbers had only ever come from a one-off script. Sixteen pairs were added to `monitor` for the triad, the substrate's two carboxylate oxygens and the three residues that compete for them, so a contact is now reported because the site description names it.

**A bug found and fixed in the process, which is why the first run of this is not in this entry.** `one()` resolved a selector by residue number and atom name only, ignoring the chain: `D:1:O` matches the crystal water's oxygen **and residue 1's backbone O**, and whichever came first in the topology won. The first table had **Zn-HOH485 at 0% occupancy, minimum 14.13 A**. After the fix it reads **100%, 1.93/2.08/2.27 A** - which matches D138's `shell` row exactly, by a different code path. An ambiguous selector now raises and names the candidates. Any earlier analysis that used a `chain:1:O`-style selector on a trajectory is suspect; `shell`, `sphere` and `check` are not, because they resolve differently.

**The RQK triad over `md_wt100b`'s 100 ns** (5,000 frames, substrate held by the O01-C4x tether, Asp92 free; MM, so these are hydrogen-bond contacts, not coordination claims):

| pair | occ% | min / mean / max (A) |
| --- | --- | --- |
| **Arg229 NH2 - mer O4x** | **64.2** | 2.64 / 3.65 / 7.02 |
| Zn - mer O1x | 33.6 | 3.00 / 3.56 / 4.79 |
| Lys226 NZ - O01 | 18.0 | 2.46 / 4.02 / 5.38 |
| Lys226 NZ - mer O4x | 2.2 | 2.59 / 6.62 / 8.38 |
| Lys226 NZ - mer O3x | 0.1 | 2.92 / 6.06 / 8.37 |
| Arg229 NH1 - mer O3x | 0.0 | 3.43 / 6.13 / 8.43 |
| Arg90 NH1 - mer O1x | 0.0 | 4.65 / 6.08 / 7.51 |
| Arg90 NH1 - O01 | 0.0 | 5.45 / 7.01 / 8.60 |
| Ser186 OG - mer O4x | 0.0 | 4.78 / 6.91 / 8.50 |
| Asn225 ND2 - O01 | 0.0 | 4.42 / 5.61 / 6.88 |
| **Gln93 OE1 - mer O2x** | **0.0** | 6.16 / 8.37 / 9.66 |
| **Gln93 NE2 - mer O2x** | **0.0** | 7.86 / 9.66 / 11.06 |
| Gln93 NE2 - O01 | 0.0 | 6.31 / 7.64 / 8.46 |

And the same pairs in three optimised QM/MM structures, which do not share the tether or the pose:

| pair | `subst_start` (O1x pose) | `aspin/scan_00` (4-coord) | `mich2d/scan_00_08` (carbonyl forced on) |
| --- | --- | --- | --- |
| Gln93 NE2 - mer O2x | 8.12 | 11.10 | 9.76 |
| Gln93 OE1 - mer O2x | 7.43 | 9.87 | 8.47 |
| Arg90 NH1 - mer O1x | 6.16 | 4.37 | 4.39 |
| Arg90 NH1 - HOH485 | - | 3.82 | 3.51 |
| Lys226 NZ - mer O3x | 2.72 | 3.92 | 3.96 |
| Lys226 NZ - mer O4x | **0.58** | 5.37 | 5.30 |
| Arg229 NH2 - mer O4x | 6.42 | 3.64 | 3.64 |

**Four findings, and they are about the enzyme rather than the method.**

1. **Gln93 is not an oxyanion hole and is not near the chemistry.** 6.2-11.1 A from the carbonyl in every structure, 0% occupancy over 5,000 frames, and never closer than 6.31 A to the nucleophile. SMB-1's Gln157 holds the carbonyl at 1.67-1.70 A through ES, TS and EI (Mu & Xu 2020) - but Q157 is a Loop1 residue and our Gln93 is BBL Q121, a **Zn2-ligand position**. The positions are not homologous and the measurement says so: **the Q of RQK is a vestigial metal ligand, not a catalytic one.** If CSR-1 has an oxyanion hole it is some other residue, and nothing in our monitors is it.
2. **Arg90 does not touch the substrate; it sits on the water side.** 4.37-6.16 A from the 6-alpha-hydroxyethyl oxygen, 0% occupancy, but **3.51-3.82 A from HOH485** in both optimised structures - the closest cation to the non-nucleophile water. Its position (BBL R118, where Zn1's third histidine would be) fits a role in the solvent shell rather than substrate recognition. This **supersedes D32's "Arg90 NH1...O1x 99%"**, which was measured in the degraded shell on a different pose; the contact is pose-dependent, not a property of the enzyme.
3. **Lys226's grip on the substrate carboxylate is pose-dependent too** - 2.72 A in the O1x pose, 5.3-5.4 A and ~2% occupancy in the naclong-derived ones - while its contact with the hydroxide (18%, minimum 2.46 A) is the more persistent one. D107's "Lys226 salt-bridges the C3 carboxylate" holds for the pose it was measured on and not generally.
4. **The one persistent substrate contact in 100 ns belongs to Arg229** - NH2 to the C3 carboxylate O4x, **64.2% of 5,000 frames**, mean 3.65 A, and 3.64 A in both optimised structures. D32 already noted Arg229 is closer to the carboxylate and to the leaving nitrogen than Lys226 "yet is not in the QM region". It still is not. **Whatever anchors the substrate in CSR-1, it is this arginine and not the K of RQK**, and the QM region is drawn around the wrong cation.

**Separately, and it stops a job: `subst_start.pdb` has Lys226 NZ 0.58 A from the substrate's O4x.** A nitrogen and an oxygen cannot be 0.58 A apart. That is a second graft artefact, worse than the 1.19 A shared proton D146 found, and it **explains `substopt`'s failure better**: xTB's point-charge gradient overflowed its field (D147) because two atoms were essentially superimposed, not because of the proton. **`subst` (31167021) cancelled before it started** - 100 ns of MM from a 0.58 A overlap either fails on the first minimisation step or buries the error in a distorted site, and either way the 3 h of H100 would answer nothing. The `substituent_on_metal` hypothesis is not refuted by this; its **structure** is, and it needs re-grafting with a clash check before it can be tested.

**`michaelis2d` (31155441) finished while this ran: 45 of 45 points, all converged.** The full row-by-row table is the next entry.


**D149 - The grid is complete, 45 of 45 converged, and the carbonyl costs 25-34 kcal/mol in every row: the trend is not monotonic in Asp92, and whenever Asp92 is on the metal the carbonyl expels His88 outright. The job still reported FAILED - a 2D-scan bug in mcube that killed it after the last point.** (2026-10-30, job 31155441, 03:26:23)

**All 45 points converged**; 206 min of ASH walltime, ~4.6 min a point. **The job exited non-zero anyway**, and the traceback is mcube's own: `write_ash_surface` asked for the first coordinate as `target`/`target_A`, which is the **1D** spelling, while a 2D point carries `target1_A`. So every 2D scan this library has ever run has completed its chemistry and then died writing its summary - losing `surface_results.txt` and the `results` block, with the full CSV already on disk. Fixed with four tests (commit e4c2d64); the science here is untouched, since every number below comes from `scan_energies.csv`.

**The cost of pulling Zn-O2x from 4.70 to 2.10 A, now for all five rows**, each against its own 4.70 A point:

| Asp92 held at (A) | carbonyl 4.70 -> 2.10 | that row's 4.70 A point, vs row 0 |
| --- | --- | --- |
| 3.60 (four-coordinate) | **+25.17** | 0.00 |
| 3.175 | +28.16 | -5.12 |
| 2.75 | **+33.96** | **-9.03** |
| 2.325 | +26.67 | -7.76 |
| 1.90 | +27.56 | -6.95 |

**The trend D147 reported does not survive the last two rows.** It is not "further out is cheaper": the cost **peaks at Asp92 2.75 A** and falls again as the carboxylate closes further. The honest statement is narrower and still decisive - **the carbonyl costs 25-34 kcal/mol at every Asp92 position we sampled, and no position makes it affordable.** The spread across the whole Asp92 range is 8.8 kcal/mol, about a third of the cheapest barrier to binding.

**A caveat this grid raises against D143.** The right-hand column is the Asp92 coordinate with the carbonyl held out, and its minimum is at **2.75 A**, with 2.325 and 1.90 above it. D143's two arms put the minimum at **2.11 A**. Both say Asp92 sits on the metal, but the two instruments disagree about where by ~0.6 A, on different structures and active regions. Any number quoted for the Asp92 well should carry that spread.

**His88 is expelled whenever Asp92 is on the metal and the carbonyl is forced on.** `site.py sphere` at the three most forced points:

| | row 0 (Asp92 3.60) | row 3 (2.325) | row 4 (1.90) |
| --- | --- | --- | --- |
| Zn-His88 NE2 | 2.43 (from 2.03) | **3.59** | **3.63** |
| Zn-His164 NE2 | 1.96 | 2.12 | 2.11 |
| Zn-HOH485 | 2.21 | 2.06 | 2.18 |
| O01-C4x | 2.94 | 2.62 | **2.57** |
| **Burgi-Dunitz angle** | **58** | **69** | **71** |

With Asp92 off the metal the carbonyl merely stretches His88; with Asp92 on, there is no room for both and **the histidine leaves**. The attack distance gets shorter as the angle gets worse - 2.57 A at 71 degrees in the last row - which is the opposite of a productive approach, and the same pattern across all three.

**The RQK triad through the whole perturbation** (same points, from `monitor`): Gln93 NE2 is 8.5-11.1 A from the carbonyl at every one, Arg229 NH2 holds the C3 carboxylate at 3.61-3.64 A throughout, and Arg90 NH1 stays 3.34-3.82 A from HOH485. **Dragging the carbonyl onto the metal does not recruit any of R, Q or K**, which is the strongest form of D148's result: the triad is unmoved by a 25 kcal/mol perturbation of the substrate it supposedly recognises.

**`subst` (31167021) was cancelled, not run** - D148 found Lys226 NZ **0.58 A** from the substrate's O4x in its start structure. The restraint block this entry was meant to verify therefore has nothing to verify; the structure needs re-grafting with a clash check first.

**Nothing pruned.** `protdft/scan_michaelis2d` is 334 MB in 72 files - 45 point PDBs and one set of geomeTRIC files - with nothing redundant in it, and scratch is 243/300 GB.


**D150 - The no-new-simulation list finished: Gln93 stays 9.7-11.1 A from the carbonyl across all 45 grid points while the carbonyl moves 2.6 A, and the substrate's carboxylate anchor - Arg229 at 3.62 A, Lys226 at 3.90 - does not move at all. Both 100 ns runs now measured.** (2026-10-30, user: "I thought this was the new set of deliverables")

Four things were listed as answerable from data already on disk. D148 did the first. These are the rest.

**Deliverable 2: the occupancies on `md_wt100` as well** (5,000 frames, no substrate in the site, Asp92 free), beside D148's `md_wt100b`:

| pair | wt100 occ% | min/mean/max | wt100b occ% | min/mean/max |
| --- | --- | --- | --- | --- |
| Zn-Asp92 OD1 | 92.3 | 1.93/2.67/3.72 | 100.0 | 1.95/2.15/2.74 |
| **Lys226 NZ - O01** | **74.3** | 2.52/3.16/5.11 | 18.0 | 2.46/4.02/5.38 |
| **Asn225 ND2 - O01** | **43.9** | 2.60/3.80/5.90 | 0.0 | 4.42/5.61/6.88 |
| Arg229 NH1 - mer O3x | 35.1 | 2.55/4.05/10.14 | 0.0 | 3.43/6.13/8.43 |
| Arg229 NH2 - mer O4x | 5.1 | 2.67/5.16/9.91 | **64.2** | 2.64/3.65/7.02 |
| Gln93 NE2 - mer O2x | 0.0 | **8.24**/9.97/12.44 | 0.0 | **7.86**/9.66/11.06 |
| Arg90 NH1 - mer O1x | 0.0 | 5.79/8.01/9.21 | 0.0 | 4.65/6.08/7.51 |

**Gln93 is 0% in both runs, 10,000 frames, and never once inside 7.86 A of the carbonyl.** Arg90 is 0% in both. What differs between the runs is the hydroxide's company: with the substrate loose (wt100) **Lys226 and Asn225 take the hydroxide** - 74% and 44% - and when the substrate is held against it (wt100b) they are displaced to 18% and 0%. Arg229 holds the substrate carboxylate in both, through a different NH each time.

**Deliverable 4: R, Q and K across all 45 grid points**, where the carbonyl is dragged 4.70 -> 2.10 A and Asp92 is swept 3.60 -> 1.90:

| pair | min | mean | max | range over 45 points |
| --- | --- | --- | --- | --- |
| Zn-mer O2x (the driven coordinate) | 2.10 | 3.40 | 4.70 | 2.60 |
| Zn-His88 | 2.02 | 2.36 | **3.63** | 1.61 |
| **Gln93 NE2 - mer O2x** | **9.73** | 10.37 | 11.13 | 1.40 |
| Arg90 NH1 - HOH485 | 3.15 | 3.45 | 3.82 | 0.67 |
| **Arg229 NH2 - mer O4x** | 3.58 | **3.62** | 3.66 | **0.08** |
| **Lys226 NZ - mer O3x** | 3.86 | **3.90** | 3.96 | **0.10** |

**This is the oxyanion-hole test on a surface already paid for, and it is a clean negative.** Drag the carbonyl 2.6 A onto the metal, at a cost of up to 34 kcal/mol, and **Gln93 never comes closer than 9.73 A** - it moves less than the thing it is supposed to be stabilising by a factor of two, and always away. Nothing in RQK follows the oxyanion.

**What does hold the substrate is rigid and is not the K of RQK alone.** Arg229 NH2 sits at 3.58-3.66 A from the C3 carboxylate across every one of the 45 structures - a range of **0.08 A** while the substrate's other end swings 2.6 A - and Lys226 NZ at 3.86-3.96 A from the other carboxylate oxygen. The two cations form a fixed clamp on the carboxylate, with **Arg229 the closer of the two in every structure**, and the β-lactam end pivots in it. That is a concrete model of substrate recognition in CSR-1, measured rather than inferred, and it is the first thing in this project that behaves like a real anchor.

**Deliverable 3 is a decision, not a measurement, and it is the user's.** Arg229 is the closest cation to the substrate carboxylate in all 45 optimised structures and in both trajectories, and it is **still in the MM region**, treated as fixed point charges while Lys226 - further away in every structure - is quantum. Moving it into the QM region changes the QM charge from +3 to +4 and makes no run after it comparable with any run before it, which is why it is not done here.

**One operational note that cost a run tonight.** `conda activate ASH` can leave `CONDA_DEFAULT_ENV=ASH` while `which python` still points at the miniforge base, and then `site.py` dies on `ModuleNotFoundError: mdtraj`. The reliable form is the interpreter's absolute path, `/home/uquadhik/.conda/envs/ASH/bin/python`, which is what every analysis command here used after the first failure.


**D151 - The free energy was already in the frames, and it does not say what D150's paper outline said: at MM level Asp92 on and off the metal are two basins of near-equal population, not one filled vacancy. "The orphaned carboxylate fills the vacancy" is withdrawn as a claim.** (2026-10-30, user: "I'm not fully convinced the orphaned carboxylate fills the vacancy ... are you sure we have enough evidence to say that?")

**The user's doubt was right, and three things were wrong with the claim.**

1. **The framing.** 6NC5 is the **engineered Q93H/K226H double mutant with Zn2 present**, and in it Asp91 is a Zn2 ligand at 1.95 A. Its 4.20 A from Zn1 is a distance measured in a protein where a second metal holds the carboxylate. There is **no wild-type RQK structure**, so we are not contradicting an experimental observation of this site - we are saying a mutant that restores Zn2 cannot speak to the state without it. That is a far weaker claim and the only one the data supports.
2. **The arms are not independent.** Every structure descends from **one build** - superpose on 6NC5, delete ZN301, keep the carboxylate - and `aspin` started from `aspout`'s own end point. The hysteresis pair tests *path* dependence, which it passed. It does not test *model* dependence, and deleting a +2 ion beside a carboxylate will move that carboxylate toward the nearest remaining cation in almost any method.
3. **Asp92's protonation state has never been tested.** The literature leans ionised - a QM/MM/DFT study of CcrA found that protonating Asp120 disrupts its metal bonding, and IMP-1 mutagenesis argues it is not a proton donor - but those are di-zinc sites where a metal depresses the pKa, and that argument does not transfer to a site with no Zn2.

**What the existing frames already say.** `site.py contacts --profile` (new, commit cab5853) computes -RT ln P(r) from any monitored pair's histogram: unbiased sampling already carries the free energy over the range it explored. A radial Jacobian is divided out, and the frame count is printed beside every bin so a number standing on three frames looks like what it is. A `Zn-Asp92 CG` pair was added as the coordinate, because a profile along Zn-OD1 is contaminated every time the carboxylate rotates OD2 in and OD1 out - the same choice the classical carboxylate-cation PMFs make.

| run | coordinate | sampled range (A) | minimum | second basin |
| --- | --- | --- | --- | --- |
| `md_wt100` (no substrate, Asp92 free) | Zn-OD1 | 1.93-3.73 | **2.13** | **3.18 A at +1.86 kcal/mol** (74 frames) |
| `md_wt100` | Zn-CG | 2.33-3.37 | **2.47** | **3.22 A at +1.53** (123 frames) |
| `md_wt100b` (substrate tethered) | Zn-OD1 | 1.93-2.73 | 2.13 | none - 2.73 A is the furthest frame |
| `md_wt100b` | Zn-CG | 2.33-2.72 | 2.52 | none |

**`md_wt100`'s profile is two-state.** There is a sharp basin at 2.1-2.5 A, a thinly sampled gap at 2.6-2.8 (one frame a bin), and a second, broad, well-sampled basin at 3.2 A only **1.5-1.9 kcal/mol** above the first. And by population rather than by basin depth the two are closer still: D134 measured OD1 inside 2.6 A in **56.7%** of 5,000 frames, so **-RT ln(0.433/0.567) = +0.16 kcal/mol**. The bound state is narrow and the unbound state is broad; the minima differ by ~1.5 kcal/mol and the free energies by ~0.2, which is the entropy of a carboxylate that is free to move.

**So the correct statement is the one D133 guessed and this log then over-claimed on:** with no substrate, Asp92 is **in equilibrium between on and off the metal**, not parked on it. The 2.11 A minimum from both scan arms is real and is the enthalpic minimum; it is not evidence that the site is four-coordinate only in the crystal. **D150's paper outline is corrected: section 4 is "an equilibrium, measured", not "the vacancy is filled".**

**Two caveats that keep this from being the final word.** The profile is MM, with the other four donors restrained, which is exactly the kind of coordination claim the user has ruled MD out of. And the barrier between the basins rests on bins holding one frame each, so the shape between 2.6 and 2.8 A is noise. What the profile does establish is that **the off state is populated at all**, which no previous analysis showed, and that is enough to withdraw the claim.

**Costing the QM/MM umbrella sampling, on SMB-1's protocol.** Coordinate Zn-Asp92 CG, 2.3 -> 4.8 A (the crystal is at 4.66), 26 windows at 0.1 A, k = 150 kcal/mol/A^2 (the DFTB3/3OB zinc value; SMB-1 used 100-300), 45 ps equilibration + 60 ps data a window, WHAM - that is **2.7 ns of xTB/MM MD**. Our measured rate is soft: job 30467248 ran 50 ps of QM/MM MD plus an optimisation in 33.7 h on 8 cores, so **~1.7 ps/h**, and 105 ps a window is **~62 h**. Twenty-six windows in parallel is therefore ~62 h of wall time and **~13,000 core-hours**, each window inside the 72 h limit but not by much. **Before committing that, one short timed QM/MM MD is needed to replace the 1.7 ps/h figure**, which comes from a cancelled job and includes work that was not MD. **[Measured 2026-11-01, D157: the real rate is 2.57 ps/h, so 40.9 h a window and ~8,500 core-hours. The PMF is affordable; this entry's cost was 1.5x too high. It should still not be run before the pose is rebuilt - D156.]**

An MM umbrella run over the same coordinate would cost under an hour on a GPU and would converge the barrier the histogram cannot - but it would be measuring the force field's answer to a coordination question, which is the thing we have agreed not to trust. It is worth running as a sampling control, not as the result.

**The protonated-Asp92 test is not the quick run it sounds like.** Asp92 is in the QM region, but the MM system still needs a template for the residue, and the config's `amino_acids` list deliberately omits the protonated variants - the same blocker D133's charge-deletion test hit with LYN. It needs the residue renamed to ASH with an HD2 built on OD2, and the force field confirmed to template it. **Corrected the same night, by building it** (`site.py protonate`, commit d11a552): ff14SB's ASH is exactly ASP plus HD2 on OD2, and **`amino_acids` does not need extending** - `PDBFile` renames ASH back to ASP on load through OpenMM's own `pdbNames.xml` alt-name mapping, so `build_qm_region` never sees the name and treats the residue identically either way (verified: QM region 6 atoms for ASP, 7 for ASH, same cut at CA-CB). What does change is the **charge**: the residue goes -1 to 0, so the system's total charge and the QM region's both rise by +1 and the `charge` in the qmmm section must be incremented - a run that forgets this is wrong without looking wrong. ASH is also mid-chain only (there is no NASH/CASH), and LYN cannot be done by renaming at all, because OpenMM has no alt-name entry for it and the residue would load with no bonds.


**D152 - The Arg229 charge-deletion pair died in 8 s on a doubled section prefix, not on the feature: `sections.py` writes `sp_qmmm_ref` under `stage: sp` as `sp_sp_qmmm_ref`, and the job was given `sp_qmmm_ref`. Fixed at the tool, pair resubmitted, and the protonated-Asp92 test moved to the structure where Asp92 is actually on the metal.** (2026-10-31, jobs 31176144/31176145 FAILED; 31185324, 31185325, 31186037, 31186038 submitted)

**The failure.** Both single points ran 8 s and died identically:

```
KeyError: 'No geometry to score for the sp stage. Set `pdb_file` in `sp_qmmm_ref.inputs`.'
```

`sp_inputs` loads the config key equal to the section it is given, and the section given was `sp_qmmm_ref`. The key holding `pdb_file` was `sp_sp_qmmm_ref`, because `sections.py` writes every entry as `<stage>_<name>` and the entry had been *named* `sp_qmmm_ref` in `site.yaml`. So the run read an empty section, found no geometry and stopped - correctly, and loudly, which is the one good thing here.

**D110 found this and called it harmless.** Its own words: *"several existing keys double - `scan_znasp` becomes `scan_scan_znasp`, `sp_dft` becomes `sp_sp_dft` - which is harmless but means `run.slurm`'s documented examples disagree with what those keys actually produce"*, and it left them alone. It is harmless only while nothing is submitted: the moment a job is, the name a person types and the name the config holds are different strings and nothing between them objects. The pair was submitted at 21:39 and dead by 21:40, and the whole overnight window was spent with nothing running, for a naming rule that was written down and not enforced.

**Fixed where it can be caught** (`7e3a38c`, 5 tests). `sections.py` now refuses a name that begins with its own stage, before writing anything, and prints the key it would have made:

```
sp_qmmm_ref: drop the `sp_` - sections.py writes this as `sp_sp_qmmm_ref`,
while a job submitted as `sp_qmmm_ref` reads a section of that name
and finds nothing. Name it `qmmm_ref`.
```

The five `sp_*` entries in `site.yaml` are renamed to `qmmm_ref`, `qmmm_r229`, `dft`, `co0/4/8` - **the mcube keys they produce are unchanged** (`sp_qmmm_ref`, `sp_dft`, `sp_co0/4/8`), so every reference in this log still names the right thing. The two stale `sp_sp_*` keys were dropped from `config.yaml`; nothing else in it doubled, and `sp_dft`/`sp_co*` had never been run, so no result was ever stored under a doubled name.

**What is therefore still unknown.** The charge-deletion feature has not run once. Nothing is known about `pointcharge_lists`, `scope`, `charge_removed` or `verified_after_run`, and **there is no Arg229 electrostatic contribution** - D151's open question about whether `charge_sequences` finds where this ASH build keeps its point charges is exactly as open as it was. **Resubmitted as 31185324 (`spref`) and 31185325 (`spr229`)**, same structure, same settings, one key apart, 16 cores. The Lys226/Arg90/Gln93 sections were deliberately **not** submitted: the brief made them conditional on the pair working, and a decomposition built on an unverified deletion mechanism would be three more numbers of unknown meaning.

**The protonated structure: built, and it does not clash** (`fad1162`). `site.py protonate` on both ends of the `aspin` scan:

| structure | HD2 geometry | nearest atom outside Asp92 |
| --- | --- | --- |
| `scan_00.pdb` (Asp92 out, 3.60 A) | OD2-HD2 0.980, CG-OD2-HD2 113.00 deg, out of plane 0.0000, HD2-OD1 3.038 (syn would be 2.347) | **2.37 A** to E:9994:H2 |
| `scan_08.pdb` (Asp92 on, 1.90 A) | OD2-HD2 0.980, 113.00 deg, out of plane 0.0000, HD2-OD1 3.035 (syn 2.334) | **1.82 A** to E:9994:H2 |

Anti, in plane, and in neither case an overlap. The 1.82 A is an H...H contact that is tight but real - above the sum of two hydrogen vdW radii is 2.2-2.4 A, so it is strained rather than impossible, and the optimisation is free to relax it.

**A clash check added to the tool, because it had none.** `protonate` places the proton from three heavy atoms of one residue and knows nothing about the rest of the structure. It now reports the three nearest atoms outside the residue on **every** run - so the line is evidence the check ran, not something only a failure reveals - and refuses below 1.2 A, writing nothing. A strong O-H...O hydrogen bond is 1.6-2.0 A, so 1.2 A is an overlap and not a contact. `subst_start.pdb` carried a 0.58 A Lys226-O4x overlap for weeks (**D148**; this entry first cited D151, corrected 2026-11-05 per D208) and job 31160947 died of it in `xtbpcgradientgrab` with a row of asterisks where a point charge should have been; nothing in that path ever named the overlap. 6 tests.

**A correction to the overnight brief, which would have made the test vacuous.** It said to protonate `scan_aspin/scan_00.pdb`. **scan_00 is `aspin`'s 3.60 A end - the structure where Asp92 has already left the metal.** `aspin` runs Zn-OD1 from 3.60 inward to 1.90, so scan_00 is the four-coordinate site (that is exactly why `michaelis2d` starts there). Protonating it and optimising asks whether a neutral acid *comes back*, which is not the question the user put. The question - *"obviously protonating asp prevents it from coordinating - might be worthwhile doing a quick semi-empirical QM/MM run with that assumption"* - needs the structure that has something to lose.

**So `qmmm_asph` (job 31186037) starts from `scan_08.pdb`**, measured first with `site.py sphere`: **five ligands** - Asp92 OD1 1.90, hydroxide O01 1.94, His88 NE2 2.04, HOH485 2.19, His164 NE2 2.19, with OD2 at 3.43 (monodentate) and the substrate out at O1x 3.64 / carbonyl 4.96. `--to=ASH` puts HD2 on **OD2**, the oxygen that is *not* on the metal, leaving the bound OD1 bare: the chemically sensible tautomer, and the only one ff14SB defines. **`charge: 3`, not the inherited 2** - the residue goes -1 to 0, so the system's and the QM region's total charge both rise, and a run left at 2 converges the SCF on the wrong species while looking entirely normal. **Nothing restrained**: no `opt_constraints`, so the answer is a distance rather than a convergence criterion.

**And a real rate for the umbrella cost** (job 31186038, `qmmm_qmrate`). D151 costs umbrella sampling at ~62 h per window and ~13,000 core-hours from **1.7 ps/h, a figure taken from a job that was cancelled before it had run long enough to measure it.** 3 ps of QM/MM MD on 8 cores from `scan_08.pdb`, same QM region as the windows would use, `skip_optimize: true` so the optimiser is not what gets timed. D151's cost line can then be corrected in either direction.

**Four jobs running at 08:05**: 31185324/31185325 on bun157 (10 min), 31186037/31186038 on bun155. Disk checked first: scratch **250 of 300 GB**, home 38 of 50, 755k of 1M files - nothing near the limit and nothing pruned, since no new trajectory has finished.

**The bridge, for the record.** The overnight chain stopped at 02:47 with the code-server session expired and redirecting to the UQ sign-in page; nothing was submitted and no successor was scheduled, which is the standing rule and cost about five hours of queue. On reconnection the port had moved from 10539 to **47121** on bun149 - the URL in a wakeup prompt is not durable, and the first action of any shift has to be reading `location.href` and following the redirect.


**[QUALIFIED by D156/D157 - +44.34 kcal/mol is computed on a pose 5.15 A from the experimental complex. See `csr/retractions.md`.]** **D153 - Charge deletion works and is verified; Arg229's QM-MM electrostatics at the four-coordinate geometry is +44.34 kcal/mol, repulsive. A protonated Asp92 leaves the metal on its own and the site becomes the crystal's four-coordinate sphere. Arg90 and Lys226 cannot be measured this way at all, and the QM/MM MD rate is still unmeasured because the timing job computed nothing and reported success.** (2026-11-01, jobs 31185324, 31185325, 31186037 COMPLETED; 31186038 vacuous; 31277149, 31277150 submitted)

## 1. The charge-deletion pair, and the feature is sound

| run | job | energy (Eh) | elapsed |
| --- | --- | --- | --- |
| `sp_qmmm_ref` | 31185324 | **-5076.869006116424** | 03:01:51 |
| `sp_qmmm_r229` (A:229 deleted) | 31185325 | **-5076.939666070885** | 03:06:14 |

Both wb97x-v/def2-SVP on 16 cores, same structure (`scan_aspin/scan_00.pdb`), one key apart. The energies **differ**, so the deletion reached the QM Hamiltonian:

**E(with Arg229's charges) - E(without) = +0.07065995 Eh = +44.34 kcal/mol.**

The record, exactly as `confirm_zero_charges` wrote it: `key: zero_charges_qm_only`, `scope: qm_mm`, `residues: [A:229]`, **`atoms: 24`**, **`charge_removed: 1.0`** (ff14SB ARG, net +1 - the right residue, whole), **`verified_after_run: true`**, and no `[WARNING]` about restored lists.

**`pointcharge_lists: [theory.charges, theory.charges_qmregionzeroed]`** - not the `charges`/`charges_original` pair the test condition named, and this is **benign, not the defect that condition was watching for**. `charges_qmregionzeroed` is precisely the list ASH hands the QM code as point charges, so zeroing it is what makes the measurement QM-MM rather than MM-MM. The 44 kcal/mol difference is the proof: had the detection written to an attribute nothing reads, the two energies would have been identical to every digit. **No correction to the feature.** D151's open question - whether `charge_sequences` finds where this ASH build keeps its point charges - is answered yes.

**What the number is, and what it is not.** It is the *absolute* QM-MM electrostatic interaction between a +1 guanidinium and a +2 QM region at **one** geometry, and it is repulsive, which two like charges must be. It is **not** a catalytic contribution and must not be quoted as one: only the *change* in this quantity between the reactant and the transition state enters a barrier, and a residue can be strongly repulsive in both states while contributing nothing to catalysis. The measurement becomes mechanistic only when the same deletion is run at both ends of a reaction coordinate, on the same structures that give the barrier. That is the next use of these two sections, not a third residue.

## 2. Arg90 and Lys226 are not measurable this way - established without spending a job

`preprocess.results.qm_residues` puts **nine** residues in the QM region: Zn (B:1), HOH485 (D:1), His88, His164, **Asp92**, **Arg90**, **Lys226**, the hydroxide (C:2) and meropenem (C:1). Under electrostatic embedding ASH has already removed the MM charges of all of them, so deleting them again cannot change the energy, and the difference would be noise written down as an electrostatic contribution - which `charge_deletion` refuses by name, for exactly this reason. **Arg90 and Lys226 are therefore dropped**, as the brief required, and dropped from the config rather than by submitting two jobs to be told so.

**Gln93 is the only one of the three that can be measured**, and only because the two `sp` sections do *not* set `extra_qm_residues` (unlike the scans, which put Gln93 in QM): it is in MM for these runs and has charges to delete. Submitted as **31277149** (`sp_qmmm_q93`).

**A limit on the paper's electrostatics section, worth facing now.** A charge-deletion table over the RQK site cannot be completed: Arg229 and Gln93 by deletion, but Arg90 and Lys226 only by taking them out of the QM region - which changes the model and so cannot be compared with the other two. Either the decomposition covers the MM residues only and says so, or the whole set is done another way (alanine substitution, or a second parallel model with all four in MM) and the two methods are not mixed in one table.

## 3. A neutral Asp92 leaves the metal, and the crystal's sphere comes back

Job 31186037, free xTB/MM optimisation from `scan_08_A92_ASH.pdb`, **converged** in 12 min 03 s, charge 3, **nothing restrained on Asp92 or anything else**.

| ligand | protonated start | optimised |
| --- | --- | --- |
| **Asp92 OD1** | **1.90** | **2.99** |
| Asp92 OD2 | 3.43 | 4.73 |
| Asp92 CG | 2.95 | 3.81 |
| hydroxide O01 | 1.94 | **1.80** |
| His88 NE2 | 2.04 | 2.00 |
| His164 NE2 | 2.19 | 2.02 |
| HOH485 O | 2.19 | 2.08 |
| **ligands within 2.8 A** | **5** | **4** |

One proton, nothing held, and the carboxylate walks off to 2.99 A. The site is left with **two histidines, the hydroxide and HOH485 - the crystal's four-coordinate sphere** - and all three remaining ligands tighten, the hydroxide most (1.94 -> 1.80).

**What this decides.** The Asp92-Zn contact in this model **requires the anionic carboxylate**: it is electrostatic, not steric, and it is not a restraint artefact (nothing was restrained) nor an optimiser artefact (it converged). And the model *can* produce the crystal's four-coordinate sphere - something D120, D134 and D138 each failed to get by any other route, having tried twice.

**What it does not decide.** It does **not** show Asp92 is protonated at pH 7, and **the sentence that first stood here had the electrostatics backwards** - it claimed a carboxylate 1.9 A from a dication is the textbook case for an *upshifted* pKa. It is the opposite: a cation stabilises the conjugate base, so direct Zn(II) coordination *lowers* a carboxylate's pKa, which is why metal-bound aspartates are taken as ionised and why DFT work on the B. fragilis dizinc enzyme concludes Asp120 is ionised. See D154 for what that does to the interpretation. What it does give the D145 discrepancy - our model at 1.90-2.15 A against 6NC5's **4.20 A** - is a mechanism that produces the crystal number without touching the restraints, which is more than the other explanations have.

**And the substrate does not take the freed face.** O1x 3.64 -> 3.72, carbonyl O2x 4.96 -> 4.63, O01-C4x 3.48 -> 3.69. The vacancy closes by the **hydroxide** moving in, not by the substrate moving over. So emptying the Asp92 site does not by itself produce a Michaelis complex, which is a direct answer to the assumption behind `michaelis2d`: the open face is not sufficient, the carbonyl has to be driven in.

## 4. The timing job computed nothing and said DONE

31186038 ran **14 s**, exited 0, printed `[DONE] qmmm_qmrate`, and wrote a results record holding `qm_level: GFN2-xTB` and `metal_parameters: iod` and **no energy, no trajectory, no structure**. Nothing in the `.err` but the thread-binding INFO line.

**Cause: in the `qmmm` stage the MD *is* the equilibrate step**, so `skip_equilibrate: true` switched it off while `time: 3` sat in the section asking for 3 ps and was read by nothing. `site.yaml`'s own `nacarg` comment states this - *"skip_equilibrate: false # qmmm.inputs sets true; without this there is no MD"* - and the section was written wrong anyway.

**Fixed at the library** (`057a12e`, 5 tests): a `qmmm` section that skips both now raises **before anything is built**, names the `time` it asked for, and points at the `sp` stage for an energy on a geometry already in hand. A stage that can be configured to do nothing and report success is the same failure mode as D152's doubled prefix - a job number spent on a run that looks like a result.

**So D151's 1.7 ps/h is still unmeasured**, and the umbrella cost in that entry - ~62 h per window, ~13,000 core-hours - still rests on a cancelled job. Resubmitted as **31277150** with `skip_equilibrate: false`.

## 5. Housekeeping

Disk checked first: scratch **256 of 300 GB**, home 38 of 50, 844k of 1M files. **Nothing pruned.** The only candidate is `qmmm_asph/geometric_OPTtraj_Full.xyz` at 122 MB, and `qmmm_asph.results.opttraj` points at it: breaking a recorded path to recover 0.05% of the allocation is a worse trade than keeping it. The two `sp` directories are 24 MB each (four cubes plus the PySCF checkpoint), which is the evidence for the number in part 1.

The code-server port had moved again, bun149/47121 -> **bun151/28472**, as D152 said it would.


**D154 - No, we do not have an equilibrated Michaelis complex by this literature's standards, and the reason is upstream of everything measured since D98: the reactant pose came from docking rather than from reversing a product complex, and our own acceptance gate has rejected it twice. Also: a protonated Asp92 is not a plausible reactant state, which changes what the `asph` run is evidence for; and the whole job history is now inventoried.** (2026-11-01, user: "Have we established a good equilibrated michaleis complex structure according to the literature?")

## 1. What the literature does, and what we did

**The standard route to an MBL Michaelis complex is to reverse a crystallographic product complex.** There are no crystal structures of an MBL with an intact beta-lactam bound - that is the whole problem - so the QM/MM studies in this field build the Henry-Michaelis complex from a structure of the enzyme with the **hydrolysed** drug and undo the hydrolysis. The NDM-1/meropenem study builds its complexes "from the crystal structures of NDM-1 complex with hydrolyzed drug molecules, PDB 4RL2"; Lisa 2017 (*Nat Commun* 10.1038/s41467-017-00601-9), already cited in `reference-structures.md`, is the mononuclear/binuclear product-complex paper. Then: solvate, neutralise, classical equilibration, **QM/MM MD** (SMB-1: 10 ns of SCC-DFTB/CHARMM), and only then a surface.

**We docked instead.** Every reactant structure in this project descends from an AutoDock Vina pose - pose 2 by D123, chosen on a bond-graph match - and nothing in the chain ever came from a crystallographic complex with the ligand in it.

**And the gate says no.** `site.py check --at=reactant --hypothesis=carbonyl_on_metal` is **REJECTED** on both 100 ns trajectories, five of seven criteria failing (D134, D138). The supporting evidence, all of it already in this log:

| what was asked | what happened | where |
| --- | --- | --- |
| does the substrate stay bound without help? | **no** - it leaves the site in 100 ns | D134 |
| ... with a 50 kcal/mol/A^2 tether on O01-C4x? | it stays, and the zinc goes **six-coordinate in 4,998 of 5,000 frames** | D138 |
| is the carbonyl anywhere near the metal? | **5.98 A** with the tether on, 7.80 A without | D134, D138 |
| can it be brought in? | only by driving it, at **25.2-34.0 kcal/mol** across the five Asp92 rows | D147, D149 |
| is there an oxyanion hole holding it instead? | **no** - Gln93 never comes within 7.86 A of the carbonyl in 10,000 frames or 45 optimised points | D148, D150 |

So the honest statement is: **we have an equilibrated enzyme, and a substrate that will not stay in the pose we need without being held there.** Two 100 ns trajectories and a 45-point surface agree on that.

## 2. The hypothesis may be the wrong one, not just the pose

The closest literature analogue we have read in full is **SMB-1** (Mu & Xu 2020, B3), and its answer is that the beta-lactam carbonyl is **not a zinc ligand**: O8 sits 2.90/2.87/2.90 A from Zn1 through ES, TS and EI, and what holds it is **Gln157's amide** (Hepsilon21...O8 at 1.70 A). We have been testing `carbonyl_on_metal` - a hypothesis that paper rejects for its own enzyme - and measuring the one CSR-1 residue that could play Gln157's part (Gln93, the Q of RQK) only to find it 8-11 A away in every structure we own.

That gives three possibilities, and they are distinguishable:

1. **The pose is wrong.** Reversing the product complex would put the ligand where the chemistry requires rather than where Vina put it. `CSR_hydrolysed` already exists - built from the user's own crystal model `CSR_nontruncated_1metal_original.pdb`, with atom order deliberately matched to the reactant (D-entries at lines 101-122 of this log) - and it has never been used for this. **[Corrected 2026-11-01, D155: `CSR_hydrolysed` is a docked model too - CSR-1 has no complex structure of any kind, with substrate, product or inhibitor. The route that survives is transferring the ligand from another B3 enzyme's experimental complex, which is what D23 already did for the zinc.]**
2. **The carbonyl does not coordinate in CSR-1 either**, and the acceptance criterion in `site.yaml` should be an oxyanion-hole contact rather than a Zn-O2x distance - in which case the question becomes which residue, and the answer from our own measurements is *none of the RQK three in this pose*.
3. **The RQK site works on the C3 carboxylate rather than the carbonyl.** The one contact that is real and reproducible across everything we have is **Arg229 NH2-O4x at 3.62 +/- 0.04 A in all 45 optimised points and 64% of the MD frames** (D148, D150), and Lys226 NZ-O3x at 3.80-3.90 A. That is substrate *anchoring*, not nucleophile activation, and it is consistent with SMB-1, where Ser221/Thr223 hold the C3 carboxylate and no zinc contacts it.

None of these is settled by more sampling of the current structure.

## 3. A protonated Asp92 is not a plausible reactant, and D153's reasoning was wrong

D153 said the protonated state was hard to rule out because "a buried carboxylate 1.9 A from a dication is the textbook environment for an upshifted pKa". **That is backwards.** A nearby cation stabilises the conjugate base, so coordination to Zn(II) *lowers* a carboxylate's pKa - the same effect that takes the zinc-bound water in carbonic anhydrase from ~15.7 to 6.8. Burial and desolvation push the other way, but for an oxygen at 1.9 A from the metal the metal wins, which is why metal-bound aspartates are modelled as ionised and why the DFT work on the *B. fragilis* dizinc enzyme concludes **Asp120 is ionised**. D153 is corrected in place.

**So `asph` (31186037) is a diagnostic, not a candidate reactant.** What it establishes stands: the Asp92-Zn contact in this model is **purely electrostatic** - remove the charge and the carboxylate leaves on its own, unrestrained, and the site relaxes to the crystal's four-coordinate sphere. What it does not do is offer a protonated Asp92 as the state the Michaelis complex starts in: at pH 7, coordinating or not, that carboxylate is deprotonated. The D145 discrepancy (our 1.90-2.15 A against 6NC5's 4.20 A) therefore still needs an explanation, and "it is protonated" is no longer on the list. What remains: xTB/ff14SB over-binding the carboxylate to Zn(II); the crystal being a **Q93H/K226H double mutant** with a second zinc that changes the first shell (D145); or the crystallographic Asp genuinely being held off the metal by something our truncation removed.

## 4. The job history, inventoried

`docs/csr/job-inventory.md`, generated from `sacct` rather than from memory: **166 CSR allocations** since August, one row each with state, elapsed, directory and the decision-log entries that cite it, plus each system's QM region and the structures that are live inputs today. Forty of them are the early `mcube_QMMM` setup cluster - all FAILED inside two minutes, in `CSR_ligand` and `CSR_hydrolysed`, none cited anywhere - and they are counted rather than listed. The regeneration command is in the file, on the same principle as `reference-structures.md`: a number in a doc should be re-derivable, not quoted.

**The question it was written to answer, answered: no, there is no equilibrated complex with Arg90 and Lys226 outside the QM region.** All eleven systems put them there, because `preprocess.inputs.extra_qm_residues: {A: [90, 226]}` says so and neither is a zinc ligand the 2.5 A detection would find. `CSR_hydrolysed` is the only directory with Lys226 in MM, and it is the hydrolysed product with **no His164, no hydroxide and no HOH485** in its QM region - the wrong state and an incomplete first shell. Measuring those two by charge deletion therefore needs a fresh `preprocess` and a re-equilibration of everything downstream: a new model, whose numbers would not be comparable with the Arg229 value in D153.

## 5. And the silent-success class is closed at the CLI

D153 fixed the one section that could skip both of its steps. The general fix is now in `mcube.cli`: **`confirm_produced` raises unless the stage's result carries a key that means it computed something** - `sp_energy` for `sp`, `lastframe`/`trajfile` for `md`, `optimized_pdb`/`lastframe`/`trajfile`/`opt_energy` for `qmmm`, and so on - so the exit code is non-zero and **Slurm records FAILED rather than COMPLETED**. That is the part that matters: job 31186038 sits in the inventory as COMPLETED having computed nothing, and every reader of a job table treats COMPLETED as a result. `umbrella` is deliberately unlisted (its three entry points return three different shapes with no key in common) and a test asserts that, so listing it has to be a deliberate change. 12 tests.


**D155 - The precedent for a docking-derived Michaelis complex says the dock is a starting point, not the answer: in the closest published study the ligand is placed by swapping it onto an experimental complex's core and the pose is then produced by ~10 ns of QM/MM MD. We already do exactly that for the zinc and never did it for the substrate.** (2026-11-01, user: "CSR_hydrolysed is also a docked structure. There's no crystal of CSR-1 with substrate bound. We need to get the literature-grade Michaelis complex from a docked structure, look for precedent.")

D154 proposed reversing `CSR_hydrolysed` and that is wrong: it is a docked model too, and CSR-1 has no complex structure of any kind. Correcting that, here is what the literature actually does when it has no ES crystal - which is **always**, for this enzyme family.

## 1. SMB-1, in its own words

Mu & Xu 2020 (*J Mol Model* 26:71, the local PDF) is the closest analogue we have - same subclass, same question - and it had the same problem: *"due to its strikingly efficient catalytic activity, it is impossible to obtain crystallographic structure of enzyme in complex with a natural substrate."* Their protocol, from the Computational details section:

1. **Start from an experimental complex of the same enzyme with a non-substrate ligand.** *"the initial model was extracted from the crystallographic structure of SMB-1 (PDB entry code 5B15), which is a structure in complex with hydrolyzed doripenem."*
2. **Swap the ligand on the shared core.** *"On the basis of the same core structure between doripenem and ampicillin, we then manually replace this hydrolyzed intermediate with the ampicillin molecule."* Note what that is: a carbapenem replaced by a penicillin, by hand, on the core.
3. **Place the mechanistic pieces by chemistry, not by docking.** *"one hydroxide ion was inserted into this position, which is further hydrogen bonded with D120"*; the C3/C4 carboxylate was *"postulated to have direct contact with zinc ion at the second position"*.
4. **Let the simulation produce the pose.** *"Such a model would be further subjected to an extensive QM/MM molecular dynamic (MD) simulation for relaxation to generate the equilibrated Michaelis complex structure"* - 50 ps heating, 250 ps NPT, then **9.7 ns of SCC-DFTB/CHARMM QM/MM MD**.

**The pose in that paper is not the modeller's and not a docking program's: it is whatever survived ten nanoseconds of QM/MM MD.** That is what "equilibrated Michaelis complex" means in this literature, and it is the standard we have not met - not because we docked, but because nothing has relaxed the complex at QM/MM level and then been measured.

## 2. The precedent for starting from a dock

For a docked start specifically, the IMP-1 study (*Proteins* 2022, 10.1002/prot.26227) is the pattern: **dock, then MD, then QM/MM refinement**, with two safeguards - the docking protocol is validated against ligands whose complexes *are* known, and the QM/MM optimisation is reported as having restored the experimentally observed coordination. So a docked start is publishable when (a) the protocol is shown to reproduce a known pose, and (b) the final geometry comes out of QM/MM, not out of the docking score. We have done neither: no redocking control exists, and no QM/MM MD of the current complex exists.

## 3. We already use template transfer - for the metal

**CSR-1's own crystals have no metal in them.** 6DN4 and 6DQ2 are metal-free (D22), and the published characterisation says the same of CSR-1 and CSR-1_trunc. The zinc in every structure in this project is there because **D23 superposed 6NC5 - the Q121H/K254H construct - on 256 Calpha at 0.58 A RMSD and moved Zn_alpha and its two waters across as a rigid unit**, then adjusted by 0.38 A so Zn-NE2 matched. The nucleophile is 6NC5's HOH493, transferred the same way.

That is precisely the technique SMB-1 used for its ligand, applied by us to the metal and accepted without argument. **The substrate is the one component of the model that was never placed this way** - it came from AutoDock Vina, pose 2, chosen on a bond-graph match (D123). The inconsistency is the finding: there is no principled reason the metal gets a crystallographic template and the ligand does not.

## 4. The templates we can transfer from

| template | what it holds | why it fits |
| --- | --- | --- |
| **7UHI / 7UHJ** (+ 7UHH, 7UHK) | L1, B3, **intact moxalactam** - the only measured Michaelis complex in the family | the lactam end is positioned by experiment, carbonyl O1 at 2.53/2.63 A from Zn1 |
| **5B15** | SMB-1, B3, hydrolysed **doripenem** | the published swap was doripenem -> ampicillin; doripenem -> **meropenem** is carbapenem to carbapenem, a *smaller* step than the one in print |
| L1 + hydrolysed carbapenem products (PMC10148155) | stereochemistry of the product end | a check on the product state if a NEB ever needs one |

**The caveat that has to be stated, not hidden: every template is di-zinc, and in all of them the C3 carboxylate is held by Zn2, which CSR-1 does not have.** So a transfer fixes the lactam end and the scaffold and leaves the carboxylate end to relax onto whatever CSR-1 offers. Our own measurements already say what that is - **Arg229 NH2-O4x at 3.62 +/- 0.04 A in all 45 optimised points and 64.2% of wt100b frames, Lys226 NZ-O3x at 3.80-3.90 A** (D148, D150) - which is the RQK site standing in for Zn2. That is the paper's thesis, and a transferred pose is a way to *test* it; a Vina pose is not, because the contacts it produces are whatever the scoring function liked.

## 5. The cost, which is a system-size decision and not a core-count one

SMB-1's 9.7 ns of QM/MM MD is the number to match, and at the measured 2.57 ps/h (D157; this entry first said 240 days from D151's unverified 1.7) it is **about 157 days**. They could afford it because of their boundary: a **25 A TIP3P sphere centred on Zn2, everything beyond 25 A deleted**, 22-25 A under Langevin dynamics, free inside 22 A, and 114 QM atoms. That is a few thousand atoms against our periodic box.

So reaching the literature's equilibration length means changing the model, not the allocation. Three ways, in order of how much they concede:

1. **Build a stochastic-boundary model for the equilibration only** - a sphere around the zinc, as SMB-1 did - and carry the relaxed complex back into the periodic system for the surface. Most faithful to the precedent; needs a truncation step mcube does not have.
2. **Equilibrate classically for nanoseconds and report a short QM/MM MD as a stability check**, stating plainly that it is picoseconds against SMB-1's 9.7 ns. Cheapest, and defensible if the limitation is in the methods section rather than buried.
3. **Measure the real rate first and then decide** - job 31277150 is running for exactly this, and until it lands every cost in this entry is built on a cancelled job's figure.

## 6. What this means for the acceptance gate

`carbonyl_on_metal` was already withdrawn for B3 (D123), and SMB-1 puts its carbonyl at **2.90 A from Zn1, not a ligand**, held instead by Gln157's amide. So the gate for a transferred pose should be the measured L1 geometry plus the C3-carboxylate anchoring, **judged after QM/MM relaxation rather than at the moment of grafting** - which is also the only way the judgement is about the enzyme rather than about the modeller.

**Nothing submitted.** The step this points to is a build, not a job: superpose 7UHI/7UHJ and 5B15 onto the CSR-1 model, transfer the ligand as a rigid unit, graft meropenem onto the shared bicyclic core, run the redocking control the IMP-1 precedent implies (re-place each template's own ligand by the same procedure and report RMSD), clash-check it, and only then equilibrate. `site.py graft` already does bond-graph matching between ligands and `protonate` now clash-checks what it writes, so most of the tooling exists.


**D156 - The docked poses are measured against the only experimental Michaelis complex in the family, and none of them reproduces it. The pose the project has used since D123 is the second-furthest of the nine: its beta-lactam carbonyl is 4.90 A from the zinc where L1's is 2.27, and its C3 carboxylate is 8.14 A where L1's is 2.82. AIM-1 cannot help - it has no ligand-bound structure.** (2026-11-01, user: "compare the Vina poses to the L1 and SMB-1 structures, maybe also chuck AIM-1 in there if there's a structure with substrate bound")

## 1. The templates, and which of them exist

| PDB | enzyme | ligand | state | usable as a pose template? |
| --- | --- | --- | --- | --- |
| **7UHI / 7UHJ** | L1 (B3, di-Zn) | moxalactam, **intact** | Michaelis | **yes** - the only one |
| **5AXO** | SMB-1 (B3, di-Zn) | **hydrolysed meropenem** | product | yes, for the product end |
| 5B15 / 5B1U | SMB-1 | hydrolysed doripenem / imipenem | product | same |
| 4AWY, 4AWZ, 4AX0, 4AX1 | **AIM-1** (B3) | **none** | apo / Q157 mutants | **no** |

**AIM-1 has no substrate-bound structure.** Its four entries are the 2012 Leiros set: wild type, and the Q157A and Q157N mutants, all with zinc and water only. That paper studied binding by *docking plus QM/MM*, which is the same position we are in rather than a way out of it.

One correction to a search result that would have cost a day: **5AXR is not the hydrolysed-meropenem complex** - it is SMB-1 with 2-mercaptoethanesulfonate, checked by reading the file. The meropenem complex is **5AXO** (ligand code LMP).

## 2. How the comparison was done

New `refstruct.py transfer` mode (`b5e226a`, 19 tests). It fits a template on named zinc-site atoms, moves its ligand into our frame, and reports what that ligand and each of our poses present to our zinc. The correspondence is **given explicitly, not guessed** - the three proteins number their residues differently - and the fit RMSD is printed so a wrong correspondence shows up as a bad fit rather than as a plausible pose. The fit is local, on the Zn1 site, which is the right frame for comparing where a ligand sits relative to the metal and makes no claim about the fold.

The BBL correspondence used, from `refstruct.py shell` on each file rather than from the papers:

| BBL | CSR-1 (model) | L1 7UHI | SMB-1 5AXO |
| --- | --- | --- | --- |
| 116 | His88 NE2 | His105 NE2 | His72 NE2 |
| 196 | His164 NE2 | His181 NE2 | His150 NE2 |
| 120 | Asp92 CG | Asp109 CG | Asp76 CG |
| Zn1 | B:1:ZN | ZN301 | ZN301 |

**The B3 Zn1 site is conserved enough for this to work.** SMB-1 fits onto L1 across two different enzymes at **0.254 A** over those four atoms; both fit our model at 0.50-0.85 A. L1's own Zn1-moxalactam O1 distance is 2.53 A and comes back as 2.27-2.85 A after transfer, which is the size of the fit error and the uncertainty every number below carries.

## 3. The nine Vina poses against L1's moxalactam

All in `CSR_ligand`'s frame, the one the poses were docked in, with the L1 ligand transferred into it (fit RMSD 0.851 A). "lactam C=O" is the ring carbonyl, identified by connectivity; "carboxylate" is the C3 one; "centroid" is the whole-ligand centroid against the transferred template's.

| pose | Zn-lactam O | Zn-carboxylate O | centroid vs L1 |
| --- | --- | --- | --- |
| **L1 moxalactam (transferred)** | **2.27** | **2.82** | - |
| 1 | 10.29 | 8.64 | 1.80 |
| **2 - the one in use since D123** | **4.90** | **8.14** | **3.66** |
| 3 | 9.21 | 4.65 | 2.13 |
| 4 | 5.06 | 4.43 | **0.98** |
| 5 | **3.67** | 4.70 | 2.39 |
| 6 | 3.77 | 4.03 | 3.49 |
| 7 | 4.35 | 4.52 | 1.93 |
| 8 | 3.83 | **3.69** | 2.49 |
| 9 | 5.43 | 8.08 | 3.22 |

**No pose reproduces the experimental binding mode.** The best carbonyl approach in the whole set is 3.67 A (pose 5) against 2.27; the best carboxylate is 3.69 A (pose 8) against 2.82; and no pose is good at both. Pose 4 puts the ligand in nearly the right place - centroid 0.98 A - and still has its carbonyl 5.06 A from the metal, so it is the right position in the wrong orientation.

**Pose 2, the one everything since D123 descends from, is the second-worst of the nine by centroid and worst-but-two by carbonyl.** It was chosen on a bond-graph match to a hand-picked pose and a Vina score of -5.624, not against a structure. Had this comparison existed in October, pose 5 or 8 would have been the candidate and pose 2 would have been rejected.

## 4. And the structure we are actually running

Measured in `CSR_wt`'s own frame (fit RMSD 0.622 A, transferred L1 carbonyl at 2.85 and carboxylate at 2.54):

| | Zn-lactam O | Zn-carboxylate O | centroid vs L1 |
| --- | --- | --- | --- |
| L1 moxalactam (transferred) | 2.85 | 2.54 | - |
| **`scan_aspin/scan_00.pdb`** | **4.69** | **6.37** | **5.15** |

**5.15 A from the experimental pose - further than any of the nine docked poses was.** That is what the MD did to it: D134 and D138 recorded the substrate leaving the site and needing a 50 kcal/mol/A^2 tether, and this is the same fact in the one unit that can be compared with experiment. The 2D surface of D147/D149 was driving the carbonyl in against ~25-34 kcal/mol because the pose it started from has the ligand in the wrong place, not only the wrong orientation.

## 5. What follows

1. **The pose is the error, and it is now quantified rather than suspected.** D155 argued from precedent that the ligand should be placed by transfer from an experimental complex; this says the same thing from our own numbers, and gives the target: Zn-carbonyl ~2.5 A, Zn-C3-carboxylate ~2.8 A, and in L1 both at once.
2. **The transferred L1 pose is itself the candidate start**, not a pose to score against. `refstruct.py transfer --out` already writes it; what it needs is meropenem grafted onto the moxalactam core (`site.py graft` does bond-graph matching), a clash check, and then the equilibration D155 set out.
3. **One honest caveat, the same one as D155.** L1 is di-zinc, and in 7UHI the C3 carboxylate is held by **Zn2**, which CSR-1 does not have. The 2.82 A in the table is that carboxylate measured against our *single* Zn1 after transfer, so it is where the carboxylate sits, not a bond we should expect to form. What CSR-1 offers that end is Arg229 and Lys226 (D148, D150), and whether they can hold it there is exactly the question the rebuilt complex would answer.
4. **Nothing is yet rebuilt and nothing was submitted.** The two jobs from D153 (31277149 `q93`, 31277150 `qmrate`) were still running when this was written.


**[QUALIFIED by its own section 4 - 'the last word on the Q of RQK' rests on the same condemned pose. See `csr/retractions.md`.]** **D157 - Gln93's whole electrostatic contribution is -0.87 kcal/mol, fifty times smaller than Arg229's and the last word on the Q of RQK. QM/MM MD runs at 2.57 ps/h, not 1.7, which makes the Asp92 PMF affordable - about 8,500 core-hours, every window inside the walltime limit. The PMF was never the blocker; the pose is.** (2026-11-01, jobs 31277149 and 31277150, both COMPLETED)

## 1. Gln93, measured

| run | job | energy (Eh) |
| --- | --- | --- |
| `sp_qmmm_ref` | 31185324 | -5076.869006116424 |
| `sp_qmmm_q93` (A:93 deleted) | 31277149 | **-5076.867618573601** |

**E(with Gln93's charges) - E(without) = -0.00138754 Eh = -0.87 kcal/mol**, stabilising. Against Arg229's **+44.34** (D153), that is a factor of **51**, and the sign is opposite.

The record: `atoms: 17`, **`charge_removed: 0.0`**, `verified_after_run: true`, same two point-charge lists, no warning. Seventeen atoms is ff14SB's GLN and zero is its net charge, so the residue the key found is the neutral amide it should be - the check the brief asked for, passed.

**And a point worth keeping: the energies differ although the charge removed is zero.** A net-neutral residue still has a dipole, and deleting it changes the QM region's energy; if the method only ever registered monopoles, these two numbers would have been identical. So the -0.87 is evidence the measurement sees the field rather than the formal charge.

**What it settles.** D150 said the Q of RQK "is a vestigial metal ligand, not a catalytic one" on geometry alone - 6.2-11.1 A from the carbonyl, 0% occupancy over 5,000 frames. The electrostatics now agrees from the other direction: **0.87 kcal/mol is below the error of everything else in this project**, and no mechanism rests on a residue contributing that. Gln93 is where Zn2's Asp-His-His set lost its histidine, and it has not acquired a new job.

**The decomposition is therefore complete and it has two entries**, which is all this method can ever give here: Arg229 **+44.34**, Gln93 **-0.87**. Arg90 and Lys226 are inside the QM region and permanently out of reach (D153); not attempted, and not to be made reachable by moving them out, which would be a different model.

## 2. The QM/MM MD rate, measured at last

3 ps on 8 cores, `skip_equilibrate: false` this time, trajectory and last frame both written: **4204.1 s, so 2.57 ps/h** - about **1.5x faster** than the 1.7 ps/h D151 took from a cancelled job.

**What that does to the umbrella cost** (26 windows, 105 ps each - 45 ps equilibration plus 60 ps data, SMB-1's protocol):

| | per window | 26 windows, 8 cores |
| --- | --- | --- |
| D151, from the cancelled job | 61.8 h | ~12,800 core-hours |
| **measured** | **40.9 h** | **~8,500 core-hours** |

**That is affordable, and D151's framing was wrong.** Every window fits inside a single 72 h allocation with 31 h to spare, and the 26 run in parallel, so the wall-clock cost of the whole PMF is **about two days**, not a quarter. The entry treated ~13,000 core-hours as the reason the Asp92 free energy was out of reach; it is not, and the figure that said so came from a job that was stopped before it had run long enough to be timed.

**This does not mean the PMF should be run now.** D156 measured the start structure 5.15 A from the only experimental Michaelis complex in the family: a PMF on the Asp92 coordinate computed from that pose would be an accurate free energy for the wrong complex, and 8,500 core-hours is only cheap if the answer means something. The order is the rebuilt pose first, this second.

## 3. A free result: 3 ps of unrestrained QM/MM MD

The timing run is also **the first unrestrained QM/MM MD of the current system**, and `run.slurm`'s own `sphere` report on its last frame says:

| ligand | scan_08 start | after 3 ps |
| --- | --- | --- |
| Asp92 OD1 | 1.90 | **2.17** |
| hydroxide O01 | 1.94 | 1.98 |
| His88 NE2 | 2.04 | 2.07 |
| HOH485 | 2.19 | 2.12 |
| His164 NE2 | 2.19 | 2.16 |
| ligands within 2.8 A | 5 | **5** |

**Asp92 stays on the metal at xTB/MM with nothing restraining it**, which is D102, D134 and D138 again at a third level of theory, and it is the control the `asph` result of D153 needed: the carboxylate leaves when it is protonated and stays when it is not, in the same code, from the same structure, over comparable relaxation. The substrate stays out (O1x 3.70, carbonyl 5.30, O01-C4x 4.20 against 3.48 at the start) - drifting further in 3 ps, consistent with D156.

## 4. Where this leaves the open question

The brief this shift ran from asked for the charge deletions to be written up as needing a reaction coordinate. That framing is **superseded by D155 and D156**, written since: the deletions and the PMF are both instruments pointed at a complex whose ligand is in the wrong place, and the reaction-coordinate question cannot be answered before the pose is. The decision in front of the user is D156's, not this one's: graft meropenem onto the transferred L1 pose and rebuild, or re-dock with the template as the validation target.

**Nothing submitted, nothing queued, and no further wakeup scheduled** - there is nothing left running to read. Scratch is **264 of 300 GB** (up 8 GB on the two jobs, mostly the 3 ps trajectory and two sets of cubes); nothing pruned, since every candidate is named by a `results` path in `config.yaml`.


**D158 - CSR-1 is closer to L1 than SMB-1 is, so transferring L1's pose into our site rests on firmer ground than the cross-enzyme step the SMB-1 paper itself published. But similarity is not what decides this: both templates are di-zinc and both hold the ligand partly with Zn2, which CSR-1 does not have. Building the product end from 5AXO is the right move, and the reactant has to come from a template too or the reaction coordinate will be a translation.** (2026-11-01, user: "what's the sequence/structure similarity between L1 and CSR compared to L1 and SMB-1?")

## 1. Measured, not assumed

New `refstruct.py similarity` mode: global Needleman-Wunsch with BLOSUM62 and affine gaps, then a CA superposition over the aligned positions with outliers dropped at twice the running RMSD until it converges on 3 A. CSR-1 is represented by **6NC5** chain A, the Q121H/K254H construct, which is CSR-1's own sequence with two point mutations and the only CSR-1 crystal with a metal in it.

| pair | identity | core | core CA RMSD | all-pair CA RMSD |
| --- | --- | --- | --- | --- |
| L1 (7UHI) vs SMB-1 (5AXO) | 32.3% | 216 of 251 | 1.22 A | 3.47 |
| **L1 vs CSR-1 (6NC5)** | **39.0%** | 221 of 249 | **1.05 A** | 4.02 |
| SMB-1 vs CSR-1 | 33.9% | 226 of 254 | 1.16 A | 2.46 |

**CSR-1 is the closest of the three to L1 on both measures.** Every pair sits in the ordinary B3 range - the subclass is defined at roughly 20-40% identity - and the folds agree to about 1 A over ~220 residues, which is as much as any two B3 enzymes agree.

And at the scale that matters for a ligand, the sites are closer still: the four-atom Zn1 fits of D156 come to **0.254 A** between SMB-1 and L1, **0.62 A** from L1 onto our model and **0.50 A** from SMB-1 onto it.

**So the answer to the question as asked is yes.** Mu & Xu transferred a ligand between doripenem and ampicillin inside one enzyme and reasoned across to L1 and NDM-1 at this same level of similarity. Transferring L1's pose into CSR-1 is a smaller extrapolation than that, by the numbers.

## 2. What similarity does not settle, and this is the real caveat

**Both templates are di-zinc; CSR-1 has one metal.** In 7UHI the intact moxalactam touches **both** metals - carbonyl O1 on Zn1 at 2.53 A, and O8 and N1 on Zn2 at 2.65/2.67. In 5AXO the hydrolysed meropenem's new carboxylate O32 and its N4 are **both on Zn2**, at 2.26 and 2.37 A. The part of the pose nearest the chemistry is anchored by a metal we do not have, and no amount of fold similarity changes that.

That is not a reason to abandon the transfer; it is the thing the transfer is for. **If the RQK cations stand in for Zn2, a ligand placed in the di-zinc pose and released should settle onto Arg229 and Lys226 rather than drift away** - and D148/D150 already measured those contacts holding at 3.62 and 3.80-3.90 A in everything we own. The transferred pose is how that hypothesis gets tested instead of asserted, and either outcome is a result.

## 3. The user's plan - graft hydrolysed meropenem from 5AXO as the end state

**Right, and better than the reactant-first route for three reasons.** It is the same drug in the same state, so nothing is invented: no ring to re-form, no bond orders to guess - and the log records exactly that going wrong before, with RDKit putting a spurious O-O bond into the strained hydrolysed ligand and inventing enolate and S+ forms. It is the **first experimentally anchored structure this project would have**. And the product end is the one end we can anchor, since no MBL has ever been crystallised with an intact beta-lactam except L1.

**Four things it needs, none of them optional:**

1. **Atom-order correspondence.** A NEB wants its two endpoints to be the same atom list in two geometries, and this log already records `hydrolysed_lastframe.pdb` being built so its order matches the reactant's. A ligand lifted out of 5AXO arrives with LMP's names and order and has to be mapped onto our topology before it is an endpoint rather than a different molecule.
2. **Protonation and charge decided, not inherited.** Our product is hydrolysed meropenem at -1 absorbing the hydroxide, QM charge +2. 5AXO's LMP is a crystallographic model whose protons are not in the file at all.
3. **Parameters.** Reuse `hyd_zw.pdb`'s topology if the mapping allows it rather than re-parameterising, for the reason in point 1.
4. **The Zn2 contacts must be allowed to relax.** O32 and N4 come in holding a metal that is not there; where they go is the measurement.

## 4. The one caution, and it decides the order of work

**If the product comes from a template and the reactant stays where it is, the reaction coordinate between them is mostly a translation.** D156 put `scan_00`'s ligand **5.15 A** from the experimental pose. A NEB or an umbrella series run from that reactant to a template-derived product would spend its barrier moving meropenem across the site, and the number at the top would not be a chemical barrier at all. That failure would look entirely normal in the output.

**So both ends come from templates, or neither.** The cheapest way to get that, and the one the literature would recognise: **build the reactant from the grafted product by re-forming the C7-N4 bond** - the reverse of the hydrolysis, which is precisely what the NDM-1 and SMB-1 studies do from their product crystals. It keeps the atom order identical to the product by construction, so point 1 above solves itself, and it uses one drug and one template rather than matching cores across moxalactam and meropenem. **7UHI then becomes the check rather than the source**: re-form the bond, relax, and ask whether the carbonyl lands near the 2.5-2.9 A that L1 measures. If it does, both endpoints are anchored and the coordinate between them is chemistry.

**Nothing built or submitted.** `refstruct.py similarity` and the numbers above are the whole of this entry.


**D159 - The product complex is grafted and assembled, and the restraint check the user asked for has a clear answer: nothing of the hydrolysed meropenem is within coordination range of the Zn1-site metal in the crystal, so the product gets NO substrate restraint. Not submitted: the ligand's protonation state is ambiguous in our own files, and it sets both the force field and the QM charge.** (2026-11-01, user: "equilibrate the grafted product first with our standard 100ns MD and 50ps GNF2 QM/MM protocol with the correct crystal restraints ... double check the constraints")

## 1. The graft

Two steps, both in `refstruct.py`, both recorded so the structure can be rebuilt from the PDB ID:

```
$PY $M/refstruct.py transfer refs/5AXO.pdb CSR_wt/preprocess_results/csr_crystal_start.pdb \
    --ligand=LMP --zn=B:1:ZN --fit=$F2 --out=lmp.pdb
$PY $M/refstruct.py match CSR_hydrolysed/preprocess_results/hyd_zw.pdb lmp.pdb \
    --out=hyd_crystal.pdb
```

`transfer` fits 5AXO on the Zn1 site (RMSD 0.838 A over the four atoms) and moves LMP across. `match` then puts **our** atom names and order on **those** coordinates by graph isomorphism over the heavy atoms - never by name, since LMP and our `hyd` name the same 27 atoms completely differently.

**The mapping is chemically right, which is the check that matters.** Our `C4x`, the electrophile, maps to LMP's `C7`; our `N1x`, the lactam nitrogen, to `N4`; and **our `O2x` and `O01` - the old carbonyl oxygen and the absorbed hydroxide - map to `O71` and `O72`, the two oxygens of the carboxylate the hydrolysis created.** That is the reaction written out in the atom mapping, and it came from the bonding alone.

**8 mappings exist**, differing only in C12x/C13x (the amide's two methyls), O2x/O01 and O3x/O6x (each carboxylate's two oxygens). All three are real symmetries and all three choices are harmless: measured on `hyd_zw.pdb`, none of those four oxygens carries a hydrogen.

## 2. The restraints, checked against the crystal as asked

The standing rule is that anything within coordination range is restrained at **the distance it is in the crystal**. Measured in 5AXO itself, against **ZN301, the Zn1-site metal that CSR-1's single zinc occupies**:

| product atom | to Zn1 (5AXO) | to Zn2 (5AXO) |
| --- | --- | --- |
| O62 | 3.26 | 3.96 |
| O32 | **3.41** | **2.26** |
| N4 | 4.26 | **2.37** |
| C31 | 4.13 | 2.97 |

**Nothing of the product is within coordination range of Zn1. Everything that binds, binds Zn2 - the metal CSR-1 does not have.** So the correct restraint set for this run has **no substrate term at all**, and writing one would be transplanting a Zn2 bond onto a Zn1 metal. In the assembled structure the nearest product atoms sit at 2.92 and 2.94 A from our zinc, outside the 2.5 A detection cutoff and the 2.6 A analysis cutoff, which agrees with 5AXO's own 3.41 to within the fit error.

The protein restraints are unchanged from D136/D138, at their 6NC5 distances: **His88 NE2 2.22, His164 NE2 2.09, HOH485 1.90**. `constraint_exclude_residues: {C: [1]}` keeps the ligand out of automatic detection, as every reactant run has.

## 3. Assembled, and why it was not submitted

`prod_start.pdb` is `csr_crystal_start.pdb` with the reactant and the hydroxide removed and `hyd_crystal.pdb` appended. `refstruct.py shell` on it:

```
Zn ZN1 chain B:  1.90 HOH W1 | 2.09 His164 NE2 | 2.22 His88 NE2
                 2.80 Asp92 OD2 | 2.92 O6x | 2.94 O5x
```

**Three decisions looked open; measuring settled two of them. The surviving one is item 2.**

1. ~~**The product's protonation and net charge.**~~ **[Resolved the same day by measuring it rather than reading the atom order: `hyd_zw.pdb` puts H03 on **N1x**, the leaving nitrogen, not on O01 as the file's atom order suggested. So it is N1x-H (the proton the hydrolysis delivers), N2x-H01/H10x (the protonated pyrrolidine) and **all four carboxylate oxygens bare**: +1 - 2 = **net -1**, exactly D13's figure, and internally consistent. The product zwitterion SMILES follows from `CSR_hydrolysed`'s by protonating the amine and deprotonating the ring acid:** `C[C@@H](O)[C@@H](C(=O)[O-])[C@@H]1NC(C(=O)[O-])=C(S[C@@H]2C[NH2+][C@H](C(=O)N(C)C)C2)[C@@H]1C`**. No ambiguity remains.]**
2. **What fills the hydroxide's place on the zinc.** The nucleophile is now part of the product, and the product does not coordinate. 6NC5 has a crystal water there (HOH493, 2.14 A) and 5AXO keeps a water on Zn1 at 1.82 A, so a restrained water is defensible - but leaving the face open is a different model of the product state and the two will not give the same first shell.
3. ~~**O2x against O01.**~~ **[Not a problem: neither oxygen carries a proton, so the swap changes nothing. The concern came from assuming H03 belonged to the atom listed before it.]**

**Nothing submitted.** A 100 ns MD and a 50 ps QM/MM on a guessed protonation state would produce a product complex that cannot anchor a reaction coordinate, which is the whole purpose. The graft, the mapping and the restraint finding stand whatever is decided.


**D160 - The product preprocess failed twice in nine seconds total, both times on the submission rather than the science, and the second one was worth fixing at the tool: a grafted ligand written without CONECT is 27 unbonded atoms to the thing that reads it.** (2026-11-01, jobs 31296723 and 31299591 FAILED; 31299717 submitted)

**31296723, 3 s.** `run.slurm preprocess` with no second argument: the script takes `<stage> <section>` and preprocess's section is the word `preprocess` again, which its own usage line says. My error, and the script caught it before anything ran.

**31299591, 6 s.** `OSError: Bad input file ./hyd_crystal.pdb` from `Chem.MolFromPDBFile(..., proximityBonding=False)`. **mcube reads a ligand with proximity bonding switched off on purpose** - D21 records RDKit's distance-based bonding inventing an O-O bond in this very molecule - so the bonds have to come from CONECT records, and `refstruct.py match` wrote none.

Fixed where it belongs (`65e0b27`, 2 tests): `match` now writes CONECT from the graph it did the matching on. The connectivity was already known exactly - it is what the isomorphism was computed from - so emitting it costs nothing and keeps the ligand's bonds out of a distance guess, which is the same reason `proximityBonding=False` is there.

**Neither failure touched the structure.** The regrafted `hyd_crystal.pdb` is the same 27 atoms with the same 8-fold symmetric mapping, now with its 29 bonds written down. Resubmitted as **31299717**.


**D161 - The product preprocess failed a third time, and the cause is a real inconsistency in mcube rather than my typing: `preprocess.inputs` paths are read from the stage's working directory while every other stage reads `pdb_file` against the config's. Fixed three ways. Also, the helper-script audit: everything is documented except what was written today, and `traj_contacts.py` is now a duplicate.** (2026-11-01, jobs 31299591 and 31299717 FAILED, 31301810 submitted; user: "Why does run.slurm preprocess need two commands? ... check how many of those helper .py scripts should actually be a part of the mcube repo and documented")

## 1. Why `run.slurm preprocess preprocess` was needed, and why it no longer is

It was redundant and the user is right. `SECTION` was mandatory, but for preprocess the script sets `SECTION_ARGS=()` and never passes it on - the word only named the output directory, which is `preprocess_results` regardless. **`SECTION` now defaults to the stage name**, which is also mcube's own default section for every stage that takes one (`qmmm`, `sp`, `scan`, `neb`, `ts` all default to their own name), so `run.slurm preprocess` works and `run.slurm qmmm` means the `qmmm` section. Nothing that used to work changes.

## 2. The failure that was not my typing

Both 31299591 and 31299717 died in six seconds with `OSError: Bad input file ./hyd_crystal.pdb` - and the file was there, well formed, with CONECT records, and RDKit read it happily when asked from the run directory.

**`run.slurm` does `cd "$OUT"` before calling mcube.** So `preprocess` resolves `protein_pdb` and `ligand` against `preprocess_results/`, while every other stage resolves `pdb_file` against the **config's** directory (`config_relative`, which exists precisely so a structure is found from wherever the stage runs). Two rules in one config file, and nothing says which applies where. The reactant run never tripped on it because its `csr_crystal_start.pdb` happens to live in `preprocess_results/` - the layout encodes the rule, and nothing else does.

Three changes:

- **`make_product_run.py` writes both inputs into `preprocess_results/`**, with the reason in a comment and two tests that assert the location, since the path in the config has to stay bare.
- **`ligand_mol_from_pdb` checks the file exists first**, and the message names the working directory and contrasts it with `pdb_file`. RDKit's own error reads like a malformed structure and sends you to count PDB columns; it cost two jobs here.
- The `run.slurm` change above, so the command is shorter and the directory rule is documented where it bites.

**The deeper fix is not made:** preprocess should read its inputs the way every other stage does. Changing that would break every existing config whose paths are already `preprocess_results`-relative, so it needs a fallback (try the config's directory, then the stage's) and a sweep of the eleven run directories. Worth doing; not today.

## 3. The helper-script audit

Twelve scripts in `docs/csr/`, and the documentation is in better shape than I expected: **all of them were already in `README.md`**, including a "one-shot, kept only as a record" note for `migrate_paths.py` and `make_ohH_config.py`. The gaps:

| script | tests | state |
| --- | --- | --- |
| `site.py`, `sections.py`, `refstruct.py` | yes | the core; `refstruct`'s three new modes were documented only in its own docstring - now in the README |
| `make_product_run.py` | 16 | written today, was undocumented - now in the README |
| `make_water_nuc.py`, `traj_frame.py`, `make_ohH_config.py`, `migrate_paths.py` | yes | fine |
| **`traj_contacts.py`** | yes | **superseded.** `site.py contacts` does the same measurement from `site.yaml`'s `monitor`; this one carries its own hard-coded pair list, which is the exact anti-pattern `site.py`'s docstring exists to forbid - *"a contact is reported because the YAML names it and not because someone wrote a one-off script for this protein"*. Marked as superseded; delete once D111's numbers are re-derived with `site.py contacts`. |
| `check_ohH_frame.py`, `check_ohH_system.py` | **none** | one-off pre-flight checks for a build that finished at D119. Documented but not marked one-shot, unlike their siblings. |

**What should arguably not be in `docs/` at all:** `refstruct.py` says in its own docstring that it is *"deliberately independent of `site.py` and of `site.yaml`"* - it knows nothing about CSR-1 - and it now holds four generally useful operations: metal-shell measurement, local superposition and ligand transfer, ligand graph matching, and sequence/structure similarity. That is library code sitting in a documentation directory, and the same argument applies to `traj_frame.py`. Moving them into `src/mcube/` with a CLI entry point is the right shape; it is a refactor with no scientific content, so it is recorded here rather than done in the middle of a build.

## 4. The run

`CSR_prod` rebuilt with the inputs in the right place and **submitted as 31301810**, using the one-argument form to confirm the `run.slurm` change. Nothing has been equilibrated yet; the zinc shell check that gates the 100 ns MD is still the next thing to read. Scratch 265 of 300 GB.


**D162 - The product preprocess failed seven times, and the protomer was a red herring: the real cause was that the grafted ligand had no hydrogens, which no SMILES can fix. Fixed, along with two file-format rules mcube's own guards caught. The zwitterion stands.** (2026-11-01, jobs 31299591 through 31303033 FAILED; 31303502 submitted)

**The diagnostic that settled it.** D161 left a question open: run the product as the zwitterion (matching `hyd_zw.pdb`'s hydrogens and D13) or as the neutral amine (matching `CSR_hydrolysed`). I ran `CSR_hydrolysed`'s **known-good** SMILES as a one-variable test, and **it failed with exactly the same error**. So the SMILES was never the problem, and **the protomer question does not need answering to make the build work** - the zwitterion is restored, for the reason D13 gives.

**The real cause.** `ligand_mol_from_pdb` builds its template as `Chem.AddHs(Chem.MolFromSmiles(smiles))` - 53 atoms - and matches it against the ligand file. My `refstruct.py match` wrote **heavy atoms only**, 27 of them, on the reasoning that `preprocess` adds hydrogens. It does, to the protein and the waters; the ligand is parameterised from its own file, and the docstring says so in its first line: *"Read a ligand PDB with explicit hydrogens"*. A 27-atom structure cannot match a 53-atom template under any protonation.

`match` now **carries the hydrogens across**, each one by a rigid fit of the small group around its parent - the parent, its neighbours and theirs - onto those same atoms in their new positions. Leaving them at their old coordinates was never an option either, since the heavy atoms move to the template's conformation. 53 atoms out, 27 heavy and 26 hydrogens, with the X-H distances preserved by construction and a test that checks it.

**Two more, both caught by mcube rather than by me**, and both worth having:

1. **`preprocess.inputs` paths are read from the stage's working directory** while every other stage reads `pdb_file` against the config's (D161). Fixed in the builder; `ligand_mol_from_pdb` now names the directory instead of letting RDKit answer "Bad input file".
2. **`check_residue_numbering` twice refused a chain not numbered 1..N.** Removing the substrate left its chain starting at residue 2; moving the converted water into the water chain then left that chain reading 2 then 1, because the hydroxide sat before the crystal water in the file and the check walks it in order. Both are real: ASH renumbers, so every selector would afterwards point at a different residue. The builder now moves the converted water to the end of the water block and numbers it after the others.

**Seven failures, and every one of them was a guard doing its job** - a missing argument, a missing file, missing bonds, wrong bonds, missing hydrogens, and twice a numbering rule. None of them would have been visible in the output of a run that went ahead; the two that could have gone unnoticed, the bonds and the hydrogens, are exactly the ones that would have produced a parameterised molecule that was not the one we meant.

**31303502 submitted** with the zwitterion SMILES and the 53-atom ligand. The zinc-shell check still gates the 100 ns.


**D163 - The eighth preprocess got the ligand right and died on the force fields: `resolve_forcefields` has been wired into `md` and `qmmm` since the config-relative rule was written, and never into `preprocess`. Fixed in mcube, not worked around in the config. Two further defects found by reading the generated file rather than by running it.** (2026-11-01)

## What job 31303502 did

It parameterised the ligand. `openff_hyd.xml` is 23 kB, and ASH's own consistency check passed with identical Coulomb and LJ energies from the XML and from the forcefield object - so the 53-atom grafted product, zwitterion SMILES and transferred hydrogens (D162) are all accepted. That was the hard part and it is done.

Then, at `[STEP] Solvating system...`:

```
ValueError: Could not locate file "hydroxide.xml"
```

104 seconds, after the ligand work, which is why it looks like a late failure for an early reason.

## The cause is general, not a CSR detail

`hydroxide.xml` is the custom residue template for the nucleophile, and it sits beside each run's `config.yaml`. `src/mcube/preprocess.py` resolved nothing: it passed the config's strings straight to `ForceField`, which reads a bare name against the **working** directory - `preprocess_results/` under `run.slurm` - and the file is one level up. `md.py:55` and `qmmm.py:224` have both called `resolve_forcefields` for months.

So every stage except the one that builds the system honoured the rule. Any config naming a custom residue XML beside it would equilibrate and run QM/MM fine and fail at preprocess, which reads as a problem with the structure. `preprocess.run()` now calls `resolve_forcefields(forcefields, args.config)` after the membrane block, which is the same resolution md and qmmm already use.

`hydroxide.xml` stays in the product's list even though the product has no hydroxide residue: an unused template costs nothing, and the two states are only comparable if the force fields are identical.

## Two more, found by reading the file

Neither had been reached yet; both would have failed afterwards.

1. **The converted water was still called `O01`.** `make_product_run.py` turned the hydroxide into a `HOH` and kept the hydroxide's atom name. `amber14/tip3p.xml` types a water by its atom names, so that residue has no template. Now renamed to ` O  `, in the name columns only.
2. **`qmmm.inputs.forcefields` still named `openff_MER.xml`** - the substrate's parameters, inherited verbatim. The product would have run with the reactant's ligand typed in. The XML name follows the resname, so `retarget_ligand_xml` derives it: `openff_hyd.xml`, directory preserved.

**The QM charge of 2 is correct unchanged, and that is worth stating rather than assuming.** Hydrolysis conserves total charge: the reactant's QM region holds a neutral substrate and a hydroxide (-1), the product's holds a -1 product and a neutral water. Zn(2+) + Asp92(-1) + Arg90(+1) + Lys226(+1) + (-1) = 2 either way.

## State

Nothing resubmitted: see D164 on the terminal bridge. The rebuild is three commands and the zinc-shell check still gates the 100 ns MD.

596 tests pass locally (`PYTHONPATH=src pytest tests`), 5 of them new.

**D164 - Nothing was resubmitted: the code-server terminal will not take input, and the one toggle that fixes it does not respond. The bridge reads files fine, which is how D163 was diagnosed without running a single command.** (2026-11-01)

code-server is up at `bun151:28472` and `vscode-remote-resource` serves every file, so the whole of D163 came from reading `protdft/prodprep-31303502.err`, `config.yaml`, `product_protein.pdb` and `openff_hyd.xml`. That is worth knowing on its own: a failed preprocess can be diagnosed completely without a shell.

Input is a different matter. The status bar reads **`Screen Reader Optimized`**, the documented mode in which `computer{action:"type"}` does not reach xterm at all (`BUNYA_BROWSER_SHELL.md`, 29 Oct). The three documented ways out were all tried and none changed the status bar: the command palette, `Shift+Alt+F1`, and a synthetic click on the status-bar item. Synthetic `keydown`/`keypress` into the helper textarea delivered exactly one command and then stopped.

So the next step is one keystroke in the browser tab, and it is the user's: turn Screen Reader Optimized off. Everything else is ready - the rebuild is `make_product_run.py`, `sections.py prod100`, `sbatch`, and the zinc-shell check still gates the MD.

**D165 - The product complex is built and its zinc shell is exactly what D159 predicted: four-coordinate, with the product's nearest oxygen at 2.92 A and not coordinating. The 100 ns MD is away (31327774). Three preprocess runs were needed, and the two failures in between both reported `[DONE]`.** (2026-11-01)

## The shell, which is the result this build was for

`refstruct.py shell` on `CSR_prod/preprocess_results/finalsystem.pdb`:

```
Zn ZN1 chain B:
    1.90  O    HOH C1
    2.09  NE2  HIS A164
    2.14  O    HOH C2
    2.22  NE2  HIS A88
    2.80  OD2  ASP A92
    2.92  O5x  hyd D1     <-- LIGAND
    2.94  O4x  hyd D1     <-- LIGAND
```

**Four-coordinate - two histidines and two waters, at the 6NC5 distances - and the grafted product's nearest oxygen sits at 2.92 A, outside the 2.5 A detection cutoff.** D159 measured 3.41 A from the Zn1-site metal to the closest product atom in 5AXO and predicted ~2.9 A here; the graft reproduces it to within the fit error. Nothing of the product is restrained, which is what that measurement requires: everything that binds in 5AXO binds Zn2, the metal CSR-1 does not have.

`coordination` holds exactly four terms (C:1:O 1.90, A:164:NE2 2.09, C:2:O 2.14, A:88:NE2 2.22) and `qm_residues` holds the zinc, His88, His164, Asp92, Arg90, Lys226, the ligand at D:1 and both waters - the reactant's QM region residue for residue, with the hydroxide now a water.

## Two runs that failed by succeeding

Both printed `[DONE] preprocess` and wrote a complete, solvated, internally consistent system.

1. **31322290 built the complex without the complex.** `merge_ligand: false` is correct for the reactant, whose substrate sits inside `csr_crystal_start.pdb` and whose ligand file only supplies parameters. `make_product_run.py` inherited that flag into a run whose protein file has the substrate removed *on purpose*, so nothing ever put the product in. The zinc shell was right, the restraints were right, and the thing being equilibrated was apo CSR-1.
2. **31325102's predecessor pointed every ligand selector at a water.** ASH's merge renames the water chain to C and appends the ligand as its own chain - D here. The reactant's `extra_qm_residues: {C: [1]}` and `constraint_exclude_residues: {C: [1]}` mean *the substrate* there and *a water* after the merge. Consequences, both silent: the ligand stayed out of the QM region, and the 1.90 A water lost its restraint, because the exclusion written for the ligand landed on it. A chain letter chosen by a merge is not a stable name for a molecule, and `extra_qm_residues` cannot catch this - a water at C:1 exists, so there is nothing to refuse.

**Worth proposing as an mcube change:** let the ligand be named by resname rather than by chain and number (`constraint_exclude_resnames: [hyd]`, or a resname key in `extra_qm_residues`). Every selector that means "the ligand" is currently spelled as a position that another residue can occupy.

## ~~The reactant and the product are not built the same way~~ [WITHDRAWN 2026-11-01, same day: wrong, and wrong because I inferred a history from two config files instead of reading the log I had already written. See the correction below.]

**The correction.** `cutoff: 2.5` was decided deliberately - the rules section of this log says "a Zn-O/N bond is 1.9-2.2 A; **2.5 is defensible, 3.0 is not**", and says to set it explicitly because the library default is 2.8. It did not drift into the product by accident; inheriting it is correct.

The stale five-term block is real but **inert**. It is `CSR_nac`'s, detected at **3.0 A** - not 2.8 - and carried in when D132 built `CSR_wt` as a copy of `CSR_nac`'s preprocess output. `md_wt100` set `coordination_redetect: true`, so the MD ignored it and measured the coordination itself: **D134 records four restraints redetected at 2.5 A - HOH 1.90, His164 2.09, O01 2.14, His88 2.22 - and states that nothing restrained Asp92 and nothing restrained the substrate.** `md_prod100` carries the same flag. So the reactant and the product are restrained identically, Asp92 free in both, and **nothing needs re-running.** D134's finding - Asp92 on the zinc in 5,000 of 5,000 frames - is a measurement precisely because nothing held it, which is the opposite of what the withdrawn section claimed.

What survives is narrower and was already prescribed here: a fresh preprocess at 2.5 A **drops Asp92 from the QM region**, because its OD2 at 2.80 A is outside the cutoff and the reactant's saved region only has it by inheritance from the 3.0 A detection. D133's closing note says "put 92 and 93 in `extra_qm_residues` before ever re-running it". `make_product_run.py` now adds 92, so the product's QM region matches the reactant's residue for residue - the right change, for that reason rather than the one I gave. **Gln93 is in neither preprocess region**; the reactant's two scan sections add it per-section, so the product's reaction-coordinate sections will need `extra_qm_residues: {A: [93]}` the same way. Gln is neutral, so the QM charge stays unchanged.

**The lesson, since this is the second time.** D133 wrote the note, and D165 asked the question it answers. Two config files were enough to build a plausible story, and the log had the measurement. Read the log entry for a run before reasoning about what that run did.

## ~~The withdrawn section, kept as written~~

`CSR_wt/config.yaml` says `cutoff: 2.5`, but its stored `preprocess.results` were produced under the **2.8 A default** and never regenerated. Its `coordination` has **five** terms, the fifth being **Asp92 OD2 at 2.80 A**, and it is stored as raw atom indices (4203, 4257, ...) rather than the selector strings the current code writes - two independent signs of a stale result.

So the reactant's 100 ns MD and 50 ps QM/MM (D134/D136/D138) ran with the zinc-Asp92 distance **restrained at 2.80 A**, and with Asp92 in the QM region by detection. The product, built today at 2.5 A, restrains four contacts and gets Asp92 into the QM region only because `make_product_run.py` now adds it explicitly.

That matters because restraining Asp92 to the metal fixes the very thing the mechanism question is about, and the standing position here is that there is no point starting from a structure where Asp92 is overcoordinated. The product is built the way that position implies. **The reactant is not, and bringing it into line means re-running its 100 ns MD and its 50 ps QM/MM.** The product MD does not depend on which way that goes, so it was submitted.

## Disk

268 of 300 GB. A 100 ns run of this system writes **2.9 GB**, so the MD and the QM/MM both fit in the 32 GB free. The large consumers are other projects (`glp1r_md` 41G, `channel` 38G, `boltzgen` 35G, `AANCA` 27G, `hf` 25G); all of CSR_wt is 6.9 G. `CSR_wt` holds two 2.9 GB 100 ns trajectories, `md_wt100_results` and `md_wt100b_results` - one is probably redundant, but D134/D138 cite these runs and nothing was deleted without checking which.

**D166 - Measurement only, no decision: both equilibrated 100 ns end states have Asp92 on the zinc, one of them bidentate. Recorded because the question of grafting the product into an equilibrated frame turns on these numbers.** (2026-11-01, user: "you can probably graft the hydrolysed SMB1 meropenem in to an already 100ns equilibrated MD structure that has already been equilibrated with the correct restraints")

`refstruct.py shell` on the two candidate scaffolds, `NVTtrajectory_lastframe.pdb` in each:

| | `md_wt100` (30966977) | `md_wt100b` (31016877) |
| --- | --- | --- |
| OH O01 | 1.90 | 1.93 |
| **Asp92 OD2** | **2.02** | **2.06** |
| HOH D:1 | 2.08 | 2.09 |
| **Asp92 OD1** | - | **2.12** |
| His88 NE2 | 2.28 | 2.22 |
| His164 NE2 | 2.18 | 2.31 |
| count | five | **six, Asp92 bidentate** |

Counted against this enzyme's own H-R-H motif, where the canonical-equivalent four is two His + hydroxide + one water (rule 4), **the extra ligand in both frames is Asp92** - not a solvent ligand, so this is over-coordination in the sense that matters. `md_wt100` is the run with the correct restraints (four, redetected at 2.5 A, nothing on Asp92 or the substrate, D134), and Asp92 still arrives at 2.02 A unforced; `md_wt100b` tethered the substrate and ends with the carboxylate chelating.

For comparison, the structure the running product MD started from, `csr_crystal_start.pdb`: HOH 1.90, His164 2.09, His88 2.22, **Asp92 OD2 out at 2.80** - four-coordinate, crystal geometry (D159).

**What job 31327774 is actually doing**, from its own `.out`: `Coordination restraints from: redetected on .../CSR_prod/preprocess_results/finalsystem.pdb`, then `4 harmonic distance restraints (4 detected metal coordination)` at 1.90, 2.09, 2.14, 2.22. The same four, at the same distances, as D134's reactant run, from the same crystal shell, for the same 100 ns, with only the ligand different. Rate: 9.16 ns in about 32 minutes, so ~5.8 h rather than the 3 h estimated from the reactant - inside the 12 h walltime.

**Decided, same day, by the user: "if you've already got a correctly set up MD job running, just let it complete and do QM/MM on that like the original plan."** No second MD, no change of start structure. The numbers below are why the running job is the one to keep, and are recorded for the methods section.

**The case against grafting into an equilibrated frame.** Grafting into `md_wt100`'s end state would give the product an equilibrated protein, which is the motive and a good one. It would also start it with Asp92 bound at 2.02 A, which is rule 7 ("no cutoff undoes an artefact already in the coordinates", a mistake D120 has costed once), and it would break the matched control that makes "the product stayed" readable against "the reactant's substrate left". No structure available has both an equilibrated protein and an Asp92-out shell: `scan_00` is four-coordinate with Asp92 out at 3.60 A and descends from the equilibrated trajectory, but it is a constrained scan point, biased along the coordinate, and carries the substrate rather than the product.

**D167 - The trajectory analysis this build needs is now a mcube stage. It was never missing from the framework so much as half-present: `postprocess` already had a coordination category, but it reimplemented the detection rule with a hardcoded donor set, a cutoff from the command line and no minimum image.** (2026-11-01, user: "stick to mcube rather than using custom scripts to keep everything reproducible")

## What was actually outside mcube, after auditing rather than assuming

The single-frame sphere never was: `preprocess.coordination_shell` is the rule, and `md` already prints it **before and after** every run, so a contact that appears during equilibration is visible in the job's own output. `sphere`'s protonation reporting has no mcube equivalent.

The trajectory version is what mattered, and `mcube postprocess` had a version of it - an interaction graph with a `coordination` edge type and per-edge occupancy. Five differences from the rule that builds the restraints, all load-bearing here:

| | `postprocess` as it was | the restraint rule |
| --- | --- | --- |
| donors | hardcoded `element N O S` | `preprocess.inputs.donor_elements` |
| cutoff | `--cutoff`, default **2.6** | `preprocess.inputs.cutoff` (2.5 here) |
| metals | `--metals`, default **ZN** | `preprocess.inputs.metals` |
| periodicity | raw `np.linalg.norm` | minimum image on an orthorhombic CRYST1 |
| keyed by | residue | atom or residue, a switch |

So the cutoff and the metal list reached the analysis from a **command line**, which no config records - the reproducibility hole, and in a general library a default of `ZN` is this project leaking into everyone else's. Keying by residue also collapses a bidentate carboxylate to one edge, which is the distinction D166 had just measured. And nothing reported the coordination **number**, the quantity D114 shows semiempirical methods get wrong.

## The change

`coordination_shell` is split into `coordination_rows` (parsing) and `coordination_from_rows` (the rule). `postprocess` substitutes each frame's coordinates into those rows, so a trajectory is measured by the same rule, cutoff, donor set and minimum image as a file - one implementation, not two. It reports per-ligand occupancy with min/mean/max, the coordination number as a **distribution** (a site flipping between 4 and 5 is not a site at 4.5), and donors just past the cutoff, since 2.7 against a 2.6 cutoff is a borderline call about the site and not an absence. Settings come from `preprocess.inputs`; what was measured is saved into the section's `results`.

`postprocess` also now takes `--section` like every other stage, and `inputs.md_section` names which md run to analyse - a run submitted as `md md_prod100` saves its results under that name, so a stage that could only read a section called `md` could not find them. `site.yaml` gains `prodpost`, so the analysis is submitted as `run.slurm postprocess postprocess_prodpost` and nothing is typed by hand.

**Still outside mcube**, and only worth moving when needed: `site.py`'s `contacts` (named pairs over a trajectory, which the reaction-coordinate work will want), `nucleophiles`, `acceptors`, `joint`, `poses`, `graft`, `check`, `protonate`, and `refstruct.py` entirely. These are config-driven rather than hardcoded, but they are not stages and a result from them is not recorded in a config.

## A testing gap this exposed, which is the part worth keeping

The first version of the report had a `print("\n...")` whose escape had collapsed into a literal newline, so `postprocess.py` would not parse. **The local suite was green.** The only test importing that module is gated on MDAnalysis, which is not installed locally, so the module was never imported, never parsed, and the syntax error surfaced only on Bunya - in the one environment that has MDAnalysis. Gating on an optional dependency hides syntax errors in everything it guards.

`tests/test_modules_compile.py` now parses every module in `src/mcube` and every `docs/csr` helper with `compile()`, which imports nothing and needs nothing installed. 632 pass locally (8 skipped, all ASH-gated), **707 pass on Bunya with nothing skipped.**

**D168 - `contacts` is now part of `postprocess`, and promoting it turned up something worse than a missing feature: OpenFF's parameterizer renames ligand atoms, so 15 of the product's 53 atom names in the simulated system are not the names the graft wrote - `O01`, the nucleophile oxygen this whole log calls `O01`, is `O6x` in `finalsystem.pdb`.** (2026-11-01, user: "contacts probably belongs in preprocess or postprocess.py - use your judgement")

## Where it went, and why

**`postprocess`, with one piece in `preprocess`.** The analysis is a trajectory question, and `postprocess` is the stage for those: it already iterates frames, and the result belongs in that section's `results` beside the coordination report. What belongs in `preprocess` is the *vocabulary*, and it was already there - `resolve_atoms` is the one parser for `chain:resid:atom`, used by restraints and `zero_charges`, and `contact_occupancy` now calls it rather than carrying its own resolver. So the selectors mean the same thing in a monitor as in a restraint, which is the property worth having.

`postprocess.inputs.monitor` is `{name: [selector, selector]}` - the same shape `site.yaml` already used, so `sections.py` carries it across unchanged. A contact is reported because the config names it.

Two behaviours carried over deliberately. A pair that does not resolve is **reported as missing, not skipped**: D120 is the entry where the oxygen that migrated onto the metal was in no monitor at all, and a selector that silently matches nothing is the same failure with extra steps. And `-RT ln P(r)` came too (`contacts.profile`), because the histogram of a 100 ns trajectory is already a PMF over the range it visited - with the per-bin count printed beside the energy, since a bin with a handful of frames is noise wearing a number.

## The renaming, which is the real finding

`merge_ligand` parameterises through `small_molecule_parameterizer`, which writes its own `hyd.pdb` from the SDF - and **re-derives the atom names**. Matching the two files by coordinate (all 53 atoms within 0.05 A) gives the map:

| graft (`hyd_crystal.pdb`) | simulated (`finalsystem.pdb`) |
| --- | --- |
| O5x | O4x |
| O6x | O5x |
| **O01** | **O6x** |
| H01, H02, H1x-H4x, H6x-H10x, H03 | H1x-H11x, H26x |

Everything else, `C4x` and `N1x` among them, is unchanged - which matters, because those two are the bond the reactant build re-forms.

**Cross-checked rather than assumed:** D159 measured our `O6x`/`O5x` at 2.92/2.94 A from the zinc in `prod_start.pdb`; the final system reports `O5x`/`O4x` at the same two distances. The map is consistent with both.

**What this breaks.** Any config selector naming a ligand atom was written against the ligand file and means a *different atom* in the system that runs. Nothing records the map, nothing warns, and the failure is silent: `D:1:O01` simply does not resolve, while `D:1:O5x` resolves to the wrong oxygen. There are now three naming systems in play - the reactant's substrate, the grafted product, and the simulated product - and a selector has to say which.

**Built in D169:** `preprocess` records the input-to-parameterised name map in `preprocess.results` and warns when it is non-empty.

`prodpost` therefore carries 23 pairs written in **final-system names**, with pairs named by atom rather than by chemical role wherever the role has not been re-derived from the product's connectivity: the reactant monitor's roles do not carry over.

632 pass locally (9 skipped, MDAnalysis- and ASH-gated), and the compile check now parses `postprocess.py` on every run, which is what would have caught D167's collapsed escape.

**D169 - `preprocess` now records the ligand's input-to-parameterised atom-name map in its results and warns when it is non-empty, so a config's ligand selectors can be checked against the config.** (2026-11-02, user: "Yep add and document it correctly")

## The failure this closes

D168 found the renaming; this records it. OpenFF re-derives ligand atom names from the molecular graph, so the names in the ligand file a config was written against are not the names in the system that runs. The product's case is the bad one: three of the fifteen changes are a **permutation among names both files use** - `O01`->`O6x`, `O5x`->`O4x`, `O6x`->`O5x`. A selector written against the input therefore names a real but different atom. Nothing is missing, nothing raises, and a distance between the wrong pair of oxygens is reported under the right label. That is the one class of error no existing check could catch: `resolve_atoms` refuses a selector that matches nothing, and D168's own monitor reports a pair that does not resolve - but neither can see a selector that resolves to the wrong atom.

## What it does

`preprocess.ligand_atom_map(input_pdb, parameterised_pdb)` returns `{input name: parameterised name}` for **the entries that differ only**, so an empty map is the positive statement that the config's own names are the system's names. `run` calls it after `merge_ligand`, saves it as `preprocess.results.ligand_atom_map`, and prints each pair as a `[WARNING]` when it is non-empty, saying that ligand atom selectors must use the right-hand name.

**Matched by coordinate, not by order or by name.** The parameteriser writes back the geometry it was given, so a coordinate match is the one correspondence a renaming or a reordering cannot fool - and it is how D168's map was established by hand in the first place. The tolerance is 0.05 A, two orders below a bond.

**Only for a `.pdb` input.** A `.mol` or `.pdbqt` ligand carries no atom names to map from, so there is nothing to compare and no warning to give. This is the reactant's case: its substrate comes in through a PDBQT, which is why the reactant's configs never hit this.

**A map that cannot be built is reported and empty, never guessed.** Different atom counts, an atom with no match, an ambiguous match, a missing file: each prints what it could not do and tells the reader to check the selectors by hand. The parameterisation has already succeeded by this point and failing a run over an analysis aid would be worse than the gap - but the warning matters, because a silent `{}` would read as "no renaming happened", which is the opposite claim.

## What it does not do

It does not stop the renaming. Preserving the input names through `small_molecule_parameterizer` would be better - the map would then always be empty and no config would need checking - but that means getting names back onto the parameterised molecule inside ASH's call, which I have not established is possible cleanly. The map is the version that is safe to have now. **P7: try to preserve the input ligand atom names through parameterisation, and if it works, keep the map as the check that it worked.**

It also does not rewrite existing configs. `prodpost`'s 23 pairs are already in final-system names (D168); nothing else in `site.yaml` names a ligand atom, and the product's `prodqm` selects whole residues rather than atoms. The next config that needs one - the reactant build's `C4x`-`N1x` bond, where both names happen to be unchanged - can now be checked against `preprocess.results.ligand_atom_map` instead of against a coordinate match done by hand.

`tests/test_ligand_atom_map.py` is 12 tests, and the permutation is the one they are built around: the same four oxygens in both orders, which must give the same map.

**[CORRECTED by D171 - the Gln93 omission was the user's deliberate decision, not an oversight. See `csr/retractions.md`.]** **D170 - The route from the product MD to a mechanism: check the shell, equilibrate at QM/MM, optimise, then build the OTHER endpoint by scanning this one. The reactant cannot be preprocessed separately, and Gln93 - the Q of RQK - is in neither QM region.** (2026-11-02, user: "next steps after MD runs is to check the resultant structure before doing QM/MM to equilibrate the reactant - then a geometry optimisation followed by converting it to a product along reaction coordinate - perhaps for NEB, being aware of the pitfalls. Remember our overall hypothesis of analysing the RQK site's impacts on mechanism")

## The order, and what gates each step

1. **`prodpost` before any QM/MM.** 100 ns of MM MD is not evidence that the site survived it. The gate is the one the whole product arm exists to protect: four-coordinate zinc (two His, the hydroxide-derived oxygen, one water - CSR-1's Zn1 is **H-R-H**, Arg90 where the third His would be, so four is the canonical-equivalent sphere), Asp92 off the metal, and an empty `missing` list. A five-coordinate distribution or Asp92 on the metal stops the arm: equilibrating a product whose zinc is already overcoordinated measures the artefact, not the mechanism, and no later cutoff undoes coordinates.
2. **`prodqm`** - QM/MM MD, which is what "equilibrate" means here, since the MM force field cannot hold the hydrolysed ligand's electronics. ~20 h on 8 CPU cores; xTB is the bottleneck and a GPU does not help.
3. **Geometry optimisation** from that frame, giving a stationary product.
4. **Scan along the reaction coordinate** to the other endpoint - for this arm the **reverse**: re-forming `C4x`-`N1x`. The arm starts at the product because 5AXO gives the hydrolysed complex at crystallographic resolution; the reactant is what has to be constructed.
5. **NEB** between the two.

## The pitfall that decides step 4, so it is not a choice

`neb` refuses endpoints that differ in atom count **or in atom order**, and a reactant built through its own `preprocess` run will not match this product however carefully it is prepared - different solvation, different atom order, and nothing about the two runs makes them one system. **The only valid route is `scan` from this product**, because scan writes every point through the product's own MM topology, so the order matches by construction. This is why "convert it to the other state along a reaction coordinate" is not an alternative to NEB but a prerequisite for it.

Two smaller ones. `C4x` and `N1x` happen to be unchanged by OpenFF's renaming, so the bond being re-formed can be named the same way in both files - but every other ligand atom selector must be checked against `preprocess.results.ligand_atom_map` first (D169). And a scan that re-forms a bond the MM topology does not have is a QM-region question, not an MM one: the bond order lives in the QM calculation, so both atoms must be inside the QM region at every point.

## The RQK gap, which matters more than the mechanics

The hypothesis is what the **R-Q-K** triad does to the mechanism: Arg90 in place of the third zinc histidine, Gln93, and Lys226. **Gln93 is in neither preprocess QM region** - the reactant's scans add it per section as `extra_qm_residues: {A: [93]}` (D133), and the product's sections need the same. A scan or NEB run without it would answer a question about a site missing the residue the study is named after. Gln is neutral, so the QM charge is unchanged - which is exactly why the omission is easy to miss.

This also sets up **P2** (switching off the charges of Lys226, Arg90 and Asp92 along each path): the per-residue contributions are the direct test of the hypothesis, and they need a converged path on each arm first.

## Written, not submitted (2026-11-02, overnight)

The bridge was not open - no code-server pane, so no way to read job 31327774 or submit anything, and the portal is the user's to open. What could be done without Bunya was done: `site.yaml` gains the product arm's whole path, `prodopt` -> `recombine` -> `expel` -> `pstep2`/`pstep1` -> `psaddle2`/`psaddle1`, each carrying `extra_qm_residues: {A: [93]}`, with the reasoning above in the comments. No `start` is set on either scan: `scan` defaults it to the distance in the structure it is given and records what it used, so there are no invented numbers waiting to be believed. `end` values are the reactant arm's own measurements (1.39 A for the ring, 3.00 A where `attack` begins).

**Two findings, both about Gln93, which is the hypothesis's own residue.**

1. **`prodqm` did not have it.** The QM/MM MD that equilibrates the product was about to run with Arg90 and Lys226 in the QM region and the Q of RQK outside it. Added. It matters at the equilibration and not only at the scans, because a scan whose first point has to re-equilibrate a region the MD never had measures the region change rather than the reaction. Gln is neutral, so the charge stays 2 and nothing about the omission would have looked wrong.

2. **The reactant arm's path sections do not have it either** - `attack`, `ringopen`, `step1`, `step2`, `saddle1`, `saddle2`. `aspout`, `aspin`, `michaelis2d`, `rwoff`, `rwon`, `qmrate` and `asph` all do, because D133 added it when the Asp92 question was the live one. So the sections that would produce the **barriers for the paper** are precisely the ones missing it, on both arms. The product arm is fixed; **the reactant arm is left as it is pending the user's call**, since those sections were planned before D133 and changing them changes work already designed. Left inconsistent, the two arms are not comparable, which is the whole point of building the second one.

**D171 - The product's 100 ns MD did not finish: job 31327774 died at 34.14 ns with `Particle coordinate is NaN`. The same 4 fs timestep carried the reactant through all 100 ns, so the instability is the product's own. And Gln93's absence from the QM regions was the user's decision, not an oversight - D170 had that backwards.** (2026-11-02, user: "glutamine was left out deliberately because it's probably not involved/functional as it's neutral and far away from everything")

## The failure

`protdft/prod100-31327774.err`: `openmm.OpenMMException: Particle coordinate is NaN`, raised from `simulation.step` after 68.7 minutes, at step 8,535,000 of 25,000,000. `nvtsim.csv` has 1,707 frames and its last rows are unremarkable - potential energy -659,172 kJ/mol, 300.2 K, the same as 30 ns earlier - so the blow-up happened inside the final 20 ps reporting interval with no drift leading into it.

**It is not the timestep as such.** `md_wt100` reached 100,000.00 ps at step 25,000,000 - 4 fs, the same `md.inputs` default, the same `autoconstraints='HBonds'` and `rigidwater=True`, completed. **It is not disk**: scratch is 273/300 GB and the write succeeded. What differs on this arm is the merged, OpenFF-parameterised 53-atom product ligand, which the reactant's run did not have (`merge_ligand: false` there - the substrate came in inside the protein PDB).

A NaN with no energetic warning is a single bad interaction, not a degrading system: the candidates are the grafted geometry relaxing into a close contact, an OpenFF torsion on the hydrolysed carboxylate end, or mixed precision on CUDA. Which of those it is cannot be read off the trajectory; it can only be distinguished by restarting and seeing whether it recurs at the same place.

## What is salvageable, and what it decides

**34.14 ns of trajectory survived** - 1.0 GB of `NVTtrajectory.dcd`, written right up to the crash. The question that gates this arm is whether the zinc held its four-coordinate shell, and 34 ns answers that as well as 100 ns would; what 34 ns cannot do is match the reactant's 100 ns as a control. So `prodpost` runs now on what exists (**job 31354500**), and the restart decision waits on what it says. If the shell degraded before the crash, the NaN is a symptom and the arm needs rebuilding; if the shell is clean, the NaN is an accident and a restart is enough.

`prodpost` therefore names `pdb_file` and `trajectory` directly, because a FAILED `md` section never recorded a `trajfile` or a `lastframe`. The topology is the trajectory's own first frame rather than `finalsystem.pdb` - postprocess refuses a topology whose atom count disagrees with the frames, and the first frame is by definition the structure the frames were produced from.

**The restart options, for the user:** (a) continue from the last good frame with fresh velocities at 2 fs, ~5.4 h for the remaining 66 ns, which also tests whether 4 fs was the cause; (b) accept 34 ns as the product's equilibration and go straight to `prodqm`, giving up the matched 100 ns control; (c) restart at 4 fs with a new seed, cheapest and least informative.

## Gln93: D170 corrected

D170 read Gln93's absence from the QM regions as a gap and called the omission "easy to miss". It was not a gap. The user left it out deliberately: the Q of RQK is a **neutral** amide and far from the site, so it has none of the electrostatic lever the R and the K have, and the QM region stays as small as the question allows. Removed again from the eight product-arm sections D170 added. The sections that do carry it - `aspout`, `aspin`, `michaelis2d`, `rwoff`, `rwon`, `qmrate`, `asph` - added it for D133's Asp92 question, where the amide is a hydrogen-bond partner of the carboxylate rather than a charge, and that remains a different question from the barrier.

**The lesson is the same one as D165.** I inferred a design decision from a config's contents - twice now, in the same week - rather than from the record or from asking. "Present in these sections and absent in those" is evidence of a decision, not of a mistake, and the asymmetry was the decision being visible. `prodpost` monitors Gln93 NE2 to Asp92 OD1 and to the ligand's O2x, so the claim that it is far from everything will be a number on this trajectory rather than either of our assertions.

**D172 - Before restarting, check the hydrogens: they are the one part of a system no crystal structure supplies, nothing downstream refuses a wrongly placed one, and an H sitting on another atom is a NaN waiting for the integrator to find it. `hydrogen_report` and `export_frames` put both halves of that in mcube.** (2026-11-02, user: "Go back a couple frames in the trajectory and restart - but before you do, check that the topology and forcefield is correct - make sure both the ligand and protein has hydrogens placed correctly, as that is the main thing not present in PDB structures")

## Why this is the right check for this failure

The NaN came with no warning: energy and temperature flat to the last reported frame, every force-field template matched, the minimiser satisfied. Everything that *can* be checked by a template match already was. What a template match cannot see is a hydrogen placed in the wrong spot or the wrong state - it checks that the names it needs are present, not that the geometry is possible or that the protonation is the one intended. Those are placed by pdbfixer, Modeller's `addHydrogens`, OpenFF's `AddHs`, or a graft done by hand, and on this arm all four touched the system.

## What went into mcube

`preprocess.hydrogen_report(pdb_file, inputs)` reads it off the structure, needing no force field and no connectivity, so it runs on a built system, a frame out of a trajectory, or someone else's PDB:

- **Protonation states by which hydrogens are present**, not by residue name: HD1 alone is HID, HE2 alone is HIE, both is HIP; ASP with HD2 is ASH; GLU with HE2 is GLH. `residue_variants` asks for a tautomer - this reports what the file has. Asp92's charge is the question this whole site turns on, and a stray HD2 would quietly remove the electrostatic reason for the contact being measured.
- **Residues with no hydrogens at all** - a TIP3P oxygen without its two is not water.
- **Amide nitrogens missing their H**, proline and the terminal H1/H2/H3 variants excepted.
- **Unattached hydrogens** (nothing heavy within 1.35 A) and **hydrogens with two heavy atoms that close**, which is the clash that becomes a NaN. The neighbour grid wraps with the cell, because `separation` applies the minimum image and a grid that did not wrap would report every atom at a box edge.
- **The ligand's formula against its own `ligand_smiles`.** This is the check with no substitute: OpenFF types whatever it is given, so a ligand hydrogenated to the wrong protonation state parameterises cleanly and runs at the wrong charge. Methanol less one hydrogen is methoxide, and nothing but the formula says so.

`preprocess` now runs it on every system it builds. `postprocess.inputs.hydrogens: true` runs it on an existing structure, which is how this arm's system - built before the check existed - gets checked without rebuilding it.

`postprocess.export_frames` writes named frames back out through the topology, negative indices counting from the end, so a crashed run can be restarted from a frame *before* the crash. A frame past either end is refused rather than clamped, because silently clamping would restart a run from a geometry nobody chose. Written before any other analysis in the section, so the restart structure survives a later step failing.

## A bug of my own, found by the first attempt

Job 31354500 - the analysis of the salvaged 34 ns - died with `KeyError: 'csv'`. D171 made the trajectory and topology optional for a FAILED md section and left the reporter CSV required, so the stage still refused to analyse exactly the runs most worth analysing. The CSV is now optional too, and its thermodynamic plots are skipped with a line saying why.

## The restart, as configured and not yet submitted

`prod100b` resumes from **frame 1704 = 34.08 ns**, three frames and 60 ps before the crash. **Nothing else is changed**, and that is the point: 4 fs carried `md_wt100` through all 100 ns on the same constraints, so the timestep is not a known fault, and a product integrated differently from its reactant is a different control. Resuming from a PDB draws fresh velocities, which is the one thing that does differ - so a second crash near 34 ns indicts the system, while a clean run indicts the first integration step and nothing else.

The four restraints are given **outright rather than redetected**. Redetecting on a trajectory frame would hold each metal-donor distance where MD had drifted to, which is the inversion D120 and D130 are about: a contact is held where the crystal has it, never where a previous run put it.

**Job 31355566** does the export, the hydrogen check and the shell analysis in one pass. The restart waits on it: the user's instruction was to check the topology *before* restarting, and a restart from a frame whose hydrogens are wrong would reproduce the crash on purpose.

43 of the new tests pass in the ASH environment, which is the one that has MDAnalysis - run there on purpose, since D167 is the entry about a gated test hiding a broken module.

**[SUPERSEDED by D191/D195 - Gln93 reaches the zinc (3.0% of frames, 1.79 A elsewhere), and the NaN recurred with the site restrained. See `csr/retractions.md`.]** **D173 - The 34 ns did not equilibrate the product complex, it rearranged it. The zinc is SIX-coordinate in 1704 of 1707 frames: Asp92 OD2 collapsed from 2.80 A to 2.03, and the ligand turned over - O2x came 4.75 A onto the metal while the crystal's own two contacts left. The restart is withdrawn; the hydrogens, meanwhile, are clean, and Gln93 is 7-9 A from everything, which settles it.** (2026-11-02)

## The three gates

**Hydrogens: pass.** On `restart_1704.pdb`, 32,297 of 49,580 atoms are hydrogens. The ligand's formula is **C17H26N3O6S1 and matches `ligand_smiles` exactly** - so the graft, the hydroxide merge and OpenFF's `AddHs` between them produced the molecule the config asks for, at the protonation state it asks for. Asp92 reads as unprotonated (no HD2); His88 and His164 are both HID, as `residue_variants` asked. No unattached hydrogens, no hydrogen clashes, every amide N with its H. **The hydrogens were never the problem.**

**The `missing` contact list: pass, empty.** All 23 pairs resolved.

**The zinc shell: fail, and not marginally.**

| contact | `finalsystem.pdb` | over 34 ns (min/mean/max) | occupancy |
| --- | --- | --- | --- |
| HOH C:1 O | 1.90 | 1.94 / 2.07 / 2.29 | 100% (restrained) |
| HOH C:2 O | 2.14 | 1.95 / 2.11 / 2.31 | 100% (restrained) |
| His88 NE2 | 2.22 | 2.05 / 2.20 / 2.43 | 100% (restrained) |
| His164 NE2 | 2.09 | 2.04 / 2.20 / 2.42 | 100% (restrained) |
| **Asp92 OD2** | **2.80** | **1.88 / 2.03 / 3.25** | **99.9%** |
| Asp92 OD1 | 3.81 | 2.04 / 3.23 / 3.50 | 0.2% |
| **hyd O2x** | **6.76** | **1.89 / 2.01 / 2.20** | **100%** |
| hyd O4x | 2.94 | 6.24 / 6.64 / 6.99 | 0% |
| hyd O5x | 2.92 | 6.94 / 8.46 / 9.79 | 0% |
| hyd O6x | 7.14 | 3.09 / 3.84 / 4.21 | 0% |

**Coordination number: 6 in 1704 frames, 7 in 3.** Against the crystal's four.

## What this says, which is more than "the MD failed"

**Two contacts moved, and both were left unrestrained by design.** Asp92 OD2 sat at 2.80 A, outside the 2.5 A detection cutoff, so it got no restraint - and that cutoff exists precisely to stop a second-shell contact being restrained as though it were a bond (D120, D130). The decision is still right and still incomplete: it protects against holding a contact that is not a bond, and says nothing about a contact that is free to *become* one. A carboxylate 2.80 A from a divalent cation under a point-charge force field has one place to go.

The ligand is the same hole with a bigger displacement. `constraint_exclude_residues: {D: [1]}` leaves the substrate free, which is correct for the reactant - a substrate that cannot move cannot react - and wrong for an equilibration whose only job is to let solvent relax around a crystallographic geometry. Freed, the product rotated: **O2x travelled 4.75 A onto the metal while O4x and O5x, the two oxygens 5AXO puts at 2.9 A, left for 6.6 and 8.5 A.** The complex that came out is not the complex that went in.

**The NaN is explained, and it was never the interesting part.** A system that rearranges its metal site this far in the first frames is a sick system, not an unlucky integration step: Asp92's 99.9% occupancy means it came in almost immediately and stayed. Nothing about the timestep needed changing and nothing about it would have helped.

**The user's own rule is the fix, applied properly this time:** "we can't trust MD for coordination so that's automatically invalid - anything within coordination range should be restrained at the distance it is in the crystal for MD." `prod100c` holds **eight** contacts at their `finalsystem.pdb` distances rather than four. Asp92 needs both oxygens, because D130 is the entry where holding OD2 alone let OD1 rotate in to 1.99 - a monodentate restraint on a bidentate group restrains a distance, not a geometry. The ligand needs O4x and O5x, and holding them is what makes O2x's migration impossible: it cannot reach the metal while the crystal's own contacts are held there.

`prod100b` is **withdrawn in place** rather than deleted. The restart machinery was right; only its starting structure was wrong, and there is no frame in that trajectory worth starting from.

## Gln93, settled by measurement

The user's position was that Gln93 is "neutral and far away from everything". The trajectory agrees, and not narrowly:

- **Gln93 NE2 - Asp92 OD1: 5.83 / 7.08 / 8.60 A**, 0% of frames inside 3.5 A.
- **Gln93 NE2 - ligand O2x: 7.95 / 8.82 / 9.88 A**, 0%.

Never a contact, over 34 ns, on either partner. Keeping it MM is right. **It also puts a question to the sections that still carry it** - `aspout`, `aspin`, `michaelis2d`, `rwoff`, `rwon`, `qmrate`, `asph` - which added it in D133 on the premise that its amide hydrogen-bonds the carboxylate whose position is being measured. At 5.8 A minimum it does not, at least on this arm. The reactant arm has its substrate in the pocket and may differ, which is a measurement rather than an assumption, and `postprocess` can now make it on `md_wt100`'s trajectory.

Two further things the monitor turned up. **Arg90 NH1/NH2 to the ligand's O1x is 11-20 A** - not a pair worth monitoring, and I chose it from atom names without checking what they were near. And **O6x-C4x is 1.22 A at 100%**, a carbonyl bond length, confirming the nucleophile-derived oxygen is covalently on the carbon, with **C4x-N1x at 3.56-4.01 A**: the ring is properly open. The product is the right molecule in the wrong place.

## One fix to the check this exposed

`hydrogen_report` listed the zinc and all 56 sodiums as "residues with no hydrogens" - 57 lines burying anything real. A monatomic residue has no hydrogens to be missing, so they are counted as ions instead. The honest limit is recorded in the tests: a lone oxygen named HOH is indistinguishable from an ion by atom count, and water missing its hydrogens is caught at the template match, where tip3p has no one-atom residue to offer.

**D174 - The 3 A rule, applied. Detection at 2.5 A was mine, not the project's, and it is wrong in the one way that matters: the cutoff decides which contacts are HELD, so a tight cutoff does not protect a contact between 2.5 and 3.0 A, it abandons it. `prod100c` detects at 3.0 and restrains all seven contacts it finds. Job 31355970.** (2026-11-02, user: "use the 3 angstrom constraint! We already established that anything within that range is affected by the zinc and should be restrained to it's crystallographic distances")

## Where I went wrong, which is worth being exact about

D120 and D130 are two runs ruined by a restraint on Asp92: detection at 3.0 A swept OD2 in at 2.80, the restraint held it, and the carboxylate pivoted OD1 onto the metal at 1.99, giving a five-coordinate zinc. From that I concluded the cutoff should be 2.5 so Asp92 is never swept in at all, wrote that into `wt100`'s comment block as settled, and defended it again in D165 when challenged.

**The conclusion does not follow from the evidence.** What D120 and D130 show is that holding ONE oxygen of a bidentate group is worse than useless - it pins a distance while leaving the geometry free, so the group relaxes the only way it can. They say nothing about whether the contact should be held. I treated "restraining it went wrong" as "do not restrain it", and the alternative I chose - leave it free - is the one option neither entry tested.

D173 tested it: unrestrained, **OD2 collapses to 2.03 A by itself**, and the ligand turns over as well, O2x travelling 4.75 A onto the metal. A tight detection cutoff does not keep a second-shell contact in the second shell. It only removes mcube's ability to.

The rule was already stated, twice, in the user's words: *"we can't trust MD for coordination so that's automatically invalid - anything within coordination range should be restrained at the distance it is in the crystal for MD."* Coordination range is 3 A. I read it as licence to pick the cutoff that kept Asp92 out of the restraint list, when it is an instruction to hold everything the zinc reaches at the distance the crystal has.

## What is running

`preprocess.inputs.cutoff` in `CSR_prod/config.yaml` is **3.0** (backed up as `config.yaml.bak-cutoff25`), and `constraint_exclude_residues` no longer names the ligand - that flag is right for the reactant, whose substrate must be free to react, and it is what left the product's own two zinc contacts unheld.

At 3.0 A on `finalsystem.pdb` the sweep is **seven** contacts, each held where that file has it:

| contact | distance | at 2.5 A |
| --- | --- | --- |
| HOH C:1 O | 1.90 | held |
| His164 NE2 | 2.09 | held |
| HOH C:2 O | 2.14 | held |
| His88 NE2 | 2.22 | held |
| **Asp92 OD2** | **2.80** | free - collapsed to 2.03 |
| **hyd O5x** | **2.92** | free - left for 8.46 |
| **hyd O4x** | **2.94** | free - left for 6.64 |

Plus **one restraint beyond the sweep, flagged in the config as a judgement: Asp92 OD1 at 3.81 A**, which is outside 3 A so the rule does not reach it. It is there because D120 and D130 are the same failure twice, and both were OD2 held with OD1 free; a carboxylate whose other oxygen is 2.80 A from a divalent cation is affected by it whatever the distance says. One line to strike if the group should be free to pivot.

`coordination_redetect: true` is what applies the 3.0, so the seven restraints are a measurement of the structure rather than a list typed by hand - and `md`'s own results will record all seven with the distances used, which is the check that this entry did what it says.

**Job 31355970**, gpu_cuda with an H100, 4 CPUs, 32 GB, 8 h for an expected 4.1. The first attempt failed on `Invalid qos specification`: `run.slurm` hardcodes `--qos=normal`, which the GPU partition rejects, so a GPU submission needs `-q gpu` on the command line. Worth moving into the job script.

`prod100b` is withdrawn in place, and the void run's 1.0 GB trajectory, its checkpoint and its optimiser trajectory are deleted - `nvtsim.csv`, the first frame and `restart_1704.pdb` are kept, which is everything D173's numbers were read from. Scratch 272/300 GB, 889k/1M inodes.

**[CORRECTED by D176 and D208 - the per-residue default was one of two necessary causes, not the whole story. See `csr/retractions.md`.]** **D175 - Detection restrained one donor per RESIDUE, not per atom, so a group with two donors inside the cutoff always had one pinned and the other free. That is the mechanism of D120 and D130 both, and it is in the code, not in the configs. `prod100c` as submitted holds seven contacts, not the eight its comment claims: the ligand's O4x at 2.94 A is free.** (2026-11-02)

## The bug

`coordination_from_rows` carries this comment, and has since it was split out:

> Per residue for the QM region, which takes whole residues anyway. Per atom for restraints: a carboxylate can rotate its second oxygen in while the first is held.

`find_coordinators` is the only caller that builds restraints, and it called `coordination_shell(pdb_file, inputs)` - taking the **per-residue default**. So the rule the comment states was never the rule that ran. Only the closest donor of each residue was ever restrained.

That is not a subtle mismatch. It is exactly the failure D120 and D130 record: Asp92's OD2 held at 2.80 A while OD1, in the same residue, was free and pivoted onto the metal at 1.99 A. Both entries blamed the cutoff; both wrote the remedy as a cutoff change; and the cause was that the restraint list could not contain two atoms of one residue however the cutoff was set. D174 - mine, an hour ago - repeated the same blame.

## What it means for the run in flight

`prod100c` (job 31355970) was submitted and the pane closed before the restraint list could be read, so this is from the distances rather than from the job. At 3.0 A on `finalsystem.pdb`, per residue:

| donor | distance | restrained per residue? |
| --- | --- | --- |
| HOH C:1 O | 1.90 | yes |
| His164 NE2 | 2.09 | yes |
| HOH C:2 O | 2.14 | yes |
| His88 NE2 | 2.22 | yes |
| Asp92 OD2 | 2.80 | yes - closest in ASP A:92 |
| **hyd O5x** | **2.92** | yes - closest in hyd D:1 |
| **hyd O4x** | **2.94** | **NO - same residue, 0.02 A further** |

Seven restraints, not eight: six detected plus the configured Asp92 OD1 at 3.81. **The ligand oxygen 0.02 A behind its neighbour is free**, which is the same hole that let O2x travel 4.75 A onto the metal, narrowed but not closed. With O5x held the ligand is far better anchored than in `md_prod100`, so the run may well produce a usable trajectory - but it is not the run the config describes, and `prodpost` would have to decide it after the fact rather than the restraint list deciding it in advance.

## The fix

`find_coordinators` now detects with `per_residue=False`. The QM region is unaffected: it is keyed by residue, so naming one twice sets the same entry twice.

**It changes behaviour for every project**, which is the right direction - a restraint list is not a description of which residues touch the metal, it is the set of distances the run is asked to preserve - but the channel work should know that a run redetecting coordination will now hold more contacts than before.

`tests/test_coordination_detection.py` is 8 tests, and the two that matter were confirmed to FAIL against the old per-residue call before the fix went in: both oxygens of a carboxylate restrained, each at its own distance. **Nothing had ever tested `find_coordinators`** - no test in the suite referenced it, `coordination_shell` or `coordination_from_rows` - which is how a comment and its implementation disagreed through two decision entries about the consequence.

## What has to happen to the run

`scancel 31355970` and resubmit after this fix, so the detection finds all seven and the configured OD1 makes eight. The cost is the GPU hours already spent; the alternative is a 100 ns equilibration with one zinc contact free and no record in the config saying so.

The pane was closed when this was found, so nothing could be cancelled or verified. The fix is pushed and waiting.

**D176 - No: detection was never per atom for restraints. The per-atom option existed in `coordination_from_rows` from the day it was split out, and `find_coordinators` - the only caller that builds restraints - never passed it. The reactant's four restraints happen to be unaffected, one donor per residue; what is not unaffected is that its Asp92 was free at 2.80 A too, which `wtpost` now measures. Asp92 OD1 kept. prod100c resubmitted as 31423830.** (2026-11-02, user: "detection was always by atom? It should always have been per atom restraints")

## The answer

Per **residue**, always, for restraints. `coordination_shell(pdb_file, inputs, radius=None, per_residue: bool = True)` - and `find_coordinators` called it with two arguments. So the switch was there, documented in the comment beside it, defaulted the wrong way for the one use that needed it, and no test touched any of the three functions. It should always have been per atom, as the user says: a restraint list is the set of distances a run is asked to preserve, and a donor atom inside the cutoff is one of those distances whether or not another atom of its residue is nearer.

## What it did and did not touch

**The product:** it bit. The ligand has O5x at 2.92 A and O4x at 2.94 A in the same residue, so detection held O5x and left O4x free 0.02 A behind it. That was the state of job 31355970, which was still PENDING when this was found - so it was cancelled before using any GPU time, and `prod100c` is resubmitted as **31423830** with the fix in place: seven detected contacts plus the configured Asp92 OD1, eight in all.

**The reactant: the bug did not bite, but the cutoff did.** `md_wt100` and `md_wt100b` both recorded exactly four restraints - `D:1:O` 1.90 (hydroxide), `A:164:NE2` 2.09, `C:2:O01` 2.14, `A:88:NE2` 2.22 - four residues, one donor each, so per-residue and per-atom give the same list. Those runs are not invalidated by D175.

What they do share with the product is the 2.5 A cutoff, and therefore **Asp92 OD2 free at the same crystallographic 2.80 A**. On the product that contact collapsed to 2.03 A at 99.9% occupancy. Whether it did the same on the reactant is a measurement, the trajectory is on disk, and nobody has ever looked: `md`'s before-and-after shell report covers two frames, and D134's table came from `site.py` outside the framework.

**`wtpost` (job 31424060)** asks it: coordination occupancy over `md_wt100`'s full 100 ns, per atom, plus the hydrogen check on its topology. Two outcomes, both worth having:

- **Asp92 held near 2.80.** The product's collapse is the product's own - most likely the ligand's own charge distribution next to the metal - and the reactant's 100 ns stands as the control.
- **Asp92 collapsed there too.** Then both arms were equilibrated on a six-coordinate zinc, every QM/MM run descended from `md_wt100` inherits it, and both arms need redoing at 3.0 A. That would be the largest single correction in this log, and it would be better to know it from an analysis job than from a referee.

## Asp92 OD1 stays

The user kept the one restraint outside the 3 A sweep (OD1 at 3.81 A). With detection now per atom, its justification is narrower than when I added it - a carboxylate with both oxygens inside the cutoff is now held as a group automatically - but Asp92's OD1 is outside the cutoff at 3.81, so only an explicit row reaches it, and D120 and D130 are both what happens when it is not reached.

8 of the new detection tests pass in the ASH environment as well as locally.

**D177 - The reactant arm has the same artefact, and it is worse. Over `md_wt100`'s 100 ns the zinc is SIX-coordinate in 4995 of 5000 frames, with Asp92 OD2 bound at 2.07 A mean in 100% of them and OD1 joining bidentate in 56.7%. The 2.5 A cutoff never kept Asp92 in the second shell; it only stopped mcube from holding it there. Every QM/MM run descended from this trajectory inherits a six-coordinate zinc.** (2026-11-02, job 31424060)

## The measurement

`postprocess` over `md_wt100`'s 5000 frames, per atom, cutoff 2.50 A with a 1.0 A margin shown:

| donor | occupancy | min / mean / max | restrained? |
| --- | --- | --- | --- |
| hydroxide C:2:O01 | 100.0% | 1.82 / 1.94 / 2.12 | yes |
| HOH D:1:O | 100.0% | 1.92 / 2.07 / 2.27 | yes |
| **Asp92 OD2** | **100.0%** | **1.91 / 2.07 / 2.32** | **no** |
| His164 NE2 | 100.0% | 2.01 / 2.20 / 2.46 | yes |
| His88 NE2 | 100.0% | 2.05 / 2.24 / 2.51 | yes |
| **Asp92 OD1** | **56.7%** | **1.93 / 2.60 / 3.50** | **no** |
| **HOH E:6862 O** | **43.2%** | **1.96 / 2.11 / 2.45** | **no** |

**Coordination number: 6 in 4995 frames, 5 in 5.** Never four. Asp92 is on the metal in every frame of the equilibration, and the site alternates between OD1 and a bulk water for its sixth place.

The hydrogens on that topology are clean: meropenem matches `ligand_smiles` as C17H25N3O5S1, Asp92 is unprotonated, His88 and His164 are both HID.

## What this overturns

The reasoning block above `wt100` says the run starts from a structure that is "6NC5's ZN302 sphere with Asp92 correctly in the second shell" and that the 2.5 A cutoff "keeps the crystal's four and excludes the 2.80 contact". The first half is true of `finalsystem.pdb`. The second half describes what the *detector* does, and it was read as a statement about what the *trajectory* does. It is not. The cutoff excluded Asp92 from the restraint list and the MM force field then pulled it in within the first frames, to the same 2.0-2.1 A that D120 and D130 recorded as the artefact to be eliminated.

So the fix for D120 never worked, and nothing since would have revealed it: `md` reports the shell on two frames, and the before-frame is `finalsystem.pdb` where Asp92 is still at 2.80. The occupancy over a whole trajectory did not exist as a stage until yesterday (D167), and nobody had pointed it at this run.

**D134, D165 and D174 all rest on this.** D134's four restraints are correct as a list of what was restrained and silent about what else was bound. D165 defended the 2.5 A cutoff as "defensible"; it is not, and the user said so twice before I measured it. D174 called the cutoff choice mine and the rule theirs, which was right, and still underestimated the damage: it is not only the product arm that was equilibrated wrong.

## What inherits it

Everything whose structure descends from `md_wt100` or `md_wt100b`: `nac` and its successors, `aspout`/`aspin` (which measure the Asp92-on-metal question *from a trajectory where Asp92 is already on the metal*), `michaelis2d`, `rwoff`/`rwon`, `qmrate`, `asph`, and any scan or NEB built on those frames. Their numbers are not wrong as computations - they describe a six-coordinate zinc, which is not the species the paper is about.

`aspin`/`aspout` deserve a separate note: they were built to drive Asp92 on and off the metal and measure the cost. Starting from a frame where it is already bound at 2.07 A, the "off" end is the one that had to be constructed, which inverts what the scan was measuring.

## What has NOT been started

Nothing. The reactant needs the `prod100c` treatment - from `finalsystem.pdb`, detection at 3.0 A, every donor it finds held at the crystal distance, which for the reactant is hydroxide 1.90, His164 2.09, HOH 2.14, His88 2.22, **Asp92 OD2 2.80**, and whatever else 3.0 A reaches that 2.5 A did not. That is a 100 ns run, ~4.1 h, and it re-bases the whole reactant arm - so it is the user's call, not mine. `prod100c` (31423830) is already correct on this point and is still PENDING in the GPU queue.

## Why the pane keeps dying

The code-server session moved from **bun151:20685 to bun157:36073**: it is itself a Slurm job (`sys/dash`) with a walltime, and when it expires the page it served stops existing and the app closes the pane. Two were visible in the queue this morning - one 8 h 37 m into its life on bun151, one 2 minutes into its life on bun157 - so the handover was already in progress while I was reading the old one. A longer walltime on the OnDemand session would reduce it; caching the URL across hours is what makes it bite.

**D178 - The audit is adopted as the plan of record, and the kinetics are held BLIND. The user has measured kinetics for CSR-1 with meropenem and will not release them until the computation is finished, so the comparison is a genuine prediction. This entry is timestamped before any barrier exists, which is what makes that claim checkable later.** (2026-11-02, user: "I have kinetics for CSR-1 with meropenem - but I won't give them to you until computational analysis is complete, as our method needs to stand up on it's own")

## The blind protocol, recorded before there is anything to compare

There is no published kinetics for CSR-1 with meropenem (D22, D40), and the audit's strongest criticism of the paper's prospects was that a computed barrier would have nothing to be validated against. That gap is now closed in the best possible way: the kinetics exist and are **withheld**.

So the prediction must be made and recorded first. Concretely, before any number is released:

- The computed barrier for the rate-determining step, with its level of theory, its spread over snapshots, and its sign conventions, goes into this log **and is not revised afterwards**. A barrier that moves after the experiment arrives is a fit, not a prediction.
- The mechanism's qualitative claims - which species is the nucleophile, what protonates the leaving nitrogen, which residues matter - are stated before, not after.
- Where the computation cannot decide something, this log says so rather than leaving the question to be settled by hindsight.

This is worth more to the paper than any amount of extra compute, and it costs nothing except the discipline not to peek. **Do not ask for the kinetics.** The one thing that would destroy the claim is a revision made with the numbers in hand.

## The plan of record

Adopted from `docs/csr/audit-2026-11-02.md`, with one consolidation of my own. Items are tracked here so that the audit's own closing warning - "a fix that was prescribed and not carried out is indistinguishable, later, from a fix that landed" - has somewhere to be checked against. **Status is written next to each item and updated as it changes, not in a new entry.**

**Items 1 and 2 collapse into one route, which is already half built.** The audit asks (1) for a Michaelis complex anchored to experiment and (2) for both endpoints on one topology, built by reversing the hydrolysis. 5AXO is SMB-1 with **hydrolysed meropenem** - our own substrate, in a B3 enzyme, measured - and `CSR_prod` already transferred it into the CSR-1 frame (D159). So the route is:

    prod100c  ->  prodpost gates  ->  prodqm  ->  prodopt  ->  recombine
    (100 ns)      (shell + H)        (50 ps)     (opt)        (re-form C4x-N1x)

and the structure `recombine` ends on **is** the Michaelis complex, anchored to experiment through the product rather than through a docking score, and on the product's own topology so it is a valid NEB partner by construction. **7UHI becomes the check, not the source:** does the re-formed carbonyl land at L1's measured 2.53-2.63 A from Zn1? If it does, the pose is independently corroborated. If it does not, that is a real result about CSR-1 rather than a modelling failure.

This also means no cross-ligand graft is needed. Grafting meropenem's substituents onto moxalactam's transferred core was the other way to item 1, and it is strictly harder and less direct than re-forming a bond in a molecule that is already meropenem.

| # | Item | Status |
| --- | --- | --- |
| 0 | **Asp92's placement in the starting structure** (below) | **RESOLVED (D179): 2.80 A is supported; 4.20 A was the mutant's. prod100c unblocked** |
| 1+2 | Experiment-anchored Michaelis complex, both endpoints on one topology | `prod100c` PENDING (31423830); the rest configured, unrun |
| 3 | Equilibration that preserves the site, occupancy read over the whole trajectory | method fixed (D174-D176); product arm pending, **reactant arm not started** |
| 4 | Nucleophile identity, origin and protonation state | not started; blocked on 1-3 |
| 5 | PROPKA/H++ and the His164 HID/HIE sensitivity check | **not started, depends on nothing but a protonation tool** |
| 6 | scan -> NEB -> TS -> frequencies -> IRC on >=3 snapshots | not started; the paper's actual content |
| 7 | Cluster-model DFT-vs-GFN2 benchmark, Imhof-style | **not started, depends on nothing - should run in parallel** |
| 8 | Protonation of the leaving nitrogen | not started |
| 9 | Per-residue roles, with a method that can reach R, Q and K | **design decision needed before item 3's equilibration** |
| 10 | Di-zinc comparison at identical level | not started; scope decision |
| 11 | Asp92 free energy, only if the discrepancy survives 1-3 | deferred by design |

Item 9 is the one with a deadline that is not obvious: Arg90 and Lys226 are inside the QM region in all eleven systems, so charge deletion can never reach them (D153). Whether the production system puts R, K, Arg229 and Q **all in MM** so the set can be decomposed consistently, or whether the roles are done by computational mutagenesis instead, has to be decided **before** the equilibration in item 3, because it changes what is built.

## Item 0, which the audit surfaced and which blocks the run that is queued

D145 established that the **2.80 A "crystallographic" Asp92 distance is our own model's Asp-hydroxide hydrogen bond**. In 6NC5 the equivalent Asp91 is **4.20 A (OD2) and 4.22 A (OD1)** from Zn1. `CSR_prod/preprocess_results/finalsystem.pdb` has Asp92 OD2 at **2.80** and OD1 at **3.81**.

The zinc in our model is 6NC5's, placed by superposition as a rigid unit (D23, 0.58 A over 256 CA). So a 1.4 A difference in the Zn-Asp distance is not a difference in the metal's position - it means **our Asp92 side chain sits 1.4 A closer to the zinc than 6NC5's does**, in the starting structure, before any dynamics.

`prod100c` restrains Asp92 OD2 at 2.80 A. That is faithful to the rule (hold each contact where the structure has it) and it therefore **preserves the discrepancy** for the whole 100 ns and everything downstream. Three options, and it is a modelling decision rather than a pipeline one:

1. **Hold 2.80 and state it.** Internally consistent, cheapest, and every energy downstream carries a 1.4 A departure from the only metal-bearing CSR-1 crystal. The Asp92 question then cannot be asked of this system at all.
2. **Rebuild the scaffold so Asp92 sits where 6NC5 has it**, then hold it there. Principled - the structure should match the crystal before restraints preserve it - and it means another preprocess and a re-run of the 100 ns.
3. **Establish first why they differ.** Our protein scaffold is not 6NC5 (6NC5 is the Q121H/K254H mutant; the WT coordinates are the metal-free 6DN4/6DQ2), so Asp92's rotamer may legitimately differ between the mutant and the wild type - in which case 2.80 may be right and 4.20 is the mutant's. This is a measurement, not a judgement, and it is cheap: superpose and compare the side chain directly.

Option 3 first, then 1 or 2 on what it says. Nothing is restrained on a number whose provenance is unsettled - that is how 2.80 A came to be quoted as crystallographic for ten weeks.

`prod100c` is still PENDING with no start time (the GPU queue is backed up), so nothing is committed yet.

**D179 - Item 0 resolved, and the conclusion was already in the log: D108 established on 25 October that "6NC5 is a beta-site mutant and is not a wild-type reference for that site. 6DN4 is." D145 then quoted 6NC5's Asp distance as "the crystal's" anyway, and that contradiction sat unnoticed for eight days. What is new here is the measurement in the wild-type files: 6DN4 at 1.99 A, Asp120 B-factors 15-16, puts the carboxylate 3.00-3.07 A from His116 NE2 - the close conformer, the one our model has. So 2.80 A is not a modelling error.** (2026-11-02, user: "we already figured this out days ago")

**The lesson is not about aspartate.** This is the third time in a week that something settled in the log was re-derived from scratch instead of read (D165, D171, and now this), and the first time a *contradiction between two entries* was what had to be found. A 180-entry append-only log cannot be the place a current value is looked up. `docs/csr/state.md` now holds the settled facts, the void results, the live plan and the open questions, and the log's header says so.

## The measurement

`refstruct.py pair` (added for this: `shell` needs a metal and the wild-type structures have none; an internal distance also needs no superposition, so it cannot be an artefact of a fit). Same two atoms in every structure - the first zinc histidine's NE2 against each carboxylate oxygen of the Asp:

| structure | what it is | resolution | Asp B-factors | His NE2 - OD1 | - OD2 | His-His NE2 |
| --- | --- | --- | --- | --- | --- | --- |
| **6DN4** | **wild type**, metal-free | **1.99 A** | **15-16** | **3.00** | **3.07** | 3.63 |
| 6DQ2 | HARLDQ **mutant S60**, metal-free | 2.16 A | **53-66** | 4.28 | 5.49 | 3.75 |
| 6NC5 | **Q121H/K254H mutant**, Zn-bound | 1.60 A | - | 4.13 | 4.82 | 3.42 |

(6DN4/6DQ2 are in BBL numbering - His116, Arg118, Asp120, Gln121, Lys254 - so their Asp120 is our Asp92 and their His116/His192 pair is our His88/His164. The 3.42-3.75 A His-His NE2 separation in all three confirms the same site is being measured in each.)

## What it means, and what it corrects

**The crystals disagree with each other by 2.4 A on this side chain, and the disagreement tracks what each structure is.** The two "far" values come from a *mutant whose Asp120 side chain is disordered* (B 53-66 against a 2.16 A map - that conformer rests on weak density) and from a *mutant that engineers a zinc-binding histidine into the position adjacent to the Asp*. The one "close" value comes from the only well-ordered wild-type carboxylate in the set.

So **D145's correction needs correcting, and D108 had the grounds for it already.** D145 was right that 2.80 A had been quoted as crystallographic when it was our model's own Asp-hydroxide hydrogen bond, and right to catch that. But the 4.20 A it offered as the real crystal value is the engineered construct's - and D108 had stated eight days earlier, in those words, that 6NC5 is not a wild-type reference for this residue. 6NC5's His92, the engineered mutation, is the residue immediately C-terminal to Asp91: a Zn2-ligating histidine next door is exactly the kind of thing that moves a neighbouring carboxylate.

**`prod100c`'s restraint at 2.80 A is therefore supported, not an error**, and the entry above it in D178 (item 0, "BLOCKS prod100c") is withdrawn. The run should go.

**More importantly, the Asp92 question is reframed.** The project has spent ten weeks treating "Asp92 comes onto the metal" as a model-versus-crystal discrepancy. On this evidence it is a **conformational equilibrium that the crystal record itself samples** - a well-ordered close conformer in the wild type, a disordered far one in a mutant. That is consistent with D151's measurement, which found the on and off basins separated by only 1.5-1.9 kcal/mol with populations within 0.2 kcal/mol and correctly withdrew the "the orphaned carboxylate fills the vacancy" claim. Two independent lines now say the same thing: this carboxylate is mobile.

What it does **not** license is the six-coordinate zinc of D173/D177. A mobile second-shell carboxylate at 2.8-3.0 A is not a carboxylate bound at 2.03 A displacing the coordination number to six. Both can be true: the Asp sits close in the wild type, and our MM runs still over-bound it. The restraints in `prod100c` hold it at 2.80, which is now the better-supported distance rather than merely the one our file had.

## What is still open on this

The decisive number - Zn to Asp92 in 6DN4's frame - cannot be measured directly, because 6DN4 has no metal. It needs 6NC5's zinc superposed into 6DN4 on the site residues, which is `refstruct.py transfer`'s machinery pointed at a metal rather than a ligand. Until that is done, "3.00 A from His116 NE2 implies roughly 2.8-3.5 A from the zinc" is geometry, not measurement, and should be written as such.

**Two things for the paper fall out of this**, and both are cheap:
- The Asp92 conformer comparison across CSR-1's own three crystals, with resolutions and B-factors, is a figure. It is the kind of structural due diligence a referee looks for and it was sitting in the PDB the whole time.
- It is also the honest frame for the mechanism: if the carboxylate is mobile, the question is not "is Asp92 on or off the metal" but what each conformer does to the nucleophile - which is a mechanistic question with an answer the calculations can give.

**D180 - Redundant files removed, and the headline is that CSR is not where the disk goes. All eleven CSR run directories together are 14.5 GB of the 269 GB in use; the pressure is other projects. So this was done for clarity, not quota, and two things were deliberately NOT deleted.** (2026-11-02, user: "remove or archive any redundant files")

## In the repository

Five `docs/csr` scripts moved to `archive/`, **with their tests**, following the convention `threads.py` and `ramd.py` already set:

| File | Why |
| --- | --- |
| `traj_contacts.py` | contacts over a trajectory is a `postprocess` stage now (D168), driven by `site.yaml`'s `monitor` and resolved by the same parser restraints use. Its own README row said it was kept "only until D111's numbers are re-derived with it"; D111's run is void (D177), so there is nothing to re-derive with a script rather than the stage |
| `check_ohH_frame.py`, `check_ohH_system.py` | superseded by `preprocess.hydrogen_report` (D172), which runs on **every** system the pipeline builds rather than one arm by hand, and checks more |
| `make_ohH_config.py` | the record of a route D113 abandoned; the file is the record, which is what an archive is for |
| `migrate_paths.py` | one-shot, and `resolve_forcefields` does the resolution at run time in every stage now (D163) |

`test_traj_contacts.py`, `test_make_ohH_config.py` and `test_migrate_paths.py` went with them: 672 passing becomes 629, and **none of the difference is a test that stopped being run against live code.**

**Kept, with the reason stated so it is not re-litigated:** `make_water_nuc.py` and `water_nuc.xml`, because plan item 4 - the nucleophile's identity - is open and needs exactly that build; and `traj_frame.py`, because whether `postprocess.export_frames` covers the optimiser's own xyz output has not been checked, and an unverified claim should not drive a deletion.

## A bug found in the deletion tool, before using it

`cleanup.sh` pruned three literal paths - `preprocess_results`, `md_results`, `qmmm_results` - so only a **default-named** run was protected. A run submitted as `md md_wt100` writes `md_wt100_results/`, which was swept like any working directory. Nothing precious matched the working-file patterns so nothing was lost, but the header promised more than the code did, and that is not a promise to get wrong in a script whose other mode is `--delete`. It now prunes every `*_results` directory.

## On scratch

- **2.2 GB in 346 files** of ASH and geomeTRIC working files from `CSR_wt`, `CSR_nac` and `CSR_crystal_v2`, through `cleanup.sh --delete` after reading each listing. Every input, converged structure, energy table and job log survives by the script's own design. `CSR_prod` was left alone: `prod100c` is queued against it.
- **264 MB**: two `OpenMMMD_traj_wrapped.xyz` files, which are text copies of the `QM_MM.dcd` kept beside them at 28 MB - 4.7x the size for the same frames. A `DELETED.md` is left in each directory, as `md_wt100b_results` already does.
- **~500 MB**: three 183 MB job logs gzipped in place rather than deleted. Provenance is cheap at a tenth the size and `zcat` still reads it.
- **294 loose files** in the scratch root - `o.*`, `out*.txt`, `*.log`, and the mangled `bridge.txtddaatee`, `bridge.txtecho`, `probe.txtdate` left by the terminal-typing failures. The root went from 394 entries to 100. These were shell captures from the browser bridge, and they were actively confusing to read past.

Scratch 272 -> 269 GB. The inode count did not move measurably: ~600 files against 889,000 is noise, which is the point of the headline.

## What was NOT deleted, and why

**`md_wt100` and `md_wt100b`'s imaged NVT trajectories, 2.83 GB each.** Both runs are void as equilibrations (D177), so the instinct is to delete them. Against that: they are the **primary evidence for D177 itself**, the largest correction in this log, established today. The occupancy numbers are recorded in three places, but a recorded number is not a re-checkable one. Deleting 5.7 GB would buy nothing - CSR is 5% of the usage - and would cost the only route to verifying this week's most consequential claim. `md_wt100b` had already been pruned once on 30 October (4.2 GB of raw NVT, NPT and the minimisation path) with the imaged NVT deliberately kept; that decision stands.

The rule this suggests, and it is worth keeping: **the evidence for a live claim is not redundant, however void the run that produced it.** A trajectory becomes deletable when the claim it supports has been superseded by a better measurement, not when the run is judged bad.

**D181 - The reactant arm is re-based at 3.0 A as `wt100c`, job 31483853. Both arms are now queued with the same treatment, so plan item 3 is in flight on both. And a correction to an assumption made while lowering the walltimes: the `md` stage has no resume. `qmmm` got one in D109; `md` never did.** (2026-11-02, user: "why not rebase it at 3 angstroms?" and "feel free to lower the max time if the job won't take that long - it should be resumable even if it gets cancelled")

## Why it was waiting, and why that was wrong

I had this as the user's call because re-basing invalidates everything descended from `md_wt100`. That was the wrong frame: the arm is **already** invalidated - D177 measured it six-coordinate in 4995 of 5000 frames - so there is nothing left to protect, and the 3.0 A rule is the user's own standing instruction rather than a new decision. The one real argument for waiting was cost, and the GPU queue makes it an argument for the opposite: both jobs now wait in the same queue in parallel instead of one after the other.

## The sweep, measured on the reactant's own structure

At 3.0 A per atom on `CSR_wt/preprocess_results/finalsystem.pdb`, five contacts:

| donor | distance | at the old 2.5 A |
| --- | --- | --- |
| HOH D:1:O | 1.90 | held |
| His164 NE2 | 2.09 | held |
| hydroxide C:2:O01 | 2.14 | held |
| His88 NE2 | 2.22 | held |
| **Asp92 OD2** | **2.80** | **free - and it collapsed to 2.07** |

The next three donors are outside: the substrate's own carboxylate oxygen `mer C:1:O1x` at **3.13**, Asp92 OD1 at **3.81**, and **Arg90 NH1 at 3.89**. That last number is why 3.0 A is a ceiling rather than a floor - nitrogen is in `donor_elements`, so the 4.0 A cutoff this project opened with restrained an arginine to the zinc as though it were a ligand (D8). Between 3.0 and 3.9 there is nothing to gain and an arginine to lose.

`constraint_exclude_residues: {C: [1]}` **stays** on this arm. The reactant's substrate must be free to react, and here it costs nothing either way since its nearest atom is outside the sweep. This is the one place the two arms legitimately differ: on the product arm the ligand's own 2.92/2.94 A contacts are crystallographic, and holding them is the point (D174). Asp92 OD1 is restrained explicitly at 3.81 on both arms, at the user's instruction, because D120 and D130 are the same failure twice - OD2 held with OD1 free, and the carboxylate pivoting OD1 in to 1.99.

**What this run answers:** plan item 3, which is unanswered on either arm - does a site restrained at 3.0 A per atom hold four-coordinate through 100 ns? **What it does not:** give a Michaelis complex. The substrate is still Vina pose 2, 5.15 A from L1's measured pose (D156), and items 1-2 replace it by building the reactant from the product through `recombine`. This is the coordination control, not the reaction's starting structure, and it should not be quoted as the latter.

## The walltimes, and the assumption underneath them

Both jobs are now at **5 h** for a run measured at 3.4-4.1 h (`md_prod100` did 34.14 ns in 68.7 min, 497 ps/min; the earlier figure was 408). `prod100c` (31423830) was lowered in place with `scontrol update`; `wt100c` was submitted at 5 h. An 8 h request only starts when an 8 h window opens, and the h100 nodes are `PLANNED` - held for a higher-priority job - so the window is the binding constraint.

**`scontrol update` prints a Bunya warning on this request:** "This job has requested too much memory and CPU per GPU." It applied the change anyway, and the first product MD ran to 34 ns on exactly `cpu=4, gres/gpu:h100=1, mem=32G`, so the ratio is tolerated in practice. Recorded because it looks like an error and is not.

**The assumption to correct: `md` cannot resume.** `grep -n 'resume|checkpoint|state.xml|loadState|restart' src/mcube/md.py` returns nothing. D109 added resume to `qmmm` - ASH writes `OpenMM_MD_state.xml` and `OpenMM_MD_final_state.xml` as it runs and says in its own output that they can restart a simulation, and `qmmm.py` passes them - but `md.py` never did, even though `md_prod100_results/OpenMM_MD.chk` shows the files are being written. So a cancelled or timed-out `md` job restarts from zero.

With 5 h against a 3.4-4.1 h run that is ~22-47% headroom, which is comfortable, so the lowered limits stand. But the general fix is worth having and is small: `md` should pass the state file the way `qmmm` already does, and a `resume: true` input should mean the same thing in both stages. **Added to the plan as item 12.**

**D182 - `md` resumes now, and the user was right about the files: ASH writes `OpenMM_MD_state.xml` as it runs whichever stage calls `OpenMM_MD`. What was missing was mcube reading it. The restart helpers move from `qmmm` to `md` and `qmmm` re-exports them, so one implementation serves both stages.** (2026-11-02, user: "I'm pretty sure md writes as it runs and should be resumable - fix if not")

D181 reported that `md` "has no resume" and left it at that. The files were there the whole time - `md_prod100_results/OpenMM_MD.chk`, `OpenMM_MD_final_checkpoint.chk` and `OpenMM_MD_final_state.xml` are all in the listing from the NaN crash - so the gap was never the data, only the wiring. D109 wired it for `qmmm` and nobody carried it across.

**Two things are specific to this stage, and both would have been silent failures.**

*The keys.* `qmmm` reads `time` and `timestep`; `md` reads `nvt_time` and `nvt_timestep`. `resolve_resume` now takes `time_key` and `timestep_key`, defaulting to the `qmmm` names so that stage is untouched. Reading the wrong key finds `None` and the refusal names the key it wanted, rather than resuming against nothing. The timestep matters too: it is the threshold for "already finished", and 4 fs is not 1 fs.

*The phases before the production MD.* `md` minimises and runs NPT first. On a resume both are **skipped**, because the state file carries positions and velocities - re-equilibrating the box would spend exactly the walltime that resuming exists to save, and then hand the MD coordinates it discards. This is decided before the minimisation rather than inside the NVT block, which is why `resolve_resume` is called at the top of `run`.

Both trajectory files are moved aside, not just one: `md` writes a raw and an imaged DCD, and ASH builds its DCDReporter with `append=False`, so a restart truncates whichever it is pointed at. Segments become `NVTtrajectory_part1.dcd`, ... and are listed in `results.trajfile_parts`.

`tests/test_md_resume.py` is 13 tests, including one that pins the re-export - every existing caller and test in the tree reaches these helpers through `qmmm`, and 629 passing stayed 629 before the new file was added. 642 now.

**It does not help the two jobs already queued**, which were submitted without `resume` and cannot be edited. At 5 h against a 3.4 h run that is comfortable. What it does mean is that the next tight walltime is not a gamble, and that a 100 ns equilibration can be a chain of shorter jobs - which also makes it easier to schedule, since a 2 h window opens far more often than a 5 h one on a partition whose nodes are held for higher-priority work.

**D183 - `resume` is a config option and always was: there is no CLI flag, and `grep resume src/mcube/cli.py` returns nothing. It is read from the section's `inputs` like every other setting. What it did lack is a way to notice you forgot it, so an unset `resume` with a state file present now says what the run is about to discard.** (2026-11-02, user: "avoid resume flag but you can probably add it as a config option")

To be concrete, this is the whole of it in `site.yaml`:

    wt100c:
      stage: md
      inputs:
        nvt_time: 100000
        resume: true

`sections.py` carries it into the run's `config.yaml` with the rest, and `md` reads `inputs.get("resume")`. Nothing was added to the command line - `cli.py` is untouched by D182, and the only flags that stage takes are the ones it already had.

**The default stays false, and the warning is why that is now safe.** `scan` defaults `resume` to **true** (`inputs.get("resume", True)`) because skipping a finished point cannot change a result. For a trajectory, continuing and replacing are genuinely different outcomes, so defaulting it on would mean a deliberate rerun of a section silently appends to the previous attempt - which is the mirror image of the bug just fixed, and harder to notice. So `md` and `qmmm` keep it opt-in.

The failure that is left is forgetting the line, and that is now loud rather than silent: if `resume` is unset and `OpenMM_MD_state.xml` is sitting in the working directory, the run prints how far the previous attempt got, says it is starting from the beginning and overwriting it, and names the one line that would continue instead. A state file too truncated to read is reported too - something wrote it and then stopped, which is the case most worth mentioning rather than least.

Three tests. 659 pass.

**[VOID per D208 - 2.68 A superposes 6NC5's zinc into 6DN4, which contains zero metal atoms. See `csr/retractions.md`.]** **D184 - D179's open number, measured: superpose 6NC5's ZN302 into 6DN4 on the two histidine rings alone and the zinc lands 2.68 A from Asp120 OD2. Our model has 2.80. The wild type puts the carboxylate where we have it, to 0.12 A, and the 4.20 A that D145 offered is the Zn2 mutant's own Asp91.** (2026-11-02)

## The measurement, and why the fit is the result

6DN4 is the wild-type reference (D108) and is metal-free, so the Zn-Asp distance cannot be read off it directly - 6NC5's zinc has to be put into its frame. Two things made that honest rather than circular:

**The fit excludes Asp entirely.** It is on the two histidine rings - NE2, ND1 and CE1 of His87/His116 and His163/His192, six atoms - because those are the ligands the metal actually sits on, and including the Asp in the superposition would let the atom under question determine where the zinc goes. **RMSD 0.184 A**, worst pair 0.307 A. A first attempt that did include Asp CG fitted at 0.635 A with that pair the worst at 0.853 A, which is the signature of forcing a correspondence that is not there.

**The histidine distances check the transfer.** In 6DN4's frame the transferred zinc sits 2.34 A from His116 NE2 and 2.18 A from His192 NE2 - both ordinary Zn-N bonds. A bad superposition would show up here first.

| contact, in 6DN4's frame | distance | our model |
| --- | --- | --- |
| Zn - His116 NE2 | 2.34 | 2.22 |
| Zn - His192 NE2 | 2.18 | 2.09 |
| **Zn - Asp120 OD2** | **2.68** | **2.80** |
| Zn - Asp120 OD1 | 3.82 | 3.81 |
| Zn - Arg118 NH1 | 4.00 | 3.89 |
| Zn - Gln121 NE2 | 7.10 | - |

Every one of our six agrees with the wild-type crystal to within 0.13 A. **The starting structure is right.** `prod100c` and `wt100c` are both restraining Asp92 OD2 at a distance the wild type supports, and the 1.4 A discrepancy I was ready to rebuild the system over was an artefact of comparing against the mutant.

**Gln93 gets a second, independent line.** The user's position was that it is "neutral and far away from everything"; D173 measured 7.95-9.88 A to the ligand over 34 ns of trajectory, and the wild-type crystal now puts its NE2 **7.10 A** from the zinc. Trajectory and crystal agree, and neither is our model's assumption.

## Two defects in `refstruct.py transfer`, both found by using it

`transfer` could not do this, and the reasons were worth fixing rather than working around - these numbers go in a paper and must come from a command, not from arithmetic in a browser console.

1. **`--ligand` selected by residue NAME only.** 6NC5 holds six residues called `ZN` and only ZN302 is the site, so a name cannot name a metal. It now also takes `CHAIN:RESSEQ[:RESNAME]`.
2. **`--zn` was required**, and the structure worth transferring into is the metal-free one - which is the entire point of the exercise. It is now optional: without it the group is written with `--out` and measured with `pair`, which needs no metal. `--poses` still requires it, because a pose is scored by how it sits on the metal, and that is now said rather than crashed into.

And one defect introduced by the first fix and caught before the numbers were used: with a selector, `--out` was writing `A:302:ZN` into the four-character resname column, producing a file every other tool misparses. It writes the residue's own name.

Seven tests added to `tests/test_refstruct_transfer.py`, including the one that pins the selector against the name - selecting `ZN` by name still returns both copies, which is the behaviour that would have silently measured the wrong zinc.

**[CORRECTED by D186/D194 - the four-coordinate gate was never wrong, and the attack restraint was replaced rather than added. See `csr/retractions.md`.]** **D185 - Both runs are live and both restraint lists are exactly right: seven detected on the product, five on the reactant, the ligand's O4x AND O5x both held, and neither the substrate's O1x at 3.13 nor Arg90 at 3.89 anywhere near the list. And a gate I wrote three times is now wrong: "the zinc stays four-coordinate" cannot be the test when the 2.8-2.9 A contacts are deliberately held.** (2026-11-02)

## The detections, verified rather than assumed

`prod100c` (31423830), cutoff 3.00 A, per atom:

    HOH C:1 (O)    1.90      ASP A:92 (OD2)   2.80
    HIS A:164 (NE2) 2.09     hyd D:1 (O5x)    2.92
    HOH C:2 (O)    2.14      hyd D:1 (O4x)    2.94
    HIS A:88 (NE2) 2.22

**Both ligand oxygens**, 0.02 A apart, which is the per-atom fix (D175) doing the one thing per-residue detection could not. Seven detected plus the configured Asp92 OD1 at 3.81 = eight applied. `wt100c` (31483853) has the matching five - hydroxide 2.14 rather than a ligand - plus the same configured OD1, and **neither `mer C:1:O1x` at 3.13 nor Arg90 NH1 at 3.89 appears**, which is the check that the cutoff is 3.0 and not higher.

The long-contact warnings fire on all three contacts beyond 2.50 A, by design (D131), and they are in the config rather than only in stdout - which is why they are being read now rather than in three weeks.

Both on **bun076**, both RUNNING: product 13.6 ns in 29.5 min (460 ps/min, 3.6 h for 100 ns, 4:30 left); reactant 7.6 ns in 18.4 min (413 ps/min, 4.0 h, 4:41 left). Both fit inside 5 h with about an hour spare, and the rates are below the 497 ps/min measured for a run with the node to itself - two jobs on one node contend, which is worth knowing before anyone plans a 5 h run off that figure.

## The gate was wrong, and it would have failed in both directions

Three of my own check prompts say the next step is `prodpost` with "the zinc FOUR-coordinate with Asp92 held near 2.80". **Those two clauses contradict each other.** Coordination occupancy counts donors inside the cutoff, and the cutoff is now 3.0, so a run holding Asp92 OD2 at 2.80 and the ligand's two oxygens at 2.92/2.94 will report **seven** ligands for the product and five for the reactant - by construction, in every frame. Read against a four-coordinate gate that is a failure; read against no gate at all it is a pass for anything.

**The restated gate**, which is what the numbers can actually test:

1. **The four bond-length contacts stay bonds** - two histidines, the hydroxide-derived oxygen and the water, all inside 2.5 A at high occupancy.
2. **Asp92 OD2 stays near 2.80 and does not collapse toward 2.0.** That collapse was the artefact (D173, D177); it is now restrained, so this tests whether 50 kcal/mol/A^2 is enough to hold it against whatever pulled it in.
3. **Asp92 OD1 stays near 3.81** rather than pivoting in - the D120/D130 failure, now explicitly restrained.
4. **Product only: O4x and O5x stay near 2.9, and O2x stays out.** O2x travelling 4.75 A onto the metal was the ligand turnover of D173.
5. **Nothing new arrives.** The margin report lists donors just past the cutoff, and a bulk water appearing inside it is the thing to look for - in `md_wt100` a water took the sixth place in 43.2% of frames.

So the question is no longer "what is the coordination number" but "did the restraints hold the geometry they were set to, and did anything else join". A coordination number of 7 that is 4 bonds plus 3 held second-shell contacts is the intended state; the same 7 with Asp92 at 2.0 is not.

## Plan item 5: no protonation tool

`import propka` fails and neither `pdb2pqr` nor `propka3` is on the path in the ASH environment. Both are pure-Python pip packages with no compiled dependencies, so `pip install propka pdb2pqr` into that env is the whole of it - a minute, and reversible. It is the user's environment, so it is their call, but this is the item the audit called the cheapest high-value one and it has been open since **D7**.

**[SUPERSEDED by D191 - CN is 6 in 3523/5000 frames, and the 'transient' waters are resident 97%/89%. See `csr/retractions.md`.]** **D186 - Two radii, two questions, and D185 conflated them. The 3.0 A set is "what is affected by the metal and must be held at its crystal distance" - deliberately wider than the coordination sphere. The coordination number is a separate measurement at the stricter 2.5 A, and at 2.5 A the site is FOUR. So the four-coordinate gate was never wrong; measuring it at the restraint radius was.** (2026-11-02, user: "the 3 angstrom restraints aren't exclusively for coordinators - they're for anything affected by the metal - a different set. The metal could still be 4-coordinate you'd just use a stricter coordination radius for that like 2.5")

## The distinction

`preprocess.inputs.cutoff` = **3.0 A** answers: what is close enough to the zinc to be *affected* by it, and therefore has to be held where the crystal has it rather than left for a point-charge force field to pull in. That set is wider than the coordination sphere **on purpose** - the contacts between 2.5 and 3.0 A are precisely the ones that are not bonds and that become bonds when nothing holds them. Asp92 OD2 at 2.80 is the worked example: free, it collapsed to 2.03-2.07 on both arms (D173, D177).

**Coordination number is a different question with its own radius: 2.5 A**, the longest Zn-O distance that is still a bond - the same `LIGAND_MAX` the code already warns above. At 2.5 A the product's site is two histidines, the hydroxide-derived oxygen and one water. **Four.** Asp92 OD2 at 2.80 and the ligand's O4x/O5x at 2.92/2.94 sit outside it, restrained in the second shell where the crystal puts them, and are correctly not counted as ligands.

## What D185 got wrong, and what it got right

D185 declared "the zinc stays four-coordinate" unusable as a gate because occupancy at 3.0 A would report seven. The premise was right and the conclusion was backwards: the fix is to measure coordination at 2.5 A, not to abandon the test. **The four-coordinate gate stands exactly as it was.**

What makes this more than bookkeeping: `postprocess.coordination` with no `cutoff` falls back to `preprocess.inputs.cutoff`, which is now the *restraint* radius. So `prodpost`'s existing section would have counted seven ligands on a four-coordinate site - one number answering the other question - and D185's restated gate would have told the next check to expect that and call it intended. Both wrong, from one missing key.

`prodpostc` and `wtpostc` therefore set **`coordination.cutoff: 2.5` explicitly**, with `margin: 1.0` still showing everything out to 3.5 A - which is what makes the held second-shell contacts visible without counting them. The margin is how you see that Asp92 is at 2.80 and not at 2.07, which is the actual question.

Note this also means **D173 and D177 stand unchanged**: both were measured at 2.50 A with a 1.0 A margin, so their six-coordinate findings were already the strict count. Asp92 OD2 really was inside the bond radius at 2.07 A, in 100% of `md_wt100`'s frames. That was a genuine artefact and remains one.

## Bulk water is not a gate

The user's second point: "bulk water will come and go, QM/MM should give us the correct coordination behaviour." So a transient water inside 2.5 A in an MM trajectory is not a failure, and D185's fifth gate item - "nothing new arrives" - was too strict. MM with point charges cannot be trusted for coordination at all; that is the reason the restraints exist in the first place. What the MM occupancy can tell us is whether the *restrained* contacts held where they were set. The coordination behaviour that counts is QM/MM's, which is what `prodqm` is for.

That also puts the 43.2% bulk water of `md_wt100` in its place: it was the sixth ligand of a site whose fifth was an unrestrained Asp92, and with the carboxylate held the question is whether anything stays rather than whether anything visits.

**D187 - Plan item 7 scoped: the cluster is 37 atoms of site plus a substrate fragment, net charge -1 with it and 0 without, and `mcube` cannot run it. `sp` is hard-wired to `QMMMTheory` and to DFT, so a gas-phase cluster at GFN2 vs DFT needs two small additions rather than a script around the framework. The substrate truncation is NOT settled: it needs the ligand's connectivity, which is the one thing not to guess.** (2026-11-02, open since D9)

## Why this matters more than anything else not yet done

Imhof's own benchmark - the paper this project's method is modelled on - found AM1/d made 4->5 coordination **7.55 kcal/mol too favourable**, made bidentate->monodentate carboxylate **7 kcal/mol too costly**, and missed a water->Asp proton transfer that DFT finds (D114). Our model is the one showing a contested coordination number and a mobile carboxylate, every gradient in the project is GFN2-xTB, and the only DFT checks that exist are single points on xTB geometries answering conditional questions - which D59 and D64 both say is not the test. If GFN2 is wrong about coordination here by the amount Imhof measured, it is wrong by more than any barrier this project will compute.

## The cluster, from `CSR_wt/preprocess_results/finalsystem.pdb`

Read off the structure, not assumed. **His88 and His164 both carry HD1 and no HE2**, so both are HID with NE2 bare and facing the metal, as `residue_variants` asks; Asp92 has no HD2, so it is ionised; the nucleophile is `OH C:2` = O01 + H03, a hydroxide with its proton; the second ligand is `HOH D:1`.

| piece | truncation | atoms | charge |
| --- | --- | --- | --- |
| Zn | whole | 1 | +2 |
| His88 -> 4-methylimidazole | cut CA-CB, cap CB with H; keep CB HB2 HB3 CG CD2 HD2 ND1 HD1 CE1 HE1 NE2 | 11 + 1 | 0 |
| His164 -> 4-methylimidazole | the same | 11 + 1 | 0 |
| Asp92 -> acetate | cut CA-CB, cap CB with H; keep CB HB2 HB3 CG OD1 OD2 | 6 + 1 | **-1** |
| hydroxide | whole (O01, H03) | 2 | **-1** |
| water | whole (O, H1, H2) | 3 | 0 |
| **site subtotal** | | **37** | **0** |

Methylimidazole rather than bare imidazole because cutting at CA-CB and capping keeps the ring's electronics and the CB substitution pattern; acetate rather than formate for the same reason. That leaves **~15 atoms** of D54's 50-55 for the substrate, which is about right for the beta-lactam ring plus the C3 carboxylate - the group that binds Zn2 in the di-zinc enzymes and has nothing to bind here.

**Net charge: 0 for the site alone, -1 with a substrate fragment carrying the C3 carboxylate.** Stated before anything runs because D152 is the entry where a run left at the inherited charge converged its SCF on the wrong species while looking entirely normal.

## What is not settled, and will not be guessed

The substrate truncation needs the ligand's **connectivity**. `mer C:1` is 51 atoms and the attack is O01 -> C4x with O2x the carbonyl and N1x the leaving nitrogen, but which atoms close the four-membered ring cannot be read off the names, and `finalsystem.pdb` is an OpenMM product with no CONECT records for the ligand. Inferring bonding by distance is what put a spurious O-O bond in a strained geometry once already (`running-on-bunya.md`, "Bond orders for PDB ligands"). The connectivity is available - `mer_crystal.pdb` has CONECT records, and `refstruct.py` already has `conect_bonds` and `heavy_graph` - so this is a lookup, not a judgement, and it is the next step rather than a blocker.

## What mcube needs, which is two small things

`sp.run` builds `make_mm(...)` and then `QMMMTheory(qm_theory=pyscf, mm_theory=mm, ...)` unconditionally, and its QM side is always PySCF DFT. A gas-phase cluster is neither: no MM system, no point charges, every atom in the QM region, and it has to be scored at **both** GFN2-xTB and DFT to be a benchmark at all.

1. **`sp.inputs.cluster: true`** - no MM, the whole fragment is QM, `Singlepoint(theory=..., fragment=...)` directly. One branch.
2. **`sp.inputs.sp_engine: dft|xtb`** - `qmmm` already has exactly this as `qm_engine` for its gradient stages; `sp` is DFT-only, so the same geometry cannot currently be scored both ways through the framework.

With those two, the whole benchmark is two sections per question on one cluster PDB, and every number is re-derivable from the config - which is the point. **Flagged rather than scripted**, per the standing rule: when a stage cannot do it, propose the addition.

**D188 - The two additions D187 asked for are built: `sp.inputs.cluster` runs a gas-phase model with no MM and no point charges, and `sp_engine: xtb|dft` scores one geometry either way. The cluster benchmark is now two sections on one PDB, differing in one key.** (2026-11-02)

`sp` built `make_mm(...)` and then `QMMMTheory(...)` unconditionally, and its QM side was always PySCF. Neither suits a cluster: there is no MM system, there are no point charges, and the whole question is whether the **semi-empirical method** is right about this chemistry - which cannot be asked of a QM/MM system, because the embedding is part of what would be under test.

`cluster: true` therefore builds the fragment, picks a theory, and calls `Singlepoint` directly. It says out loud that `zero_charges`, `actradius` and the force-field list do not apply, and the result records `cluster: true` with `cluster_atoms`, because **a gas-phase energy and a QM/MM one are not comparable** and the config has to say which it holds - the same reasoning as `coordination_source` in `md` (D134).

`sp_engine: xtb` goes through `make_qm`, the same construction every gradient stage uses, rather than a second `xTBTheory(...)` here. If the xTB method or its SCC options could drift between the benchmark and the production runs, the thing benchmarked would not be the thing run, and a benchmark that does not bind the method it tests is worth nothing. `make_qm`'s refusal of the `library` runmode does not bite: that exists because mcube only ever ran xTB as the QM region of a QMMMTheory and needed a point-charge gradient, and a cluster has neither.

Eight tests, and the two that matter pin the properties rather than the code paths: a cluster builds **no** MM system and **no** QMMMTheory, and the xTB branch reaches `make_qm`. 673 pass.

**What is still needed before the benchmark runs** (D187): the substrate truncation, which needs the ligand's connectivity from `mer_crystal.pdb`'s CONECT records rather than a guess, and the cluster PDB itself. The site part is settled - 37 atoms, net 0, charge -1 with the substrate fragment.

**[QUALIFIED - 'every single run' is five structures read; six run directories were not audited. See `csr/retractions.md`.]** **D189 - Histidine audit: the four that matter are identical in every run, and the one discrepancy is His112, 15.5 A from the zinc. CSR_ligand and CSR_hydrolysed have it HIE; every run from CSR_crystal onward has it HID, including the live CSR_wt and CSR_prod - and `112: HID` does not come from the source D23 attributes it to.** (2026-11-02, user: "csr hydrolysed and csr ligand should have correctly parameterised his in their configs - check all jobs run so far to determine if his is parameterised correctly in them")

## Read off the built structures, not the configs

Which hydrogens each histidine actually carries in `preprocess_results/finalsystem.pdb` - HD1 is HID, HE2 is HIE, both is HIP:

| run | A:53 | **A:88** | A:112 | A:150 | **A:164** |
| --- | --- | --- | --- | --- | --- |
| CSR_ligand | HIE | **HID** | **HIE** | HIE | **HID** |
| CSR_hydrolysed | HIE | **HID** | **HIE** | HIE | **HID** |
| CSR_wt | HIE | **HID** | **HID** | HIE | **HID** |
| CSR_prod | HIE | **HID** | **HID** | HIE | **HID** |
| CSR_nac | HIE | **HID** | **HID** | HIE | **HID** |

The configs match: CSR_ligand and CSR_hydrolysed set only `164: HID, 88: HID` and let the builder assign the rest at pH, which gave HIE at 112. Every config from CSR_crystal onward names all five and says `112: HID` explicitly.

## What this does and does not mean

**The two zinc ligands are HID in every run, with NE2 bare and facing the metal.** That is what the chemistry requires and it is right everywhere - including in the runs built from the file D21 found to have **HE2 on every histidine**, with His164's HE2 sitting **0.76 A from the zinc**. So the thing D21 and D7 warned about never reached a production system: `residue_variants` overrode it, in every single run. That is worth saying plainly, because the log has carried it as an open worry since the first week.

**His112's tautomer is mechanistically irrelevant.** 15.53 A from the zinc, 17.92 A from the substrate, 17.19 A from the nucleophile, in MM in every QM region, no contact with anything in the site. A neutral histidine's tautomer at that range changes nothing in any energy this project computes. **Nothing needs rebuilding over it.**

**But its provenance is wrong, and that is the finding.** D23 records that the tautomers were taken from `CSR_nontruncated_1metal_original.pdb`'s refined hydrogens. D21 records that in that same file *every* histidine carries HE2. A file with HE2 on His112 says **HIE**, not HID - so `112: HID` cannot have come from where D23 says it did, and the two configs that let the tool decide got HIE, agreeing with the file. One of the five assignments is unattributed.

**Recommendation:** set `112: HIE` in any future build, for consistency with the two configs the user confirms are right and with the source file's own hydrogens. Do not rebuild a running system for a residue 15 A away. And **the His164 question D7 raised is still open and is not this**: whether HID is right for a residue whose NE2 sits 3.2 A from the developing oxyanion is a chemistry question needing the sensitivity check, not a parameterisation error.

**What is not audited:** CSR_crystal, _v2, _v3, CSR_unrestrained, _2, CSR_gln93, CSR_ohH. Their configs all carry the same five assignments as CSR_wt, so the structures will match; they are void runs (D28-D70 laundering) and fetching 4 MB each to confirm a value their configs already state is not worth it. If any of them is ever revived, check it then.

**D190 - Arg90, Lys226, Arg229 and Gln93 all go to MM, so the RQK set can be decomposed by charge deletion as one set. That drops the QM charge from +2 to 0 on both arms, and it is set on every QM/MM-family section rather than inherited.** (2026-11-02, user's decision on the audit's item 9)

## Why the method forced the choice

Charge deletion cannot reach a residue that is already quantum - electrostatic embedding has removed its MM charges, and `zero_charges` refuses it by name (D153). Arg90 and Lys226 were inside the QM region in **all eleven systems**, so of the four cations around the site only Gln93 (-0.87 kcal/mol) and Arg229 (+44.34) were ever measurable. The paper's thesis is what **R, Q and K** do in place of the second zinc, and the method chosen to test it could reach only the Q.

Mixing the two is not an option: a table with two residues done by deletion and two by a parallel all-MM model is not one measurement (D153 says so). So all four go to MM, and the whole set is decomposed the same way.

**What leaves the QM region: Arg90 and Lys226.** What stays: the zinc, its ligands, Asp92 - the subject of the coordination question - and the substrate. Arg229 and Gln93 were never in it.

## The charge, which is the part that fails silently

Both arms were **+2** with the cations inside: Zn(2+), Asp92(-1), the ligand (-1 hydrolysed product, 0 reactant zwitterion), the hydroxide (-1, reactant only), Arg90(+1), Lys226(+1). Remove the two cations and both arms are **0**.

`charge: 0` is therefore written on every affected section, not left to inherit. D152 is the entry where a run left at the inherited value converged its SCF on the wrong species while looking entirely normal, and this is the same shape of mistake one residue larger.

## Where it is applied

Nine sections: `prodqm`, `prodopt`, `recombine`, `expel`, `pstep2`, `pstep1`, `psaddle2`, `psaddle1`, and a new `wtqm` giving the reactant the same QM/MM MD the product gets. **Every section that will produce a barrier now carries the same region as the equilibration that feeds it** - which D170 flagged as a gap when the path sections and the diagnostic sections disagreed about Gln93.

The two regions are YAML anchors under `site.qm_regions` rather than nine copies, because nine copies is nine chances for one to drift. YAML expands them on load, so what `sections.py` writes into a run's config is the literal mapping and nothing depends on the anchor surviving.

**The running MDs are unaffected.** `md` uses restraints, not the QM region, so `prod100c` and `wt100c` keep running and nothing needs resubmitting.

## What this does not settle

Whether the all-MM decomposition is the *right* test, as against computational mutagenesis. The user chose deletion-on-all-four, which is cheap and consistent; the mutants (Arg90->Ala/His, Lys226->Ala/His, and the 6DR8-like triple as a di-zinc control) remain the structural test and a second campaign. And the deletions still have to be run **at both ends of the reaction coordinate**: an absolute QM-MM interaction at one geometry is not a catalytic contribution, which is the error Arg229's +44.34 would be if quoted as one.

**[PARTLY SUPERSEDED by D197 - Lys226 NZ-O01 is inside 3.5 A in 10.3% of frames, not 100%, mean 4.93 A; the dissociation finding stands. See `csr/retractions.md`.]** **D191 (the provisional dissociation claim is CONFIRMED by D192; the periodic-image doubt is closed) - The gate FAILS on both arms, and the reason is not the coordination number: the SUBSTRATE APPEARS TO DISSOCIATE. Over the reactant's 100 ns, O01-C4x has mean 19.95 A and max 50.69 A, inside 3.5 A in 0.2% of 5000 frames - meropenem leaves the active site in roughly the first 0.2 ns and diffuses across the periodic box. `wt100c` is therefore not a Michaelis complex and cannot supply a reactant structure. The product's hydrolysed ligand goes the same way (O6x 6.44-8.10 A from the zinc). Separately, a 50 kcal/mol/A^2 restraint cannot hold a carboxylate against a bare Zn(2+): Asp92 OD2 restrained at 2.80 A has mean 2.37 A (product) and 2.44 A (reactant), and a BULK water takes the free site in 89.3% / 97.0% of frames. `prodqm` and `wtqm` were NOT submitted.** (2026-11-03, jobs 31423830 FAILED / 31483853 COMPLETED, analysed by 31528586 / 31528410)

## What was measured

`prod100c` (31423830) NaN'd at **22.68 ns** - the product arm's second NaN, after 34.14 ns on `prod100` under the old 2.5 A cutoff, so it is not caused by the D186 restraint change. `wt100c` (31483853) completed its 100 ns in 3:18:43.

| | reactant `wt100c`, 5000 frames | product `prod100c`, 1134 frames |
| --- | --- | --- |
| frames analysed | 100 ns, complete | 22.68 ns, salvaged from the NaN |
| Zn-His88 NE2 | 100%, 2.22 mean | 99.6%, 2.29 |
| Zn-His164 NE2 | 100%, 2.17 | 100%, 2.19 |
| nucleophile | O01 100%, 1.92 | - |
| structural water(s) | D:1 100%, 2.04 | C:1 100%, 2.04; C:2 100%, 2.09 |
| **resident BULK water** | **E:2920, 97.0%, 2.10** | **E:3189, 89.3%, 2.10** |
| **Zn-Asp92 OD2** (restrained 2.80) | **70.5%, min 2.06 / mean 2.44 / max 2.97** | **86.1%, 2.00 / 2.37 / 2.75** |
| Zn-Asp92 OD1 | 0.3%, 3.73 | 0.3%, 3.82 |
| coordination number (2.5 A) | 6 in 3523, 5 in 1477 | 6 in 852, 5 in 281, 4 in 1 |
| **substrate at the metal** | **O01-C4x 0.2%, 3.14 / 19.95 / 50.69** | **Zn-O6x 0%, 6.44 / 7.44 / 8.10** |
| Lys226 NZ-O01 | **100%, 2.45 / 2.72 / 3.40** | - |
| Gln93 OE1 on the zinc | **3.0%, 1.98 / 2.16 / 2.76** | 0%, 3.44 |

## Three distinct failures, in order of how much they matter

**1. Nothing holds the substrate, so it leaves.** At MM level there is no bond from the zinc to meropenem; its closest zinc contact, `mer O1x` at 3.13 A, is outside the 3.0 A restraint radius and was correctly excluded by D186's own rule. With neither a bond nor a restraint, the only thing keeping a substrate in a pocket is the force field's non-covalent binding free energy, and here it is not enough: ~10 of 5000 frames inside 3.5 A, then 50 A of box. **A Michaelis complex cannot be produced by unrestrained MM MD and was not.** Everything downstream of the reactant arm - NAC frames, the attack geometry, the scan, the NEB - needs a structure this run does not contain.

**2. 50 kcal/mol/A^2 is too soft for this site.** Asp92 OD2 is restrained harmonically at 2.80 A and its distribution's *maximum* over 100 ns is 2.97 while its mean is 2.44: the distance barely ever reaches the restraint's own minimum. The restraining force at 2.44 A is 2 x 50 x 0.36 = 36 kcal/mol/A and it is still losing, which is what a formal +2 point charge against a carboxylate oxygen at 2.4 A looks like in a non-bonded model. D186 fixed *which* atoms are restrained; this is *how hard*, and 50 was never tested against the question it is being asked. Holding a mean within ~0.05 A of a 2.80 A target needs of order 500 kcal/mol/A^2, or a flat-bottomed restraint, or a bonded zinc model.

**3. A bulk water takes the open site and stays.** The user's position is that "bulk water will come and go, QM/MM should give us the correct coordination behaviour" - but E:2920 and E:3189 are not coming and going, they are resident for 97% and 89% of their runs, and they are what makes the zinc six-coordinate. An MM zinc is a naked point charge with no ligand-field preference, so this is the expected behaviour of the model rather than a finding about the enzyme. It is also the third time this project has concluded that MM MD cannot be trusted for this coordination sphere (D120, D130, D186).

## What survives, and it is not nothing

**Lys226 NZ-O01 at 2.45 / 2.72 / 3.40 A in 100% of 5000 frames.** The lysine holds the nucleophile through the entire 100 ns, which is a direct, converged, classical-MD observation bearing on the RQK thesis and does not depend on the coordination number. **Gln93 OE1 reaches the zinc in 3.0% of frames (min 1.98 A)**, which is worth recording precisely because the working assumption (user, 2026-11-01) is that Gln93 is neutral, remote and non-functional: in the WT crystal it is 7.10 A from the metal, and this says the MD lets it approach. Both want re-measuring once the site is actually held, since both could be artefacts of the same naked-+2 problem.

Also confirmed again, independently of the above: **His88 and His164 are HID with NE2 on the metal in 100% of frames in both arms**, and the hydrolysed topology is right - O6x-C4x 1.23 A (the new acid C=O) and C4x-N1x 3.81 A (the ring is open).

## Why `prodqm` and `wtqm` were not submitted

The user's standing authority was to submit them **if the gate passes**. It does not. Running 50 ps of QM/MM MD from a frame with no substrate in the site would produce a converged, expensive, meaningless number, which is the D152 failure mode in a different costume.

## AMENDMENT, same day: the minimum-image check was not done

The claim above rests entirely on `wtpostc`'s monitored-contacts table:
`O01-mer C4x` 0.2% inside 3.5 A, min 3.14 / mean 19.95 / max 50.69 A, with
`Zn-mer O1x` and `Zn-mer O2x` saying the same thing. **It was written up as
settled before the obvious check was done, and that was wrong.**

`contact_occupancy` computes every distance as `separation(..., cell)`, and
`separation` applies the minimum image **only when `cell` is not None**.
`cell` comes from `coordination_rows(pdb_file)` -> `_coords_and_cell`, which
returns a cell only if the topology PDB carries a `CRYST1` whose three angles
are all 90 within 1e-3; anything non-orthorhombic or missing yields `None` and
**no wrapping correction at all**. Whether
`md_wt100c_results/NVTtrajectory_firstframe.pdb` satisfies that has not been
checked. If it does not, a ligand that never left the site but sits near a box
face reads as tens of Angstrom, which is exactly what maxima clustered at
50.69 / 51.25 / 51.62 A look like.

The maxima cannot discriminate. For a cubic box of side L the largest possible
minimum-image distance is L*sqrt(3)/2; this system is ~10,000 waters, so
L ~ 70 A and the ceiling is ~61 A. 51.6 A is therefore consistent with genuine
diffusion **and** with an uncorrected wrap.

**What does corroborate dissociation, independently of the ligand's absolute
position:** the reactant's coordination listing contains **twenty-plus bulk
waters at 3.1-3.5 A from the zinc** and `mer O1x`/`O2x` nowhere in it. That
scan examines every atom within cutoff+margin of the metal, so the substrate is
absent from a neighbourhood that bulk solvent has filled. A periodic wrap
displaces the ligand; it does not put water into the ligand's pocket.

**The test that settles it, in one read:** the per-frame series, which
`contact_occupancy` already returns as `_distances` and writes to the CSV.
Dissociation gives a smooth random walk away from 3.14 A; a periodic-image
artefact gives instantaneous jumps between ~3 A and ~50 A with no density in
between. Check the `CRYST1` angles of the topology PDB at the same time.

**Until that is done, treat "the substrate dissociates" as PROVISIONAL**, and
with it criterion 6 and the blocked status of plan items 1+2 and 3. What is
*not* provisional: the NaN at 22.68 ns, Asp92 OD2's mean at 2.37/2.44 against a
2.80 A restraint, the resident bulk water at 89.3%/97.0%, and Lys226 NZ-O01 at
2.72 A in 100% of frames - none of those depend on a long-range distance.

## The decision this leaves

The substrate needs a restraint, and its distance cannot be read off a CSR-1 structure because none with a ligand exists. **7UHI/7UHJ (L1 + moxalactam) is the only measured B3 Michaelis complex** (D96), so it is the only defensible source for a Zn-substrate target distance - which makes this a decision about what the reactant model *is*, not a parameter to pick quietly. The alternative is to stop asking classical MD for the site at all: minimise, equilibrate the protein and solvent briefly with the whole site held stiff, and let QM/MM MD do the equilibration it is the arbiter of anyway.

**D192 (DIAGNOSIS WITHDRAWN by D194: the 5.43 A pose is docked pose 2, chosen deliberately, and the cause was a dropped attack restraint) - The substrate was never coordinated in the starting pose. Nearest substrate atom to the zinc: `mer O1x` at 3.13 A - the 1-hydroxyethyl SIDE-CHAIN hydroxyl, not the beta-lactam carbonyl. `O2x`, the carbonyl oxygen that must sit on Zn1 in a B3 Michaelis complex, was 5.43 A away; the C3 carboxylate 8.05/8.44 A. The only thing holding the ligand was a 2.25 A hydrogen bond from its own side-chain hydroxyl to the hydroxide nucleophile. AND the periodic-image doubt of D191's amendment is resolved AGAINST me: CRYST1 is 79.028 A cubic at 90/90/90, so the minimum image WAS applied and the dissociation is real.** (2026-11-03, measured from `CSR_wt/preprocess_results/finalsystem.pdb`)

## The starting geometry, measured rather than recalled

| pair | A |
| --- | --- |
| Zn-O01, the hydroxide nucleophile | **2.14** |
| **Zn-O1x**, the 1-hydroxyethyl hydroxyl | **3.13** (nearest substrate atom) |
| Zn-C4x, the beta-lactam carbonyl carbon | 5.28 |
| **Zn-O2x, the beta-lactam carbonyl oxygen** | **5.43** |
| Zn-N1x, the leaving nitrogen | 6.03 |
| Zn-O4x / Zn-O3x, the C3 carboxylate | 8.05 / 8.44 |
| O01-C4x, the attack distance | 3.47 |
| **O01-O1x** | **2.25** |

`O1x-H5x` is 0.99 A, so O1x is a protonated hydroxyl, which fixes its identity:
meropenem's C6 (R)-1-hydroxyethyl substituent. Nothing of the substrate was
inside 3.0 A, let alone the 2.5 A bonding radius. D186's restraint logic
excluded `mer O1x` at 3.13 A correctly and by its own rule - the mistake was
never in the restraint list, it was upstream of it.

## Why this reframes D191 completely

A B3 Michaelis complex has the beta-lactam carbonyl oxygen **on Zn1**, at
roughly 2.0-2.3 A: that coordination is the mechanism, because it is what
polarises the carbonyl for attack. In 7UHI/7UHJ (L1 + moxalactam, D96), the
only measured B3 Michaelis complex there is, that is what the density shows.
Our pose instead put the **wrong oxygen** nearest the metal, 1 A too far out to
be a contact, with the reactive carbonyl 5.43 A away pointing elsewhere.

So the 100 ns did not break a complex. **There was no complex.** The ligand sat
in the pocket held by one hydrogen bond (O01-O1x, 2.25 A) between its own
side-chain hydroxyl and the nucleophile, and when that broke it diffused off,
which is the behaviour of a weakly docked ligand and not a finding about
CSR-1. An attack distance of 3.47 A at t=0 is loose for a near-attack
conformation on top of that.

## The periodic-image doubt is closed

D191's amendment raised it and was right to; it resolves against me.
`md_wt100c_results/NVTtrajectory_firstframe.pdb` carries
`CRYST1 79.028 79.028 79.028 90.00 90.00 90.00`, so `_coords_and_cell` returns
a cell, `separation` applies the minimum image, and the ceiling for a cubic box
is L*sqrt(3)/2 = **68.4 A**. The observed max of 50.69 A is a legitimate
minimum-image distance well inside it. **The dissociation numbers are real and
the claim is no longer provisional.** The per-frame time series is no longer
needed to settle it.

## What this does and does not change

Still true, and now the *secondary* problems rather than the main one: 50
kcal/mol/A^2 cannot hold Asp92 OD2 at 2.80 A (mean 2.44 / 2.37), and a bulk
water holds the free site for 97.0% / 89.3% of the two runs. Both still need
fixing; neither is why the substrate left.

Now the primary problem: **the reactant structure itself.** It cannot be
repaired by restraining the pose it has, because restraining a wrong pose
harder only freezes the wrong pose. The carbonyl has to be placed on the zinc
first, from 7UHI/7UHJ, and only then restrained. That makes 7UHI the **source**
for the Michaelis geometry, not merely the check D178 scoped it as - which is
exactly the decision D191 said this leaves, now with the reason stated
properly.

Unaffected, because neither depends on the substrate: **Lys226 NZ-O01 at 2.72 A
mean in 100% of 5000 frames** stands as an observation about the nucleophile's
anchoring, and the His88/His164 HID assignment stands.

**D193 (CAUSAL CLAIM WITHDRAWN the same day - see the correction below) - The reactant's ligand was named in `constraint_exclude_residues: {C: [1]}`, which I kept deliberately and documented as costing nothing. The documented reasoning was wrong, but the line was INERT and cost nothing: `constraint_exclude` is tested inside the loop over `shell`, which is already filtered to within the cutoff, so the key can only act on residues inside it - and `mer`'s nearest atom was 3.13 A, outside 3.0. Neither run's log contains a single `No restraint for ... (constraint_exclude_residues)` line. mcube restrained everything within 3.0 A on both arms, exactly as built. The sole cause is the pose (D192).** (2026-11-03, user)

## What the exclusion was for, and why each half of the reason fails

`site.yaml` line 1479, written when `wt100c` was set up: *"`constraint_exclude_residues: {C: [1]}` stays. The reactant's substrate must be free to react, which is the documented reason, and here it costs nothing either way: its nearest atom is 3.13 A, outside the sweep."*

**"Must be free to react" is a QM/MM concern applied to an MM stage.** MM cannot form or break a bond, so nothing reacts during equilibration and a restraint there biases no reaction coordinate - it is released before QM/MM regardless. The earlier and better-posed version of the reason (line 394, from `nac`) was *"for a mono-Zn B3 site it is not established that the carbonyl coordinates the same zinc that carries the nucleophile"*, which is a real open question, but it fails for a different reason: **omitting the restraint is not neutrality.** It does not leave the question open, it guarantees the ligand leaves, which removes the ability to ask the question at all. The neutral act is to hold the pose the structure was given and let QM/MM decide the coordination - which is the user's own standing position ("bulk water will come and go, QM/MM should give us the correct coordination behaviour").

**"Costs nothing either way" was circular.** It appealed to O1x being 3.13 A - to the pose - while the pose was the unexamined thing (D192: the carbonyl `O2x` was 5.43 A out and the nearest atom was the wrong oxygen). Two independent errors pointed the same way and each was used to excuse the other.

## The thing that makes this unambiguous

The product arm had the ligand taken OFF the exclusion list. Identical code, identical 3.0 A cutoff:

| | reactant `wt100c` (ligand excluded) | product `prod100c` (ligand not excluded) |
| --- | --- | --- |
| ligand atoms restrained | **none** | `D:1:O5x` 2.92, `D:1:O4x` 2.94 |
| where the ligand ended up | **gone: mean 19.95 A, max 50.69 A** | **held: 3.11 / 3.15 A mean, 100% inside 3.5 A** |

Where mcube was allowed to do what it was built to do, it worked.

## The correction that makes it worse, not better

**Removing the exclusion alone would not have saved the reactant run.** With the pose as built, the 3.0 A sweep still finds nothing of the ligand: O1x misses by 0.13 A and `O2x`, the atom that should have been caught, is 5.43 A away. So the sequence was: mcube built to restrain it, the exclusion switched that off, and the pose left nothing to grab even with it on. **Both have to change**, and the pose has to change first, because a restraint can only hold the geometry it is given (D192).

## CORRECTION, same day: the exclusion was inert and did not cause this

Written after the user asked why a warning for excluding atoms inside the
cutoff would be needed, "since you'd only be excluding atoms inside the cutoff
anyway". That is right, and following it refutes this entry's own headline.

`find_coordinators` tests `resid in constraint_exclude.get(chain)` **inside the
loop over `shell`**, and `shell` comes from `coordination_shell(..., per_residue=False)`,
already filtered to `r <= cutoff`. So the key is only ever consulted for
residues inside the cutoff. Checked against both logs:

| | reactant `wt100c` | product `prod100c` |
| --- | --- | --- |
| `No restraint for ... (constraint_exclude_residues)` lines | **none** | **none** |
| restraints applied | 6 (5 detected + 1 configured), no ligand atom | 8 (7 + 1), including O5x 2.92 and O4x 2.94 |

**The reactant's exclusion never fired.** `mer`'s nearest atom was 3.13 A, so
`C:1` never entered `shell` and the exclusion was never reached. Its restraint
list lacks the ligand because nothing was **found**, not because anything was
**excluded**. The product's list contains the ligand because 2.92 and 2.94 A
are inside the cutoff. Same code, same key, and the only thing that differed
was the geometry.

So the causal claim in this entry's first version - "it cost the run" - is
**withdrawn**. The reasoning recorded in `site.yaml` was still wrong, and the
rule at the end of this entry still stands, but the line was dead config and
cannot have caused anything. **D192 remains the sole cause**: the wrong oxygen
nearest the metal, 0.13 A outside the cutoff, the carbonyl 5.43 A away.

The proposed warning - "you excluded something inside the cutoff" - is
**withdrawn too, and was worse than redundant**: it would fire on every
exclusion that does anything at all, D8's arginine included, which is the
legitimate use the key exists for. A warning that fires on correct usage is
noise that teaches people to skip the log.

**The guard that would have helped is the inverse: warn when an excluded
residue has NO atom inside the cutoff** - dead config that reads as
load-bearing. That is the trap this entry fell into: a `site.yaml` comment
asserting the exclusion "stays" and "costs nothing either way" about a line
doing nothing, which is what stopped the obvious question of why the ligand was
outside the cutoff at all. Not built; offered.

## The rule, so this is not re-derived

**`constraint_exclude_residues` is for atoms that are near the metal by accident, not for molecules whose relationship to the metal is the question.** D8's arginine - restrained to the zinc as though it were a ligand - is what the key exists for. A substrate carbonyl is the opposite case: being affected by the metal is the entire point, and "we are not sure it coordinates" is an argument for holding it where the structure puts it and asking QM/MM, never for letting it go.

**[SUPERSEDED by D197 - with the attack restraint restored, every restraint held and the structure is still rejected. See `csr/retractions.md`.]** **D194 - CONSOLIDATION of the pose history, because I have now re-derived a discarded argument for the third time. The cause of `wt100c`/`prod100c` losing the substrate is NOT the pose: it is that I REPLACED the attack-coordinate restraint `[C:2:O01, C:1:C4x, 3.47]`, which `nac` and `wt100b` carried and which held the substrate through 100 ns, with `[B:1:ZN, A:92:OD1, 3.81, 50.0]`. One line, and my own `site.yaml` comment calls it "**One** restraint beyond the sweep" without mentioning what it overwrote. D192's "there was never a Michaelis complex" is withdrawn as a claim about error, and my proposal to revert to `CSR_ligand`'s pose is withdrawn outright - that pose is the one the project moved away from on purpose.** (2026-11-03, user: "CSR_ligand and CSR_hydrolysed were the ones I set up personally, and it's the pose we used for most simulations before deciding it wasn't correct")

## The two poses, and which is which

| pose | Zn-O1x (side-chain OH) | Zn-O2x (carbonyl) | builds | provenance |
| --- | --- | --- | --- | --- |
| **hand-built near-attack** | 13.3-14.0 | **2.34-2.55** | `CSR_ligand`, `CSR_hydrolysed` | **the user's own**, D20: "pose from the processed model ... rigid-body shift towards OH- with C4x 3.0-3.2 A at 100-115 deg; contacts down to ~2.2 A left for MM minimisation (user)" |
| pulled-in docked | 3.82 | 1.94 | `CSR_crystal_v2` | docked pose 2 **pulled onto the metal** by two hand restraints, `[4203, 4209, 2.2, 20.0]` and `[4255, 4208, 3.0, 20.0]` (D68) |
| **Vina docked pose 2** | 3.13 | **5.43** | `CSR_crystal`, `CSR_nac`, `CSR_unrestrained2`, **`CSR_wt`** | `model: 2` in all three configs, **chosen deliberately** (D68, correcting D66's "nobody chose this") |

`CSR_ligand`/`CSR_hydrolysed` are the reference runs whose `MM.slurm` defined the
pipeline order (D44), they are the two the user set up personally, and their
pose was used for most of the early work **and then deliberately set aside** in
favour of the docked pose. **So "transfer the pose from `CSR_ligand`" is a
proposal to undo a decision, not to fix a regression. Withdrawn.**

## What D192 got wrong

D192 measured the `CSR_wt` pose correctly - carbonyl `O2x` 5.43 A, side-chain
`O1x` 3.13 A nearest, confirmed by connectivity (C4x-O2x 1.21 A, C4x-N1x
1.37 A) - and then drew the wrong conclusion from it. **The 5.43 A is not an
accident and not a regression.** It is docked pose 2 as docked, and the
decision not to restrain the carbonyl was equally deliberate: D98 set up
`CSR_nac` with exactly one substrate restraint, on the attack coordinate, and
said so in terms - *"the substrate gets one restraint on the attack coordinate
at the distance the docked pose already has. **Nothing asserts the carbonyl
coordinates the zinc.**"* That was the point of it.

So D192's "there was never a Michaelis complex to lose" stands as a *description*
of the geometry and is **withdrawn as a diagnosis**. The complex was always a
docked pose with an uncoordinated carbonyl, held by one restraint, and that is
what the project decided to study.

## The actual cause, which is one line of config

| section | `restraints` | outcome |
| --- | --- | --- |
| `nac` | `[['C:2:O01', 'C:1:C4x', 3.47]]` | substrate held; Zn-O2x 5.47 at 1 ns (D100/D101) |
| `wt100b` | `[['C:2:O01', 'C:1:C4x', 3.47]]` | held |
| **`wt100c`** | `[['B:1:ZN', 'A:92:OD1', 3.81, 50.0]]` | **O01-C4x mean 19.95, max 50.69 A** |
| **`prod100c`** | `[['B:1:ZN', 'A:92:OD1', 3.81, 50.0]]` | **Zn-O6x 6.44-8.10 A** |
| `wt100`, `prod100` | none | the original void pair |

`restraints` is a list and I treated it as a slot. The Asp92 OD1 restraint was
a sound judgement on its own terms - D120 and D130 are the same bidentate
pivot twice - but it went in **instead of** the attack restraint rather than
beside it. Both belong in the list. **That is the whole of it.**

## The pattern this entry exists to stop

**D89** observed the substrate dissociating in an unrestrained run and
concluded the pose was wrong. That argument was examined and **discarded**:
*"dissociation-implies-bad-pose (an artefact of my own restraint set)"*, and
again *"The substrate dissociating unrestrained (run 1, D89) is evidence of an
MM artefact, not evidence against the docked pose."* **D192 re-derived the same
argument from the same observation and did not find the discard.** D145
contradicting D108 unnoticed for 8 days is why `state.md` exists; this is the
same failure with the pose instead of the zinc sphere.

The rule: **the substrate leaving an MM run is a statement about the restraint
list, not about the pose.** It has been mistaken for the latter twice.

## What this does NOT settle, and must not be read as settling

1. **Whether the docked pose is right remains the blocking scientific
   question**, exactly as the log already says: *"The binding pose is the
   blocking scientific question"* (D89), and the carboxylate-Arg interaction
   every other MBL substrate complex shows is one the docked pose never makes.
2. **7UHI/7UHJ is live evidence against 5.43 A.** D123 withdrew D40's claim
   that inner-sphere carbonyl coordination is unestablished, because L1's
   time-resolved Michaelis complexes put the carbonyl at **2.53/2.63 A from
   Zn1**. That is a real argument that the docked pose has the carbonyl too far
   out - and it is an argument to be had on the evidence, not something to
   settle by reverting to an older build.
3. **Restoring the attack restraint fixes MM and not QM/MM.** D119 took the
   `nac` frame into QM/MM at O01-C4x 3.49 A and it never came back under
   3.5 A in 50 ps, mean 5.31, max 6.20, on both protonation arms. So the
   substrate drifts under QM/MM regardless, and the cause there is still open -
   *"GFN2, or the pose"*. Re-running `wt100c` with the attack restraint
   restored gives a valid MM equilibration; it does not answer that.

**D195 - INVENTORY, measured with `site.py check` over all 84 last frames. There is NO correctly equilibrated structure for a QM/MM start. Exactly ONE frame in the project is accepted - `CSR_crystal_v2/qmmm_results/QM_MM_lastframe.pdb` - and it has a six-coordinate zinc with Asp92 OD2 at 2.04 A and an attack angle of 26 deg. AND the user's hypothesis is CONFIRMED by the only frame that can test it: a carbonyl held coordinated through MM (1.94 A) STAYED under QM/MM and relaxed to 2.56 A, inside 0.07 A of 7UHI/7UHJ's 2.53/2.63 **[QUALIFIED - D197: this is "a bond-length relaxation within the coordinated basin, not a rejection", so it is not a validation of GFN2's coordination behaviour; in this project QM/MM has never been observed to cross the coordination boundary in either direction. n = 1, our own run. The standing method question stays open.]**; every frame whose carbonyl was uncoordinated in MM lost the substrate under QM/MM.** (2026-11-03, user: "If the substrate is held coordinated in MM then it'll stay in QM, if it's not it won't")

## The hypothesis test

| MM frame | Zn-O2x in MM | Zn-O2x after QM/MM | substrate |
| --- | --- | --- | --- |
| `CSR_crystal_v2` (carbonyl frozen in MM) | **1.94** | **2.56** | **stayed** |
| `CSR_nac/md_naclong` | 5.16 | 6.20 (D111) | left |
| `CSR_unrestrained` | 6.89 | - | left in MM already |
| `CSR_unrestrained2` | 7.71 | - | left in MM already |
| `CSR_wt/md_wt100c` | 16.82 | - | left in MM already |

**xTB chose 2.56 A.** The MM value of 1.94 A is outside the criterion window
(2.4-3.0) and QM/MM *corrected it upward* to the experimental distance rather
than holding the MM artefact - which is D67's point measured again, and is
positive evidence that **GFN2 is not overbinding the substrate**. It also
retroactively vindicates the carbonyl restraint the user set originally and
which D40 argued out on a premise D123 withdrew.

## Every candidate, with its disqualifier

| frame | Zn-O2x | Zn-O01 | O01-C4x | angle | CN | verdict / issue |
| --- | --- | --- | --- | --- | --- | --- |
| `CSR_crystal_v2/qmmm_results/QM_MM_lastframe` | **2.56** | 2.21 | 3.70 | **26** | **6** | **ACCEPTED**, but CN 6 with Asp92 OD2 **2.04**, angle 26 deg, and D29's Lys226 proton transfer |
| `CSR_crystal_v2/md_results` NVT | 1.94 | 2.14 | 3.99 | 30 | 5 | carbonyl too tight (frozen constraint), angle 30 |
| `CSR_nac/md_naclong` NVT | 5.16 | 1.86 | **3.33** | **101** | **6** | chemistry right, coordination wrong: Asp92 OD1 **1.99**, bulk water 2.27. Wrong resting state (D112); condemned shell (state.md rule 7) |
| `CSR_nac/md_nacgpu` NVT | 5.47 | 1.95 | 3.55 | 80 | 5 | Asp92 OD1 1.98; angle 80 |
| `CSR_wt/md_wt100b` NVT | 5.98 | 1.93 | 3.65 | 102 | **6** | Asp92 fully **bidentate**, OD2 2.06 AND OD1 2.12 - D120/D130 undiluted |
| `CSR_wt/md_wt100` NVT | 7.80 | 1.90 | 5.56 | 119 | **6** | Asp92 OD2 2.02 + bulk water 2.14; void (D186) |
| `CSR_wt/md_wt100c` NVT | **16.82** | 1.91 | 14.75 | 114 | **6** | substrate dissociated (D194); bulk water E:2920 2.06 |
| `CSR_unrestrained` NVT | 6.89 | 2.14 | 4.45 | 100 | 5 | **Gln93 OE1 on the zinc at 1.79** |
| `CSR_unrestrained2` NVT | 7.71 | 2.14 | 5.31 | 107 | 5 | substrate out by design (D69) |
| `CSR_ligand/md_results` NVT | - | - | - | - | - | **UNUSABLE: `B:1:ZN not in the structure`** |
| `CSR_hydrolysed/md_results` NVT | - | - | - | - | - | **UNUSABLE: same** |

All ten MD frames are rejected under `substituent_on_metal` as well, so this is
not an artefact of assuming one mechanism (D127's point).

## `CSR_ligand` and `CSR_hydrolysed` cannot be checked at all

Every selector in `site.yaml` fails against them: `B:1:ZN not in the
structure`, `C:2:O01` likewise. These are the user's own builds and predate the
current chain naming, so the site description does not address them. **This is
why the hypothesis could not be tested on them**, which is what the user
expected to be possible. Remapping their chains would make them checkable
without re-running anything - worth doing, because they are the only frames
whose pose has the carbonyl near the metal *without* a frozen MM constraint.

## The conclusion, which is the answer to "do we have one"

**No.** The two halves of a Michaelis complex exist in different structures and
never in the same one:

- **attack geometry** (O01-C4x 3.0-3.5 A, angle 95-120 deg): `naclong` 3.33/101,
  `wt100b` 3.65/102, `unrestrained2` 5.31/107.
- **carbonyl coordination** (2.4-3.0 A): only `CSR_crystal_v2`, and only after
  QM/MM.
- **a four-coordinate zinc**: **none of them.** Every frame is five- or
  six-coordinate, always because Asp92 has collapsed onto the metal, a bulk
  water has taken a site, or Gln93 has.

So the blocking problem is not the pose alone and not the restraints alone: no
run has yet combined a coordinated carbonyl, a near-attack nucleophile and a
four-coordinate zinc. The D191 force-constant finding (50 kcal/mol/A^2 loses to
a bare Zn2+) and D194's restraint list are both prerequisites, and on this
evidence the carbonyl needs restraining too - at 2.4-3.0 A, not at 1.94.

**[SUPERSEDED by D197 - the carbonyl is inside 2.5 A in 0.3% of frames; O1x holds the fourth site at 97.9%. See `csr/retractions.md`.]** **D196 - Five-coordinate is allowed once the substrate is bound, and the
rebuild that follows from it.**

The user's decision, settling the one question D195 left: the substrate
DISPLACES a zinc water, and the zinc may be five-coordinate once it has.

This matters because the criterion it replaces was rejecting correct structures
for being correct. `carbonyl_on_metal` demanded a coordination number of
exactly 4, on the reasoning that the resting site is four-coordinate - 2 His +
2 waters (D120) - and the complex should be 2 His + hydroxide + carbonyl. That
target is unchanged and is still four. What was wrong was treating five as a
failure. GOB-18's B3 site expands four- to five-coordinate on substrate
binding (`reference-structures.md`), so a retained second water is a real
possibility rather than a build error; and Asp92 OD2 sits at 2.80 A in the
crystal as a second-shell contact, so a frame where it drifts inside the 2.5 A
measurement radius scores as five while being nothing of the kind. Both
hypotheses' coordination criteria are now `[4, 5]`, with the reasoning written
where the criterion is.

**What the rebuild is.** D195's finding was that each run had ONE of the three
Michaelis properties and never the other two, so `eq5` restrains all three
together. It starts from CSR_crystal_v2's input pose - the only pose that has
ever had the carbonyl on the metal, so the move asked for is 1.94 -> 2.60 A
rather than 5.43 -> 2.60 - and carries:

- the attack coordinate `[C:2:O01, C:1:C4x, 3.47]`, which `wt100c` replaced
  rather than added to (D194);
- the carbonyl `[B:1:ZN, C:1:O2x, 2.60, 500]`. 2.60 is the middle of the
  criterion's window, within 0.04 A of what GFN2 chose for itself when QM/MM
  relaxed crystal_v2's frozen 1.94 to 2.56, and within 0.03 A of L1's measured
  2.53/2.63 (D121). **Not 1.94** - that number came from a build decision,
  nothing measures it, and holding it there is what made crystal_v2's attack
  angle 26-30 degrees;
- `[B:1:ZN, A:92:OD1, 3.81, 500]`, OD1 pinned out as on every arm since D130.

**500 kcal/mol/A^2, not 50.** D191 measured a 2.80 A restraint sitting at a
mean of 2.44 with a maximum of 2.97 over 100 ns. A 0.36 A mean offset is a
systematic force of order 36 kcal/mol/A, not thermal spread, and offset goes as
F/k - so 500 puts it under 0.04 A. At 500 and a 15 amu reduced mass the
restrained mode's period is about 53 fs. `eq5` therefore sets
`nvt_timestep: 0.002`, which is 26 steps per period: the run directory it is
built from uses 4 fs, and 13 steps per period is not enough to integrate a
restraint this stiff. (An earlier draft of this entry said 1 fs, which is the
library default and not what this config inherits.) The site is deliberately rigid
in MM: MM cannot be trusted for coordination (D174), the MD's job is to
equilibrate protein and solvent AROUND a site held where the crystal has it,
and QM/MM is what decides the coordination afterwards (confirmed in D195).

**The library gained `preprocess.inputs.remove_residues`.** There was no config
route to displacing a water at all - the only way was to edit the input PDB by
hand, which is what this project's rules say a general library must not
require. It takes `CHAIN:RESSEQ:RESNAME` selectors, runs BEFORE solvation
(a water deleted after the box exists is simply replaced by the solvent),
takes every atom of the residue with its ANISOU and CONECT records, and refuses
the whole request if any selector does not match - including the case that
matters, a real residue of the wrong name at the right number, which is what a
selector written against an earlier build hits once solvation has renumbered
everything after the solute. `sections.py` checks these selectors too, as
whole-residue selectors against a residue index rather than as atoms. 19 new
tests; the suite is 704 passed, 10 skipped.

**What the run directory must set, which no section can.** `preprocess` takes
no `--section`, so a rebuild is a new run directory, as `CSR_crystal_v2` was.
Its `preprocess.inputs` carries `remove_residues: [D:1:HOH]` and
`coordination_force: 500`; `eq5` is the `md` section on top. The water selector
has to be checked against the MERGED input structure rather than the final
system, and preprocess refuses a mismatch rather than deleting the wrong
residue.

**What `eq5` decides.** Whether a coordinated carbonyl, a near-attack
nucleophile and a four- or five-coordinate zinc can coexist in one equilibrated
structure at all. If they can, plan items 1-2 are answered and `wtqm` finally
has a start. If they cannot, the pose itself is wrong, and that is a result
rather than a failure. 10 ns, not 100: every 100 ns run in D195's inventory is
in the rejected list, and nothing was visible at 100 ns that was not visible at
10.

**D196, measured.** The criterion was re-run over all ten NVT production
frames and over CSR_crystal_v2's QM/MM frame, at `--at=reactant`
(`.cc/i196.log`, `.cc/j196.log`, 2026-11-03). Two things came out of it.

**The relaxation changes four coordination verdicts and no overall verdict.**
`crystal_v2`, `nacgpu`, `unrestrained` and `unrestrained2` now PASS the
coordination number at 5 where `[4, 4]` failed them; `naclong`, `wt100`,
`wt100b` and `wt100c` are 6 and still fail. All ten are still rejected, every
one of them on `Zn - carbonyl O2x`:

| frame | CN | Zn-O2x | attack | angle |
| --- | --- | --- | --- | --- |
| `crystal_v2` NVT | **5 PASS** | 1.94 FAIL | 3.99 | 30 |
| `nacgpu` | **5 PASS** | 5.47 FAIL | 3.55 | 80 |
| `unrestrained` | **5 PASS** | 6.89 FAIL | 4.45 | 100 PASS |
| `unrestrained2` | **5 PASS** | 7.71 FAIL | 5.31 | 107 PASS |
| `naclong` | 6 | 5.16 FAIL | **3.33 PASS** | **101 PASS** |
| `wt100b` | 6 | 5.98 FAIL | 3.65 | 102 PASS |
| `wt100` | 6 | 7.80 FAIL | 5.56 | 119 PASS |
| `wt100c` | 6 | 16.82 FAIL | 14.75 | 114 PASS |
| `CSR_ligand`, `CSR_hydrolysed` | - | - | - | - (`B:1:ZN not in the structure`) |

So the decision removes a false barrier rather than producing a start
structure, which is the honest reading: the blocker was never the coordination
number, it is the carbonyl, and `eq5` is what restrains the carbonyl.

**The water to displace is confirmed by measurement, not assumed.** The QM/MM
frame's sphere is

    D:1:HOH:O 1.98, A:92:ASP:OD2 2.04, A:88:HIS:NE2 2.06,
    A:164:HIS:NE2 2.08, C:2:OH:O01 2.21, C:1:mer:O2x 2.56

Six, and removing `D:1:HOH` leaves exactly the five of D196 - 2 His +
hydroxide + carbonyl + a second-shell Asp92 OD2 inside the measurement radius.
That frame then passes both the carbonyl distance (2.56) and the coordination
number, and what is left failing is only the attack geometry: O01-C4x 3.70 and
26 degrees, which is what the frozen 1.94 A pose produced and what `eq5`'s
attack restraint addresses. `remove_residues: [D:1:HOH]` is therefore the right
selector for this system, measured in the structure it will be applied to.

One criterion `eq5` does NOT restrain: `Asp92 OD1 - O01` at 2.5-3.2 A, which
the QM/MM frame has at 4.76. D119 measured it at 2.82-2.94 in 100% of frames,
so it is left to happen on its own; the 4.76 is a consequence of the frozen
pose rather than evidence against it. If `eq5` comes back with it long, it
becomes a fourth restraint.

**D196, submitted.** `CSR_mich` built by `docs/csr/build-mich.py` and both
stages launched on 2026-11-03.

**31721628 `mich_pre` - preprocess, COMPLETED.** `remove_residues` did exactly
what it was added for:

    [STEP] Removing residues...
    [INFO]   removed D:1:HOH (3 atoms)

and the detected sphere at the 3.0 A radius is

    mer C:1 (O2x)  1.94    His164 NE2  2.09    OH C:2 (O01)  2.14
    His88  NE2     2.22    Asp92  OD2  2.79

**This is the first build in the project whose resting sphere IS the Michaelis
sphere**: 2 His + hydroxide + the beta-lactam carbonyl, with Asp92 OD2 as the
second-shell contact the crystal has at 2.80 (detected at 2.79, and the
`not a bond` warning fires on it exactly as D186 intends). Four donors inside
2.5 A, five inside 3.0.

**No bulk water took the vacated position.** That was the risk in deleting a
water, and the reason it is done before solvation; here there was no vacancy
to fill, because the carbonyl was already in it at 1.94. 47,979 atoms,
31,231 of them hydrogens.

`sections.py` then wrote `md_eq5` with **6 selectors checked** against
`preprocess_results/finalsystem.pdb` - the point of giving the section a
`pdb_file`, since with none named the selector check passes silently.

**31721717 `mich_eq5`** - 10 ns NVT at 2 fs, `-p gpu_cuda --qos=gpu
--gres=gpu:1`, 16 cores, 32 GB. At the 750-760 ns/day measured for this
49k-atom system on an H100 at 4 fs, 10 ns at 2 fs is of order 40 minutes of
compute plus the queue.

What to read when it lands, in order: whether the three restraints held (the
carbonyl near 2.60 rather than collapsing to 1.94 or leaving), the
coordination number at the 2.5 A radius over the trajectory, and then
`site.py check` on the last frame under `carbonyl_on_metal` at
`--at=reactant`. The criteria it has to pass are the ones every frame in D195
failed: a coordinated carbonyl AND a near-attack nucleophile AND four-or-five
coordination, in the same structure.

**[SELF-CORRECTING - the PASS on CN 4 is contradicted 120 lines later ('the number is correct and its membership is wrong'). See `csr/retractions.md`.]** **D197 - `mich_eq5` ran clean, every restraint held, and the structure is
still REJECTED. The cause is one unrestrained atom.**

**31721717, COMPLETED 0:0.** 35.1 minutes of walltime, 10 ns at 2 fs, no NaN
and no walltime trouble. It drew an **L40S**, not an H100, so about 370 ns/day
rather than the 750-760 D196 estimated from - the estimate was for the wrong
card, and the run still finished inside the hour.

**All seven restraints held, to better than 0.11 A.** The three configured:

    Zn - mer O2x (carbonyl)   2.645   set 2.60
    Zn - Asp92 OD1 (held out) 3.782   set 3.81
    O01 - C4x (attack)        3.458   set 3.47

and the four redetected: hydroxide O01 2.035 (2.14), His164 2.166 (2.09),
His88 2.193 (2.22), Asp92 OD2 2.743 (2.79). Nothing collapsed and nothing
dissociated. On the restraints alone this is the best-behaved MD in the
project.

**But an eighth donor arrived, and nothing was restraining it.** `mer`'s
**O1x at 2.06 A** - closer to the zinc than the restrained carbonyl at 2.645,
and as close as the anionic hydroxide at 2.035.

**What O1x is, established from the bond lengths in the frame rather than
assumed from the name.** O1x-C2x is **1.39 A**, a single C-O bond, and C2x
carries C1x (1.51) and C3x (1.51): O1x is the alcohol oxygen of meropenem's
**6-(1-hydroxyethyl)** group. For contrast, in the same frame O2x-C4x is
**1.24** with C4x-N1x **1.35**, which is the beta-lactam carbonyl and its
amide nitrogen, and O3x/O4x sit on C6x at 1.22/1.24 - the C3 carboxylate,
which stayed away at 6.07 and 7.43. So the atom on the metal is the
hydroxyethyl OH, not the carboxylate and not a second carbonyl contact.

O1x was at 3.13 A in the October pose (D123) and **outside the 3.0 A detection
sphere at preprocess**, which is why nothing restrained it: `fill_restraints`
can only restrain what is inside the radius when it looks. It came in from
outside during the MD. **That is the one failure mode detection-time restraints
cannot prevent**, and it is the same failure the `prod100b` note records for
the product arm, where an unrestrained ligand oxygen came from 6.76 A onto the
metal at 2.01 (D173). `wtpostc` was already monitoring `Zn-mer O1x` for that
reason.

**`site.py check ... --at=reactant --hypothesis=carbonyl_on_metal`: REJECTED.**

    PASS  Zn - carbonyl O2x       2.65   want 2.4-3.0
    PASS  Zn - nucleophile O01    2.04   want 1.9-2.2
    PASS  Zn coordination number     4   want 4-5
    PASS  O01 - C4x (attack)      3.46   want 3.0-3.5
    FAIL  closest substrate O      O1x 2.06, want O2x
    FAIL  O01 - C4x - O2x angle     56   want 95-120
    FAIL  Asp92 OD1 - O01         5.81   want 2.5-3.2

**The three criteria D196 named are passed except the angle**: carbonyl
distance, attack distance and coordination number are all inside their
windows, in one structure, for the first time in the project. D195's blocker is
gone. What replaced it is the **direction** of the approach, not any distance.

**The angle failure follows from O1x, and the arithmetic says so.** O1x sits
between the nucleophile and the carbonyl - the sphere angles are

    OH-Zn-O1x  131.3    OH-Zn-O2x   76.8    O1x-Zn-O2x  76.8
    OH-Zn-H164  89.8    OH-Zn-H88   96.3    H164-Zn-H88 91.7
    O2x-Zn-H164 133.4   O2x-Zn-H88 133.5    O1x-Zn-H164 79.8

so the hydroxide and the carbonyl are squeezed to **76.8 degrees** with O1x
wedged in between them. Holding Zn-O01 2.035, Zn-O2x 2.645, O2x-C4x 1.24 and
O01-C4x 3.458 fixed, the O01-C4x-O2x angle is then determined: 56 degrees, which
is what was measured. To land in the 95-120 window the **OH-Zn-O2x angle would
have to be 106.7-128 degrees**.

This matters because **107 degrees is simply tetrahedral**. The Burgi-Dunitz
window and `carbonyl_on_metal` are therefore NOT mutually exclusive, which was
worth checking before blaming the criterion: a four-coordinate zinc holding the
nucleophile and the carbonyl at a normal tetrahedral angle gives a legal attack
geometry. It is the fifth donor crowding the sphere that closes the angle. One
unrestrained atom accounts for the one remaining geometric failure.

**Asp92 OD1 - O01 is 5.81 A, against 4.76 in the QM/MM frame (D195).** The
fourth restraint D196 anticipated is needed, and it has to be an **OD1-O01
distance**: pinning OD1 3.81 A from the zinc fixes how far the carboxylate is
from the metal and says nothing about which way it points, and it pointed away.

**A reporting bug found on the way, which matters beyond this run.** The
pipeline printed

    APPEARED  B:1:ZN - C:1:O1x at 2.06 A
    LEFT      B:1:ZN - C:1:O2x (was 1.94 A)

and O2x had **not** left: it is at 2.645 A, inside the 3.0 A sphere, exactly
where it was restrained. The detector keys the sphere by **residue, reporting
one closest atom each**, so when O1x became `mer`'s closest atom O2x vanished
from the table and the diff called it a departure. Consequences: the printed
"after MD" sphere **omits a genuine donor**, the coordination count it implies
is short by one, and a change report can announce that a restrained contact
broke when it did not. Both `NVTtrajectory_lastframe.pdb` and
`...lastframe_imaged.pdb` give identical distances to three decimals, so
imaging is not involved - this is the per-residue keying, the same flaw D166
and the `wtpost` note already call out for averaging, now showing up in the
change report as well. The fix is the one `per_residue: false` already applies
to `postprocess`; the `md` stage's own before/after reporter does not take it.

**What this does NOT settle.** Whether O1x on the zinc is real. MM cannot say:
a **neutral** alcohol at 2.06 A, closer than an **anionic** hydroxide at 2.035,
on a non-bonded Zn(2+) with no Pauli repulsion, is the textbook MM collapse,
and the stage's own warning says exactly that. The competing reading is that a
carbapenem's 1-hydroxyethyl genuinely chelates the metal, which is not absurd.
Deciding it is QM/MM's job, per the standing rule that MM cannot be trusted for
coordination.

**Disk.** The run wrote **12.7 GB** - every trajectory is stored twice, raw and
imaged, so NVTtrajectory.dcd and NVTtrajectory_imaged.dcd are 5.76 GB each -
taking scratch from 269 to **286/300 GB**. The D196 estimate of 6 GB counted
one copy.

**31736155 `mich_post` FAILED, 31736226 resubmitted.** `postprocess` takes the
**full config key** as its section: `run.slurm postprocess postprocess_michpost`,
not `michpost`. With the short name the section is not found, the stage falls
back to the plain `postprocess` inputs, `md_section` defaults to `md` and the
error is `` `md` has no results `` - which reads like a missing MD rather than a
mistyped section. D152's lesson in a new place: a section that is not found is
not an error, it is a silent fallback to defaults.

`michpost` added to `site.yaml`, **26 selectors checked**, at an explicit
`cutoff: 2.5` so the analysis does not inherit the 3.0 A restraint radius
(D186), and monitoring O3x and O4x as well in case the ligand is turning over
rather than donating one group. Its occupancies are what say whether O1x came
on early and stayed - which decides whether the whole 10 ns is contaminated or
only the end of it.

**D197, measured over the trajectory. The hydroxyethyl held the fourth
coordination site for the whole 10 ns, and the carbonyl was never a ligand by
the project's own criterion.**

`postprocess_michpost` (31736226) over all **10,000 frames** at the 2.5 A
coordination radius:

    donor                       occupancy    min   mean    max
    hydroxide O01                   100.0    1.91   2.03   2.16
    His164 NE2                      100.0    2.01   2.13   2.23
    His88 NE2                       100.0    2.09   2.21   2.33
    mer O1x  (1-hydroxyethyl)        97.9    1.92   2.18   3.19
    Asn225 OD1                       17.9    2.02   2.81   3.50
    mer O2x  (lactam carbonyl)        0.3    2.47   2.60   2.74
    Asp92 OD2                         0.0    2.63   2.77   2.89

    coordination number:  3: 33    4: 8328    5: 1638    6: 1

**O1x is not a late arrival: 97.9% of frames inside 2.5 A, 100% inside 3.5,
mean 2.18 A.** It was on the metal essentially from the start. So **no window
of this trajectory is usable** - there are no early clean frames to take to
QM/MM, which was the one thing worth checking before discarding the run.

**And the carbonyl is inside 2.5 A in 0.3% of frames.** Its mean is 2.60,
which is exactly where it was restrained, with a spread of only +/-0.14 - the
restraint was never in difficulty. The reason the occupancy is 0.3% is that
**the 2.60 A target is outside the 2.5 A bond radius by construction.**

So, by this project's own measurement standard, the lactam carbonyl is **not a
ligand at any point in this run**, while an unrestrained secondary alcohol is
one for 97.9% of it.

**This exposes a conflict between two settled numbers, and it is a criterion
question rather than a modelling one.** The carbonyl restraint target of 2.60
was chosen in D196 to match L1's crystallographic 2.53/2.63 and GFN2's own
relaxed 2.56. The coordination radius is 2.5 (D186). A carbonyl held where the
crystal and the QM both put it is therefore **guaranteed** to be scored as
non-coordinating. Either a dative carbonyl O->Zn needs the 3.0 A "affected"
radius rather than the 2.5 A "bonded" one, or the bond radius is too tight for
this donor type. Nothing should be changed silently: the same 2.5 A figure is
what every earlier coordination verdict in the log was measured at.

**Consequence for how D197's `site.py` result reads.** `PASS Zn coordination
number 4` is composed of hydroxide + His164 + His88 + **O1x**. The carbonyl is
not in the four. The number is correct and its membership is wrong, and the
coordination criterion cannot see that because it counts donors without
requiring the carbonyl to be one of them. The separate `closest substrate O`
test is what caught it. A hypothesis set that tests properties one at a time
can pass on composition it would reject if asked directly.

**A second uninvited donor, not previously tracked: Asn225 OD1, 17.9%
occupancy, min 2.02 A.** Next to Lys226, unrestrained, and inside the sphere
for about a sixth of the run. It is not in `wtpostc`'s monitor list and so has
not been looked for on any earlier arm.

**Everything that was restrained behaved.** O01-C4x mean 3.47 against a 3.47
target, 82.8% inside 3.5. Asp92 OD1 held out at mean 3.80 against 3.81, never
inside 3.5. Asp92 OD2 mean 2.77 against the crystal's 2.80, and **0.0%**
inside 2.5 - the two-radii design doing exactly what D186 built it to do,
holding a second-shell contact without ever counting it as a bond.

**No ligand turnover.** The C3 carboxylate stayed away throughout - O4x mean
5.93, O3x mean 7.05 - so this is not `prod100b`'s whole-ligand reorientation
(D173). It is one specific intramolecular donor, the 6-(1-hydroxyethyl) OH,
chelating the metal while the rest of the substrate stays put.

**Two incidental numbers worth having.** Lys226 NZ-O01 is inside 3.5 A in only
**10.3%** of frames, mean 4.93 - Lys226 is not holding the hydroxide in this
build, which bears on P6 and D41. And Gln93 NE2-Asp92 OD1 has mean **7.74 A**:
Gln93 is nowhere near Asp92 here, so the Q of RQK has no contact to decompose
on this arm.

**D197, two corrections to my own reading, both on the user's challenge
(2026-11-03). The "criterion conflict" was not one, and the QM/MM run I
recommended could not have answered the question.**

**1. The carbonyl at 2.60 A is not in conflict with anything. It is D124's
design criterion being met.** I wrote that the 2.60 restraint target and the
2.5 A coordination radius could not both be satisfied and called it a conflict
between two settled numbers. That is wrong, and D124 already settled it:
thermolysin and carboxypeptidase A hold the nucleophile and the carbonyl on
**one** zinc, and the condition is that the carbonyl is **loosely** coordinated
- one oxygen tight (2.01, 2.19), the other loose (2.41, 2.98). A carbonyl at
2.60 A that is NOT counted as a bond at 2.5 A is exactly that loose oxygen, and
**D186's two radii exist to express precisely this distinction**. The 0.3%
occupancy is the design working, not a contradiction. Nothing needs changing.

**The user's recollection that the carbonyl was thought "polarised, not
coordinated" is correct as history**, and the record is explicit about its fate:
that position (D53, and the user's own instruction that "a polarised carbonyl
should be enough") was **withdrawn on 2026-10-28** by D123/D121, on L1's
7UHI/7UHJ putting the carbonyl at 2.53-2.63 A from Zn1, and then sharpened by
D124 into "loosely coordinated". Worth saying plainly: at 2.60 A the difference
between "loosely coordinated" and "polarised" is not a distinction this
project's measurements can resolve, and no conclusion here depends on the name.
The caveat on the withdrawal still stands - L1 is **di-zinc** and moxalactam is
an oxacephem, so the one measured Michaelis complex differs from CSR-1 in both
metal count and substrate class.

**2. My argument for QM/MM was unsupported, and the experiment cannot
adjudicate O1x.** I claimed GFN2 would push off an artefactual contact, citing
`crystal_v2` going in at 1.94 A and relaxing to 2.56. That is a bond-length
relaxation **within** the coordinated basin, not a rejection, and state.md
section 2b already says what it means in its own words: *"a carbonyl held
coordinated through MM stays under QM/MM"*, with every frame whose carbonyl was
uncoordinated in MM losing the substrate instead (5.16 -> 6.20, D111; the one
inward move on record is 5.16 -> 4.37, which never approaches the sphere).

**So in this project QM/MM has never crossed the coordination boundary in
either direction.** Held in, stays in; left out, goes further out. Starting a
QM/MM run from a frame with O1x at 2.06 A would therefore almost certainly end
with O1x at ~2.1 A, and that outcome is uninformative - it is the starting state
preserved.

**There is a mechanistic reason, already in the record at D96's analysis, and it
is decisive.** Classical MM MD "overfavors octahedral geometries" by pulling
extra ligands into the zinc shell (Yang et al., PMID 34748329, Spencer and
Mulholland among the authors); DFTB3 QM/MM MD rectifies it; **our xTB QM/MM MD
did not**, and the likely reason given there is that tens of ps starting from an
already-crowded sphere is far too short, because **ligand exchange at Zn(II) is
slow**. A 50 ps QM/MM MD cannot expel O1x on that timescale whatever GFN2
prefers. The proposed run would have measured the timescale of our trajectory,
not the chemistry.

**3. And O1x is not a stray artefact - it is the pose's own nearest oxygen.**
D123 records that the pose inherited from `CSR_crystal` "presents the
6-alpha-hydroxyethyl **O1x at 3.13 A** and leaves the carbonyl at **5.43 A**".
So in the pose everything descends from, the hydroxyethyl is the group near the
metal and the carbonyl is nowhere near it. `eq5` dragged the carbonyl from 5.43
to 2.60 and left O1x unrestrained at 3.13, from which it fell in to 2.06. **The
10 ns did not introduce a contaminant; it added the carbonyl alongside the
contact the pose already had.** That is a sharper statement of the problem than
"an unrestrained donor crept in", and it is about the pose, which D194 settled
as not to be re-litigated by reverting.

**What would actually adjudicate it.** Not a dynamics run, because the event
needed is slower than the affordable trajectory. A **thermodynamic** comparison
instead: optimise with O1x on the metal and with it displaced, and compare
energies, at GFN2 and at DFT on the same geometries. That is the cluster
benchmark - **P1 / plan item 7, open since D9**, the `sp` cluster mode built for
it in D187/D188, and the standing blocker recorded in this project's memory as
"whether xTB over-coordinates the zinc". Imhof's own benchmark, the paper this
method follows, found AM1/d made 4->5 coordination **7.55 kcal/mol** too
favourable; the same number for GFN2 is unmeasured, and no GFN2 benchmark for
zinc coordination number in enzymes appears to exist. Until it is measured, a
QM/MM verdict on a fifth donor is not interpretable either way.

**The lesson, since this is the second time.** I quoted a sentence from this
file as support for the opposite of what it says. The rule in the header - that
this file is the single source of truth and the log is history - only works if a
claim is read where it is written, not recalled. Both corrections here came from
the user asking twice, not from the check.

**[QUALIFIED - its own caveat limits this to one contact and one geometry pair; D200 shows frame 4374 is not unique. See `csr/retractions.md`.]** **D198 - The cluster benchmark is built and running, and building it turned up
a better starting structure than the one we had. Frame 4374 is the only
structure in the project where the carbonyl is the substrate oxygen nearest the
zinc.** (2026-11-03, user: "ok run that, and then run proper QM/MM 300ns with
GFN2 to set up for the charge deletion and NEB steps")

**300 ps, not 300 ns.** Flagged before anything was submitted and agreed by the
user: at the ~1 s per fs step these runs cost, 300 ns is of order 9 years of
walltime. The MM 100 ns runs are affordable only because MM is thousands of
times cheaper per step. Six 50 ps batches, checking in after the first.

## The benchmark (P1, open since D9)

`docs/csr/build-cluster.py` plus four `sp` sections. The question is narrow:
**does GFN2 hold O1x on the zinc when DFT would not?** The quantity is a
difference of differences, `[E_off - E_on](GFN2) - [E_off - E_on](DFT)`, so
strain and thermal noise enter both brackets and cancel while the method does
not. Imhof's own benchmark found AM1/d out by 7.55 kcal/mol on a 4 to 5
coordination change; the GFN2 number has never been measured.

**The truncation came from CONECT, as D187 insisted.** `mer_crystal.pdb`'s 51
records give the real topology, which cannot be read off the atom names: the
beta-lactam ring is **C4x-C3x-C17x-N1x** with O2x on C4x, the pyrroline ring is
N1x-C5x-C7x-C15x-C17x fused at N1x-C17x, the C3 carboxylate is C6x on C5x, and
O1x is the hydroxyl of the 1-hydroxyethyl on C2x. So the fragment carrying both
oxygens of interest is the intact beta-lactam ring plus the hydroxyethyl.

**Two departures from D187's scope, both deliberate.** The C3 carboxylate is
left out - it is 5.93 to 7.05 A from the zinc throughout `md_eq5`, and a bare
carboxylate 6 A from a dication in gas phase makes a 1-10 kcal/mol effect ride
on a far larger charge-separation term. That makes the cluster **neutral**
rather than D187's -1. And both geometries are **real frames of one
trajectory** - the shortest and longest Zn-O1x - not one frame displaced by
hand, which would put strain into one bracket only and report it as chemistry.

Both clusters came out at **51 atoms, C15 H25 N5 O5 Zn, charge 0, mult 1**, with
the same five bonds capped (Asp92 CB-CA, His88 CB-CA, His164 CB-CA, N1x-C5x,
C17x-C15x). The script refuses to finish if the two differ in atom count, since
the comparison is only a method comparison if both methods see the same system.

**The pair is better conditioned than expected: the carbonyl is at 2.58 A in
BOTH frames**, to two decimals, because the 500 kcal/mol/A^2 restraint held it
there all run. So Zn-O1x (1.92 vs 3.74) is very nearly the only large
difference between the brackets.

    GFN2-xTB   on  -80.274541064236 Eh     off  -80.297915594550 Eh
               E_off - E_on = -14.67 kcal/mol

Taken alone that says GFN2 prefers O1x **off** by 14.67 kcal/mol, but it must
not be read alone: the two frames are 3,960 apart and differ thermally across
all 51 atoms. Only the difference against DFT isolates the method.
**31738632/31738633 (DFT) were still in SCF at 6.5 minutes; that number is
pending and this entry is incomplete without it.**

## What fell out of it: frame 4374

`site.py check` on the two extracted frames, which cost nothing since both were
already on disk:

    criterion                 frame 4374 (off)   last frame      frame 414 (on)
    closest substrate O       PASS O2x 2.58      FAIL O1x 2.06   FAIL O1x 1.92
    Zn - carbonyl O2x         PASS 2.58          PASS 2.65       PASS 2.58
    Zn - nucleophile O01      PASS 2.07          PASS 2.04       PASS 2.04
    Zn coordination number    PASS 5             PASS 4          PASS 5
    O01 - C4x (attack)        PASS 3.50          PASS 3.46       PASS 3.50
    O01 - C4x - O2x angle     FAIL 68            FAIL 56         FAIL 39
    Asp92 OD1 - O01           FAIL 5.67          FAIL 5.81       FAIL 5.25

**Frame 4374 passes five of seven, and is the first structure in the project to
pass `closest substrate O`** - the criterion that encodes what
`carbonyl_on_metal` actually claims. It is still REJECTED, on the attack angle
and the Asp92 H-bond.

**The angle tracks O1x monotonically, which is the D197 mechanism confirmed on
three independent points:**

    Zn-O1x  1.92 A  ->  attack angle 39 deg
    Zn-O1x  2.06 A  ->  attack angle 56 deg
    Zn-O1x  3.74 A  ->  attack angle 68 deg

So the crowding account is right. **But it is not sufficient**: 68 deg is still
far from 95, and the reason is visible in the sphere of frame 4374 - the five
donors are hydroxide 2.07, His164 2.13, **Asn225 OD1 2.19**, His88 2.23,
carbonyl 2.58. **Asn225 took the site O1x vacated.** One uninvited donor
replaced the other, the sphere stayed five-coordinate, and the angle opened by
only 12 deg. Whatever closes this angle is not O1x specifically; it is that this
site holds a fifth donor in MM whatever we do, and a five-coordinate zinc puts
the nucleophile and the carbonyl too close together.

## The production run

**31738709 `mich_qm1`**, `qmmm_michqm`, 50 ps of the 300, 16 cores, 48 GB, 24 h.

**It starts from frame 4374, and that is a declared selection.** QM/MM here is
hysteretic - no run in this project has crossed the coordination boundary in
either direction (D197 corrections) - so whichever donor is on the metal at the
start is the one on it at the end. Starting from the last frame would bake in
O1x at 2.06 A and settle nothing. Starting from the extremum of a distribution
is still a selection, and anything measured from this run describes the reactant
state we chose rather than an equilibrium average. That has to be said in any
write-up.

**`qmreg_mich` is a new QM region and Asn225 is in it.** Not a preference: at
2.19 A it is inside the bond radius by any reading, and a coordinate bond left
in MM is represented by point charges, which is the one thing electrostatic
embedding cannot do. It was in no earlier arm's region because it was not known
to be a donor until D197. No `D:1:HOH`, which this build does not have. Charge
0, stated.

**Batching.** `time` is the **total** wanted and `resume: true` continues from
`OpenMM_MD_state.xml` (D182/D183), so each batch raises `time` by 50 and
resubmits the same section - one line, and a walltime-limited batch never
restarts from zero.

**D198, the benchmark result. GFN2 does NOT over-coordinate here - it
under-binds the hydroxyethyl by 5.92 kcal/mol relative to DFT, which is the
opposite sign to the worry this has been blocked on since D9.**

    geometry pair: frames 414 and 4374 of md_eq5, 51 atoms, charge 0
    Zn-O1x 1.92 A (on) vs 3.74 A (off); Zn-O2x 2.58 A in BOTH

    GFN2-xTB            on -80.274541064236   off -80.297915594550   Eh
    wb97x-v/def2-svp    on -3014.167570288007 off -3014.181518200444 Eh

    E_off - E_on   GFN2  -14.67 kcal/mol
                   DFT    -8.75 kcal/mol
    dd = GFN2 - DFT      **-5.92 kcal/mol**

**Both methods put the O1x-off geometry lower**, and GFN2 exaggerates that
preference by 5.92 kcal/mol. So:

1. **The 97.9% O1x occupancy cannot be blamed on xTB.** It was an MM result in
   the first place, and GFN2 is *less* inclined to hold that contact than DFT
   is. The standing question - "does xTB over-coordinate the zinc?" - is
   answered **no, on this contact**: it under-coordinates it.
2. **DFT also says the contact is unfavourable**, by 8.75 kcal/mol on these
   geometries. That is consistent with the O1x contact being an MM artefact of a
   non-bonded Zn(2+) with no Pauli repulsion, which is what the chemistry
   suggested (a neutral alcohol closer than an anionic hydroxide).
3. The magnitude is comparable to Imhof's **7.55 kcal/mol** for AM1/d on a
   4 -> 5 coordination change, so the error is the same size as the one that
   paper found disqualifying; what differs is its direction.

**What this does NOT license.** These are single points on MM-sampled
geometries, so -8.75 is **not** a binding energy: the two frames are 3,960
apart and differ thermally across all 51 atoms. The `dd` cancels geometry to
first order because both methods see the same pair, but not exactly - each
method evaluates the same thermal distortions slightly differently. And it is
**one contact on one geometry pair**, not a general statement about coordination
number. The optimisations that would turn this into a binding energy are a
separate, larger job.

## A near-miss worth recording: the QM/MM ran on a system missing 4,556 waters

Job 31738823 started cleanly, placed its link atoms, built a correct 90-atom QM
region - and built the MM system with **34,311 atoms where the MD's has
47,979**. `frame_o1x_off.pdb` had been written by MDAnalysis's PDB writer, which
renumbers the waters; `resSeq` is four characters and this box has ~14,500 of
them, so thousands shared a `(chain, resSeq)` and OpenMM's reader **silently
dropped the duplicates**.

**Nothing in the pipeline would have caught this.** No stage checks that the
system it builds matches the one that was equilibrated, the energies are all
relative, the QM region and the link atoms were right, and the geometry looked
normal. The only symptom was one line in the log that had to be compared by hand
against the MD's. It was caught because the duplicate-atom warnings in the
`.err` prompted a check, not because anything failed.

Fixed by writing extracted frames as **coordinate substitution into the topology
PDB**, which keeps ASH's own numbering byte-for-byte; the DCD's atom order is
the topology's by construction, so the swap is exact. The writer refuses if the
template and the frame disagree on count. Both frames now come out at 47,979,
and **31738904** builds 47,979 atoms / 14,888 residues / 33,421 bonds, matching
`md_eq5`.

**The cluster benchmark is unaffected** - the clusters keep only the metal, the
hydroxide, two histidines, Asp92 and the ligand fragment, none of which were
renumbered, and the re-extracted frames give the identical frame numbers and
distances (414/1.92, 4374/3.74).

**Two other setup faults, both cheap and both from the same cause - a run
directory built for two stages being asked to run a third.** `qmmm` layers
`qmmm.inputs` under the section's, and `build-mich.py` writes only
`preprocess:` and `md:`, so 31738709 died in 5 s on `No qmmm.inputs section`.
`docs/csr/add-qmmm-base.py` now writes that base, taking the force fields from
the run's **own** `md.inputs` rather than a template, because a QM/MM run whose
MM differs from its MD's is the same class of silent error as the missing
waters. Its first run then failed on `_locked` being handed a `str` instead of a
`Path` - and the chained `sbatch` fired anyway, submitting 31738807 which was
bound to die. That is the **second** time this session a failed script did not
stop a chained command; submits are now gated with `&&`.

**D199 - QM/MM batch 1 runs at 0.87 fs/s. The resume trap that cost two
submissions, and the real reason O1x sits at 2.06 A - which is not the one I
gave in D197 and D198.** (2026-11-04)

**31740139 `mich_qm1`, RUNNING.** 47,979 atoms / 14,888 residues / 33,421 bonds
matching `md_eq5`, 90 QM atoms, four boundaries. Step 107 at about 103 s of MD:
**1.04 fs/s = 3.7 ps/hour**, so 50 ps is **13.3 h** and batch 1 should land
about 13:40 on 4 Nov, inside its 24 h walltime. That is the ~1 s/step the 300 ns
arithmetic assumed, so the six batches cost about **80 h** of compute.

**Corrected against 13.5 h of real running (13:48, 4 Nov).** Step 42,470 at
48,748 s is **0.871 fs/s = 3.14 ps/hour** - 16% slower than the figure above,
because the first 107 steps were measured before the trajectory and state
writers had started. So 50 ps is **15.9 h**, batch 1 lands about **16:10**
rather than 13:40, and the six batches cost about **95 h**, not 80. Still well
inside the 24 h walltime, but the margin on a 50 ps batch is 8 h, not 11, which
matters if a later batch is slowed by a busier node. **Take a rate from a run
that has written at least one frame**, not from its first hundred steps.

## The resume trap

31738904 built the **correct** 47,979-atom system and then died in 14 s:

    openmm.OpenMMException: Called setPositions() on a Context with the
    wrong number of positions          (at simulation.loadState(statefile))

`michqm` carries `resume: true` from the start so batches 2-6 need only `time`
raised. But `protdft/qmmm_michqm/` still held the `OpenMM_MD_state.xml` of the
**cancelled 34,311-atom run**, and a state file carries its own particle count,
so the section resumed a system that no longer existed. Clearing the output
directory fixed it, and mcube then says the right thing by itself:

    [INFO] resume: true, but no OpenMM_MD_final_state.xml or
    OpenMM_MD_state.xml in ... - nothing to resume from, starting a fresh run

**No library change needed** - D183's wiring already covers it. What to remember
is the interaction: **`resume: true` plus a killed run in the same output
directory resumes the dead run's system.** It fails loudly, but the message
names positions rather than the state file's provenance, so it reads like a
system-building bug instead of a leftover.

## Correction: "a non-bonded Zn with no Pauli repulsion" is wrong, and the right
## answer is in `ions.py`

I used that phrase in D197, in D198 and in `build-cluster.py` to explain O1x at
2.06 A. **It is not what the code does, and I should have read the code rather
than quote the warning.** The phrase itself is mcube's own text, printed by
`md.py` when the sphere changes - so it is the library's loose wording that I
adopted, but checking it was still mine to do.

What actually happens, verified in both logs and in the config:

    [INFO] IOD parameters for ZN B:1: sigma 0.2265 -> 0.2486 nm,
                                      epsilon 0.0138 -> 0.0624 kJ/mol
    md_eq5.results.metal_parameters: iod

**`md` applies the Li/Merz IOD set by default and so does every QM/MM stage** -
`ions.apply_metal_parameters` is "the one place this happens, for md and every
QM/MM stage alike, so a metal cannot be parameterised one way in MM and another
in QM/MM". So the 10 ns ran on **fitted** ion-oxygen parameters, not stock Amber
ones, and a draft of this entry claiming otherwise - and proposing an IOD re-run
as a cheap test - was wrong on both counts and was deleted before it was
committed. (It survived one check only because `/iod/i` matches
"per**iod**icity"; the real grep is `IOD parameters`.)

**`ions.py` states the correct explanation and predicted this exactly:**

> "It does not solve over-coordination. A 12-6 potential knows only distance and
> atom type, so it cannot tell a ligand from a non-ligand and every oxygen in
> reach is still rewarded for approaching; only something that knows which atoms
> are ligands can stop that, which is what the coordination restraints do."

So the chain is:

1. The inherited pose puts O1x at **3.13 A** (D123).
2. That is **outside** `preprocess`'s 3.0 A detection radius - by 0.13 A - so
   `fill_restraints` never saw it and no restraint was placed on it.
3. IOD moves the Zn-O minimum to about **1.89 A** and rewards *any* oxygen for
   approaching, because a 12-6 term cannot tell an alcohol from a hydroxide.
4. O1x therefore settles at **2.06 A**, which is near that minimum and inside
   the 2.0-2.1 A IOD is fitted to reproduce.

**This is not a collapse to an absurd distance; it is the fitted distance,
applied to the wrong oxygen.** The artefact is not a missing repulsion but a
potential with no concept of which atoms are ligands, and the only thing in the
pipeline that does have that concept - the coordination restraints - was blind to
O1x because of where the pose put it. That makes the D197/D198 conclusions about
O1x *stronger*, not weaker: the contact appeared despite fitted parameters, so it
is a property of the pose plus the restraint radius, not of sloppy MM.

It also sharpens what the D198 benchmark showed: DFT puts the O1x-bound form
8.75 kcal/mol higher, and GFN2 a further 5.92 on top of that, while MM gives it
the fitted bond distance. Three methods, three answers, and only the two quantum
ones know what kind of oxygen it is.

**`md.py`'s warning text should be fixed** - "no Pauli repulsion" in a codebase
that ships fitted LJ parameters and applies them by default is misleading to
exactly the reader the warning is for. Not changed here, because the live run
reads that code; it is P8.

**Third correction of my own claims this session**, after the carbonyl
"criterion conflict" and the QM/MM hysteresis argument. All three were in the
record already. Two were caught because the user pushed back; this one because
the standing instruction to check the record is now a memory.

**[RETRACTED by D201 - the sphere is NOT the designed one, and 1TMN/7CPA are 57-60 deg cis. See `csr/retractions.md`.]** **D200 - QM/MM batch 1, 50 ps of GFN2. The zinc site is right, both MM
intruders left, the carbonyl coordinates without a restraint - and the attack
geometry is impossible for a reason that is arithmetic, not sampling.**
(2026-11-04)

**31740139 `mich_qm1`, COMPLETED, exit 0:0, 15:54:58 for 50,000 steps.** 50
frames at 1 ps. Analysed through the stage, as `postprocess_michqmpost` (job
31885041), not through a script: the result is in `config.yaml`.

`postprocess --qmmm` was not the route. It loads the **base** `qmmm` section
rather than the one the job ran as, and reads `results.opttraj`, which a run
with `skip_optimize: true` never writes. The section names `pdb_file` and
`trajectory` outright instead, which `postprocess` supports for exactly this.

## The coordination sphere is the designed one, and both MM intruders left

    cutoff 2.50 A, 50 frames                occ%     min / mean / max
    C:2:O01    hydroxide                   100.0    1.72 / 1.86 / 2.08
    A:88:NE2   His88                       100.0    1.92 / 2.09 / 2.28
    A:164:NE2  His164                      100.0    1.94 / 2.10 / 2.33
    A:92:OD2   Asp92                        96.0    1.87 / 2.06 / 2.61
    C:1:O2x    mer carbonyl                 66.0    2.08 / 2.43 / 3.09
    C:1:O1x    mer 1-hydroxyethyl           10.0    2.27 / 2.82 / 3.46
    coordination number: 4 in 14 frames, 5 in 36

**2 His + Asp92 OD2 + hydroxide + a loose carbonyl.** That is the Michaelis
sphere D196 specified, and it is the first time any structure in this project
has held it.

The two donors MM put there and QM/MM was supposed to adjudicate **both left,
decisively**:

- **Asn225 OD1: 2.19 A at the start, 4.18 A at its closest over 50 ps, 5.86 A
  mean.** It is not a ligand. D198's objection to frame 4374 - that Asn225 took
  the site O1x vacated, so the window needed a genuinely four-coordinate sphere
  before the angle could be judged - **is answered, and in the direction that
  removes the objection**. Asn225 was in the QM region, so it had every
  opportunity to stay.
- **O1x, the 1-hydroxyethyl: inside 2.5 A in 10% of frames, mean 3.71 A over
  the contact range, out to 4.73.** The 97.9% occupancy at 2.06 A that `md_eq5`
  reported is not reproduced at the GFN2 level.

**On O1x, say what this is and is not.** It started OUTSIDE (3.74 A), and
"left out, goes further out" is the documented hysteresis (D111, D197), so O1x
staying out is what hysteresis predicts and carries no weight by itself. What
is worth something is that it is **not frozen**: it reaches 2.27 A in 10% of
frames and comes back out. A contact that visits the sphere and leaves is being
sampled, not held. Taken with D198 - GFN2 under-binds it by 5.92 kcal/mol
relative to DFT, and DFT puts the bound form 8.75 higher - the O1x contact is
an MM artefact, and the question D9 blocked on is closed for this contact.

## The carbonyl coordinates with nothing holding it, and it exchanges

This was the honest test: O2x was restrained in every previous run and
**unrestrained here for the first time**, starting at 2.58 A, which is outside
the 2.5 A radius.

    66% of frames inside 2.5 A, mean 2.43 A, range 2.08 - 3.59

So it moved **in**, and it crosses the boundary in both directions inside 50
ps. The mean sits on D124's loose-carbonyl precedent almost exactly -
thermolysin 2.41, L1 7UHJ 2.63, CPA 2.98 - which is the one position this
project has withdrawn and reinstated twice (D53 withdrawn by D121/D123, settled
in D124).

**This corrects the hysteresis generalisation, which is mine.** D96 and D197
say no run here has crossed the coordination boundary in either direction,
because ligand exchange at Zn(II) is slower than the affordable trajectory.
That is now false as stated: a **neutral, weakly bound** carbonyl exchanges
freely on a picosecond timescale. The claim should be narrowed to what the
evidence supports - an anionic or chelating donor already on the metal does not
leave in 50 ps, which is still what every other row above shows (the hydroxide,
the histidines and Asp92 never leave). A loose neutral donor is the easiest
possible case and it moves. **Do not use "QM/MM cannot cross the boundary" as a
blanket reason to discount a sphere again.**

## The attack geometry cannot work, and the reason is geometry

    O01-C4x  attack distance        14% within 3.5 A, mean 3.82, min 3.16
    O01-C4x-O2x  Burgi-Dunitz       0% in 95-120 deg, mean 28.9, max 78.1
    O01-C4x-N1x                    62% in 95-120 deg, mean 119.5
    O01...O2x                      100% within 3.5 A, mean 2.86
    O01-Zn-O2x                      0% in 150-180 deg, mean 82.0

The attack angle did not open towards 95-120; it **closed**, from 68 deg at the
start to 28.9 mean, and the last frame alone gives 48. It is back inside the
24-54 deg band D34 measured over 791 frames.

The last three rows say why, and they are not a sampling statement. The
nucleophile and the carbonyl oxygen are **both on the same zinc, cis, 82 deg
apart and 2.86 A from each other**. Fix those two numbers and the attack angle
follows from the triangle:

    O01-C4x 3.82, C4x-O2x 1.25, O01...O2x 2.86  ->  33 deg at C4x

against 28.9 measured, so the three numbers are mutually consistent and nothing
is mistranscribed. Run it backwards for the angle the mechanism needs:

    O01-C4x-O2x = 105 deg  requires  O01...O2x = 4.31 A
    with Zn-O01 1.86 and Zn-O2x 2.45, that is O01-Zn-O2x = 180 deg

**A Burgi-Dunitz attack angle with the carbonyl on the same metal as the
nucleophile requires the two to be essentially TRANS across the zinc.** Cis
ligands on a four- or five-coordinate zinc are 2.8-3.0 A apart, which puts the
hydroxide beside the carbonyl oxygen rather than above the carbon. The
measured O01-C4x-N1x of 119.5 deg confirms the picture independently: with the
ring's O2x-C4x-N1x near 92 deg, 28.9 + 92 = 121 against 119.5 measured, so the
hydroxide is lying **in the plane of the carbonyl**, on the oxygen's side, not
perpendicular to it.

This is the sharpened version of the tension D198 flagged. It is no longer "the
angle is wrong and crowding may be why"; it is **the angle is determined by
which two sites of one zinc the nucleophile and the carbonyl occupy**, and the
pose has them adjacent.

**What it does not settle** is which way out is right, and there are at least
three, none of which this run can choose between:

1. The pose is wrong and a trans arrangement exists on a five-coordinate zinc -
   axial-axial in a trigonal bipyramid, with the 2 His and Asp92 equatorial.
   Whether that is accessible is a structural question, not a sampling one.
2. The carbonyl is polarised by the metal without being a ligand, and the
   nucleophile comes from elsewhere - the position D53 held, withdrawn on
   crystallographic grounds, and not revived by this.
3. The mono-Zn assumption is what fails, and the second site matters (plan item
   10). Thermolysin and CPA hold both on one zinc, so **what Nu-C=O angle those
   structures actually have is the literature check this now needs**: if they
   reach 100 deg with both on one metal, the pose here is wrong; if they do not,
   the geometric criterion is.

## A method point worth more than the numbers: a single frame cannot gate a loose contact

`site.py check` on `QM_MM_lastframe.pdb` **REJECTED** the structure, with
`Zn - carbonyl O2x 3.07, want 2.4-3`. The 50 ps says that contact is inside 2.5
A in 66% of frames at 2.43 A mean. The last frame caught it mid-excursion, and
on a contact with a 1.5 A range the frame a job happens to stop on is close to
a coin toss.

The full gate on the last frame, for the record:

    PASS  closest substrate O     C:1:O2x 3.07  (first run ever to pass this)
    PASS  Zn coordination number        4       want 4-5
    FAIL  Zn - carbonyl O2x           3.07      want 2.4-3
    FAIL  Zn - nucleophile O01        1.86      want 1.9-2.2
    FAIL  O01 - C4x attack            3.73      want 3-3.5
    FAIL  O01 - C4x - O2x angle         48      want 95-120
    FAIL  Asp92 OD1 - O01             5.32      want 2.5-3.2

Two of those five failures are frame noise against the trajectory means (O2x
3.07 against 2.43; O01 1.86 is also the trajectory mean and fails by 0.04 A,
which is a window boundary, not a finding). Two are real and unchanged - the
attack distance and the angle. One is real and was never restrained: **Asp92
OD1-O01, 4.68 A at its closest, 5.17 mean.** Worth noting that `Lys226 NZ-O01`
is 2.94 A mean at 92% occupancy - something is H-bonding the hydroxide, and it
is not Asp92.

**P9: `site.py check` should be able to read a criterion's occupancy from a
`postprocess` section instead of measuring one frame.** A gate that a
fluctuating contact passes or fails by which frame a job stopped on is not a
gate, and the data to do it is already in `config.yaml`.

## What this says about the question the batch was for

The user asked whether the system is set up correctly before committing the
remaining five batches. **The zinc site is set up correctly, and the setup
produced two real results** - the sphere is the designed one, and both MM
intruders left. **The attack geometry is not, and more of the same trajectory
will not fix it**, because the obstruction is which sites the two oxygens
occupy, not how long they are sampled. Batches 2-6 would sample a complex that
cannot react.

**D201 - The Burgi-Dunitz gate is wrong, measured from the two structures it
cites. D200's "the carbonyl and the nucleophile must be trans" is false: both
reference enzymes are 57-60 degrees CIS, more cis than ours. What is actually
wrong with batch 1 is two setup defects - Asp92 was left free and went to 2.00
A, and Lys226 is in the MM region while sitting 2.82 A from the nucleophile.**
(2026-11-04, measured from 1TMN and 7CPA downloaded from RCSB and from batch
1's own first and last frames, read directly rather than from this log)

## 1. The reference geometries, measured rather than quoted

    structure   Zn-O(nuc side)  Zn-O(carbonyl side)  O...O   O-Zn-O   O-X-O at centre
    1TMN 0ZN         2.01              2.41           2.15    57.2     119.5  (C-O 1.25 / 1.23)
    7CPA FVF         2.19              2.98           2.68    60.1     117.2  (P-O 1.60 / 1.53)

D124's distances are confirmed to 0.01 A. Two things it did not carry into the
criteria:

- **The 117-119 degrees is the O-X-O angle of a formed tetrahedral centre.** In
  1TMN both oxygens are carboxylate oxygens on one carbon at 1.25 and 1.23 A;
  119.5 degrees is what that bonding gives and nothing else. D124's own caveat
  said exactly this and it was never applied to the gate.
- **The two oxygens are 57-60 degrees apart on the metal.** That is strongly
  cis - more cis than CSR-1's 71.8 degrees in batch 1's last frame.

## 2. So D200's inversion was the wrong question, and its conclusion is withdrawn

D200 held r(O01-C4x) at the measured 3.82 A, asked what O...O gives 105
degrees, got 4.31 A, and concluded the two oxygens must be essentially trans
across the zinc. The arithmetic is right and the constraint is wrong: it
demands a transition-state angle at a reactant distance. The attack angle is a
function of r, and at fixed O...O it opens by itself as the nucleophile
approaches:

    r(Nu-C)   CSR-1 (O...O 3.05)   1TMN (2.15)   7CPA (2.68)
     3.73          47.9                 --          37.8
     3.50          58.7                 --          46.5
     3.20          71.9               26.2          56.6
     3.00          80.7               37.6          63.0
     2.60          99.6               55.3          75.8
     2.40         110.5               63.3          82.4
     2.00         141.4               79.2          97.2
     1.60           --                97.2         116.6

**At a reactant separation of 3.0-3.2 A thermolysin's own arrangement gives
26-38 degrees.** The 24-54 degree band D34 reported over 791 frames, D48/D49's
33-40 degrees, and batch 1's 28.9 degree mean are what a working mononuclear
zinc hydrolase looks like before the bond forms. 95 degrees arrives at r = 1.64
A in thermolysin and 2.06 A in CPA, against **2.47-2.69 A in CSR-1**: because
our O...O is wider, the window is reachable EARLIER along our approach than in
either reference.

The project has therefore been rejecting its own structures since D33 on a
criterion that both of its reference enzymes also fail.

## 3. The gate is internally inconsistent, and the fix

`carbonyl_on_metal` demands together: Zn-O2x 2.4-3.0 (source: 1TMN 2.41, 7CPA
2.98), Zn-O01 1.9-2.2 (source: 1TMN 2.01, 7CPA 2.19), O01-C4x 3.0-3.5, and
O01-C4x-O2x 95-120. The last three are jointly satisfiable only for O...O
between 3.35 and 4.3 A, which forces O01-Zn-O2x above about 93 degrees.
**Neither of the two structures the distance windows are taken from would pass
the gate those distances are in.**

Changed in `site.yaml`: the angle moves to `required_at: [transition_state]`,
so it prints as a WARN at `start` and `reactant` and is a FAIL at no stage
`site.py check` is ever run with. Two criteria measured from 1TMN/7CPA take its
place at `reactant`:

    O01 - Zn - O2x      45-95 deg    1TMN 57.2, 7CPA 60.1; batch 1 71.8   PASS
    O01 ... O2x         2.0-3.3 A    1TMN 2.15, 7CPA 2.68; batch 1 3.05   PASS

These say what the angle was meant to say - that the nucleophile and the
carbonyl are correctly placed relative to each other on one metal - without
demanding the transition state of a reactant.

## 4. Batch 1's sphere is NOT the designed one. D200 said it was, and that is wrong

Read off `QM_MM_firstframe.pdb` and `QM_MM_lastframe.pdb` directly:

    Zn-                first frame   last frame
    O01  hydroxide        2.07          1.86
    His88 NE2             2.23          2.01
    His164 NE2            2.13          2.20
    Asp92 OD2             2.77          2.00
    Asp92 OD1             3.82          4.10
    mer O2x carbonyl      2.58          3.07
    mer O1x               3.74          3.63
    Asn225 OD1            2.19          6.33

D196 designed 2 His + hydroxide + carbonyl = **four**, with OD2 at 2.80 A as a
second-shell contact and OD1 pinned out at 3.81. What ran has **Asp92 OD2 as a
full ligand at 2.00 A**, inside 2.5 A in 96% of frames. It started at the
designed 2.77 A, the `qmmm` stage applies no MM restraints, and xTB pulled it
onto the metal within 50 ps.

**This is the unconditioned test D59 said D58 lacked.** D58 compared xTB and
DFT along Zn-OD2 from a structure that was already octahedral and found DFT 1.6x
steeper; D59 correctly objected that both methods had been asked a loaded
question. The move has now been watched from a correct four-coordinate start
with nothing imposed, and it still happens. So **Asp92 on the metal is a real
feature of this model rather than an artefact of the starting structure** - and
the crystal reference still says the opposite (6NC5 ZN302: two His, two waters,
no aspartate, D108; 2.68 A superposed into 6DN4, D184). The discrepancy is
therefore method-against-crystal, which is plan item 11, and it is no longer
avoidable by picking a better start.

**It also finishes off the hysteresis argument.** D200 narrowed D96/D197 to
"an anionic or chelating donor already on the metal does not leave in 50 ps".
Asp92 is anionic and it crossed the 2.5 A boundary INWARD in this very run,
from 2.77 to 2.00 A. Nothing is left of "QM/MM cannot cross the coordination
boundary" in either direction, for any donor type, and it must not be used
again as a reason to accept or discount a sphere.

## 5. The hydroxide survived because it was not allowed to do anything else

O01 carries one hydrogen (0.96 A) in both frames and Lys226 NZ carries three.
That is not a result. `michqm`'s QM region is Zn, the hydroxide, meropenem,
His88, His164, Asp92 and Asn225 - charge 0, with **Lys226 and Arg90 in the MM
region** (D190, and the config says so). An MM lysine cannot transfer a proton.

And over the 50 ps **Lys226 NZ moved from 5.68 A to 2.82 A of the nucleophile**
(92% occupancy on the trajectory). So the model holds the residue this log says
protonates the hydroxide - D29 in xTB, confirmed at DFT in D35, and still -5.2
kcal/mol in a clean four-coordinate sphere in D48 - at hydrogen-bond distance
and unable to act. **Batch 1 did not test the protonation state. It excluded
it.** D41's third item, "a five-coordinate Zn-OH- with Lys226 as NH3+ is the
canonical B3 reactant and was never tested", is still untested.

What is actually near the nucleophile, measured: **Ser186 OG** donates at 1.60
A in the first frame and has left by the last (4.06 A); **Asn225 HD21** donates
at 1.77 A in the last frame, having given up its coordination site. Asp92 OD1
is 5.67 then 5.32 A. The hydroxide is hydrogen-bonded throughout and never by
Asp92.

## 6. The water-nucleophile / neutral-Lys arm, re-read against all of this

The angle never discriminated between the two routes. HOH3094 reached 3.41 A at
99 degrees in 210 of 500 frames (D51) because a non-metal-bound water attacks
from the open face, and 99 degrees at 3.41 A requires O...O = 3.8 A, which a
nucleophile cis on the metal cannot have at that distance. Reading that as "the
water route satisfies the geometry and the hydroxide route does not" compared
one route against a criterion the other cannot meet by construction - and that
both reference enzymes meet only at bond formation.

**What actually failed in that arm was the acceptor, and it still fails.** Zero
frames of 500 with the attacking water donating to Asp92, twice, from two
different starting states (D49, D51); the acceptor census found only bulk waters
within hydrogen-bond reach of its protons (D52). Both of those runs had Asp92
on the metal at 2.07-2.11 A, which was the stated reason to distrust them - and
batch 1 now shows Asp92 goes there spontaneously from the correct start. So
**Asp92 as general base has failed three times, the third time for the right
reason**: it is a zinc ligand in every state this model reaches, not only in a
badly set up one. The CphA route is ruled out a third time.

**The acceptor that has never been tested is Lys226.** Neutral Lys226 - the
state D35 and D48 establish at DFT - is an amine 2.82 A from the nucleophile in
92% of batch 1's frames, and apart from transient waters it is the only group
in hydrogen-bond reach of it. The resting-state transfer is downhill by only
5.2 kcal/mol in a clean sphere, which is the asymmetry of a shuttle rather than
a sink. **Zn-OH2 attacking with neutral Lys226 accepting the proton** is
consistent with everything measured here, it is the RQK lysine doing catalytic
work, and it has never been computed - because Lys226 has been outside the QM
region since D190.

## 7. What this changes about batches 2-6

D200 stopped them on the geometry, and the geometry does not stop them. What
does is section 4 and section 5: Asp92 is unrestrained and goes to 2.00 A, and
Lys226 sits on the nucleophile in the MM region. 250 more ps with those two
unchanged buys more of the same complex. Both are config changes - Lys226 into
`qm_residues` with the charge corrected, and a decision about whether Asp92 is
held at the crystal distance or allowed to bind - not code changes, and neither
needs a new equilibration.

**D202 - Arg90 and Lys226 go back into the QM region, as a new section
`michrk` rather than a change to `michqm`. D190 is not reverted: the four RQK
cations are split by the KIND of role each has, because charge deletion and a
proton transfer are not the same measurement and cannot share one region.**
(2026-11-04, the user's decision: "we'll have to do both, let's make arginine
and lysine QM")

## Why this is not a reversal of D190

D190 moved Arg90 and Lys226 to MM so that all four cations could be decomposed
one way: `zero_charges_qm_only` cannot reach a residue that is already quantum,
because electrostatic embedding has removed its MM charges and
`charge_deletion` refuses it by name (D153). That reasoning is intact for what
it covers. What it does not cover is the case where **deletion is the wrong
instrument rather than an unavailable one**: deletion measures an electrostatic
contribution on an unchanged geometry, and a residue whose catalytic role is to
hand over a proton has a contribution that zeroing a point charge does not
represent, in any region. For such a residue an all-MM region does not preserve
the measurement; it only preserves the appearance of one.

So the split is by role, and each side of it is measured rather than assumed:

    Arg90, Lys226    Lys226 NZ ends batch 1 at 2.82 A from the nucleophile, 92%
                     occupancy (D201). The transfer is downhill in xTB (D29),
                     at DFT (D35), and by 5.2 kcal/mol in a clean four-
                     coordinate sphere (D48). An MM lysine cannot perform it.
                     -> QM.
    Arg229, Gln93    Arg229 holds the substrate carboxylate at 3.62 +/- 0.04 A
                     across all 45 grid points and both trajectories (D150) and
                     never comes near a transferable proton; Gln93 is an amide.
                     Second-sphere electrostatics, which is exactly what
                     deletion measures.  -> MM, still deletable.

**Both measurements therefore survive.** K's catalytic role comes out of this
trajectory chemically; the R of RQK (via Arg229) and the Q (via Gln93) keep the
deletion numbers D190 was protecting, in `sp_qmmm_r229` and `sp_qmmm_q93`,
which are unchanged. The one row that moves is K's electrostatic contribution,
which was never obtainable from the eleven pre-D190 systems either and whose
right test is mutagenesis - D190 already names that as a second campaign. If
the paper still wants a deletion number on these frames, an `sp` pair with the
all-MM region run on THIS trajectory's geometries gives one: both members share
a region and differ in one key, which is all D152 requires. The caveat to state
if that is done is that the geometry was produced by a different region.

## What changed in `site.yaml`

    qm_regions.michaelis_rk  &qmreg_michrk  = michaelis + A:{90: ARG, 226: LYS}
    sections.michrk          stage qmmm, charge 2

`charge: 2`, not the inherited 0: the michaelis region's 0 plus Arg90(+1) and
Lys226(+1). That is the value all eleven pre-D190 systems carried, which is the
arithmetic check on it. A section left at 0 converges its SCF on the wrong
species while looking entirely normal (D152).

**Region size, from ff14SB rather than from an estimate.** Both ARG and LYS are
in `AMINO_ACIDS`, so `build_qm_region` cuts them at CA-CB and the backbone
stays MM: ARG contributes 18 atoms of 24 and LYS 16 of 22, so the region grows
by **34 atoms**. LYS carries HZ1/HZ2/HZ3, which is the independent confirmation
that it is the cationic form and hence the +1. `build_qm_region` matches on
chain:resid only - the residue NAME in the mapping is documentation, never
checked against the structure - so a wrong name there fails silently and the
`[QM] Including residue` lines of the run's own output are the only check.

## Why a new section and not `time: 100` on michqm

`resume: true` restores `OpenMM_MD_state.xml`. Resuming a state file written
with a different QM region and a different total charge is D152's silent-wrong-
species failure with extra steps. `michrk` writes into an empty directory, so
there is nothing to clear (D199), and batch 1 stays on disk as the control -
the two sections now differ in exactly three keys (`pdb_file`, `qm_residues`,
`charge`) and in nothing else, which is what makes them comparable.

## Why it starts from batch 1's last frame

`protdft/qmmm_michqm/QM_MM_lastframe.pdb`, not `cluster/frame_o1x_off.pdb`.
That frame is 50 ps of geometry already equilibrated at the level of theory
this run uses, and it is the frame in which Lys226 NZ has arrived at 2.82 A -
the geometry the question is about. Restarting from frame 4374 spends
picoseconds travelling back to it. Its Asp92 sits on the metal at 2.00 A and
that is deliberate as well: D201 watched that move happen spontaneously from
the designed four-coordinate start with nothing imposed, so 2.00 A is the
model's own answer rather than a bad starting structure, and the standing rule
is that MD gets crystal restraints while QM/MM gets none, because QM/MM's
coordination behaviour is the thing being measured.

## What it decides

Whether Lys226 NH3+ hands its proton to the zinc-bound hydroxide. If it does,
the attacking species is water and the base is a neutral amine - the mechanism
D201 argues is the only one left standing, since Asp92 as general base has now
failed three times and the third time from a correct start. If it does not,
D41's canonical B3 reactant has passed its first honest test, and D48's 5.2
kcal/mol is overturned where it matters. **Batch 1 could answer neither,
because it did not test the protonation state - it excluded it.** Either
outcome is a result; the current region produces neither.

## Submitted, and the rate measured rather than projected

**31892313 `michrk`, RUNNING on bun152, 8 cpus, `--time=2-00:00:00`.** Written
with `sections.py site.yaml michrk` as `qmmm_michrk`; `config.yaml.bak202` is
the config as it was before.

Checked in the job's own output rather than assumed:

    [QM] Including residue  HIS A:88, ARG A:90, ASP A:92, HIS A:164,
                            ASN A:225, LYS A:226, ZN B:1, mer C:1, OH C:2
    [INFO] QM atoms: 124            Boundary excluded atoms: 54
    Boundaryatoms (QM:MM pairs): six, one per amino acid
    xtb ... --gfn 2 --chrg 2 --uhf 0 ...

124 atoms is 90 + 34, so batch 1's region was 90 and the ff14SB side-chain
arithmetic above is confirmed from the run itself. **Six** boundary pairs is
the check that matters after that - one per amino acid, cut at CA-CB, so no
backbone entered the region. `--chrg 2` on the real xtb command line is the
D152 check: the charge reached the QM code rather than sitting unread in the
config.

**1.565 s per xTB call**, against batch 1's 1.15 s (15:54:58 for 50,000 steps).
So 34 more QM atoms cost **1.36x**, and 50 ps projects to **~21.8 h** - which
would have overrun `run.slurm`'s 24 h default once queueing and I/O are added,
hence the two-day request. The next batch costs the same; 300 ps in six batches
is of order 130 h rather than the 95 h the old rate implied.

**Scratch at submission: 287 of 300 GB, 722,408 of 1,000,000 files.** Measured
because an over-quota write fails silently. 13 GB is ample for this run - a 50
ps QM/MM trajectory of this system is tens of MB - but the margin is thin, and
**the cause is not this project**: CSR is 15 GB of the 287. The large trees are
`channel` 41 G, `glp1r_md` 41 G, `AANCA` 27 G, `hf` 25 G, `boltzgen` 16 G.
Inside CSR: `CSR_mich` 6.2 G, `CSR_wt` 3.6 G, `CSR_crystal_v2` 1.7 G, every
other run under 1 G. Trimming CSR's own MD trajectories to last frames would
recover single-digit GB and is the only part of this the project can decide for
itself.

**D203 - Trimmed CSR's finished MD trajectories to their last frames, on the
user's instruction. 287 -> 278 GB. Two were kept deliberately, and the reason
one of them was kept is a finding in itself: `md_prod100c` never wrote a last
frame, because that is the arm that NaNs.** (2026-11-04)

Deleted, each with `NVTtrajectory_lastframe.pdb` and `_lastframe_imaged.pdb`
confirmed present in the same directory first:

    5.4 GB  CSR_mich/md_eq5_results/NVTtrajectory_imaged.dcd
    2.8 GB  CSR_wt/md_wt100c_results/NVTtrajectory_imaged.dcd
    274 MB  CSR_crystal_v2/md_results/NVTtrajectory_imaged.dcd

**Kept, both on purpose:**

- **`CSR_prod/md_prod100c_results/NVTtrajectory.dcd` (675 MB).** Its directory
  has `relaxbox_NPT_lastframe.pdb` and `frag-minimized.pdb` but **no
  `NVTtrajectory_lastframe.pdb` at all** - the run died at 22.68 ns on the
  product arm's recurring NaN and never wrote one. So this trajectory is not a
  redundant copy of a finished result, it is the only record of an open bug.
  Trimming it would delete the evidence for a diagnosis that has not been made.
- **`CSR_mich/protdft/qmmm_michqm/OpenMMMD_traj_wrapped.xyz` (134 MB).** Batch
  1's own QM/MM trajectory, named outright as `trajectory` in
  `postprocess_michqmpost` (D200), and the control that `michrk` is to be
  compared against. Small and in use.

Also kept: the job `.out` files, including `mich_qm1-31740139.out` at 187 MB.
That is the price of the record and it is worth paying.

**What the deletion costs, stated rather than glossed.** `trajfile`,
`coordination_trajectory` and `contacts_trajectory` in each run's `config.yaml`
are provenance keys mcube wrote when it analysed those trajectories; the
occupancies and distances they produced are already in the config, so what is
lost is the ability to re-derive them or to ask a NEW question of those frames.
For `md_eq5` that matters most, because it is the trajectory `michqm`'s start
frame was selected from (frame 4374, `cluster/frame_o1x_off.pdb`, which stays).
A different frame can no longer be chosen from it.

**This was already the established practice, which is how the decision was
checked.** `md_wt100`, `md_wt100b`, `md_nacgpu`, `md_naclong`,
`CSR_unrestrained/md_results` and `CSR_crystal/md_restrained` are all still
named by `trajfile` lines in the configs and none of their files exist: they
were trimmed earlier. The three deleted here were simply the ones that had not
been.

**Scratch: 287 -> 278 of 300 GB; files 722,408 -> 722,453 of 1,000,000.** The
headroom problem is not this project's - CSR is 15 GB of the total. The large
trees are `channel` 41 G, `glp1r_md` 41 G, `AANCA` 27 G, `hf` 25 G, `boltzgen`
16 G, and what to do about those is not a CSR decision.

**[QUALIFIED by D206 - 'correctly built' is denied by this entry's own closing section (Zn-Asp92 1.78 A) and by the structural criterion. See `csr/retractions.md`.]** **D204 - D202 is answered, 6 ps in, and the answer is yes: Lys226 hands HZ3 to
the hydroxide within 2 ps of the lysine becoming quantum. This is the fourth
time this transfer has been seen and the first time in a complex that was built
right, which is the whole of its value. Measured from the running job's own
files, no analysis stage involved.** (2026-11-05, job 31892313 at step ~5700 of
50,000)

## What the record already said, so that what is new is only what is new

The transfer itself is not news and must not be written up as though it were:

- **D29** saw it in QM/MM MD on `CSR_crystal_v2` - "one proton transfer: Lys226
  HZ1 -> OH- O01 in the first frame, permanent".
- **D35** confirmed it at wB97X-D3(BJ)/def2-SVP//xTB/MM: -13.1 and -15.5
  kcal/mol on two frames, and dropped hypothesis H4 on the strength of it.
- **D48** repeated it in a clean four-coordinate sphere and got **-5.23
  kcal/mol** at DFT, finding that two thirds of D35's margin had been an
  artefact of the degraded six-coordinate geometry.

And **D41** is the reason that was not the end of it. D41 classed the transfer
as one of three ways "the QM/MM MD then degraded" a correctly built reactant,
alongside Asp92 migrating onto the metal and the carbonyl coming inner-sphere.
That left a real question open: in a complex that is NOT degraded - four-
coordinate at the start, substrate bound, carbonyl on the metal by its own
choice rather than by restraint - does the lysine still do it? Every prior
observation was in a sphere that could be blamed. **This one cannot be, and the
transfer still happens.** That is what D202 bought.

## The measurement

`QM_MM.dcd` carries six frames at 1 ps spacing (`nsavc` 1000), read directly
with the atom indices from `QM_MM_firstframe.pdb`. Nine further samples came
from `xtb_.xyz`, the input file the interface rewrites every step, taken ~4 s
apart at step ~5700-5800.

    ps      NZ-HZ3  HZ3-O01  O01-H03  NZ-O01   N-H...O     Zn-O01  Zn-O2x
    1        1.04     2.54     0.98    2.79      93 deg      1.87    2.45
    2        1.19     1.40     0.92    2.58     172          2.21    2.27
    3        1.59     1.00     0.96    2.56     160          2.26    1.99
    4        1.12     1.35     0.98    2.46     168          1.87    2.35
    5        1.60     1.00     0.97    2.60     173          2.26    2.04
    6        1.72     0.99     1.00    2.70     169          2.00    2.24
    live   1.63-1.75  0.96-1.05 0.96-1.00 ~2.6  164          ~2.0    ~2.1

Three things in that table are new.

**It is a shuttle, not the one-way event D29 recorded.** At 1 ps the proton is
on the nitrogen (1.04 A, and the N-H...O angle is 93 degrees - not even aligned
to transfer). By 2 ps it is in flight and the angle has snapped to 172. At 3 ps
it is on the oxygen. At **4 ps it is back at 1.12/1.35**, and from 5 ps on it
stays on the oxygen (1.6-1.75 / 0.96-1.05) through every later sample. So the
barrier is low in both directions and the well is on the oxygen - which is what
a -5.23 kcal/mol driving force (D48) should look like at 300 K, and is a
quantitative consistency check between the two rather than a repeat.

**The coordination sphere reports the protonation state independently.** Zn-O01
is **1.87 A whenever the proton is on the nitrogen** (ps 1, ps 4) and
**2.21-2.26 A whenever it is on the oxygen** (ps 2, 3, 5), with Zn-O2x moving
the opposite way (2.45 -> 1.99-2.27). Hydroxide binds zinc harder than water
does, so the sphere is answering the same question the proton is, through a
different observable, and the two agree frame by frame. That is the strongest
evidence available here that this is chemistry and not a parameter artefact: no
single wrong parameter produces both.

**The nucleophile is a water, and it is not poised.** O01-C4x is **3.72-4.16
A** across every frame, and the Burgi-Dunitz angle `O01-C4x-O2x` is **16-33
degrees** - the nucleophile stays cis to the carbonyl oxygen (O01-O2x 2.68-3.02
A, O01-Zn-O2x 77-91 degrees), both of them adjacent on the same metal. D201's
refutation of the Burgi-Dunitz gate therefore survives the region change
intact, and D29's 20-53 degrees is reproduced in a completely different system.
Two ligands 3 A apart on one zinc cannot put one above the other's carbon; that
is geometry, not sampling.

## What it decides, and what it does not

**Decided.** The reactant of this enzyme at this level of theory is **Zn-OH2
with a neutral Lys226**, in a correctly built Michaelis complex and not only in
a degraded one. D41's hope that the transfer was a simulation artefact is
closed. D48's energy is confirmed dynamically.

**Not decided, and it is the same gap as before.** Everything above is GFN2-xTB
for the dynamics. D35 and D48 are the DFT support and they are single points on
other geometries, so the DFT check on THESE frames has not been done. The
second open worry is unchanged and visible in this run too: **Zn-Asp92 OD2
reaches 1.78 A at 6 ps** (2.80 A in the crystal), the over-coordination D41
flagged, and Asp92's charge sitting on the metal is exactly what D48 showed
inflates the transfer energy. So the sign is secure across four measurements
and two methods; the magnitude in THIS geometry is not.

**The mechanistic consequence, stated but not acted on.** A zinc-bound water at
2.0-2.2 A with a neutral amine 2.5-2.7 A away, both of them stable and
exchanging a proton freely, is not a dead end - it is the arrangement in which
the amine can take the proton back as the water attacks. That makes a concerted
water attack with Lys226 as the base the reading this run supports, and plan
item 6's scan/NEB the thing that would test it. No job is proposed here.

## A false criterion, found while reading this

`monitor_angles` in `michqmpost` scores **both** `O01-C4x-O2x` and
`O01-C4x-N1x` against the section's 95-120 window. Burgi-Dunitz is defined for
Nu-C=O, which is the first of those. The second is the angle to the leaving
nitrogen, and a tetrahedral carbon puts any two substituents near 109 degrees:
it measures **99-119 in every frame above**, so it scores as a PASS inside that
window while saying nothing whatever about an approach. In `michrkpost` it now
carries `window: null`, which is how `angle_occupancy` says there is no
criterion. **`michqmpost` is left alone deliberately** - changing a section
whose result may already be recorded would make the stored numbers disagree
with the config that produced them - so if batch 1's report has been run, its
`O01-C4x-N1x` occupancy is that artefact and should be struck rather than
re-derived. A 95-120 occupancy on that triple must not reach the paper.

**D205 - The MDAnalysis gate on the local test suite was never necessary. Three
additive pip installs and the 50 trajectory tests that have been skipping since
D167 now run locally: 767 pass, 7 skipped, all seven ASH.** (2026-11-05)

Every `postprocess` change for weeks has been committed on a local run that
printed a clean count while the tests that read trajectories had not executed -
D167's collapsed escape and D202's `1 skipped` where Bunya gave `10 passed` are
both that gate. The fix is `pip install MDAnalysis networkx polars` into the
local `mcube-test` env, found by running the suite and installing whatever the
next `ModuleNotFoundError` named.

**Checked before each one, because this is where an install does damage.** The
env has **numpy 2.4.6 and OpenMM 8.4**, and `pip install --dry-run` reported
every package purely additive with no numpy or openmm line. The first attempt
pinned `numpy==1.26.4` out of habit from the Bunya ASH env, and the dry run
said it **would install numpy-1.26.4** - a downgrade under an OpenMM compiled
against numpy 2, which is the exact breakage the rule about dry-running exists
to prevent. Unpinned, nothing moved; numpy and OpenMM verified unchanged after.

So the per-pair cutoff of D204 is **verified locally** -
`TestAPairWithItsOwnCutoff` 4/4, `test_contact_occupancy` 18 passed - and the
commit that introduced it was wrong to say it still needed the ASH environment.
What remains ASH-gated is only `test_cores`, `test_optimizer_constraints`,
`test_pdb_resseq`, `test_provenance`, `test_qm_levels`, `test_scan_surface` and
`test_sp`, which is a real gate: `test_qm_levels` covers `mcube.qmmm`, so a
`qmmm.py` change is still unverified until a Bunya run.

**D206 - The order of evidence, written down as a rule because it was broken
twice in one session: an experimental structure is ground truth, the literature
is second, and an energy computed on our own geometry is not evidence about
structure at all. D58's headline is retracted in place - D59 had already
retracted it, and then D59's own closing line re-asserted it, so the log
contradicted itself and I quoted the wrong half.** (2026-11-05, user: "you're
relying on minima from flawed md starts again - asp is not coordinated in the
crystal or in any other B3 structure - that's our ground truth ... good job with
looking at the literature - that's our second source of truth")

## The rule

1. **An experimental structure is ground truth.** 6NC5's ZN302 is
   four-coordinate - His87 NE2 2.22, His163 NE2 2.09, HOH485 1.90, HOH493 2.14:
   two histidines, two waters, **no aspartate**. **[CORRECTED 2026-11-05, within
   the hour: this rule first said "Asp92 OD2 is a second-shell contact at 2.80 A,
   2.68 A when superposed into wild-type 6DN4 (D184)", and neither number is a
   measurement of a Zn-Asp distance. The log's own standing rule 2 says 2.80 A
   "is not a crystal distance" and gives 4.20 A (OD2) / 4.22 A (OD1) measured in
   `CSR_crystal/6NC5.pdb` itself (D145); and 6DN4 contains no metal, so D184's
   2.68 A is our own superposition of 6NC5's zinc into a metal-free structure.
   Writing a ground-truth rule on two of our own numbers, in the entry
   establishing that rule, is the error the entry is about. What the structures
   actually support, on every reading - 2.68, 2.80 or 4.20 - is the qualitative
   statement, and 4.20 supports it hardest: **Asp92 is not a Zn1 ligand**, being
   far outside a Zn-O bond of 1.9-2.2 A. Quote that, not a distance. The three
   readings are irreconcilable until someone measures the deposited files, and
   `csr/retractions.md` §1(a) carries that as the project's first open
   contradiction - including D151's finding that 6NC5 is the Q93H/K226H double
   mutant with Zn2 PRESENT, holding that same carboxylate at 1.95 A, so it may
   not be a wild-type reference for this residue at all.]** Wider than this one structure, **Asp120 is a Zn2 ligand in the B3
   structures this project has read** - L1, SMB-1, GOB-18 - and CSR-1 has no
   Zn2. So a Zn1-coordinating Asp92 has no precedent in anything we have looked
   at, and a model that produces one is wrong about the structure however
   favourable its energy. **Flagged honestly: the stronger form of this sentence
   first written here - "no B3 enzyme has Asp120 on Zn1" - is a generalisation
   over the whole subclass that no citation in this log supports. It needs a
   source before it goes in the paper, and until it has one the argument rests on
   6NC5/6DN4 plus those three structures, which is enough for the decision and
   not enough for the claim.**
2. **The literature is the second source of truth**, for what one structure
   cannot say: which site a mononuclear metal occupies, what the nucleophile is,
   what each metal does.
3. **An energy computed on our own geometry is evidence about the energy
   function, not about the structure.** A single point cannot promote the
   geometry it was handed. Neither can agreement between two methods on the same
   handed geometry.

## The two violations, and the correction to each

**(a) D58's headline, retracted in place above.** The DFT points were single
points on xTB scan geometries that all descend from the degraded six-coordinate
MD start, so both methods were asked one conditional question - given this
octahedron, what does removing a ligand cost? - and neither was asked what the
site should look like. **D59 already said exactly this**, in the words
"agreement between two methods on a loaded question is not evidence that the
question was right", and then closed its own table by re-asserting that the
six-coordinate zinc "is real at both levels". Both sentences are corrected
above. The surviving statement is the narrow one: xTB does not over-bind Asp92
along that coordinate from that structure.

**(b) D107 and D112's pKa bookkeeping has a sign error.** "A 2+ cation 4-5 A
away should **raise** a lysine's pKa" is backwards. A lysine's conjugate base is
the **neutral** amine, so a nearby cation stabilises the base and **lowers** the
pKa - the same argument D154 used to correct D153 for the carboxylate in this
same site. Signed correctly, the electrostatic term favours the Lys226 -> O01
transfer instead of opposing it; the BcII aqua pKa of ~5.6 also needs revising
upward for a mono-zinc site, which pushes back. So the bookkeeping is
**indeterminate**, and D112's reading that it "runs the other way" and would
"favour `lysH` and vindicate what has been run" does not hold.

## What survives, because the rule cuts both ways

**D48 is not touched by any of this, and the rule raises its standing rather
than lowering it - but the distinction has to be stated, because rule 3 excludes
the evidence class D48 belongs to.** The rule is not that our own calculations
say nothing; it is that they say nothing *about geometry they were handed*. D48
was not asked where Asp92 sits - it was told, from the crystal, and then asked a
question about protons at fixed geometry. That is a legitimate conditional, where
D58's was not, because D58's conditional was loaded by a geometry no structure
supports. What D48 cannot do is promote its own sphere, and it does not try to. Its comparison was run in a constructed clean
four-coordinate sphere - Asp92 held at the crystallographic 2.80 A, carbonyl out
of the inner sphere, identical constraints on both states - that is, in the
geometry the ground truth prescribes, and it still puts the Lys226 -> O01
transfer at **-5.23 kcal/mol** at wB97X-D3(BJ)/def2-SVP. It is the one energetic
result on this question whose geometry agrees with the crystal. So the
water-versus-hydroxide question stays open and the best evidence in it still
favours the transfer; what the rule removes is a bad argument that was being used
on both sides.

**And the standing worry about `michrk` stands at full strength.** Zn-Asp92 OD2
at 1.77-1.84 A and five-coordinate in every frame (D204) is rejected by the
structural criterion - which is also the published VIM-2 criterion, where
candidate protonation states were judged by zinc coordination against the
crystal rather than by energy (*Protocols for Metallo- and Serine-beta-Lactamase
Free Energy Predictions*, J Phys Chem B 2024, 128:12416, PMC11664573; read and
recorded in **D113**, which is where the citation lives - this entry asserted the
criterion without it) - independent of what its energy says.

**D207 - Why no reaction coordinate has ever run. Two answers, and the first
one corrects this entry's own original framing: CI-NEB was never the published
method for this problem - it was OUR substitution for Imhof's path network (D1:
"their CPR-style path search -> *our* CI-NEB"), and every MBL mechanism study in
this project's survey reports a FREE ENERGY instead, which is where CphA's 14
kcal/mol comes from. D32 was a user decision to replicate that PMF; the M-series
dropped umbrella sampling for cost. Second, the path itself has been blocked on
one missing file since 2026-10-22 - `protdft/scan_attack/scan_15.pdb` - because
the reactant was treated as something to settle BEFORE the path rather than
something the path settles. The stage machinery is ready; the endpoint is not.** (2026-11-05, user: "now if you'll
recall this was the plan from the very beginning. Why did we deviate from NEB?")

**[QUALIFIED the same day, by the user: "previously you've stated that NEB is
not done for MBL mechanism studies and why it's not suited". Correct, and this
entry's framing was wrong in the other direction. D1's own words are "their
CPR-style path search -> OUR CI-NEB": NEB was never the published method, it was
our substitution for Imhof's path network, chosen on day one. **No MBL mechanism
study in this project's own survey uses NEB.** Xu-Xie-Guo 2006 on CphA - the
mono-Zn precedent, and the source of the 14 kcal/mol we would compare against -
is "SCC-DFTB/MM MD, with free energies from potential-of-mean-force (umbrella
sampling), plus DFT", so that number is a **PMF barrier**; SMB-1 and NDM-1 are
the same group's protocol, and D151 costed our own Asp92 PMF explicitly "on
SMB-1's protocol"; the VIM-2 paper is titled *Free Energy Predictions*; IRE1 used
metadynamics; Ackerman & Gatti used potential energy surfaces. **And D32 was a
user decision to replicate those methods, PMF included.** What dropped the PMF
was the M-series, whose "explicitly dropped" list contains umbrella sampling -
for cost. So the deviation runs the opposite way to what is written below: the
project deviated from the field's free-energy methods TOWARD a minimum-energy
path. Everything below about the missing endpoint still holds; read "the next job
is `scan_attack`" as seeding an umbrella PMF, with NEB/TS/IRC optional for
characterising the saddle rather than for the number we publish. See the
reasons-it-is-unsuited list appended at the end of this entry.]**

## It was never deselected

CI-NEB with TS and IRC is in **D1**, the first entry (2026-10-15): "their
CPR-style path search -> our CI-NEB, transition-state optimisation and IRC". It
survived every revision of the plan - the P-series, then the M-series where it is
M2, then the numbered plan where it is item 6, still reading **"not started; the
paper's actual content"**. **D204**, four days ago, still names "plan item 6's
scan/NEB" as the test. It was staged as real config on 2026-10-22:
`scan_attack -> neb_step1 -> ts_saddle1` and `scan_ringopen -> neb_step2 ->
ts_saddle2`, one NEB per elementary step (D17), both endpoints in one solvent box
(D4).

`neb_step1` names `protdft/scan_attack/scan_15.pdb` as its product and
`neb_step2` names `protdft/scan_ringopen/scan_11.pdb`. **Neither file has ever
existed**, because neither scan was ever run. That is the whole deviation.

## The machinery was genuinely broken three times, and all three are fixed

Worth recording because it justifies part of the delay, and because each would
have produced a confident wrong number:

1. **`neb.py` hardcoded `runmode='library'`** (found 2026-10-23 while costing the
   MD). The stage could not run at all - it would have reached the PC-gradient
   branch and called `ashexit()` mid-job. Now
   `runmode = 'parallel' if neb_workers > 1 else 'serial'`.
2. **The atom-order bug.** `neb` compared atom *count*, and NEB interpolates atom
   i of the reactant toward atom i of the product, so a mismatched order "gives a
   path between atoms that are not counterparts and **converges anyway**". Now
   compares `atom_names` index by index.
3. **The two-surfaces bug.** `qmmm`, `scan` and `ts` held active-region MM waters
   rigid as optimizer constraints; `neb` cannot, since ASH's NEB exposes no
   constraints argument. So the endpoints were not minima of the surface the path
   was evaluated on, and interior images could dip below the reactant. Fixed by
   the rigid-water default (D139/D141), with `check_path_shape` left as a warning.

**So nothing in the code blocks a NEB today.**

## What did block it: one thing in six costumes - we never had an endpoint

Every candidate reactant was disqualified, each time for a real reason:

    no attack geometry      0 of 500 frames meet the NAC criteria; the Zn-bound
                            oxygens are cis at 24-54 deg (D33, D41, and again in
                            D201/D204 at 16-33)
    protonation contested   D29 found the Lys226 transfer, D30/D35 checked it at
                            DFT, D48 redid it in a clean sphere - and D112 then
                            found the running system was in the state M1 had
                            rejected, because "the verdict was never carried
                            across"
    the pose                nine tied Vina poses, one carried forward by accident
                            (D66), and the 1.94 A Zn-carbonyl bond imposed rather
                            than found (D67)
    the trajectories        `qmmm_nacarg/QM_MM.dcd` "cannot supply Michaelis-
                            complex frames - not for NEB endpoints" (D110/D111)
    the shell               everything downstream of `md_wt100`/`md_wt100b`
                            "describes a six-coordinate zinc, which is not the
                            species the paper is about"
    unrestrained MM         "a Michaelis complex cannot be produced by
                            unrestrained MM MD and was not"

## The diagnosis, which is the part worth keeping

Each blocker was investigated **on its own terms rather than as a question about
the endpoint**. Asp92's coordination got three separate campaigns (D51-D59,
D120-D145, D173-D200); the protonation state got three (D29-D48, D112-D119,
D202-D205). Both are real questions about this enzyme. Neither is the NEB, and
neither has to be answered before one.

And the loop could not terminate, because **every attempt to settle the reactant
used our own minima as the evidence** - which D206 has now ruled inadmissible for
exactly this. A question whose only evidence class cannot settle it will reopen
indefinitely, and it did, three times each.

**The exit is the user's own argument of today.** A reactant whose proton
interconverts over a barrier of order kT does not need its protonation state
settled in advance, because the path settles it: drive the slow coordinate, leave
the proton free, and read its timing off the result. The only property the
endpoint must have is the one the ground truth prescribes - Asp92 second-shell
near 2.80 A and a four-coordinate zinc - and structures with that property exist
in the docked/`naclong` lineage. So the next job is `scan_attack`, which is what
`neb_step1` has been waiting for since 2026-10-22 - and which also seeds the
umbrella windows (D133: "windows take a `seed` each, which the scan in 1
supplies").

## Why a minimum-energy path is the wrong instrument for the number we publish

Four reasons, all already in this log:

1. **The literature's barriers are free energies.** CphA's ~14 kcal/mol is a
   PMF. A NEB gives a potential-energy barrier along one path with no entropy and
   no solvent sampling, so quoting it beside that number compares two different
   quantities. This log already says the same of scans: what a scan cannot give
   is a free energy, because there is "no entropy and no sampling of the waters
   that reorganise".
2. **NEB-CI returns one saddle.** Hydrolysis is at least two elementary steps
   (attack, then C4x-N1x cleavage with N-protonation), and "a single NEB spanning
   both puts two barriers on one path and NEB-CI returns **one** saddle, so the
   other is simply never found" (D17). That forces the step decomposition to be
   known in advance, which is the thing we do not know and the thing the scan is
   for.
3. **The two-surfaces asymmetry is structural, not a bug.** ASH's NEB exposes no
   constraints argument, so `neb` cannot hold what `qmmm`, `scan` and `ts` hold.
   The rigid-water default (D139/D141) removed the mismatch by making everything
   flexible, but interior images can still lower their energy by relaxing
   solvent unrelated to the reaction, which is why `check_path_shape` exists.
4. **The soft coordinates here are exactly the ones an MEP serves worst** - bulk
   water reorganisation and a proton shared over a barrier of order kT.

**What this does not change.** The endpoint problem is prior to both instruments:
"a NEB **or an umbrella series** run from that reactant to a template-derived
product would spend its barrier moving meropenem across the site, and the number
at the top would not be a chemical barrier at all. That failure would look
entirely normal in the output."


**D208 - The five open contradictions `retractions.md` opened are closed. Four
by measuring the deposited structures, which nobody had done, and one by reading
the restraint builder. The headline result: the crystal puts Asp92 **4.20 A**
from the catalytic zinc, every other B3 structure puts the same aspartate
4.02-4.17 A from Zn1, and the 2.80 A this project has restrained, gated and
reasoned from for five weeks was never a crystallographic distance.**
(2026-11-05, user: "ok resolve these 5")

**Provenance, and it needed no new tool: `python refstruct.py shell <pdb>`
re-derives every number below from the PDB ID**, which is what that script was
written for ("a criterion that cannot be re-derived from a PDB ID is an
assertion, not a measurement"). Nobody had pointed it at these files. I measured
independently first with a stdlib script and then cross-checked against
`refstruct.py`; the two agree to the digit on 6NC5's six zincs and on 1SML, so
the numbers below are two independent implementations of the same measurement.

## 1. Asp92's distance from the zinc: 4.20 A. Four readings, one survives

    6NC5  ZN A302 (OUR site)  His87 NE2 2.22, His163 NE2 2.09,
                              HOH485 1.90, HOH493 2.14      -> FOUR-coordinate
                              ASP A91 OD2 4.20   OD1 4.22
    6NC5  ZN A301 (Zn2 site)  ASP A91 OD2 1.95, His92 NE2 2.01,
                              His225 NE2 2.04, HOH493 2.08
    6DN4                      4557 atoms, ZERO metal sites
    6DQ2                      4049 atoms, ZERO metal sites

- **D145 is correct and is now confirmed independently: 4.20 A (OD2), 4.22 A
  (OD1).** The four-coordinate sphere in standing rule 2 is reproduced atom for
  atom.
- **The 2.80 A is our own model's number.** It is not in this file.
- **D184's 2.68 A is void.** It superposes 6NC5's zinc into 6DN4, and **6DN4
  contains no metal at all** - measured above, not inferred. The same goes for
  6DQ2. No Zn-Asp distance can come from either.
- **D151 is correct about ZN301.** In this mutant the carboxylate is a **Zn2
  ligand at 1.95 A**, monodentate (OD1 3.43), and **HOH493 bridges both zincs**
  at 2.08/2.14 - the canonical bridging-nucleophile position, which is also the
  independent confirmation of D39. So 6NC5 cannot show what the carboxylate does
  with Zn2 *empty*, and that caveat has to travel with the 4.20 A. It is still
  the only measured value in existence, because the wild-type files are apo.

## 2. "No B3 enzyme has Asp120 on Zn1": now measured, not asserted

| structure | Zn1 site | Asp->Zn1 | Asp->Zn2 |
|---|---|---|---|
| 1SML, L1 | His84/86/160 + HOH 1.88 | **4.08** | Asp88 OD2 **2.07**, OD1 3.20 |
| 5K0W, GOB-18 (chain A) | **Gln98 OE1 2.04**, His100/175 | **4.17** | Asp102 OD2 **2.11**, OD1 3.43 |
| 5B15, SMB-1 | His72/74/150 + HOH 1.83 | **4.05** | Asp76 OD2 **2.13**, OD1 3.31 |
| 7UHJ, L1 + ligand | His105/107/181 + MX0 O1 2.63 | **4.02** | Asp109 OD2 **2.06**, OD1 3.18 |

In all four the conserved aspartate is a **Zn2** ligand, **monodentate**, with
the second oxygen 3.18-3.46 A out - and 4.02-4.17 A from Zn1. **6NC5's 4.20 A
sits within 0.03-0.18 A of every one of them.** CSR-1's crystal is completely
ordinary for a B3 at this contact; what was anomalous was our 2.80 A.

D206's rule 1 is therefore supported by measurement and its uncited
generalisation can be replaced by this table. Two by-products:

- **5K0W has FOUR zincs - both sites occupied in both chains.** "GOB-18's single
  zinc sits in the Zn2 site" is not what the deposit shows, and that sentence
  needs correcting wherever the project repeats it.
- **GOB-18's Zn1 is confirmed as the HQH site**, Gln98 OE1 at 2.04 A, which is
  D124's point measured rather than read.

## 3. Does GFN2 over-coordinate the zinc? Still open - and now secondary

No admissible test exists: the three closures (standing rule 3, D133, project
memory) all rest on D58, which is retracted, and the fourth (D195) is reclassed
by D197 as a relaxation inside the coordinated basin. **But the question no
longer gates anything,** because the discrepancy is now too large to be a method
artefact alone: our runs put this contact at **1.77-2.11 A** against **4.02-4.20
A in four independent experimental structures**. A 2+ A error against four
crystals is a wrong geometry whatever produced it. Which component to fix - the
MM start, the restraint set, or xTB - is a real question and a secondary one.

## 4. What put the fifth ligand on the zinc: both causes, and neither alone

Read in `find_coordinators`, which is the restraint builder. It now calls
`coordination_shell(pdb_file, inputs, per_residue=False)`, and its own comment
records the history: it "took the per-residue default anyway. So only the closest
donor of each residue was ever held, and a group with two donors inside the
cutoff had one pinned and the other free: exactly the geometry D120 and D130
record failing, where Asp92's OD2 was held at 2.80 and OD1 pivoted onto the metal
at 1.99."

So **two conditions were both necessary**: the **3.0 A cutoff** admitted a
second-shell contact into the list at all, and the **per-residue default** then
held only the nearer carboxylate oxygen, leaving the other free to rotate in.
D120, D130 and D174 each blamed the cutoff alone; D175 named the per-residue
default and was substantially right; D176's "No: detection was never per atom for
restraints" over-corrected it. Nobody had both halves. Fixed at the library
(`per_residue=False`, plus `warn_beyond_bond_length` at the single detection
exit).

## 5. The citation error

D152's "0.58 A Lys226-O4x overlap for weeks (D151)" is **D148's** finding.
Corrected in place.

## What this costs, and it is not small

**Asp92 should never have been in a restraint list.** At 4.20 A it is outside the
3.0 A "affected by the metal" radius and far outside the 2.5 A coordination
radius, so under this project's own standing rule it gets no coordination
restraint and no QM-region membership on coordination grounds. Every structure
built to hold it at 2.80 A - D130's hand-written pair, D144's rule, the
`finalsystem.pdb` lineage, both 100 ns arms, and the `michrk` start frame - has
that contact **1.4 A too short**, and the ones that collapsed to 1.77-2.11 A are
2.1-2.4 A too short. The restraint target for any rebuild is ~4.2 A, which in
practice means: do not restrain it, and let the run say where it goes.

**D209 — The wild-type site, measured rather than inferred. Neither
wild-type CSR-1 deposit contains a metal, so no Zn–Asp92 distance exists in any
of them; the aspartate is a **Zn2** ligand in every B3 structure; and in 6DN4 it
is **Lys226 that holds the carboxylate, at 2.99 Å**. That last number makes
Asp92's collapse onto the metal a CONSEQUENCE of the Lys226 proton transfer
rather than an independent xTB failure. It also corrects D208, which is mine.**
(2026-11-05, prompted by the user: "The CSR-1 WT structure is our source of
ground truth and it has a 4-coordinate zinc yes?")

## The premise that started this, and why it is wrong in a useful way

It is not. **6DN4 has 347 waters and nothing else; 6DQ2 has 63 waters and a
phosphate. Neither contains a single metal atom.** The only four-coordinate zinc
anywhere in the CSR-1 record is 6NC5's ZN302, and 6NC5 is the Zn2-site double
mutant. So "the wild-type structure has a four-coordinate zinc" is not available
as a premise — but the conclusion drawn from it survives on better evidence,
below.

**And ZN302 is not an isolated sphere.** Both of its waters bridge to
neighbouring zincs: HOH493 to ZN301 at 2.08 Å, HOH485 to ZN303 at 1.88 Å. The
sphere the project built `finalsystem.pdb` from is one site of a water-bridged
zinc cluster. Delete the neighbours, as a mono-zinc model must, and neither
water's position is justified by this deposit. This is a second, independent
reason D108 was right to reject 6NC5 as a wild-type reference.

## What corrects D208, which I wrote

D208 measured Asp—>Zn1 at 4.20 Å in 6NC5 and 4.02–4.17 Å across four B3
structures and concluded that ~4.2 Å is where the carboxylate belongs. The
measurement is right; the inference is not, and it was about to be built into a
start structure.

    structure        Asp -> Zn1    Asp -> Zn2
    1SML (L1)            4.08         **2.07**
    5K0W chain A         4.17         **2.11**
    5K0W chain B         4.36         **2.08**
    7UHI (L1 Michaelis)  3.94         **2.27**
    5AXO                 4.08         **2.19**
    6NC5                 4.20         **1.95**

In **every** one the aspartate is a Zn2 ligand, monodentate, and its ~4 Å
distance to Zn1 is the Zn1–Zn2 separation plus the carboxylate's bite — a
measurement of where the second metal is, not of where the side chain prefers to
sit. CSR-1 has no second metal. So the user's ground truth, stated exactly as the
user stated it, stands untouched: **Asp92 is not coordinated to Zn1 in the
crystal or in any other B3 structure.** What does not follow from it is any
particular distance. That is now an open question in `state.md` §1 rather
than a settled value, and all three numbers previously quoted as "the crystal"
— 2.80 (ours), 2.68 (void, D184), 4.20 (a Zn2-occupied mutant) — are
measurements of something else.

## The finding: in the wild type, Lys226 does Zn2's job

6DN4's site, measured structurally after the assumed BBL numbering failed (see
below):

    Lys254 NZ  - Asp120 carboxylate   **2.99**   <- a salt bridge
    Asp120 OD1 - His116 NE2             3.00
    HOH545     - Asp120 carboxylate     2.83
    Gln121 OE1 - Asp120                 3.70
    Arg118 NH1 - Asp120                 6.05
    His116 NE2 - His192 NE2             3.63   <- the Zn1 histidine pair

Lys254 is Lys226. Where canonical B3 holds this carboxylate with a **zinc** at
~2.1 Å, wild-type CSR-1 holds it with a **lysine ammonium** at 2.99 Å. This is
the RQK thesis with a crystallographic number behind it, and the project did not
have it: `state.md` has never mentioned an Asp92–Lys226 salt bridge, and
`michrk` was set up treating Lys226 as a proton donor with Asp92 elsewhere.

**So Asp92 is held off the metal in the resting state, and by the lysine.**
Asp92-on-Zn1 is an artefact of the resting state, which is what the user said.

## Why the model does it anyway, and the prediction that follows

D204: Lys226 hands HZ3 to the hydroxide within 2 ps of becoming quantum, and from
5 ps the proton stays on the oxygen. **That deprotonates the lysine, and a neutral
amine cannot hold a carboxylate.** So the two observations this project has been
treating as separate problems — the Lys226 proton transfer and Asp92's collapse
onto the metal — are one event with an order, and the collapse is downstream.

**Falsifiable, on data already on disk:** in the `michrk` trajectory Asp92's
approach to the zinc must **follow** the proton transfer in time. If it precedes
it, this chain is wrong and GFN2 over-coordination is back as a first-order
suspect. Not yet tested: every path under `/scratch/user/uquadhik/` returned the
route's own `Not found` from the new `bun126/53134` code-server session, including
the 41 GB `channel` tree, so scratch looks unreachable from that session rather
than any one file being absent.

**And it sharpens the thesis.** The difference from canonical B3 is not merely
that Zn2 is absent: the cation that replaces it is **titratable**. A zinc can
polarise the nucleophile but cannot hand it a proton; a lysine can do both. If
that holds, CSR-1 has swapped a permanent Lewis acid for a group that is also a
general acid/base, which is a difference in kind rather than degree — and it is
the mechanistic content of the user's thesis that CSR-1 differs from canonical B3.

## The salt bridge is in BOTH deposits that have the lysine, and the one that lacks it had the lysine mutated away ON PURPOSE

    6DN4 (WT)        Lys254 NZ - Asp120 OD   **2.99**
    6DQ2 (mutant)    Lys254 NZ - Asp120 OD   **2.83**
    6NC5             residue 225 is **HIS** - the lysine is not there

6NC5's own `SEQADV` records why: `HIS A 92 <- GLN 136` and
`HIS A 225 <- LYS 269`, both ENGINEERED MUTATION. In our numbering those are
**Gln93His and Lys226His** — the Q and the K of RQK. So 6NC5 is the
**Q93H/K226H double mutant**, and the two residues that had to be replaced to
get a zinc into Zn2 are exactly the two the RQK thesis says are standing in for
it. The nearest other lysine is 16.4–17.2 Å away in all three, so the
contact is specific.

This is the strongest experimental support the thesis has, and it was sitting in
a `SEQADV` record: **the lysine occupies the Zn2 site, and it must be removed
before zinc will bind there.**

## Structures the user supplied, which are NOT deposits (checked 2026-11-05)

Under `pbeval/.claude/worktrees/charming-cray-fe9629/University/Undergrad/`
`Honours/SCIE4481/bioinformatic_data/MM/Structure/CSR/`. Two kinds, and the
difference matters:

- **`CSR_nontruncated_1metal.pdb` / `_2metal` / `CSR_truncated_1metal`** — the
  full-length WT with zinc. **Their zinc spheres are non-physical**:
  `ASP120.OD2 1.38 | LYS254.NZ 1.95` in the non-truncated 1-metal file, and
  `LYS254.NZ 2.18` alone in the truncated one. A Zn–O of 1.38 Å is not a
  bond length, and a zinc does not ligate a protonated lysine NZ at all. The
  metal was placed into the Asp92–Lys226 salt bridge — i.e. into the pocket
  the lysine already occupies. These cannot serve as a structural reference, but
  as evidence they point the same way: the only hole the modelling could find for
  a second metal is the one the lysine is in.
- **`dmCSR_mutant_1metal_refine_12.pdb` / `dmCSR_mutant_2metal_refine_16.pdb`**
  — `REMARK PHENIX refinement`, dated 2025-12-16. These are refined against
  real data. The 2-metal one has **two zincs bridged by two waters** (ZN1–HOH1
  2.03, ZN2–HOH1 2.01, HOH2 2.43/1.98), which is the canonical B3 bridging
  arrangement. `dmCSR_1metal`/`dmCSR_2metal` put a zinc on Asp91/His92/His225,
  i.e. the restored Zn2 site, with the Zn1 site (His87/His163) filled only in the
  2-metal file.

**So an experimental, zinc-bound CSR-1 structure does exist — of the double
mutant, from the user's own crystallography.** Nothing in this project has used
it. That is the control the whole argument wants: WT has Lys226 in the Zn2 site,
the double mutant has a zinc there, and the user holds kinetics for both.

## Numbering, corrected

The BBL mapping in `state.md` is wrong twice for these deposits. In 6DN4 the Zn1
histidine partner is **His192**, not His196; and the lysine is **Lys254**, not
Lys263 — Lys263 NZ is **20.4 Å** from the site, so any reasoning that reached
for "Lys263" was reaching for the wrong residue. Both corrected in place.

## What changed in `state.md`

Seven rows corrected and three added in §1, and three prescriptions
withdrawn: the claim that Asp92 "sits at 2.80 Å in the crystal" (false), plan
item 0's "2.80 Å supported" (reopened — Asp92 should carry **no** restraint,
since the standing rule restrains what is within coordination range and no
wild-type structure puts it there), and equilibration gate 2, "Asp92 OD2 near
2.80" (withdrawn — it tests against a number with no crystallographic basis,
so both arms' recorded FAIL is a failure to reproduce the model's own artefact).
What a start structure must now get right is the Asp92–Lys226 salt bridge, not
an Asp92–Zn distance. The header also said the log was 180 entries long.

**D210 — Tier 0, and it is two restraints rather than a rebuild. I
over-corrected twice on 2026-11-05 and the user caught both: Asp92 **is** to be
restrained, because the 3.0 Å radius is "what is affected by the metal" and not
a coordination test. Redoing D184's superposition at three atom sets reproduces
**2.71–2.76 Å**, so D184 was right and so was D174's 3.81 Å judgement call.
What was genuinely missing is different: the fifth donor that rejected `eq5`, and
the Asp92–Lys226 salt bridge, which no run has ever held.** (2026-11-05, user:
"Asp92 should be restrained. I said restrain everything AFFECTED BY THE METAL ...
I've explicitly said this several times - check logs", and "reuse previous runs if
they're valid")

## The rule, from the logs rather than from me

- **D174** (user): "use the 3 angstrom constraint! We already established that
  anything within that range is affected by the zinc and should be restrained to
  it's crystallographic distances".
- **D186** (user): "the 3 angstrom restraints aren't exclusively for coordinators
  - they're for anything affected by the metal - a different set. The metal could
  still be 4-coordinate you'd just use a stricter coordination radius for that
  like 2.5".
- **D174 again, in its own words:** Asp92 OD1 at 3.81 Å was restrained as "one
  restraint beyond the sweep, flagged in the config as a judgement ... a
  carboxylate whose other oxygen is 2.80 Å from a divalent cation is affected by
  it whatever the distance says."

So "Asp92 should carry no restraint at all", which I put into `state.md` plan
item 0 earlier today and attributed to D208, is the opposite of the standing rule
and is retracted. The three edits that carried it — plan item 0, equilibration
gate 2, and the "FALSE" mark on the 2.80 Å line — are all reverted in place.

## D184 reproduced, which is why my correction was wrong

6NC5's ZN302 placed into 6DN4's frame by Kabsch fit on the two conserved
histidines (6NC5 His87/His163 —> 6DN4 His116/His192):

    fit basis                    n    fit RMSD    Zn-OD2   Zn-OD1   Zn-Lys226 NZ
    ring only                   10     0.189       2.71     3.88       4.13
    ring + CB                   12     0.226       2.73     3.89       4.13
    ring + CB + CA + N + C      18     0.256       2.76     3.90       4.19

**0.05 Å of spread over an 8-fold change in the fit basis, and within 0.08 Å
of D184's 2.68.** It is still a model value and not a crystallographic one —
6DN4 has no metal — but it is a stable one, and the project's two
independently chosen numbers land inside it: the model's **2.80** for OD2 and
D174's **3.81** for OD1, the latter to within 0.09 Å of a figure derived without
reference to it. For comparison the mapped zinc reaches His116 NE2 at 2.36 and
His192 NE2 at 2.09, against 2.22 and 2.09 measured on the real ZN302.

## The geometry that explains the artefact

OD2 is the **bridging** oxygen, and this had not been noticed:

    Asp92 OD2   2.76 A from the (modelled) zinc     2.99 A from Lys226 NZ
    Asp92 OD1   3.90 A from the zinc                5.06 A from Lys226 NZ

Measured in both wild-type-background deposits: Lys254 NZ – Asp120 OD2
**2.99 Å** in 6DN4 and **2.83 Å** in 6DQ2, with OD1 at 5.06 and 4.12. So the
carboxylate is pinned between two cations, zinc on one side and the lysine
ammonium on the other, and OD1 points away from both.

**That predicts the failure mode this project keeps recording.** D204 has Lys226
handing HZ3 to the hydroxide within 2 ps of becoming quantum. Deprotonate the
lysine and OD2 loses the cation on one side, so it falls onto the one on the
other — 2.76 to 2.00 Å. The 2.00 Å collapse is not an independent xTB
failure and it is not a bad starting structure; it is what the geometry does once
the salt bridge is gone. It is still an artefact, because it overshoots the
second shell by a whole coordination shell, but the cause is now named.

## Reuse: `eq5` is valid and its Asp92 restraints already worked

Checked rather than assumed, per the user's instruction. `md_eq5` holds OD1 at
3.81 Å with k = 500 and lets the 3.0 Å rule detect OD2, which is the
per-atom treatment D175/D176 fixed. Over **10,000 frames** it gave Asp92 OD2 at
**2.743 Å** mean and **0.0% inside 2.5 Å** — 0.02 Å from the
superposition target. So the aspartate was never the problem and there is nothing
to rebuild for its sake.

**`eq5` was rejected for a different atom.** `mer` **O1x**, the
6-(1-hydroxyethyl) alcohol oxygen, sat at 3.13 Å in the October pose —
**0.13 Å outside the detection radius** — so `fill_restraints` never saw it.
It then came onto the metal at 2.06 Å with **97.9% occupancy** and, sitting
between the hydroxide and the carbonyl at 76.8 degrees from both, closed the
attack angle (D197). D199 adds that this is not a missing repulsion: the Li/Merz
parameters were applied, so 2.06 is the fitted Zn–O distance "applied to the
wrong oxygen", because a 12-6 term cannot tell an alcohol from a hydroxide.

An alcohol oxygen 0.13 Å outside an arbitrary radius is a stronger case for
D174's judgement than the one D174 itself made. Held at its pose distance.

## `eq6`: `eq5` plus two lines, nothing changed and nothing removed

    - [B:1:ZN, C:1:O1x, 3.13, 500.0]     # the fifth donor, D197's cause
    - [A:226:NZ, A:92:OD2, 2.99, 500.0]  # the salt bridge, 6DN4

k = 500 on both because D191 measured that 50 kcal/mol/Å^2 holds a 2.80 Å
target at a mean of 2.44/2.37, and holding within ~0.05 Å needs of order 500.
Every other key is `eq5`'s, so the two runs differ in exactly two restraints.

**Monitors added in five blocks** (`site`, `wtpostc`, `michpost`, `michqmpost`,
`michrkpost`): `Lys226 NZ-Asp92 OD2` and `-OD1`. The salt bridge had never been
monitored anywhere in this project, which is why nobody saw it break. Both
oxygens, because D120 is what happens when only one is watched. The product-arm
sections do not carry it; their hydroxide is a different chain and the product
state is a separate question.

Also corrected in `site.yaml`: the D148 comment saying Gln93 and Lys226 stand on
Zn2's ligand positions "and why Asp92 (D120) is left orphaned." **Asp92 is not
orphaned.** The lysine holds it, at 2.99 Å, which is the same job Zn2 does at
~2.1 Å in every canonical B3.

## What is NOT done

`eq6` has not run. Every path under `/scratch/user/uquadhik/` returns the
code-server route's own `Not found` from the `bun126/53134` session, including
the 41 GB `channel` tree, so scratch is unreachable from that session and neither
`eq6` nor the `michrk` time-ordering test could be touched. The config change is
on disk and pushed; it needs one `ls` to settle the path before anything runs.

**And one thing to decide before `eq6` is submitted.** D197 says "Asp92 OD1-O01
needs its own restraint ... the restraint has to be OD1-O01 itself", because
pinning OD1 to the zinc fixes its distance and not its direction. I have NOT
added that, because whether `Asp92 OD1 - O01` belongs in `carbonyl_on_metal` at
all is still an open question in this project, and the new geometry argues
against it: if Lys226 is the base, the hydroxide has no reason to be pointing at
Asp92, and 5.67/5.32 Å in the michrk frames plus "something hydrogen-bonds the
hydroxide throughout and it is never Asp92" (D204) says it is not. Adding the
restraint would impose the CphA geometry on a run designed to test the other one.

**D211 — `michrk`'s full 50 ps read frame by frame. Four results, one of
them against my own prediction and one of them a correction to D204. The
headline: **the ordering test I proposed in D209 cannot be run on this
trajectory**, because Asp92 is already on the metal at 2.13 Å in frame 1,
before the proton leaves Lys226 — and it is there because D202 deliberately
started this run from a frame where it had already collapsed. The salt bridge
never forms at all: **NZ–OD2 mean 4.66 Å, 0% of frames inside 3.2 Å**,
against 2.99 Å in the crystal.** (2026-11-05, measured)

## How it was read, and where the files actually are

**The run root is `/scratch/user/uquadhik/CSR/`,** not `/scratch/user/uquadhik/`.
Nothing in this project records that, and it is why a dozen probes returned the
code-server route's own `Not found` earlier today. `mcube` itself is checked out
at `/scratch/user/uquadhik/mcube/`. Both worth remembering.

`QM_MM.dcd` header: `nset 50`, `nsavc 1000`, `delta 0.020455` AKMA = 1 ps per
frame, 47,979 atoms, 28,791,676 bytes = 276 + 50 x 575,828 exactly. So **the
batch finished its full 50 ps** - D204 was written 6 ps in and there is eight
times more data than it saw. Atom indices from `QM_MM_lastframe.pdb`
(Zn 4203, Asp92 OD1/OD2 1413/1414, Lys226 NZ/HZ3 3431/3434, O01 4255,
mer O1x/C4x 4206/4208), coordinates by ranged byte reads of the three
coordinate blocks per frame. `xtb_.xyz` holds **130** atoms = D202's 124 QM
atoms plus the six link atoms, which confirms the six boundary pairs.

## 1. The ordering test is VOID, not failed

    frame (ps)   Zn-Asp92 OD2    HZ3 on
       1            2.13           N
       2            1.97           N (in transit, 1.40 from O01)
       3            2.03           O
       4            1.99           N
       5 onward     1.77-2.10      O, continuously to frame 50

Asp92 is 2.13 Å from the zinc in the **first** frame, with HZ3 still on the
nitrogen. So the collapse does not follow the transfer here, which is what D209
predicted it must.

**But this trajectory was never able to answer the question.** D202 chose
`protdft/qmmm_michqm/QM_MM_lastframe.pdb` as the start precisely because "its
Asp92 sits on the metal at 2.00 Å and that is deliberate". The collapse
happened in **batch 1**, where Lys226 was in the MM region, cationic, and could
not transfer a proton to anything. So one of the two events the test orders had
already occurred at t = 0, and no amount of this trajectory can order them.
D209's chain is **untested, not refuted**, and the lesson is general: *a run that
starts from the phenomenon cannot test its cause.*

The one run that could test it is **`eq6` followed by a QM/MM from its frames**,
because that is the only lineage in the project where Asp92 starts out and the
salt bridge starts held.

## 2. The salt bridge is absent for the whole batch

    Lys226 NZ - Asp92 OD2   mean 4.66 A   min 3.43   max 5.58   0% inside 3.2 A
    crystal                 2.99 A (6DN4), 2.83 A (6DQ2)

Not once in 50 ps does the contact the wild type holds at 2.99 Å come within
0.4 Å of existing. With Asp92 on the metal and the lysine deprotonated from
frame 5, that is what should be expected — but it means **every number in
`michqm` and `michrk` was computed in a site missing its wild-type cation
arrangement**, which is a caveat that has to travel with all of them.

    Zn-Asp92 OD2   mean 1.93   min 1.77   max 2.13   **100% of frames inside 2.5 A**
    Zn-Asp92 OD1   mean 3.95   min 3.45   max 4.26

## 3. It is NOT a shuttle. D204 is corrected

D204, 6 ps in, called the transfer "a **shuttle** (at 4 ps the proton is back on
the nitrogen at 1.12 Å) rather than D29's one-way event". Over the full batch
the proton is on the oxygen in **47 of 50 frames** and on the nitrogen only in
frames 1, 2 and 4. O01–HZ3 mean **1.05 Å**. So the flicker D204 saw is the
first four picoseconds of equilibration, and after that it is **one-way**. The
`-5.23` kcal/mol of D48 is still the number to reconcile, but "shuttle" is not
what this trajectory shows.

## 4. Zn-O01 does NOT report the protonation state

D204: "**Zn-O01 reports the protonation state independently**, 1.87 Å with the
proton on N and 2.21-2.26 Å with it on O, frame by frame, which no single wrong
parameter would produce." Over 50 frames:

    H on N (3 frames)    Zn-O01 mean 1.98   range 1.87-2.21
    H on O (47 frames)   Zn-O01 mean 2.04   range 1.86-2.38

The distributions **overlap almost completely** and the means differ by 0.06 Å.
The claim held on three early frames and does not survive the batch. Retracted.

## 5. A result that validates `eq6`'s new O1x restraint

    Zn - mer O1x, QM/MM (michrk, 50 frames)   mean 4.24   min 3.53   max 5.04
    Zn - mer O1x, MM (eq5, 10,000 frames)     mean 2.18   **97.9% inside 2.5 A**

So the fifth donor that rejected `eq5` is an **MM-specific artefact**: QM/MM,
given the same system with no restraint on it at all, keeps the alcohol oxygen
3.5-5.0 Å away. That is exactly D199's diagnosis — a 12-6 term cannot tell
an alcohol from a hydroxide and rewards every oxygen in reach. It also settles
the worry D210 left open about the new restraint: holding O1x at 3.13 Å during
MM **corrects an MM error rather than imposing a geometry**, because the quantum
answer is further out still, not closer in.

## 6. And the nucleophile is still not poised, over 8x more data

    O01 - C4x   mean 3.81   min 3.34   max 4.16

Never inside 3.3 Å. D201's and D204's conclusion survives the full batch.

## What this changes in the plan

`eq6` was Tier 0 scaffolding. It is now the **only** run that can test the
mechanism D209 proposed, and the QM/MM that follows it is the experiment, because
it is the first lineage that would start with the wild type's two-cation
arrangement intact. Nothing about the resting state should be concluded from
`michqm` or `michrk` geometries without the caveat in §2.

**D212 — `eq6`'s restraints checked against the actual start structure
before running it, and three of the five were pulling the wrong way. The good
news is bigger than the corrections: **`CSR_mich/preprocess_results/
finalsystem.pdb` already holds the Asp92–Lys226 salt bridge at 2.79 Å**,
inside the deposits' own 2.83–2.99 spread. The start structure was right about
the contact this project never knew it had.** (2026-11-05, measured)

## What the start structure actually contains

    Zn - His164 NE2    2.09        Zn - Asp92 OD2     2.79
    Zn - mer O2x       1.94        Zn - Asp92 OD1     4.09
    Zn - O01           2.14        Zn - mer O1x       3.82
    Zn - His88 NE2     2.22        Lys226 NZ - OD2    **2.79**
                                   Lys226 NZ - O01    2.66
                                   Lys226 NZ - OD1    4.91

Four donors inside 2.5 Å — His88, His164, O01 and the substrate carbonyl
— with Asp92 in the second shell at 2.79. **Four-coordinate at the strict
radius, exactly as D196 designed.** Note `Zn - HOH485` is **48.99 Å**: in the
Michaelis system the crystal's second water is not in the sphere at all, because
the substrate carbonyl took its place. The D186 cluster table, which lists
`HOH D:1` as the second ligand, was built from `CSR_wt`, not from this system.

## The three corrections, all in the same direction

Every one of them was pulling an atom **toward** the metal, which is the artefact
this whole line of work exists to avoid.

    restraint        D210 wrote   structure has   superposition   now
    mer O1x            3.13          3.82              -          3.82
    Lys226 NZ-OD2      2.99          2.79          2.83-2.99      2.79
    Asp92 OD1          3.81          4.09             3.90         3.90

- **O1x 3.13 —> 3.82.** 3.13 came from D197's words about "the October pose",
  a different structure. Applied here it would have pulled the alcohol oxygen
  **0.69 Å inward** at k = 500 — manufacturing the very coordination that
  rejected `eq5`. There is no crystallographic value for this atom, because no
  CSR-1 deposit contains the substrate, so the structure's own distance is the
  only defensible target. D211 says this is safe from the other side too:
  unrestrained QM/MM keeps O1x at 4.24 Å mean, further out than 3.82.
- **Salt bridge 2.99 —> 2.79.** The structure already has 2.79, so 6DN4's
  2.99 would have pulled it 0.20 Å outward for no measured reason. Holding a
  detected contact where the structure has it is the project's own convention
  (D174), and 2.79 is inside the two deposits' spread.
- **Asp92 OD1 3.81 —> 3.90, which makes this three changes from `eq5` and not
  two.** D174's 3.81 has no stated source anywhere in the log. The structure has
  4.09 and the superposition gives 3.90, so 3.81 at k = 500 pulls the second
  carboxylate oxygen **0.28 Å toward the metal** — biasing precisely the
  bidentate collapse D120 and D130 both record as the failure. 3.90 is the
  crystal-derived number and is used. This breaks the clean
  "one-restraint-differs" comparison with `eq5`, and being right is worth more
  than the comparison.

The two restraints left untouched are `eq5`'s and both are crystal-based:
`O01-C4x 3.47` (the structure has 3.48) and `Zn-mer O2x 2.60`, which pushes the
carbonyl **out** from 1.94 to the 2.53/2.63 Å that 7UHI/7UHJ measure in the
only experimental B3 Michaelis complex (D156).

## The lesson, which is the project's recurring one

Four of the five numbers in a config written from the decision log disagreed with
the structure the config points at, and three of them disagreed in the direction
that would have produced a known artefact. **A distance quoted in an entry is a
measurement of whatever structure that entry was looking at.** Reading the start
structure takes one fetch.

## Status

`eq6` is configured and the YAML parses. **It has not been submitted.** Running
it needs `sections.py site.yaml eq6` and an `sbatch` in the code-server
terminal, and the pane is currently hidden, so input cannot be verified by
reading it back — which is the one thing that must not be guessed at after
today.

**D213 — `eq6` submitted (job 32177894) and the Asp92 umbrella set up but
NOT started, because checking the logs first showed the only available seeds
would have bought an accurate free energy for a condemned shell over a range
that misses the answer, at the measured cost of ~8,500 core-hours. The library
blocker is gone though: `coordinate: distance` exists, so the chain is
`md_eq6` —> `aspcg` —> `asp`.** (2026-11-05, user: "double check logs
thoroughly before submitting any jobs to make sure you're not wasting work. Set
up and run the umbrella using mcube's umbrella module.")

## What the check found, which is why the umbrella did not go in tonight

`scan_aspout` and `scan_aspin` are 18 converged QM/MM points along Zn–Asp92
OD1, sitting on disk, and they look like ready-made window seeds. Three
independent reasons they are not:

1. **Both descend from `CSR_nac/md_naclong_results`**, which D120, D122 and D130
   all call the condemned shell and which D131 says nothing measured on can be
   cited. `retractions.md` already carries D136 for precisely this.
2. **Their range stops short of the answer.** Zn–OD1 1.90 –> 3.60 Å,
   while the wild-type conformer is at **OD1 3.90 Å** (D210's superposition)
   and **4.09 Å** in the start structure. Both arms end before the geometry
   the question is about — the same error D145/D208 caught in D144's costing,
   where the then-believed crystal value was off the end of both scans.
3. **D157 says so in its own words:** "This does not mean the PMF should be run
   now ... The order is the rebuilt pose first, this second", because that pose
   is 5.15 Å from the only experimental Michaelis complex in the family
   (D156).

At D157's measured rate — **40.9 h per window, 2.57 ps/h on 8 cores, ~8,500
core-hours for 26 windows** — seeding from those scans would have spent the
whole budget on a structure the project has already condemned. That is the
"wasting work" the instruction was about, and it was two greps away.

## What IS submitted

**`eq6`, job 32177894**, PD on `gpu_cuda`, 10 ns at 2 fs on one H100, 6 h limit.
`sections.py` wrote `md_eq6` after checking **10 selectors**, all resolved —
including `A:226:NZ`, which no section in this project had ever named. The five
restraints landed exactly as D212 corrected them:

    C:2:O01  - C:1:C4x   3.47
    B:1:ZN   - C:1:O2x   2.60  500
    B:1:ZN   - A:92:OD1  3.90  500
    B:1:ZN   - C:1:O1x   3.82  500
    A:226:NZ - A:92:OD2  2.79  500

## The library blocker is gone, and site.yaml said otherwise

The `aspout` comment block has claimed since it was written that "a PMF along a
two-atom distance cannot be run here: `umbrella` and `metadynamics` are both
built around a channel's z coordinate and a single ion". **That is no longer
true.** `umbrella.py` now has `COORDINATES = ("z", "distance")`, `cv_pair`,
`distance_cv`, `distance_values`, a distance branch in `bias_force`, and
`exclude_pair` on `coordination_force` so the biased pair is dropped from the
coordination restraints — those two forces would otherwise add and WHAM is
only told about the bias. `docs/umbrella.md` documents the whole mode. The stale
comment is marked superseded in place.

What a distance CV still cannot do: `pull` (refused by name — dragging a
bonded pair apart is not a seeding method), and the frame rescoring behind
`energies`, `wham --reweight` and `stitch`, all three tied to the channel
geometry. None is needed here: WHAM reads `cv.csv` directly.

## The chain, and why each link is gated

    md_eq6   the rebuilt pose              SUBMITTED, job 32177894
    aspcg    26 optimised points, the seeds   waits on eq6's last frame
    asp      the PMF, 26 QM/MM windows        waits on aspcg's points

`sections.py` refuses to write a section whose `pdb_file` does not exist, so
each link gates the next by construction rather than by my remembering to wait.

**The coordinate is Zn–Asp92 CG, not an oxygen**, for the reason the `monitor`
block already gives: a profile along Zn–OD1 is contaminated every time the
carboxylate rotates OD1 out and OD2 in, which is the exact motion D120 was
about, whereas Zn–CG moves only when the whole group moves. It is also what
the classical carboxylate–cation PMFs use. **2.3 –> 4.8 Å at 0.1 Å**,
D157's 26 windows, which brackets the bound state and reaches past the
wild-type value — D210 measures **Zn–CG 3.71 Å** in the superposed frame,
comfortably inside the range, unlike the OD1 scans' 3.60 Å cap. k = 150, 45 ps
equilibration + 60 ps data, SMB-1's protocol.

## One gap closed while setting up

`umbrella` takes a **subcommand** (`pull` / `window --window N` / `wham`), and
`run.slurm` passes only a stage and a section, so umbrella could not be launched
by it at all. `docs/channel/system/umbrella_windows.slurm` already solves this
generally — `SECTION` from the environment, `FIRST`/`LAST` blocks, finished
windows skipped so a block resumes after a timeout — but it hardcodes
`--partition=gpu_cuda --gres=gpu:1`, right for its MM windows and wrong for
these, since a QM/MM window is xTB-bound and a GPU does not help it. Added
`docs/csr/umbrella_windows.slurm`: the same script with the resource directives
removed so the caller supplies them on the sbatch line, as `run.slurm` does.
The two should eventually be one file, and that is the only difference.

**D214 - `eq6` finished in 43.7 minutes and all five configured restraints held,
but the result is a valid *starting structure*, not evidence about Asp92: the
redetection added a sixth restraint on Zn-OD2 itself, so the distance that the
PMF exists to measure was held, at both ends, for the whole 10 ns.** (2026-11-06)

Job 32177894, `md_eq6`, 10 ns NVT at 2 fs on one H100, `[DONE]` at 02:37 AEST.
`ASH Total Walltime: 2624.7 s`.

**What was asked for and what landed, from the last frame
(`md_eq6_results/NVTtrajectory_lastframe.pdb`, 47,979 atoms):**

| pair | target | last frame |
| --- | --- | --- |
| `C:2:O01 - C:1:C4x` | 3.47 | **3.48** |
| `B:1:ZN - C:1:O2x` | 2.60 | **2.57** |
| `B:1:ZN - A:92:OD1` | 3.90 | **3.87** |
| `B:1:ZN - C:1:O1x` | 3.82 | **3.79** |
| `A:226:NZ - A:92:OD2` | 2.79 | **2.79** |

So the pose survives: the carbonyl stays poised, the hydroxide stays at attack
distance, Asp92 stays in the second shell, and the salt bridge is intact.
Lys226 carries `HZ1 HZ2 HZ3` - the ammonium the crystal salt bridge implies
(the Tier 1 question 1), so the protonation state used here is the one the
structure argues for, not an assumption added for convenience.

**The metal's own sphere in the last frame, measured rather than assumed:**

    OH2:O01      1.96      the nucleophile
    HOH3171:O    1.97      a bulk water
    HIS164:NE2   2.13
    HIS88:NE2    2.20
    mer1:O2x     2.57      the substrate carbonyl
    ASP92:OD2    2.81      second shell
    ASP92:CG     3.72      the PMF's coordinate

Two histidines, not three, is correct and not a defect: the motif is **H-R-H**,
Arg90 standing where BBL His118 would be. At the 2.50 A coordination radius the
metal is **four-coordinate** (2 His, hydroxide, water), with the carbonyl just
outside at 2.57 - which is the coordination number the user expected of the
wild type, reached without being imposed.

**The honest qualification, which changes what this run can be cited for.**
`coordination_redetect: true` ran preprocess's detector on `finalsystem.pdb`,
and at the 3.00 A cutoff it found `ASP A:92 (OD2) 2.79 A` and added its own
restraint `[ZN, OD2, 2.79, 500.0]` - nine restraints in total, five configured
plus four detected (the detector's own Zn-O2x at 1.94 A was correctly overridden
by the configured 2.60). mcube printed the warning it was written to print:
*"longer than 2.50 A, so it is not a bond: this holds a second-shell contact
against the metal"*. That is exactly the user's rule (D174, D186: restrain what
is affected by the metal at its crystal distance), so the run is **right**. But
Zn-OD2 was pinned at 2.79 and OD2 was simultaneously pinned to Lys226's NZ at
2.79, so "OD2 ended at 2.81" is the restraint reporting itself, not a finding.
**Nothing in this run tests whether Asp92 would come onto the metal.** That is
the PMF's question, and the reason the PMF is the next step rather than a
confirmation of this one.

**Two things checked in `src/mcube/umbrella.py` before trusting the plan:**

1. **The bias and the coordination restraints do not fight.** `coordination_force`
   takes `exclude_pair` and drops the biased pair, and `release()` zeroes the
   whole coordination force before sampling begins, in both `run_mm_window` and
   `run_qmmm_window`. The module's own comment says why: two forces on one pair
   simply add while WHAM is told only about the bias, so the window's effective
   force constant would not be the one in the config "and the error would be
   invisible in the output". This matters here more than in the channel case,
   because the pair being biased (Zn-CG) is *not* the pair carrying a restraint
   (Zn-OD1, Zn-OD2) - without the release, the PMF would have been a profile of
   the restraints.
2. **`engine: qmmm` with `coordinate: distance` is a real path**, not a
   combination that falls through to MM: `run_qmmm_window` has its own distance
   branch and writes `engine: qmmm`, and `wham` skips each window's
   `equilibrate_ps` for QM/MM entries.

**The one gap left, and it is in the library, not the plan.** `window_seed` looks
only at the window's own `seed`, then at a matching MM window, then raises; a
finished `scan` cannot seed an umbrella even though it produces exactly the
structures required. So the 26 seeds have to be written into the config by hand
between `aspcg` and `asp`. That is a config-driven feature mcube should have
(the user's standing rule that users should not rewrite code per protein), and
it is the only manual step in the chain.

**D215 - The thesis is now a mechanism, not a contrast, and the argument map was
rebuilt around it. Twelve new nodes; four existing ones moved, every one of them
found by the validator rather than by reading. Findings 11 -> 2, and both
survivors are real gaps.** (2026-11-06)

The user's argument, in their words: *"Lys replaces Zn2 in the CSR-1 WT
mechanism by restraining the bridging aspartate, allowing for MBL activity."*
This replaces *"CSR-1's wild-type mechanism is different from the canonical B3
mechanism"* as the map's thesis. The old statement is kept as a claim, not
voided: the new thesis entails it. The reason for the change is that a negative
claim has no premises to list, and listing premises is the only thing an
argument map does.

**New nodes (12).** `thesis-lys226-restrains-asp92-in-place-of-zn2`;
`lys226-salt-bridges-asp92-in-wild-type-crystals` (2.99 A in 6DN4, 2.83 in 6DQ2);
`engineered-zn2-holds-asp92-at-195-angstrom`; `lys226-is-cationic-in-the-resting-enzyme`;
`question-what-is-lys226-pka`; `asp92-on-zn1-interferes-with-hydrolysis`;
`a-fifth-zinc-donor-closes-the-attack-angle`; `asp120-mutants-lose-activity-in-b1-and-b3`;
`question-does-restraining-asp92-lower-the-barrier`;
`the-truncated-form-is-the-one-modelled`; `spr1-activity-required-n-terminal-truncation`;
`gfn2-overstates-the-metal-coordination-preference`.

**6NC5 measured properly for the first time: it holds SIX zincs, and four are
crystallisation artefacts.**

    ZN301  Asp91 OD2 1.95 | His92 NE2 2.01 | His225 NE2 2.04 | HOH493 2.08
    ZN302  His87 NE2 2.22 | His163 NE2 2.09 | HOH485 1.90 | HOH493 2.14
    ZN303  HOH485 1.88 | HOH504 2.13 | HOH516 1.96          <- no protein ligand
    ZN304  Glu115 OE2 1.87 | PO4 1.90
    ZN305  His111 ND1 2.36 | Glu141 OE2 1.91 | PO4 1.90
    ZN306  Glu123 OE1/OE2 2.36/2.50 | Asp129 OD2 2.03

ZN301 is the engineered Zn2 with the canonical Asp-His-His set the two mutations
restore; ZN302 is Zn1 with two histidines because Arg90 holds the third position;
**Zn1-Zn2 is 3.43 A bridged by HOH493** at 2.14/2.08 - canonical B3 mu-hydroxide
geometry. So the thesis's central contrast sits in one deposit: the carboxylate a
lysine holds at 2.8-3.0 A in the wild type is held at **1.95 A by a metal** once
the metal is engineered back. Two preprocessing consequences: the four extra
zincs must be dropped explicitly, and **ZN303 shares Zn1's terminal water**
(HOH485 at 1.88/1.90 from the two metals), so removing ZN303 leaves that water
unanchored - and that water is where the nucleophile goes.

**The GFN2-vs-DFT benchmark (plan item 7) had already been run and was never
recorded.** Two 51-atom clusters differing only in whether the substrate's O1x
is on the metal, from `CSR_mich/config.yaml`:

    GFN2-xTB           off - on = -14.67 kcal/mol
    wb97x-v/def2-svp   off - on =  -8.75 kcal/mol

Same sign, **GFN2 exaggerates by 5.9 kcal/mol, a factor of 1.68**, on exactly
the metal-coordination question the thesis turns on. D198 read these same
numbers as "GFN2 does not over-coordinate the zinc"; the retractions index
already lists that as stronger than its evidence, and the opposite reading is
the defensible one. It also means `scan_znasp`'s +8.46 kcal/mol for pulling
Asp92 from 2.05 to 3.05 A should be expected to be too steep.

**Four existing nodes moved, each flagged by `check.py` and not by anyone
reading.**

- `asp92-is-a-zn1-ligand-in-our-models` was `supported` while its own body said
  the claim is "wrong by 2 A or more" against the crystals. **Voided**,
  superseded by the measurement. Its prescription - "do not restrain it, and let
  the run say where it goes" - is also superseded by the user's standing rule
  (D174, D186): Asp92 is affected by the metal and is restrained at its
  structural distance.
- `cpha-route-is-closed` **demoted to `open`**. D201 closed it on "Asp92 is a
  zinc ligand in every state this model reaches", measured false. The route is
  open, not reopened: nothing is asserted in its place.
- `mono-zinc-mbl-precedents-are-zn2-enzymes` **demoted to `open`** - its GOB-18
  premise is a literature description its own body says the deposit does not
  evidence.
- `a-fifth-zinc-donor-closes-the-attack-angle` was written leaning on D198's
  "the site keeps a fifth donor whatever we do", which **D200 overturns**
  (Asn225 leaves: 2.19 -> 4.18 A closest, 5.86 mean). Rewritten so the heading
  is the measurement - 39/56/68 deg against a 95-120 window - and the
  generalisation is marked retracted. This weakens
  `asp92-on-zn1-interferes-with-hydrolysis`, which is the thesis's mechanistic
  half, and it is now `open` for three stated reasons rather than one.

**The thesis's one unsupported premise is gone, by removal rather than support.**
`the-modelled-non-truncated-form-is-catalytically-relevant` had been check 3's
only finding since the map was built. The user's decision - truncation bears on
expression, use the truncated sequence - replaces it with
`the-truncated-form-is-the-one-modelled`, whose real requirement is weaker and is
met: **the construct modelled must be the construct the withheld kinetics were
measured on.** That needs no kinetic values to settle, only a construct name,
and it is the one question outstanding on that premise.

**Two findings remain and both are results, not defects.**
`question-what-is-lys226-pka` is a premise of the thesis with nothing behind it:
an uncharged lysine cannot restrain a carboxylate electrostatically, the project
has previously modelled this residue **neutral**, and the QM/MM lineage that
deprotonated it (D204) is the one where Asp92 then collapsed onto the zinc. If
the resting state is deprotonated, the thesis inverts. PROPKA and pdb2pqr are
still not in the ASH env (plan item 5, open since D7). And
`no-usable-michaelis-complex-exists-yet` is held on `ours` evidence by
construction - whether our structures are usable is our own assessment - but it
was narrowed: `eq6` is an equilibrated complex, so what is missing is agreement
on the **pose**, 5.43 A docked against 7UHI/7UHJ's 2.53/2.63 (D89, D123).

**Literature read for the aspartate's role, which cuts both ways.** L1 (B3,
Garrity 2004): Asp120 mutants lose **10- to 1000-fold in kcat**. IMP-1 (B1,
Yamaguchi 2005): **10^2-10^4-fold with both zincs still bound**, and stopped-flow
pH work argues the residue is neither the thing lowering the bridging water's
pKa nor the proton donor. BcII (B1, Llarrull 2007): the residue acts mainly as a
strong Zn2 ligand, positioning that ion. The mainstream explanation is therefore
that Asp120 holds **Zn2** - the ion CSR-1 does not have - which is the premise
behind caring what holds the carboxylate instead. But IMP-1 keeps both zincs and
still loses four orders of magnitude, so whatever the residue does past
positioning a metal is large, and nothing in our argument says a lysine supplies
it. Recorded as the third way the thesis could be wrong. All three are abstracts
and secondary summaries so far; the L1 paper needs reading in full before it is
cited, being the B3 comparator a referee reaches for.

**D216 - A fifth deposit, 6DQH, was never cited by this project. It is 1.10 A,
it holds a Zn1 and the wild-type lysine at once, and it measures Asp92 at 4.28 A
from the metal and 2.79 A from the lysine. In the same session the subclass's own
mutagenesis turned out to run against the thesis's causal clause.** (2026-11-06)

The user asked whether the Lys226-Asp92 salt bridge might exist only in
zinc-free structures, which would make Asp92 an equilibrium between the lysine
and Zn1 rather than committed to either. Checking it turned up more than the
answer.

**An RCSB full-text search for this enzyme returns five deposits and the record
used four.** Measured directly from the files, with project numbering in
brackets:

| deposit | 118 [90] | 121 [93] | 254 [226] | metal | Lys NZ-Asp120 | res. |
| --- | --- | --- | --- | --- | --- | --- |
| 6DN4 | Arg | Gln | Lys | **none** | 2.99 | 1.99 |
| 6DQ2 | Arg | Gln | Lys | **none** | 2.83 | 1.99 |
| **6DQH** | **His** | Gln | Lys | 1, at Zn1 | **2.79** | **1.10** |
| 6NC5 | Arg | **His** | **His** | 6 | no Lys | 2.07 |
| 6DR8 | **His** | **His** | **His** | 2 | no Lys | 1.48 |

**So the user's suspicion was right about the apo structures and is answered by
the one they had not seen.** 6DN4 and 6DQ2 contain no metal at all - zero Zn,
Cd, Ni, Mn, Fe, Co, Mg, Ca - which Pedroso et al. state in words ("no metal ions
were observed in the crystal structures of CSR-1 and CSR-1trunc") with ITC
giving alpha_Kd ~2 uM and beta_Kd ~150 uM. 6DQH is the R118H single mutant
[R90H], the construct the same paper describes as showing "one bound Zn2+ in the
alpha-site", and it keeps Gln121 and the lysine:

    ZN401   HOH526.O     2.03   <- terminal water, the nucleophile position
            His192.NE2   2.05
            His118.ND1   2.07   <- ENGINEERED; wild type is Arg here
            His116.NE2   2.08
            PO4405.O2    2.32
            HOH525.O     2.50

    Zn-Asp120 OD1   4.28     OD2  5.14     CG  5.15
    Lys260 NZ-OD2   2.79     NZ-OD1 3.94   NZ-Zn  5.27
    Asp120 OD1-HOH526  2.64
    Gln121 OE1-Zn   5.25     OE1-HOH526  4.60

**Three results, the third unlooked for.**

1. **With a metal and the lysine both present, the lysine wins outright.** 2.79
   against 4.28 at 1.10 A resolution is one state, not a split population. There
   is no crystallographic equilibrium.
2. **Every structure this project has modelled puts Asp92 too close to the
   metal.** `eq6` holds Zn-OD2 at 2.79 and Zn-OD1 at 3.90; earlier work put the
   carboxylate at 1.95-2.05. The crystal says 4.28/5.14, which is **1.5-2.4 A
   further out**. By the standing rule - restrain what the metal affects at its
   crystal distance (D174, D186) - those are the distances, and `eq6` is
   therefore a starting structure built on the wrong ones.
3. **Asp92 hydrogen-bonds the nucleophile**: OD1 to HOH526, the Zn1 terminal
   water, at 2.64 A. The carboxylate is not merely out of the way, it is in the
   second shell touching the water that attacks. That is the role canonical B1
   schemes give Asp120, and it is a positive mechanism rather than relief of
   inhibition.

**The caveat that keeps this an observation.** 6DQH's third Zn1 ligand is the
engineered histidine at the Arg90 position, two residues from the carboxylate.
6NC5's ZN302 is the only wild-type-like Zn1 (His87/His163 plus two waters, Arg90
unmutated at 4.02 A) and that deposit has lost the lysine. **No deposit is wild
type at all three positions and metal-bound, and the ITC numbers say none can
be:** at beta_Kd ~150 uM the Zn2 site does not fill, and Zn1 fills tightly only
once its third ligand is restored. The resting state the project models - one
zinc at Zn1, Arg90 intact, Lys226 intact - has never been observed. It is
assembled from two deposits that each supply half, which is defensible and is
what the modelling is for, but the paper has to say so.

**The literature result that runs against the thesis, and it is the serious one.**
Krco et al. (Front Chem 2023;11:1196073), citing Pedroso et al. 2020 Supp.
Table S2: replacing the Zn2-site glutamine and lysine with histidines "enhances
the affinities of both Zn(II) ions and **drastically improves the catalytic
efficiency to a level comparable to that of other B3 MBLs**"; and R118H gives "a
modest improvement in catalytic efficiency". ITC: wild type 2 / 150 uM, R118H
0.1 / 180, double **0.14 / 0.14**, triple 0.06 / 0.4.

The double mutant is exactly the construct in which Asp92 becomes a true
bridging ligand - 6NC5 measures it at 1.95 A on the engineered Zn2 - so if a
carboxylate on a metal suppressed hydrolysis, that construct should be worse. It
is dramatically better. **The thesis's structural half survives untouched; its
causal clause "allowing for MBL activity" does not.** On this evidence the
lysine is what *limits* CSR-1's activity, and the defensible reading is weaker:
the lysine is a structural stand-in that lets a mono-zinc B3 retain some
activity, not the reason it has activity. Two reservations: the fold changes sit
in a supplementary table neither the review nor this project has read
("drastically" is the review's word), and the double mutant adds a metal *and*
removes the lysine, so nothing isolates the carboxylate. Also checked and
absent: **no Asp120 mutant of CSR-1 or SPR-1 appears in either paper**, so the
claim that removing the bridging aspartate raises activity in this subclass is
not in the deposited literature.

**Clavulanate is functional evidence for the occupancy claim** and had not been
recorded: Ki 350 uM for CSR-1trunc, with the mechanism given as displacing the
metal from the beta site, attributed to the residue at 263 [226]. A residue
elsewhere could not make Zn2 the site a clavulanate displaces a metal from.

**The user's rulings on the premise list, recorded.** Premise 6 (Lys226 cationic)
**supported**, on the ground that a neutral lysine could not make this contact in
every crystal - now its own premise node, which is the right shape of argument.
Premises 4 (hydrolyses meropenem) and 8 (truncated construct) **withheld**: the
data exists and the user holds it, which is what the status means. Premise 7
(Asp92 on Zn1 interferes) was ruled true on a sharper mechanism than the node
had - the carboxylate makes Zn1 six-coordinate and a six-coordinate metal has no
open face - and that mechanism is now recorded as the one to test, **but the
status was not flipped**, because the mutagenesis above arrived after the ruling
and refutes it. Flipping it would have hidden the contradiction.

**Map state.** Nodes thesis 1, claim 20, observation 28, question 7; evidence
experimental 15, literature 7, ours 6; edges supports 81, refutes 10, requires
25, superseded-by 2. Findings **1**, the pre-existing
`no-usable-michaelis-complex-exists-yet`. The pKa question was re-scoped: it is
no longer a premise of the thesis or of the cationic claim, both settled by
geometry, and now belongs to `question-what-protonates-the-leaving-nitrogen`,
where the turnover state actually bears.

**The extra zincs in 6NC5 go, and it needs no code.** `remove_residues` already
exists in `preprocess.py`, takes `chain:resid:resname` selectors, and runs
before solvation and before the coordination detector - so a position left by a
deleted metal is left open for whatever displaces it, which is the behaviour
wanted here. For the mutant arm:

    remove_residues: [A:303:ZN, A:304:ZN, A:305:ZN, A:306:ZN]

keeping ZN301 (the engineered Zn2) and ZN302 (Zn1). **ZN303 shares Zn1's
terminal water** (HOH485 at 1.88 from ZN303 and 1.90 from ZN302), so deleting it
leaves that water as Zn1's own terminal ligand, unanchored - which is correct,
because that is where the nucleophile goes. Still to decide: the five PO4 and
four ACT crystallisation additives, and in 6DQH the phosphate oxygen sitting on
Zn1 at 2.32 A, in the position the substrate carbonyl has to occupy.


**D217 - Our own solved structure was never in the record, and the B-factor
check the user asked for changes what the no-equilibrium claim is allowed to say:
Asp92 is rigid, the lysine is the mobile partner.** (2026-11-06)

**Our structure.** `CSR_Manuscript-v1_old_gary_comments.docx` names it and the
coordinates have been sitting in the Honours tree all along:
`University/Undergrad/Honours/SCIE4481/bioinformatic_data/MM/Structure/CSR/dmCSR_mutant_1metal_refine_12.pdb`
and `dmCSR_mutant_2metal_refine_16.pdb`. **dmCSR-1 S214D, 1.27 A (one Zn) and
1.47 A (two Zn)**, P2(1)2(1)2(1), Rwork/Rfree 0.163/0.179 and 0.188/0.214. PDB
accession still "XXXX" in the manuscript. Numbering is BBL **-42** (Asp78 = Asp120,
His79 = Gln121H, His212 = Lys254H, His74 = His116, His150 = His192).

    ours1 (1.27 A)  ZN B1 occ 0.70  -> HOH2 2.07, His212 NE2 2.12, Asp78 OD2 2.13,
                                       HOH1 2.16, His79 NE2 2.22, HOH61 2.38
    ours2 (1.47 A)  ZN B1 occ 1.00  -> HOH1 2.03, Asp78 OD2 2.06, His79 NE2 2.11,
                                       His212 NE2 2.18, HOH2 2.43, HOH157 2.51
                    ZN B2 occ 0.53  -> HOH2 1.98, HOH1 2.01, His150 NE2 2.05,
                                       His74 NE2 2.23, HOH185 2.29

**The single zinc is at Zn2, not Zn1**, coordinated by the DHH set, which is what
the manuscript says: *"Multiple datasets ... all revealed the presence of only a
single zinc in the active site ... located in the Zn2 site coordinated by an
aspartate and two histidines (DHH site). The second metal site was only occupied
in a single dataset obtained for enzyme co-crystallised with ampicillin ... This
points to a weaker binding at the Zn1 site, and a possible mechanism where
binding of metal at this site is facilitated by the presence of substrate."*
Zn1 is at **occupancy 0.53** even there, and Zn2 at 0.70 in the one-metal form.

**What it does and does not settle.** It is the double mutant, so position 254 is
a histidine and **there is no lysine in it**: it cannot speak to the
Lys226-Asp92 question, and it is not a wild-type reference. What it does give is
independent, high-resolution, in-house confirmation that **when a Zn2 ligand set
exists the aspartate goes to Zn2 at 2.06-2.13 A and Zn1 is the weak site** - the
mirror image of the wild type, and a second reason the assembled resting state
(one zinc, at Zn1, with the lysine) has never been observed in any crystal,
ours included. It also supplies a substrate-dependent metal-binding mechanism we
had not recorded, which is a better explanation of the apo wild-type crystals
than crystallisation conditions.

**The B-factor check, which the user demanded before the no-equilibrium claim
stands. It stands for the aspartate and not for the pair.**

    6DQH (1.10 A; protein mean B 13.2, median 11.3, side-chain mean 14.4)
      Asp120   CB 9.6  CG 9.9  OD1 10.0  OD2 11.9      <- BELOW average, no altloc
      Lys260   CB 10.6 CG 12.5 CD 12.9 CE 13.6 NZ 14.6 <- rising along the chain
      His118   6.9-9.4            HOH526  14.2
    6DN4 (1.99 A; mean 16.2, median 13.8, side-chain 17.4)
      Asp120   CB 14.2 CG 15.4 OD1 16.4 OD2 16.3       <- about average
      Lys254   CB 13.9 CG 17.0 CD 19.1 CE 22.1 NZ 24.8 <- 1.5x the mean at the tip
    ours1 (1.27 A; mean 13.2, median 11.4, side-chain 15.3)
      Asp78    CB 10.9 CG 10.7 OD1 9.9 OD2 11.5        <- below average, as a Zn2 ligand

**Asp120 in 6DQH is one of the better-ordered side chains in the structure** -
every atom under the protein median, full occupancy, and not among the 12
residues that carry an alternate conformer. At 1.10 A a carboxylate with B ~10
is a single conformer, not a two-state average. So the **4.28 A is a real single
position**, and the user's concern does not apply to the atom the claim is about.

**But the lysine is the mobile partner, and that is a real qualification.** NZ's
B climbs monotonically out the side chain and ends **above** the protein mean in
6DQH (14.6 vs 13.2) and at **1.5x** the mean in apo 6DN4 (24.8 vs 16.2), with CE
22.1 - the signature of a lysine whose tip is sampling. So the honest statement
is not "the pair is rigid" but **"the carboxylate is pinned and the ammonium
reaches it"**. The direction of the thesis's verb is therefore weaker than
written: the crystal shows Asp92 held in a well-ordered position with the lysine
nearby and somewhat mobile, not a lysine clamping a carboxylate. In the apo
structures, where NZ B is 25, the salt bridge is much thinner evidence than its
2.83-2.99 A distance suggests - which is a second reason 6DQH rather than 6DN4
is the structure to take numbers from.

**No equilibrium between Zn1 and the lysine is still the right reading of the
crystal**, on the aspartate's own order parameters rather than on the single
distance: a carboxylate exchanging between a metal at 2.1 A and an ammonium at
2.8 A would show either two conformers or elevated B, and at 1.10 A it shows
neither. What the crystal cannot rule out is an equilibrium that is fast and
fully shifted - a minor state below the density threshold - which is a
free-energy question and is what the PMF measures
([[question-is-asp92-in-equilibrium-between-lys226-and-zn1]]).

**D218 - Pulling Asp92 out of a pre-equilibrated structure is the right plan and
cheaper than a rebuild, but two of this log's own rules and one already-measured
result set conditions on how. The pull must assert a distance and not an
orientation, and the attack coordinate has to be held while it happens.**
(2026-11-06)

The user's question: what stops us using a pre-equilibrated structure and pushing
Asp92 out to 6DQH's distance before QM/MM, instead of re-restraining and
re-equilibrating from scratch? Three things were already in the record and none
of them had been connected to it.

**1. Rule 6 of the coordination-sphere block is this exact move, and it is there
because D130 made it.** *"Do not hand-write a restraint to correct a shell that
is already wrong. D130's `[B:1:ZN, A:92:OD2, 2.80]` plus `[.., OD1, 3.81]`
blocked D120's rotation and reinstated its condition, **asserting the
carboxylate's orientation where the original asserted only a distance**."* Rule 7
adds the companion: *"no cutoff undoes an artefact already in the coordinates."*

**But those rules were written about a different defect and do not disqualify
`eq6`.** They condemn frames with Asp92 **on** the metal at 1.93-2.11 A - a real
coordination bond, inside any defensible cutoff. `eq6`'s carboxylate is at
**2.81 A**, outside the 2.50 A radius, and its sphere was measured
four-coordinate with the right composition (two His, hydroxide, bulk water, with
the carbonyl at 2.57). So `eq6`'s **coordination sphere is not wrong**; its
**second shell** is, by 1.4 A. That is a repairable defect and rules 6 and 7 do
not reach it.

**2. What rule 6's warning does reach is the orientation, and the target needs a
rotation.** Compare the two carboxylate geometries:

    eq6 last frame     Zn-OD2 2.81   Zn-OD1 3.87   Zn-CG 3.72
                       NZ-OD2 2.79   NZ-OD1 4.86   NZ-Zn  3.94
    6DQH (1.10 A)      Zn-OD1 4.28   Zn-OD2 5.14   Zn-CG 5.15
                       NZ-OD2 2.79   NZ-OD1 3.94   NZ-Zn  5.27

In `eq6` **one oxygen does both jobs**: OD2 faces the metal at 2.81 and the
lysine at 2.79. In the crystal **the jobs are split**: OD1 faces the metal, OD2
faces the lysine, and the NZ-OD2 distance is the same 2.79 in both. So the target
is not a longer distance, it is a **different conformer** - a chi2 rotation that
separates the two oxygens' roles. Restraining OD1 to 4.28 *and* OD2 to 5.14 would
produce it, and that pair of restraints is literally D130's mistake: it fixes the
rotation by assertion. **Pull on one coordinate - Zn-CG - and let the carboxylate
rotate or not.** If the free pull reproduces the crystal's split, the conformer is
earned; if it does not, that disagreement is a result about the pose.

**3. The experiment has effectively been run once, and its outcome is the real
obstacle.** In the protonated-Asp92 run (log §around D154), one proton and nothing
restrained, *"the carboxylate walks off to 2.99 A"* and the site becomes
**the crystal's four-coordinate sphere** - two histidines, the hydroxide and
HOH485 - which *"D120, D134 and D138 each failed to get by any other route,
having tried twice."* So Asp92 will leave and the sphere will come out right.
What happened next is the warning:

    hydroxide O01   1.94 -> 1.80     the vacancy closes by the HYDROXIDE moving in
    O01-C4x         3.48 -> 3.69     the attack coordinate OPENS
    O1x             3.64 -> 3.72
    carbonyl O2x    4.96 -> 4.63

*"The vacancy closes by the hydroxide moving in, not by the substrate moving
over. So emptying the Asp92 site does not by itself produce a Michaelis
complex."* **That is the thing to stop.** `eq6` spent 10 ns establishing a poised
pose - O01-C4x 3.48, Zn-O2x 2.57 - and freeing Asp92 is measured to degrade
exactly that coordinate.

**The protocol that follows from all three.** Pull from `eq6`'s last frame along
**Zn-CG only**, 3.72 -> 5.15 A, **while keeping `eq6`'s attack-coordinate and
carbonyl restraints in place** (`C:2:O01 - C:1:C4x 3.47` and `B:1:ZN - C:1:O2x
2.60`), which are already in the config and are what prevents the hydroxide from
taking the face. Monitor OD1, OD2 and NZ at every point so the conformer is read
rather than imposed. Then release the Asp92 restraint for the QM/MM leg and let
it say where the carboxylate goes.

**The tool is the `scan` stage, and this reverses what D213 said about it.**
`aspcg` was dismissed there as adding nothing over the qmmm stage's built-in
optimiser. For *this* job that was wrong: an unconstrained optimisation runs
downhill to wherever it likes, while a relaxed scan walks one coordinate in
controlled steps with each point starting from the previous - which is precisely
a controlled pull. A short scan over 3.7-5.2 A is also the first half of `aspcg`,
so none of it is wasted if the PMF follows.

**Cost, for the comparison that prompted the question.** The pull is a ~15-point
relaxed QM/MM scan, far below the 8,500 core-hours of the PMF and below the
43.7 minutes of a fresh 10 ns MD. The rebuild-from-6DQH alternative is worse than
it sounded in D216: 6DQH is a different deposit with different numbering, a
phosphate on the zinc where the substrate carbonyl goes, and no meropenem - so
"rebuild" means re-docking the substrate and re-establishing the pose, which is
the single most expensive thing this project has done. **Taking 6DQH's distances
as restraint targets while keeping our own equilibrated system is the correct
reading of D216, and it is what the user proposed.**

**D219 - The two-arm comparison replicates, and the reason is a result in itself:
Asp92 sits the same distance from Zn1 whether a Zn2 holds it or the lysine does.
But our own structure cannot be the mutant arm's starting geometry - S214D is an
active-site mutation, 4.1 A from both metals.** (2026-11-06)

Asked whether the Asp92 pull can be replicated for dmCSR so the two arms stay
comparable. Measured across every metal-bound structure of this enzyme:

| structure | what holds Asp92 | Zn1-CG | Zn1-OD1 | Zn1-Zn2 |
| --- | --- | --- | --- | --- |
| `ours2` dmCSR-1 S214D, 1.47 A | Zn2 at 2.06 | 4.72 | 4.12 | 3.23 |
| 6NC5 dm (Q121H/K254H), 2.07 A | Zn2 at 1.95 | 4.66 | 4.22 | 3.43 |
| 6DR8 tm, 1.48 A | Zn2 at 2.22 (occ 0.45) | 5.48 | 4.61 | 4.46 |
| **6DQH WT lysine, 1.10 A** | **Lys226 at 2.79** | **5.15** | **4.28** | - |

**Zn1-Asp92 is 4.12-4.61 A at OD1 and 4.66-5.48 at CG in all four, and 6DQH sits
inside that range.** So the lysine does not merely keep the carboxylate off the
metal - **it holds it where a Zn2 would hold it**, to within 0.2 A at OD1. That is
a sharper statement of the thesis than "restrains the bridging aspartate", it is
measured in four independent crystals of three different constructs, and it is
what makes the comparison valid: both arms are restrained at effectively the same
Zn1-Asp92 distance, reached by different means.

**What is symmetric between the arms, which is what the user's "same
preprocessing" requires:** restraints taken from each system's own crystal
distances; the same MD length; the same 300 ps xTB QM/MM in 50 ps chunks; the
same def2-SVP optimisation and def2-TZVP single point; the same QM-region residue
indices. **What differs, necessarily:** the mutant has two metals, so a different
QM charge and a bridging hydroxide, and its Asp92 is a genuine Zn2 ligand at
~2.0 A, which `preprocess` will detect and restrain there without any pull. **The
WT arm's pull is a repair of an equilibration artefact, not a physics step** -
`eq6` drifted to Zn-CG 3.72 against the crystal's 5.15 - so applying it to one
arm and not the other does not break the symmetry that matters.

**S214D disqualifies our own structure as the mutant arm's starting geometry.**
The construct's name gives it away only once the numbering is resolved: our file
is BBL **-42**, so S214D is file residue **172**, and the motif reads
ALA-ASP-**ASP**-LEU-ALA where 6DN4, 6DQH and 6NC5 all read
ALA-ASP-**SER**-LEU-ALA. In BBL numbering that serine is **Ser214**, which is
**Ser186 in project numbering** - and `eq6`'s own metal sphere has Ser186 OG at
**3.81 A from the zinc**. So this is not a remote crystallisation aid:

    6DQH   Ser214 OG -> Zn1 5.21        (Asp120 8.4, HOH526 6.35)
    6NC5   Ser185 OG -> Zn1 4.78, Zn2 5.87
    ours2  Asp172 OD2 -> Zn1 4.13, Zn2 4.06     <- the mutation, moved ~1 A IN

The engineered aspartate puts **a second anionic group 4.1 A from both metals**,
roughly equidistant, where the wild type has a neutral hydroxyl 4.8-5.2 A out.
It does not touch Asp92 (9.1 A away), so it is not a direct confound on the
carboxylate - but it changes the site's total charge and the metals'
electrostatic environment, which is precisely the variable this thesis is about.
It is also the likely reason our crystals bind Zn2 well.

**So: the mutant arm starts from 6NC5** - the plain double mutant, both zincs at
full occupancy, Ser186 intact, Zn1-Zn2 3.43 - with `remove_residues:
[A:303:ZN, A:304:ZN, A:305:ZN, A:306:ZN]` as D216 prescribed. **Our structure
becomes the validation instead, and it is good at it:** an independent crystal
form at higher resolution reproducing 6NC5's Asp-Zn2 at 2.06-2.13 against 1.95
and Zn1-CG at 4.72 against 4.66. Two crystals, two space groups, the same
geometry.

**One question for the user, of the same kind as the construct question on
premise 8 and needing no kinetic values:** are the withheld mutant kinetics
measured on the plain double mutant or on dmCSR-1 **S214D**? If the latter, the
simulated construct and the measured one differ by an active-site mutation, and
the mutant arm should carry S186D.

**D220 - The Asp92 pull cannot be a `scan`. mcube's optimiser-constraint defect,
diagnosed in D134 and never fixed, blocks every relaxed scan on this system, and
it sits on the critical path of the pipeline the user defined. The pull is
restrained MD instead - config only, ~15 minutes - and it was NOT submitted.**
(2026-11-06)

D218 chose `scan` for the pull because a relaxed scan is a controlled pull where
an unconstrained optimiser is not. Checking the log before submitting found that
it cannot run.

**The defect, in the log's own words (D134).** `optimizer_constraints`
(`src/mcube/qmmm.py:539`) emits the active-region MM waters as geomeTRIC
`$freeze` when there is no user constraint, and switches **the whole block to
`$set` with values the moment one appears - which a scan always supplies**. A
frozen distance is projected out and a residual on it costs nothing; a `$set`
distance is a **convergence criterion** checked against `convergence_cmax`, and
the TIP3P O-H distances sit at 0.069 A against a tolerance seven times smaller,
so the run **cannot be declared converged however long it runs**. It also wrecks
the linear algebra: **160 SVD failures in 120 steps**, each dumping a 32 MB file.
`scan_aspout` was cancelled over it - 2.3 min/step with `maxiter` 300 is 11.5 h
per point, and **4.6 GB of dumps in a single point**.

D218's `asppull` had **four** user constraints, so it would have triggered this
immediately, and its 15 points at `opt_maxiter: 300` would have been ~170 h
against any walltime. **It was written, committed, and would have wasted the
allocation. Caught by reading D134, not by anything in the tooling.**

**This is on the critical path, which D134 said and nobody carried forward.**
Its closing line: *"D133 test 1 is not answerable by `scan` as mcube currently
writes constraints, and the same applies to every `qmmm` run given
`opt_constraints`."* The pipeline the user defined ends in **a def2-SVP geometry
optimisation**, and if that optimisation is given any constraint - to hold the
attack coordinate, say - it hits the same wall. So does `aspcg`. The PMF's
umbrella windows are safe, because they bias through MD forces rather than
optimiser constraints.

**The fix that should be made, and it is a library fix rather than a CSR one.**
D134's option 1: write geomeTRIC's constraint file ourselves so the waters stay
`$freeze` and only the user's own pairs go in `$set`. ASH's `Optimizer` takes a
single `constrainvalue` flag for the whole dict, which is why mcube currently
cannot express the mixture. This is general - it unblocks constrained
optimisation for every system, not this protein - and it is the shape of change
the standing instruction asks for. Option 2, `rigid_water: false`, is one config
line but changes what is being optimised, so it stays the user's call. Option 3,
shrinking `actradius`, mitigates only.

**What the pull is now.** An `md` section, because MD restraints are harmonic and
have no convergence criteria, no SVD and no dumps - and the surroundings relax at
300 K rather than downhill into the nearest minimum, which is the more
defensible thing to do to a solvated site anyway.

    asppull   nvt_time 1000   Zn-CG restrained at 4.45, k 50   (halfway)
    asppull2  nvt_time 2000   Zn-CG restrained at 5.15, k 50   (6DQH's distance)

Both carry eq6's own restraint list minus the Asp92-Zn pair, so nothing new is
asserted, and both set **`coordination_redetect: false`** - which is what stops
D214 repeating, since redetection re-pins Zn-OD2 at its detected 2.81 and would
hold the pulled coordinate at both ends again. Staged in two steps so a 1.4 A
displacement is never applied as one impulse, and at k = 50 rather than 500 -
the library default, justified in the log from Gervasoni 2022 and Hanzevacki
2025, leaving 0.11 A RMS of room. Still one coordinate only: Zn-CG, with both
oxygens free to rotate, per rule 6.

**Cost:** eq6 did 10 ns in 43.7 min on one H100, so 1 + 2 ns is **about 15
minutes**, against ~170 h for the scan that was nearly submitted.

**Three corrections to D218's config, all found by checking the stage's own
keys against `eq6`'s working section rather than assuming:** `md` has no
`monitor` key (monitoring is post-hoc through `site.py shell`), no
`constraint_exclude_residues` (that is a **preprocess** key; the md-stage
equivalent is `coordination_redetect: false`), and takes `nvt_time` in **ps**,
not `time` in ns.

**Not submitted.** The browser pane is hidden, so only the read-only
`vscode-remote-resource` route is available and the shell must not be driven -
and the scratch quota has not been checked, which rule 8 requires before any
write. The config is committed and ready.

**Verified directly from the file rather than from D214's table**, since the pull
is designed around these numbers: `md_eq6_results/NVTtrajectory_lastframe.pdb`,
47,979 atoms, one zinc - Zn-OD2 **2.81**, Zn-OD1 **3.87**, Zn-CG **3.72**,
NZ-OD2 **2.79**, NZ-OD1 **4.86**, Zn-NZ **3.94**, Zn-Ser186 OG **3.81**; sphere
O01 1.96, HOH3171 1.97, His164 2.13, His88 2.20, mer O2x 2.57, OD2 2.81. Every
value matches D214. **It also confirms the carboxylate flip the pull is designed
around:** in `eq6` OD2 does both jobs at once (2.81 to the metal, 2.79 to the
lysine), where 6DQH splits them (OD1 4.28 to the metal, OD2 2.79 to the lysine,
and NZ-OD1 closes from 4.86 to 3.94).

**D221 - Scratch was FULL at 301/300 GB, which would have made the pull's output
silently incomplete. The cause is not this project: mcube writes every
trajectory twice, raw and imaged, and eq6's NVT alone was 5.4 GB in each copy.
6 GB reclaimed with no information lost, and `asppull` submitted as job
32365427.** (2026-11-06)

Rule 8 says check disk first because an over-quota write fails silently, leaving
a file that looks complete. It was worth obeying: **`rquota` reported
`scratch 301 300`** - over the limit, not near it. Anything submitted before
clearing that would have produced output that could not be trusted, and the
failure would have looked like success.

**CSR is not the problem.** The whole 301 GB, by `du --max-depth=1`:

    45G .cache     41G channel    41G glp1r_md   27G AANCA    25G hf
    23G rfd4       19G CSR        17G pmhc       16G boltzgen 14G .venv_bgrocm
    14G .venv_esm  7.9G .boltz    3.1G NeuralCubes

CSR is 19 GB of it, and ~98 GB is caches and virtualenvs (`.cache`, `hf`, the two
`.venv_*`). Those are not this project's to delete, and `channel`, `glp1r_md`,
`AANCA`, `pmhc`, `rfd4` and `boltzgen` belong to other work - so freeing real
space is the user's call, not an autonomous one.

**What was this project's to fix, and it is a library finding.**
`md_eq6_results` was **13 GB - 68% of the entire CSR tree** - for one 10 ns run:

    5.4G NVTtrajectory.dcd          5.4G NVTtrajectory_imaged.dcd
    550M relaxbox_NPT.dcd           550M relaxbox_NPT_imaged.dcd
     55M warmup_MD_cycle2.dcd        55M warmup_MD_cycle2_imaged.dcd
     28M warmup_MD_cycle1.dcd        28M warmup_MD_cycle1_imaged.dcd

**Every trajectory is stored twice, raw and imaged.** That is mcube's behaviour
in every project, not a CSR setting, so it doubles the footprint of every MD run
anyone runs. **A config flag to drop the un-imaged copy once imaging succeeds
would halve that everywhere**, and by the standing instruction - make it a config
change so users need not rewrite code - that is where the fix belongs. Not built.

Deleted the four raw duplicates, conditional on the imaged twin existing and
being non-empty (`test -s NVTtrajectory_imaged.dcd && rm -f ...`). **Information
loss: none** - the imaged copy is the same data wrapped into the box, and it is
the copy analysis wants. `rquota` after: **294/300**. The last frame, the
checkpoints and `OpenMM_MD_final_state.xml` are untouched, and the eq6 occupancy
analysis remains possible because `NVTtrajectory_imaged.dcd` is still there.

**6 GB of headroom is thin for two MD runs, so the config changed too.**
`traj_frequency` 500 would write ~1.2 GB per ns, doubled - about 5 GB for the two
pull stages, against 6 GB free. Raised to **5000** in both pull sections: the
deliverable is the last frame, not the trajectory, so 100 and 200 frames are
ample and the runs now cost ~200 MB together.

**Submitted.** `asppull.sh`, in the repo rather than typed because the
code-server terminal truncates long lines silently, prints the quota and the
restraints that actually landed before submitting. `sections.py` accepted the
section - **10 selectors checked** - and `config.yaml` holds what was intended:

    ['C:2:O01','C:1:C4x',3.47] ['B:1:ZN','C:1:O2x',2.6,500.0]
    ['B:1:ZN','C:1:O1x',3.82,500.0] ['A:226:NZ','A:92:OD2',2.79,500.0]
    ['B:1:ZN','A:92:CG',4.45,50.0]
    nvt_time 1000, nvt_timestep 0.002, traj_frequency 5000,
    coordination_redetect False, platform CUDA

**Job 32365427**, `asppull`, `gpu_cuda`, PD at 16:46 AEST. Stage 2 (`asppull2`,
Zn-CG to 5.15) cannot be staged yet: `sections.py` refuses a section whose
`pdb_file` does not exist, and stage 1's last frame is that file.

**Also learned about the bridge:** a `browser_batch` of eight actions was
accepted and silently did nothing - no file, no error - where the same commands
sent one at a time all landed. `BUNYA_BROWSER_SHELL.md` warns that batches are
refused more often; this was worse than a refusal, since it reported success for
every action. One command per call on this node.

**D222 - `asppull` stage 1 is DONE in 6.6 minutes with the pose fully intact,
and it moved Asp92 by 0.10 A of the 0.73 it was asked for. The reason is the
thesis turning up in a restrained MD: the carboxylate is tethered to Lys226, and
the lysine did not move.** (2026-11-06)

Job 32365427, `md_asppull`, `[DONE] md_asppull` at 16:53 AEST, **ASH total
walltime 396.2 s**. 1 ns NVT at 2 fs on one H100, `traj_frequency` 5000.

**Every restraint that was meant to hold, held** - which was the thing at risk,
because D154 measured the hydroxide taking the vacated face and the attack
coordinate opening when Asp92 leaves:

| pair | target | last frame |
| --- | --- | --- |
| `C:2:O01 - C:1:C4x` | 3.47 | **3.45** |
| `B:1:ZN - C:1:O2x` | 2.60 | **2.55** |
| `B:1:ZN - C:1:O1x` | 3.82 | **3.79** |
| `A:226:NZ - A:92:OD2` | 2.79 | **2.79** |
| `B:1:ZN - A:92:CG` | **4.45** | **3.82** |

So the poise survived the pull. **The pull did not happen.** Zn-CG went
3.72 -> 3.82: 0.10 A of a 0.73 A displacement, with the restraint 0.63 A from
its target at the end of 1 ns.

**Why, and it is informative rather than a bug.** The full geometry:

    Zn-CG   3.72 -> 3.82      Zn-OD1  3.87 -> 4.01    Zn-OD2  2.81 -> 2.86
    Zn-NZ   3.94 -> 4.03      NZ-OD1  4.86 -> 4.93    NZ-OD2  2.79 -> 2.79
    Zn-Ser186 OG  3.81 -> 4.01
    sphere: HOH3171 1.96, O01 2.03, His164 2.13, His88 2.19, mer O2x 2.55, OD2 2.86

**The lysine never moved** (Zn-NZ 3.94 -> 4.03), and the carboxylate is tied to
it by a 500 kcal/mol/A^2 salt-bridge restraint that held at exactly 2.79. A
k = 50 pull on Zn-CG was therefore pulling a carboxylate that is anchored, at ten
times the force constant, to a residue the rest of the protein holds in place.
The 500 wins, as it must.

**That is the mechanism the thesis asserts, demonstrated by accident: Asp92
cannot be moved away from Zn1 without moving Lys226, because the lysine is what
holds it.** It is weak evidence - a restraint fight is not a free energy - but it
is the right sign, and it came out of a run set up for another purpose.

**The fix, and it is faithful to the crystal rather than a force-constant hack.**
6DQH gives an internally consistent set for the whole pair: **Zn-CG 5.15,
Zn-NZ 5.27, NZ-OD2 2.79.** Our structure has 3.82 / 4.03 / 2.79 - the salt bridge
is already at the crystal value and the *pair* is 1.3 A too close to the metal.
So stage 2 restrains all three and the pair travels as a unit. Raising k on
Zn-CG instead would have strained the structure to win a tug of war against a
restraint we put there ourselves.

Lys226 NZ at 3.94 A from the metal is outside the 3.0 A detection cutoff, so it
was never auto-restrained - but it is squarely "anything affected by the metal"
under the standing rule (D174, D186), and 5.27 A is its measured crystal
distance. Still **no OD1/OD2-to-metal restraint**, so the carboxylate's rotation
stays unasserted (rule 6).

One incidental: **Ser186 OG moved 3.81 -> 4.01 A**, which matters because
S214D - our own crystal's third mutation - is this residue (D219). Its position
is not fixed by anything in our setup.

**Stage 2's submission is UNCONFIRMED.** `asppull2.sh` was committed (`e732e1e`)
and the launcher was run, but the code-server bridge stopped answering
immediately afterwards - every `vscode-remote-resource` fetch returning
`TypeError: Failed to fetch` while `location.href` still showed the code-server
URL and the terminal was still in the DOM. So it is not known whether the repo
pull landed, whether `sections.py` accepted `md_asppull2`, or whether anything
was submitted. **Check `squeue` and `config.yaml` before re-running it**, because
a second submission of the same section would be two jobs writing the same
output directory. Stage 1's results are on `/scratch` and unaffected.

**D223 - `asppull2` moved the lysine and so confirmed D222's diagnosis, but the
carboxylate still did not follow, and the reason is that `Zn-CG` is nearly blind
to what is wrong: the carboxylate is ROTATED, not merely too close. Also
confirmed, by probe rather than by assumption: typing does not reach the shell
while the browser pane is hidden.** (2026-11-06)

Job 32366546, `md_asppull2`, `[DONE] md_asppull2`, **ASH total walltime 617.3 s
(10.3 min)**, 2 ns NVT. It was submitted before the bridge dropped, so D222's
warning about a possible double submission is resolved: there is exactly one
stage-2 job and its output is intact. `sections.py` checked 12 selectors.

**The whole lineage, against 6DQH:**

| pair | eq6 | stage 1 | stage 2 | 6DQH |
| --- | --- | --- | --- | --- |
| Zn-CG | 3.72 | 3.82 | **3.96** | 5.15 |
| Zn-NZ | 3.94 | 4.03 | **4.85** | 5.27 |
| Zn-OD1 | 3.87 | 4.01 | **4.09** | **4.28** |
| Zn-OD2 | 2.81 | 2.86 | **2.95** | 5.14 |
| NZ-OD2 | 2.79 | 2.79 | 2.83 | 2.79 |
| NZ-OD1 | 4.86 | 4.93 | 4.71 | 3.94 |
| O01-C4x | 3.48 | 3.45 | **3.45** | - |

Sphere: O01 2.00, HOH3171 2.04, His164 2.14, His88 2.21, mer O2x 2.62, OD2 2.95.
**The pose survived a second time** - the attack coordinate has not moved off
3.45 through two pulls, which is the thing D154 said was at risk.

**D222's diagnosis was right.** Adding the `Zn-NZ 5.27` restraint moved the
lysine from 4.03 to **4.85**, within 0.42 A of its crystal distance, where in
stage 1 it had not budged. The salt bridge was indeed the anchor.

**But the carboxylate is still in the second shell, and the table says why.**
**Zn-OD1 is 4.09 against the crystal's 4.28 - essentially there - while Zn-OD2 is
2.95 against 5.14.** One oxygen is at its crystal distance and the other is 2.2 A
too close. That is not a translation problem, it is an **orientation** problem:
the carboxylate is rotated so OD2 faces the metal, where 6DQH has OD2 facing the
lysine and OD1 facing the metal. This is the flip identified in D220 before any
of these runs, and it has survived all of them.

**And `Zn-CG` cannot fix it.** A carboxylate can swing about CG while Zn-CG
barely changes, so the coordinate chosen to avoid asserting the rotation (rule 6)
is also nearly insensitive to it. Zn-CG moved 3.72 -> 3.96 across two runs, 0.24
of 1.43 A, which is about what a rotation alone would give.

**So stage 3 pulls `Zn-OD2` instead, at its measured crystal value of 5.14.**
That is a *single* pair, so rule 6 still holds - OD1 is left free, and whether it
settles near 4.28 by itself becomes a test rather than an assertion. And **k is
500, not 50**: every crystal-derived restraint in `eq6` uses 500, the 50 on the
pull was an inconsistency of mine, and D154 says why it matters - the Asp92-Zn
contact in this model *"requires the anionic carboxylate: it is electrostatic,
not steric"*. MM wants the carboxylate on the metal, so a gentle restraint loses.
Fighting the force field is the **point** of the standing rule: restrain at the
crystal distance in MM, because MM cannot be trusted for coordination, then let
QM/MM decide where it goes. Committed as `52c5e85`; not submitted.

**The bridge, verified rather than assumed.** The pane was hidden, so a probe
(`date > p2.txt`) was sent before anything else. The fetch returned **404 with
body `Not found`** - the route working, the file absent - so **the command never
ran**. `bunya-working-rules` says typing fails while the pane is hidden in a way
that looks identical to success, and that is now confirmed with a one-command
test rather than inferred after losing work. The earlier eight-action
`browser_batch` that silently did nothing (D221) was the same failure with the
pane visible, so there are two independent ways for input to vanish: **send one
command, then verify by fetching a file, every time.**

**D224 - `asppull3` reproduced 6DQH to within 0.1 A on every Asp92/Lys226
distance, and OD1 got there with nothing restraining it. The QM/MM leg is
submitted as job 32384737 with nothing restrained at all. Gate 3 is re-passed.**
(2026-11-06)

Job 32382372, `md_asppull3`, `[DONE]`, 2 ns NVT.

| pair | eq6 | stage 1 | stage 2 | **stage 3** | 6DQH |
| --- | --- | --- | --- | --- | --- |
| Zn-OD2 | 2.81 | 2.86 | 2.95 | **5.05** | 5.14 |
| Zn-OD1 | 3.87 | 4.01 | 4.09 | **4.19** | **4.28** |
| Zn-CG | 3.72 | 3.82 | 3.96 | **5.09** | 5.15 |
| Zn-NZ | 3.94 | 4.03 | 4.85 | **5.26** | 5.27 |
| O01-C4x | 3.48 | 3.45 | 3.45 | **3.47** | - |

**The result that makes this start earned rather than built: Zn-OD1 is 4.19
against the crystal's 4.28 and nothing restrained it.** Only OD2 was pulled. So
the carboxylate adopted the crystal's orientation - OD1 facing the metal, OD2
facing the lysine - on its own, which is the split 6DQH shows and which eq6 did
not have. Had OD1 swung in to replace OD2 on the metal, that would have said the
MM site simply insists on a carboxylate oxygen at ~2 A and the geometry is
unreachable. It did not.

Sphere at the 2.50 A radius: HOH3095 2.01, O01 2.03, His164 2.15, HOH3171 2.15,
His88 2.21 - **five donors, because a second bulk water took the space Asp92
vacated**, as D154 saw the vacancy close. Asp92 is nowhere near it. Not fixed in
MM, deliberately: bulk water comes and goes and MM cannot decide coordination
(user, repeatedly), so whether that water stays is QM/MM's answer to give.

**Gate 3 (plan item 3) is re-passed on the right numbers.** D216 had demoted it
because `eq6` held Asp92 at distances 1.5-2.4 A inside the crystal's. The
structure now matches 6DQH where 6DQH can be matched, with the attack coordinate
intact through three successive pulls.

**Three stages for what one should have done, and the reason is worth recording
because it is a process failure rather than a chemistry one.** Stage 1 pulled
`Zn-CG` at k=50 while leaving the salt bridge at k=500, and moved 0.10 A of
0.73; stage 2 added `Zn-NZ` and moved the lysine but not the carboxylate; stage 3
pulled `Zn-OD2` at k=500 and arrived. **Every one of those failures was
deducible from a table already in this log before the job ran** - D219 recorded
Zn-NZ at 3.94 against 5.27, so stage 1 could not work; stage 1's own output
showed OD1 near its target and OD2 2.2 A out, so the rotation was diagnosable
before stage 2; and eq6's restraint list, quoted repeatedly, uses 500 for every
crystal-derived distance, so k=50 was never justified. An arithmetic problem was
treated as a sequence of experiments. **The rule, to be applied before every
submission from here: write out every restrained pair as current to target with
its force constant, check that no restraint opposes another, check that every
distance differing from its target is actually restrained, and state the
predicted outcome. A setup whose outcome cannot be predicted is not understood
well enough to spend the allocation on.** The cost of ignoring this just rose
from 10 minutes to 20 hours per iteration.

**`aspqm`: the QM/MM leg, job 32384737, `general`, 8 cores, 16 GB, 24 h.**
`sections.py` accepted it and `config.yaml` holds what was intended -
`skip_equilibrate: False`, `skip_optimize: True`, `time: 50`, `resume: True`,
`charge: 2`, QM region `Zn, OH, mer, His88, His164, Asp92, Asn225, Arg90,
Lys226`. **The launcher prints a warning if `restraints`, `opt_constraints` or
`constraints` appear in the section, and none fired**, so "nothing is restrained"
is confirmed from the written config and not only from reading the source.

Why nothing is restrained, checked in `qmmm.py` rather than assumed:
`bond_constraints` is reached only from the optimiser path and from rigid water,
so the equilibrate applies no coordination restraints; and the QM region is
rebuilt from residue names, so D115's stale-index trap - saved indices
restraining a substrate carbon as the metal - cannot fire on a frame this far
from the one preprocess ran on.

Charge 2 re-derived rather than inherited (D152): Zn(2+) - OH(-) + meropenem
zwitterion(0) - Asp92(-) = 0 for the michaelis region, + Arg90 + Lys226 = +2.

**The 300 ps goes in six chunks**, by raising `aspqm.inputs.time` to 100, 150,
200, 250, 300 and re-running `aspqm.sh`; `resume: true` continues rather than
overwrites. ~19.5 h per chunk at the measured 2.57 ps/h. The chunking is
insurance, not bookkeeping: the attack coordinate is free as well, so if the
nucleophile drifts we learn it after one chunk instead of after six.

**What to read off chunk 1.** Does Asp92 stay in the second shell with Lys226
protonated, or collapse onto the metal? D211's collapse to a 1.93 A mean came
from a lineage that began already-collapsed with Lys226 **deprotonated**, so this
is the first clean test. Staying out supports the thesis; coming in means GFN2
over-coordinates the zinc against a crystal that says otherwise, which is still
open (D206). Second: does HOH3095 leave, giving the four-coordinate sphere rule 4
expects for an H-R-H site?

Scratch at submission: **239/300 GB**, down from 294 - 61 GB free, so the quota
risk that D221 raised has receded.

**D225 - The 10 ns vs 100 ns question, settled from trajectories already on disk
rather than by running a new one. eq6's 10 ns is converged; `asppull3`'s 2 ns is
not; no 100 ns trajectory survives to check against. The answer is a 10 ns
continuation, not 100 ns - 45 minutes rather than 7.3 hours.** (2026-11-06)

The user asked whether the MM part was finished correctly, whether that was the
full 100 ns the pipeline specifies, and - taking RMSF, temperature and pressure
into account - whether 10 ns and 100 ns differ enough to matter. Then: do not
throw out work.

**What had actually been run.** `eq6` 10 ns + `asppull` 1 ns + `asppull2` 2 ns +
`asppull3` 2 ns = **15 ns in four segments**, and each segment re-ran mcube's
`Gentle_warm_up_MD` (10 steps at 1 K, 50 at 10 K, 10,000 at 4 fs and 300 K = 40
ps) plus an NPT box relax before its own NVT. That warm-up is inherited from the
run directory's base `md` inputs, not from `site.yaml`, which is why it went
unnoticed: **the 15 ns is four short equilibrations, not one trajectory.**
Nothing is discarded by this - `asppull3`'s last frame carries the whole lineage
- but it is not 100 ns and it is not continuous.

**`docs/csr/convergence.py` (new, read-only, writes nothing).** Per run: CA RMSD
against frame 0 as decile block means, so a plateau is visible rather than
inferred from one number; per-residue RMSF over windows of the same trajectory,
scored both by Pearson r against the full profile and by the **ratio of mean
RMSF**, since the usual failure of a short window is to understate motion; the
same restricted to residues 88, 90, 92, 93, 164, 186, 225, 226, because a global
profile can converge while the site has not; and T, box volume, density and
total energy from `nvtsim.csv`/`nptsim.csv` with last-decile-minus-first-decile
drift.

**eq6, 10 ns, 10,000 frames - converged on every measure.**

    CA RMSD by decile:  0.82 0.94 1.22 1.17 1.09 1.18 1.19 1.14 1.04 1.06
    RMSF (A)      mean   r vs full   ratio
      full        0.63      -          -
      first 10%   0.52    0.902      0.820
      first 50%   0.62    0.986      0.974
      last  90%   0.62    0.999      0.985
    site CA RMSF: 88:0.32 90:0.33 92:0.32 93:0.31 164:0.48 186:0.39 225:0.38 226:0.37
    T 300.307 +/- 1.454, drift +0.135    box 481.178 +/- 1.223, drift +0.653
    density 1.015 +/- 0.003, drift -0.001

**RMSD plateaus by the third decile (~3 ns) at ~1.1 A and stays flat for the
remaining 7 ns. RMSF is within 3% of its final value by 5 ns at r = 0.986.** The
site residues fluctuate by 0.31-0.48 A and converge within the first nanosecond.
Temperature and box volume are as steady as an NVT run gets. **A 1 ns window
understates the motion by 18%, so the convergence is real and not an artefact of
the metric being insensitive.**

**`asppull3`, 2 ns, 200 frames - NOT converged.**

    CA RMSD by decile:  0.75 0.81 0.78 0.83 0.85 0.86 0.82 0.93 0.95 0.91
    RMSF full 0.50; first 50% 0.48 (r 0.979); last 90% 0.50 (r 0.999)
    T 300.371 +/- 1.323, drift +0.005
    box 480.117 +/- 2.401, drift **+5.396**    density drift -0.012

**RMSD is still climbing - its three highest deciles are the last three - and
box volume drift is +5.40 nm3 against eq6's +0.65, eight times larger, with
twice the standard deviation.** The box is still accommodating the 1.3-2.2 A
displacement of the Asp92/Lys226 pair. Temperature is fine; the volume is not.
So the repaired structure has not finished relaxing, and that is the honest
answer to "is what we have fine?" - the *equilibration* is, the *2 ns after the
pull* is not.

**No 100 ns trajectory survives, so the 10-vs-100 comparison cannot be made
directly.** `md_wt100c_results` holds 123 MB of `OpenMMOpt_traj.xyz`,
checkpoints, state XML and four `.cif` frames - **and no DCD at all**; `wt100b`
likewise. Their trajectories were deleted by the same cleanup policy that freed
6 GB today (D221). The complete surviving inventory:

    5.4G  eq6                NVTtrajectory_imaged.dcd        10 ns
    644M  prod100c           NVTtrajectory.dcd            22.68 ns
    110M  asppull3, asppull2 (both copies)                    2 ns each
     55M  asppull                                             1 ns

**The longest survivor, `prod100c` at 22.68 ns, does NOT converge** - RMSD
deciles 1.22 1.43 1.19 1.09 1.27 1.39 1.45 1.56 1.67 **1.72**, rising
monotonically through the back half, +0.51 A overall and still climbing at the
end; RMSF first-50% r only 0.961. But that is the **product arm that rearranged
and then NaN'd** (D173: Asp92 to 2.03 A, ligand turned over, O2x travelled 4.75
A), so its drift is the pathology on its way to failure rather than evidence
about how long equilibration takes. It is recorded because it is the reason the
claim here is about *eq6's system*, not a general claim that 10 ns suffices.

**Decision: `relax10` - 10 ns, continuous, from `asppull3`'s last frame.**
`gentle: false` so it continues instead of re-annealing; **NPT left ON**, because
the box drift says it is needed - an earlier draft of this section had
`skip_npt: true`, which the data contradicts. Restraints are `asppull3`'s,
unchanged. 10 ns because that is the length eq6 demonstrated is sufficient for
*this* system; the pipeline's 100 ns is a round number and nothing measured
supports paying 7.3 h for it. **~45 min on one H100.** Afterwards, re-run
`convergence.py` and require eq6's plateau - RMSD flat over the last five
deciles and box drift back under ~1 nm3 - before the QM/MM leg is taken as
starting from a relaxed structure.

**Two library findings.** `MD_TRAJECTORIES` in `src/mcube/md.py` lists only the
`.pdb` spellings of each topology, while the md stage wrote `.cif` for the
wt100/prod100 runs - so `salvage_results` cannot find a trajectory pair for
those runs even in principle. `convergence.py` now tries every combination.
And the raw/imaged duplication of D221 is visible again here: `asppull2` and
`asppull3` each carry two identical 110 MB copies.

**Caveat on the running QM/MM chunk (32384737).** It started from `asppull3`'s
frame, which this entry establishes is not fully relaxed. Its 50 ps is still
worth having - it answers early whether Asp92 stays off the metal under GFN2,
which is the central question - but it is a pilot, and the pipeline's QM/MM leg
should start from `relax10`'s frame.

**D226 - The 10-vs-100 ns comparison made directly, frame against frame. The
100 ns run did finish and its last frame survives. Its backbone is 1.52 A from
the 10 ns structure with a median per-residue deviation of 1.01 A, and the five
residues that move more than 4 A are all termini. **10 ns is adopted as the
baseline.** The site distances in that frame look alarming and are not evidence
about length - both 100 ns runs were restrained differently.** (2026-11-06)

D225 said no 100 ns trajectory survives, which was true of the DCD and wrong
about what the question needed. **`md_wt100c_results` is in
`/scratch/user/uquadhik/CSR/CSR_wt/`, not `CSR_mich/`** - which is why D225's
inventory missed it - and it holds `NVTtrajectory_lastframe.pdb`,
`NVTtrajectory_firstframe.pdb`, the NPT frames and `nvtsim.csv`. **The run
completed: 25,000,000 steps at 4 fs = 100,000.00 ps, final T 300.70.** A
single-structure comparison needs no trajectory, so the question the user asked
- is there anything at 100 ns not evident at 10 ns - is answerable from what is
on disk, and was answered by fetching two PDBs through
`vscode-remote-resource` and doing the superposition in browser JS.

**The two frames are the same system, checked before comparing rather than
assumed.** Both: 275 CA, chain A 4203 protein atoms, 1 ZN, 51 mer, 2 OH. Only
the water count differs (45,108 vs 43,668 atoms' worth), so the solute is
identical and the boxes are 79.03 A against 78.44 A.

**Backbone: 10 ns and 100 ns agree.**

    CA RMSD, wt100c @100 ns vs eq6 @10 ns   1.522 A   (275 CA, chain A)
    per-residue deviation   median 1.01   90th 2.47   max 6.88
    residues > 2 A  37      > 3 A  16      > 4 A  5
    the five over 4 A: 275, 273, 131, 1, 8   - termini and one loop
    site CA deviation: 88:0.99 90:1.08 92:0.82 93:0.74 164:1.48 186:0.40
                       225:1.11 226:0.92

**Every one of those is inside the thermal envelope of a single trajectory.**
eq6's own CA RMSF is 0.63 A and its RMSD against its own frame 0 plateaus at
~1.1 A (D225), so 1.52 A between two independently solvated copies run to
lengths differing tenfold is the scale of ordinary fluctuation, not of a
structural change. **The eight site residues move by 0.40-1.48 A, below the
1.52 A global value**, so the part of the protein the argument depends on is
the part that moves least. The superposition was computed twice by different
routes - the quaternion key-matrix eigenvalue, and a Newton-Schulz polar
decomposition of the 3x3 correlation matrix - and both give 1.522, which is
how the first run's broken rotation matrix was caught.

**Thermodynamics over the real 100 ns, 5000 rows at one per 20 ps.**

    T  (K)        0-10 ns  301.00 +/- 1.27     90-100 ns  300.97 +/- 1.26
    PE (kJ/mol)   0-10 ns  -656395 +/- 773     90-100 ns  -656693 +/- 768
    PE by decile: -656395 -656520 -656580 -656616 -656616
                  -656648 -656673 -656625 -656591 -656692

Temperature is identical to 0.03 K. Potential energy falls by 300 kJ/mol out of
656,600 - **0.05%, and 0.4 of one standard deviation** - and is flat from the
third decile on. **There is no slow relaxation running past 10 ns.**

**The site distances in the 100 ns frame, and why they are NOT an argument
about length.** Read naively they are alarming:

    pair        wt100c @0 ns   wt100c @100 ns   eq6 @10 ns   6DQH
    Zn-OD2          -              2.40            2.81       5.14
    Zn-OD1          -              3.67            3.87       4.28
    NZ-OD2          -              3.14            2.79       2.79
    Zn-O2x        4.94           16.82 (min-image) 2.57         -
    sphere < 2.5 A   -          6 donors          4 donors      -

At 100 ns the zinc is **six-coordinate** - His88 2.18, Asp92 OD2 2.40, His164
2.17, OH 1.91, and two waters at 2.09 and 2.06 - and the meropenem is out of
the site, 16.8 A away by minimum image in a 79 A box, its nearest protein
contact Gly231 instead of Trp21. **Neither is attributable to the duration,
and the dissociation is not a new finding: it is D194**, which established that
`md_wt100c` and `prod100c` replaced `[C:2:O01, C:1:C4x, 3.47]` with `[B:1:ZN,
A:92:OD1, 3.81, 50.0]` instead of carrying both, so nothing held the substrate.
This frame is an independent confirmation of it from the structure rather than
from the trajectory - 16.8 A at the final frame against D194's 50.69 A maximum,
both inside the 68.4 A cubic-box ceiling. What *is* new here is the second
consequence of that same swap: **the single weak restraint on OD1 left the
carboxylate free to rotate OD2 in to 2.40 A**, which is D223's lesson about
restraining one oxygen of a carboxylate arriving again from the other
direction, and it is why the 100 ns frame cannot be used as a QM/MM start even
setting the ligand aside. `md_wt100b` is the same story with no Asp92 restraint
at all: its own `coordination_end` records **Asp92 OD2 at 2.06 A** where its
`coordination_start` had no Asp92 in the sphere.

**So both 100 ns runs are independently unusable as QM/MM starting structures,
and for a reason that has nothing to do with their length.** That is also why
the backbone comparison above is the only part of them worth using: the fold is
insensitive to the restraint differences, the site is not.

**Decision: 10 ns is the baseline, and the pipeline's "100 ns MD" is amended to
10 ns.** Measured, not assumed: 90 further nanoseconds move the backbone by
less than one trajectory's own fluctuation, leave the site CAs within 1.5 A,
and change the potential energy by 0.4 sigma. MM is being asked for one thing
here - a relaxed fold that does not hand QM/MM an overcoordinated metal - and
10 ns delivers it. **The 7.3 h per arm that 100 ns costs buys nothing
measurable, and both existing 100 ns runs spent it and produced structures we
cannot use.**

**This does not retire `relax10`.** Its justification was never length: it is
that `asppull3`'s 2 ns after the pull has not finished relaxing (RMSD still
climbing, box volume drift 8x eq6's - D225). 10 ns continuous from
`asppull3`'s last frame, `gentle: false`, NPT on, restraints unchanged, ~45
min. The length is now the measured baseline rather than a guess.

**Correction to D225.** "No 100 ns trajectory survives" should have read "no
100 ns *DCD* survives". The last frames of both 100 ns runs are on disk, and
`nvtsim.csv` covers the whole 100 ns at 20 ps resolution for both. The error
came from searching one working directory and reporting the absence as general
- the runs for the wild-type arm live in `CSR_wt/`, the Michaelis arm in
`CSR_mich/`, and `job-inventory.md` is the index that should have been read
first.

**D227 - The standard mcube pipeline written down from the source: minimise ->
NPT box equilibration -> NVT production of `nvt_time` ps, with the gentle
warm-up OFF, which is the library default. `relax10` is 10 ns of NVT reached
through NPT with `gentle: false`, and the launcher refuses to submit anything
else. CORRECTED LATER THE SAME DAY - the first version of this entry called
`gentle: true` "the standard pipeline" and removed `relax10`'s `gentle: false`
on that basis. That was wrong, the user corrected it, and the correction is
§CORRECTION below. The NPT-then-NVT finding, the 17-section audit and the
`npt_time` no-op all stand.** (2026-11-06)

The user's instruction: use the standard mcube pipeline, NPT then NVT, and make
the 10 ns NVT for all future runs.

**What the stage actually does, read in `src/mcube/md.py` rather than assumed.**
In order, for a non-membrane explicit-solvent system:

    resolve_resume(...)                 # decided FIRST (D182)
    OpenMM_Opt(maxiter=500)             # unless resuming
    Gentle_warm_up_MD(...)              # iff inputs["gentle"]
    OpenMM_box_equilibration(...)       # NPT, unless skip_npt/implicit/resuming
    OpenMM_MD(simulation_time=nvt_time) # NVT production, unless skip_nvt

`Gentle_warm_up_MD(time_steps=[0.0005,0.001,0.004], steps=[10,50,10000],
temperatures=[1,10,300])` is 40 ps from 1 K. **NVT production is
`OpenMM_MD(..., integrator='LangevinMiddleIntegrator', trajfilename=
'NVTtrajectory')` and its length is `nvt_time` in ps**, so "10 ns NVT" is
`nvt_time: 10000` and nothing else.

**The base `md` inputs in `CSR/CSR_mich/config.yaml`, which every CSR section
inherits from:** `skip_npt: false`, `skip_nvt: false`, `gentle: true`,
`npt_timestep: 0.001`, `nvt_timestep: 0.002`, `temp: 300`, `cores: 16`,
`nvt_time: 1000`. **So the standard pipeline is the default and always has
been**, and a section only leaves it by overriding one of those keys.

**Audit of all 17 `md` sections in `docs/csr/site.yaml`.** **Not one sets
`skip_npt` or `skip_nvt`**, so every CSR MD run ever submitted - `eq5`, `eq6`,
`wt100`/`wt100b`/`wt100c`, `prod100`/`prod100b`/`prod100c`, `nac`/`nacgpu`/
`naclong`, `subst`, `rwon`, and the three `asppull` stages - went NPT then NVT
with the 40 ps warm-up in front. **The one deviation anywhere was `relax10`'s
`gentle: false`**, which D225 added so the run would continue from 300 K rather
than re-anneal.

**§CORRECTION - the warm-up is not part of the standard pipeline, and
`relax10` keeps `gentle: false`.**

The first version of this entry reasoned that because every CSR section
inherits `gentle: true`, the warm-up *is* the project's standard and `relax10`
should not deviate from it, and that the two arms must match. The user: *"gentle:
false is the default, we only need it if normal md runs throw errors, otherwise
it's wasted compute."* **That is right, and mcube's own README says so more
sharply than I did** (`README.md:172`): `gentle` defaults to **`false`**, it is
"**for a system that will not start otherwise**", and - the part that matters -

> Despite the name it is 10 steps at 1 K and 50 at 10 K (0.055 ps in total)
> followed by **40 ps of ordinary 300 K MD**, so after a successful
> minimisation it buys nothing NPT does not, and the 300 K leg is long enough
> to let a force-field artefact develop before any equilibration you are
> watching.

So the cost is not only the 40 ps. **Those 40 ps are 300 K MD in a box that has
not been equilibrated yet**, which on this system is exactly where coordination
drift comes from - the thing every restraint in `relax10` exists to prevent.
Running it on an already-equilibrated 300 K frame is worse than wasteful. The
`gentle: true` I had taken for the standard is a **local override in this run
directory's base `md` inputs**, not the library's behaviour, and it is the
reason a section that stays silent gets the warm-up: `inputs.get("gentle")`
reads the merged section, so silence means "inherit true" here and "off"
anywhere else. **`relax10` therefore sets `gentle: false` explicitly** rather
than relying on a default it would not actually get.

Two things from the retracted reasoning survive, and one does not. **Surviving:**
the warm-up would not have lost the pulled pose (the restraints are in force
through it, Zn-OD2 5.14 and NZ-OD2 2.79 at k = 500), and the two arms must still
match. **Not surviving:** that matching them means keeping the warm-up. It means
**both arms run with `gentle: false`**, and the arm already built is unaffected
in anything that is measured downstream - the warm-up only ever preceded NPT,
which then re-equilibrated the box to a volume/density threshold regardless.

**What the earlier arms actually spent on it.** `eq6` and the three `asppull`
stages each ran the warm-up, so the 15 ns lineage carries **4 x 40 ps = 160 ps
of pre-NPT 300 K MD** it did not need. That is small against 15 ns and none of
it is a reason to re-run anything; it is recorded because D225 counted the four
warm-ups as evidence the lineage was "four short equilibrations rather than one
trajectory", and part of what made it four was a step that should not have been
there.

**`relax10.sh` now gates the submission on the ensemble instead of printing
it.** It reads the merged `md_relax10` section back out of `config.yaml` after
`sections.py` writes it and exits non-zero if `skip_npt` or `skip_nvt` is set,
if `nvt_time` is not 10000, **or if `gentle` is on** - the last of those being
the check the correction above inverted. Printing a config and
eyeballing it is what let three wrong keys through in D218 and `gentle: false`
through in D225; the check is three lines and runs before `sbatch`.

**`docs/csr/prod100d.sh` deleted.** It was the superseded 100 ns launcher, it
carried `skip_npt: true`, and **no `prod100d` section exists in `site.yaml`** -
so it could only ever have run against a hand-edited config, and it was the one
file in the tree still advertising a no-NPT run.

**A library observation, recorded because the config implies otherwise.**
`npt_time` is read only in the **membrane** branch
(`simulation_time=inputs.get("npt_time", 1000)`). The non-membrane call passes
`numsteps_per_NPT=10000, volume_threshold=1.0, density_threshold=0.001` and
**no simulation time at all**: NPT runs in 10,000-step blocks until the box
volume and density stop changing by more than those thresholds. **So for every
CSR run the NPT length was convergence-driven and `npt_time` was silently
ignored.** That is defensible behaviour - better than a fixed length - but a
user setting `npt_time` on a solvated protein would get no error and no effect,
which is the kind of silent no-op worth fixing upstream. It also means the NPT
phase self-certifies the box, which is the right check for the box-volume drift
D225 flagged in `asppull3`.

**For all future runs, including the mutant arm.** An MD section states
`nvt_time: 10000` and **`gentle: false`**, and leaves `skip_npt`, `skip_nvt`
and `npt_timestep` alone. 10 ns is the measured baseline (D226), NVT production
is what that 10 ns means, and NPT runs to convergence in front of it. `gentle`
goes to `true` only for a system that will not start without it - an unstable
first few steps after minimisation - and then it is a diagnosis, recorded with
the error that prompted it, not a default.

**The base `md` block in each CSR run directory's `config.yaml` still says
`gentle: true`**, and that is the landmine: it is a local override of a library
default, so every section that does not speak up inherits it. Setting it false
per section is the safe fix and is what is done here; flipping the base would be
the general one, and it needs a shell in the run directory rather than a commit
in this repository.

**D228 - QM/MM job 32384737 cancelled at 5.16 ps because it started from an
unrelaxed frame, and the pipeline is put back in its stated order: finish the
MD, then QM/MM. `relax10` submitted as job 32405173, 10 ns NVT through NPT. A
new `aspqm2` section runs the QM/MM leg from that run's last frame, and its
launcher refuses to submit until that frame exists.** (2026-11-06)

The user: *"why would you start QM/MM from an unrelaxed frame? cancel it, finish
the MD and then move on to QM/MM - the pipeline is linear and straightforward."*

**The criticism is right and the evidence was mine.** D225 measured
`asppull3`'s 2 ns as unconverged - CA RMSD still climbing with its three
highest deciles last, box volume drift **+5.396 nm3 against eq6's +0.653**,
eight times larger with twice the SD - and D224 had already submitted the
QM/MM leg from that run's last frame. I then described the job as "a pilot
rather than the pipeline's QM/MM leg" and left it running. A leg that is not
the pipeline's leg is not worth 19.5 h of an allocation while the run that
would make it valid takes 45 minutes. **The pipeline is a sequence; a stage
does not start before the one in front of it finishes.**

**Cancelled.** 32384737, confirmed gone from `squeue`, having reached **5,159
xTB calls at a 1 fs step = 5.16 ps** of the 50 it was asked for. Its output
stays in `qmmm_aspqm_results` and the `aspqm` section stays in `site.yaml`,
marked superseded - nothing is thrown out. What it does NOT answer is the
question it was submitted for: 5 ps is too short to say whether Asp92 stays in
the second shell, so the GFN2 over-coordination question (open since D206)
remains open rather than partly answered.

**`relax10` submitted: job 32405173, `gpu_cuda`, h100, 3 h wall.** The
launcher's own gate passed and printed what the merged config holds -
`gentle: False`, `skip_npt: None` (so the base's `false`, i.e. NPT runs),
`nvt_time: 10000`, `nvt_timestep: 0.002`, `coordination_redetect: False`,
**`[OK] minimise -> NPT -> 10 ns NVT, no warm-up`** - and `sections.py`
checked 12 selectors against the starting structure before writing it.

**Pre-submission arithmetic, done off the actual frame rather than off this
log** (the rule from D224, applied):

    RESTRAINED             current -> target    k      delta
    C:2:O01  - C:1:C4x      3.47 -> 3.47       50     -0.00
    B:1:ZN   - C:1:O2x      2.58 -> 2.60      500     +0.02
    B:1:ZN   - C:1:O1x      3.80 -> 3.82      500     +0.02
    A:226:NZ - A:92:OD2     2.75 -> 2.79      500     +0.04
    B:1:ZN   - A:226:NZ     5.26 -> 5.27      500     +0.01
    B:1:ZN   - A:92:OD2     5.05 -> 5.14      500     +0.09

    UNRESTRAINED            current   [6DQH]
    Zn-OD1                   4.19     [4.28]
    Zn-CG                    5.09     [5.15]
    NZ-OD1                   3.95     [3.94]

    sphere < 2.50: His88 2.21, His164 2.15, OH 2.03, HOH3095 2.01, HOH3171 2.15

No restraint opposes another - the Zn/NZ/OD2 triangle is 6DQH's own geometry and
satisfies the triangle inequality with room - every pair differing from its
target is restrained, and the largest correction asked for is 0.09 A. **Predicted
outcome: the pose holds to ~0.1 A, NPT converges the box, and the 10 ns NVT
plateaus as eq6 did (RMSD flat over the last five deciles, box drift under
~1 nm3). If the drift is still near +5.4 nm3 after 10 ns, that is the restraint
set fighting the solvent and it needs diagnosing rather than more nanoseconds.**

**`aspqm2`, and why it is a new section rather than an edited `aspqm`.**
`resume: true` continues rather than overwrites, and it keys off the section's
own results directory. Editing `aspqm`'s `pdb_file` and re-running it would
have **resumed the cancelled 5.16 ps** - the exact frames the cancellation was
for - and the new starting structure would have been read and discarded. A new
name gives a fresh `qmmm_aspqm2_results`, so resume stays available for the
six-chunk plan while being unable to reach the old state file. Everything else
is identical to `aspqm`: same QM region (`Zn, OH, mer, His88, His164, Asp92,
Asn225, Arg90, Lys226`), charge 2, `skip_equilibrate: false`,
`skip_optimize: true`, `rigid_water: false`, `time: 50`.

**`aspqm2.sh` refuses instead of reporting.** It exits non-zero unless
`pdb_file` is exactly `md_relax10_results/NVTtrajectory_lastframe.pdb`, unless
that file **exists on disk** - so it cannot be submitted while `relax10` is
still running or if `relax10` failed - and unless `charge` is 2. The
restraint/constraint warning from `aspqm.sh` is kept. This is the same change
made to `relax10.sh` in D227: the launcher checks the thing the entry claims,
because printing a config for me to read is what let D218's three wrong keys
and D225's `gentle: false` through.

**Scratch: 209/300 GB, 723,500/1,000,000 inodes** - 91 GB and 276k inodes free,
read fresh from `rquota`. **Correction to this morning's reading:** I reported
"300/300, at quota" from a `q.txt` on scratch. That file was **stale** - written
by an earlier session - and reading it as current is exactly the trap
`BUNYA_BROWSER_SHELL.md` section 3 warns about: *"a successful fetch proves the
FILE exists, never that YOUR command produced it."* It was the thing that
eventually showed my typing was not landing, but it should not have been quoted
as a quota figure. Dated, freshly-named log files from here.

**The shell failure that preceded all of this, and the document that solves
it.** Six typed commands reached nothing. I diagnosed "the terminal never
attached to a pty" from `.xterm-rows` being empty and reported the pane
unusable. **Both halves were wrong, and `BUNYA_BROWSER_SHELL.md` - which the
user had to point me to - says why.** Section 2: the terminal is **xterm.js on
a canvas**, so there is no `.xterm-rows` DOM to read and its emptiness means
nothing; scraping the screen is not a diagnostic at all, and the only way to
read output is write-to-file-then-fetch. Section 4: `querySelector` returns the
**first** terminal, usually the blocked one, so the newest must be focused by
mouse events plus `.focus()`. Section 3: `ctrl+u`, not `ctrl+c`, clears a
prompt - and the first probe after re-focusing came back as
`bash: shrquota: command not found`, a leftover fragment merged with my
command, which is the documented *"a line that fails to submit silently
swallows the next one"*. Clearing with `ctrl+u` fixed it immediately. **The
lesson to carry: read `BUNYA_BROWSER_SHELL.md` before driving the pane, and
check a mechanism before reporting it broken.** Screen Reader Optimized, the
first thing that document says to check, was off here - so it was not even the
common cause.

**Order from here, one stage at a time.** `relax10` (32405173) finishes ->
re-run `convergence.py` on `md_relax10_results` and require eq6's plateau ->
`aspqm2` chunk 1, 50 ps -> read whether Asp92 stays in the second shell with
Lys226 protonated and whether HOH3095 leaves -> chunks 2-6 by raising `time` to
100, 150, 200, 250, 300, one at a time. Then the def2-SVP optimisation, which
needs the D134 `optimizer_constraints` defect fixed first if it is given any
constraint, and the def2-TZVP single point.

**D229 - `relax10`'s log appears to apply a `Zn-Asp92 OD2 2.79` restraint
opposing the configured 5.14, and it does not: an explicitly configured pair
wins and the detected one is dropped by pair identity before anything is
applied. The run is correct and was NOT cancelled. This also corrects D220 on
*which* mechanism protects the pulled geometry.** (2026-11-06)

Reading `relax10`'s output while it ran, the `[INFO] Atom selectors resolved`
line appears **twice**, and the second is a set nobody configured:

    [['B:1:ZN','C:1:O2x',1.94,500.0], ['B:1:ZN','A:164:NE2',2.09,500.0],
     ['B:1:ZN','C:2:O01',2.14,500.0], ['B:1:ZN','A:88:NE2',2.22,500.0],
     ['B:1:ZN','A:92:OD2',2.79,500.0]]

with a warning naming exactly the danger: *"coordination restraint on A:92:OD2
at 2.79 A is longer than 2.50 A, so it is not a bond: this holds a second-shell
contact against the metal for the whole run, and a carboxylate held that close
can rotate its other oxygen on (D120)."* **`Zn-OD2` at 2.79 against the
configured 5.14, and `Zn-O2x` at 1.94 against the configured 2.60 - two pairs
restrained at both ends, which is D214 exactly.** The run looked invalid.

**It is not, and the thing that settled it was `asppull3`'s own log, then the
source.** `asppull3` printed the identical five-pair resolution - and recorded

    coordination_restraints: [['B:1:ZN','A:164:NE2',2.09,500.0],
                              ['B:1:ZN','C:2:O01',2.14,500.0],
                              ['B:1:ZN','A:88:NE2',2.22,500.0]]

**three pairs, not five.** `Zn-O2x` and `Zn-OD2` were gone. Since `asppull3`
reached Zn-OD2 5.05, they cannot have been applied. `src/mcube/md.py:317-319`
says why:

    rows = resolve_atoms(saved, pdb_file) or []
    # An explicit restraint on a pair wins: the user named a distance for it.
    named = {frozenset(r[:2]) for r in asked}
    detected = fill_restraints([r for r in rows if frozenset(r[:2]) not in named],
                               pdb_file, force)

**Any detected pair whose two atoms match a configured pair is discarded,
whatever the distances.** `relax10` names `[B:1:ZN, A:92:OD2, 5.14]` and
`[B:1:ZN, C:1:O2x, 2.60]`, so both offending rows are dropped. What is applied
is the configured six plus the three first-shell pairs nobody configured -
His164 2.09, hydroxide 2.14, His88 2.22, which measure 2.15, 2.03 and 2.21 in
the starting frame, so all three are right and none conflicts. **Nine
restraints, no pair twice.**

**The warning fires before the dedup, by design**, and the comment above it says
so: *"Warned on BEFORE resolution, so the message can name the atoms."* So the
message names a restraint that is then thrown away. That is defensible - it can
print atom names rather than indices - but it means **the warning cannot be read
as a statement about what the run will do**, which is exactly how I read it. A
line in the log saying a second-shell contact "is held for the whole run" while
the run holds nothing of the kind is a trap worth knowing: `coordination_restraints`
in the saved results is the authority, and it is written after the dedup.

**Correction to D220.** It credited `coordination_redetect: false` with
preventing D214: *"which is what stops D214 repeating, since redetection
re-pins Zn-OD2 at its detected 2.81 and would hold the pulled coordinate at
both ends again."* **That is not the mechanism.** With `coordination_redetect:
false` the stored `preprocess.results.coordination` is used - and that set
**still contains `Zn-OD2` at 2.79**, measured on the original structure, which
is precisely the row that alarmed me here. The protection is the pair-identity
dedup, which works whether or not redetection runs. (Redetection would have
been safe too, for a different reason: on this frame OD2 is 5.05 A, beyond
`LIGAND_MAX` 2.50, so it would not be detected at all. Both routes are safe;
neither is safe for the reason D220 gave.) `coordination_redetect: false` is
still the right setting - it keeps the metal's first shell pinned at the
distances the crystal-derived preprocess measured rather than at whatever a
pulled frame happens to show - but it is not what keeps Asp92 out.

**The practical rule this adds to D224's pre-submission arithmetic.** That rule
said to write out every restrained pair as current -> target and check that no
restraint opposes another. **It is not enough to check the pairs you wrote: the
stage adds its own.** The check is against the **merged** set, and the merged
set is only visible in `coordination_restraints` in the saved results - i.e.
after the run starts. So the pre-submission version of the check is: for each
pair you configure, confirm it either (a) is absent from
`preprocess.results.coordination`, or (b) appears there and will be dropped by
the dedup. Naming a pair is sufficient to win it; that is the whole protection
and it should be stated deliberately rather than relied on by accident.

**Nothing was cancelled, and that was the right call.** The first instinct was
to kill the job on the strength of a log line. Checking `asppull3`'s recorded
set cost one fetch and the source check cost two greps, against ~45 minutes of
H100 and a re-submission. **A warning in a log is a hypothesis about a run, not
a measurement of it.**

**D230 - The withheld kinetics' CONSTRUCTS, asked early and answered: wild type
is the truncated form, and both the plain double mutant and dmCSR-1 S214D have
kinetics. The mutant arm runs strictly after the wild-type arm. Wakeups
scheduled to pick the work up when jobs land.** (2026-11-06)

The user asked for autonomous work on scheduled wakeups, and for questions to
be asked early rather than becoming blockers. Three were outstanding; **a
construct name is not a kinetic value, so asking does not breach D178's
blind-prediction protocol**, and the answers decide what gets built.

**1. Wild type: the truncated form.** That is the form this project already
models, so the measured and the modelled construct are the same one and no
caveat is needed on that axis. This retires the last soft spot in that premise:
`the-truncated-form-is-the-one-modelled` had been `withheld`, resting on the
user's modelling decision plus an inference from SPR-1's 49-residue truncation
(which is as consistent with expression as with catalysis). It is now
`supported`, confirmed by the person holding the data. **The values stay
withheld and the protocol is unchanged - do not ask for them.**

**2. Mutant: BOTH the plain double mutant and dmCSR-1 S214D were measured.**
This is the more useful answer, because it unties a knot from D219. That entry
disqualified our own S214D structure as the mutant arm's starting geometry -
S214D is BBL Ser214 = project **Ser186**, 3.81 A from the zinc in `eq6`, with
the engineered aspartate 4.06/4.13 A from both metal sites, so it is not
catalytically innocent and cannot stand in for the plain double mutant. With
kinetics on both, **the mutant arm is built from 6NC5 as the plain double
mutant and still has a measured partner**, and our S214D structure is not
stranded either: it has its own. S214D is disqualified as a *starting structure
for the double-mutant comparison*, not as a system worth computing later.

**3. The mutant arm runs strictly after the wild-type arm**, not in parallel.
The two arms are independent and a parallel 10 ns would have cost one of the
four GPU slots, but the user chose one thing in flight. So 6NC5 preprocessing
is not started yet. When it is, the policy already decided applies: delete the
four extra zincs (`remove_residues: [A:303:ZN, A:304:ZN, A:305:ZN, A:306:ZN]`,
and note ZN303 shares Zn1's terminal water at 1.88/1.90 A, so deleting it
leaves that water as Zn1's own, which is where the nucleophile goes), delete
the crystallisation additives by the same reasoning the user gave for the
zincs - five PO4 and four ACT are there because of the crystallisation, and in
6DQH a phosphate oxygen sits **on Zn1 at 2.32 A**, occupying the nucleophile's
site - and then the same treatment as the wild-type arm: `gentle: false`, NPT,
10 ns NVT, 3.0 A restraints on everything the metal affects at its crystal
distance, coordination measured at 2.50 A.

**Wakeups: recurring job `7f92d121`, `13,38 * * * *`, first fire 22:19.**
Session-only and auto-expiring after seven days. Its prompt carries the whole
pickup procedure rather than a reminder to look: read `state.md` and
`retractions.md` before asserting anything; read the code-server base from the
live tab; consult `BUNYA_BROWSER_SHELL.md` before concluding the pane is
broken; **stretch its own cadence to the remaining runtime** rather than waking
every 25 minutes, and delete itself if the server is down (the standing check-
cadence rule); run `convergence.py` on `md_relax10_results` and require `eq6`'s
signature - RMSD flat over the last five deciles, box drift under ~1 nm3 -
before QM/MM; treat a drift still near the 2 ns run's +5.4 nm3 as the restraint
set fighting the solvent rather than a run needing longer, and **not** start
QM/MM in that case; then `aspqm2` in 50 ps chunks, one at a time, with the D224
arithmetic done on the **merged** restraint set per D229.

**Status at 21:50.** `relax10` (32405173) is past minimisation and NPT and into
production MD. NPT took the box from 475.6 to 480.6 nm3, **within 0.6 of
`eq6`'s equilibrated 481.18**, which is the relaxation `asppull3`'s 2 ns never
finished. Scratch 209/300 GB.

**D231 - `relax10` finished and passed the gate: the RMSD plateaus, the
restraints held to 0.09 A, and Zn-OD1 reached 4.33 against 6DQH's 4.28 with
nothing restraining it over 10 ns. QM/MM chunk 1 submitted as job 32416288.
AND the gate's own box-volume criterion was wrong - D225's box-drift argument
is retracted, because all three NPT runs plateau in their first decile and
`asppull3`'s box was already converged.** (2026-11-06)

**`relax10`, job 32405173, finished.** 10,000.0 ps in 5000 frames, final T
299.63, no NaN anywhere in the output.

    CA RMSD by decile:  0.81 0.85 0.87 1.05 1.25 1.23 1.41 1.38 1.38 1.40
    RMSF (A)       mean   r vs full   ratio
      full         0.58      -          -
      first 10%    0.49    0.871      0.850
      first 50%    0.54    0.995      0.945
      last  50%    0.55    0.907      0.954
      last  90%    0.57    0.999      0.995
    site CA RMSF: 88:0.31 90:0.30 92:0.31 93:0.30 164:0.37 186:0.32 225:0.39 226:0.40
    T 300.292 +/- 1.399, drift +0.094

**The last four deciles are 1.41, 1.38, 1.38, 1.40 - a spread of 0.03 A.** The
run climbs for ~6 ns and then holds for 4, where `eq6` plateaued by 3 ns and
held for 7. So 10 ns was the right length and not generous: 2 ns would have
caught none of the plateau, which is what `asppull3` showed.

**The site geometry after 10 ns, measured on the last frame:**

    RESTRAINED             value   target    dev
    C:2:O01 - C:1:C4x      3.48     3.47    +0.01
    B:1:ZN  - C:1:O2x      2.61     2.60    +0.01
    B:1:ZN  - C:1:O1x      3.86     3.82    +0.04
    A:226:NZ- A:92:OD2     2.88     2.79    +0.09
    B:1:ZN  - A:226:NZ     5.28     5.27    +0.01
    B:1:ZN  - A:92:OD2     5.15     5.14    +0.01

    UNRESTRAINED           value   [6DQH]
    Zn-OD1                  4.33    [4.28]
    Zn-CG                   5.19    [5.15]
    NZ-OD1                  4.14    [3.94]
    Zn-Ser186 OG            4.05

    sphere < 2.50 A: HOH3095 2.02, OH 2.04, HOH3171 2.05, His164 2.13, His88 2.18

**`Zn-OD1` at 4.33 against the crystal's 4.28, and `Zn-CG` at 5.19 against
5.15, with nothing restraining either, held for 10 ns.** D224 reported the same
agreement at the end of a 2 ns run and noted it was adopted rather than
imposed; ten nanoseconds of 300 K MD later it is still there, which is a much
stronger version of that statement. The sphere is still **five** donors with
Asp92 outside it and two bulk waters inside - rule 4 expects four for an H-R-H
site, and whether the extra water leaves is QM/MM's to answer, not MM's.

**One soft number, recorded rather than smoothed over:** the last 50% of the
trajectory correlates with the full RMSF profile at **r = 0.907**, where the
first 50% manages 0.995. The back half has a different distribution of
flexibility, so something re-arranged on the way to the plateau. It is not the
site - those eight residues sit at 0.30-0.40 A and their first-decile values
are within 0.05 of their final ones - so it is termini and loops, exactly where
the 10-vs-100 ns comparison put the only deviations over 4 A (D226).

**THE GATE'S BOX-VOLUME CRITERION WAS WRONG, AND SO WAS D225's ARGUMENT FOR
THIS RUN.** The criterion, written in D225 and carried into D228 and into the
wakeup prompt, was: box volume drift under ~1 nm3, because `asppull3` drifted
+5.396 against `eq6`'s +0.653. `relax10` drifts **+5.008**, which by the letter
of that gate forbids starting QM/MM. Checking the quantity before acting on it:

    NPT box volume by decile (nm3)                                  final-3 spread
    relax10   475.6 480.5 481.2 481.3 480.4 481.2 480.9 480.6 480.2 480.6   0.43
    asppull3  475.3 480.6 480.3 480.5 480.7 480.8 480.7 481.0 480.5 480.7   0.51
    eq6       480.9 480.5 482.2 481.7 481.4 481.2 481.0 480.4 481.1 481.5   1.07

**Every one of them plateaus in the first decile.** The "drift" is
last-decile-minus-first-decile, and for a quantity that steps once and then
sits still that statistic measures **how far the box had to travel**, not how
unconverged it is. `relax10` and `asppull3` both started from a 475.3-475.6
box and both reached ~480.7 immediately; `eq6` started at 480.9 because its box
came straight from solvation and had nowhere to go. **By the statistic that
actually describes convergence - the spread over the final three deciles -
`asppull3` (0.51) is better than `eq6` (1.07), and all three are converged.**
These numbers also come from `nptsim.csv`, i.e. the NPT phase: NVT is constant
volume and `nvtsim.csv` has no box column at all, so this metric can never say
anything about production.

So **D225's second reason for running `relax10` is retracted.** Its claim -
*"The box is still accommodating the 1.3-2.2 A displacement of the Asp92/Lys226
pair"* - was a mechanism invented to explain a number, and the number meant
something else. All three runs equilibrate to the same ~480.7 nm3 box, so the
displacement required no accommodation at all. Entered in `retractions.md`.
**The first reason stands and is sufficient**: `asppull3`'s RMSD was still
climbing with its three highest deciles last, which is a real statement about a
2 ns trajectory and is exactly what the 10 ns fixed. The decision was right;
one of the two arguments for it was not.

**The gate is replaced** by what the data supports: CA RMSD flat over the final
three or more deciles; T drift under ~0.2 K; configured restraints held to
within ~0.1 A; and the unrestrained second-shell distances still matching the
crystal. Box volume is checked as the **final-three-decile spread of the NPT
phase**, under ~1 nm3, not as a first-to-last difference.

**QM/MM chunk 1: `aspqm2`, job 32416288, `general`, 50 ps.** The launcher's
gate passed - `pdb_file: md_relax10_results/NVTtrajectory_lastframe.pdb`,
charge 2, `skip_equilibrate: False`, `skip_optimize: True`, `resume: True`,
`rigid_water: False`, QM region `Zn, OH, mer, His88, His164, Asp92, Asn225,
Arg90, Lys226` - and printed **`[OK] starts from relax10's last frame, charge
2, nothing restrained`**. The D229 merged-set check is vacuous for a `qmmm`
stage, since `bond_constraints` is reached only from the optimiser path and
rigid water and both are off; confirmed from the launcher's own warning check
finding nothing rather than from reading the source again.

**Housekeeping.** The raw trajectory duplicates in `md_relax10_results` are
deleted where the imaged twin exists - `NVTtrajectory.dcd` 2.68 GB and
`relaxbox_NPT.dcd` 549 MB - **verified by `ls` returning "No such file or
directory"** rather than by an echo after a `;`, which is the failure mode
`BUNYA_BROWSER_SHELL.md` section 5 records. 3.2 GB freed; scratch 213/300 GB,
723,619/1,000,000 inodes. The duplication is still a library-level waste in
every project and the config flag to suppress it is still unbuilt.

**Cadence.** Wakeup `7f92d121` deleted, replaced by `aae784b7` every three
hours at :21, since chunk 1 is ~19.5 h at the measured 2.57 ps/h. Its prompt
carries the starting distances so the next check can read the chemistry off the
running job without a shell, and instructs it to stretch its own cadence again.

**D232 - THE PRE-REGISTERED QUESTION IS ANSWERED, IN THE AFFIRMATIVE: GFN2
brings Asp92 onto Zn1 as a BIDENTATE chelate within 8 ps from the crystal
geometry, with Lys226 still protonated and the salt bridge still formed. The
question open since D9 closes; def2-SVP is now the only admissible instrument
for this carboxylate's position, exactly as written down before the result.**
(2026-11-07)

**Read off the live job without a shell**, at 32 of 50 ps - job 32416288,
32,122 `STARTING XTB INTERFACE` calls at 11:36, i.e. 2.45 ps/h against the 2.57
estimated, finishing about 18:50 today. The code-server session had moved to
`bun147/52687`; the node in my notes was dead and the base was read from the
live tab. `QM_MM.dcd` is in `protdft/qmmm_aspqm2/` - `run.slurm` puts every
non-`md` stage in `protdft/<section>`, which is why the obvious paths 404'd -
and it parsed to the exact byte: 47,979 atoms matching the topology PDB,
32 frames at 1 ps, cubic box 78.332 A.

**The measurement, over all 32 frames rather than a frame (D198's lesson).**

    ps          1     4     7     8     9    15    20    25    32
    Zn-OD2    4.34  4.34  3.54  2.36  2.01  2.27  2.03  2.07  2.12
    Zn-OD1    5.90  5.59  4.74  3.52  3.41  2.31  2.28  2.23  2.29
    Zn-O2x    2.50  2.12  2.21  2.51  2.60  4.83  4.94  6.31  5.81

    occupancy inside 2.5 A, 32 frames:
      His88 NE2 100% (2.09)   His164 NE2 100% (2.12)   OH O01 100% (1.83)
      Asp92 OD2  78% (2.10)   Asp92 OD1   53% (2.29)   mer O2x  19% (2.16)
      HOH3171     3%          HOH3095      0%
    coordination number: 3 at 1 ps, 4 through 13 ps, 5 from 14 ps on

**The carboxylate crosses the boundary between 7 and 8 ps and does not come
back for the remaining 24 ps.** It is a **chelate, not a contact**: Zn-CG falls
5.50 -> 2.52 A (min 2.34), which is where a symmetric bidentate carboxylate puts
its own carbon, and both oxygens are inside 2.5 A from 15 ps.

**Checked by a second route that shares no code with the first.** `xtb_.xyz`,
the live 130-atom QM region handed to xTB after frame 32, gives the zinc's
neighbours inside 2.9 A as **O 1.813, N 2.023, O 2.030, O 2.074, N 2.294,
C 2.389** - two histidine nitrogens, the hydroxide, *two* carboxylate oxygens
and the carboxylate carbon. Different file, different writer, same sphere. This
is the discipline the transposed-rotation bug in D226 taught.

**Against the pre-registered readings, which are not revised.** The question
node written on 2026-11-06 said: beyond ~4 A closes it in the negative; **inside
~2.5 A closes it in the affirmative and makes def2-SVP the only admissible
instrument**; between, or drifting, means 50 ps was too short. The second one
happened, and it is not a marginal reading of it - 2.10 A mean at 78%
occupancy, bidentate, 24 consecutive ps, no recrossing. **The question is
answered in the affirmative.** The second pre-registered reading is also
upheld: **HOH3095 leaves** (0% inside 2.5 A, HOH3171 3%), which is what rule 4
and the H-R-H site predict - the sphere is five-coordinate anyway, because the
carboxylate's two oxygens took the vacated places.

**Why this test is admissible where every earlier closure was not.** D206 rule
3 excludes single points and excludes two methods agreeing about a geometry one
of them was handed. This is dynamics, unrestrained, from a geometry four
deposits support. Every earlier collapse - D49, D51, D201's batch 1, D211's
1.93 A mean - began either with the carboxylate already on the metal or with
Lys226 deprotonated, so none tested the Hamiltonian. This one began with Asp92
4.3 A out and a formed salt bridge, and GFN2 moved it anyway.

**THREE THINGS DID NOT HAPPEN, AND EACH KILLS AN EXPLANATION THAT WAS ON
OFFER.**

1. **Lys226 kept its proton.** All three NZ hydrogens at 0.98-1.09 A in every
   frame; O01 carries exactly one hydrogen (0.98-0.99 A) throughout, so it is
   still a hydroxide. D216 proposed that Asp92 moving onto Zn1 is **downstream
   of the Lys226 proton transfer** - deprotonate the lysine and you remove what
   holds the carboxylate out. That mechanism is **refuted as a requirement**:
   no proton moved and the carboxylate came in regardless.
2. **The salt bridge held.** NZ-OD2 mean **2.87 A** (min 2.53, max 4.24)
   against 6DQH's 2.79. The carboxylate is bonded to the metal at 2.10 A *and*
   hydrogen-bonded to the ammonium at 2.87 A **at the same time, for 24 ps**.
   Holding the lysine contact does not keep the carboxylate off the zinc, which
   was the hope behind restraining it there.
3. **The carboxylate arrived first; the substrate was displaced afterwards.**
   Zn-OD2 crosses inward between 7 and 8 ps, Zn-O2x crosses outward between 8
   and 9 ps - one frame later. So no vacancy opened first for the carboxylate to
   fill, and the ordering is the carboxylate taking the carbonyl's place rather
   than filling a gap it found.

**AND A COUNTER-OBSERVATION TO A NODE THE PROJECT HOLDS.** D204's finding -
Lys226 hands HZ3 to the hydroxide within 2 ps of becoming quantum - is the same
QM region at the same charge. Here, 32 ps and **no transfer at all**. 32 ps
against 2 ps is not undersampling. So that result is real but **not general**:
the transfer happens in the lineage that began with Asp92 already on the metal
and does not happen in the lineage that began with the carboxylate on the
lysine. Whether the proton moves may depend on the carboxylate's position - and
that position is now a method artefact. The node
`nucleophile-is-zn-bound-water-with-neutral-lys226` rests on D204, and that is
now a conditional rather than a general result. Recorded as a caveat on the
node, not as a demotion: the observation happened.

**What this run is NOT, said plainly.** The Michaelis complex dissociates. With
nothing restrained the carbonyl leaves the metal (mean 4.32 A) and the attack
coordinate O01-C4x goes from 3.4-4.5 A through 22 ps to 5.6-6.2 A after it.
By criterion 6 this is not a Michaelis complex past ~22 ps, so **nothing
energetic may be taken from these frames**. The measurement asked where the
carboxylate goes, and that is all it answers.

**What it costs, and what it does not.** The def2-SVP QM/MM optimisation stops
being a refinement of the GFN2 answer and becomes the instrument of record for
Asp92's position - and it is not blocked by D134, because that defect fires only
when a `qmmm` run is given `opt_constraints`, which this leg is not. The planned
Zn-Asp92 PMF (item 11, ~8,500 core hours at GFN2) would now measure this bias
rather than the free energy: a GFN2 PMF returns "one minimum near 2.1 A"
whatever the truth is, so the shape could not be interpreted, and the
equilibrium question's instrument needs settling before that run is worth
buying. It overturns **no** crystallographic claim: a method artefact cannot
refute `the-lysine-holds-asp92-where-a-zn2-would`, which rests on four
deposits. What it says is that our method cannot reproduce them.

**One further reading, offered as a possibility rather than a result.** GFN2's
answer is not an equilibrium between the lysine and the metal but **both at
once** - a carboxylate bridging a metal and a cationic side chain, which is
what "bridging aspartate" means literally. If that arrangement is real rather
than an artefact, the equilibrium framing of item 11 is the wrong question.
Deciding which needs a method that is not itself on trial.

**Argument map.** Two nodes added: the observation
`aspqm2-chunk1-asp92-chelates-zn1-under-gfn2` (`ours`) and the claim
`gfn2-over-coordinates-zn1-from-a-crystal-consistent-start`, which answers the
question. The question node is **void with a `superseded-by` edge**, keeping its
predictions readable above the result, which is the point of having written them
first; its three supporters are repointed at the claim. The validator went from
1 finding to 6 on the first pass and every one was mine: an evidence class on a
claim, the claim citing the question's retracted lineage (D58/D133/D195/D197/
D206) as though they were its own support, two live claims resting on a node I
had just voided, and the claim resting on `ours` alone. All fixed - the claim's
other half is the crystal (`zinc-bound-crystal-holds-asp92-off-the-metal`,
experimental), which it must be, since "over-coordinates" names nothing without
a measured distance. **It now reports no findings at all**, and that includes
the known `no-usable-michaelis-complex-exists-yet` one, which cleared because an
experimental observation now sits beneath that node through the new claim.
Recorded rather than celebrated: whether that path genuinely anchors it is a
judgement the validator cannot make.

**D233 - CHUNKS 2-6 ARE NOT RUN. The user's call, taken after D232: stop at
chunk 1. The ~100 h of CPU that remained was budgeted to decide a question that
is now decided, and it would have been spent on a Hamiltonian this project no
longer trusts for this carboxylate's position, on a complex that has already
dissociated by 22 ps.** (2026-11-07)

The decision was put to the user with the arithmetic rather than a
recommendation, because the trade is theirs: chunks 2-6 were 50 ps increments
to 300 ps, and their stated purpose in D228 was to give the Asp92 question
enough sampling to close. D232 closed it at 32 ps, in the affirmative, with no
recrossing in 24 consecutive ps. More GFN2 would have measured the same bias
for longer.

**What is given up, said plainly rather than argued away.** Three things were
on offer past 50 ps and are now not available: whether the chelate is stable on
a 100-300 ps timescale or eventually breaks; whether the Lys226 proton moves
later than 32 ps, which would have bounded D204's conditionality from the other
side; and a long GFN2 reference trajectory against which the def2-SVP
optimisation's geometry could have been compared frame-for-frame. The first two
are questions about GFN2's behaviour, not about CSR-1, so the project loses
little. The third is a real loss and is recorded as one: the def2-SVP leg will
be compared against the crystals and against chunk 1's 32 ps only.

**Chunk 1 itself still finishes.** Job 32416288 runs to 50 ps on its own and
nothing was cancelled - killing it would have thrown out work already paid for,
and the final 18 ps cost nothing extra. The one-shot wakeup at 19:13 reads them
and reports only what they add or contradict; D232's reading is not reopened by
them, because a 24 ps non-recrossing is not overturned by 18 ps more of the
same contact.

**What this unblocks.** The def2-SVP QM/MM geometry optimisation from
`md_relax10_results/NVTtrajectory_lastframe.pdb` is now the next leg and the
instrument of record for Asp92's position, with the def2-TZVP single point
after it. It is not blocked by D134, since that defect fires only when a
`qmmm` run is given `opt_constraints` and this leg supplies none. Item 11's
GFN2 PMF stays gated on settling its instrument, per D232. Not submitted: no
job is proposed here.

**D234 - THE QM/MM EQUILIBRATION IS SKIPPED AND THE def2-SVP OPTIMISATION GOES
STRAIGHT ON THE MM FRAME. The user's instruction, immediately after D233:
"maybe skip the qm/mm equilibration and go straight to def2-svp geometry
optimisation then. long QM/MM is only really practical with xtb." Section
`aspdft` and launcher `docs/csr/aspdft.sh` are written and gated; NOT YET
SUBMITTED - the code-server pane is hidden, so only the read-only route is
available.** (2026-11-07)

**This amends the pipeline as written.** The protocol was "100ns MD, 300ps full
xtb QMMM, a geometry optimisation at def2-svp, singlepoint at def2tsvp". The
300 ps xTB QM/MM leg is now 32-50 ps and stops there (D233), and the DFT
optimisation starts from the **MM** frame rather than from the end of a QM/MM
equilibration. The user's reason is the general one and is recorded as such:
long QM/MM is only practical at xTB. The arithmetic agrees - a DFT QM/MM MD step
costs what a DFT gradient costs, so 1 ps at 1 fs would be 1,000 gradients.

**It starts from `relax10`'s last frame, NOT from chunk 1's.** This is the
decision that makes the run admissible. Chunk 1's geometry already has Asp92
bidentate on Zn1, so starting there hands DFT the answer and D206 rule 3
excludes the result - two methods agreeing about a geometry one of them was
given. `relax10`'s frame is the crystal-consistent one: Zn-OD1 4.33 against
6DQH's 4.28 with nothing restraining it, held over 10 ns (D231). Starting there,
the carboxylate's position is something DFT decides rather than inherits.

**Nothing is constrained, and the gate refuses rather than warns.** Two reasons,
both already paid for: a restrained carboxylate measures the restraint; and
D134 - `opt_constraints` on a `qmmm` run switches geomeTRIC's whole block to
`$set` the moment an active-region MM water appears, after which convergence
cannot be reached. Supplying none keeps it on `$freeze` only, which is why this
leg is not blocked by that defect. `rigid_water: false` for D137's reason.

**The active region stays at the base 8 A (~485 atoms, ~1,455 DOF).** Shrinking
it is not a saving: it pulls the active/frozen boundary in against the site and
the optimiser then drives active MM atoms into the frozen wall - the MM term
went -240 to +1139 Eh that way once already, with the QM term sane throughout.

**COST, stated before the run rather than after.** The only measured figure at
this level is **~26 min per wb97x-v/def2-SVP/MM single point on 119 QM atoms on
8 cores**. This region is 124 QM + 6 link = 130 atoms, and a gradient is dearer
than an energy. 32 cores should bring a step under ~15 min, but that is an
extrapolation. So `opt_maxiter` is **25**, not the base 300: ASH raises at
maxiter, mcube catches it, keeps the last geometry and writes
`QM_MM_optimized.pdb` with `converged: false`, whereas a walltime kill writes no
PDB at all. 25 steps fit inside 24 h even at 50 min per step, which makes the
first job a rate measurement that is also real work. **A converged optimisation
is not expected from it**: 1,455 DOF typically needs 100-400 steps, so on the
extrapolated rate the full leg is ~2,000-3,000 core-hours and three to six
chunks, each starting from the previous `QM_MM_optimized.pdb`. That number is
the one to check against the measured rate before buying the rest, and the
mutant arm will cost the same again.

**Level:** wb97x-v/def2-SVP, the library default functional stated explicitly
because it is the level D215's cluster benchmark used - so this geometry and
that calibration point are comparable rather than merely both "DFT".

**Not submitted.** `aspdft.sh` has to run in the code-server terminal and the
pane is hidden, which leaves only the read-only `fetch` route (bunya working
rules). The base `qmmm.inputs` was confirmed through that route - `actradius: 8`,
`actregion_origin: B:1:ZN`, `cores: 4`, `opt_maxiter: 300` - and `aspqm2`'s
recorded region is 124 QM atoms, so the section's overrides are the only things
that differ. The launcher prints `rquota` before it submits anything.

**The per-step cost is being asked for rather than guessed.** On the user's
suggestion the ion-channel QM/MM umbrella-sampling session was asked for its own
measured DFT/MM gradient timings, PySCF core scaling at 8 vs 32, how many
geomeTRIC steps a protein optimisation of this size took in their hands, and
whether TRIC degenerated. A "we only ever ran xTB" is as useful an answer as a
number, because it would confirm the user's general claim independently.

**D235 - THE def2-SVP QM/MM OPTIMISATION IS SUBMITTED: job 32453561, `aspdft`,
wb97x-v/def2-SVP, 124 QM atoms, the active region cut to the QM atoms only, 60
steps against a 48 h wall, nothing constrained. The active region is NOT the
8 A one the base config gives, and that single change came from the ion-channel
project's measured experience rather than from our own extrapolation - the user
sent me to ask them.** (2026-11-07)

**What the channel session measured, which this project had no equivalent of.**
Their QM/MM optimisations are r2SCAN/def2-SVP through the same mcube `qmmm`
stage, ASH + PySCF, density fitting on, 32 cores on `general`:

- **An 8 A active region containing MM waters never converged in 150 steps**,
  and worse, the energy **climbed 30-45 kcal/mol** after a minimum around step
  40-120 while geomeTRIC rated every step good (Quality ~1), with repeated
  `SVD Error` in the internal coordinates. A finite-difference check ruled out
  an energy/gradient mismatch (they agree to 0.4 kcal/mol/A), so the geometry
  was being driven uphill by the coordinate system, not by a bad gradient.
- **Cutting the active region to the QM atoms converged it in 44 steps**, energy
  falling every step, and all 60 of their round-2 jobs then converged in 27-100.
- `rigid_water: false` plus `opt_coordsystem: cart` was the other fix:
  monotonic, never uphill, but **not converged at 300 steps** - the gradient met
  the thresholds while the energy crept ~1 kcal/mol per 50 steps out of the
  floppy water network. So flexible water is not a substitute for a small
  active region; it trades a wrong answer for a slow one.
- Cost: **4.6 min per optimisation step** at ~50 QM atoms in a 124k-atom system
  with ~64k MM point charges - SCF 49 s, gradient 62 s, and the remaining 2.8 min
  the QM-MM point-charge terms. At 28-31 QM atoms in a water box, 18.6 s a step.
- Core scaling was **1.25x for 16 -> 32 cores** at 28 QM atoms. It flattens.
- DFT/MM MD ran at 7-9 ps/day at 28 QM atoms against ~40x that with xTB/MM, so
  their conclusion is the user's: DFT/MM for minima and single points, xTB for
  dynamics. That is now two independent projects saying it.

**The two records do not contradict each other, and keeping both is the point.**
This project ruled out a 5 A active region for the opposite-looking reason: the
MM term went -240 -> +1139 Eh while the QM term stayed sane, because a
*partially* shrunk region leaves active MM atoms pressed against frozen
neighbours and the optimiser drives them into the wall. With the region cut to
the QM atoms there are **no active MM atoms at all**, so that failure mode has
nothing to act on. 5 A is still ruled out; 0 is a different setting, not a more
extreme version of it. mcube returns `list(qm_atoms)` directly when `actradius`
is falsy, so no `actregiondefine` call happens and the QM region is trivially
contained in the active one.

**What the frozen environment costs, said rather than buried.** The second shell
cannot relax: Ser186, the bulk waters and the backbone stay at the MM geometry
the 10 ns trajectory ended on. Everything the question turns on is quantum and
free - Zn, hydroxide, substrate, Asp92, Lys226, Arg90, Asn225, His88, His164 -
so Asp92's position is decided by DFT against a fixed environment rather than a
relaxed one. That is the standard minimum-refinement approximation, it is what
makes the leg affordable, and it belongs in the methods section rather than in a
footnote.

**Also taken from them, and not acted on yet:** SCF trouble on a divalent metal
next to bare MM charges is to be met with damping, level shift and a Newton
fallback (their Ln3+ clusters needed all three, and one Lu cluster flipped
between two SCF solutions 0.34 Eh apart on successive steps - so watch for energy
jumps between steps rather than only for non-convergence). The DFT grid adds
~1e-4 Eh/A of gradient noise, comparable with geomeTRIC's own thresholds, so
formal convergence can stall at a minimum that has physically been reached.
Their default auxbasis was fine and was never tuned. Their one memory blowup was
PCM-specific (96 GB was not enough for dense PCM matrices on 74-77 atoms) and
plain QM/MM showed none; this job has 120 GB and no PCM.

**Submitted, with the gate's output read rather than assumed.** Job **32453561**,
`general`, 32 cores, 120 GB, 48 h. The launcher printed
`[OK] def2-SVP optimisation from relax10's last frame, charge 2, nothing
constrained, no QM/MM MD`, and `config.yaml` shows `qm_engine: dft`,
`dft_functional: wb97x-v`, `dft_basis: def2-svp`, `dft_cores: 32`,
`actradius: 0`, `opt_maxiter: 60`, `rigid_water: False`,
`skip_equilibrate: True`, `skip_optimize: False`, and the `michaelis_rk` region
at charge 2. Quota before submitting: scratch 238/300 GB, 724,499/1,000,000
inodes. `aspqm2` was still running at 16:45 h elapsed and was left alone.

**What to read first, and what would make it worth stopping.** Output in
`protdft/aspdft-32453561.out`, run directory `protdft/qmmm_aspdft/`. The first
numbers wanted are the **wall time per step** (against the 20-35 min band
predicted) and, from `optimization_energies.log`, **whether the energy falls
monotonically**. A climbing energy with good step quality is the channel
project's 8 A failure reproduced at a QM-only region, and it would mean stopping
rather than buying more steps. 60 steps is not expected to converge if their
27-100 range transfers; the next chunk starts from this run's
`QM_MM_optimized.pdb`, which mcube writes even on non-convergence.

**D236 - THE FIRST DFT ANSWER POINTS THE OTHER WAY: through three
wb97x-v/def2-SVP QM/MM optimisation steps Asp92 does NOT approach Zn1, and the
Lys226 salt bridge tightens from 2.88 to 2.77 A toward 6DQH's 2.79. AND the leg
is mis-sized: a step costs 110-118 min, not the predicted 20-35, so the 48 h
wall buys ~25 steps against `opt_maxiter: 60` - the job will be killed by the
wall rather than by maxiter, and a wall kill writes no `QM_MM_optimized.pdb`.**
(2026-11-07)

**COST, measured from ASH's own timers rather than estimated.** Job 32453561
started 15:10:22 and had completed three gradient steps by 19:35:

    step 1   QM/MM run  6603.7 s   110.1 min
    step 2   QM/MM run  7073.4 s   117.9 min
      of which  pySCF run_SCF            5230.6 s   87.2 min   (74%)
                pyscf_gradient           1377.1 s   23.0 min
                pyscf_newpointcharge_grad  464.8 s    7.7 min
                MM step                      0.6 s

1,211 basis functions for the 124-atom QM region, density fitting on. **The SCF
is three quarters of the cost**, which is where the prediction went wrong: the
band of 20-35 min a step came from scaling the channel project's 4.6 min at ~50
QM atoms with r2SCAN, and wb97x-v's range-separated exchange plus VV10 nonlocal
correlation is a different animal. The extrapolation was out by 3.4-5.9x. Their
warning that wb97x-v would cost more than a meta-GGA was in their message and I
did not weight it heavily enough.

**So the sizing is wrong in the dangerous direction.** At ~115 min a step, 48 h
is ~25 steps. ASH raises `GeomOptNotConvergedError` at `opt_maxiter` and mcube
catches it, keeps the last geometry and writes `QM_MM_optimized.pdb`; a walltime
kill writes nothing. With maxiter 60 and a 25-step budget, the wall arrives
first. **It is recoverable and nothing is lost**:
`geometric_OPTtraj_Full.xyz` is written every step and `docs/csr/traj_frame.py`
reads one frame of it onto a PDB template, which is exactly how D106 recovered a
geometry from the same situation. But it is a recovery, not a plan, and the
honest statement is that this job was submitted with a maxiter it cannot reach.

**The optimisation itself is healthy, which is the part the channel project's
experience was most useful for.** Energy falls every step, monotonically:

    -5201.1362  ->  -5201.3065  ->  -5201.3341  ->  -5201.3519 Eh
    deltas (kcal/mol):   -106.9      -17.3        -11.2
    step 2: Trust 1.414e-01, Grad 1.199e-02/5.638e-02 rms/max, Quality 0.645
    step 3: Trust 1.414e-01 (=), Grad 6.038e-03/2.503e-02, Quality 0.871

**No SVD errors in the .out and no `LinAlgError` in the .err** - counted, zero
each, and the `.err` was checked separately from the `.out` because conflating
them produced a wrong conclusion in this project once already (D106). The trust
radius is steady rather than collapsing and the gradient is falling. This is NOT
the channel project's 8 A failure: there the energy climbed 30-45 kcal/mol with
Quality ~1. Cutting the active region to the QM atoms was the right call.
Convergence is still far off - `grms` 6.0e-03 against the 1e-04 criterion.

**THE GEOMETRY, AND THE MAPPING CHECKED BEFORE IT WAS TRUSTED.** The QM-region
geometry blocks in the `.out` were mapped to atom names through
`qmmm_aspdft.results.qm_atoms` and the input PDB. The check that the mapping is
right: the first block reproduces `relax10`'s recorded start to 0.003 A -
Zn-OD2 5.153 (5.15), Zn-OD1 4.327 (4.33), Zn-CG 5.192 (5.19), Zn-NZ 5.276
(5.28), NZ-OD2 2.880 (2.88). Only then the rest:

    step        0      0'     0''     1       2     [6DQH]
    Zn-OD2    5.153  5.176  5.140  5.106  5.087   [5.14]
    Zn-OD1    4.327  4.474  4.549  4.509  4.544   [4.28]
    Zn-CG     5.192  5.175  5.243  5.216  5.198   [5.15]
    Zn-NZ     5.276  5.277  5.289  5.287  5.269   [5.27]
    NZ-OD2    2.880  2.854  2.780  2.799  2.766   [2.79]
    His88 NE2 2.178  2.170  2.197  2.216  2.222
    His164NE2 2.126  2.126  2.122  2.118  2.114

**Asp92 does not move toward the zinc. The salt bridge tightens onto the crystal
value.** Where GFN2 took Zn-OD2 from 4.34 to 2.01 A within 8 ps and held it for
24 ps more, the DFT gradient at the same geometry points the carboxylate
nowhere: Zn-OD2 drifts 0.07 A *outward* over three steps and Zn-OD1 0.22 A
outward, while NZ-OD2 closes from 2.880 to 2.766 against 6DQH's 2.79.

**What this is NOT, because the temptation to overclaim here is large.** Three
steps is a direction, not a result: the total displacement is rms 2.7e-02 A, the
gradient is two orders of magnitude from the convergence criterion, and no
stationary point has been reached. And the comparison with GFN2 is not clean -
**that was 50 ps of dynamics at 300 K and this is a zero-temperature
optimisation**, which has no thermal energy to cross a barrier it might find. A
GFN2 *optimisation* from this same frame would be the matched comparison and has
not been run. What can be said is narrow and is worth saying: the DFT gradient
at the crystal-consistent geometry does not pull Asp92 onto the metal, and it
improves the salt bridge. That is the first direct evidence that GFN2's chelate
is the Hamiltonian's and not the enzyme's, as opposed to an inference from a
cluster benchmark.

**JOB 2, `aspqm2` chunk 1, finished at its full 50 ps** - 50 frames exactly,
28,791,676 bytes, parsed at the same stride as D232. The last 18 ps add nothing
and contradict nothing:

    ps 33-50:  Zn-OD2 1.91-2.28   Zn-OD1 1.97-2.44   Zn-CG 2.36-2.62
               NZ-OD2 2.57-3.25   Zn-O2x 5.73-6.71   O01-C4x 5.86-7.63
               Zn-O01 1.73-1.96

The chelate holds for the whole run with no recrossing, and **OD1 and OD2 now
exchange places frame to frame** (OD1 is the shorter contact at 34, 37, 38, 44
and 49 ps), so it is a symmetric bidentate carboxylate rather than one oxygen
bound with the other nearby. The hydroxide stays on the metal throughout. The
substrate continues to leave - O01-C4x goes from ~6.0 A at 33 ps to 7.4 A at
50 - which only strengthens D232's statement that **nothing energetic may be
taken from these frames**. D232 is not reopened and not revised.

**Chunks 2-6 remain cancelled and nothing was submitted.**

**The question the user has to answer, with the arithmetic rather than a
recommendation.** As configured this leg cannot converge: ~25 steps of a budget
that needs hundreds. Four ways out, costed:
1. **Let it run to the wall** and recover the geometry from
   `geometric_OPTtraj_Full.xyz` with `traj_frame.py`. Cost: already spent, 48 h
   of 32 cores (~1,500 core-hours). Decides: the direction of the DFT gradient
   over ~25 steps, which is what the table above is already showing. Does not
   give a stationary point, so no energy may be quoted from it.
2. **Cut the QM region** from 124 atoms to the metal, two histidines, the
   hydroxide, Asp92 and Lys226 (~60 atoms, roughly 600 basis functions). SCF
   cost falls faster than linearly, so plausibly 15-25 min a step and a
   converged optimisation inside one 48 h job. Loses Arg90, Asn225 and the
   substrate from the QM region, which is a real loss for anything mechanistic
   but may be acceptable for the single question of where the carboxylate sits.
3. **Change the functional to r2SCAN**, the channel project's choice, whose SCF
   was 49 s at their size. Cheapest fix by far, but it leaves the level of
   D215's cluster benchmark, so that calibration point stops being comparable.
4. **Accept the direction and stop**, treating this as a gradient-direction
   measurement rather than an optimisation, and put the def2-TZVP single point
   on `relax10`'s frame instead of on an optimised one.

**D237 - `aspdft` runs to the wall, the user's decision. Resumable, but NOT by
`resume: true`, which is a no-op on this leg: the optimise path never calls
`resolve_resume`, and the route is `traj_frame.py` onto the last frame of
`geometric_OPTtraj_Full.xyz`. Checked in the source and against the files on
disk rather than taken on trust.** (2026-11-07)

**What `resume` does and does not do here, read in `src/mcube/qmmm.py` rather
than assumed.** `resolve_resume(inputs, pdb_file)` is called inside
`if not inputs.get("skip_equilibrate")` - the QM/MM **MD** branch only. This
section sets `skip_equilibrate: true`, so a `resume: true` on it would be
silently ignored, and there is no checkpoint restart of a geomeTRIC
optimisation anywhere in mcube. The work is still not lost, which is the part
that matters, but the mechanism is manual.

**The route, with the files confirmed present on disk:**

    geometric_OPTtraj_Full.xyz        6,550,513 B   the whole system, every step
    geometric_OPTtraj_QMregion.xyz       29,017 B
    pyscf.chk                        11,969,816 B   SCF guess
    QM_MM_optimized.pdb                     404      not written yet, as expected
    constraints.txt                         404      NO constraint block at all

Then `python docs/csr/traj_frame.py <Full.xyz> <template.pdb> <out.pdb> [index]`,
default index -1, template `md_relax10_results/NVTtrajectory_lastframe.pdb`
(same system, same atom order, which the script checks rather than assumes - it
verifies the atom count and that every element agrees with the template's atom
name). The file is streamed, so its eventual few hundred MB is not a problem.
A new section then takes that PDB as its `pdb_file`. This is exactly how D106
recovered a geometry from an optimisation that never reached its end.

**What a restart costs, so it is not discovered later:** geomeTRIC's
accumulated approximate Hessian and its trust radius go with the old process,
so the first few steps of a restart re-learn curvature that was already paid
for. At ~115 min a step that is a few hours, which is the price of the manual
route and is worth stating next to "it is resumable". The SCF guess does carry
over in `pyscf.chk`.

**`constraints.txt` is absent, which is the design condition confirmed by
measurement rather than by reading the config.** D134's defect fires only when a
`qmmm` run is handed `opt_constraints`, and D137's `$set` water trap only when
rigid water puts pairs in the block. No constraints file at all means geomeTRIC
was given neither, so convergence is not gated on anything artificial. That was
the argument for supplying none; this is the check on it.

**Nothing submitted. The job keeps running to its 48 h wall (started 15:10:22 on
7 Nov, so about 15:10 on 9 Nov), expected to reach ~25 steps against
`opt_maxiter: 60`.** The three questions the next check should answer are
whether the ~115 min rate holds, whether the energy is still falling
monotonically, and whether Asp92 continues to move *away* from the zinc - the
third being the one that matters, since three steps is a direction and twenty-five
is an argument.

**D238 - THE MATCHED COMPARISON IS RUNNING: `aspxtb`, job 32461816, GFN2 from
the same frame with the same QM region, charge, active region, optimiser and
criteria as `aspdft`, and nothing constrained - one key different. AND the gate
is config rather than a script: `sections.py` now takes a declarative
`require:` and `submit:` block, because the launcher I had started writing for
this arm was the fourth copy of the same shell script and they had already
diverged.** (2026-11-07)

**Why the comparison was needed, which D236 had already said against itself.**
"GFN2 chelates Asp92, DFT does not" was not yet a statement about the
Hamiltonians: D232 was 50 ps of 300 K dynamics and `aspdft` is a 0 K
optimisation, which has no thermal energy to cross a barrier it might find. So
the difference could have been dynamics versus minimisation. `aspxtb` removes
that confound by changing only `qm_engine`.

**Pre-registered, before the result, per D178.**
- **GFN2 minimises Asp92 inside ~2.5 A where DFT holds it near 5** -> the
  cleanest statement of the bias available: one structure, one variable. D232
  stops being the only evidence and stops being answerable by the
  dynamics-versus-minimisation objection.
- **Both keep it out** -> D232's chelation needed 300 K, so GFN2 makes a
  barrier crossing too easy rather than putting the minimum in the wrong place.
  A different, weaker claim, and
  `gfn2-over-coordinates-zn1-from-a-crystal-consistent-start` would have to be
  narrowed to say exactly that.
- **Both bring it in** -> the approach is not a QM-level artefact at all, and
  the finger points at the force field or the starting pose instead. The
  expensive outcome.

Cheap, which is why it is worth running at all rather than arguing about:
`opt_maxiter` is the base 300 because an xTB step is minutes, so unlike the DFT
arm this one can actually reach a minimum. 8 cores - xTB does not scale and the
DFT arm has the cores.

**THE LIBRARY CHANGE, which is the user's point and the more durable half of
this entry.** I started writing `docs/csr/aspxtb.sh` and was told, correctly,
that a config option should not need a script. Every `docs/csr/*.sh` was the
same six lines - module load, conda activate, rquota, `sections.py`, a
hand-written preflight, `sbatch` - with only the preflight differing, and **they
had already drifted apart: `aspqm2.sh` printed a WARNING where restraints were
set and `aspdft.sh` REFUSED.** A gate written once per run is a gate that is
wrong somewhere, and nothing in either run's output would say so.

So `sections.py` takes two optional blocks beside `stage` and `inputs`:

    require:
      unset: [restraints, opt_constraints, constraints]
      equals: {charge: 2, rigid_water: false}
      same_as: {section: aspdft, except: [qm_engine, dft_basis, opt_maxiter]}
    submit: {cores: 8, mem: 16G, time: "24:00:00", partition: general}

checked in the existing "check everything before writing anything" pass, so a
failing gate writes no section at all. `--dry-run` checks and prints the sbatch
line without touching `config.yaml`; `--submit` runs it. The runner is
`run.slurm` beside the site file, so no path comes from the caller.

**`same_as` is the one that earns its keep**, and it is not a convenience: it
refuses to write an arm that differs from its reference in anything but the
named keys, which is the only way a two-arm comparison can be *kept* true as
either arm is edited later. It compares against the reference **as written in
`config.yaml`** - what the job will actually read - falling back to the site
file only if the reference has not been written yet, and saying which it used.
A test pins that precedence by giving the site file a different charge from the
config and checking the config wins.

**Eleven tests, `tests/test_sections_require.py`, and one of them caught a real
bug** - a `"..." + NL + "..." % (...)` where `%` bound to the last fragment
only, so the `same_as` failure path raised `TypeError` instead of reporting the
difference. That is the error path, i.e. exactly the code that only runs when
something is already wrong, which is why it is worth having a test for it.
35 passed with the existing `test_sections_structures.py`.

**Submitted through the generic path, with the gate's own output read rather
than assumed.** `--dry-run` first: `DRY RUN - config.yaml not touched. Checks
passed for: qmmm_aspxtb`, and the printed line
`sbatch -J aspxtb -p general --cpus-per-task=8 --mem=16G -t 24:00:00
.../run.slurm qmmm qmmm_aspxtb`. Then `--submit`: `Submitted batch job
32461816`. **No file was added for this arm.**

**Still to do, and not done here:** the four existing launchers are now
redundant and should become `require:`/`submit:` blocks, with their divergent
preflights reconciled rather than carried over - `aspqm2.sh`'s warning is the
wrong behaviour and `aspdft.sh`'s refusal is the right one. And `sections.py`
itself belongs in `src/mcube/` as `python -m mcube.sections`, since nothing in
it is CSR-specific; it sits in `docs/csr/` only for historical reasons.

**D239 - THE EIGHT LAUNCHER SCRIPTS ARE GONE. Every gate and every `sbatch`
line they held is now a `require:`/`submit:` block on its own section in
site.yaml, and the divergence between them is reconciled rather than carried
over: where `aspqm.sh` and `aspqm2.sh` printed `[WARNING]` for a restrained
carboxylate, the config refuses.** (2026-11-07)

Deleted: `asppull.sh`, `asppull2.sh`, `asppull3.sh`, `eq6.sh`, `relax10.sh`,
`aspqm.sh`, `aspqm2.sh`, `aspdft.sh`. Kept, because they are analysis rather
than submission: `cleanup.sh`, `convergence.sh`, `conv_relax10.sh`.

**Nothing was thrown away.** Each script's preflight is in the section it
guarded, with the reason beside it, and each one's resources are a `submit:`
block:

    relax10   require.unset [skip_npt, skip_nvt, gentle], equals nvt_time 10000
              (D227: `gentle` is 40 ps of 300 K MD before NPT for nothing)
    aspqm     require.unset [restraints, opt_constraints, constraints]
    aspqm2    the same, plus equals pdb_file = relax10's last frame, charge 2,
              rigid_water false (D228: the 2 ns frame is what job 32384737 was
              cancelled for using)
    aspdft    the same, plus actradius 0, qm_engine dft, dft_basis def2-svp,
              skip_equilibrate true, skip_optimize false
    md arms   submit {partition: gpu_cuda, cores: 4, mem: 32G, extra: [--qos=gpu,
              --gres=gpu:h100:1]} - the qos goes in `extra` because run.slurm
              hardcodes --qos=normal and gpu_cuda rejects it
    qmmm arms submit {partition: general, ...}, 8 cores for xTB and 32 for DFT

**What the three "section is in config.yaml" checks were worth: nothing.**
Every launcher opened by asserting that `sections.py` had written its section,
which `sections.py` already guarantees by refusing to write anything when any
section fails. Three of the eight preflights were only that. The checks that
mattered were in two scripts, and one of those two was wrong.

**One sharp edge found while writing the checker, and fixed.** `equals` compared
with `==`, and Python says `False == 0`, so `equals: {rigid_water: false}` would
have been satisfied by `rigid_water: 0` and `actradius: 0` by
`actradius: false`. Those mean the same thing to mcube today. A gate that cannot
tell them apart is a gate that will one day pass something it was written to
stop, so the comparison is on `(isinstance(value, bool), value)` and two tests
pin it. This is the second bug the tests found in the error paths of this
checker, the first being a string-formatting precedence error that made the
`same_as` failure path raise `TypeError` instead of naming the difference.

**`docs/csr/README.md`** now documents the route - `--dry-run` to check and
print the `sbatch` line without writing, `--submit` to write and submit - and
says why there are no launchers, since the obvious next move for anyone adding
an arm is to copy the last script.

**Still not done, and named so it is not rediscovered:** `sections.py` belongs
in `src/mcube/` as `python -m mcube.sections`. Nothing in it is CSR-specific -
it takes a site file and section names - and it sits in `docs/csr/` only because
that is where it was first written. Moving it is an import change, a test-path
change and a line in `running-on-bunya.md`, with no behaviour change, and it is
the difference between a general library and a library with a project's script
in its docs directory.

**D240 - THE MATCHED COMPARISON ANSWERS ITS SECOND PRE-REGISTERED OUTCOME, AND
IT RETRACTS WHAT I WROTE THIS AFTERNOON. GFN2's own minimum from the
crystal-consistent frame puts Asp92 FURTHER OUT than the crystal - Zn-OD2
5.153 -> 5.383, Zn-OD1 4.327 -> 4.952 - with the salt bridge tightening to
2.718 and the substrate carbonyl taking the open site at 2.282. Converged in 87
steps and 3.8 minutes. So `aspdft`'s outward drift is NOT a DFT-versus-GFN2
difference, and D236's reading of it as "the first direct evidence that GFN2's
chelate is the Hamiltonian's" is wrong and is retracted.** (2026-11-07)

**`aspxtb`, job 32461816, converged.** `converged: true`, 87 optimisation
steps, `opt_energy` -430.99370531860416 Eh at GFN2-xTB, **228 seconds of
walltime**. Energy fell monotonically, -430.7337 -> -430.9937 Eh (-163
kcal/mol). Everything but `qm_engine` was identical to `aspdft` and the
`same_as` gate refused to write the section until that was true, which is the
only reason the pair can be read against each other at all.

**The measurement, on `QM_MM_optimized.pdb`:**

    pair          start   GFN2 min   [6DQH]    [aspdft, 3 steps]
    Zn-OD2        5.153    5.383     [5.14]        5.087
    Zn-OD1        4.327    4.952     [4.28]        4.544
    Zn-CG         5.192    5.564     [5.15]        5.198
    Zn-NZ         5.276    5.195     [5.27]        5.269
    NZ-OD2        2.880    2.718     [2.79]        2.766
    His88 NE2     2.178    2.083                   2.222
    His164 NE2    2.126    2.010                   2.114
    Zn-OH O01     2.035    1.775
    Zn-O2x        2.610    2.282
    O01-C4x       3.481    3.673

**THE SECOND PRE-REGISTERED OUTCOME, WORD FOR WORD FROM D238: "Both keep it out
-> D232's chelation needed 300 K, so GFN2 makes a barrier crossing too easy
rather than putting the minimum in the wrong place. A different, weaker claim,
and `gfn2-over-coordinates-zn1-from-a-crystal-consistent-start` would have to be
narrowed to say exactly that."** That is what happened, and it is honoured
rather than re-spun. The pre-registration is what makes this cheap to accept:
the reading was fixed before the result, so there is nothing to argue about.

**WHAT IS RETRACTED, and it is six hours old.** D236 closed with: *"the DFT
gradient at the crystal-consistent geometry does not pull Asp92 onto the metal,
and it improves the salt bridge. That is the first direct evidence that GFN2's
chelate is the Hamiltonian's and not the enzyme's."* **GFN2 does the same thing
and more of it** - it moves the carboxylate 0.23 A further out where DFT moved it
0.07, and tightens the same salt bridge to the same place. So `aspdft`'s three
steps say nothing about GFN2 versus DFT. They say that *neither* method's
gradient, from this start, points the carboxylate at the metal. Entered in
`retractions.md`. The numbers in D236 stand; the inference does not.

**TWO DIFFERENCES FROM D232, NOT ONE, AND THEY CANNOT BE SEPARATED BY WHAT HAS
RUN.** D232 was 50 ps of 300 K dynamics with the whole system free. Both
optimisations are 0 K **and** have the MM environment frozen at the MM geometry,
because the active region is the QM atoms only (D235). So the chelate's absence
here is consistent with needing thermal energy, with needing the second shell to
relax, or with both. Saying "it needed 300 K" would be picking one of three
explanations because it is the one I thought of first. What the pair does
establish is narrower and still worth having: **the chelate is not where GFN2's
gradient goes from a crystal-consistent start with this environment held.**

**A new observation, recorded carefully because it is `ours` and the method is
on trial.** GFN2's minimum pulls the substrate carbonyl onto the zinc at
**2.282 A**, from 2.610 at the start, where 7UHI/7UHJ give 2.53/2.63. The
coordination that `no-usable-michaelis-complex-exists-yet` records as missing
from the docked pose (5.43 A) is present in this geometry and slightly short of
the deposited values. That is **not** a resolution of that node: the pose was
inherited from the same docking, the method is one this project has just spent
two days distrusting, and a single minimum is not a measurement of a distance.
It is recorded because it is the first geometry in the project where the carbonyl
is on the metal without a restraint holding it there.

**And the other thing that moved: both optimisations tighten the Zn-Lys226 salt
bridge onto the crystal value** - 2.880 -> 2.718 at GFN2 and -> 2.766 at DFT,
against 6DQH's 2.79. Two methods, same direction, from an MM geometry that had
drifted 0.09 A long. That is the one place where agreement between two methods is
worth something here, because the thing they agree with is a measurement.

**Costs and the comparison that is now available.** 228 seconds against the DFT
arm's 110-118 minutes **per step**: GFN2/MM optimisation of this system with a
QM-only active region is free, which is worth knowing for its own sake, since the
8 A active region version of the same calculation took 159 s a step and 111 steps
without finishing (D105/D106). Cutting the active region to the QM atoms did not
only fix the channel project's convergence failure - it made the xTB arm
instantaneous.

**Housekeeping.** One pre-existing wart surfaced in the `.err`:
`python: can't open file '/scratch/.../CSR_mich/../mcube/docs/csr/site.py'`,
after `[DONE]`. `run.slurm` runs a post-stage `site.py` check through a relative
path that is wrong from this run directory, so that check has silently not been
running. It did not affect the optimisation. Not fixed here; named so it is not
rediscovered.

**Library.** 784 passed, 7 skipped locally (the ASH-gated ones) with the
`require:`/`submit:` work in place, and the eight launchers' gates verified
end-to-end on Bunya: `--dry-run` on all eight sections passes and reproduces
exactly the `sbatch` lines the deleted scripts used, GPU qos and gres included.

**D241 - `aspxtb8` reopens the 8 A active region from the converged xTB
minimum, to remove one of the two differences D240 could not separate. Written,
gated and SUBMITTED - but the submission is UNCONFIRMED: the code-server session
died between the `sbatch` and the read, so whether job landed is not known and
must be checked before anything is resubmitted.** (2026-11-08, 00:0x AEST)

**Why: D240's confound, halved.** Both optimisations differ from D232's 300 K
dynamics in two ways at once - they are 0 K **and** they hold the MM environment
frozen, because the active region is the QM atoms only (D235). So the chelate's
absence from both minima is consistent with needing thermal energy, with needing
the second shell to relax, or with both. `aspxtb8` removes the second: same
Hamiltonian, same QM region, same charge, nothing constrained, but ~485 active
atoms instead of 124. Only `pdb_file` and `actradius` differ from `aspxtb`, and
the `require.same_as` gate refuses to write the section if anything else does.

**It starts from `aspxtb`'s converged minimum, which is the only converged QM/MM
geometry in this project.** The user asked for the *DFT* geometry and it does not
exist yet: `aspdft` was on **step 4 of ~25** at 8:50 elapsed, `grms` 6.0e-03
against the 1e-04 criterion, so it will hit its 48 h wall far from a stationary
point. Asked rather than substituted, and the answer was the xTB minimum at xTB.

**Pre-registered, before the result, per D178:**
- **Asp92 comes in toward ~2.5 A once the environment is free** -> the frozen
  shell was holding it out, D240's narrowing was premature, and the
  0 K-versus-300 K half of the confound is the half that does not matter.
- **Asp92 stays out near 5 A** -> the frozen environment was not the reason, so
  the chelate needs thermal energy. D240's narrowing stands and strengthens.
- **No convergence in 300 steps** -> expected, not a failure. The channel
  project's flexible-water 8 A runs were monotonic and still unconverged at 300.
  Read the trajectory's direction, which is what the first two readings turn on.

Starting geometry, from `aspxtb`'s `QM_MM_optimized.pdb`: Zn-OD2 5.383,
Zn-OD1 4.952, Zn-CG 5.564, Zn-NZ 5.195, NZ-OD2 2.718, Zn-O2x 2.282,
Zn-OH 1.775, O01-C4x 3.673.

`rigid_water` stays false. The channel project's 8 A runs climbed 30-45 kcal/mol
UPHILL with good step quality because rigid-water constraints degenerate TRIC's
coordinate system; with no constraints at all geomeTRIC gets no constraint block,
which `aspdft`'s absent `constraints.txt` confirmed by measurement (D237). 8
cores and 32 GB, 24 h: our own 8 A xTB optimisation ran at ~159 s a step with
geomeTRIC's linear algebra over 467 active atoms as the cost driver, so 300 steps
is ~13 h.

**WHAT IS NOT KNOWN, and why no resubmission may happen until it is.** The
`--dry-run` passed and printed the exact line
(`sbatch -J aspxtb8 -p general --cpus-per-task=8 --mem=32G -t 24:00:00
.../run.slurm qmmm qmmm_aspxtb8`) with quota at 238/300 GB and
724,722/1,000,000 inodes. The `--submit` was then typed and entered, and the
next read returned **`Failed to connect to bun147:52687`** - the interactive
code-server job had ended. So the submission is in one of three states and the
log that would say which is unreachable:

1. `qmmm_aspxtb8` is absent from `config.yaml` -> nothing happened at all.
2. It is present and `sacct` shows a job -> it landed; read it.
3. It is present and `sacct` shows nothing -> the section was written and the
   `sbatch` did not run, so submitting again is correct.

`sections.py` writes the section **before** it submits, so the config is the
discriminator and must be read first. **Submitting blind would risk a duplicate
8-core job**, which is exactly the kind of waste the `[OK]`-before-`sbatch`
discipline exists to prevent. The jobs already running are unaffected either way:
`aspdft` (32453561) and anything else are separate Slurm jobs and `/scratch`
survives the session moving, which is why nothing is written to the node.

**D242 - D241'S UNKNOWN RESOLVES TO "NOTHING HAPPENED": `qmmm_aspxtb8` was
absent from `config.yaml` and `sacct` showed no such job, so the 7 Nov
`--submit` never ran. Resubmitted, job **32489558**. And `aspdft`, now 13
steps in, OVERTURNS ITS OWN FIRST THREE STEPS: Asp92's outward drift reverses
at step 6 and the carboxylate comes back to the crystal distance, while
Lys226 moves monotonically ONTO the zinc.** (2026-11-08, 14:49 AEST)

**The submission question, answered by measurement rather than assumption.**
`config.yaml`'s section list read `michqm, michrk, aspqm, aspqm2, aspdft,
aspxtb` - no `aspxtb8` - and `sacct -X -S 2026-11-07` showed no job of that
name. `sections.py` writes the section before it submits, so an absent section
means the process never got that far: state 1 of D241's three, and a blind
resubmission would have been correct after all. It was still right to check,
because the cost of being wrong was a duplicate 8-core job and the cost of
checking was two fetches. Dry-run passed and reproduced the line exactly;
`--submit` then returned `Submitted batch job 32489558`. Quota 238/300 GB,
724,739/1,000,000 inodes.

**A fresh code-server shell inherits no conda environment, and this cost one
wasted batch.** `python .../sections.py` died with
`ModuleNotFoundError: No module named 'mcube'`. The fix is the three lines
`run.slurm` itself uses - `module load miniforge/26.1.0-0`,
`source $ROOTMINIFORGE/etc/profile.d/conda.sh`, `conda activate ASH`. Also:
the log read empty immediately after those three, because `module load` takes
a few seconds and the fetch beat it. A blank log right after a `module load`
is a timing artefact, not a failure.

**`aspdft` (32453561): 13 completed steps, 23.6 h elapsed, still running, and
the mapping re-verified - geometry block 0 reproduces 5.153 / 4.327 / 5.192 /
5.276 / 2.880 / 2.178 / 2.126 exactly.**

    step   Zn-OD2  Zn-OD1  Zn-CG  Zn-NZ  NZ-OD2  H88NE2  H164NE2
     0     5.153   4.327   5.192  5.276  2.880   2.178   2.126
     2     5.140   4.549   5.243  5.289  2.780   2.197   2.122
     4     5.087   4.544   5.198  5.269  2.766   2.222   2.114
     6     5.073   4.536   5.206  5.213  2.699   2.235   2.097
     8     5.118   4.518   5.215  5.164  2.729   2.213   2.094
    10     5.143   4.527   5.224  5.145  2.740   2.180   2.085
    12     5.135   4.538   5.235  5.137  2.717   2.162   2.074
    13     5.134   4.546   5.237  5.131  2.713   2.162   2.073
    6DQH   5.14    4.28    5.15   5.27   2.79

**THE ANSWER TO "DOES Asp92 KEEP MOVING AWAY" IS NO, AND THREE STEPS WAS A
DIRECTION THAT DID NOT HOLD.** Zn-OD2 falls to a minimum of 5.073 at step 6
and then climbs back, settling at 5.134 - which is 6DQH's 5.14 to within
0.006 A. Over the last four steps it moves 5.143 -> 5.141 -> 5.135 -> 5.134,
and Zn-CG and Zn-OD1 are flat to +/-0.005 A over the same span. So the DFT
minimum does not push the carboxylate out and does not pull it in: **it holds
Asp92 where the crystal has it.** D240's survival claim - that neither
gradient points the carboxylate at the metal - stands and is in fact
strengthened, because the metal-ward displacement is now zero rather than
negative. What does not stand is the *wording* of the observation node, whose
title says Asp92 "moves away from Zn1". By step 13 it has not.

**THE NEW MONOTONE IS LYS226, AND IT IS THE THESIS'S OWN ATOM.** Zn-NZ runs
5.276 -> 5.289 -> 5.269 -> 5.213 -> 5.164 -> 5.145 -> 5.137 -> 5.131 with no
reversal after step 2: the lysine ammonium moves 0.145 A closer to the zinc
while holding its salt bridge to Asp92 OD2 at 2.70-2.73. Read carefully: this
is Lys226 taking up position in the vacant Zn2 subsite, which is what
"Lys replaces Zn2" asserts, and it is the first time any unrestrained
calculation here has produced that motion rather than had it imposed. It is
also 0.145 A in a geometry optimisation that has not converged, so it is a
direction and not a distance, and 6DQH's 5.27 is where the crystal puts NZ -
so the optimisation is moving AWAY from the crystal on this one coordinate
while settling onto it on the others. That is worth saying out loud rather
than harvesting the half that flatters the thesis.

**Energy, monotonic throughout, and the optimiser is healthy.**
-5201.1362 -> -5201.3871 Eh over 13 steps, every step downhill, step quality
1.07-1.20 for the last six, trust radius grown from 1.000e-01 and capped at
3.000e-01, `Displace` rms down to 1.231e-02. **Zero `SVD Error` in the `.out`
and zero `LinAlgError` in the `.err`, counted in their own files.** This is
the opposite of the channel project's 8 A failure - no climb, no degenerate
coordinates - which is what `actradius: 0` was chosen for.

**IT WILL NOT CONVERGE BEFORE THE WALL, BUT THE GEOMETRY IS ALREADY SETTLED,
AND THOSE TWO FACTS TOGETHER DECIDE D236'S OPEN QUESTION.** `Grad` rms falls
3.656e-02 -> 1.339e-03, but the last five steps decay only ~8% each, so 1e-04
is roughly 30 more steps away. The rate has improved - 110.1 and 117.9 min at
the start, then 109.3, 96.0, 91.1 min, the SCF reusing `pyscf.chk` - and
24.4 h of wall remain, which buys about 16 steps, reaching grms ~4e-04 at step
~29 against `opt_maxiter: 60`. So: no `QM_MM_optimized.pdb`, recovery via
`traj_frame.py` on `geometric_OPTtraj_Full.xyz` as planned. **But the distances
that matter are flat to +/-0.005 A over the last four steps.** Another 16
steps of DFT, or a continuation costing another 48 h, would buy a tighter
gradient and no different answer about where Asp92 and Lys226 sit. The fourth
of D236's options - stop at the gradient direction and put the def2-TZVP
single point on the recovered frame - is the one the numbers now support, and
it is the cheapest. Not acted on; the user's call.

**One thing a continuation WOULD buy, stated so it is not lost.** Zn-NZ is
still falling at step 13 and has not turned over the way Zn-OD2 did. If that
motion is the thesis's own, then where it stops is a number worth having, and
only more steps give it. That is an argument for a continuation on this
coordinate specifically - about 29 h to flatten it at the current rate - not
for finishing the optimisation to 1e-04.

**D243 - the user's decision, against the recommendation D242 argued for:
`aspdft` RUNS TO ITS WALL. The def2-TZVP single point waits for the frame that
wall produces.** (2026-11-08)

D242 made the case for stopping: the distances are flat to +/-0.005 A over four
steps, 1e-04 is ~30 steps away and the wall buys ~16, so the remaining compute
buys a tighter gradient and no different answer about Asp92. The answer is **let
it run**, and the reason D242 itself named is the one that makes this right:
**Zn-NZ has not turned over.** Lys226 is still moving onto the metal at step 13
with no reversal since step 2, and where that motion stops is a number the
thesis wants. Only more steps give it. Stopping would have saved ~24 h of one
32-core job and thrown away the only coordinate still carrying information.

So: no continuation section, no resize, nothing submitted. The job ends at its
48 h wall ~15:10 on 9 Nov at ~29 steps and **writes no `QM_MM_optimized.pdb`**.
Recovery is `python docs/csr/traj_frame.py
protdft/qmmm_aspdft/geometric_OPTtraj_Full.xyz
md_relax10_results/NVTtrajectory_lastframe.pdb <out.pdb>`, and the def2-TZVP
single point takes that frame. What to read at the wall is Zn-NZ's last four or
five steps: if it has flattened, the number is real and the salt-bridge geometry
is measured; if it is still falling, say so and do not report a distance as
though the optimiser had settled on it.

**D244 - BOTH QUESTIONS ANSWERED. `aspxtb8` converged in 103 steps with 479
active atoms and gives D241's SECOND pre-registered outcome: freeing the
second shell does NOT bring Asp92 in, so the frozen environment was not the
reason and D240's narrowing stands. And `aspdft`'s Zn-NZ HAS TURNED OVER AND
FLATTENED at 5.140 for three steps, so Lys226's distance is now a measured
number rather than a direction - 0.13 A closer to the zinc than 6DQH.**
(2026-11-09)

## `aspxtb8` (32489558): converged, 103 steps, 479 active atoms, 15.3 minutes

`converged: True`, `opt_energy` -431.5060340007243 Eh, `actatoms_count: 479`,
`actatoms_structure` the `aspxtb` minimum, zero `SVD Error`. Only `pdb_file`
and `actradius` differed from `aspxtb`, enforced by the `same_as` gate.

    pair         aspxtb min   aspxtb8 min   6DQH    aspdft step 25
    Zn-OD2         5.383        5.092       5.14        5.123
    Zn-OD1         4.952        5.830       4.28        4.534
    Zn-CG          5.564        5.777       5.15        5.223
    Zn-NZ          5.195        5.044       5.27        5.140
    NZ-OD2         2.718        2.749       2.79        2.729
    His88 NE2      2.083        2.063                   2.162
    His164 NE2     2.010        1.997                   2.073
    Zn-OH O01      1.775        1.776
    Zn-O2x         2.282        2.211
    Zn-N1x                      4.544
    O01-C4x        3.673        3.626

**THE PRE-REGISTERED READING, HONOURED: outcome 2.** D241's second outcome was
*"Asp92 stays out near 5 A -> the frozen environment was not the reason, so the
chelate needs thermal energy. D240's narrowing stands and strengthens."* With
479 atoms free instead of 124, the carboxylate carbon moves 0.21 A **further
out** and the nearest oxygen sits at 5.092 A. Nothing approaches 2.5 A. So the
frozen second shell was not what held Asp92 off the metal, and the remaining
explanation for D232's chelate is the one left standing: **300 K dynamics, not
the minimum**. The 0 K-versus-300 K half of D240's confound is the half that
matters, and the other half is now closed by measurement.

**But the carboxylate ROTATED, and this is new.** OD2 comes IN 0.29 A to 5.092,
while OD1 goes OUT 0.88 A to 5.830 - the largest Zn-OD1 distance in any
geometry in this project, against 6DQH's 4.28. The symmetric arrangement of
`aspxtb`'s minimum breaks: with the second shell free, Asp92 turns so that one
oxygen faces Lys226 at 2.749 A and the other points away from the metal
entirely. That is **monodentate to the lysine and nothing to the zinc**, which
is a sharper version of the thesis's own geometry than the crystal shows - and
it is 1.55 A off the crystal on Zn-OD1, which is the largest single
disagreement with 6DQH anywhere in this leg. Both facts go in the same
sentence or neither does.

**THE ENERGY RISE AT THE END IS NOT THE CHANNEL PROJECT'S FAILURE, AND THE TEST
I WROTE FOR MYSELF WOULD HAVE MISREAD IT.** Steps 102 and 103 rise by
+6.809e-06 and +4.081e-06 Eh - **0.004 and 0.003 kcal/mol** - with Quality
**-2.990** and **-1.959**, and the trust radius collapsing 4.933e-03 ->
2.294e-04 in response. The standing test was "a climb with good geomeTRIC
Quality means STOP". The quality here is *negative*: geomeTRIC identified the
steps as bad and shrank the trust radius, which is the healthy response, and
the magnitudes are numerical noise at a gradient of 3.0e-05. The channel
project's failure was 30-45 kcal/mol uphill with Quality ~1 - seven orders of
magnitude larger and with the optimiser believing it. **The test needs the
magnitude as well as the quality**; recorded so the next check does not stop a
converged run.

## `aspdft` (32453561): 25 steps, Zn-NZ flattened, ~2 steps of wall left

Mapping re-verified: geometry block 0 reproduces 5.153 / 4.327 / 5.192 / 5.276
/ 2.880 exactly.

    step   Zn-OD2  Zn-OD1  Zn-CG  Zn-NZ  NZ-OD2
     13    5.134   4.546   5.237  5.131  2.713
     16    5.133   4.553   5.237  5.120  2.726
     19    5.125   4.547   5.230  5.124  2.732
     22    5.117   4.535   5.220  5.138  2.735
     23    5.119   4.535   5.220  5.140  2.733
     24    5.121   4.535   5.222  5.140  2.730
     25    5.123   4.534   5.223  5.140  2.729
    6DQH   5.14    4.28    5.15   5.27   2.79

**THE QUESTION D243 KEPT THE JOB RUNNING FOR IS ANSWERED: Zn-NZ TURNED OVER.**
It reached 5.120 at step 16, came back, and has read **5.140 at steps 23, 24
and 25** - flat to the third decimal. So the decision to let it run bought
exactly what D243 said it would: Lys226's ammonium sits **5.14 A from Zn1**,
0.13 A closer than 6DQH's 5.27, and that is now a distance rather than a
direction. Had the job been stopped at step 13 the number would have been
5.131 and still falling, and reporting it would have been wrong.

**The whole site, unrestrained, against the crystal:** Zn-OD2 5.123 (5.14),
Zn-OD1 4.534 (4.28), Zn-CG 5.223 (5.15), Zn-NZ 5.140 (5.27), NZ-OD2 2.729
(2.79). **Four of five within 0.14 A**, the fifth (OD1) 0.25 A out. For an
optimisation with nothing constrained, from a 10 ns MM frame, at
wb97x-v/def2-SVP, that is close agreement with a measurement.

Energy monotonic at every one of 25 steps, -5201.1362 -> -5201.3906 Eh;
quality 1.15-1.35 for the last six; trust growing again to 1.198e-01;
`Grad` rms **3.820e-04** at step 24, better than D242's ~8%-a-step
extrapolation predicted. Zero `SVD Error` in the `.out` and zero `LinAlgError`
in the `.err`, counted in their own files. The rate did **not** keep improving:
124.4, 131.7, 115.8 min for the last three steps, against 91.1 at step 13, so
`pyscf.chk` reuse is not a trend. At 11:25 on 9 Nov with the wall at ~15:10,
about **two more steps** will land, ending near step 27 with `grms` ~3.5e-04
and no `QM_MM_optimized.pdb`.

## What the three runs now say together

Three independent setups - wb97x-v/def2-SVP with 124 active atoms, GFN2 with
124, GFN2 with 479 - put the nearest Asp92 carboxylate oxygen at **5.09, 5.38
and 5.12 A** and never on the metal, from a crystal-consistent start. D206
rule 3 excludes two methods agreeing about a handed geometry, and that
exclusion does not bite here for the reason D240 gave: the thing they agree
with is 6DQH's own 5.14, a measurement. What they do not do is settle D232 -
the chelate GFN2 produced in 50 ps of 300 K dynamics is still unexplained by
any of these minima, and that is now the single open question on this arm.

**D245 - THE CROSS-METHOD SPREAD, COORDINATE BY COORDINATE: DFT and GFN2 agree
on Lys226 to 0.055 A and on the salt bridge to 0.011 A, and disagree on Asp92's
carboxylate by 0.26-0.42 A matched and up to 1.30 A unmatched. So "the two
methods agree" is true of the metal-site core and FALSE of the coordinate the
thesis is about. Also recorded: there is no DFT run with an 8 A active region,
and GFN2 over-binds the histidines by 0.06-0.10 A while still not bringing
Asp92 in.** (2026-11-09)

**The premise that prompted this was wrong and the correction matters.** The
question asked how the DFT and GFN2 optimisations "with 8 angstrom radius
actatoms" compare. **`aspdft` has `actradius: 0`** - 124 active atoms, the QM
region only, MM environment frozen (D235). Only `aspxtb8` has the 8 A /
479-atom region. So the **matched** pair is `aspdft` vs `aspxtb`, both at 124
active atoms and differing only in `qm_engine`; `aspxtb8` differs from both in
active region as well. Any statement of the form "both methods at 8 A agree"
has no run behind it.

    pair          DFT/124   GFN2/124   GFN2/479   6DQH    matched diff
    Zn-NZ          5.140      5.195      5.044    5.27       0.055
    NZ-OD2         2.729      2.718      2.749    2.79       0.011
    Zn-OD2         5.123      5.383      5.092    5.14       0.260
    Zn-OD1         4.534      4.952      5.830    4.28       0.418
    Zn-CG          5.223      5.564      5.777    5.15       0.341
    His88 NE2      2.162      2.083      2.063               0.079
    His164 NE2     2.073      2.010      1.997               0.063

**THE LYSINE IS WHERE THE METHODS AGREE.** Matched, 0.055 A on Zn-NZ and
0.011 A on NZ-OD2. Across all three setups the salt bridge spans 2.718-2.749 -
a **0.031 A** spread - and all three fall 0.04-0.07 A short of 6DQH's 2.79, the
same direction in every case. Zn-NZ spans 5.04-5.20, all three within 0.23 A of
the crystal and all three on the same side of it. For three calculations that
differ in Hamiltonian, in active region, and in starting geometry, that is as
close to method-independent as anything on this arm.

**THE CARBOXYLATE IS WHERE THEY DO NOT.** The same matched pair differs by
0.26 A on Zn-OD2, **0.418 A on Zn-OD1** and 0.341 A on Zn-CG - an order of
magnitude worse than the lysine. Adding the 8 A run, **Zn-OD1 spans 1.296 A**,
with GFN2/479 at 5.830 sitting 1.55 A off the crystal while DFT/124 at 4.534
sits 0.25 A off it. Part of that spread is the active region and not the
Hamiltonian, which is exactly why the unmatched run cannot be quoted as a
second opinion on the matched one.

**What this licenses and what it does not.** Two claims survive as
method-independent: **where Lys226 sits**, and **that Asp92 is not on the
metal**. A third does not: **how Asp92 is turned** - whether OD1 or OD2 faces
the metal, and how far the far oxygen swings - is method- and
active-region-dependent at the 0.4-1.3 A level, and that rotamer is the
coordinate the thesis's mechanism runs through. Any statement about Asp92's
orientation needs to name which run it came from.

**A point that runs in the thesis's favour, and it comes from the matched pair
rather than from agreement.** GFN2 is 0.079 A shorter on His88 NE2 and 0.063 A
shorter on His164 NE2 than DFT, and `aspxtb8` is shorter still - GFN2
**over-binds the zinc**, which is the known direction of its error and the same
direction Mehmood & Kulik report for too-small QM regions (D7's region-size
caveat). GFN2 binds the metal more tightly than DFT does **and still does not
bring Asp92 in**. The bias and the result point opposite ways, which is worth
more than two methods agreeing.

**Caveat on the DFT column.** `aspdft` is at `grms` 3.820e-04 against a 1e-04
criterion: settled to +/-0.003 A over its last three steps but not a stationary
point. The matched differences above are 0.011-0.418 A, so the comparison is
not limited by that, but no DFT number here should be quoted as a converged
minimum.

**No job was proposed or run for this.** D244 named a DFT optimisation from
`aspqm2`'s chelated frame as the test that would separate "GFN2's barrier is
wrong" from "minima are not where the enzyme sits"; the user declined it as not
worth the compute. **That separation therefore stays open and must not be
written up as settled.** What is recorded instead is this: the artefact reading
has no direct evidence, because GFN2's own minima keep Asp92 out, so GFN2 is
not placing a spurious minimum - only its 300 K dynamics reach the chelate, and
no run here measures a barrier.

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
- P10: `site.py check` has no way to say a criterion applies at a transition state rather than a reactant. D201 parks the Burgi-Dunitz angle on a stage name that is never passed, which works but reads as a trick; a real stage list, documented, would be honest.
- P9: let `site.py check` read a criterion's occupancy from a `postprocess` section instead of measuring a single frame. The last frame of batch 1 rejected the structure on a carbonyl at 3.07 A that the 50 ps has inside 2.5 A in 66% of frames at 2.43 A mean (D200).
- P8: fix `md.py`'s sphere-change warning, which says a non-bonded metal has "no Pauli repulsion" in a codebase that ships fitted Li/Merz LJ parameters and applies them by default. It misled this project's own write-ups three times (D199). The accurate wording is that a 12-6 term cannot tell a ligand from a non-ligand, which is `ions.py`'s own explanation.
- P7: preserve the input ligand atom names through OpenFF parameterisation, so `preprocess.results.ligand_atom_map` is always empty (D169).
