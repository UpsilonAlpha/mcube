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

**D5 – Starting structures come from unknown earlier runs.** (open)
- Evidence: in `ligand_lastframe.pdb` the hydroxide O–H was stretched to 1.54 Å, and H03 sat 2.83 Å from Asp92 OD1. That suggests it came from an earlier QM/MM run caught mid-proton-transfer, whose history isn't recorded. `hydrolysed_lastframe.pdb` was built so its atom order matches the reactant's.
- Decision: use them as they are, but re-equilibrate with MM before any QM/MM (D14) and select snapshots from that MD (D16).
- Revisit: **yes**. The write-up should state how the product structure was generated.

## Force field and protonation

**D6 – Force fields.** (inherited, 2026-10-15)
- Protein: Amber ff14SB. Water: TIP3P.
- Ligand: OpenFF (ASH's `small_molecule_parameterizer` default, `openff-2.0.0.offxml`).
- Hydroxide: custom `hydroxide.xml`, with charges O −1.02 / H +0.02, an O–H of 0.975 Å, and the LJ parameters listed in that file.
- Zn²⁺: nonbonded ion, with geometry held by distance constraints (D15).
- Revisit: **yes**. The hydroxide parameters have no stated source. A nonbonded Zn model with constraints is crude; consider a bonded model (MCPB.py) or 12-6-4 LJ for production MD.

**D7 – Residue protonation.** (inherited)
- His88 and His164 are HID (`residue_variants`). Everything else comes from pdbfixer at pH 8.0 (mcube's explicit-solvent default).
- Measured: Lys226 is protonated (3 HZ). Asp92 is deprotonated.
- Revisit: **yes**. The pH 8.0 default and the His tautomers haven't been checked; a pKa prediction (e.g. PROPKA) would help.

**D13 – Meropenem is modelled as the zwitterion.** (2026-10-15, user decision)
- Evidence: the ring carboxylate pKa is about 3 and the pyrrolidine amine pKa about 7.5–9, so at pH 7 the molecule is mostly COO⁻ / NH₂⁺. The earlier model had a neutral COOH and a neutral amine. Metallo-β-lactamases typically bind the C3 carboxylate with a Lys, and Lys226 may be that residue (NZ is 5.4 Å from C6x in the current frame).
- Decision: move H01 from O3x to N2x in both the reactant and the product, so atom order is unchanged. Re-parameterise with zwitterion SMILES. Reactant ligand net charge 0; hydrolysed product −1. The QM charge is unchanged.
- Alternatives: neutral, as before; anionic (COO⁻ with a neutral amine).
- Revisit: possibly, as a sensitivity check.

## Metal site

**D8 – Zn coordination cutoff of 3.0 Å.** (2026-10-15, user decision)
- Evidence: typical Zn–ligand bonds are 2.0–2.3 Å, and coordination cutoffs are usually 2.5–3.0 Å. In this frame His88 NE2 is at 2.77 Å and Arg90 NH1 at 3.50 Å. The old 4 Å cutoff wrongly constrained Arg90 to the Zn.
- Result: the coordinating atoms are OH⁻ O01 (2.34), Asp92 OD2 (2.40; OD1 at 2.60), meropenem carbonyl O2x (2.55) and His88 NE2 (2.77).
- Revisit: no.

**D15 – During MD, Zn–ligand distances are fixed except the substrate carbonyl's.** (2026-10-15, user decision)
- Evidence: fixing every Zn–ligand distance pins the substrate carbonyl at 2.55 Å and biases the Michaelis complex.
- Decision: keep the Zn–OH⁻, Zn–Asp92 and Zn–His88 distances fixed. Leave Zn–O2x free (`constraint_exclude_residues: {C: [1]}`).
- Revisit: **yes**. Check whether the carbonyl stays near the Zn. A free hydroxide or a bonded model would be alternatives.

## QM/MM model

**D9 – QM method: GFN2-xTB with electrostatic embedding.**
- Run through ASH. QM/MM boundaries cut at CA–CB with link atoms. Non-amino-acid residues (Zn, ligand, OH) are fully QM, with no link atoms.
- Evidence of concern: in the first reactant optimisation, Zn–OH shortened to 1.75 Å, against the usual ~1.9–2.0 Å.
- Revisit: **yes, required**. Benchmark against DFT on a cluster model, and recompute DFT/MM energies along the xTB paths (P5).

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
- Decision: `extra_qm_residues: {A: [90, 226]}`, about 105 QM atoms.
- Revisit: after sampling, check which other residues contact the reacting atoms (Asn, other Lys/His, waters).

**D11 – Active region of 8 Å around the Zn.**
- Built from atom 4203 (431 atoms in the first test). The same active atom list is used for every endpoint, intermediate and NEB (`--actatoms`).
- Revisit: after the QM region change, check that the converged results don't depend on the radius.

## Mechanism and path protocol

**D16 – Snapshots must be lined up for attack.** (2026-10-15, user decision)
- Evidence: in the current reactant, O01–C4x is 3.81 Å and the O01–C4x–O2x angle is 51°. The hydroxide lies beside the carbonyl O rather than on the Bürgi–Dunitz trajectory.
- Decision: minimise, warm up, equilibrate at constant pressure, then run 2 ns at constant volume (4 fs steps), saving a frame every 2 ps. Keep frames with O01–C4x ≤ 3.2 Å and O01–C4x–O2x at 95–120°. If fewer than 5 qualify, consult before adding restraints.
- Alternative: a flat-bottom restraint, which is biased.
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
