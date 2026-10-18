# CSR-1 meropenem hydrolysis: QM/MM pathway plan

How we go from the restrained MM model to reaction barriers. The decisions behind each step are in [csr-decision-log.md](csr-decision-log.md).

## Questions

CSR-1 is an unusual B3 with an RQK-HRH active site. The paths must answer:
1. **Lys226 (K254):** does it stabilise meropenem hydrolysis, as Lys92 does in EcoRV?
2. **Arg90 (R118) and Arg229:** are the arginines involved?
3. **Asp92 (D120), the "bridging" aspartate:** is it involved?

How each question is answered:
- **Who moves which proton:** along every NEB path, monitor Lys226 NZ, Asp92 OD1/OD2, the Zn water and O01 protons.
- **Per-residue contributions** (as Imhof did for Lys92 in EcoRV): recompute each path's energies with the charges of Lys226, Arg90, Arg229 and Asp92 zeroed in turn, then optimise the path again with the residue removed or neutralised for the important cases.
- **Where they sit at the start (D29):**
  - Arg90 and Arg229 hold the 6α-hydroxyethyl and C2 side chain;
  - Asp92 binds the Zn and H-bonds both water oxygens;
  - Lys226 sits 2.7 Å from the O01 nucleophile.

## Methods we replicate (checked against the abstracts, 2026-10-18)

**Imhof, Fischer & Smith, *Biochemistry* 2009 (EcoRV; PMID 19678693)**
- **Approach:** AM1/d/MM. They computed a whole *network* of reaction pathways exhaustively, rather than one guessed path. The full network method is reviewed in Imhof, *Methods Enzymol.* 2016 (PMID 27497170).
- **Validation:** AM1/d against B3LYP/6-31++G(d,p) on small model systems.
- **Mechanism they found:**
  - a **water nucleophile**, whose proton goes to a carboxylate oxygen of **Asp90**;
  - **Mg²⁺ activates the water**;
  - **Lys92 positions the phosphate and stabilises the intermediate**;
  - the protein electrostatically stabilises the transition states.

**Xu, Xie & Guo, *J. Biol. Chem.* 2006 (mono-Zn B2 CphA, biapenem; PMID 16423823)**, and **Wu, Xu & Guo, *JACS* 2010 (PMID 21138257)**
- **Approach:** SCC-DFTB/MM MD, with free energies from potential-of-mean-force (umbrella sampling), plus DFT.
- **Mechanism they found:**
  - a **non-metal-bound water** nucleophile;
  - the **general base is the Zn-bound Asp120**, not His118;
  - a single, near-concerted addition/elimination transition state, with a PMF barrier of ≈14 kcal/mol;
  - protonating Asp120 breaks its O–Zn bond, and the leaving N replaces it as the Zn ligand.
- The 2010 paper found that the crystallised enzyme–intermediate complex lies on a *minor* pathway.

**Other MBLs, same group:**
- NDM-1 (B1, PMID 24025144) and SMB-1 (B3, PMID 32146530): the *Zn-bridging hydroxide* is the nucleophile, with Asp120/124 assisting. That is the di-Zn case.

**Relevance to CSR-1:**
- Our **Asp92 is D120** (crystal 120 = model 92). As in CphA, it binds the Zn, and in our QM/MM MD it H-bonds both water oxygens.
- CSR-1 has one Zn, like CphA. So the precedent favours **water + Asp general base** over a pre-formed hydroxide.

## Replication protocol for CSR-1 with modern methods

| | Imhof 2009 / Xu–Guo 2006 | CSR-1 (mcube) |
| --- | --- | --- |
| QM level for sampling and paths | AM1/d (Imhof); SCC-DFTB (Guo) | GFN2-xTB/MM (electrostatic embedding) |
| Energies | B3LYP model systems; DFT corrections | ωB97X-D3(BJ)/def2-SVP/MM single points along paths; def2-TZVPD on stationary points (D31) |
| Reactant sampling | minimised structures (Imhof); QM/MM MD (Guo) | 50 ps QM/MM MD (D28/D29), several snapshots |
| Paths | exhaustive path network (Imhof) | relaxed scans → NEB-CI → TS optimisation + IRC for every candidate step, assembled into a network |
| Free energies | PMF by umbrella sampling (Guo) | xTB/MM umbrella sampling along the best path's coordinate(s), with WHAM/MBAR. Needs a new mcube stage. |
| Validation | AM1/d vs B3LYP on model systems | xTB vs ωB97X on the QM/MM structures (D30), plus a cluster-model check |
| Residue roles | Lys92 positioning and stabilisation | charge deletion / mutation for Lys226, Arg90, Arg229, Asp92 along the best paths |

**Step 1 – enumerate mechanisms (Imhof-style network).** Candidate first steps, each a scan → NEB → TS/IRC on 2–3 snapshots:
- **H1 (CphA-like):** a **non-metal-bound water** attacks C4x and **Asp92** accepts its proton. This needs the second-shell waters near C4x, which we haven't analysed yet: the QM/MM MD analysis covered only O01 and HOH485.
- **H2:** the Zn-bound water (O01, after the Lys226 → O01 transfer) attacks, and Asp92 or Lys226 NZ accepts the proton.
- **H3 (classic B3):** Zn-bound OH⁻ attacks, with Lys226 protonated. This depends on D30 confirming that this protonation state is reasonable.
- **H4:** Lys226 acts as the base.

Then, from the tetrahedral intermediate or directly (concerted, as in CphA): C4x–N1x cleavage, and N1x protonation by Asp92-H, Lys226 or water. The network is the set of intermediates joined by the computed steps. The answer is the path with the lowest highest barrier, plus any competitive alternatives (cf. the "minor pathway" result of Wu *et al.*).

**Step 1 results – nucleophile candidates in the 50 ps QM/MM MD** (`scripts/nucleophile_candidates.py`, 1,000 frames; `qmmm_analysis/nucleophiles*`).

Criteria: attack = d(O···C4x) ≤ 3.5 Å and O–C4x=O2x angle 95–120°; loose attack = ≤ 4.0 Å and 90–125°; H-bond = O···acceptor ≤ 3.2 Å with an H ≤ 2.5 Å.

| Candidate | In QM? | Zn-bound | d(O···C4x) min | Angle | Loose attack frames | H-bond partners (occupancy) |
| --- | --- | --- | --- | --- | --- | --- |
| O01 (former OH⁻ + Lys226 H) | yes | 100% | 2.91 Å | ≤ 87° | 0 | meropenem C3-carboxylate O4x 88%; Lys226 NZ 2.7 Å, Asp92 OD2 2.9 Å (heavy atoms) |
| HOH485 (Zn water) | yes | 100% | 3.20 Å | ≤ 87° | 0 | **Asp92 OD1 97%**, OD2 45% |
| **HOH E3094** (second shell) | MM | 0% | 3.60 Å | mean 105° | 4.7% (47 frames) | meropenem O1x (6α-OH, held by Arg90 NH1 99%) 97%, HOH E2493 99%, E2155 94%, **E3166 94%** |
| HOH E3166 | MM | 0% | 3.32 Å | ≤ 75° | 0 | **Asp92 OD1 67%**, HOH485 (mean 3.1 Å) |
| HOH E2493 | MM | 0% | 3.93 Å | mean 123° | 0.7% | E3094 |

- The Zn-bound oxygens are close to C4x but beside the carbonyl O; they never reach the attack angle.
- Only the second-shell water **E3094** sits on the Bürgi–Dunitz trajectory. Its route to a base is the relay **E3094 → E3166 → Asp92 OD1**.

**H1–H4, revised with these data.** Each is a relaxed xTB/MM scan of d(O_nuc–C4x) from its starting value to 1.45 Å (0.1 Å steps) from 3 frames, with protons free. Monitored throughout: O_nuc–H···acceptor, C4x–N1x, Zn–O2x and Zn–O_nuc.

| | Nucleophile → proton acceptor | Precedent | QM region change | Start frames |
| --- | --- | --- | --- | --- |
| **H1** | second-shell **E3094 → (E3166) → Asp92 OD1** | CphA (Xu–Guo 2006): non-metal water, Zn-bound D120 base | **add E3094 and E3166** (the scans need them in QM) | E3094 loose-attack frames spread in time, e.g. 105, 280, 472 |
| **H2** | Zn-bound **HOH485 → Asp92 OD1** (direct, 97% H-bond) | EcoRV (water → Asp90, metal-activated) | none | frames with the best HOH485 angle |
| **H3** | Zn-bound **O01 (water) → meropenem C3-carboxylate O4x** (88% H-bond), or → Lys226 NZ | substrate-assisted / Lys as base (Q1) | none | frames 575, 244, 192 (closest O01 approach) |
| **H4** | classic **Zn-OH⁻** (proton held on Lys226) attacks; Lys226-NH3⁺ stabilises | di-Zn B1/B3 bridging hydroxide | none (Lys226 NZ–H constrained during the scan) | same as H3; **only if D30 says Lys-NH3⁺/OH⁻ is reasonable at DFT** |

Then:
- If a scan's end point is a tetrahedral intermediate, optimise it unconstrained. Run NEB R→TI and TS/IRC, then the C4x–N1x cleavage and N1x-protonation steps (from Asp92-H, Lys226, water or the carboxylate).
- If the attack and ring opening are concerted, as in CphA, a single NEB R→EI⁻ covers both.
- The paths connect into the Imhof-style network.

**Step 2 – DFT energies along every path** (ωB97X-D3(BJ)/def2-SVP; TZVPD at stationary points). Paths whose ordering changes between xTB and DFT are re-examined.

**Step 3 – free energy of the winning mechanism (Guo-style).**
- QM/MM umbrella sampling along the reaction coordinate(s), e.g. d(O–C4x) − d(C4x–N1x), and a proton-transfer coordinate to Asp92.
- PMF by WHAM/MBAR, with DFT corrections.
- At ~1.4 s per xTB/MM MD step, a 30-window, 5 ps-per-window PMF is ≈ 60 h of compute per coordinate, run in parallel.

**Step 4 – residue roles (the three questions).**
- **Q1:** compare paths with Lys226 protonated vs neutral; charge deletion along the best path; whether Lys226 stabilises the oxyanion or the leaving N.
- **Q2:** zero the MM charges of Arg229 (it is in the MM region) along the paths; for Arg90 (QM) use a neutralised or truncated QM model.
- **Q3:** Asp92 as base (H1/H2) vs a spectator; protonated-Asp92 intermediates; whether Asp92 leaves the Zn when protonated (CphA behaviour).

**mcube work this needs:**
- an umbrella-sampling stage plus WHAM/MBAR analysis;
- a charge-deletion single point (zero chosen MM residues' charges);
- second-shell water / nucleophile candidate analysis (`contact_report.py` and `site_report.py` cover most of it);
- network bookkeeping (one config section per step; a small script to assemble barriers).

## Status (2026-10-17)

- **QM/MM MD:** 50 ps (job 29937641; first 39.5 ps analysed in D29). Lys226 protonates the Zn-OH⁻ at once, and no frame reaches attack geometry.
- **DFT check:** the protonation-state check (D30, job 29981678) decides how reactants are built.
- **Then:** relaxed O01–C4x scans from QM/MM frames, with protons free and monitored, give TI endpoints for NEB R→TI, then TS/IRC, per snapshot.

## Principle

**The metal site is never left free in MM.** Every classical stage keeps the Zn coordination geometry held by constraints and restraints. That is the purpose of the mcube pipeline. The Zn site is only free in QM/MM, where the Zn, its ligands and the substrate are all in the QM region and no constraints apply.

So the sequence is:

**preprocess → restrained MM equilibration/NVT → QM/MM** (optimisation, scans, NEB, TS/IRC)

There is no unrestrained MM MD at any point.

## System

`/scratch/user/uquadhik/CSR_crystal`. Model numbering = crystal numbering − 28. Atom indices are 0-based in `preprocess_results/finalsystem.pdb`.

| Item | Value |
| --- | --- |
| Zn ligands | His88, His164 (NE2), OH⁻ (O01), Asp92 OD2, meropenem β-lactam O2x (pulled in) |
| Key atoms | Zn 4203, O2x 4209, C4x 4208, OH⁻ O01 4255 (H03 4256), Zn-bound water O 4257 |
| QM region | 119 atoms: Zn, OH⁻, water, meropenem, His88, His164, Asp92, Arg90, Lys226; charge +2, multiplicity 1 |
| QM method | GFN2-xTB with electrostatic embedding (exploration). DFT/MM for final energies. |
| Active region | 8 Å around the Zn (atom 4203). One `active_atoms` list per snapshot, reused by every later step for that snapshot. |

## Steps

### 1. Restrained MM equilibration and NVT (running: job 29848762)

- mcube md stage with:
  - **restraints:** Zn–O2x 2.2 Å and OH⁻–C4x 3.0 Å, k = 20 kcal/mol/Å²;
  - **constraints:** Zn to His88, His164, OH⁻ and Asp92.
- Sequence: minimisation, gentle warm-up, NPT, then 1 ns NVT.
- Checks:
  - Zn–O2x reaches ~2.2 Å;
  - OH⁻–C4x ~3.0 Å with an attack angle of 95–120°;
  - the Zn-bound water (HOH485) leaves the Zn;
  - the protein and the rest of the ligand stay stable (RMSD).

### 2. Snapshot selection (from the restrained NVT)

- Use `scripts/select_nac_frames.py` on `md_restrained/NVTtrajectory_imaged.dcd`.
- Criteria: OH⁻–C4x ≤ 3.2 Å, attack angle 95–120°, Zn–O2x ≤ 2.5 Å, frames ≥ 50 apart (100 ps). Take ≥ 5 snapshots.
- The frames carry the restraint bias by design. Step 3 removes it: the QM/MM optimisation has no restraints.

### 3. QM/MM reactant per snapshot (mcube qmmm)

```bash
mcube qmmm config.yaml --pdb nac_frameXXXXX.pdb --section qmmm_snapN
```

- Settings:
  - `skip_equilibrate: true` (no QM/MM MD), `singlepoint: false`;
  - `actradius: 8`, `actregion_origin: 4203`.
- Outputs: `QM_MM_optimized.pdb` (the reactant minimum, R) and `active_atoms`.
- Check: does the carbonyl stay on the Zn, and the OH⁻ near C4x, once the restraints are gone?

### 4. Intermediates and product by relaxed QM/MM scans (new mcube option)

All structures are built inside the same box from R, so every endpoint has identical atoms.

| Step | Coordinate scanned | Result |
| --- | --- | --- |
| 4a | O01–C4x, ~3.0 → 1.45 Å | Tetrahedral intermediate (TI). Optimise it without constraints: if it isn't a minimum, attack and ring opening are **concerted**. |
| 4b | C4x–N1x, 1.4 → ~2.6 Å | Ring-opened anionic intermediate (EI⁻), with the negative charge on N. |
| 4c | Proton donor–H…N1x | Δ2 product (P). One scan per candidate donor. |

Candidate proton donors for 4c (decision log D18):
- H03 from the original OH⁻, directly or relayed through Asp92;
- Lys226 NZ–H;
- a water molecule, if one reaches N1x.

The Δ1 tautomer (protonation at C2) comes later (D19).

The scans give a first energy profile, show which steps are separate, and supply the endpoints for NEB.

### 5. NEB per elementary step (mcube neb)

- Climbing-image NEB between optimised minima: R→TI, TI→EI⁻, EI⁻→P (or R→EI⁻ if step 4a shows a concerted path).
- `neb.inputs`:
  - `images` (default 12);
  - CI on;
  - convergence settings;
  - the same `--actatoms` as the reactant;
  - xTB `runmode: library`;
  - cores split between images and xTB threads.

### 6. TS refinement and IRC

See "TS refinement and IRC" below. Run on each climbing image:
1. A transition-state optimisation with a QM-region Hessian; check for exactly one imaginary frequency along the reaction coordinate.
2. An IRC in both directions; check it ends in the intended minima.

### 7. Energies

- DFT/MM single points on the optimised minima, the transition states and the NEB images:
  - ωB97M-V/def2-SVP/MM, then ωB97M-V/def2-TZVP(D)/MM;
  - PySCF, with gpu4pyscf on Bunya's NVIDIA H100/A100 nodes (D9).
- If GFN2 and DFT disagree badly, re-optimise the key stationary points at DFT/MM.
- Check GFN2-xTB on a cluster model of the Zn site: Zn, its ligands' side chains, OH⁻ and a β-lactam fragment.

### 8. Sampling

- Repeat steps 3–6 for ≥ 5 snapshots, at least for the rate-determining step.
- Report barriers as mean ± spread.

### 9. Analysis

- **Per-residue contributions:** zero the charges of Lys226, Arg90 and Asp92 along each path (as Imhof et al. did for Lys92 in EcoRV).
- **Asp92:** does it H-bond the OH⁻ or coordinate the Zn (D24)?
- **Later:** one Zn vs two Zn (P3), and the Δ1 product (P4).

## mcube support for steps 3–6

Implemented (decision log D26; see the README for every key):

1. **Constrained QM/MM optimisation:** `qmmm.inputs.opt_constraints`, `[[i, j]]` (hold) or `[[i, j, d]]` (set).
2. **Relaxed scan:** `mcube scan config.yaml --pdb R.pdb --section scan_attack --actatoms active_atoms`, with `atoms`, `start`, `end` and `points` in the section. Writes `scan_XX.pdb` and `scan_energies.csv`.
3. **NEB settings from the config:** `neb.inputs` (or a per-step section given with `--section`): `images` (12), `CI`, `maxiter`, tolerances, interpolation, springs, `method: nebts`. Writes `saddlepoint.pdb`.
4. **TS refinement and IRC:** `mcube ts config.yaml --pdb saddlepoint.pdb --section ts_attack --actatoms active_atoms`. Runs the TS optimisation with a partial Hessian over `hessian_atoms`, a frequency check, and a quasi-IRC along the imaginary mode (or Sella), with `monitor` distances reported.

Still to do: **DFT/MM single points along paths**, PySCF on GPU (needs gpu4pyscf in a separate env).

Example for step 4a on snapshot 1:

```yaml
scan_attack_s1:
  inputs: {atoms: [4255, 4208], start: 3.0, end: 1.45, points: 16}
ts_attack_s1:
  inputs:
    hessian_atoms: [4203, 4208, 4209, 4255, 4256]   # Zn, C4x, O2x, O01, H03; add N1x and neighbours (look up indices)
    monitor: [[4255, 4208], [4203, 4209]]   # OH–C4x, Zn–O2x
```

## Cost estimate (to be measured)

- GFN2-xTB QM/MM gradient with 119 QM atoms: ~1–2 min. That is extrapolated from ~60–70 s per step with 71 QM atoms (`running-on-bunya.md`).
- One NEB (12 images × ~100 iterations) ≈ 1 day on 16–32 cores.
- The full step-by-step path for one snapshot ≈ 3–4 days.
- Run the rate-determining step on all snapshots and the other steps on 1–2.

## TS refinement and IRC

**TS refinement.**
- NEB represents the path as a chain of images. The climbing image is pushed uphill towards the saddle point, but it only converges to the NEB tolerance on a discretised path, so it is an approximate transition state.
- TS refinement starts from the climbing image and runs a saddle-point optimisation, which follows the one direction of negative curvature uphill and everything else downhill. geomeTRIC does this with partitioned rational-function steps, using a Hessian computed for the QM region.
- Success means a precisely converged first-order saddle point, with **exactly one imaginary vibrational frequency**, whose mode is the reaction coordinate (e.g. the O01···C4x bond forming).

**IRC (intrinsic reaction coordinate).**
- Starting from the refined TS, the IRC follows the steepest-descent path, in mass-weighted coordinates, downhill in both directions along the imaginary mode.
- It shows which two minima the TS actually connects, which may not be the intended endpoints: for example, a side-chain rotation instead of the chemical step.
- It also gives a clean minimum-energy path, used for energy profiles and DFT single points.
- Together, TS refinement and IRC turn "the highest NEB image" into a verified transition state for that elementary step.

## Minimal protocol (M-series) — supersedes the P-series below

Decided 2026-10-18 (D47): CPU only, and as little simulation as will test the three
questions, staying with what the source papers actually did.

**This reverses the tier B decision (D44).** Re-reading the methods settles it: both
source papers used *semi-empirical* QM/MM as the primary method, with DFT only for
validation on small models.

| | Imhof 2009 (EcoRV) | Xu & Guo 2006 (CphA) | CSR-1, here |
| --- | --- | --- | --- |
| Primary QM | AM1/d/MM | SCC-DFTB/MM | **GFN2-xTB/MM** |
| DFT used for | B3LYP on small model systems, as validation | DFT corrections | cluster-model validation + single points on the final path |
| Stationary points optimised at DFT/MM? | no | no | **no** |

GFN2-xTB is a later, better-parameterised method than either AM1/d or SCC-DFTB, so
using it as the primary method is faithful to these papers rather than a compromise.
Re-optimising every stationary point at DFT/MM (tier B) went beyond both of them, and
its cost is what forced the GPU detour. `qm_engine: dft` stays available in mcube for
a final check on one structure, but is not the route.

### M1. Settle the reactant - in flight

Job 30053828 (D42): proton on Lys226 versus on O01, in a clean four-coordinate
sphere, xTB/MM geometries with DFT/MM single points. One number decides whether we
have a hydroxide. Nothing else starts until it lands.

### M2. One mechanism, one structure

The precedent from both papers and our own D33/D41 geometry point the same way: a
**non-metal-bound water attacking C4x with Asp92 (= Asp120) as the general base**,
exactly the CphA case. Compute that one first, on the single best reactant:

scan (2D: O-C4x x H...Asp92 OD1) -> NEB-CI -> TS optimisation -> IRC -> DFT/MM single
points along the path. Everything xTB/MM except the last step.

Only if it fails do the alternatives get the same treatment, in this order: Zn-bound
water with Asp92 as base, then concerted N1x protonation.

### M3. Answer the three questions by charge deletion

The cheap part, and what the questions actually require. On the one converged path,
recompute the energies with selected residues' charges zeroed:

| Question | Residue | Region |
| --- | --- | --- |
| 1. Does Lys226 stabilise hydrolysis, as Lys92 does in EcoRV? | Lys226 | QM |
| 2. Are the arginines involved? | Arg90, Arg229 | Arg90 QM, Arg229 MM |
| 3. Is the bridging Asp92 involved? | Asp92 | QM (and M2 tests it as the base directly) |

Imhof's EcoRV conclusion - Mg(II) activates the water, Asp90 takes the proton, Lys92
positions the substrate and stabilises the intermediate - is the direct comparison for
question 1.

### M4. Validation, as the papers did it

One cluster model (QM region only, no MM), xTB versus wB97X-D3(BJ)/def2-TZVPD, for the
rate-determining step. This is Imhof's AM1/d-vs-B3LYP check, not a QM-region
convergence sweep.

### Deliberately not doing

Dropped from the P-series as more than the questions need: DFT/MM optimisation and
frequencies of every stationary point; QM-region size convergence sweeps; three
structures per mechanism; umbrella-sampling free energies (Guo's PMF - only if a
barrier must be compared with k_cat); and the 6DR8 triple-mutant control. Each is
defensible, none is necessary to answer the three questions, and any can be added
later if a referee asks.

## Proposed protocol (P-series) — publication grade

Proposed 2026-10-18 (D43), awaiting a decision on the level of theory. The point is
that the three questions (Lys226, the arginines, Asp92) are answered by *differences*
between computed paths, so every path must be computed to the same, defensible standard.

### P-A. Establish the Michaelis complex before any path

The degraded reactant (D41) cost us every scan so far, so this is settled first and
stated explicitly in the paper. Candidate states, all from the same NVT frame, all
constrained identically, compared on energy at DFT/MM:

| State | Zn sphere | Proton |
| --- | --- | --- |
| 1 | 2 His + OH- + H2O | Lys226-NH3+ |
| 2 | 2 His + H2O + H2O | Lys226 neutral |
| 3 | 2 His + OH- (no HOH485) | Lys226-NH3+ |
| 4 | 2 His + OH- + H2O + Asp92 OD2 | Lys226-NH3+ |

States 1 and 2 are job 30053828 (D42). States 3 and 4 test whether HOH485 belongs in
the sphere and whether Asp92 is a ligand or a base — the question D41 identified as
defining the mechanism. The chosen state is the one lowest in energy, not the one the
setup happened to produce. Report the carbonyl distance as an *outcome*, not a restraint:
re-run the winner with no Zn-O2x constraint to see where it settles.

### P-B. Level of theory (decision required)

| Tier | Paths | Stationary points | Comment |
| --- | --- | --- | --- |
| A | xTB/MM | xTB/MM + DFT/MM single points | What we have been doing. A referee will object that barriers from semi-empirical geometries are not quantitative. |
| B | xTB/MM (scans, NEB) | **re-optimised at DFT/MM**, DFT/MM frequencies and IRC | Standard for MBL QM/MM papers. **Recommended.** |
| C | DFT/MM throughout | DFT/MM | Defensible but likely unaffordable for 3 mechanisms x 3 structures. |

Under tier B, xTB/MM is used only to find approximate paths; every number that appears
in the paper comes from DFT/MM. Final energies on stationary points at
wB97X-D3(BJ)/def2-TZVPD (D31).

Tier B was chosen (D44). mcube gained `qm_engine: dft` for this; the cost figures and
the three ways to make it affordable are in D44, and the level split below follows from them:

| Purpose | Level |
| --- | --- |
| Paths (scans, NEB) | GFN2-xTB/MM |
| Stationary-point geometries and frequencies | DFT/MM, meta-GGA + dispersion (e.g. r2SCAN-D4/def2-SVP), `dft_*` keys |
| All reported energies | wB97X-D3(BJ)/def2-TZVPD on those geometries, `sp_*` keys |

Frequencies use `partial_hessian_atoms` on the reacting subset: a full 119-atom numerical
Hessian is 714 gradients and is not affordable at any DFT level we can reach.

### P-C. Convergence and validation (both expected by referees)

- **QM region size**: R and TS at ~119, ~150 and ~200 atoms; show the barrier is converged. Also vary the active-region radius (8 -> 10 A).
- **xTB vs DFT**: for each elementary step, compare on a cluster model (P1). Report the deviation rather than assuming xTB geometries are adequate.
- **Basis/dispersion**: def2-SVP -> def2-TZVPD effect on the barrier; D3(BJ) contribution reported separately.

### P-D. Free energies, not just potential energies

Barriers compared with k_cat need DG. Under tier B, DFT/MM frequencies on the QM region
give ZPE and thermal corrections almost free. For the rate-determining step only, add
umbrella sampling or metadynamics at xTB/MM along the reaction coordinate to capture the
condensed-phase contribution, and state clearly which terms are included.

### P-E. Sampling

At least three independent reactant structures per mechanism, with the spread in barriers
reported. One path from one frame is not a result.

### P-F. Mechanism discrimination

All three hypotheses computed to the same level and compared, not just the favoured one:
Zn-OH- attack; second-shell water with Asp92 as general base (CphA-like, Xu & Guo 2006);
concerted attack with N1x protonation. The comparison *is* the answer to question 3.

### P-G. Residue contributions — the three questions

On the converged paths: charge deletion / in-silico mutation for **Lys226** (question 1),
**Arg90 and Arg229** (question 2) and **Asp92** (question 3). Arg229 in MM, Arg90 and
Asp92 in QM (P2).

### P-H. An experimental control: the triple mutant

**PDB 6DR8** (1.48 A) is the R118H/Q121H/K254H triple mutant and, unlike the wild type,
**is crystallographically metal-loaded with two zincs**. It gives us (i) a structure with
an experimentally defined zinc site to validate the whole protocol against, and (ii) the
natural comparison: the canonical di-zinc site versus the RQK variant. Computing the same
mechanism for both turns "what does CSR-1 do?" into "why does the RQK substitution change
the mechanism?", which is a stronger result and uses an existing structure.

### P-I. Limitations to state explicitly

- No wild-type CSR-1 structure retains zinc (6DN4 is apo), so the mono-zinc Michaelis
  complex is constructed; its justification is P-A, and the zinc position comes from the
  Q121H/K254H mutant 6NC5, at a site whose own ligands are wild type (D41).
- The non-truncated form is modelled, as the more faithful model of the natural enzyme;
  residues 1-14 have no experimental coordinates and sit 9-13 A from the zinc.
- Semi-empirical sampling underlies the paths even under tier B; P-C quantifies that.
