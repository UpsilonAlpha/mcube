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
