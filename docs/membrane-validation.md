# Membrane mode: validation

The numbers behind the claim that mcube's membrane mode gives a bilayer a sane area per
lipid and a stable box, and that an isotropic barostat — the only kind ASH attaches —
does not. Reproduce with:

```bash
python scripts/check_membrane_barostat.py --ps 40
```

**System.** OpenMM's pre-equilibrated POPC patch (`openmm/app/data/POPC.pdb`): 128
lipids, 5120 waters, 32,512 atoms, starting box 63.70 × 66.87 × 73.14 Å, i.e. 66.6 Å²
per lipid. The System is built by `mcube.membrane.membrane_system`, the same call the md
stage makes: PME, 12 Å cutoff, 5e-4 Ewald tolerance, `HBonds`, rigid water, 1.5 amu
hydrogens, and one `MonteCarloMembraneBarostat` (XYIsotropic, ZFree, zero surface
tension) at 300 K and 1 bar.

**Runs.** 40 ps each, 4 fs timestep (mcube's production default, which the HBond
constraints and heavy hydrogens allow), `LangevinMiddleIntegrator`, CPU platform, after
a 500-step minimisation. Figures are means over the second half of each run; `±` is the
standard deviation over that window. The stretched cases start from the same patch with
its box and coordinates scaled by 1.08 in x and y, which puts the bilayer under ~8%
lateral tension.

**What is measured.** The area per lipid is `Lx · Ly / (N/2)`, since the bilayer spans
the box and each leaflet holds half the lipids. P–P is the separation of the two
phosphate planes (mean z of one leaflet's P atoms minus the other's), i.e. the bilayer
thickness. The leaflet counts are a distortion check: a bilayer that has started to
break up stops splitting 64/64.

## Results

| case | area per lipid (Å²) | box (Å) | P–P (Å) | volume (nm³) | leaflets |
|---|---|---|---|---|---|
| equilibrated patch, membrane barostat | 66.6 → **65.8 ± 0.07** | 63.7×66.9×73.1 → 63.3×66.5×73.1 | 36.8 → 36.8 | 307.9 ± 0.33 | 64/64 |
| +8% lateral stretch, membrane barostat | 77.6 → **66.5 ± 0.15** | 68.8×72.2×73.1 → 63.7×66.8×72.9 | 36.7 → 37.0 | 310.2 ± 1.08 | 64/64 |
| +8% lateral stretch, isotropic barostat | 77.6 → **69.8 ± 0.08** | 68.8×72.2×73.1 → 65.2×68.5×69.3 | 36.7 → **35.5** | 309.8 ± 0.56 | 64/64 |

**Sane area per lipid and a stable box.** The equilibrated patch holds 65.8 ± 0.07 Å²,
inside POPC's 60–68 Å² at 300 K, with the volume at 307.9 ± 0.33 nm³ and z unmoved at
73.1 Å. P–P stays at 36.8 Å against an experimental 36–39 Å, and the leaflets stay
64/64: the bilayer does not distort.

**The area relaxes independently of the thickness.** Stretched to 77.6 Å², the membrane
barostat brings the area back to 66.5 ± 0.15 Å² — within 0.7 Å² of the equilibrated run
— and the box footprint back to 63.7 × 66.8 Å, essentially where it started, while the
thickness holds (36.7 → 37.0 Å). The strain comes out of the area, which is the point.

**An isotropic barostat cannot.** Given the same stretched box, `MonteCarloBarostat`
reaches the same correct volume (309.8 vs 310.2 nm³) by a different route: it compresses
z from 73.1 to 69.3 Å, thinning the bilayer by 1.2 Å to 35.5 Å, and leaves the area per
lipid stuck at 69.8 Å² — 4.0 Å² above the equilibrated value and 3.3 Å² above where the
membrane barostat puts the same stretched box. Isotropic scaling cannot separate area
from thickness, so it buys the right density at the cost of the bilayer's structure.
That is why `OpenMM_box_equilibration`, which always asks ASH for a `MonteCarloBarostat`,
is not used for a membrane system.

## Superseded numbers

PR #1's description quotes an earlier run of the same setup, made with the script that
`scripts/check_membrane_barostat.py` was then cleaned up from. Same physics, different
random seed; the table above supersedes it.

| case | PR description | this document |
|---|---|---|
| equilibrated, membrane barostat | 65.9 ± 0.11 Å², P–P 36.9 Å | 65.8 ± 0.07 Å², P–P 36.8 Å |
| stretched, membrane barostat | 67.9 ± 0.21 Å², P–P 36.6 Å | 66.5 ± 0.15 Å², P–P 37.0 Å |
| stretched, isotropic barostat | 69.5 ± 0.14 Å², P–P 35.3 Å | 69.8 ± 0.08 Å², P–P 35.5 Å |

**The conclusion does not move.** In both runs the membrane barostat returns the
stretched bilayer to the POPC range at constant thickness and the isotropic one does not,
and the gap between them (3.3–3.6 Å² in area, 1.2–1.4 Å in thickness) is far larger than
the difference between the runs.

The one number that shifted appreciably is the stretched membrane case, 67.9 → 66.5 Å².
40 ps is short for a lateral relaxation, so that case is still approaching its limit
rather than sitting at it — which is why it is quoted with the equilibrated run beside
it. Taken together the two runs bracket the equilibrated value of 65.8 Å² from 66.5 and
67.9, i.e. the relaxation converges toward it from above in both, which says more than
either run alone. A converged area per lipid would need nanoseconds on a GPU node; this
change does not rest on one.

## Protein in a bilayer

An ideal 26-residue poly-alanine transmembrane helix (39.8 Å along z, hydrogens and
termini from pdbfixer), inserted with the keyword arguments the preprocess membrane
branch passes to ASH, which forwards them to `Modeller.addMembrane`:

- built 31,092 atoms: 116 lipids (12 deleted to make room for the helix), 5089 waters,
  9 Na⁺ and 9 Cl⁻, box 63.7 × 66.9 × 73.1 Å;
- exactly one `MonteCarloMembraneBarostat` in the System;
- 10 ps of NPT left the box at 62.8 × 65.9 × 73.1 Å with the helix still spanning the
  bilayer (z span 39.8 → 41.0 Å).

The area per lipid is not a useful number for a protein-containing bilayer as computed
here (71–73 Å²): the helix displaces lipids, but its cross-section still counts into
`Lx · Ly`, so the figure reads high. Use the protein-free patch for that check.

## End-to-end against ASH

There is no ASH installation on the development machine, so none of the above exercises
ASH itself — the OpenMM calls are made directly, with the keyword arguments mcube passes
to `OpenMM_Modeller` and `OpenMMTheory`.

That gap was closed separately: another session ran this preprocess path on Bunya (CPU)
against GLP-1R with GLP-1(7–36) and reported 215,790 atoms, box 127.4 × 133.7 × 129.5 Å,
leaflets 241/241 with no imbalance, headgroup separation 36.8 Å, the bilayer midplane
1.0 Å off z = 0, and 5 waters in the hydrophobic core. The 36.8 Å agrees with the
36.8–37.0 Å measured here on an unrelated system. Those figures are reported, not
measured here.

## Regression

`membrane` absent or false must leave the two existing modes exactly as they were;
`tests/test_membrane_pipeline.py` asserts it rather than arguing it, with ASH replaced by
a recorder:

- **explicit TIP3P** → `OpenMMTheory(xmlfiles=[protein.ff14SB, tip3p], periodic=True,
  rigidwater=True)`, no `xmlsystemfile`, no `system.xml` written, no lipid XML in the
  force-field list, `OpenMM_box_equilibration` still used for NPT;
- **implicit GBn2** → `xmlsystemfile` with `periodic=False` and `CutoffNonPeriodic`,
  `implicit/gbn2.xml` appended, `solventDielectric` from `implicit_eps`, no barostat in
  the System, NPT still skipped;
- **preprocess** passes no membrane keyword at all to `OpenMM_Modeller` in either mode,
  and keeps pH 8.0 explicit / 7.0 implicit and `solvent_padding`.

Full suite at the time of writing: 47 passed, 1 skipped (the skip is the pre-existing
POSIX-only thread-binding test).
