# Ion-channel permeation PMF: decision log

A running record of the modelling decisions for the umbrella-sampling study of Ca²⁺ and two
lanthanide ions passing through the channel in `docs/channel/sequence.fasta`. Each entry gives
the decision, why, the evidence, and whether to revisit it. Same format as
[csr-decision-log.md](csr-decision-log.md). The protocol is in
[channel-us-plan.md](channel-us-plan.md).

**Aim.** Potentials of mean force for a light and a heavy lanthanide and for Ca²⁺ crossing the
channel, built with mcube's existing stages as closely as possible.

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

**C3 – Bilayer: POPE, single lipid type, at 310 K.** (2026-10-22; open)
- mcube already accepts `membrane_lipidtype: POPE` (`membrane.py:LIPID_TYPES`), but only POPC
  has been validated (`docs/membrane-validation.md`). A *mixed* POPE/POPC bilayer is not
  supported: `addMembrane` takes one lipid type.
- POPE's main phase transition is near 25 °C, so 300 K sits on it; run at 310 K.
- The outer membrane's LPS outer leaflet cannot be built by mcube; POPE is the usual stand-in.
- Revisit: whether to add a mixed bilayer (user: "may have to add POPE in addition to POPC"),
  and a POPE area-per-lipid check (literature ~56–60 Å² at 310 K) before production.

**C4 – Ions: Ca²⁺, La³⁺ and a heavy lanthanide; heavy partner not yet chosen.** (2026-10-22; open)
- Gd³⁺ is 4f⁷, multiplicity 8. Dy³⁺ does not fix this: it is 4f⁹ (multiplicity 6) with orbital
  angular momentum, which is harder, not easier. Lu³⁺ (4f¹⁴) is closed shell and is the heaviest
  lanthanide, so La³⁺ vs Lu³⁺ spans the whole lanthanide contraction with two singlets.
- Pending: literature check of ground terms, hydration numbers and GFN2-xTB lanthanide
  treatment; the user's choice of partner.

## Method

**C5 – Hybrid: MM umbrella sampling over the full path, QM/MM umbrella sampling over the
constriction, MM→QM/MM reweighting, then xTB→DFT single-point reweighting.** (user, 2026-10-22)
- Why: QM/MM windows can afford ~10–20 ps at ~1.4–2 s per step (`running-on-bunya.md`), while
  first-shell water exchange for these ions is on the ns scale, so QM/MM alone cannot sample
  dehydration. MM supplies the sampling; QM/MM corrects the chemistry where it matters.
- Primary literature reference: pending (research in progress).

**C6 – One independent simulation per window.** (user, 2026-10-22)
- Each window is its own job and its own config section, seeded from a frame of the steered pull.
  Embarrassingly parallel, restartable per window, and a failed window costs only itself.

**C7 – A new `umbrella` stage in `src/mcube/umbrella.py`, beside `neb`.** (user, 2026-10-22)
- One transparent, minimal config section, like every other stage: inputs written by the user,
  results written back by mcube. Built out from ASH's biased-sampling example.
- Whether the bias is an OpenMM force added to the system or a PLUMED input: pending the ASH
  source check.

**C8 – Metals stay set in `preprocess`.** (user, 2026-10-22)
- The ion is listed under `metals`, so preprocess detects its coordination sphere and defines the
  QM region. Coordination restraints apply only while each window equilibrates; they are released
  for sampling and are never applied in QM/MM. This matches the rule "restrained for the initial
  MD equilibration, not for QM/MM".
- Consequence: the `md` stage cannot be reused unchanged for the MM windows, because it keeps
  the coordination restraints through production (`md.py:run`, restraints passed to
  `OpenMMTheory` for every stage). The umbrella stage holds them for its equilibration leg only.

**C9 – The QM region is rebuilt per window.** (2026-10-22)
- The ion moves ~70 Å along the path, so one fixed QM atom list cannot serve every window. Each
  QM/MM window builds the metal + coordination-sphere region from its own starting frame with
  the existing `build_qm_region` logic.

## Findings

**F1 – `qmmm` does not apply `metal_parameters`.** (2026-10-22)
- `md.py:apply_metal_parameters` swaps in IOD LJ terms, but `qmmm.py:make_mm` builds the MM
  system from the force-field XMLs alone. A QM metal therefore interacts with MM atoms through
  the stock Amber LJ in QM/MM and through IOD in MM. The umbrella stage must use one set in
  both; whether to fix `qmmm` too (it affects the CSR runs) is the user's call.

## Planned checks, not yet done

- Q1: GPU. OpenMM platforms in the ASH env on Bunya (CUDA on `gpu_cuda`, HIP on MI355X); OpenMM
  MM part of QM/MM on GPU vs CPU; whether openmm-plumed is installed.
- Q2: do `amber14/tip3p.xml` templates exist for Ca²⁺, La³⁺, Lu³⁺ in the installed OpenMM?
- Q3: one QM/MM window in the membrane system: timing, energy drift, box stability (QM/MM in a
  membrane has never been run in mcube).
- Q4: pore radius profile of the AF2 model before planning windows.
- Q5: 12-6-4 as a `metal_parameters` option (user asked why not; formula and tables pending).
