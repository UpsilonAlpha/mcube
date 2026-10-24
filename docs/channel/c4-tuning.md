# Pair-specific 12-6-4 C4 tuning: what the Merz group does, and how it maps to the channel

Written 2026-10-25 from a literature pass (an agent read the papers' methods sections; items
marked *unverified* could not be checked). Why it matters here: decision log C17. 12-6-4
over-favours Ln³⁺ against Ca²⁺ by 11–13 kcal/mol where r2SCAN binds Ca to acetate through one
oxygen, so the Ca-vs-Ln answer of an MM PMF through the Glu300 constriction is suspect until the
ion–carboxylate C4 terms are checked.

## The precedent

| Paper | Tuned | Target | Free energy method | Result |
|---|---|---|---|---|
| Panteva, Giambaşu & York, JPCB 2015, 119, 15460 (m12-6-4) | ion–site C4 only (phosphate O, N7), per ion | site-specific log K from titrations | TI with flat-bottom restraint + analytic standard-state term; PMFs only for kinetics | error 4–5 → ≤ 0.1 kcal/mol |
| **Jafari et al., JPCB 2024, 128, 684 (metal–acetate)** | acetate-O polarizability → ion–O C4, **per ion and per water model** (TIP3P, SPC/E, OPC); 11 ions incl. Ca²⁺, no Ln | experimental ΔG (Ca: −1.27 kcal/mol, log K ≈ 0.93) | PMF on M–carboxyl-C distance, 0.05 Å windows < 5 Å, 0.1 Å above, 6 ns/window, WHAM, 3 replicates | within 0.35 kcal/mol; hydration-fitted 12-6-4 over-binds acetate; Ca still bidentate (~52 % in TIP3P) |
| Fındık et al., JCTC 2024, 20, 4298 (Ca/Mg–phosphates) | M²⁺–phosphate-O C4 per ion **and per protonation state**; negative C4 allowed | experimental association constants | PMF on P–M distance, 2 + 8 ns/window, WHAM, 3 replicates | ±0.3 kcal/mol |
| Zhang, Jafari … Merz, Hu, Nat. Commun. 2024 (BbZIP); Jafari et al., JCTC 2025, 21, 4221 | used the tuned C4 sets for Zn²⁺ with Asp/Glu/His | — | **well-tempered volume-based metadynamics** (PLUMED; spherical position + coordination number), ≥ 500 ns per replica; **umbrella sampling failed** (high barriers, helix unfolding) | validated experimentally (cross-linking, accessibility, transport) |

Sources: PMC11137831, PMC12203879 (DOI 10.1021/acs.jpcb.3c06567), PMC11549095, PMC12020355,
PMC4762653, PMC10372143 (ZIP8, m-12-6-4 with TI, Zn over Cd by 3.8 kcal/mol).

Added 2026-10-25 (read from the papers by hand):

- **Kantakevičius, Mathiah, Johannissen & Hay, JCTC 2022, 18, 2367** ("Ch-BE", doi
  10.1021/acs.jctc.1c00898): pair C4 for carboxylate O and chelator amine N, fitted by TI to
  EDTA/NTA stability constants (TIP3P, ff14SB). Carboxylate-O C4 (default → Ch-BE): Ca²⁺ 34.4 →
  29.2, Mg²⁺ 52.4 → 12.3, Y³⁺ 85.1 → 14.9, **La³⁺ 59.9 → 9.5**. Tested on lanmodulin and calmodulin
  EF-hand peptides: LanM order La > Y > Ca reproduced (La − Ca −40 vs −46 kJ/mol exp.), absolute
  binding still too strong for every ion. No Lu³⁺.
- **Jafari et al. TIP3P values (Table 2):** Ca²⁺ acetate-O polarisability 0.600 Å³ → C4 36.27
  (default 34.4); the full table for 11 ions is in `src/mcube/data/c4_pairs.yaml`.
- Our default values are exactly theirs (Ca 34.4, La 59.9), so their corrections apply directly.

## In mcube

`c4_set: kantakevicius2022_chbe_tip3p` or `c4_pairs: {La3+: {O2: 9.5}}` in md (and any later
stage) — see the README rows and `src/mcube/data/c4_pairs.yaml`. Decision log C20 has the plan.

## Transfer to the channel

1. **What to tune:** only the ion–carboxylate-O C4, separately for Ca²⁺, La³⁺ and Lu³⁺, with the
   acetate model that has the Glu side-chain charges. Ion–water terms stay (Li/Song/Merz TIP3P),
   so hydration free energies are untouched. In mcube this is a per-pair override of the C4
   table that `ions.apply_metal_parameters` builds from the polarizability ratio.
2. **Targets:** first stability constants of the ion–acetate complex at I = 0. Ca: −1.27
   kcal/mol (Jafari). La³⁺/Lu³⁺: log β₁ is roughly 1.8–2.0 and nearly flat across the series
   (*unverified*; take NIST SRD 46 or Martell–Smith and correct to I = 0, which is large for a
   3+/1− pair). Express the targets as ΔΔG(Ln − Ca) and ΔΔG(La − Lu), since selectivity depends on those.
3. **PMF protocol (Jafari's, plus what it leaves out):** M–carboxyl-C distance 1.8–12 Å, 0.05 Å
   windows below 5 Å and 0.1 Å above, ≥ 6 ns each (longer for Lu³⁺), 3 replicates, MBAR or WHAM; add
   the Jacobian 2kT ln r; K° = 4π∫₀^rc r² e^(−W/kT) dr / 1661 Å³ with rc at the first barrier;
   handle the 3+/1− Coulomb tail (≈1.3 kcal/mol at 8 Å in TIP3P, *estimate*): large box or
   subtract the screened Coulomb. Record contact/solvent-shared and mono/bidentate populations.
4. **Cross-check against our QM:** after the fit, compare with the r2SCAN water→acetate exchange
   energies (C17; TZVP + counterpoise running) and check that Ca²⁺ moves toward monodentate.
   Matching ΔG with the wrong denticity is the trap.
5. **Pitfalls:** C4 absorbs missing charge transfer, so negative values are allowed; values fitted in
   bulk water may over-correct inside the low-dielectric constriction (check against
   Ca/Ln–carboxylate crystal geometries or lanmodulin sites); hydrolysis and acetate protonation are
   not in the model; and umbrella sampling failed for ZIP, so check overlap and convergence near
   Glu300, and consider a coordination-number CV or metadynamics as a cross-check.
