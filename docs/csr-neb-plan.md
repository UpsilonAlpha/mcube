# CSR-1 meropenem hydrolysis: QM/MM pathway plan

How we go from the restrained MM model to reaction barriers. The decisions behind each step are in [csr-decision-log.md](csr-decision-log.md).

## Questions

CSR-1 is an unusual B3 with an RQK-HRH active site, and the substitutions are at
metal-ligand positions: **Arg90 (BBL R118) occupies a Zn1 ligand position** and
**Lys226 (BBL K263) occupies the Zn2 ligand position**. A cation cannot ligate
zinc, so the canonical second metal cannot bind and the site carries one zinc at
Zn1.

**Central hypothesis (D107).** Arg90 and Lys226 together supply, as two
monovalent cations, the **+2 that the absent Zn2 would have carried**. They do
Zn2's substrate-anchoring job electrostatically - Lys226 salt-bridges the C3
carboxylate (2.61 A at closest over `nacqm`'s 50 ps - but **this opens to 4.26 A over the following 50 ps and no cation holds the substrate at all; see D111**), Arg90 holds the
6-alpha-hydroxyethyl in 99% of frames - while being unable to do its catalytic
jobs: lowering the nucleophile's pKa and stabilising the anionic nitrogen. That
would account for CSR-1 binding broadly and turning over poorly.

**Where the two cations actually sit, measured in `QM_MM_opt_last.pdb`** - they
act at opposite ends of the site, and only one of them is near the metal:

    Arg90 NH1   Zn 4.41  O01 5.13  N1x 7.03  O3x  9.66
    Lys226 NZ   Zn 8.73  O01 7.13  N1x 6.93  O3x  4.80
    Arg229 NH2  Zn 9.34  O01 7.76  N1x 6.61  O3x  4.38

So **Arg90** is the cation positioned to affect the nucleophile, while
**Lys226** points away from the metal at the substrate carboxylate. A Zn2 metal
would sit ~3.5-4 A from Zn1, not 8.7, so "Lys226 supplies the charge Zn2 would
have carried" is geometrically loose; what it does plainly inherit is Zn2's
**carboxylate-binding** role, since the C3 carboxylate binds Zn2 in canonical
MBLs. **And Arg229 NH2 is closer to both the carboxylate and the leaving nitrogen
than Lys226 is, while sitting outside the QM region** - on geometry alone,
including Lys226 and excluding Arg229 is hard to justify if the question is
cationic stabilisation of the anionic nitrogen.

The paths must answer:
1. **Lys226 (K254 / BBL K263):** does its charge contribute to the barriers, and
   to which step? The hypothesis predicts it matters more for **ring opening**,
   where the anionic nitrogen appears, than for the attack step. That is the
   falsifiable part.
2. **Arg90 (R118) and Arg229:** the same question for the arginines. Arg90 sits
   at a Zn1 ligand position, beside the metal; **Arg229 is cationic and contacts
   the substrate in 89% of frames but is not in the QM region.**
3. **Asp92 (D120):** is it the general base? Its proton-accepting geometry says
   yes (OD1 2.50 A from the nucleophile's proton). Note it is **monodentate to
   the single zinc** (OD1 2.19, OD2 3.47), so it bridges nothing - the di-zinc
   "bridging aspartate" label does not transfer to a mono-zinc site.

How each question is answered:
- **Who moves which proton:** along every NEB path, monitor Lys226 NZ, Asp92 OD1/OD2, the Zn water and O01 protons.
- **Per-residue contributions** (as Imhof did for Lys92 in EcoRV): recompute each path's energies with the charges of Lys226, Arg90, Arg229 and Asp92 zeroed in turn, then optimise the path again with the residue removed or neutralised for the important cases.
- **Where they sit at the start.** The two models disagree about Lys226, which is
  why D107 measured it over a trajectory rather than in one frame:
  - Arg90 and Arg229 hold the 6α-hydroxyethyl and C2 side chain (D29, 99%/89%)
    — **in classical MD. Under QM/MM neither keeps a 3.5 Å contact with the
    substrate in any frame** (D111), so this bullet describes the force field;
  - Asp92 binds the Zn and H-bonds both water oxygens (D29, 99.5%/100%);
  - Lys226 is **2.7 Å from O01 in `CSR_crystal_v2`** (D29) but **5.6-7.3 Å from it
    in `CSR_nac`**, where it salt-bridges the C3 carboxylate at 2.61 Å over the
    first 50 ps, loses that contact over the next 50 ps (D111), and
    never comes within 4.7 Å of the leaving nitrogen (D107). Any claim about its
    role has to say which model it refers to.

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


## How long should the QM/MM equilibration be? (checked 2026-10-25)

Asked because the length has to be defended in the write-up. Everything here was
read in a Methods section unless marked otherwise.

| Study | QM level | classical MD | QM/MM MD | DOI |
| --- | --- | --- | --- | --- |
| Chudyk *et al.* 2014, class A + **meropenem** | SCC-DFTB/ff12SB | min/heat | **50 ps equilibration, then 300 ps production** | 10.1039/C4CC06495J |
| Chudyk *et al.* 2022, same group | SCC-DFTB/ff12SB | - | **300 ps before reaction calcs; 3 x 1 ns per system** | 10.1021/acsinfecdis.2c00152 |
| Dal Peraro *et al.* 2007 | Car-Parrinello BLYP | 4 ns | **20 ps equilibration, ~70 ps per path** | 10.1021/ja0657556 |
| Das & Nair 2017, NDM-1 | CP, PBE | 2 ns NPT + 10 ns NVT | **~10 ps before metadynamics** | 10.1039/C6CP08769H |
| Simona *et al.* 2009, CphA | Car-Parrinello BLYP | - | **3 ps**, ~3 ps per constrained step | 10.1074/jbc.M109.049502 |
| Saito & Takano 2022, urease - **the only enzyme GFN2-xTB/MM MD paper found** | **GFN2-xTB** | 10 ns NPT | **2.0 ps** | 10.1021/acs.jpcb.1c10200 |
| Ackerman & Gatti 2013 | B3LYP scans | 100 ps | **none** | 10.1371/journal.pone.0055136 |
| Lisa *et al.* 2017 | SIESTA/PBE opt | 10 ns NVT | **none** | 10.1038/s41467-017-00601-9 |
| Gervasoni *et al.*, IMP-1 | DFTB3 then B3LYP opt | 20 ps heat + 240 ps NPT + 100 ns NVT | **none**, "starting from the last frame" | - |

**Where our 50 ps sits.** Generous against DFT-level QM/MM MD, which runs 3-20 ps
in this literature. Twenty-five times longer than the one published GFN2-xTB/MM
enzyme protocol (n = 1). But **at the floor against the closest methodological
neighbours**: 50 ps is exactly Chudyk 2014's *equilibration* figure for a
meropenem complex, and there it is followed by 300 ps of production before
structures are chosen - by 2022 the same group runs 300 ps before reaction
calculations and 1 ns per system. The people who use a cheap Hamiltonian spend
the savings on sampling, so that is where a reviewer would push.

**Planned, not yet run (user, 2026-10-25): about 300 ps of production** after the
equilibration, which is what Chudyk 2014 does and what brings us level with the
closest comparable protocol. At the rate measured in `running-on-bunya.md`
(~1.08 s/step at 1 fs with 4 threads) that is roughly 90 h of wall clock per run,
so it wants either a longer walltime request on `general` or splitting across
jobs. **Splitting is now the easier of the two:** `resume: true` on the qmmm
stage continues from ASH's state file and runs only the remainder of `time`, so
a 300 ps target can be a chain of jobs that resubmit the same section unchanged
(D109).

**There is no convention for GFN2-xTB/MM, and no guidance on sufficiency.** No
methodological paper, benchmark or review gives a target length or a convergence
criterion for QM/MM MD equilibration. Convergence statements in these papers
attach to the *classical* stage only. Dal Peraro is the only one that reasons
about QM/MM MD length at all, and negatively: it argues a conformational flip is
"very unlikely to occur spontaneously" on a ~70 ps timescale and so models both
pathways separately rather than trusting sampling. Chudyk tests convergence of
the *free energy* by repeating umbrella sampling three times, not by lengthening
equilibration. **So the honest write-up line is that the field states its length
and does not defend it.**

**Running no QM/MM MD is about half of the MBL set** (Ackerman & Gatti, Lisa,
Gervasoni), and it correlates with DFT/static-path methodology. Also, across the
whole set **no paper reports how many structures entered its path calculations** -
Gervasoni says "starting from the last frame", one structure. Reporting a snapshot
count would put us ahead of the field on that point.

**Generic length is the wrong yardstick for our blocking question.** Zn(II) ligand
exchange is slow next to every timescale above, so none of these numbers - not
even Chudyk's 1 ns - is evidence that any equilibration length escapes an imposed
coordination geometry.

**But equilibration length is not what we judge the zinc against: the reference is
the crystal structure** (user, 2026-10-25). CSR-1 is deposited - 6DN4 (mature wild
type, 1.99 A), 6DQH/6NC5 (2.07 A) - and D39/D40 identify ZN302 as the Zn1 H-R-H
site we keep and ZN301 as the canonical Zn2 site that only the R118H/Q121H/K254H
mutant restores. So the question "is our zinc sphere right?" is answered by
comparison with that density, not by how long anyone's MD ran. This also means the
inaccessible Methods of Yang *et al.* 2021 (10.1021/acs.jcim.1c01109) do not
undermine the zinc argument: that paper was only ever a secondary prop.

**One thing to read off the PDB rather than this log.** D40 lists ZN302's ligands
as "His87/His163 + two waters", but the mapping recorded in D22 is model =
deposited - 28 (model His88 = PDB His116, model His164 = PDB His192), under which
87 and 163 would be model 59 and 135 - not zinc ligands. One of the two numbering
statements is wrong, and since the crystal sphere is now the reference standard,
**the ligand list and its distances should be measured directly in the deposited
file before being quoted.** Our optimised model currently has a five-coordinate
zinc (O01 1.88, His88 2.13, His164 2.14, HOH 2.17, Asp92 OD1 2.19), and whether
that matches the crystallographic sphere is exactly the comparison to make.

**That comparison was made in D120 and the answer is no.** Measured from 6NC5
(ZN302), the crystal sphere is four — HOH485 1.90, His163 NE2 2.09, HOH493 2.14,
His87 NE2 2.22 — with **no carboxylate**. Asp92 OD1 at 2.19 Å above is the fifth
ligand the crystal does not have, and a restraint this pipeline generated put it
there.

**Could not be read:** Imhof 2009 Methods (ACS paywall - so "minimised structures"
is still unverified), Xu-Xie-Guo *JBC* 2006 (not OA; its QM/MM MD length is
unknown, and the 500 ps figure that search surfaces belongs to the separate 2005
binding paper and is snippet-level, so it must not be quoted as a length), VIM-1
(10.1021/acscatal.1c04786, ACS 403), Yang 2021 Methods, Wu-Xu-Guo 2010 SI.
**Feng *et al.* 2017 (10.1038/s41467-017-02339-w) has no computational component
at all** - crystallography and NMR only - so it can be cited for structure, as
D107 does, but never for a protocol.

## Replication protocol for CSR-1 with modern methods

| | Imhof 2009 / Xu–Guo 2006 | CSR-1 (mcube) |
| --- | --- | --- |
| QM level for sampling and paths | AM1/d (Imhof); SCC-DFTB (Guo) | GFN2-xTB/MM (electrostatic embedding) |
| Energies | B3LYP model systems; DFT corrections | ωB97X-D3(BJ)/def2-SVP/MM single points along paths; def2-TZVPD on stationary points (D31) |
| Reactant sampling | minimised structures (Imhof - **Methods not accessible, unverified**); QM/MM MD (Guo - **length unconfirmed**) | 50 ps QM/MM MD (D28/D29), several snapshots |
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

## Status (2026-10-29)

> **Superseding note, 2026-10-29.** Everything in the 2026-10-26 status below
> was written before D120, and one of its bullets is now known false: the zinc
> sphere was listed as **five-coordinate and crystal-referenced, "established
> and not in doubt"**. The crystal sphere is **four** — HOH485 1.90, His163 NE2
> 2.09, HOH493 2.14, His87 NE2 2.22 (6NC5 ZN302) — **and Asp92 is not on it.**
> The fifth ligand was put there by a coordination restraint this pipeline
> generated from a 3.0 Å detection cutoff, which held Asp92 OD2 at its 2.80 Å
> second-shell distance until the carboxylate rotated OD1 onto the metal (D120,
> D121).
>
> **What that costs the rest of the bullet list.** Every trajectory measurement
> below was made on that five-coordinate zinc, so the ones that involve Asp92 or
> the metal are not independent of the artefact. In particular, *"the hydroxide
> stays on the zinc and Asp92 OD1 holds it, 100% of frames"* was measured with
> Asp92 restrained against the metal — it has to be re-established on a shell
> that matches the crystal before it can be cited. The cation positions (Arg90,
> Lys226, Arg229) are not affected.
>
> **Two further blockers found since.** The pose everything ran on from October
> was the wrong one, and 2.55 Å was its restraint's own output rather than a
> measurement (D123, D126); and the start structures built in D126/D128 graft
> onto `md_naclong`'s frame, which carries the condemned shell forward, so
> `md_subst` (job 30923516) was cancelled an hour in (D130).
>
> **The current blocker is therefore the shell, not the pose**, and D120 already
> costed the route through it: preprocess from the crystal with `cutoff: 2.5`,
> no Asp92 restraint, then re-run everything from `md_nac` onward. Read
> `## Standing constraints` at the top of the decision log before setting any of
> it up. M0 below is still the right shape; its starting structure is not.

**The blocker is the substrate pose, and it sits upstream of everything else in
the M-series.** Two poses are now excluded. The pinned one - every structure
before D102 descended from a meropenem whose carbonyl was held on the zinc at
1.94 A through all of equilibration - is an artefact of its own restraint (D56,
D59, D60). The free one has now been sampled, across 100 ps of continuous
GFN2-xTB QM/MM MD in two consecutive segments (D104, D110), and it does not
form: no cation holds the substrate in a single frame, the Lys226-carboxylate
contact opens from 2.61 to 4.26 A, and the O01-C4x attack distance wanders over
3.63-7.23 A (D111). **There is no Michaelis complex to compute a path from**,
which is what M0 below is for.

**What is established and not in doubt:**

- **The nucleophile and its base.** Over the same 100 ps the hydroxide stays on
  the zinc (1.73-1.97 A, 100% of frames) and Asp92 OD1 holds it (2.63-3.20 A,
  mean 2.94, 100%) - D111. Asp92's proton-accepting geometry is the one part of
  the mechanism every calculation here agrees on.
- ~~**The zinc sphere**, five-coordinate under GFN2 and referenced to the
  crystal rather than to any simulation length (D108; user, 2026-10-25): O01
  1.90, His164 NE2 2.10, Asp92 OD1 2.11, HOH 2.16, His88 NE2 2.18. Asp92 is
  **monodentate** (OD2 4.22), so it bridges nothing.~~
  **[FALSE — D120.]** The crystal sphere is **four** and carries no carboxylate;
  "referenced to the crystal" was never checked against the deposited file.
  Asp92 OD1 at 2.11 Å is a restraint artefact, not a ligand. Nothing that
  depends on this bullet survives it.
- **Both cations are where the hypothesis needs them**, and cationic throughout:
  Arg90 at a Zn1 ligand position, Lys226 at the Zn2 position (D107, D110). The
  +2-substitution hypothesis stands as a hypothesis; what it does not have is
  evidence from any trajectory (D111).
- **100 ns of classical MD** reached a near-attack conformation (3.33 A, 101
  degrees; D102) that 50 ps of QM/MM destroyed. The D29 contact fractions are a
  force-field result QM/MM does not reproduce (D111), so they cannot be cited as
  structure.

**Infrastructure, so it is not re-derived:** QM/MM MD costs 1.08 s/step at 1 fs
on 4 threads, 1.72 s/step with Arg229 in the QM region; xTB stays on CPU; and
`resume: true` lets a long run be a chain of jobs (D109, `running-on-bunya.md`).

## Principle

**The metal site is never left free in MM.** Every classical stage keeps the Zn coordination geometry held by constraints and restraints. That is the purpose of the mcube pipeline. The Zn site is only free in QM/MM, where the Zn, its ligands and the substrate are all in the QM region and no constraints apply.

So the sequence is:

**preprocess → restrained MM equilibration/NVT → QM/MM** (optimisation, scans, NEB, TS/IRC)

There is no unrestrained MM MD at any point.

## System

**The live system is `/scratch/user/uquadhik/CSR_nac`** (2026-10-26). The table
below describes `CSR_crystal`, which every structure up to D102 came from and
whose atom indices the older entries in the log quote; **its substrate pose is
the pinned one M0 rejects**, so it is kept for reading the history, not for
building anything new. Model numbering = crystal numbering − 28. Atom indices are
0-based in `preprocess_results/finalsystem.pdb` and are **not** transferable
between the two systems.

| Item | Value (`CSR_crystal`) |
| --- | --- |
| Zn ligands | His88, His164 (NE2), OH⁻ (O01), Asp92 OD2, meropenem β-lactam O2x (pulled in) |
| Key atoms | Zn 4203, O2x 4209, C4x 4208, OH⁻ O01 4255 (H03 4256), Zn-bound water O 4257 |
| QM region | 119 atoms: Zn, OH⁻, water, meropenem, His88, His164, Asp92, Arg90, Lys226; charge +2, multiplicity 1 |
| QM method | GFN2-xTB with electrostatic embedding (exploration). DFT/MM for final energies. |
| Active region | 8 Å around the Zn (atom 4203). One `active_atoms` list per snapshot, reused by every later step for that snapshot. |

**`CSR_nac`, as actually run** (D102–D111): 49,424 atoms; zinc five-coordinate on
His88 NE2, His164 NE2, Asp92 OD1, OH⁻ (O01) and one water — **the wrong sphere,
and the reason this run has to be redone: the crystal has four and no
carboxylate (D120)** — with the substrate carbonyl **not** coordinated (Zn–O2x
never closer than 4.77 Å, which a saturated metal cannot avoid); QM region 119
atoms at charge +2, or **137 atoms at charge +3** with Arg229's side chain added
(`qmmm_nacarg`, D110); active region 5 Å around the Zn, unioned with the QM atoms.
Since D105 `resolve_actatoms` unions the QM atoms in unconditionally and reports
which it had to add, so a small `actradius` can no longer freeze a QM residue.

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
# qmmm_snapN.inputs.pdb_file: nac_frameXXXXX.pdb   (paths read against the config)
mcube qmmm config.yaml --section qmmm_snapN
```

- Settings:
  - `skip_equilibrate: true` (no QM/MM MD); the DFT score is a separate `sp` run;
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
  - the same `actradius` and `actregion_origin` as the endpoints, each run
    deriving the region on its own structure (compare the recorded
    `actatoms_count`; a difference means the radius cuts through something that
    moved);
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
2. **Relaxed scan:** `mcube scan config.yaml --section scan_attack`, with `pdb_file`, `atoms`, `start`, `end` and `points` in the section. Writes `scan_XX.pdb` and `scan_energies.csv`.
3. **NEB settings from the config:** `neb.inputs` (or a per-step section given with `--section`): `images` (12), `CI`, `maxiter`, tolerances, interpolation, springs, `method: nebts`. Writes `saddlepoint.pdb`.
4. **TS refinement and IRC:** `mcube ts config.yaml --section ts_attack`, with `pdb_file` (the NEB's `saddlepoint.pdb`) in the section. Runs the TS optimisation with a partial Hessian over `hessian_atoms`, a frequency check, and a quasi-IRC along the imaginary mode (or Sella), with `monitor` distances reported.

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

### M0. Establish a Michaelis complex - the blocker

Added 2026-10-26. The M-series had no such step because when it was written the
pose was assumed; D111 removed that assumption and both candidate poses with it
(see Status). Nothing downstream is computable without a bound complex. D111's
three readings are the candidates, and they want different experiments:

| Reading | Test | Cost |
| --- | --- | --- |
| The reactant protonation state is wrong, so what we are holding together is not what binds - **and the live model is in the state D48 rejected, see M1** | the `protH_lys`/`protH_oh` comparison: xTB/MM geometries, DFT/MM single points | cheap; the sections are already written |
| GFN2-xTB cannot hold this complex | DFT/MM optimisation of the complex - the retest this log has specified since D52 | moderate, one structure |
| The complex is real but needs a restraint to stay in the well through equilibration | restrained QM/MM equilibration with the carboxylate on the cations | cheap, but makes the pose an assumption again |

**The third is not a fix and should not be reached for first.** Restraining the
pose we want is how the 1.94 A artefact happened; doing it again with a different
restraint reproduces the same methodological error with the opposite bias. It is
defensible only as a *starting* device whose restraint is then released and the
complex shown to hold - and that is the test which has never passed.

**A fourth possibility has to be stated, because the data are consistent with
it: CSR-1 may bind meropenem weakly.** That would fit broad-spectrum binding
with poor turnover, and it is the +2 hypothesis's own prediction - a site that
anchors electrostatically but cannot do Zn2's catalytic work. But a drifting
simulation cannot establish it (D111), and **nothing in the protocol as written
could tell a weak complex from a bad model**, because every step from M2 onward
presumes a well to compute in. Separating them needs an affinity-like quantity
rather than a path, which is the one circumstance that would force umbrella
sampling back in from "Deliberately not doing".

### M1. Settle the reactant - open, with three incompatible answers

The proton on Lys226 versus on O01. Job 30053828 (D42) ran; its answer has not
survived:

| Source | Lys226 -> O01 transfer | Status |
| --- | --- | --- |
| D29/D30 | -13 to -16 kcal/mol | computed in the six-coordinate site, so it does not settle the question for the correct one (D38) |
| D48, clean geometry (job 30053828) | **-5.23 kcal/mol: resting state is Zn-OH2 with neutral Lys226** | about a third of D29/D30; the rest was geometric artefact. D48 calls it "not comfortably excluded, only disfavoured" |
| D107, pKa bookkeeping | **+6 to +7, i.e. uphill** | opposite in sign |

GFN2-xTB's worst benchmarked functional group is exactly the protonated amine,
MUE 5.3-7.7 kcal/mol (Arantes & Rezac 2025, doi:10.1021/acs.jctc.5c00690),
which **spans the entire disagreement** - so xTB cannot arbitrate this and the
DFT/MM single points are not optional. The `protH_lys`/`protH_oh` sections exist
but have never run in their fixed form, and they now have a clean geometry to
run in. Being also the cheapest of M0's three tests, this is the obvious first
move: it settles M0's first reading and M1 together.

**The state M1 settled on was never applied to the system we are simulating.**
D48's verdict was that the resting state is Zn-OH2 with neutral Lys226 - "the
model genuinely has no hydroxide". But `CSR_nac`, the system behind D102-D111,
carries **O01 as a hydroxide (one hydrogen, H03) and Lys226 cationic (HZ1, HZ2,
HZ3)** - the `lysH` state, the one D48 disfavoured by 5.2 kcal/mol. D98 built it
from `CSR_unrestrained2` with "the same `preprocess.inputs`", and nothing in the
log revisits the protonation at that point. So **every QM/MM result from D102
onward was computed in the state M1 rejected**, which promotes M0's first
reading from one possibility of three to the leading one (D112).

Two things stop this being a simple error to correct. D48 itself called
-5.2 kcal/mol "not comfortably excluded, only disfavoured". And D107's pKa
bookkeeping points the other way by +6 to +7, which would *favour* `lysH`. **The
two best arguments disagree in sign**, which is exactly what a DFT/MM comparison
in a clean geometry is for.

### M2. One mechanism, one structure - blocked on M0

A **water attacking C4x with Asp92 (= Asp120) as the general base**, computed
first, on the single best reactant:

scan (2D: O-C4x x H...Asp92 OD1) -> NEB-CI -> TS optimisation -> IRC -> DFT/MM single
points along the path. Everything xTB/MM except the last step. **The scan has to
drive the attack coordinate:** two independent 50 ps QM/MM runs agree the MD will
not hand us a near-attack frame (D104, D110).

**No published mechanism transfers whole, and this step should stop implying one
does.** It used to read "exactly the CphA case". It is not: CphA is a B2 enzyme
with a four-coordinate Cys-containing site and no Zn-carbonyl contact at all, so
what carries over is the principle - water nucleophile, carboxylate base - and
not a geometry or a number (D60). The closest structural analogue is the only
published mono-zinc B3 mechanism, GOB-18 (Lisa 2017), and **its metal sits at Zn2
while ours sits at Zn1** (D107), so the species that stabilises the anionic
nitrogen there is absent here. Asp92 as the base rests on our own geometry, held
in 100% of frames across 100 ps (D111), not on a borrowed scheme.

Only if it fails do the alternatives get the same treatment, in this order: Zn-bound
water with Asp92 as base, then concerted N1x protonation.

### M3. Answer the three questions by charge deletion - blocked on M0

The cheap part, and what the questions actually require. On the one converged path,
recompute the energies with selected residues' charges zeroed:

| Question | Residue | Region |
| --- | --- | --- |
| 1. Does Lys226's charge contribute to the barriers, and to which step? | Lys226 | QM |
| 2. Do the arginines? | Arg90, Arg229 | **both QM** |
| 3. Is Asp92 the general base? | Asp92 | QM (and M2 tests it as the base directly) |

Two entries changed on 2026-10-26. Arg229 was MM here until D110 put its 18
side-chain atoms in the QM region: leaving the nearest cation to the substrate
outside the quantum region while keeping a further one inside it could not be
defended. And question 3 no longer calls Asp92 "bridging" - it is monodentate to
the single zinc (D108).

**The hypothesis predicts an asymmetry, and that is the falsifiable part.** If
Arg90 and Lys226 supply the +2 the absent Zn2 would have carried, their charges
should matter more for **ring opening**, where the anionic nitrogen appears, than
for the attack step. Imhof's EcoRV conclusion - Mg(II) activates the water, Asp90
takes the proton, Lys92 positions the substrate and stabilises the intermediate -
is the **methodological** comparison; the structural analogy to Lys92 is retired
(D107).

**Why this waits on M0.** Charge deletion on a substrate sitting 5 A from every
cation would measure the drift, not the catalysis (D111).

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

**M0 may force one of them back (2026-10-26).** If the complex turns out to be
genuinely weakly bound rather than badly modelled, no path calculation can show
it, and separating the two needs a binding free energy - which puts umbrella
sampling back on the table for a reason the D47 reasoning did not anticipate:
not to compare a barrier with k_cat, but to establish that there is a well.

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
