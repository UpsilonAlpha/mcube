# Reference structures, and how each number in `michaelis_hypotheses` was got

Every range in `site.yaml`'s `michaelis_hypotheses` comes from a deposited
structure measured here, not from a review's prose. This file records the PDB
ID, the ligand code, the exact command and the output, so a reviewer — or we,
in six months — can re-derive the number rather than take it.

Download a structure into an empty directory and measure it:

```bash
curl -o 1tmn.pdb https://files.rcsb.org/download/1TMN.pdb
python refstruct.py shell 1tmn.pdb
python refstruct.py centre 1tmn.pdb 0ZN 317
```

`refstruct.py` reads only blank and `A` alternate conformations, and reports
donors within 3.2 Å of every zinc — wider than our own 2.6 Å cutoff on purpose,
so a borderline contact is visible instead of hidden by the cutoff.

## CSR-1 itself — the sphere every other number is referenced against

This belongs first and was missing from this file until 2026-10-29, which is
part of why D120 took three weeks to find: the number that decides whether a
model's metal site is right was recorded only in a decision entry.

```bash
curl -o 6nc5.pdb https://files.rcsb.org/download/6NC5.pdb
python refstruct.py shell 6nc5.pdb
```

6NC5 holds **six** zincs. Only one of them is the site we model.

| zinc | donors within 3.2 Å | what it is |
| --- | --- | --- |
| **ZN302** | HOH485 **1.90**, His163 NE2 **2.09**, HOH493 **2.14**, His87 NE2 **2.22** | **the site we model.** Wild-type ligands |
| ZN301 | Asp91 OD2 1.95, His92 NE2 2.01, His225 NE2 2.04, HOH493 2.08 | the Zn2 site — **present only because of the mutations**, see below. **Asp91 here is our Asp92**: the numbering is ours minus one, so the aspartate whose coordination we keep arguing about is a **Zn2 ligand** (D133) |
| ZN303 | HOH485 1.88, HOH516 1.96, HOH504 2.13, acetate 2.79 | adventitious (3.39 Å from ZN302) |
| ZN304–306 | Glu115/Glu141/His111/Asp129/Glu123, phosphate | adventitious surface sites |

**ZN302 has four donors: two histidines and two waters. No carboxylate.** Asp92
(BBL D120) sits at **2.80 Å** on OD2 — a second-shell contact, not a ligand.

**Why ZN301 does not contradict "CSR-1 is mononuclear", and in fact proves it.**
6NC5 is an engineered construct, titled the *HARLDQ motif*, and its `SEQADV`
records say what was engineered:

```
SEQADV 6NC5 HIS A   92  UNP  A0A384EQ6 GLN   136 ENGINEERED MUTATION
SEQADV 6NC5 HIS A  225  UNP  A0A384EQ6 LYS   269 ENGINEERED MUTATION
```

Those are **our Gln93 and Lys226** (BBL Q121 and K263) — the Q and the K of the
RQK motif — mutated to histidine. That restores a canonical Asp-His-His Zn2
ligand set, and a zinc duly binds it 3.43 Å from ZN302, sharing a bridging water
(HOH493). In the wild-type enzyme those positions are a glutamine and a lysine,
**which cannot ligate a zinc at all.** So the only way anyone has put a metal in
CSR-1's Zn2 site was to engineer two of its ligands back, and ZN302's own
ligands are untouched by either mutation.

**Counting coordination for this enzyme.** CSR-1's Zn1 site is **H-R-H, with two
coordinating residues, not the canonical three**: Arg90 (BBL R118) stands where
the third histidine would be. So the canonical-equivalent four-coordinate sphere
here is **two histidines + hydroxide + one water**, and an extra solvent ligand
must be counted against two protein ligands rather than three before it is
called over-coordination.

Three consequences, all of them load-bearing:

1. **Any model whose zinc carries five ligands, or a carboxylate, is wrong.** A
   saturated metal has no site for the substrate carbonyl, which is a sufficient
   explanation of substrate loss without blaming the QM method (D120).
2. **2.80 Å is where Asp92 belongs, and it is not a bond.** Restraining it there
   is what rotated OD1 onto the metal. See `## Standing constraints` in the
   decision log and `LIGAND_MAX` in `preprocess.py`.
3. **CSR-1 is a Zn1-site enzyme, and every mononuclear MBL in the literature is
   Zn2-site** (GOB-18, Sfh-I, CphA), so none of those analogies transfers
   cleanly — which is why the thermolysin and carboxypeptidase A rows below,
   mononuclear hydrolases with one zinc doing both jobs, carry the criteria
   instead (D124).

`6DN4` is the apo/second CSR-1 entry used as a cross-check in D120.

## The two mononuclear zinc hydrolases with transition-state analogues bound

These are the structures that decide the criteria, because they are the only
well-characterised cases of the arrangement we are proposing for CSR-1: a
zinc-bound water attacking a carbonyl that is itself on the same zinc, with a
carboxylate as general base.

| structure | analogue | Zn–O (nucleophile side) | Zn–O (carbonyl side) | O–X–O | O···O |
| --- | --- | --- | --- | --- | --- |
| thermolysin **1TMN** | `0ZN` 317, N-[(1R)-1-carboxy-3-phenylpropyl]-L-leucyl-L-tryptophan | OG1 **2.01** | OG2 **2.41** | 119° | 2.15 |
| carboxypeptidase A **7CPA** | `FVF` 309, phosphonate transition-state analogue | OPV **2.19** | OPW **2.98** | 117° | 2.68 |

`1TMN` zinc sphere, five-coordinate: His142 NE2 1.95, Glu166 OE1 1.95, **0ZN
OG1 2.01**, His146 NE2 2.04, **0ZN OG2 2.41** (Glu166 OE2 at 2.87 is the second
carboxylate oxygen, outside the sphere).

`7CPA` zinc sphere: His69 2.08, His196 2.09, Glu72 bidentate 2.24/2.29, **FVF
OPV 2.19** and **OPW 2.98**.

`refstruct.py centre` finds the tetrahedral centre by connectivity, so it
reports every P or C in the named residue carrying two or more oxygens. In both
files that means extra hits — 1TMN's `C1` carboxylate 7–8 Å from the metal,
7CPA's `C1` and `C4` — which are not the catalytic centre. The one to read is
the centre whose oxygens are the ones near the zinc.

**What this fixes, and what it does not.** It establishes that one metal can
carry both oxygens, which is what D123 wrongly denied, and it sets the shape of
the criterion: one oxygen tight (2.01, 2.19), the other **loose** (2.41, 2.98).
It does **not** validate an approach angle. Both are post-attack geometries —
the two oxygens are bonded to the same tetrahedral centre, so 117–119° is a
consequence of that bonding. Note also that 7CPA's O···O of 2.68 Å is the same
as D34's 2.67–2.77 Å, which failed: **oxygen separation does not discriminate**.
The approach angle has to be measured on our own trajectory.

## The mononuclear metallo-β-lactamases, and why none of them transfers

| structure | enzyme | site | bearing on CSR-1 |
| --- | --- | --- | --- |
| **5K0W** | GOB-18 (B3, mono-Zn) | Zn2 (DHH); Gln116 OE1 on the metal at 2.04/2.21 Å | closest by subclass, wrong site |
| **1X8I**, **1X8G** | CphA (B2, mono-Zn) | Zn2 (Asp-Cys-His) | the route ruled out twice |
| **7UHI**, **7UHJ** | L1 (B3, di-Zn), intact moxalactam | Zn1 + Zn2 | the only measured Michaelis complex: carbonyl O1 at 2.53/2.63 Å from Zn1 |
| **7UHH**, **7UHK** | L1, the 20 ms and 80 ms snapshots of the same series | Zn1 + Zn2 | the intact β-lactam is seen to 100 ms; cleavage at 150 ms |
| **5B15** | SMB-1 (B3, di-Zn) + **hydrolysed doripenem** | Zn1 + Zn2 | the structure Mu & Xu built their ampicillin Michaelis model from, by swapping the ligand on its core (D155). A carbapenem, so the core match to meropenem is closer than the swap they published |

**Every mononuclear MBL in the literature carries its zinc in the Zn2 site.**
CSR-1's is in Zn1 (His87/His163 = BBL H116/H196). It appears to be the only
mono-zinc MBL whose metal sits where the bridging hydroxide normally lives,
which is why the CphA analogy failed twice: it is not a mechanism for a
Zn1-site enzyme. GOB-18 goes the other way too — a non-metal-bound water
activated by a hydrogen-bond network, with the metal taking the substrate's C-3
**carboxylate** rather than its carbonyl, and the sphere expanding four- to
five-coordinate on binding. That is a Zn2-site description throughout, so it
describes what CSR-1's *vacant* site would have done.

On glutamine: SMB-1's Gln157 is a QM/MM proposal, not a crystallographic metal
ligand. The real structural precedent for a B3 glutamine on a zinc is GOB-18's
Gln116 in 5K0W (2.04 and 2.21 Å in the two chains, giving an HQH site in place
of the canonical 3-His) — but that is BBL position **116**, a Zn1 ligand, while
CSR-1's Gln93 is BBL **121**, a Zn2 position. Same residue type, different site.

## Papers

- Wu, Xu, Guo et al., **QM/MM of carboxypeptidase A**, on the general acid–base
  role of Glu270 — <https://pmc.ncbi.nlm.nih.gov/articles/PMC2741171/>
- **Designing hydrolytic zinc metalloenzymes**, for the mononuclear-site
  coordination conventions — <https://pmc.ncbi.nlm.nih.gov/articles/PMC3985962/>
- Lisa, Palacios, Vila et al., **A general reaction mechanism for carbapenem
  hydrolysis by mononuclear and binuclear metallo-β-lactamases**, Nat Commun
  2017;8:538 — GOB-18 — <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5599593/>
- Morán-Barrio, Vila et al., **The metallo-β-lactamase GOB is a mono-Zn(II)
  enzyme with a novel active site**, JBC 2007;282:18286 —
  <https://www.jbc.org/content/282/25/18286.long>
- Horsfall et al., on GOB Met221 —
  <https://www.pure.ed.ac.uk/ws/files/8432406/FEBS_2010_Horsfall.pdf>
- Das & Nair, **NDM-1 QM/MM**, arXiv:1702.01239 — cited in D40 for weak
  carbonyl coordination, and **withdrawn for B3** in D123: L1's time-resolved
  snapshots measure 2.53–2.63 Å.

## Putting a template's ligand in our frame (D156)

```
R=/scratch/user/uquadhik/refs          # 7UHI.pdb, 5AXO.pdb, 5B15.pdb from RCSB
M=/scratch/user/uquadhik/mcube/docs/csr
F1=A:105:NE2=A:88:NE2,A:181:NE2=A:164:NE2,A:109:CG=A:92:CG,A:301:ZN=B:1:ZN
F2=A:72:NE2=A:88:NE2,A:150:NE2=A:164:NE2,A:76:CG=A:92:CG,A:301:ZN=B:1:ZN

python $M/refstruct.py transfer $R/7UHI.pdb CSR_wt/protdft/scan_aspin/scan_00.pdb \
    --ligand=MX0 --zn=B:1:ZN --fit=$F1 --poses=... --poseres=mer
python $M/refstruct.py transfer $R/5AXO.pdb CSR_wt/protdft/scan_aspin/scan_00.pdb \
    --ligand=LMP --zn=B:1:ZN --fit=$F2
```

The nine docked poses are in `CSR_ligand`'s frame, not `CSR_wt`'s, so they are
compared against a template transferred into **that** frame
(`CSR_ligand/CSR_nontruncated_1metal.pdb`, whose zinc is `B:2:ZN`) - the same
command with `$F3`, which is `$F1` with that selector. A template transferred
into `preprocess_results/finalsystem.pdb` instead puts the poses 70 A away:
that file is the re-centred solvated box and is not the docking frame.

The correspondences above were read off `refstruct.py shell` on each deposited
file, not taken from the papers. **5AXR is not the hydrolysed-meropenem
structure** (it is the 2-mercaptoethanesulfonate complex); 5AXO is.

## How alike are the B3 enzymes we borrow from (D158)

```
python refstruct.py similarity 7UHI.pdb:A 5AXO.pdb:A 6nc5.pdb:A
```

| pair | identity | core | core CA RMSD |
| --- | --- | --- | --- |
| L1 vs SMB-1 | 32.3% | 216/251 | 1.22 A |
| **L1 vs CSR-1** | **39.0%** | 221/249 | **1.05 A** |
| SMB-1 vs CSR-1 | 33.9% | 226/254 | 1.16 A |

CSR-1 is represented by 6NC5 chain A: its own sequence with the Q121H/K254H
mutations, and the only CSR-1 crystal carrying a metal. The identities are the
ordinary B3 spread; what the table is for is the comparison *between* pairs,
and L1 is closer to CSR-1 than it is to SMB-1.

## Where Asp92 sits, across CSR-1's own three crystals (D179)

`shell` needs a metal and the wild-type structures have none, so this is
`pair`: the same two atoms in each file, with no superposition, so the answer
cannot be an artefact of a fit.

```
for P in 6NC5 6DN4 6DQ2; do curl -sO https://files.rcsb.org/download/$P.pdb; done
python refstruct.py pair 6DN4.pdb A:116:NE2-A:192:NE2 A:116:NE2-A:120:OD1 A:116:NE2-A:120:OD2
python refstruct.py pair 6DQ2.pdb A:116:NE2-A:192:NE2 A:116:NE2-A:120:OD1 A:116:NE2-A:120:OD2
python refstruct.py pair 6NC5.pdb A:87:NE2-A:163:NE2  A:87:NE2-A:91:OD1   A:87:NE2-A:91:OD2
```

| structure | what it is | resolution | Asp B | NE2–OD1 | NE2–OD2 | His–His NE2 |
| --- | --- | --- | --- | --- | --- | --- |
| **6DN4** | **wild type**, metal-free | **1.99 Å** | **15–16** | **3.00** | **3.07** | 3.63 |
| 6DQ2 | HARLDQ **mutant S60**, metal-free | 2.16 Å | 53–66 | 4.28 | 5.49 | 3.75 |
| 6NC5 | **Q121H/K254H mutant**, Zn-bound | 1.60 Å | — | 4.13 | 4.82 | 3.42 |

6DN4 and 6DQ2 are in BBL numbering (His116, Arg118, Asp120, Gln121, Lys254),
so their Asp120 is our Asp92 and their His116/His192 pair is our
His88/His164. The 3.42–3.75 Å His–His NE2 separation in all three is the check
that the same site is being measured in each.

**The carboxylate is in two conformers and the split tracks what each structure
is:** close in the only well-ordered wild-type copy, far in a mutant whose
side chain is disordered (B 53–66) and in the construct that engineers a
zinc-binding histidine into the position next to it. The 4.20 Å figure quoted
from D145 onward as "the crystal's Asp92 distance" is 6NC5's — the engineered
one.
