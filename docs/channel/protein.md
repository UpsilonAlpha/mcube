# What the channel protein is

One-page summary of what is known about the protein in `sequence.fasta`, with where each fact
came from. The reasoning and the decisions it led to are in
[../channel-decision-log.md](../channel-decision-log.md) (entries C1, C2 and their follow-ups).
Last updated 2026-10-23.

## Identity

| | |
|---|---|
| Protein | **BD122_14570** (GenBank `APO51500.1`), *Bradyrhizobium diazoefficiens* USDA 122 |
| Match | blastp, 98.9 % identity over all 373 aa |
| Family | **DUF6662** (RefSeq name of every close relative); "MetA-pathway of phenol degradation" on a few UniRef90 relatives |
| Locus | `CP013127.1:3185355..3186476` (+) |
| True start | **Met73**: RefSeq relatives are 301 aa and align from query 73. Residues 1–72 are not part of the protein |
| Signal peptide | 73–109 (`MANAVEAVVVFGMRARAQ` n, `TFWIMAACTFALA` h, `APTCVHA│EE`): an outer-membrane protein |
| Mature chain | **110–373** (264 aa), the chain folded |
| Relatives | all Alphaproteobacteria, mostly *Bradyrhizobium*, also *Methylobacterium*, *Bosea*, *Rhodopseudomonas* |

## Genomic context: a lanthanide methylotrophy cluster

Within 20 kb, RefSeq names ([identity/neighbourhood_CP013127.tsv](identity/neighbourhood_CP013127.tsv)):

| Locus | RefSeq name | Position relative to BD122_14570 |
|---|---|---|
| BD122_14565 | **lanmodulin LanM** | 61 bp upstream, same strand: probably the same transcript |
| BD122_14530 | **lanthanide-dependent methanol dehydrogenase XoxF5** | 7 kb upstream |
| BD122_14640–14655 | **lanthanide ABC transporter PedA1 / PedB (LutE) / PedC** | 10 kb downstream |
| BD122_14560, 14500 | methanol/ethanol family PQQ-dependent dehydrogenases | |
| BD122_14515–14525 | glutathione-linked formaldehyde oxidation, methanol-metabolism cytochrome c | |

That covers a periplasmic Ln-binding protein, an inner-membrane Ln importer and a Ln-dependent
enzyme. The outer-membrane step is the one missing from the cluster, so this barrel is a
**candidate for Ln uptake across the outer membrane**. That is an inference from genomic
context; it does not say whether free Ln³⁺ or a chelated complex crosses.

LanM check by eye: 130 aa, signal peptide, four EF-hand loops `DTDNDG`, `DPDHDG`, `DPDRDG`,
`NPDKDG`. Three have the Pro at loop position 2 that marks lanmodulin.

## Structure (AF2 2.3.2)

| | Monomer (use this) | Trimer |
|---|---|---|
| Confidence | pLDDT 96.7, pTM 0.91 | ipTM 0.16, pTM 0.43: no confident interface |
| Fold | closed 12-stranded β-barrel, Cα radius 12.2 Å | same barrel in each chain; chains barely touch |
| Lumen | residues 110–132 (coil) thread into the lumen | same |
| Narrowest point | 0.8–1.1 Å radius over z −14…−5 Å, lined by Phe119, Met122, Ala286, Glu300 | 0.57 Å |

Structural relatives (AF2's HHsearch against PDB70): **5O65 FapF**, the *Pseudomonas*
functional-amyloid secretion channel (probability 98.2 %, query 121–350). It is a trimeric
β-barrel and its PDB entry is tagged "plug". Next is **4RL8**, a *P. putida* outer-membrane barrel
from the meta-pathway of phenol degradation (96.3 %).

## Conservation (AF2's MSAs: 2173 sequences, ~10 % mean identity to the query)

- **Glu300 is the most conserved lumen-facing residue in the barrel.** Identity is 0.59,
  against 0.12 for the strands on average. It lines the constriction.
- The lumen is net negative: 21 Asp/Glu face it, against 11 Arg/Lys.
- The N-terminal lumen segment is aligned in only 35 % of relatives, climbing from 1 % at
  Glu110 to 58 % at Glu132. Where present, Phe119 (at the constriction) is conserved (0.63).
  Its acidic residues are not (Asp126 0.16, Glu132 0.17).

## Where the data are

| What | Where |
|---|---|
| BLAST hits (top 100) | [identity/blast_clusterednr_2026-10-23.tsv](identity/blast_clusterednr_2026-10-23.tsv): blastp, ClusteredNR, default parameters, RID `DNGYGYKX016` (expires 2026-10-24) |
| Gene neighbourhood | [identity/neighbourhood.sh](identity/neighbourhood.sh) → `neighbourhood_CP013127.tsv` (rerun to refresh) |
| AF2 models, MSAs, `pore.json`, `conservation.json` | Bunya `/scratch/user/uquadhik/channel/af2/{monomer/mature,trimer/mature_trimer}/` |
| Pore and conservation scripts | [af2/pore.py](af2/pore.py), [af2/conservation.py](af2/conservation.py) |

## Open questions

- Does the N-terminal segment take part in ion transport? It forms the constriction together
  with Glu300. If it doesn't, crop residues 110–130 for the PMFs (user, 2026-10-23).
- Does free Ln³⁺ cross, and more easily than Ca²⁺? The Ca/La/Lu PMFs test this.
