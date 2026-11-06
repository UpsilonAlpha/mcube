---
type: observation
evidence: experimental
status: supported
entries: [D217]
checked: [D217]
---

# Our own dmCSR-1 S214D structure, at 1.27 A, puts the single zinc at Zn2 with Asp92 as its ligand, and Zn1 at half occupancy

**This project had never used its own crystal structure.** It is named in
`CSR_Manuscript-v1_old_gary_comments.docx` and the coordinates are in the
Honours tree: `.../SCIE4481/bioinformatic_data/MM/Structure/CSR/dmCSR_mutant_1metal_refine_12.pdb`
(1.27 A, one Zn) and `dmCSR_mutant_2metal_refine_16.pdb` (1.47 A, two Zn),
P2(1)2(1)2(1), Rwork/Rfree 0.163/0.179 and 0.188/0.214, accession still pending.
Numbering is BBL **-42**: Asp78 [92], His79 [93, from Gln], His212 [226, from
Lys], His74 [88], His150 [164].

    ours1  ZN B1 occ 0.70  HOH2 2.07, His212 2.12, Asp78 OD2 2.13, HOH1 2.16,
                           His79 2.22, HOH61 2.38          <- the Zn2 (DHH) site
    ours2  ZN B1 occ 1.00  Asp78 OD2 2.06, His79 2.11, His212 2.18, + 3 waters
           ZN B2 occ 0.53  His150 2.05, His74 2.23, + 3 waters   <- Zn1, HALF occupied

The manuscript's own reading: only one zinc in the active site across multiple
crystallisation conditions, at Zn2; the second site occupied only in the
ampicillin co-crystal; *"a possible mechanism where binding of metal at this site
is facilitated by the presence of substrate."*

**What it cannot do.** It is the double mutant - position 226 is a histidine -
so **there is no lysine in it** and it is not a wild-type reference for anything
about Lys226. Nor does S214D have a recorded effect at the site.

**What it does do, twice over.** It confirms in-house and at high resolution that
**an available Zn2 ligand set takes the carboxylate at ~2.1 A and leaves Zn1 as
the weak site** - the exact mirror of the wild type, where there is no Zn2 set
and the carboxylate goes to the lysine instead
([[zinc-bound-crystal-holds-asp92-off-the-metal]]). And it supplies a better
explanation for the metal-free wild-type crystals than crystallisation
conditions: metal binding at Zn1 may need the substrate
([[wild-type-csr1-binds-no-zinc-in-any-crystal]]).

**It also sharpens the assembly problem.** Five deposits and our own structure,
and not one of them shows the state the project models - a single zinc at Zn1,
Arg90 intact, Lys226 intact. Ours puts the metal at the other site entirely.

## Supports
- [[csr1-is-mono-zinc-with-zn2-vacant]]
- [[wild-type-csr1-binds-no-zinc-in-any-crystal]]
