---
type: observation
status: supported
evidence: ours
entries: [¶281]
checked: [¶281]
---

# Ser186's hydroxyl donates a 1.70 angstrom hydrogen bond to the zinc-bound nucleophile

Measured in `aspxtb8`, an unrestrained GFN2 optimisation with the second shell
free (479 active atoms, nothing constrained):

    Ser186 OG  - O01 (nucleophile)    2.672      HG...O01   1.696
    Ser186 OG  - N1x (leaving N)      4.698
    Ser186 OG  - O4x (substrate COO)  6.464
    Asp185 OD1 - O01                  4.367

**1.696 Å H...O is a strong, well-formed hydrogen bond**, with the serine
hydroxyl proton pointing at the hydroxide. Nothing in the setup put it there -
no restraint on either atom, and the shell was free to relax.

**It is a fork, not just a contact.** In `l1_in_relax10.pdb`, the
experiment-derived pose (7UHI's ligand transferred at fit RMSD 0.310 Å),
**Ser186 OG sits 2.72-3.00 Å from the substrate C3 carboxylate**. In `aspxtb8`
it is 6.464 Å from that carboxylate and 2.672 Å from the nucleophile. The two
poses give the residue two mutually exclusive jobs:

| pose | Ser186's role |
| --- | --- |
| experiment-derived, 7UHI transferred | anchors the substrate C3 carboxylate, as Ser221/Thr223 do in SMB-1 and L1 |
| `aspxtb8` | hydrogen-bonds the nucleophile |

**Why it matters mechanistically.** A zinc-bound hydroxide attacking a carbonyl
has to put its proton somewhere, and this project has had no answer that did not
need Asp92 as a general base - a route whose premise is refuted
([[cpha-route-is-closed]]) - or a second metal-bound water. Ser186 at 1.70 Å is
a third option needing neither, and it is **what the only published mono-zinc
MBL mechanism predicts in kind**
([[the-published-mono-zinc-mechanism-has-no-metal-activated-nucleophile]]):
nucleophile activation by a hydrogen-bond network rather than by the metal. Ours
differs in being on the metal as well.

**Caveats carried, not buried.** GFN2 at 0 K, and GFN2 runs ~0.13 Å short on
Zn-O against DFT (¶245), so this hydrogen bond needs re-measuring at def2-SVP
before it is quoted. It also sits inside a pose 6.11 Å from the experimental
mode ([[aspxtb8-is-substrate-bound-but-not-a-near-attack-complex]]), which is
why this is an observation and not yet a claim.

## Supports
- [[the-likely-route-is-a-zn1-hydroxide-with-ser186-as-proton-partner]]
