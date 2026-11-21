---
type: observation
status: supported
evidence: ours
entries: [¶281, ¶282, ¶303]
checked: [¶281, ¶282, ¶303]
---

# Ser186's hydroxyl donates a 1.70 angstrom hydrogen bond to the zinc-bound nucleophile

**RESTORED (¶303). The test this node pre-registered has been run, and it
passes.** ¶282 voided the contact because Ser186 was an MM atom at the QM/MM
boundary and said in as many words: *"add Ser186 to the QM region and re-run. If
the contact survives at 1.8-2.0 Å with the residue treated quantum-mechanically,
this node can be restored."*

**It survives.** No new run was needed - `scan_posepull` uses
`qm_residues: *qmreg_michrks`, the region that includes Ser186 (¶287), so its
frames already satisfy the condition. Measured on `scan_00.pdb`:

    Ser186 OG - O01          2.767        (MM-boundary value was 2.672)
    Ser186 HG ... O01        **1.843**    (MM-boundary value was 1.696)
    OG - HG ... O01 angle    155.1 deg
    Asp185 OD1 - O01         4.400
    Ser186 OG - N1x          4.667

**1.843 Å is inside the 1.8-2.0 Å window this node set as its own restoration
criterion**, and the shift from 1.696 is in exactly the predicted direction:
Pauli repulsion from the QM density lengthens a contact that electrostatic
embedding had let collapse. The angle at 155.1 deg is a real if slightly bent
hydrogen bond.

**And the fork this node described is now visible inside ONE trajectory.** Along
`posepull`, as the substrate carboxylate is pulled onto the metal, Ser186 lets go
of the nucleophile and follows the carboxylate:

    scan_00   HG...O01 1.843   OG-O4x 6.54
    scan_10   HG...O01 3.628   OG-O4x 5.12   (OG-HG...O01 angle 108.2 deg)

So the residue does not do both jobs at once in our system either - which is the
table below, measured rather than inferred across two poses of one enzyme.

**Caveats that still stand, and they matter before this is quoted.** GFN2 at 0 K,
running ~0.13 Å short on Zn-O against DFT (¶245), so the number needs
re-measuring at def2-SVP. `scan_00` is the scan's first point, so it is a
constrained optimisation - though the constrained coordinate is Zn-O4x at its
existing 5.42 Å, far from either atom here. And the pose remains 6.11 Å from the
experimental mode.

**Why this now matters more than when it was voided.** ¶302 withdrew the argument
that Asp92 at 7.2 Å condemns our pose - but it also confirmed that Asp92 cannot
be the nucleophile's partner in a mono-zinc model, because its closeness in every
reference enzyme comes from ligating the absent Zn2. That leaves the question of
what positions the nucleophile here **open and unanswered**, and this contact is
the only candidate our own structures contain. It matches in kind the one
published mono-zinc MBL mechanism
([[the-published-mono-zinc-mechanism-has-no-metal-activated-nucleophile]]).

**Superseded, below, is now too strong**: the crystal-derived pose gives Ser186
the carboxylate job, but that pose is itself a transfer whose key distance is an
artefact to ~1 Å, so the two are rival readings rather than one replacing the
other.

**The original void reasoning, kept for the record:**

**VOID (¶282). SER186 IS NOT A QM RESIDUE.** `aspxtb8`-s QM region is nine
residues - Zn, OH, mer, His88, His164, Asp92, Asn225, Arg90, Lys226 - and Ser186
is not among them, nor in `extra_qm_residues`. So this contact is **a
fixed-charge MM hydroxyl hydrogen pointing at a QM hydroxide across the QM/MM
boundary**, where electrostatic embedding lets the MM charge feel the anion-s
field with no Pauli repulsion from the QM density.

**The number is the tell.** A real O-H...O hydrogen bond runs **1.8-2.0 Å**.
1.696 Å is shorter than that - not a strong hydrogen bond, an unphysical one. I
quoted its shortness as evidence of strength; it is evidence of the artefact.

**Kept rather than deleted, because the test it implies is cheap and worth
running:** add Ser186 to `extra_qm_residues` and re-run `aspxtb8`. If the contact
survives at 1.8-2.0 Å with the residue treated quantum-mechanically, this node
can be restored. **Not submitted.**

**And in the pose with experimental backing, Ser186 does something else
entirely** - it anchors the substrate carboxylate at 2.72 / 3.00 Å
([[the-crystal-derived-pose-puts-the-carboxylate-on-the-metal]]), which is the
Ser221 role in SMB-1 and L1 and needs no QM region to be believed.


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

## Superseded by
- [[the-crystal-derived-pose-puts-the-carboxylate-on-the-metal]]
