---
type: claim
status: supported
entries: [¶285]
---

# Ser186 takes over Ser221's job of positioning the substrate carboxylate

Three measurements, none of which needs the ligand pose:

1. **The residue is in the canonical place.** L1 Ser206 OG lands **1.44 Å** from
   our Ser186 OG (core 221/258, RMSD 1.22) and SMB-1 Ser175 OG **1.64 Å** (core
   219/254, RMSD 1.25) - inside each fit's own error
   ([[ser186-is-the-canonical-carboxylate-anchor]]). This settles the S214 /
   S221 numbering question from coordinates rather than from an offset.
2. **The job exists to be taken over.** Canonical di-zinc B3 binds the substrate
   carboxylate on a metal and anchors it with that serine
   ([[canonical-b3-puts-the-carboxylate-and-the-leaving-nitrogen-on-the-same-metal]]).
3. **Removing it costs activity, by occupying the site rather than the metal.**
   The S214D aspartate lands on the anchor position itself and never comes
   within 4.06 Å of any metal ([[s214d-occupies-the-substrate-carboxylate-site]],
   [[neither-aspartate-coordinates-the-metal-in-our-own-crystals]]).

**So the anchoring function is conserved and CSR-1 keeps it**, while the second
canonical anchor is not: L1 Ser208 / SMB-1 Thr177 have no counterpart, that
position being backbone in our frame (Ala188 N 2.44 Å, Ala188 O 2.78). CSR-1
has **one of the two**.

**What "takes over" does and does not mean.** It means the serine does in CSR-1
the positioning job it does in canonical B3 - same residue, same place, same
consequence when mutated. It does **not** mean the serine substitutes for the
missing second zinc; that role is split between the metal and Arg90
([[the-metal-can-both-carry-the-nucleophile-and-stabilise-the-leaving-nitrogen]]).

**The weak point is the distance, not the identity.** Ser186 OG to the
substrate's carboxylate oxygens is **2.72 / 3.00 Å** in the crystal-derived
pose, and that pose is transferred, not equilibrated
([[no-crystal-derived-pose-has-ever-been-equilibrated]]). Ser186 is also **not a
QM residue**, so its contacts in `aspxtb8` are MM-only and were already ruled
untrustworthy. The claim rests on where the residue is, which is solid, and not
on how long the hydrogen bond is, which is not.

## Supports
- [[thesis-mechanism-differs-from-canonical-b3]]

## Supported by
- [[ser186-is-the-canonical-carboxylate-anchor]]
- [[s214d-occupies-the-substrate-carboxylate-site]]
- [[canonical-b3-puts-the-carboxylate-and-the-leaving-nitrogen-on-the-same-metal]]

## Requires
- [[ser186-is-the-canonical-carboxylate-anchor]]
