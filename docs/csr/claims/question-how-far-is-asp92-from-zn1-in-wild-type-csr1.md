---
type: question
status: open
entries: [¶145, ¶151, ¶179, ¶184, ¶208, ¶284]
---

# How far is Asp92 from Zn1 in wild-type CSR-1?

**¶284 ANSWERS IT: 3.00 Å (OD2), 4.12 (OD1).** 6NC5 superposed onto
metal-free wild-type 6DN4 - the same protein, **core 257/258, RMSD 0.55 Å** -
places the M1 metal 3.00 Å from Asp120 OD2, with Lys254 NZ at 4.39 and
Arg118 NH1 at 3.66. It is a long contact, not a bond, and it is the number
this question asked for
([[neither-aspartate-coordinates-the-metal-in-our-own-crystals]]).

**Why this one is not ¶184 again.** ¶184 was voided because it superposed into
an apo structure and then quoted the result as a measurement. This is the same
construction and carries the same caveat - the metal is placed, not observed -
but three things have changed: the fit is reported (0.55 Å over 257 of 258
residues, so it is not a free parameter), the geometry is **not** our own model
but one crystal into another, and the answer is used only for the comparison it
supports, namely that the aspartate is CLOSER in the wild type (3.00) than in
the lysine-free mutant (4.12). **It is still not an observed wild-type
Zn-Asp distance, and no such observation exists.** The question therefore stays
`open` with a number in it rather than closed, which is the honest state: the
comparison it supports (3.00 in the wild type against 4.12 without the lysine) is
robust to the placement, because both ends of it are placed the same way.

**Unmeasured, and no deposited structure constrains it.** Three readings have
been quoted as "the crystal's" value and none of them is one:

- **2.80 A** - our own model's number. It appears nowhere in any PDB file
  (¶208). It was nonetheless the restraint target, the gate value and the
  reasoning baseline for five weeks.
- **2.68 A** - ¶184's superposition of 6NC5's zinc into 6DN4. **Void**: 6DN4
  contains zero metal atoms (¶208), so this is our own geometry, not a
  measurement.
- **4.20 A** - real, measured in 6NC5, but 6NC5 is the Q93H/K226H double mutant
  **with Zn2 occupied**, and that Zn2 ion is what holds the carboxylate
  ([[asp-zn1-distances-are-set-by-the-zn2-ion]]).

The structural reason the question cannot be closed from the record: every
structure that pins the number has a Zn2 doing the pinning, and the only
Zn2-free wild-type CSR-1 structures, 6DN4 and 6DQ2, are **apo** - 4557 and 4049
atoms, zero metal sites.

**What would settle it.** A wild-type CSR-1 structure with Zn1 occupied and Zn2
empty, which does not exist; or, failing that, a thermodynamic calculation
whose geometry is not handed to it (plan item 11, deferred by design). Note
that ¶206 rule 3 excludes our own minima from answering it, which is exactly
why three campaigns on this residue (¶51-¶59, ¶120-¶145, ¶173-¶200) each
reopened.

**What is settled, and is enough for the decision:**
[[asp92-is-never-a-zn1-ligand-in-b3-structures]]. The qualitative statement
does not need the distance.

## Supports
- [[gfn2-over-coordinates-zn1-from-a-crystal-consistent-start]]

## Supported by
- [[wild-type-csr1-crystals-contain-no-metal]]
