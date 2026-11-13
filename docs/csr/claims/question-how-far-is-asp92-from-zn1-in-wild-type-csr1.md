---
type: question
status: open
entries: [¶145, ¶151, ¶179, ¶184, ¶208]
---

# How far is Asp92 from Zn1 in wild-type CSR-1?

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
