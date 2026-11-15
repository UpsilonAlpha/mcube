---
type: question
status: supported
entries: [¶145, ¶151, ¶179, ¶184, ¶208, ¶284, ¶285]
checked: [¶145, ¶151, ¶179, ¶184, ¶284]
---

# How far is Asp92 from Zn1 in wild-type CSR-1?

**¶285 ANSWERS IT BY DIRECT MEASUREMENT: ~4.2 Å, AND THE LYSINE DOES NOT
CHANGE IT.** Measured in place, in the structures that actually contain metals,
with no superposition of any kind:

    6NC5 ZN302 (M1, no lysine)      Asp91  OD2  **4.20**
    refine_16 ZN B2 (M1, no lysine) Asp78  OD1  **4.12**
    6DQH ZN A401 (M1, LYSINE)       Asp120 OD1  **4.28**

and where the aspartate does hold a metal it holds **its own**, at the M2 site,
at 1.95-2.13 Å
([[the-zn2-aspartate-is-42-angstroms-from-the-m1-metal-with-or-without-the-lysine]]).

**¶284's answer of 3.00 Å is withdrawn.** It came from superposing
metal-bearing 6NC5 onto metal-free 6DN4 and carrying a metal across. 6NC5 has
its own metals, so that fit was never needed, and in place the distance is
**4.20 Å**. The comparison ¶284 built on it - that the aspartate is *closer*
with the lysine present - goes with it: 4.28 with the lysine against 4.12-4.20
without is no difference at all.

**So this closes as `supported`, on experimental evidence and without a
placed metal.** What ¶206 rule 3 excluded was our own minima answering it; this
answer is four crystals measured in their own frames. The caveat that remains
is not about the number but about its *cause*: in 6NC5 and `refine_16` the
aspartate is held out because it is busy coordinating an M2 metal, and no
structure exists with the lysine absent **and** the M2 site empty, which is the
one case that would isolate the lysine's contribution
([[thesis-lys226-restrains-asp92-in-place-of-zn2]]).

**The three readings previously quoted as "the crystal's" value, none of which
was one.** Kept because each was load-bearing for weeks and the way each
failed is the reason ¶285 insisted on measuring in place. The heading here used
to read "unmeasured, and no deposited structure constrains it" - that is now
wrong: 6NC5 and 6DQH constrain it directly, which is what closed this
question.

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
- [[the-zn2-aspartate-is-42-angstroms-from-the-m1-metal-with-or-without-the-lysine]]
