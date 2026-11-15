---
type: observation
status: supported
evidence: ours
entries: [¶285]
---

# Our Arg90 reproduces the crystal arginine to within 0.65 angstroms, and it is the metal that has moved

**The worry.** Measured against each structure's own metal, the HRH arginine is
**3.67-4.05 Å** from the M1 metal in every crystal that has one (6NC5 Arg89 NH1
**4.02**; `refine_16` Arg76 NH1 **3.67** to the M1-site metal), while in ours it
is **5.50** (`relax10`) and **5.76** (`aspxtb8`). A 1.5-1.8 Å discrepancy in the
one side chain the carbonyl-polarisation argument rests on.

**It is not the arginine.** Fitting 6NC5 onto `relax10` (core **258/259**, RMSD
**0.79 Å**) and carrying its Arg89 across:

    atom   crystal Arg in our frame     our own Arg90      difference
           (to Zn / to ligand O1)       (to Zn / to O1)
    NE          6.00 / 4.05                6.04 / 4.27        **0.38**
    CZ          6.23 / 3.89                6.19 / 4.01        **0.39**
    NH1         5.43 / **2.88**            5.50 / **2.98**    **0.33**
    NH2         7.51 / 5.09                7.38 / 5.19        **0.65**

Every atom lands **within 0.65 Å**, and NH1 - the atom that does the work -
within **0.33 Å**, well inside the fit's own 0.79 Å RMSD. **Our arginine is
already in the crystal position.** The 5.50-against-4.02 difference is a
statement about the **metal**, which sits ~1.4 Å from where 6NC5's M1 metal
sits, consistent with
[[our-single-metal-sits-between-the-canonical-zn1-and-zn2-positions]].

**Its partners are conserved, which is the independent check:**

| contact | 6NC5 | `refine_12` | `refine_16` | ours (`aspxtb8`) |
| --- | --- | --- | --- | --- |
| Arg NE - Asp OD2 | 2.77 | 2.83 | 2.88 | **2.91** |
| Arg NH2 - Asp OD1 | 2.90 | 3.00 | 2.99 | **2.84** |
| Arg NH1 - Tyr OH | 2.95 | 2.83 | 2.72 | **3.00** |

The aspartate is our **Asp125** (6NC5 124, refine 111, 6DQH 153 - an aspartate
in all four) and the tyrosine our **Tyr192**. A buried bidentate salt bridge plus
a tyrosine hydrogen bond, reproduced to 0.2 Å, is why the rotamer did not
drift. Note this is **not Asp92**: the arginine's salt-bridge partner is a
different, conserved aspartate, and that is what pins it.

**One contact does differ, and it is the histidine's doing.** Arg NH1 - His(M1)
NE2 is **3.35 Å** in 6NC5 and **4.37-4.46** in ours. Since the arginine's NH1
matches to 0.33 Å, the atom that moved is **His164 NE2**, following the metal.
The single discrepancy in the site is downstream of the metal position, not of
the arginine.

**What this does not establish.** No crystal has a metal at M1 and none at M2,
so our one-metal arrangement is still unverified - the metal's position is
exactly what differs, and nothing here tests it. And the carbonyl contact (2.88
with the crystal rotamer, 2.98 with ours) inherits the ligand pose
([[the-crystal-derived-pose-puts-the-carboxylate-on-the-metal]]).

**No precedent for arginine as a direct zinc ligand.** A literature search found
arginine described only as second shell - hydrogen bonding to zinc-bound ligands
and influencing site geometry - and never in a catalytic zinc's donor set, which
is histidine, glutamate, aspartate and cysteine plus solvent. So the expectation
that the arginine should *approach* the metal was never well founded: 3.67-4.05
Å is a second-shell distance in the crystals too.

## Supports
- [[arg90-is-on-the-carbonyl-and-cannot-reach-the-leaving-nitrogen]]
- [[our-single-metal-sits-between-the-canonical-zn1-and-zn2-positions]]
