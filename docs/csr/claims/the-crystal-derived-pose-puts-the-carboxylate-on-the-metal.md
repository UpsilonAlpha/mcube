---
type: observation
status: supported
evidence: ours
entries: [¶247, ¶282, ¶298]
checked: [¶282]
---

# In the crystal-derived pose the metal binds the substrate carboxylate, Arg90 polarises the carbonyl, and the metal-bound hydroxide is 4.7 angstroms from the scissile carbon

7UHI-s moxalactam transferred into `relax10`-s frame at fit RMSD 0.310 Å (¶247),
measured against `relax10`-s protein:

    Zn - carboxylate O8        **1.847**     the metal binds the CARBOXYLATE
    Zn - lactam carbonyl O1      3.153       NOT coordinated
    Zn - lactam N1               3.254
    Zn - our hydroxide O01       2.035

    Ser186 OG - O7 **2.72** / O8 **3.00**    anchors the carboxylate
    Arg90 NH1 - carbonyl O1 **2.98**         POLARISES THE CARBONYL
    Asp92 OD1 - Zn **4.33**  N1 **4.52**     room for a water between them
    Asp185    - Zn   5.78    N1   8.30       too far to be a candidate
    Lys226 NZ - N1 **5.830**                 further than in aspxtb8
    our hydroxide O01 - MX0 C1 **4.702**     the metal-bound OH cannot attack

**This is a different enzyme from the one `aspxtb8` describes**, built from the
same protein and an experimentally placed ligand.

**The substrate is placed and polarised without the carbonyl touching the
metal:** the carboxylate-metal bond at 1.847 Å places it, Ser186 holds the same
carboxylate, and **Arg90 NH1 hydrogen-bonds the lactam carbonyl at 2.98 Å** -
a cationic polariser where canonical B3 uses a second zinc. **Arg90 is a QM
residue** in every section that matters, so this contact does not have the
boundary problem that voided
[[ser186-hydrogen-bonds-the-nucleophile-in-our-reactant]]. It is a direct
structural statement of what the RQK arginine is for.

**Two routes are excluded by it, one is opened.** The metal-bound hydroxide sits
**4.702 Å** from the scissile carbon, so a Zn-OH nucleophile is not available
here. Lys226 is **5.830 Å** from the leaving nitrogen - slightly further than
in `aspxtb8` - so the canonical route-s substitute is not available either. What
is available is a water in the Asp92 gap
([[the-cpha-conditions-are-met-in-the-crystal-derived-pose]]).

**Qualification on Lys226, which is the one measurement here that could move.**
Its NZ position comes from the `asppull` actuator lineage and a lysine side chain
is mobile, so 5.830 Å is **rotamer-dependent, not a hard limit**. What has not
been tested is whether another rotamer reaches N1 while keeping the Asp92 salt
bridge; a side-chain scan would settle it cheaply.

**Nothing in this pose has been solvated or equilibrated** -
[[no-crystal-derived-pose-has-ever-been-equilibrated]] - so these are
single-structure distances, not populations.

**¶284 ADDS THE CAVEAT THAT MATTERS: this contact is partly inherited from the
template.** Our single metal sits **1.77 Å** from where 7UHI's **Zn2** lands in
our frame ([[our-single-metal-sits-between-the-canonical-zn1-and-zn2-positions]]),
and in 7UHI it is Zn2 that holds the carboxylate at 2.65 Å
([[canonical-b3-puts-the-carboxylate-and-the-leaving-nitrogen-on-the-same-metal]]).
So O8 arriving at 1.847 Å from our metal follows from the transfer rather than
from a discovery. Carrying MX0 across by a global CA fit instead of ¶156's
four-atom fit moves it by mean **1.03 Å** (max 1.44) and gives **Zn-O8 2.75**,
Zn-O1 3.82, against the stored 1.847 / 3.153 - the same arrangement, a looser
metal contact. The pose is still the right thing to equilibrate; it is not
independent evidence that the carboxylate binds.

**¶286: THE LIGAND IN THIS POSE IS NOT MEROPENEM.** It is MX0, an oxacephem
(C13 H16 N2 O8, no sulfur, 23 heavy atoms), and ¶156-s transfer script simply
renamed it `mer`
([[the-crystal-derived-pose-is-an-oxacephem-not-meropenem]]). The pharmacophore
atoms every argument here used - O8, O1, N1, C1 - are shared between the
classes, so the measurements stand as a **template** for where a bound
substrate-s carboxylate, carbonyl and leaving nitrogen sit. What does not
stand is the phrase "the crystal-derived pose of our substrate": there is no
such pose yet, and making one is a construction.

**So what IS independent about the carboxylate site?** The question is worth
separating, because three different things have been argued from this pose and
they do not have the same standing:

| statement | standing |
| --- | --- |
| canonical B3 binds the substrate carboxylate **on a metal** | **experimental** - 7UHI 2.65 Å, 5B15 2.23 Å, measured in place |
| **Ser186** is the residue that anchors it, and it is in the canonical place | **experimental** - protein-only fits, L1 Ser206 1.44 Å, SMB-1 Ser175 1.64 Å; no ligand involved |
| the carboxylate is **1.847 Å from our metal** | **inherited** - a transfer artefact to within ~1 Å; a global fit gives 2.75 |

The first two do not depend on this pose at all. **Only the contact distance
does**, and it is the one number that should never be quoted as a result.

**Which also answers how the substrate carboxylate-s position is known at all:
it is not, to better than about 1 Å.** What is known to better than that is
where the *anchoring residue* sits, from two deposits and the protein backbone
alone. So arguments should be run through Ser186-s position, which is measured,
rather than through the ligand oxygen-s, which is placed.

**¶298: THE ARRANGEMENT HAS NOW BEEN DRIVEN IN OUR OWN SYSTEM, AND THE CARBONYL
DOES NOT LEAVE.** `scan_posepull` (job 33244885, 11 points, all converged) drove
**Zn-O4x from 5.42 to 2.20 Å** at GFN2 QM/MM from the `aspxtb8` minimum, with
the carbonyl left completely free and **nothing restraining it** - the test ¶288
pre-registered, on the reasoning that a two-histidine zinc cannot hold the
carboxylate, the carbonyl and a nucleophile at once.

It held all three. At the endpoint: hydroxide O01 **1.78**, His164 **2.17**,
carboxylate O4x **2.20**, His88 **2.22**, carbonyl O2x **2.24** - five donors.
**The carbonyl moved 0.007 Å across the whole scan** (2.231 to 2.238) while the
carboxylate came in 3.2 Å. The attack geometry went the wrong way: O01-C4x
**3.457 -> 3.992 Å**, O01-C4x-O2x **45 -> 23 deg**.

So the failure is not that the carboxylate cannot reach our metal - it reaches it
at the template's own distance. **It is that its arrival displaces nothing.**
This bears on the pose's viability in our system and not on the first two rows
of the table above, which are experimental and do not depend on it.

**What this is not.** It is a path of constrained minima along a coordinate we
invented in order to build the structure, with the environment beyond 8 Å frozen,
so its energies are not a barrier or a binding penalty (¶297). The geometric
readings are the result.

## Supports
- [[the-cpha-conditions-are-met-in-the-crystal-derived-pose]]
- [[arg90-stands-at-the-third-zn1-ligand-position]]
