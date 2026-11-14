---
type: observation
status: supported
evidence: ours
entries: [¶247, ¶282]
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

## Supports
- [[the-cpha-conditions-are-met-in-the-crystal-derived-pose]]
- [[arg90-stands-at-the-third-zn1-ligand-position]]
