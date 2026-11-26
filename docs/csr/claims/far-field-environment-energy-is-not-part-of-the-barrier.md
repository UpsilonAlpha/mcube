---
type: observation
status: supported
evidence: literature
entries: [¶331, ¶334, ¶341]
checked: [¶344]
---

# The methodology literature judges a far-from-site environment term irrelevant to the barrier, and prices freezing it

Read rather than re-derived. This is the criterion by which ¶344 settled ¶332,
and it matters that it is **a criterion about distance**, not about magnitude: a
large environment term is not by itself an error, and a small one is not by
itself safe.

**Ryde, U.** *J. Chem. Theory Comput.* 2017, 13, 5745. DOI
10.1021/acs.jctc.7b00826. The distance criterion, verbatim: *"the addition or
removal of a single hydrogen bond far from the active site affects the energy by
~20 kJ/mol in a way that most likely is irrelevant for the reaction."* And the
mechanism behind the scatter: *"it is very hard to ensure that the surrounding
enzyme remain in the same local minima throughout a reaction sequence"* - the
**local-minima problem**. Table 1 compiles the standard deviation of QM/MM
energies across independently minimised snapshots from 24 published studies:
*"s ranges from 0.6 to 97 kJ/mol and 73% of the studies have at least one series
of energies with s > 10 kJ/mol."*

**Calixto, A. R.; Ramos, M. J.; Fernandes, P. A.** *J. Chem. Theory Comput.*
2017, 13, 5486. Prices the remedy: freezing everything outside a free region
costs **0.8 kcal/mol on average and 1.6 kcal/mol maximum** once that free region
reaches **6.00 Å**. Below 6 Å the structure is overconstrained; at or above it,
freezing is cheap in error terms. (Read from the abstract; the full text is
closed.)

**Lonsdale, R.; Harvey, J. N.; Mulholland, A. J.** *Chem. Soc. Rev.* 2012, 41,
3025. The prohibition, for the case where relaxation carries the profile onto a
conformationally distinct path: *"When this occurs, the energies obtained should
not be used."* The same review prescribes repeating adiabatic mapping in forward
and reverse directions as the diagnostic.

**Sun & Ryde 2025** report ΔE‡ of **44 kJ/mol with a relaxed environment against
70 kJ/mol with a fixed one**, and state that an MM term dominating the barrier
*"indicates that the QM region is too small or that the MM region is not enough
relaxed."*

**Kästner, J.; Sherwood, P.** *Mol. Phys.* 2010, 108, 293. DOI
10.1080/00268970903446764. The direct precedent for an environment term
exceeding the barrier it sits in: the environment-interaction term **83.3 kJ/mol
exceeds the entire ribosome-catalysed barrier of 73.1 kJ/mol**, the barrier
spread across 10 snapshots is 51.5 and 71.6 kJ/mol, and single MM residues
contribute +39.4 and −14.4 kJ/mol. Their sentence closest to ours: *"The
geometric data of the atoms directly involved in the transition are uncorrelated
with the activation barrier... the differences in the activation barriers between
the snapshots must be caused by interactions with the environment."*

**The field's four remedies**, from the same body of work: report free energies;
keep the surroundings fixed; run the coordinate forward *and* backward; average
over MD snapshots.

## What this does not say

It does not license dropping an MM term because it is large. Every spread
verified above is **across independently minimised snapshots drawn nanoseconds
apart**, which is a different object from a single optimisation at fixed
chemistry — that is why ¶344 had to measure the radial distribution rather than
cite a magnitude. And no source here says the MM term cancels between states:
the citation once used in this project for that (Ganguly, Boulanger & Thiel 2017)
says something else and is listed in `retractions.md`.

It is also not evidence about anyone's protocol. ¶334's claim that Khrenova &
Nemukhin 2018 and Levina 2020 use our scan → saddle → frequencies protocol is
**unverified and withdrawn**: Levina defers its whole protocol to the closed 2018
paper, which states no MD length, equilibration, restraints, frozen atoms or
active radius.

## Supports
- [[the-barrier-is-scored-along-the-scan]]
- [[our-stage-energies-are-not-on-a-common-surface]]
