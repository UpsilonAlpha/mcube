---
type: claim
status: supported
entries: [¶235, ¶308, ¶332, ¶341, ¶344, ¶346]
checked: [¶344, ¶346]
---

# Leg 1's barrier is scored along the scan, and no freely optimised structure enters it

¶332 asked what environment the published barrier is scored in, and carried it as
an open decision for ten entries. **¶344 settled it by measuring where the
contaminating energy lives**, which is the one thing nobody had done: every
argument for freezing an outer region is an argument about *distance* from the
active site, and the distance had never been measured.

## The measurement

`docs/csr/shell_shift.py` on the two pairs that bracket the problem, 47,979 atoms
each, binned by distance from Zn1. `sum d²` is where the strain is:

    shell/Å       scan_00 → scan_11      scan_11 → saddle1 final
                  (chemistry moves)      (chemistry FIXED at 2.02 Å)
     0-3                  10.32                     0.00
     3-6                   5.11                     0.01
     6-8                   1.76                     0.18
    8-12                   0.75                   **1.55**
    12-16                  0.28                     0.03
    >16                   0.000                    0.000
    -------------------------------------------------------
    inside 6 Å              85%                     0.6%
    in 8-12 Å                4%                      88%

**Radially inverted.** The reaction moves the active site and decays outward; the
free optimisation moves a hydrophobic cluster 8.8-10.4 Å out - Trp21 0.517 Å,
Leu42 0.136, Phe18 0.105, Asn41 0.062 - while **every reacting atom holds under
0.031 Å**: the nucleophile 0.014, the zinc 0.010, His164 0.010, Ser186 0.017,
Asp92 0.027, His88 0.006. The 46,160 atoms beyond 16 Å moved exactly 0.000 Å in
both comparisons, which is the frozen set and an independent check that the atom
ordering is consistent.

## Why that settles it

The far field **does not move with the reaction**: at most 0.246 Å in the 8-12 Å
shell across the whole 1.60 Å change in the coordinate, against 1.657 Å inside
3 Å. So whatever MM strain it carries is the same at both ends of the scan and
**cancels in the difference by construction** - a displacement measurement, not
an assumption about cancellation. A free optimisation releases that strain at one
point only, which is the entire defect behind `saddle1` finishing 6.38 kcal/mol
below its own reactant ([[our-stage-energies-are-not-on-a-common-surface]]).

So the fix is to stop differencing a freely optimised point against a scan point,
not to freeze anything. A frozen-splice construction is what ¶332 itself flagged
as introducing its own bias, and it would be solving a problem the scan does not
have. ¶235 reached the same conclusion one stage down, choosing `actradius: 0`
because "the frozen term is a constant".

It also meets the literature on the literature's own terms. Ryde (JCTC 2017, 13,
5745) is explicit that the irrelevant term is the distant one; ours is 9-12 Å out.
Calixto (JCTC 2017, 13, 5486) measures the cost of freezing outside a 6.00 Å free
region at 0.8 kcal/mol average, 1.6 maximum, and the chemistry here is inside
6 Å while we free 8.

## What it specifies

    reactant      `scan_01`, -0.293 kcal/mol - the profile's own lowest point,
                  with `scan_00` scored too so the reference is visible
    maximum       `scan_11`, 14.237 kcal/mol at O01-C4x 2.0302 Å
    intermediate  `scan_15`, 0.561 kcal/mol at 1.45 Å, also `ringopen`'s input,
                  verified to 0.065 kcal/mol
    barrier       **14.53 kcal/mol at GFN2**; the DFT value comes from
                  jobs 33876755/56/57/58, submitted and running at ¶346,
                  as `bar11 - bar01`
    level         identical across all four and therefore INHERITED, never
                  stated per section: a barrier is a difference

## What it costs, and what it retires

The barrier is a potential-energy barrier from **adiabatic mapping**, not a
saddle verified by frequencies: ¶308's "the maximum of a driven scan is not a
transition state" stands and becomes a limitation to state
([[no-reaction-coordinate-has-ever-been-computed]]). ¶308's "it is a total, not a
QM energy" is answered - that was ¶332's question. ¶308's "not quotable at this
level of theory" is what the single points are for.

Retired for leg 1: the `ts` stage, as the wrong instrument rather than an
unconverged run; the frozen-splice construction; and `qm-levels.md` step 3's DFT
re-optimisation of each stationary point, which is itself a free optimisation and
inherits the same defect. **That last is a departure from the documented sequence
and is stated as one.** Also retired, having already been withdrawn for other
reasons: ¶340's free reactant minimum and ¶342's reverse scan.

Unchanged: report a potential-energy barrier as one and **never set it beside a
free energy**. CphA's ~14 kcal/mol is a PMF, and Melia 2015 measures PE − FE at
6.3-13.2 kcal/mol in that same enzyme.

## What this does not say

It does not say the scan's points are relaxed - they are not, and `saddle1`
measured by how much. It says the unrelaxed part is far-field and common to them,
so it cancels. It says nothing about leg 2, which has no barrier to score:
`ringopen` climbs monotonically to 53.03 kcal/mol with no interior maximum.

## Supported by
- [[far-field-environment-energy-is-not-part-of-the-barrier]]
- [[our-stage-energies-are-not-on-a-common-surface]]

## Supports
- [[no-reaction-coordinate-has-ever-been-computed]]
