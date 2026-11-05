---
type: question
status: open
entries: [D107, D191, D204, D207]
---

# What protonates the leaving nitrogen, with no Zn2 to stabilise it?

**The half of the canonical mechanism CSR-1 most clearly lacks.** In di-zinc
B3 enzymes Zn2 stabilises the anionic nitrogen after ring opening; CSR-1's Zn2
is vacant. That intermediate is long-lived and its decay is rate-determining
(Lisa 2017; Das & Nair 2017), so whatever substitutes for Zn2 here is
load-bearing for the rate, not a detail.

**What is already excluded.** Lys226 cannot be the donor to N1x in any model
measured: NZ-N1x 4.71-6.92 A, mean 5.76, over all 50 frames of the `nacqm`
QM/MM trajectory, and it never comes near the leaving nitrogen (D107). The one
MBL Lys proton-transfer study (Das & Nair 2017) is cephalexin only, runs through
a two-water chain that **meropenem's R2 group blocks**, and found no direct
Lys-to-N route even under a 30 kcal/mol bias.

**And nothing being positioned to protonate N1x in the Michaelis complex is the
expected result, not a defect** (D107): three papers converge on a water that
arrives *after* ring opening, and on the anionic intermediate lasting long
enough for it to.

**What would settle it** (`state.md` section 5): step 3 computed as its own
scan/NEB, with the proton source treated as a hypothesis set rather than
chosen in advance. Blocked on plan item 6, which is blocked on the endpoint
([[no-usable-michaelis-complex-exists-yet]]).

**Why it is a question and not a claim.** The RQK hypothesis predicts a
*larger* effect here than on nucleophilic attack
([[rqk-cations-replace-zn2-electrostatically]]), which makes this the step
where the thesis is most testable - and the step with no computation behind it
at all.

## Supports
- [[thesis-mechanism-differs-from-canonical-b3]]
