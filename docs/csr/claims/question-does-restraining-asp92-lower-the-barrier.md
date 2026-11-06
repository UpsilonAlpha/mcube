---
type: question
status: open
entries: [D157, D213, D214]
---

# Does holding Asp92 off Zn1 actually lower the hydrolysis barrier, and by how much?

Open, and it is the question the thesis turns on. Everything structural is
measured in crystals; the claim that the restraint **allows activity** is causal
and nothing computed has tested it.

**Why no existing result answers it.** No reaction coordinate has ever been run
in this project ([[no-reaction-coordinate-has-ever-been-computed]]). The one
relevant profile that does exist, `scan_znasp` in `CSR_crystal_v2`, is along the
Zn-Asp92 distance rather than along the reaction, and it is on the lineage whose
pose was pulled onto the metal by hand (D68): 2.05 -> 3.05 A costs
**+8.46 kcal/mol, monotonically uphill, minimum at 2.175 A**. That says the
energy function prefers the carboxylate on the metal; it says nothing about
hydrolysis.

**What would answer it, in increasing order of cost.**

1. **Two barriers in one system at one level** - the hydrolysis step computed
   with Asp92 restrained off the metal, and again with it on. The difference is
   the answer, and the shared system cancels most of the systematic error.
2. **The double mutant as the experimental arm.** Q93H/K226H puts a real zinc
   where the lysine was and holds Asp92 at 1.95 A
   ([[engineered-zn2-holds-asp92-at-195-angstrom]]). Its barrier, computed on the
   same protocol, is the comparison the kinetics can be checked against - and
   unlike route 1 it does not require inventing a restrained reference state.
3. **A PMF along the Asp92 approach** in the equilibrated wild-type complex,
   which gives the free-energy cost of the collapse itself rather than its
   effect on the barrier. Scoped as `aspcg` -> `asp`, 26 QM/MM windows, ~8,500
   core-hours at the measured rate.

Route 2 is the user's stated priority (2026-11-06): both systems through the
same pipeline first, barriers after, and the PMFs as additional evidence.

## Required by
- [[asp92-on-zn1-interferes-with-hydrolysis]]
- [[thesis-lys226-restrains-asp92-in-place-of-zn2]]
