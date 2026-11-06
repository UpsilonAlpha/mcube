---
type: question
status: open
entries: [D7, D185, D191]
---

# What is Lys226's pKa in the wild-type site, and is it protonated at the pH the kinetics were measured at?

Open. It gates [[lys226-is-cationic-in-the-resting-enzyme]], which the thesis
requires.

**Why it is not obvious.** The residue sits in the Zn2 pocket - a site built to
hold a dication - alongside Asp92's carboxylate and within reach of the
nucleophile (NZ-O01 at 2.72 A mean in 100% of 5000 frames, D191). A buried
lysine next to a carboxylate is the classic case of a shifted pKa, and the shift
can go either way: the anion raises it, the metal lowers it.

**Three routes, cheapest first.**

1. **PROPKA or pdb2pqr** on both wild-type deposits. Both are pure Python and
   neither is in the ASH env - plan item 5, open since D7, flagged as the user's
   call in D185. An estimate, not an answer, but it is hours not days.
2. **Constant-pH MD.** Gives a titration curve in the modelled site. Not
   currently wired into mcube.
3. **Thermodynamic integration** of the deprotonation, which is the rigorous
   version and the expensive one.

**What makes it urgent rather than tidy:** D204 transferred this proton in QM/MM
and Asp92 collapsed onto the zinc afterwards. If the resting state is in fact
deprotonated, the collapse is not an artefact at all and the thesis inverts.

**Re-scoped 2026-11-06 (D216).** This was a premise of the thesis and of
[[lys226-is-cationic-in-the-resting-enzyme]]. It is neither any longer: the
user ruled the resting state settled by the crystal geometry
([[a-neutral-lysine-would-not-sit-at-28-angstroms-from-a-carboxylate]]), which
does not need a pKa. What the pKa is still needed for is the state **during
turnover** - whether this residue can be the proton source for the leaving
nitrogen, which is where D204's deprotonation actually bears. The question is
unchanged; only what depends on it has moved.

## Required by
- [[question-what-protonates-the-leaving-nitrogen]]
