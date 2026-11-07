---
type: claim
status: supported
entries: [D9, D231, D232, D238, D240]
checked: [D232, D240]
---

# GFN2-xTB's 300 K dynamics bring Asp92 onto CSR-1's Zn1 as a bidentate chelate from a crystal-consistent start, while GFN2's own minimum from that start does not: the bias is in the barrier, not in the minimum

**NARROWED 2026-11-07 (D240), exactly as pre-registered.** D238 wrote down, before running the matched GFN2 optimisation, that if both methods kept Asp92 out then "D232's chelation needed 300 K, so GFN2 makes a barrier crossing too easy rather than putting the minimum in the wrong place. A different, weaker claim, and this node would have to be narrowed to say exactly that." That is what happened: GFN2's converged minimum puts Asp92 **further out than the crystal** at Zn-OD2 5.383 ([[gfn2-minimum-puts-asp92-further-out-than-the-crystal]]), so the heading no longer says "over-coordinates" without qualification.

**What is still true, and what is given up.** D232's measurement stands - 50 ps of unrestrained GFN2 at 300 K takes the carboxylate from 4.34 to 2.01 A within 8 ps and holds it bidentate for 24 ps more, from a start four crystals support ([[aspqm2-chunk1-asp92-chelates-zn1-under-gfn2]]). What is given up is the stronger reading: GFN2 does not place its *minimum* on the wrong coordination, so a GFN2 geometry optimisation of this site is not discredited by this node, and the def2-SVP leg is a check on GFN2's dynamics rather than a replacement for a GFN2 minimum.

**One confound that cannot be removed by what has run.** Both optimisations are 0 K **and** hold the MM environment frozen (active region = the QM atoms, D235), while D232 was free at 300 K. So the chelate may need thermal energy, or the second shell's relaxation, or both; "it needed 300 K" picks one of three. The claim is written to survive all three, because each of them puts the difference between GFN2 and the crystals in the dynamics rather than in the stationary point.

**This answers [[question-does-gfn2-over-coordinate-the-zinc]], open since D9,
in the affirmative — and the answer was written down before the result
existed.** The two readings of `aspqm2` chunk 1 were pre-registered in that
question node on 2026-11-06 under the D178 protocol: *Asp92 staying beyond
~4 Å* closes it in the negative, *collapsing inside ~2.5 Å* closes it in the
affirmative, anything between means 50 ps was too short. The measurement is
2.10 Å mean at 78% occupancy, bidentate, for 24 consecutive picoseconds, with
no recrossing ([[aspqm2-chunk1-asp92-chelates-zn1-under-gfn2]]).

**Why this test is admissible where every earlier closure was not.** D206
rule 3 excludes single points and excludes two methods agreeing about a
geometry one of them was handed. This is neither: it is dynamics, unrestrained,
begun from a geometry four deposited structures support. The earlier collapses
([[our-models-put-asp92-on-the-zinc]]) all began either with the carboxylate
already on the metal or with the lysine deprotonated, so none of them tested
the Hamiltonian — D211's 1.93 Å mean is the clearest case. This one starts with
Asp92 at 4.3 Å from the metal, 2.8 Å from an intact ammonium, and GFN2 moves it
anyway.

**What it costs the project.** Every GFN2 geometry in which Asp92 sits on the
zinc is now suspect *as a geometry*, not merely as an energy — that is the
pre-registered consequence and it is not being softened here. Concretely:

- The def2-SVP QM/MM optimisation stops being a refinement of the GFN2 answer
  and becomes **the only admissible instrument** for Asp92's position.
- The planned Zn–Asp92 PMF (plan item 11, ~8,500 core hours at GFN2) would
  measure this bias rather than the free energy it is meant to measure. Its
  instrument needs settling first — see
  [[question-is-asp92-in-equilibrium-between-lys226-and-zn1]].
- It does **not** overturn any crystallographic claim. A method artefact cannot
  refute [[the-lysine-holds-asp92-where-a-zn2-would]], which rests on four
  deposits; what it says is that our method cannot reproduce it.

**What is still open.** The magnitude, and whether def2-SVP does the same thing.
One calibration point exists and points the same way:
[[gfn2-overstates-the-metal-coordination-preference]] measured GFN2 exaggerating
a coordination preference by 5.9 kcal/mol (×1.68) against wb97x-v/def2-svp on
51-atom clusters. A factor is not an error bar, and a cluster is not this site.

**The other half of the evidence is a crystal, and it has to be, or
"over-coordinates" names nothing.** This is a comparative claim: GFN2's 2.10 A
is only wrong against a measured distance, and the measured distance is 6DQH's
**4.28 A with the lysine present at 1.10 A resolution**
([[zinc-bound-crystal-holds-asp92-off-the-metal]]). Our own trajectory supplies
the first number and the deposit supplies the second, so the claim rests on both
classes rather than on ours alone - which is the one configuration in which a
calculation is legitimate evidence, because the claim is *about the
calculation*.

## Supports
- [[no-usable-michaelis-complex-exists-yet]]

## Supersedes
- [[question-does-gfn2-over-coordinate-the-zinc]]

## Supported by
- [[zinc-bound-crystal-holds-asp92-off-the-metal]]
