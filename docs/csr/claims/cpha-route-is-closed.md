---
type: claim
status: open
entries: [¶49, ¶51, ¶201, ¶284]
---

# The CphA route - Asp92 as the general base - is ruled out

**¶284 CLOSES THIS ON GEOMETRY, which the two earlier closures lacked.**
Asp92 OD1 is **4.56 Å** from the scissile carbon and OD2 is 5.22, and a search
of four hydration shells around both oxygens finds 2,362-5,147 clash-free
positions per shell with **none within 3.4 Å of that carbon** - closest
approach **3.87 Å**
([[no-water-on-asp92-can-reach-the-scissile-carbon]]). The route failed before
on a subclass argument the user rejected (2026-11-14); it fails now because the
base cannot reach, which is a reason that does not depend on what family the
enzyme is in.

**It stays `open`, and the validator is right to insist.** The 3.87 Å figure is a
search over empty space in a pose that has never been solvated or equilibrated
([[no-crystal-derived-pose-has-ever-been-equilibrated]]), so it is evidence about
our model, not about the enzyme; and the standing refutations are experimental and
literature - [[asp120-mutants-lose-activity-in-b1-and-b3]] and
[[asp92-is-never-a-zn1-ligand-in-b3-structures]] both outrank it. The premise
[[asp92-is-a-zn1-ligand-in-our-models]] is also void, and ¶284 is why
([[neither-aspartate-coordinates-the-metal-in-our-own-crystals]]) - so this claim
no longer rests on it, and nothing in the current argument needs the route closed:
the metal route wins on its own numbers.

¶201 closed this route for the third time (after ¶49 and ¶51). Its stated
reason, quoted: Asp92-as-general-base "has now failed three times - ¶49, ¶51,
and batch 1 - always **because it is a zinc ligand**, and the third time from a
correct start".

The route matters because it is the only mono-zinc precedent with a published
mechanism and a barrier to compare against: B2 CphA, where a non-metal-bound
water is activated by Asp120 as general base, concerted, dG# ~ 14 kcal/mol by
PMF (Xu, Xie & Guo 2006; ¶40, ¶207). With Zn2 vacant, Asp92 would be free to be
that base.

**The premise this closure rests on is refuted.** Asp92 is a zinc ligand only
in our own models ([[asp92-is-a-zn1-ligand-in-our-models]],
[[our-models-put-asp92-on-the-zinc]]), which ¶206 rule 3 excludes as evidence
about structure, and in no B3 structure is the residue a Zn1 ligand
([[asp92-is-never-a-zn1-ligand-in-b3-structures]]). ¶201's own bolded sentence
"Asp92 on the metal is a real feature of this model" is listed as overturned in
`retractions.md` section 3 - and that sentence is what justified starting
`michrk` from a 2.00 A Asp92 frame **and** ruling out this route a third time.

**This is the linkage failure the map exists to make visible.** The route was
closed on a premise the project had already measured its way out of; the index
records it, and nothing structurally connected the retraction to the closure.

A separate structural objection does survive, and is **not** touched by the
above: CphA is a Zn2-site enzyme and CSR-1 is not
([[mono-zinc-mbl-precedents-are-zn2-enzymes]]), so the analogy does not
transfer even if Asp92 is free (¶124). **That is now the only live reason to
set the route aside, and it is a different reason from the one on record.**

## Requires
- [[asp92-is-a-zn1-ligand-in-our-models]]

## Refuted by
- [[asp92-is-never-a-zn1-ligand-in-b3-structures]]

## Refuted by
- [[asp92-is-never-a-zn1-ligand-in-b3-structures]]
- [[asp120-mutants-lose-activity-in-b1-and-b3]]

**Demoted from `supported` to `open`, 2026-11-06.** ¶201 closed this route on the
premise that "Asp92 is a zinc ligand in every state this model reaches". That
premise has since been measured false in every deposited B3 structure
([[asp92-is-never-a-zn1-ligand-in-b3-structures]]) and in CSR-1's own wild-type
geometry, where the carboxylate sits in the second shell at 2.7-2.8 A held by
Lys226. The CphA precedent - a mono-zinc enzyme whose aspartate is not a Zn1
ligand - therefore applies better than ¶201 allowed, not worse.

Nothing here is asserted in its place: the route is **open**, not reopened. What
closes or confirms it is a barrier, and the aspartate's catalytic role in the
literature is read differently again
([[asp120-mutants-lose-activity-in-b1-and-b3]]).
