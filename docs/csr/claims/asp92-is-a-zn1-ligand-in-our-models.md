---
type: claim
status: void
checked: [¶51, ¶57, ¶58, ¶102, ¶133, ¶134, ¶138, ¶157, ¶201, ¶204]
entries: [¶51, ¶57, ¶58, ¶102, ¶133, ¶134, ¶138, ¶157, ¶201, ¶204]
---

# Asp92 coordinates the zinc in every model this project has built

Standing rule 3 records it as holding "at every level of theory we have run":
MM at 100 ns (¶102), unconstrained GFN2-xTB at 25 ps (¶51), DFT single points
(¶58/¶60), and QM/MM with no MM restraints at all (¶201).

**It is carried here as a `claim`, not an observation, and its status is a
deliberate flag.** The measurement is real and lives at
[[our-models-put-asp92-on-the-zinc]] with evidence class `ours`. What this node
holds is the *conclusion* drawn from it - that the coordination is a real
feature of the enzyme - and almost every entry behind that conclusion is
retracted, reclassed or inadmissible:

- ¶58's inference is retracted (¶59 at the time, ¶206 in place).
- ¶102, ¶134 and ¶138 are "inadmissible, not merely caveated" (¶144).
- ¶157 section 3 cites those three as corroboration, so it cannot stand either.
- ¶201 section 4's conclusion is overturned in `retractions.md` section 3.
- ¶57's curve has our own model on both axes.

**Against the crystals it is wrong by 2 A or more.** The cost is not small
(¶208): every structure built to hold Asp92 at 2.80 A has that contact too
short, including ¶130's hand-written pair, the `finalsystem.pdb` lineage, both
100 ns arms and the `michrk` start frame. The restraint target for any rebuild
is "do not restrain it, and let the run say where it goes".

## Supported by
- [[our-models-put-asp92-on-the-zinc]]

## Refuted by
- [[asp92-is-never-a-zn1-ligand-in-b3-structures]]
- [[asp92-is-a-zn2-ligand-in-b3-structures]]

## Superseded by
- [[asp92-is-never-a-zn1-ligand-in-b3-structures]]

**Voided 2026-11-06.** The status was `supported` while the body said the claim
is "wrong by 2 A or more" against the crystals - which is the exact move the
validator's check 2 exists to find, and it sat here through the map's whole
first week. The conclusion is retracted; the measurement it was drawn from is
untouched at [[our-models-put-asp92-on-the-zinc]].

Two things have changed since the node was written. The crystals now say
positively where the carboxylate sits in the wild type - the second shell, held
by Lys226 ([[lys226-salt-bridges-asp92-in-wild-type-crystals]]) - so there is a
replacement, not just a refutation. And the body's prescription, "do not restrain
it, and let the run say where it goes", is **superseded by the user's standing
rule** (¶174, ¶186, restated 2026-11-05): Asp92 is affected by the metal and is
to be restrained at its structural distance, the 3.0 A radius being "what is
affected by the metal" rather than a coordination test. `eq6` was built that way
(¶214).
