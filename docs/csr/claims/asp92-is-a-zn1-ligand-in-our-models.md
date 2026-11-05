---
type: claim
status: supported
entries: [D51, D57, D58, D102, D133, D134, D138, D157, D201, D204]
---

# Asp92 coordinates the zinc in every model this project has built

Standing rule 3 records it as holding "at every level of theory we have run":
MM at 100 ns (D102), unconstrained GFN2-xTB at 25 ps (D51), DFT single points
(D58/D60), and QM/MM with no MM restraints at all (D201).

**It is carried here as a `claim`, not an observation, and its status is a
deliberate flag.** The measurement is real and lives at
[[our-models-put-asp92-on-the-zinc]] with evidence class `ours`. What this node
holds is the *conclusion* drawn from it - that the coordination is a real
feature of the enzyme - and almost every entry behind that conclusion is
retracted, reclassed or inadmissible:

- D58's inference is retracted (D59 at the time, D206 in place).
- D102, D134 and D138 are "inadmissible, not merely caveated" (D144).
- D157 section 3 cites those three as corroboration, so it cannot stand either.
- D201 section 4's conclusion is overturned in `retractions.md` section 3.
- D57's curve has our own model on both axes.

**Against the crystals it is wrong by 2 A or more.** The cost is not small
(D208): every structure built to hold Asp92 at 2.80 A has that contact too
short, including D130's hand-written pair, the `finalsystem.pdb` lineage, both
100 ns arms and the `michrk` start frame. The restraint target for any rebuild
is "do not restrain it, and let the run say where it goes".

## Supported by
- [[our-models-put-asp92-on-the-zinc]]

## Refuted by
- [[asp92-is-never-a-zn1-ligand-in-b3-structures]]
- [[asp92-is-a-zn2-ligand-in-b3-structures]]
