---
type: claim
status: supported
entries: [¶125, ¶126, ¶196, ¶201, ¶278]
checked: [¶125, ¶196, ¶201, ¶278]
---

# Every acceptance gate encodes a mechanism, so every gate belongs on this map

**The user's principle, 2026-11-14: "each simulation is a hypothesis, and should
fit on the argmap".** This node is the parent of the gates that were previously
only thresholds in `site.yaml`, `site.py` and `state.md` - numbers that could
fail a structure without any node on this map stating what they assume.

The project already learned this once, inside the config and not here. ¶125 wrote
"the closest substrate oxygen must be the beta-lactam carbonyl" as a fact about
Michaelis complexes; it is a fact about the **carbonyl-on-metal mechanism**,
which is the thing under test, and ¶126 then found a candidate complex that
violates it while satisfying the attack geometry at 3.17 Å and 106°.
`site.py`'s own docstring draws the conclusion: *"A gate that rejects a
competing hypothesis by assuming the first one is not a gate, it is a way of
never testing anything."*

So `site.py check` takes `--hypothesis=<name>` and names the set in its output,
and `site.yaml` holds two named sets - `carbonyl_on_metal` and
`substituent_on_metal`. **What was missing is the other half: a node here for
each gate, with a status, so a gate can be argued with and refuted instead of
only passed or failed.** `docs/csr/gates.md` is the index from each criterion to
its node.

**Two corrections to my own account of the gates (¶278), made by reading the
files rather than remembering them.** There is no third criteria set: the
top-level `michaelis_criteria` key is gone from `site.yaml` and `site.py` falls
back to it only when no hypothesis is defined. And `convergence.py` carries no
science criteria at all - it reports CA RMSD by decile, RMSF correlations over
windows, and T/P/E drift, writes nothing, and is legitimate run-quality work.
The science criterion wearing a convergence badge is in `state.md`'s narrative
gate, not in the code.

The experimental case is sharper than our own: **two measured Michaelis
complexes fail our gates.** 7UHI has no water within 4.0 Å of either zinc, and
the L1 complex has no metal-bound water at all
([[the-l1-michaelis-complex-has-no-metal-bound-water]]) - so a gate requiring a
metal-bound hydroxide would reject both of the only structures that are
certainly Michaelis complexes of this enzyme class.

## Supports
- [[thesis-mechanism-differs-from-canonical-b3]]

## Supported by
- [[the-l1-michaelis-complex-has-no-metal-bound-water]]
