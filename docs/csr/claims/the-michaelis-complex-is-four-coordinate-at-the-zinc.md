---
type: claim
status: open
entries: [¶120, ¶196, ¶256, ¶278]
---

# The Michaelis complex has a four-coordinate zinc: two histidines, hydroxide and carbonyl

**This is the gate `site.yaml` applies as `Zn coordination number, range [4, 5],
required_at: [reactant]`, and the user asked where it came from, twice - the
second time asking whether it is four-coordinate in the crystals at all.**

**¶280 MEASURED IT, and the answer splits the gate in two
([[the-m1-site-is-four-coordinate-at-rest-and-not-when-substrate-bound]]).**
Across eleven crystallographic sites the coordination number runs **3 to 6**.
Four IS the resting M1 value - 6NC5 M1 (2 His + 2 waters), 1SML M1 and 5AXO M1
(3 His + water) all read exactly 4 - so the premise is a crystal value and not
invented. But **the substrate-bound sites are not four**: 5AXO`s ligand-bound M2
is **six**, and 7UHI`s M1, the one real Michaelis complex, is **three** at 2.6 Å
with the substrate oxygen at 2.53 and no water.

**So the window belongs at `start`, not at `reactant`.** Applied to a reactant it
would reject the experimental Michaelis geometry we are trying to reproduce. That
is a sharper objection than the one below and supersedes it as the reason this
node is open.

The argument for it: the resting site is four-coordinate - two histidines and two
waters (¶120). On binding, the substrate displaces one water and the other
becomes the nucleophilic hydroxide, so **2 His + hydroxide + carbonyl = four**.

**Why five is tolerated rather than rejected** (the user's decision, ¶196):
GOB-18's B3 site expands four- to five-coordinate on substrate binding, so a
retained second water is a real possibility rather than a build error; and a
second-shell Asp92 drifting inside the 2.5 Å measurement radius counts as five
here while not being five-coordinate in any chemical sense.

**Why it stays open, and this is the part the gate hid.** The count follows from
the carbonyl being on the metal, which is itself the hypothesis under test
([[carbonyl-coordinates-zn1-in-the-michaelis-complex]]), and two independent
experimental Michaelis complexes do not behave this way: **7UHI has no water
within 4.0 Å of either zinc**, and the L1 complex has **no metal-bound water at
all** ([[the-l1-michaelis-complex-has-no-metal-bound-water]]) - the substrate
displaces the nucleophile. If that happens here, the Michaelis complex is not
"four donors including a hydroxide" and the gate fails a correct structure.

**It is also not a gate MM can be asked to pass.** The count is measured at a
2.5 Å radius on a force field that keeps three waters resident at 100% (¶191)
and put the hydrolysed product's own O4x on the metal in 200 of 200 frames
(¶256) - so on that trajectory the gate fails for reasons that are about the
force field, not the mechanism. Related:
[[mm-cannot-be-gated-on-absolute-zn-donor-bond-lengths]].

## Supports
- [[every-acceptance-gate-encodes-a-mechanism]]

## Required by
- [[carbonyl-coordinates-zn1-in-the-michaelis-complex]]
