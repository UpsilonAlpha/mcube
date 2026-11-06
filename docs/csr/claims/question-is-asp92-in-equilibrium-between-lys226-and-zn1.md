---
type: question
status: open
entries: [D216]
---

# Is Asp92 in equilibrium between Lys226 and Zn1, rather than committed to one of them?

The user's alternative to the thesis, and it is the sharper statement of it: the
residue is called a *bridging* aspartate because in canonical MBLs it is shared
between two metals, so in CSR-1 it might be shared between the lysine and Zn1
rather than held by either.

**What the crystals say.** They cannot answer it, and for a structural reason
rather than a resolution one: the comparison needs both partners present at once,
and only 6DQH has them. There the carboxylate is **2.79 A from the lysine and
4.28 A from the metal** ([[zinc-bound-crystal-holds-asp92-off-the-metal]]) - one
state, at 1.10 A, not a split population - and that is now checked on the
carboxylate's own order parameters rather than on the distance alone:
**B-factors below the protein median, occupancy 1.00, no alternate conformer**
([[asp92-is-rigid-in-the-crystal-and-the-lysine-is-the-mobile-partner]]). An
exchanging carboxylate would show two conformers or raised B. What the crystal
still cannot exclude is a minor state below the density threshold, which is
exactly a free-energy question. But 6DQH carries an engineered
histidine at the Arg90 position, and
[[wild-type-csr1-binds-no-zinc-in-any-crystal]] shows no wild-type deposit holds
a metal at all, so there is no crystallographic route to the wild-type answer.

**So this is a free-energy question, and it is the one the planned PMF measures.**
A Zn-Asp92 PMF along the distance coordinate returns exactly the shape that
answers it: one minimum near 4.3 A means the lysine is committed; two minima
within a few kT means an equilibrium; one minimum near 2.1 A means our models
were right and the crystal is the outlier. This promotes the PMF from additional
evidence to the test of a stated alternative - which is a better justification
for the ~8,500 core hours than the one the plan currently carries.

**It also changes what the PMF must start from.** A window set seeded from `eq6`
starts with the carboxylate pinned at 2.79 A from the metal, 1.5 A inside the
crystal distance, so the sampling would begin inside the basin the question is
about. The scan should span at least 2.0-5.5 A to cover both candidate minima.

## Requires
- [[zinc-bound-crystal-holds-asp92-off-the-metal]]
