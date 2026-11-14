# Every gate, and the argmap node that argues it (¶278)

A gate is a number that can reject a structure or stop a chain. The user's
principle, 2026-11-14: **"each simulation is a hypothesis, and should fit on the
argmap"** — so no gate is allowed to fail a structure without a node here
stating what it assumes and carrying a status that can be argued with.

The parent node is [`every-acceptance-gate-encodes-a-mechanism`](claims/every-acceptance-gate-encodes-a-mechanism.md).
`claims/check.py` validates the graph; this file is the index from the config to
it. When a criterion is added to `site.yaml`, add its row here.

## `site.py check --hypothesis=carbonyl_on_metal`

| criterion | window | `required_at` | node | status |
| --- | --- | --- | --- | --- |
| closest substrate O | must be `O2x` | start, reactant | `carbonyl-coordinates-zn1-in-the-michaelis-complex` | the hypothesis itself |
| Zn – carbonyl O2x | 2.4–3.0 Å | start, reactant | same | the hypothesis itself |
| Zn – nucleophile O01 | 1.9–2.2 Å | reactant | `mm-cannot-be-gated-on-absolute-zn-donor-bond-lengths` | **supported** — crystal values, not MM ones |
| Zn coordination number | 4–5 | reactant | `the-michaelis-complex-is-four-coordinate-at-the-zinc` | **open** |
| O01 – C4x (attack) | 3.0–3.5 Å | reactant | `question-can-burgi-dunitz-coexist-with-carbonyl-on-one-zinc` | open |
| O01 – C4x – O2x | 95–120° | *transition_state* | `burgi-dunitz-gate-at-reactant-was-wrong` | open; prints, cannot fail (¶201) |
| O01 – Zn – O2x | 45–95° | reactant | `the-nucleophile-and-carbonyl-must-sit-cis-at-45-to-95-degrees` | **UNSUPPORTED** |
| O01 … O2x separation | 2.0–3.3 Å | reactant | same | **UNSUPPORTED** |
| Asp92 OD1 – O01 | 2.5–3.2 Å | reactant | `asp92-od1-hydrogen-bonds-the-nucleophile-in-the-reactant` | **open** |

## `site.py check --hypothesis=substituent_on_metal`

Every criterion in this set belongs to one node,
`question-does-the-substituent-on-metal-route-survive` (**open**). It had no node
on the map at all until ¶278 — a mechanism the pipeline could accept a structure
under, with nothing here to argue about it. Its `Zn – O1x 2.0–2.8 Å` window is
calibrated on ¶126's own 2.66 Å, which `site.yaml` already marks as "NOT a
reference value".

## Restraint targets, which are not `site.py check` criteria

| target | value | status |
| --- | --- | --- |
| `Zn – Asp92 OD2` | 2.80 Å | **VOID** (¶216) — a superposition into an apo site |
| `Zn – Asp92 OD1` | 3.81 Å | **VOID** (¶216) — imposed by ¶68's hand-pull |
| `Lys226 NZ – Asp92 OD2` | 2.79 Å | crystal, but phosphate-confounded (¶276); 6DN4 is the only phosphate-free deposit and gives 2.99 |

The live question is `question-how-far-is-asp92-from-zn1-in-wild-type-csr1`.
¶277 adds our own wild-type frame to it: **OD1 4.67 / OD2 5.14 / CG 5.36**,
agreeing with 6DQH to within 0.4 Å.

## Convergence

`convergence.py` is **legitimate and carries no science criteria** — CA RMSD by
decile, RMSF correlation between windows and between independent halves, T/P/E
mean, SD and drift. It writes nothing. Run-quality tests of this kind may gate a
run.

The one criterion that must not: *"unrestrained second-shell distances still
matching the crystal"*, which appears in `state.md`'s narrative gate for item 3.
That is a result, not a convergence test — see
`convergence-gates-must-not-encode-science-criteria`. Marked `unsupported` there
on purpose: it is a methodological norm, and no observation can hold up a norm.

## Two corrections to the gate list I gave on 2026-11-14 (¶278)

Both found by reading the files rather than remembering them:

1. **There is no third criteria set.** The top-level `michaelis_criteria` key is
   gone from `site.yaml`; `site.py` falls back to it only when no hypothesis is
   defined. The `reactant:` and `product:` keys I named are **QM-region
   anchors** (`&qmreg_wt`, `&qmreg_prod`), not criteria.
2. **`convergence.py` has no science criteria in it.** The science criterion
   wearing a convergence badge is in the prose of `state.md`, not in the code.

## The consequence not yet acted on

Two criteria are now marked **unsupported** and still `required_at: [reactant]`,
so they can still reject a structure: `O01 – Zn – O2x 45–95°` and
`O01 … O2x 2.0–3.3 Å`. ¶201's own precedent for a gate in this position was to
move it to a stage `site.py check` is never called with, so it prints and warns
but cannot fail. **Not done here** — changing what may reject a structure is a
scientific decision and needs the user's say.
