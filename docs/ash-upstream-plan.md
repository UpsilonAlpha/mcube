# Upstream: QM/MM point charges for xTB library mode, and g-xTB embedding

Two defects found while costing a 50 ps QM/MM equilibration. Both were submitted
upstream on 2026-10-23. This file is now the record, not the plan.

| | where | state |
| --- | --- | --- |
| g-xTB accepts point charges and ignores them | RagnarB83/ash#529 | open, +10/-0 |
| xTB library mode cannot return a PC gradient | RagnarB83/ash#530 | open, +13/-16, closes #157 |
| xtb-python does not expose the PC gradient | grimme-lab/xtb-python#130 | issue, patch held back |

## What is actually broken

**1. `gxTBTheory` accepts point charges and silently ignores them.** Its `run()`
signature takes `current_MM_coords`, `MMcharges`, `mm_elems` and `PC` — the whole
QM/MM contract — but the body references none of them and branches only on
`Grad`. Two consequences, traced through `module_QMMM.py`: at `:1219-1227` the
caller unpacks three values from a `PC=True` call while g-xTB returns two, so a
QM/MM gradient dies on an obscure unpack error; at `:1233` the energy-only path
assigns a single value, which unpacks cleanly and returns a **gas-phase energy as
though it were embedded**. That second case is wrong physics with no diagnostic.

It cannot be implemented instead: the g-xTB project lists point charge embedding
(`$pcem`) under its limitations.

**2. `xTBTheory(runmode='library')` cannot return a point-charge gradient.** Note
this is narrower than it first looks, and an earlier version of this file got it
wrong. The *setter* has always been wired up, so a library-mode **energy** with
point charges is correctly embedded. Only the gradient read was dead: the
`if PC == True:` branch printed "pc grad is not yet implemented", said to use
`runmode='inputfile'`, and called `ashexit()` — so it **fails loudly**, it does
not silently drop the embedding. Everything after that, including
`res.get_PCgradient()`, was unreachable.

Both were latent rather than historical: no result of ours came through either
path, because ASH's default is `inputfile` and we have never run g-xTB.

## What the APIs support

Checked against the sources, not assumed:

| call | status |
| --- | --- |
| `xtb_setExternalCharges(env, calc, n, numbers[n], charges[n], positions[n][3])` | in upstream `xtb.h` since 2.0.0 |
| `xtb_releaseExternalCharges(env, calc)` | in upstream `xtb.h` since 2.0.0 |
| `xtb_getPCGradient(env, res, gradient[n][3])` | in upstream `xtb.h` since 2.0.0 |
| `Calculator.set_external_charges(...)` | **wrapped** in xtb-python |
| `Calculator.release_external_charges()` | **wrapped** in xtb-python |
| `Results.get_pc_gradient(...)` | **absent** |
| `xtb_getPCGradient` in xtb-python's *vendored* `include/xtb.h` | **absent** |

That last row is the one that reshaped the work. `ffibuilder.py` `cdef`s the
vendored header wholesale, so a symbol missing from it is missing from the
compiled `lib` even when the linked `libxtb` exports it. Hence
`lib.xtb_getPCGradient` raises `AttributeError` and there is **no** workaround
from the Python side — not via the CFFI handle either. An earlier version of this
file called the missing `Results` method "the only missing piece"; it was not.

## Why xtb-python got an issue and not a patch

A 22-line patch exists locally (`Results.get_pc_gradient(npc)` plus the header
declaration) but was deliberately not submitted, because three decisions belong
to the maintainer and the patch is worthless if they go the other way:

1. **Header.** Add one declaration, or refresh the vendored header wholesale?
   `meson.build:36` already carries `# TODO: replace with include/_xtb.h`. Adding
   one declaration has no correct version suffix macro available — the vendored
   header defines only `XTB_API_SUFFIX__VERSION_6_3_0` and pins
   `XTB_API_VERSION 10000`.
2. **Minimum xtb.** `xtb_getPCGradient` arrived in API 2.0.0 (xtb 6.6.0), but
   `meson.build:32` falls back to `xtb-6.5.1`, which predates it, so a
   declaration would likely fail at link time without bumping the wrap.
3. **Signature.** `Results` cannot know the point-charge count and the C API has
   no `xtb_getNumPC`, so the count is either an explicit argument or cached at
   `set_external_charges` time. Caching is a trap: `singlepoint()`'s `copy=False`
   path reuses the `Results` and bypasses `__init__`, so a stale count would
   size the output buffer wrongly — an out-of-bounds write by the library, with
   no Python-level error.

## What was not verified

Nothing in either PR was executed against xtb; there is no libxtb on the machine
they were written on. Both are syntax-checked only, so no number in them is
verified. Specifically unverified: the shape, sign convention and units of the
returned PC gradient (`interface_DFTB.py:404-406` negates for DFTB+, so a sign
flip is a live risk), and whether releasing external charges after `singlepoint()`
disturbs the already-computed results object.

No test accompanies #530, for a concrete reason rather than squeamishness. ASH's
`environment.yml` installs `xtb-python`, so `pytest.importorskip("xtb")` would
**not** skip in CI; the test would reach the new code, find no wrapper, hit
`ashexit()` — which is `sys.exit`, a hard pytest failure rather than a skip — and
turn CI red on merge, staying red until xtb-python ships. The right check is a
finite-difference test modelled on `ash/tests/test_pyscf_pointcharge_gradient.py`,
which needs no reference numbers; it belongs in the PR conversation until the
wrapper lands. Do not copy that file's `1e-6` tolerance across without measuring
it — it is an SCF number at `conv_tol=1e-12`, not a GFN2 one.

## What this does not fix for us

mcube's own QM/MM stays on `inputfile` regardless, so the 50 ps equilibration
costs what it costs. Library mode only becomes an option once xtb-python ships a
wrapper *and* #530 is merged. This work was worth doing because g-xTB's
silent-vacuum energy is a trap for anyone using ASH, not because it unblocks the
CSR run.
