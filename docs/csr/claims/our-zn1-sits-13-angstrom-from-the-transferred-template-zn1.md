---
type: observation
evidence: ours
status: supported
entries: [¶156, ¶159, ¶257]
checked: [¶159]
---

# Our Zn1 sits 1.31 Å from where 5AXO's Zn1 lands after a site fit, so any transfer-derived distance measured against OUR metal carries that offset

`refstruct.py transfer` of 5AXO (SMB-1 + hydrolysed meropenem) onto
`relax10`'s last frame, fitted on three site atoms — His72 NE2→His88 NE2,
His150 NE2→His164 NE2, Asp76 CG→Asp92 CG — **RMSD 0.232 Å**, the closest fit
this project has achieved (per-atom 0.218 / 0.228 / 0.250).

**The fit is good and the metal still lands 1.31 Å away.** Transferring 5AXO's
ZN301 and ZN302 through the same fit and measuring against `relax10`'s own
zinc:

| | value |
| --- | --- |
| our Zn1 → transferred template **Zn1** | **1.309** |
| our Zn1 → transferred template Zn2 | 2.831 |
| transferred Zn1–Zn2 separation | 3.544 |

And the displacement is visible as a chemically impossible contact: after the
fit, the **template's Zn1 sits 1.857 Å from our His88 NE2** (2.187 from our
His164), where our own zinc sits at 2.178 and 2.126. 1.857 Å is not a Zn–N
distance.

**Why this matters, and it is a method result rather than a structural one.**
The fit aligns the *site atoms*; **the metal is not a fit atom.** So a distance
from the transferred ligand to OUR metal is the crystal distance plus the metal
offset, and the offset here is larger than the fit residual by a factor of five.
Concretely: 5AXO's own Zn1–O32 is **3.41 Å** (¶159, measured in 5AXO), and the
transferred ligand's O32 reads **2.13 Å** against our zinc — the difference is
the 1.31 Å displacement, not a prediction.

**The reading this retires before it was recorded anywhere.** During the
2026-11-11 product-arm review it looked as though two independent routes agreed
that the hydrolysed product's C3 carboxylate binds Zn1 at ~2.1 Å — the MM
trajectory (2.06) and this transfer (2.13). **They do not agree; they share an
artefact.** The MM number is MM's, and the transfer number is our metal's
position. Nothing crystallographic says that carboxylate is 2.1 Å from CSR-1's
zinc. Caught by the check ¶156 prescribes for exactly this — *does the
template's own distance come back?* — and not by inspection of the fit RMSD,
which looked excellent throughout.

**One positive result from the same transfer.** `relax10`'s Asp92 OD2 is
**2.667 Å** from the transferred **Zn2** position, against Asp76 OD2–Zn2
**2.19 Å** in 5AXO: Asp92 sits near ligating distance from where a second metal
would be, in a frame where nothing placed it there
([[the-lysine-holds-asp92-where-a-zn2-would]]). And our zinc is 2.831 Å from
that Zn2 position, so it is not at the Zn2 site
([[csr1-zinc-sits-at-zn1]]).

**What this does not say.** It does not say our zinc is misplaced. CSR-1's Zn1
has **two** protein ligands where 5AXO has three — His74 ND1 has no counterpart,
Arg90 standing at that position instead ([[csr1-zn1-site-is-h-r-h]],
[[arg90-stands-at-the-third-zn1-ligand-position]]) — so a metal with one fewer
ligand relaxing to a different position is expected, and `relax10`'s 10 ns of MM
is a second possible cause. Which of the two it is has not been measured.

## Requires
- [[relax10-holds-the-crystal-geometry-for-10-ns]]

## Supports
- [[no-usable-michaelis-complex-exists-yet]]
