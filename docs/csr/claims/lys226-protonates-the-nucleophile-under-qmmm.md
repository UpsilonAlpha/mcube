---
type: observation
status: supported
evidence: ours
entries: [D29, D48, D202, D204]
checked: [D29, D48, D204]
---

# Lys226 hands a proton to the zinc-bound oxygen within 2 ps of becoming quantum

`michrk` (D202: Arg90 and Lys226 added to the QM region, charge 2, started from
batch 1's last frame). Measured 6 ps in (D204):

- **Lys226 hands HZ3 to the hydroxide within 2 ps**, and from 5 ps on the
  proton stays on the oxygen.
- It is a **shuttle**, not a one-way event: at 4 ps the proton is back on the
  nitrogen at 1.12 A - consistent with D48's -5.23 kcal/mol at 300 K, and
  correcting D29's one-way reading.
- **Zn-O01 reports the protonation state independently**, 1.87 A with the
  proton on N and 2.21-2.26 A with it on O, frame by frame, which no single
  wrong parameter would produce.

Also measured and new: **Ser186 OG** donates to the nucleophile at 1.60 A in
the first frame and leaves (4.06 A by the last); **Asn225 HD21** donates at
1.77 A in the last frame. Something hydrogen-bonds the nucleophile throughout
and it is **never Asp92** (OD1 5.67 then 5.32 A).

**Evidence class `ours`, with two live caveats.** The dynamics are GFN2 only,
so the magnitude in this geometry is open even if the sign is not. And
**Zn-Asp92 OD2 reaches 1.78 A** in this trajectory, which is the
over-coordination D48 showed inflates this very energy and which the structural
criterion rejects outright (D206). D204's claim that the coordination sphere
corroborates the protonation state "independently" is overturned: both
observables come from the same trajectory on the same geometry
(`retractions.md` section 3).

**Conditional on the starting geometry, measured 2026-11-07 (D232).** The same
QM region and the same charge, run for 32 ps from a start where Asp92 is 4.3 A
from the metal and the salt bridge is formed, shows **no transfer at all**: all
three Lys226 NZ hydrogens stay on the nitrogen at 0.98-1.09 A for every frame
and O01 keeps exactly one hydrogen
([[aspqm2-chunk1-asp92-chelates-zn1-under-gfn2]]). 32 ps against the 2 ps this
entry reports is not undersampling. So what this node records is real but
**not general**: the transfer happens in the lineage that began with Asp92
already on the metal, and does not happen in the lineage that began with the
carboxylate on the lysine. Whether the proton moves may depend on the
carboxylate's position, and that position is itself now a method artefact
([[gfn2-over-coordinates-zn1-from-a-crystal-consistent-start]]).

## Supports
- [[nucleophile-is-zn-bound-water-with-neutral-lys226]]
