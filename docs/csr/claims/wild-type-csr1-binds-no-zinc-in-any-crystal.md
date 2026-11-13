---
type: observation
evidence: experimental
status: supported
entries: [¶216]
checked: [¶216]
---

# No crystal of wild-type CSR-1 contains a metal, and the published zinc affinity is micromolar-to-weak

Measured directly from the deposits: **6DN4 (wild type, 1.99 A) and 6DQ2 (S60,
1.99 A) contain zero metal atoms** - 347 and 63 waters, a phosphate in 6DQ2, no
Zn, Cd, Ni, Mn, Fe, Co, Mg or Ca. Pedroso et al. say the same in words: *"no
metal ions were observed in the crystal structures of CSR-1 and CSR-1trunc"*,
with ITC giving two sites at **alpha_Kd ~2 uM and beta_Kd ~150 uM**.

Zinc appears only where a ligand has been engineered back:

| deposit | 118 [90] | 121 [93] | 254 [226] | metal | Lys NZ-Asp120 | res. |
| --- | --- | --- | --- | --- | --- | --- |
| 6DN4 | Arg | Gln | Lys | **none** | 2.99 | 1.99 |
| 6DQ2 | Arg | Gln | Lys | **none** | 2.83 | 1.99 |
| 6DQH | **His** | Gln | Lys | 1, at Zn1 | 2.79 | **1.10** |
| 6NC5 | Arg | **His** | **His** | 6 (Zn1 + engineered Zn2) | no Lys | 2.07 |
| 6DR8 | **His** | **His** | **His** | 2 (Zn1, Zn2 at occ 0.45) | no Lys | 1.48 |

Two consequences for the argument.

**The one site that is wild-type-like at Zn1 is 6NC5's ZN302** - His87/His163
[88/164] plus two waters, Arg90 unmutated at 4.02 A - which is what
[[csr1-zinc-sits-at-zn1]] rests on, correctly. But that deposit has lost the
lysine. And the one deposit that keeps the lysine with a zinc, 6DQH, has an
engineered histidine at the Arg90 position. **There is no wild-type,
metal-bound, lysine-bearing structure, and the ITC numbers say there cannot be
one:** at beta_Kd ~150 uM the Zn2 site does not fill, and at alpha_Kd ~2 uM Zn1
fills only when the third ligand is restored.

**So the resting state the whole project models - one zinc at Zn1, Arg90 intact,
Lys226 intact - has never been observed.** It is assembled from two deposits that
each supply half of it. That is defensible and it is what the modelling is for,
but it is an assembly and the paper has to say so.

## Supports
- [[csr1-is-mono-zinc-with-zn2-vacant]]

## Requires
- [[csr1-zinc-sits-at-zn1]]
