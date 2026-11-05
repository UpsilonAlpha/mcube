---
type: claim
status: supported
entries: [D22, D39, D40, D107, D131, D208]
checked: [D40, D107, D131]
---

# CSR-1 is mono-zinc: the Zn2 site is vacant in the wild type, with Lys226 in it

The Zn2 ligand set in wild-type CSR-1 is **D-Q-K** (Asp92, Gln93, Lys226)
against the canonical D-H-H, so only one of the three can ligate a metal
(D40). The positive evidence is the deposited construct itself: 6NC5's
`SEQADV` records mutate Gln93->His and Lys226->His, and only then does a metal
appear in the Zn2 pocket (ZN301, Asp91 OD2 1.95, His92 NE2 2.01, His225 NE2
2.04, HOH493 2.08; D39). **The only way anyone put a metal in Zn2 was to
engineer two of its ligands back.**

Vacancy in the wild type is still an inference from that, not a direct
observation - the wild-type files are apo and so cannot show an empty site
either (retractions.md section 3, on D131). The ITC affinity and the apo
wild-type crystals corroborate it from two further directions.

**Consequence for the mechanism:** there is no second metal to stabilise the
anionic leaving nitrogen, and no bridging position for the hydroxide. That is
the half of the canonical mechanism this enzyme does not have.

## Supports
- [[thesis-mechanism-differs-from-canonical-b3]]
- [[lys226-occupies-the-zn2-site]]

## Supported by
- [[zn2-ligand-set-restored-by-q93h-k226h]]
- [[wild-type-csr1-crystals-contain-no-metal]]
- [[csr1-zn2-metal-affinity-is-weak]]
