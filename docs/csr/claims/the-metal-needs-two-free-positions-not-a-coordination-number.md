---
type: claim
status: open
entries: [¶284]
---

# The requirement is two non-protein coordination positions, not a coordination number

The user's proposition was that the M1 metal must be four- or five-coordinate to
catalyse. Measured across every substrate-bound or resting metal available,
counted at 2.7 Å:

    7UHI  L1 Zn1           **4**   3 His + lactam carbonyl 2.53
    7UHI  L1 Zn2           **5**   Asp,His,His + carboxylate 2.65 + N1 2.67
    5B15  SMB-1 Zn1        **4**   3 His + bridging water 1.83
    5B15  SMB-1 Zn2        **6**   Asp,His,His + water + OAE 2.23 + NAN 2.33
    1X8I  CphA single Zn   **5**   Asp,Cys,His + N4 2.22 + O32 2.39
    5K0W  GOB-18 Zn1       **5**   Gln98 OE1 2.04, His, His, 2 waters
    5K0W  GOB-18 Zn2       **4**   Asp,His,His + water
    refine_16 M1           **5**   2 His + 3 waters
    refine_16 M2           **6**   Asp,His,His + 3 waters
    6NC5 M1 (ZN302)        **4**   2 His + 2 waters

**Substrate-bound catalytic metals run 4, 5, 5, 6.** A [4, 5] window rejects
5B15's Zn2, which is a real enzyme-substrate complex, so the count cannot be the
requirement. It fails because it mixes protein ligands with the two positions
that have to be free:

    canonical Zn1   3 His + substrate + nucleophile                       = 5
    CSR-1 M1        2 His (Arg90 replaces the third) + substrate + nucleophile = 4

**Same chemistry, different number.** Our crystal-derived arrangement (2 His +
carboxylate + hydroxide) and `aspxtb8` (2 His + carbonyl + hydroxide) both have
exactly two non-protein positions filled, and both read 4 while the canonical
case reads 5. The gate should count **open positions**, not donors.

**Open, because the right form of the count is not yet written down.** Counting
"non-protein donors == 2" would also pass a structure with the substrate on the
metal and a stray water instead of a hydroxide, and would fail 5B15's Zn2 again
(water + carboxylate + nitrogen = 3). So the criterion probably has to name the
two occupants rather than count them, which makes it a mechanism statement and
not a geometric gate - exactly the problem
[[every-acceptance-gate-encodes-a-mechanism]] records.

## Supports
- [[the-m1-site-is-four-coordinate-at-rest-and-not-when-substrate-bound]]
- [[every-acceptance-gate-encodes-a-mechanism]]

## Refutes
- [[the-michaelis-complex-is-four-coordinate-at-the-zinc]]
