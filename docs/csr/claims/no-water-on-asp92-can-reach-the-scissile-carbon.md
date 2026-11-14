---
type: observation
status: supported
evidence: ours
entries: [¶284]
checked: [¶284]
---

# No water hydrogen-bonded to Asp92 can reach the scissile carbon, so the B2 route fails on geometry

The two conditions the user set for the CphA route (2026-11-14) are met in the
crystal-derived pose - the metal binds the substrate carboxylate and Asp92 is
nearby ([[the-cpha-conditions-are-met-in-the-crystal-derived-pose]]). The route
fails anyway, for a reason that has nothing to do with subclass:

    Asp92 OD1 - C1  **4.56 Å**      OD2 - C1  **5.22 Å**

20,000 probe points on each of four shells (2.7, 2.9, 3.1, 3.3 Å) around each
carboxylate oxygen, keeping only points at least 2.55 Å from every other heavy
atom of protein and ligand:

    shell   OD1 clash-free / nearest to C1      OD2 clash-free / nearest to C1
    2.7       3162   **3.90 Å**                   5147   **4.56 Å**
    2.9       2875   3.87                         4471   4.66
    3.1       2553   4.13                         3886   4.74
    3.3       2362   4.32                         3217   4.76

**Thousands of positions are free and not one of them is within 3.4 Å of the
scissile carbon.** The closest approach any Asp92-bound water can make is
**3.87 Å**. There is room for the water and room for the base, but not for both
plus the carbon - against the metal route's 2.63 Å
([[a-free-zinc-site-puts-a-hydroxide-263-angstroms-over-the-scissile-carbon]])
this is not a close call.

So [[cpha-route-is-closed]] is now closed on measured geometry rather than on a
family label, which is what the previous two closures lacked.

**What survives from CphA** is the metal's job description, one metal holding
the substrate's carboxylate and its nitrogen - and that is also L1's Zn2
([[canonical-b3-puts-the-carboxylate-and-the-leaving-nitrogen-on-the-same-metal]]),
so it is not a B2 feature at all.

## Supports
- [[cpha-route-is-closed]]
- [[the-metal-can-both-carry-the-nucleophile-and-stabilise-the-leaving-nitrogen]]

## Refutes
- [[asp92-od1-hydrogen-bonds-the-nucleophile-in-the-reactant]]
