---
type: observation
status: supported
evidence: ours
entries: [¶272]
checked: [¶272]
---

# In our converged QM/MM minimum, everything the missing Zn2 would have held has drifted, including the Ser186 contact

Measured in `docs/csr/structures/aspxtb8_nowat.pdb`, the converged GFN2 QM/MM
minimum, internal distances only:

    Zn - O2x   (beta-lactam carbonyl)   2.211    7UHI 2.53  -- HELD
    Zn - N1x   (lactam nitrogen)        4.544    7UHI 2.67, on Zn2
    Zn - O4x / O3x  (C3 carboxylate)    5.422 / 6.697
    Ser186 OG - O4x                     6.464    record: contact in 100% of frames

Ligand topology derived rather than assumed: the beta-lactam ring is
C3x-C4x-N1x-C17x with **C4x=O2x at 1.21 A** and C4x-N1x at 1.37 A, so O2x is the
carbonyl oxygen and the amide is intact.

**The carbonyl is the only part of the substrate that is held**, and it is held
well — more tightly than in the only experimental B3 Michaelis complex. Every
other contact at the far end of the molecule is lost, and the **Ser186 OG–O4x
contact the record describes as present in 100% of frames is absent by more than
3 A.**

**Why this is structural rather than a failed equilibration.** 7UHI anchors the
lactam nitrogen *and* the C3 carboxylate on **Zn2**
([[l1-michaelis-complex-uses-both-zincs]]), and CSR-1 has no Zn2
([[csr1-is-mono-zinc-with-zn2-vacant]]). So in this site that end of the
substrate has nothing to hold it and it swings away. This is the pose problem
([[no-usable-michaelis-complex-exists-yet]]) stated as a measurement on the
current structure instead of as a comparison to a template.

**Consequence for the plan:** re-equilibrating this lineage cannot fix it, since
the drift is what the force field correctly does to an unanchored fragment. What
would fix it is an anchor — which is what Arg90 is positioned to be
([[arg90-and-lys226-act-at-opposite-ends-of-the-site]]), and Arg90 is itself
2 A out of position in this frame
([[our-his-zn-his-angle-is-13-19-degrees-too-closed]]).

## Supports
- [[no-usable-michaelis-complex-exists-yet]]

## Supported by
- [[l1-michaelis-complex-uses-both-zincs]]
