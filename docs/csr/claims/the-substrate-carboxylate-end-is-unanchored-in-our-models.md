---
type: observation
status: supported
evidence: ours
entries: [¶272, ¶273]
checked: [¶273]
---

# The substrate's far end has no metal anchor in CSR-1, so its contacts are run-dependent rather than fixed

In CSR-1 the only part of the substrate held by a metal is the **β-lactam
carbonyl**. In `aspxtb8`, the converged GFN2 QM/MM minimum:

    Zn - O2x   (beta-lactam carbonyl)   2.211    7UHI 2.53  -- HELD
    Zn - N1x   (lactam nitrogen)        4.544    7UHI 2.67, on Zn2
    Zn - O4x / O3x  (C3 carboxylate)    5.422 / 6.697    7UHI 2.65, on Zn2

Ligand topology derived rather than assumed: the ring is C3x-C4x-N1x-C17x with
**C4x=O2x at 1.21 Å**, amide intact, so O2x is the carbonyl oxygen.

**The difference from 7UHI is expected, not a defect** (user, 2026-11-14). L1
anchors the lactam nitrogen *and* the C3 carboxylate on **Zn2**
([[l1-michaelis-complex-uses-both-zincs]]); CSR-1 has no Zn2
([[csr1-is-mono-zinc-with-zn2-vacant]]), so there is nothing in this enzyme for
those atoms to bind. A mono-zinc site is *supposed* to differ from a di-zinc
template here, and grading it against one imports the canonical mechanism the
project is testing.

**What follows is that the far end's contacts vary between runs, and the record
shows both extremes for the same pair:**

    Ser186 OG - mer O4x     100%          791-frame QM/MM run (log line 515)
                            0.0%          md_wt100b, 100 ns: 4.78 / 6.91 / 8.50
                            6.464 Å       aspxtb8 -- on the 100 ns run's mean

**And the 100% figure is the less trustworthy of the two**, because that run's
site was in the over-coordinated GFN2 state the project has been fighting —
Asp92 on the metal at 2.11 Å, coordination number 5.9. So `aspxtb8` is
**typical of this system rather than degraded**, and a claim that it "lost" the
contact quotes one half of a contradictory record (corrected in ¶273; ¶272 had
it wrong).

Two other numbers confirm `aspxtb8` is unremarkable for this site: that run put
O01···O2x at **2.77 Å** (aspxtb8: 2.781) and found attack angles of **20–53°**
over 791 frames (aspxtb8: 38.5°).

**So this node makes no claim that the pose is wrong.** It records that the
substrate's far end is held by nothing in this enzyme, which is why its
contacts are not reproducible between runs, and why a single frame's value for
any of them carries little weight.

## Supported by
- [[l1-michaelis-complex-uses-both-zincs]]
