---
type: observation
evidence: ours
status: supported
entries: [¶238, ¶240, ¶244]
checked: [¶240, ¶244]
---

# GFN2's own minimum from the crystal-consistent frame puts Asp92 further from Zn1 than the crystal does, and pulls the substrate carbonyl onto the metal instead

`aspxtb`, job 32461816: GFN2-xTB/MM optimisation, **converged in 87 steps** and
228 s, `opt_energy` −430.99370531860416 Eh, energy falling monotonically
(−430.7337 → −430.9937 Eh). Identical to [[dft-gradient-does-not-pull-asp92-onto-zn1]]
in every input but `qm_engine` — enforced by a `same_as` gate that refuses to
write the section otherwise, not by inspection.

| pair | start | GFN2 min | 6DQH | DFT, 3 steps |
| --- | --- | --- | --- | --- |
| Zn–OD2 | 5.153 | **5.383** | 5.14 | 5.087 |
| Zn–OD1 | 4.327 | **4.952** | 4.28 | 4.544 |
| Zn–CG | 5.192 | **5.564** | 5.15 | 5.198 |
| Zn–NZ | 5.276 | 5.195 | 5.27 | 5.269 |
| NZ–OD2 | 2.880 | **2.718** | 2.79 | 2.766 |
| Zn–OH | 2.035 | 1.775 | | |
| Zn–O2x | 2.610 | **2.282** | | |
| O01–C4x | 3.481 | 3.673 | | |

**This was pre-registered as the second of three outcomes (¶238), and it is the
one that happened.** The written reading was: both methods keeping Asp92 out
means ¶232's chelation is not a minimum GFN2 puts in the wrong place, and the
claim about over-coordination has to be narrowed to say so.

**Two differences from ¶232, and ONE OF THEM IS NOW CLOSED (¶244).** That was
50 ps at 300 K with the whole system free; this is 0 K **and** had the MM
environment frozen, because the active region was the QM atoms only. `aspxtb8`
(job 32489558) repeated it with **479 active atoms** instead of 124 — same
Hamiltonian, same QM region, nothing constrained, started from this minimum,
converged in 103 steps and 15.3 minutes. **Freeing the second shell does not
bring Asp92 in.**

| pair | this minimum | 8 Å active region | 6DQH |
| --- | --- | --- | --- |
| Zn–OD2 | 5.383 | **5.092** | 5.14 |
| Zn–OD1 | 4.952 | **5.830** | 4.28 |
| Zn–CG | 5.564 | **5.777** | 5.15 |
| Zn–NZ | 5.195 | **5.044** | 5.27 |
| NZ–OD2 | 2.718 | **2.749** | 2.79 |
| Zn–O2x | 2.282 | 2.211 | |

This was **pre-registered as the second of three outcomes in ¶241**, and it is
the one that happened: the frozen shell was not what held the carboxylate off
the metal, so of the two explanations for ¶232's chelate only **thermal
energy** is left. The carboxylate also *rotates* once the shell is free — OD2
in 0.29 Å to the crystal distance, OD1 out 0.88 Å to the largest Zn–OD1
distance anywhere in this project, 1.55 Å off 6DQH. Monodentate to Lys226 and
nothing to the zinc: a sharper version of the thesis's geometry than the
crystal shows, and the largest single disagreement with the crystal on this
arm. Both halves of that stand together.

**The end-of-run energy rise is not a failure.** Steps 102–103 climb
4–7e−06 Eh — 0.003–0.004 kcal/mol — with Quality −2.99 and −1.96 and the
trust radius collapsing to 2.294e−04. geomeTRIC rejected those steps, which is
the healthy response; the channel project's 8 Å failure was 30–45 kcal/mol
uphill *with* Quality ~1. The screening test needs the magnitude as well as
the quality.

**The carbonyl coordination is new and is not a resolution of anything.** Zn–O2x
2.282 Å against 7UHI/7UHJ's 2.53/2.63 is the first geometry in this project with
the carbonyl on the metal and nothing restraining it there — but the pose is
inherited from the same docking that
[[no-usable-michaelis-complex-exists-yet]] flags, the method is one this project
has been distrusting for two days, and one minimum is not a measurement.

## Supports
- [[gfn2-over-coordinates-zn1-from-a-crystal-consistent-start]]

## Requires
- [[relax10-holds-the-crystal-geometry-for-10-ns]]
- [[question-is-the-reactant-pose-experiment-anchored]]
