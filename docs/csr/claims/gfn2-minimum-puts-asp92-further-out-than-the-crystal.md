---
type: observation
evidence: ours
status: supported
entries: [D238, D240]
checked: [D240]
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

**This was pre-registered as the second of three outcomes (D238), and it is the
one that happened.** The written reading was: both methods keeping Asp92 out
means D232's chelation is not a minimum GFN2 puts in the wrong place, and the
claim about over-coordination has to be narrowed to say so.

**Two differences from D232, not one.** That was 50 ps at 300 K with the whole
system free; this is 0 K **and** has the MM environment frozen, because the
active region is the QM atoms only. So the chelate's absence is consistent with
needing thermal energy, with needing the second shell to relax, or with both,
and these runs cannot separate them. What they do establish: the chelate is not
where GFN2's gradient goes from a crystal-consistent start with this environment
held.

**The carbonyl coordination is new and is not a resolution of anything.** Zn–O2x
2.282 Å against 7UHI/7UHJ's 2.53/2.63 is the first geometry in this project with
the carbonyl on the metal and nothing restraining it there — but the pose is
inherited from the same docking that
[[no-usable-michaelis-complex-exists-yet]] flags, the method is one this project
has been distrusting for two days, and one minimum is not a measurement.

## Supports
- [[gfn2-over-coordinates-zn1-from-a-crystal-consistent-start]]
