---
type: observation
evidence: experimental
status: supported
entries: [D216]
checked: [D216]
---

# In 6DQH - the only deposit holding both a Zn1 and the wild-type lysine - Asp92 is 4.28 A from the metal and 2.79 A from Lys226, at 1.10 A resolution

**6DQH had never been cited in this project.** An RCSB full-text search for this
enzyme returns exactly five deposits, and the record used four of them. 6DQH is
the fifth, it is **1.10 A, R = 0.179** - the best data in the system by a wide
margin (6DN4 1.99, 6NC5 2.07) - and it is the single structure that holds a zinc
at Zn1 while keeping the wild-type lysine and glutamine.

Deposit numbering, with the project's in brackets: Asp120 [92], Lys260 [226],
His116/His192 [88/164], Gln121 [93], position 118 [90].

    ZN401   HOH526.O     2.03   <- the terminal water, where the nucleophile goes
            His192.NE2   2.05
            His118.ND1   2.07   <- ENGINEERED (wild type is Arg at this position)
            His116.NE2   2.08
            PO4405.O2    2.32
            HOH525.O     2.50

    Zn-Asp120 OD1   4.28        Zn-OD2   5.14        Zn-CG   5.15
    Lys260 NZ-OD2   2.79        NZ-OD1   3.94        NZ-Zn   5.27
    Asp120 OD1-HOH526  2.64     <- H-bond to the nucleophilic water
    Gln121 OE1-Zn   5.25                             OE1-HOH526  4.60

Three things follow, and the third was not expected.

1. **The lysine wins outright.** With a metal present and the lysine present, the
   carboxylate sits with the lysine, not the metal. There is no equilibrium at
   the crystal level: 2.79 against 4.28 is not a two-state population, it is one
   state. Whether the two compete on a free-energy scale is a separate question
   that no crystal can answer, because the comparison needs both partners in one
   structure and only this one has them.
2. **Every model this project has built puts Asp92 too close to the metal.**
   `eq6` holds Zn-OD2 at 2.79 and Zn-OD1 at 3.90; earlier work placed the
   carboxylate at 1.95-2.05. The crystal says 4.28/5.14. By the standing rule -
   restrain what the metal affects at its crystal distance (D174, D186) - these
   are the distances to restrain, and they are 1.5-2.4 A further out than
   anything run so far.
3. **Asp92 hydrogen-bonds the nucleophile.** OD1 is 2.64 A from the Zn1 terminal
   water. So the carboxylate is not merely kept out of the way; it is parked in
   the second shell touching the water that attacks. That is a positive role, and
   it is the role Asp120 is given in canonical B1 schemes.

**The caveat that keeps this an observation and not the whole argument:** 6DQH is
the R118H single mutant [R90H], and the engineered histidine is itself a Zn1
ligand at 2.07 A. Asp120 is two residues away from the substitution. So this
measures the carboxylate's position in a site whose third Zn1 ligand is not the
wild-type one. No deposit is wild type at all three positions *and* metal-bound,
and [[wild-type-csr1-binds-no-zinc-in-any-crystal]] says why none can be.

## Supports
- [[thesis-lys226-restrains-asp92-in-place-of-zn2]]
- [[lys226-salt-bridges-asp92-in-wild-type-crystals]]

## Refutes
- [[asp92-is-a-zn1-ligand-in-our-models]]
