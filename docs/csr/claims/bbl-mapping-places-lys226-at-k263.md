---
type: observation
status: supported
evidence: literature
entries: [D22, D40, D107, D108, D179]
checked: [D40, D107, D179]
---

# The BBL mapping puts Lys226 at K263, a Zn2 metal-ligand position

The mapping the whole argument depends on, from Pedroso 2020 and Krco 2023
(D22, D40, re-settled in D108/D179):

    BBL 116 -> His88      BBL 120 -> Asp92      BBL 196 -> His164
    BBL 118 -> Arg90      BBL 121 -> Gln93      BBL 254/263 -> Lys226

So the Zn1 site is **H-R-H** (BBL 116/118/196) and the Zn2 site is **D-Q-K**
(BBL 120/121/263), against canonical HHH/DHH.

**Lys226 is the Zn2 metal ligand, not a carboxylate-binding lysine.** Confirmed
three ways: D22's own table maps 226 -> 254 -> 263; D40 calls the Zn2 set
"D-Q-K ... one usable ligand"; and the 6DR8 mutant R118H/Q121H/**K254H**
restores exactly that histidine. It is **not** BBL Lys224, the B1 carboxylate
anchor - B3 enzymes have Ser221/Ser223 there. D12/D13's reading of it as a
carboxylate lysine is retracted (D107), and the EcoRV Lys92 analogy is retired
with it.

**Numbering hygiene, which has cost time twice.** The Lys254/Lys263 discrepancy
is one residue under two schemes - depositor sequential vs BBL. 6NC5 is our
numbering **minus 1**; 6DN4 and 6DQ2 are in **BBL** numbering (D108, D179).

## Supports
- [[lys226-occupies-the-zn2-site]]
- [[csr1-zn1-site-is-h-r-h]]
