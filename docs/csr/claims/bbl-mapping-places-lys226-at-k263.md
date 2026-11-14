---
type: observation
status: open
evidence: literature
entries: [¶22, ¶40, ¶107, ¶108, ¶179, ¶277]
checked: [¶40, ¶107, ¶179]
---

# The BBL mapping puts Lys226 at K263, a Zn2 metal-ligand position

**DEMOTED to `open`, ¶277, and the demotion is narrow.** What this node calls
"numbering hygiene" is wrong: our own wild-type crystal contains **LYS at 254 AND
LYS at 263 simultaneously**, so they are two different residues in one numbering,
not "one residue under two schemes". By alignment, the lysine the double mutant
converts to histidine - the catalytic one, our Lys226 - is **254**
([[our-crystals-number-the-catalytic-lysine-254-not-263]]).

**What survives untouched:** the site assignment. Zn1 is H-R-H and Zn2 is D-Q-K,
Lys226 is the Zn2 ligand position rather than a carboxylate-binding lysine, and
the 6DR8 mutant's K254H restores exactly that histidine - all of which this node
argues from independent evidence and none of which depends on the 254/263
identification.

**What is open:** whether BBL numbers that lysine 263. ¶269 already established
that no CSR deposit is in BBL numbering, so the question cannot be settled from
any file we hold - it needs the Garau 2004 alignment or Pedroso 2020. Until then
the title's "K263" is a literature quotation, not a measurement, and the
manuscript`s `K263H` does not name our crystals` lysine.


The mapping the whole argument depends on, from Pedroso 2020 and Krco 2023
(¶22, ¶40, re-settled in ¶108/¶179):

    BBL 116 -> His88      BBL 120 -> Asp92      BBL 196 -> His164
    BBL 118 -> Arg90      BBL 121 -> Gln93      BBL 254/263 -> Lys226

So the Zn1 site is **H-R-H** (BBL 116/118/196) and the Zn2 site is **D-Q-K**
(BBL 120/121/263), against canonical HHH/DHH.

**Lys226 is the Zn2 metal ligand, not a carboxylate-binding lysine.** Confirmed
three ways: ¶22's own table maps 226 -> 254 -> 263; ¶40 calls the Zn2 set
"D-Q-K ... one usable ligand"; and the 6DR8 mutant R118H/Q121H/**K254H**
restores exactly that histidine. It is **not** BBL Lys224, the B1 carboxylate
anchor - B3 enzymes have Ser221/Ser223 there. ¶12/¶13's reading of it as a
carboxylate lysine is retracted (¶107), and the EcoRV Lys92 analogy is retired
with it.

**Numbering hygiene, which has cost time twice.** The Lys254/Lys263 discrepancy
is one residue under two schemes - depositor sequential vs BBL. 6NC5 is our
numbering **minus 1**; 6DN4 and 6DQ2 are in **BBL** numbering (¶108, ¶179).

## Supports
- [[lys226-occupies-the-zn2-site]]
- [[csr1-zn1-site-is-h-r-h]]
