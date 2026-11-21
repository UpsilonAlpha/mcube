---
type: claim
status: open
entries: [¶301, ¶304, ¶305]
---

# A mono-zinc B3 enzyme carries its single metal in the DHH (M2) site, not the 3H (M1) site

**Why this matters more than anything else on this map right now.** Every
structural argument in this project places our single metal so that the
**beta-lactam carbonyl** coordinates it, which is an **M1** arrangement -
[[carbonyl-coordinates-zn1-in-the-michaelis-complex]], now measured at 2.53-2.63 Å
in the L1 time-resolved series. If a mono-zinc B3 instead holds its metal at
**M2**, the metal's substrate contacts are the **C-3 carboxylate and the leaving
nitrogen**, and the carbonyl is 4 Å away - the opposite assignment.

**The evidence, and its weakness.** Mono-zinc **GOB-18** is reported to bind its
single metal in the DHH/M2 site and to be fully active against penicillins,
cephalosporins and carbapenems (Moran-Barrio et al., AAC 2016,
<https://pmc.ncbi.nlm.nih.gov/articles/PMC5038331>, doi 10.1128/AAC.01067-16,
with earlier JBC spectroscopy and mutagenesis). The mono-zinc **B2** enzymes
also carry their metal in the DHH site, and there it is crystallographic.

**But there is no crystal structure of any mono-zinc B3.** The deposited GOB-18
(5K0W) is **di**-zinc, so the mono-Zn/DHH assignment rests on spectroscopy and
mutagenesis rather than coordinates, and it reached this project through a
search summary rather than the primary text. That is why this is `open` and not
`supported`.

**Three things it would explain at once if true**, each currently unexplained:

1. **Asp92 at 7.2 Å from the nucleophile** instead of the experimental 2.8-3.9 Å
   ([[asp92-od1-hydrogen-bonds-the-nucleophile-in-the-reactant]]). At M2 the
   aspartate ligates the metal, as it does in every reference enzyme, and the
   nucleophile lands in its second shell by construction.
2. **The substrate carboxylate dangling at 5.4 Å with no partner**
   ([[the-substrate-carboxylate-end-is-unanchored-in-our-models]]). In every B3
   complex measured it is pinned at 2.3-2.7 Å - by Zn2.
3. **The substrate leaving in every unrestrained QM/MM MD we have run** (¶111,
   ¶119, `aspqm2`): a complex anchored only through the carbonyl, with the
   carboxylate end free, is one contact away from nothing.

**What would settle it:** the primary GOB-18 spectroscopy read in full; a
crystal structure of any mono-zinc B3; or, internally, the metal placed at M2
and the resulting complex tested against the same gates. Our own crystals cannot
help - the only originals are the S214D mutant and no wild-type CSR-1 structure
has an experimental catalytic metal (¶277, ¶279).

**Our metal is already nearer M2 than the arrangement assumes**
([[our-single-metal-sits-between-the-canonical-zn1-and-zn2-positions]]): 1.77 Å
from where 7UHI's Zn2 lands in our frame.

**¶304: CHECKED AGAINST CSR-1'S OWN PAPER, AND IT POINTS THE OTHER WAY - but not
as cleanly as it first appeared.** Read directly from Pedroso et al.
(PMC7381538): *"in CSR-1trunc,sm electron density indicates the presence of one
bound Zn2+ in the α-site, while in the double and triple mutants both the α- and
β-sites are occupied"*. So when a CSR-1 variant holds **one** zinc, it is in α -
the 3H-equivalent site - which is the opposite of the GOB-18/B2 pattern this node
describes, and would make the metal's substrate contact the **carbonyl** after
all.

**Three qualifications that stop this closing the question.**

1. **That single mutant is the one whose α site was mutated back to canonical
   HHH** (the mutants were made *"to create variants where either the α, the β or
   both metal binding sites are identical to the canonical B3 MBL motif"*). In
   wild-type CSR-1, α is **HRH** - Arg90 in place of the third histidine - and
   *"no metal ions were observed in the crystal structures of CSR-1 and
   CSR-1trunc"*. So one-zinc-in-α is observed for an α-RESTORED variant, not for
   the wild type this project models.
2. **The ITC affinity split that was offered as the decisive evidence is
   unverified.** The figures reported to this project (α ~2 µM, β ~150 µM) are not
   in the paper's accessible text, which cites only *"Fig. S5, Table S4"*. Not
   quotable until the supplement is read.
3. **Our own one-metal crystal was read as disagreeing. ¶305 shows it does not.**
   `dmCSR_mutant_1metal_refine_12.pdb` does place its single zinc (occ **0.70**)
   at **β** - coordinated by Asp78 2.13, His79 2.22, His212 2.12 - but residues
   79 and 212 are **HIS in that file, not Gln and Lys**: `dmCSR-1` already carries
   Q121H and K263H, so its β site is canonical DHH **by engineering**. Their
   single mutant restored α and a metal appeared in α; ours restored β and the
   metal appears in β. **Both observations say only that a metal occupies
   whichever site has been repaired.**

4. **What our occupancies do show.** In `refine_16`, with both sites present in
   one crystal, **β is 1.00 and α is 0.53** - and that α site still has its
   arginine, coordinating through His74 and His150 alone. So an **HRH α site is
   half-occupied even against a canonical β**: weak, but real.

So the node stays `open`, with its content now the weaker half of a live
two-sided question rather than a standing rival. **The ligand-count argument
favours α** - α has lost one of three ligands, β two of three - and that is
reasoning, not measurement. **Settling this needs the Pedroso supplement
(Table S4), not compute**, and it gates the pose: α keeps the carbonyl contact
and `aspxtb8`, β makes the contact the carboxylate and leaving nitrogen and
discards it.

## Refutes
- [[carbonyl-coordinates-zn1-in-the-michaelis-complex]]

## Requires
- [[csr1-is-mono-zinc-with-zn2-vacant]]

## Supports
- [[question-how-far-is-asp92-from-zn1-in-wild-type-csr1]]
