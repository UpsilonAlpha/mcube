---
type: claim
status: supported
entries: [¶277]
---

# dmCSR-1, the plain double mutant, is the mutant the comparison should use

The manuscript's design is three templates - AIM-1 (HHH/DHH, canonical), CSR-1
(HRH/DQK, ours) and **dmCSR-1 (CSR-1 Q121H,K263H)**, introduced as a mutant
"to have an HRH/DHH metal-binding motif, as a third template for mutagenesis",
where "the modification of the M2 site represents a halfway-point between the
active sites of AIM-1 and CSR-1" - then the aspartate reversions on all three.

**So the wild-type-versus-mutant axis the paper is built on is CSR-1 versus
dmCSR-1, with the M2 site as the variable.** That is the same question this
project has been asking, so the computed comparison answers the paper's question
rather than a neighbouring one.

Four reasons it is the right arm, in order of weight:

1. **It isolates the variable.** `dmCSR_2metal.pdb` is **two substitutions** from
   our own wild-type file and nothing else (Q121H, K254H, verified by alignment).
   **Ser214 is intact**, so the intersite serine is held constant and the M2 site
   is the only thing that moves.
2. **Its metals are the real ones.** The zinc pair in our directory is native to
   this crystal ([[our-wt-crystal-zincs-are-grafted-from-the-double-mutant]]), so
   the mutant arm needs no metal placement at all - the one thing mcube cannot do.
3. **The substitution is remote from the substrate**, so a Michaelis pose
   transfers between the arms instead of each needing its own.
4. **The phenotype is a gain, not a confounded loss**
   ([[restoring-the-zn2-site-drastically-improves-catalysis]]), so the prediction
   is falsifiable in the informative direction: if the mono-zinc wild type runs a
   different route and the M2-restored mutant runs the canonical bridging-
   hydroxide route, that is the paper's claim turned into a mechanism.

**Not S214D**, although it is the paper's other mutant and although we hold its
crystal: [[s214d-cannot-isolate-the-serine-role]].

**What this claim does not settle.** The wild-type arm still has no experimental
metal ([[wild-type-csr1-binds-no-zinc-in-any-crystal]]) and no usable Michaelis
complex ([[no-usable-michaelis-complex-exists-yet]]); choosing a good mutant does
not fix either. It does mean the mutant arm is unblocked first, which makes it a
control on the wild-type protocol rather than merely the cheaper half.

## Supports
- [[thesis-mechanism-differs-from-canonical-b3]]

**¶279 REVISES REASON 2 AND THE SOURCE FILE.** `dmCSR_2metal.pdb` is **not an
original** ([[only-two-of-our-in-house-files-are-originals]]), so the mutant arm
must come from the deposit **6NC5** instead. And 6NC5 is **di-zinc**, both sites
at full occupancy, bridged by a water at Zn-Zn 3.43 - so "its metals are the real
ones" survives only by changing source, and the comparison is **one zinc against
two**, not one against one.

Reasons 1, 3 and 4 are untouched: the sequence ladder is established by alignment,
which a superposition cannot alter
([[our-crystals-number-the-catalytic-lysine-254-not-263]]).

**A one-zinc mutant at M1 does not exist**
([[whenever-a-canonical-m2-exists-the-metal-takes-m2]]): 6NC5 fills both sites,
the one original one-metal crystal is S214D with its zinc at **M2**, and 6DQH's
single full-occupancy M1 zinc comes with R118H, which deletes the arginine the
novelty depends on
([[the-novelty-is-a-monometallic-m1-site-with-two-histidines-and-an-arginine]]).

## Requires
- [[zn2-ligand-set-restored-by-q93h-k226h]]
- [[only-two-of-our-in-house-files-are-originals]]
