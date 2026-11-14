---
type: observation
status: supported
evidence: experimental
entries: [¶284]
checked: [¶284]
---

# Neither aspartate coordinates the M1 metal in either of our original crystals, and in the wild-type frame the Zn2 aspartate is CLOSER when the lysine is present

Our two PHENIX originals ([[only-two-of-our-in-house-files-are-originals]]),
numbered ours minus 14, counted at a 3.4 Å listing radius:

    dmCSR_mutant_2metal_refine_16.pdb
      M1 (ZN B2)  His150 NE2 2.05, His74 NE2 2.23, HOH 1.98 / 2.01 / 2.29
                  **Asp78 OD1 4.12**   (our Asp92, the Zn2 aspartate)
                  **Asp172 OD2 4.13**  (our Ser186, mutated to Asp)
                  Arg76 NH1 3.67
      M2 (ZN B1)  Asp78 OD2 2.06, His79 2.11, His212 2.18, HOH 2.03 / 2.43 / 2.51
                  Asp172 OD2 4.06

    dmCSR_mutant_1metal_refine_12.pdb
      single metal at **M2**: Asp78 OD2 2.13, His79 2.22, His212 2.12, HOH 2.07 / 2.16 / 2.38
                  **Asp172 OD2 4.23**, His74 NE2 4.03

**These are lysine-free structures** (K226H), and in them the Zn2 aspartate is
4.12 Å from the M1 metal. It does not coordinate it, approach it, or sit in its
second shell; where a canonical M2 exists it takes a metal of its own at 2.13 Å
([[whenever-a-canonical-m2-exists-the-metal-takes-m2]]).

**In the wild type it is closer, not further.** 6NC5 superposed onto metal-free
wild-type **6DN4** - the same protein, **core 257/258, RMSD 0.55 Å** - places
the M1 metal at:

    His192 NE2 1.87, His116 NE2 2.61, **Asp120 OD2 3.00**, Arg118 NH1 3.66,
    Asp120 OD1 4.12, Lys254 NZ 4.39, Ser214 OG 5.03

**3.00 Å with the lysine present, 4.12 Å with it gone.** So the proposition that
the aspartate swings onto the metal in the lysine's absence is refuted in its own
terms, and this number also answers
[[question-how-far-is-asp92-from-zn1-in-wild-type-csr1]].

**The confound, which no existing structure removes.** The second metal in
`refine_16` is itself engineered - `Q93H`/`K226H` restore the DHH set
([[zn2-ligand-set-restored-by-q93h-k226h]]) - so what holds Asp78 at 4.12 Å may
be that metal rather than the absence of the lysine. A lysine-free, single-M1-
metal CSR-1 crystal does not exist.

## Refutes
- [[asp92-is-a-zn1-ligand-in-our-models]]

## Supports
- [[asp92-is-never-a-zn1-ligand-in-b3-structures]]
- [[zinc-bound-crystal-holds-asp92-off-the-metal]]
