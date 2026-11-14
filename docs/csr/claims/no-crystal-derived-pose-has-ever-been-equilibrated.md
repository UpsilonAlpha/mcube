---
type: observation
status: supported
evidence: ours
entries: [¶196, ¶247, ¶282]
checked: [¶196, ¶282]
---

# No experimentally derived substrate pose has ever been solvated, equilibrated or run in this project

Checked against `site.yaml`, every launcher and every script:

- **`l1_in_relax10.pdb` is referenced nowhere as an input.** It contains **only
  the 23-atom MX0 ligand** - no protein, no metal. It is a measurement artefact
  from ¶247: never preprocessed, never solvated, never equilibrated, never the
  `pdb_file` of any section.
- **`CSR_mich`, the "Michaelis rebuild" of ¶196, is not it.** `build-mich.py`-s
  own instructions say to read the coordination report for *"Asp92 OD2 near 2.80,
  the carbonyl near 2.60"* - **both values voided by ¶216** - and its lineage is
  `CSR_crystal_v2`, the docked pose. The rebuild re-equilibrated the **docked**
  pose under the voided restraints.

**So every equilibrated structure in this project descends from a docked pose**,
and the one arrangement with experimental backing has never been given a solvent
box. That is why [[no-usable-michaelis-complex-exists-yet]] is still open after
this much compute, and it is a build, not a method problem: the transferred
ligand already exists at fit RMSD 0.310 Å, which was the hard part.

**It is also the experiment that would decide the mechanism**, because the two
candidate reactants support different chemistry
([[the-crystal-derived-pose-puts-the-carboxylate-on-the-metal]]) and only one of
them has ever been relaxed.

## Supports
- [[no-usable-michaelis-complex-exists-yet]]
- [[question-is-the-reactant-pose-experiment-anchored]]
