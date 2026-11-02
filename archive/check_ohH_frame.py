#!/usr/bin/env python3
"""Compare the ohH frame against the lysH frame it was built from.

    python check_ohH_frame.py lysH.pdb ohH.pdb ff1.xml ff2.xml ...

Pass the lysH frame first, then the ohH frame, then every force-field file the
qmmm stage will use - for lysH, with `hydroxide.xml`; the ohH build is then
checked with `water_nuc.xml` substituted for it automatically.

Both systems are built in full, which is what the qmmm stage itself does, and
the point of comparison is that **the two differ by where one hydrogen sits and
in nothing else**: same atom count, same total charge, a neutral lysine against
a cationic one, and a neutral water against a hydroxide. If any of those is
wrong the comparison measures something other than protonation.
"""

from __future__ import annotations

import sys
from pathlib import Path


def residue_charge(system, residue, nb):
    import openmm
    total = 0.0
    for atom in residue.atoms():
        total += (nb.getParticleParameters(atom.index)[0]
                  .value_in_unit(openmm.unit.elementary_charge))
    return total


def build(pdb_path, forcefields):
    import openmm
    from openmm.app import ForceField, PDBFile

    pdb = PDBFile(str(pdb_path))
    field = ForceField(*[str(f) for f in forcefields])
    system = field.createSystem(pdb.topology, rigidWater=False)
    nb = next(f for f in system.getForces()
              if isinstance(f, openmm.NonbondedForce))
    charges = [nb.getParticleParameters(i)[0]
               .value_in_unit(openmm.unit.elementary_charge)
               for i in range(system.getNumParticles())]
    return pdb, system, nb, charges


def find(topology, names):
    return [r for r in topology.residues() if r.name in names]


def main(argv):
    if len(argv) < 3:
        raise SystemExit(__doc__)
    lysh_pdb, ohh_pdb = Path(argv[0]), Path(argv[1])
    lysh_ff = [Path(a) if Path(a).exists() else a for a in argv[2:]]
    ohh_ff = []
    swapped = False
    for item in lysh_ff:
        text = str(item)
        if text.endswith("hydroxide.xml"):
            ohh_ff.append(text.replace("hydroxide.xml", "water_nuc.xml"))
            swapped = True
        else:
            ohh_ff.append(text)
    if not swapped:
        raise SystemExit(
            "no hydroxide.xml among the force fields, so there is nothing to "
            "substitute. Pass the lysH force-field list.")

    for path in (lysh_pdb, ohh_pdb):
        if not path.is_file():
            raise SystemExit(f"no such file: {path}")

    print(f"lysH {lysh_pdb.name}")
    lysh, lysh_sys, lysh_nb, lysh_q = build(lysh_pdb, lysh_ff)
    print(f"  {lysh_sys.getNumParticles()} particles, "
          f"total charge {sum(lysh_q):+.3f}")

    print(f"ohH  {ohh_pdb.name}")
    ohh, ohh_sys, ohh_nb, ohh_q = build(ohh_pdb, ohh_ff)
    print(f"  {ohh_sys.getNumParticles()} particles, "
          f"total charge {sum(ohh_q):+.3f}")

    if ohh_sys.getNumParticles() != lysh_sys.getNumParticles():
        raise SystemExit(
            "the two systems have different atom counts. One hydrogen should "
            "have moved, not been added or lost.")
    if abs(sum(ohh_q) - sum(lysh_q)) > 1e-4:
        raise SystemExit(
            f"total charge changed, {sum(lysh_q):+.3f} -> {sum(ohh_q):+.3f}. "
            "The proton moves inside the system, so it must not.")

    # The nucleophile: hydroxide (-1) becomes water (0).
    oh = find(lysh.topology, {"OH"})
    nuw = find(ohh.topology, {"NUW"})
    if len(oh) != 1 or len(nuw) != 1:
        raise SystemExit(f"found {len(oh)} OH and {len(nuw)} NUW residues")
    oh_q = residue_charge(lysh_sys, oh[0], lysh_nb)
    nuw_q = residue_charge(ohh_sys, nuw[0], ohh_nb)
    print(f"  nucleophile: {oh_q:+.3f} (OH, {len(list(oh[0].atoms()))} atoms) "
          f"-> {nuw_q:+.3f} (NUW, {len(list(nuw[0].atoms()))} atoms)")
    if abs(oh_q + 1.0) > 0.05 or abs(nuw_q) > 1e-6:
        raise SystemExit("the nucleophile did not go from -1 to neutral")

    # The lysine named on the command line is found by its missing HZ3.
    lysh_lys = {(r.chain.id, r.id): r for r in lysh.topology.residues()
                if r.name == "LYS"}
    changed = []
    for key, residue in lysh_lys.items():
        names = {a.name for a in residue.atoms()}
        other = next((r for r in ohh.topology.residues()
                      if r.name == "LYS" and (r.chain.id, r.id) == key), None)
        if other is None:
            continue
        if names - {a.name for a in other.atoms()} == {"HZ3"}:
            changed.append((key, residue, other))
    if len(changed) != 1:
        raise SystemExit(
            f"expected exactly one lysine to lose HZ3; {len(changed)} did")
    key, before, after = changed[0]
    before_q = residue_charge(lysh_sys, before, lysh_nb)
    after_q = residue_charge(ohh_sys, after, ohh_nb)
    print(f"  Lys {key[0]}:{key[1]}: {before_q:+.3f} -> {after_q:+.3f} "
          f"({len(list(before.atoms()))} -> {len(list(after.atoms()))} atoms)")
    if abs(before_q - 1.0) > 0.05 or abs(after_q) > 0.05:
        raise SystemExit(
            "the lysine did not go from +1 to neutral. It keeps the name LYS "
            "and should have matched the LYN template by bond graph; if it "
            "matched something else the charges would be wrong.")

    print("OK - one hydrogen moved; atom count and total charge unchanged, "
          "hydroxide -> water, cationic lysine -> neutral")


if __name__ == "__main__":
    main(sys.argv[1:])
