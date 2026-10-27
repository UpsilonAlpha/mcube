#!/usr/bin/env python3
"""Does the ohH input parameterise? Run before submitting anything.

    python check_ohH_system.py ohH_start.pdb water_nuc.xml

Loads the structure with the same force field list preprocess will use and
reports what it finds about the nucleophile. A preprocess job that fails on a
template match fails minutes in, but it fails after queueing, and the whole
point of the ohH arm is that it differ from the lysH arm in one hydrogen and
nothing else - so the things worth checking are mechanical: that the residue is
bonded, that it is neutral, that its charges are TIP3P's, and that rigidWater
leaves it flexible (which is why the equilibration uses 2 fs).

This does not solvate. It checks the one residue that is new.
"""

from __future__ import annotations

import sys
from pathlib import Path

TIP3P = (-0.834, 0.417, 0.417)


def main(argv):
    if len(argv) != 2:
        raise SystemExit(__doc__)
    pdb_path, xml_path = (Path(a) for a in argv)
    for path in (pdb_path, xml_path):
        if not path.is_file():
            raise SystemExit(f"no such file: {path}")

    import openmm
    from openmm.app import ForceField, Modeller, PDBFile

    pdb = PDBFile(str(pdb_path))
    nuc = [r for r in pdb.topology.residues() if r.name == "NUW"]
    if len(nuc) != 1:
        raise SystemExit(f"expected exactly one NUW residue, found {len(nuc)}")
    nuc = nuc[0]
    names = sorted(a.name for a in nuc.atoms())
    print(f"NUW: chain {nuc.chain.id} resSeq {nuc.id}, atoms {' '.join(names)}")
    if names != ["H03", "H04", "O01"]:
        raise SystemExit(f"NUW should hold O01, H03 and H04; it holds {names}")

    indices = {a.index for a in nuc.atoms()}
    bonds = [b for b in pdb.topology.bonds()
             if b[0].index in indices and b[1].index in indices]
    print(f"  bonds within NUW: {len(bonds)}")
    if len(bonds) != 2:
        raise SystemExit(
            "NUW needs two bonds and has "
            f"{len(bonds)}. Bonds for a non-standard residue come from CONECT "
            "records; PDBFile infers a water's bonds from the names O/H1/H2, "
            "which these are not.")

    forcefields = ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml",
                   str(xml_path)]
    print(f"  force fields: {', '.join(forcefields)}")
    field = ForceField(*forcefields)

    # Parameterise the nucleophile ALONE. This runs on the preprocess *input*,
    # whose protein still has no hydrogens - Modeller adds them during
    # preprocess - so building the whole system here is impossible and would
    # fail on the first alanine. Deleting everything else also keeps the one
    # thing being checked from hiding behind an unrelated template error.
    model = Modeller(pdb.topology, pdb.positions)
    model.delete([r for r in model.topology.residues() if r.name != "NUW"])
    kept = list(model.topology.residues())
    if len(kept) != 1 or kept[0].name != "NUW":
        raise SystemExit(f"isolating NUW left {[r.name for r in kept]}")
    indices = {a.index for a in kept[0].atoms()}
    system = field.createSystem(model.topology, rigidWater=True)
    print(f"  NUW parameterised alone: {system.getNumParticles()} particles")
    if system.getNumParticles() != 3:
        raise SystemExit("the isolated nucleophile is not three particles")

    nb = next(f for f in system.getForces()
              if isinstance(f, openmm.NonbondedForce))
    charges = [nb.getParticleParameters(i)[0]
               .value_in_unit(openmm.unit.elementary_charge)
               for i in range(system.getNumParticles())]
    print(f"  charges: {[round(c, 3) for c in charges]} (sum "
          f"{sum(charges):+.3f})")
    if abs(sum(charges)) > 1e-6:
        raise SystemExit("NUW is not neutral, so the total charge moved")
    if sorted(round(c, 3) for c in charges) != sorted(TIP3P):
        raise SystemExit(f"NUW's charges are not TIP3P's {TIP3P}")

    constrained = set()
    for i in range(system.getNumConstraints()):
        a, b, _ = system.getConstraintParameters(i)
        constrained.update((a, b))
    overlap = indices & constrained
    print(f"  rigidWater constrains NUW: {bool(overlap)}")
    if overlap:
        raise SystemExit(
            "rigidWater constrained the nucleophile. That would freeze its "
            "internal geometry, and the proton could not move - check whether "
            "the residue has been renamed to something OpenMM reads as water.")

    angles = 0
    for force in system.getForces():
        if isinstance(force, openmm.HarmonicAngleForce):
            for i in range(force.getNumAngles()):
                p1, p2, p3, _, _ = force.getAngleParameters(i)
                if {p1, p2, p3} <= indices:
                    angles += 1
    print(f"  H-O-H angle terms: {angles}")
    if angles != 1:
        raise SystemExit(
            f"NUW has {angles} angle terms, expected 1. Without it the water "
            "is floppy; with it, 2 fs is still needed because the angle is not "
            "constrained.")

    print("OK - the nucleophile is a neutral, bonded, flexible TIP3P water")


if __name__ == "__main__":
    main(sys.argv[1:])
