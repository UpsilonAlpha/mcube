#!/usr/bin/env python3
"""
Explicit lipid-bilayer support (OpenMM only, no ASH): force-field selection, the
membrane barostat, the system the md stage runs, and box/area-per-lipid reporting.

The model is the one OpenMM's `Modeller.addMembrane` implements: the membrane normal
is +z, the bilayer is centred on `membrane_center_z`, and the protein is assumed to be
oriented and positioned already. Nothing here (and nothing in ASH) rotates or shifts
the protein, so a crooked input gives a crooked bilayer. Align the transmembrane
bundle with z and centre its hydrophobic belt on `membrane_center_z` first (the OPM
database, http://opm.phar.umich.edu, publishes membrane-aligned structures).

Why the barostat lives here rather than coming from ASH: ASH only ever adds
`openmm.MonteCarloBarostat`, which scales x, y and z by the same factor. A bilayer
has to be able to change its area per lipid independently of its thickness, so the
md stage builds the OpenMM System itself and attaches a
`MonteCarloMembraneBarostat` (XYIsotropic, ZFree, zero surface tension) before
serialising it for `OpenMMTheory(xmlsystemfile=...)`.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence

import openmm
from openmm import unit
from openmm.app import ForceField, HBonds, PME
from openmm.app import forcefield as _openmm_forcefield

# Amber lipid force field (Lipid17); the one meant to be used with protein.ff14SB,
# which is mcube's default protein force field.
LIPID_XML = "amber14/lipid17.xml"

# Lipid types `Modeller.addMembrane` has pre-equilibrated patches for
LIPID_TYPES = ("POPC", "POPE", "DLPC", "DLPE", "DMPC", "DOPC", "DPPC")

# Residue names those patches use (Lipid17 keeps each lipid as one residue, truncated
# to three characters: POPC and POPE are both "POP")
LIPID_RESNAMES = frozenset({"POP", "DLP", "DMP", "DOP", "DPP"})

# ASH's own defaults for a periodic system, so a membrane system differs from mcube's
# explicit-solvent path only in the lipids and the barostat
NONBONDED_CUTOFF = 12.0       # Å
EWALD_ERROR_TOLERANCE = 5e-4
HYDROGEN_MASS = 1.5           # amu

PRESSURE = 1.0                # bar
SURFACE_TENSION = 0.0         # bar·nm (a tensionless bilayer)
BAROSTAT_FREQUENCY = 25       # steps

MEMBRANE_BAROSTAT = "MonteCarloMembraneBarostat"


# ==============================
# Force field
# ==============================

def check_lipid_forcefield(xml: str = LIPID_XML) -> str:
    """
    Check that `xml` can actually be found, and return it as given.

    mcube's `default_forcefields` (protein.ff14SB + tip3p) carry no lipid parameters,
    so `addMembrane` cannot type a lipid without this file. It ships with OpenMM, but
    it is checked against the installed OpenMM rather than assumed: a missing file
    otherwise surfaces as a residue-template error in the middle of the build.
    """
    if os.path.isfile(xml):
        return xml  # a lipid XML of the user's own
    data_dir = Path(_openmm_forcefield.__file__).parent / "data"
    if (data_dir / xml).is_file():
        return xml
    try:
        # Also picks up XML files installed by plugins, which are not under data/
        ForceField(xml)
    except Exception as exc:
        raise FileNotFoundError(
            f"The installed OpenMM ({openmm.version.version}) does not provide {xml}, which the "
            f"membrane path needs to parameterise lipids (looked in {data_dir}). Install a newer "
            "openmm from conda-forge, or point `membrane_lipid_xml` at a lipid force-field XML."
        ) from exc
    return xml


def lipid_forcefields(forcefields: Iterable[str], xml: str = LIPID_XML) -> List[str]:
    """`forcefields` with the lipid force field appended, unless it is already there."""
    forcefields = list(forcefields)
    if any(Path(f).name == Path(xml).name for f in forcefields):
        return forcefields
    return forcefields + [check_lipid_forcefield(xml)]


def check_lipid_type(lipidtype: str) -> str:
    """Validate a lipid type against the patches `Modeller.addMembrane` knows."""
    if lipidtype is None:
        return "POPC"
    if lipidtype.upper() not in LIPID_TYPES:
        raise ValueError(
            f"membrane_lipidtype {lipidtype!r} is not one of OpenMM's built-in membrane patches "
            f"({', '.join(LIPID_TYPES)})")
    return lipidtype.upper()


# ==============================
# Barostat
# ==============================

def membrane_barostat(temperature: float, pressure: float = PRESSURE,
                      surface_tension: float = SURFACE_TENSION,
                      frequency: int = BAROSTAT_FREQUENCY) -> openmm.MonteCarloMembraneBarostat:
    """
    A `MonteCarloMembraneBarostat` for a bilayer in the xy plane.

    XYIsotropic scales x and y together (the box stays square in the membrane plane)
    and ZFree lets the thickness move on its own, so the area per lipid and the
    bilayer thickness relax independently. `surface_tension` is 0 by default, i.e. a
    tensionless bilayer, which is what reproduces experimental areas per lipid.
    """
    if temperature is None:
        raise ValueError("A membrane system needs `temp` in the md inputs: the barostat is "
                         "constructed at the simulation temperature.")
    return openmm.MonteCarloMembraneBarostat(
        pressure * unit.bar,
        surface_tension * unit.bar * unit.nanometer,
        temperature * unit.kelvin,
        openmm.MonteCarloMembraneBarostat.XYIsotropic,
        openmm.MonteCarloMembraneBarostat.ZFree,
        int(frequency),
    )


def barostat_names(system: openmm.System) -> List[str]:
    """Class names of the barostat forces in an OpenMM system."""
    return [f.__class__.__name__ for f in system.getForces() if "Barostat" in f.__class__.__name__]


def check_single_barostat(system: openmm.System, expected: str = MEMBRANE_BAROSTAT) -> None:
    """
    Fail unless the system carries exactly one barostat, and that the expected one.

    ASH adds `openmm.MonteCarloBarostat` whenever `OpenMM_MD`/`OpenMM_box_equilibration`
    is given a `barostat`, and its "already present" check tests for that exact class
    name. `MonteCarloMembraneBarostat` does not match it, so ASH would happily stack an
    isotropic barostat on top of the membrane one — and isotropic scaling ties the box
    area to its thickness, which is exactly what a bilayer must not have. Hence the check,
    before and after every stage that could have added one.
    """
    found = barostat_names(system)
    if found != [expected]:
        raise RuntimeError(
            f"Expected exactly one {expected} in the OpenMM system, found "
            f"{', '.join(found) if found else 'none'}. Two barostats (or an isotropic one) would "
            "scale the bilayer area and thickness together.")
    print(f"[INFO] Barostat: one {found[0]} (XYIsotropic/ZFree)")


# ==============================
# System
# ==============================

def membrane_system(topology, forcefields: Sequence[str], temperature: float,
                    nonbonded_cutoff: float = NONBONDED_CUTOFF,
                    pressure: float = PRESSURE,
                    surface_tension: float = SURFACE_TENSION,
                    barostat_frequency: int = BAROSTAT_FREQUENCY,
                    hydrogen_mass: float = HYDROGEN_MASS) -> openmm.System:
    """
    The periodic System used for membrane MD: PME, HBond constraints, rigid water,
    heavy hydrogens, and one membrane barostat.

    Apart from the lipid force field and the barostat these are the same settings ASH
    would give the explicit-solvent path (PME, 12 Å cutoff, 5e-4 Ewald tolerance,
    HBonds, rigid water, 1.5 amu hydrogens), so the two modes stay comparable.
    """
    if topology.getPeriodicBoxVectors() is None:
        raise ValueError("The membrane PDB file has no periodic box (no CRYST1 record). "
                         "Run the preprocess stage with `membrane: true` to build one.")
    forcefield = ForceField(*forcefields)
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=PME,
        nonbondedCutoff=nonbonded_cutoff * unit.angstrom,
        ewaldErrorTolerance=EWALD_ERROR_TOLERANCE,
        constraints=HBonds,
        rigidWater=True,
        hydrogenMass=hydrogen_mass * unit.amu,
    )
    system.addForce(membrane_barostat(temperature, pressure=pressure,
                                      surface_tension=surface_tension,
                                      frequency=barostat_frequency))
    check_single_barostat(system)
    return system


# ==============================
# Reporting
# ==============================

def count_lipids(topology, resnames: Iterable[str] = LIPID_RESNAMES) -> int:
    """Number of lipid residues in a topology."""
    resnames = {r.upper() for r in resnames}
    return sum(1 for res in topology.residues() if res.name.upper() in resnames)


def area_per_lipid(lx: float, ly: float, n_lipids: int) -> float:
    """
    Area per lipid in Å², from the box side lengths in Å.

    The bilayer spans the periodic box in x and y, and each leaflet holds half the
    lipids, so the area per lipid is Lx·Ly / (N/2). POPC is roughly 60–68 Å² per lipid
    at 300 K.
    """
    if n_lipids < 2:
        raise ValueError(f"Need at least two lipids to get an area per lipid, found {n_lipids}")
    return lx * ly / (n_lipids / 2.0)


def box_table(lengths_angstrom: Sequence[Sequence[float]], n_lipids: int) -> List[Dict[str, float]]:
    """
    Per-frame box dimensions and area per lipid, one row per frame.

    Plain floats, not numpy scalars: these end up in the config file, and
    `yaml.safe_dump` cannot represent a numpy float.
    """
    rows = []
    for frame, (lx, ly, lz) in enumerate(lengths_angstrom):
        lx, ly, lz = float(lx), float(ly), float(lz)
        rows.append({
            "frame": frame,
            "box_x": lx,
            "box_y": ly,
            "box_z": lz,
            "area_per_lipid": area_per_lipid(lx, ly, n_lipids),
        })
    return rows


def membrane_report(pdb_file: str, trajfile: str, csvfile: str = "membrane_box.csv") -> Dict[str, Any]:
    """
    Box dimensions and area per lipid over a trajectory, written to `csvfile`.

    Diagnostics for a membrane NPT run: with no `barostat=` passed to ASH (the membrane
    barostat is already in the system) ASH's own CSV leaves out volume and density, and
    the area per lipid is the quantity that says whether the bilayer is equilibrated.
    """
    import csv

    import mdtraj

    traj = mdtraj.load(trajfile, top=pdb_file)
    if traj.unitcell_lengths is None:
        raise ValueError(f"{trajfile} carries no unit-cell information")
    n_lipids = count_lipids(traj.topology.to_openmm())
    rows = box_table(traj.unitcell_lengths * 10.0, n_lipids)  # mdtraj uses nm

    with open(csvfile, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    last = rows[-1]
    half = rows[len(rows) // 2:]
    mean_apl = sum(r["area_per_lipid"] for r in half) / len(half)
    print(f"[INFO] Membrane: {n_lipids} lipids, box "
          f"{last['box_x']:.1f} x {last['box_y']:.1f} x {last['box_z']:.1f} Å, "
          f"area per lipid {last['area_per_lipid']:.1f} Å² "
          f"(mean over the second half of the run: {mean_apl:.1f} Å²)")
    return {
        "lipids": n_lipids,
        "area_per_lipid": last["area_per_lipid"],
        "area_per_lipid_mean": mean_apl,
        "box": [last["box_x"], last["box_y"], last["box_z"]],
        "csv": Path(csvfile).resolve().as_posix(),
    }
