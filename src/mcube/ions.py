"""Li/Merz 12-6 IOD nonbonded parameters for metal ions, TIP3P water.

Amber's stock metal parameters put the metal-oxygen pair minimum too close: for
Zn(II) in `amber14/tip3p.xml` it sits at about 1.69 A, against a real Zn-O bond
of 2.0-2.1 A, which is why a non-coordinating donor can collapse onto the metal
during MM. The IOD set is fitted to reproduce experimental ion-oxygen distances
and moves that minimum to about 1.89 A.

It does not solve over-coordination. A 12-6 potential knows only distance and
atom type, so it cannot tell a ligand from a non-ligand and every oxygen in
reach is still rewarded for approaching; only something that knows which atoms
are ligands can stop that, which is what the coordination restraints do. These
are
simply better metal radii than the defaults.

Values are Amber Rmin/2 (A) and epsilon (kcal/mol), transcribed from
`frcmod.ions234lm_iod_tip3p` as corrected by its authors:
https://ambermd.org/tutorials/advanced/tutorial20/12_6.php - the copy shipped
in AmberTools 2016-2019 carries the 12-6-4 values for +3 and +4 ions under an
IOD label, which are not valid without the accompanying C4 term. The +2 rows
are unaffected in both. Checked against a local AmberTools install: Zn2+,
La3+, Ce3+ and Pr3+ agree exactly.

OpenMM carries the same defect in a force field: `amber14/tip3p.xml` (checked
in OpenMM 8.4; built from AmberTools 17 `frcmod.ions234lm_126_tip3p`) gives
La3+ sigma 0.31324 nm and Lu3+ 0.28954 nm, i.e. Rmin/2 1.758 and 1.625 A -
exactly the 12-6-4 rows below (LM_1264_TIP3P), with no C4 term. IOD has 1.718
and 1.588. So for the lanthanides the stock values are not a valid 12-6 model.
That makes `none` asymmetric: for Zn2+ it is a legitimate control (Amber's
stock radii are only suboptimal), but for a lanthanide it is a known-wrong
parameterisation and cannot be the baseline arm of a comparison.

Sources, as given per row in the file itself:
  +2 ions   Li, Roberts, Chakravorty & Merz, J. Chem. Theory Comput. 2013, 9, 2733
  +3/+4     Li, Song & Merz, J. Phys. Chem. B 2015, 119, 883

The 12-6-4 model (`metal_parameters: 1264`) does better for highly charged
ions: it adds an ion-induced-dipole term -C4/r^4 between the ion and every
other atom, with its own refitted 12-6 terms. It is built the way Amber's
ParmEd `add12_6_4` builds it (read: parmed/tools/add1264.py):

  C4(ion, j) = C4(ion) * alpha_j / alpha_OW,     alpha_OW = 1.444 A^3

with C4(ion) the ion-water value and alpha_j the polarisability of j's Amber
atom type, from `data/lj_1264_pol.dat`. Between two such ions the two terms add.

  12-6 terms     frcmod.ions234lm_1264_tip3p: +2 rows from Li & Merz, J. Chem.
                 Theory Comput. 2014, 10, 289; +3/+4 rows from Li, Song & Merz,
                 J. Phys. Chem. B 2015, 119, 883
  C4 values      ParmEd DEFAULT_C4_PARAMS['TIP3P'] (the same paper tables)
  polarisability lj_1264_pol.dat, copied unmodified from the AmberTools in the
                 ASH environment (sha256 9a64b54f...); it covers the protein,
                 GAFF, Lipid17, water and ion types

Every row was checked against that AmberTools install (dat/leap/parm).

C4 is not fitted to protein ligands, and the alpha-ratio rule tends to
over-bind carboxylates (Jafari et al., J. Phys. Chem. B 2024, 128, 684).
"""

from __future__ import annotations

from pathlib import Path

# element -> charge -> (Rmin/2 in Angstrom, epsilon in kcal/mol)
IOD_TIP3P = {
    "BE": {2: (1.168, 0.00063064)},
    "CU": {2: (1.409, 0.01721000)},
    "NI": {2: (1.373, 0.01179373)},
    "ZN": {2: (1.395, 0.01491700)},
    "CO": {2: (1.404, 0.01636246)},
    "CR": {2: (1.388, 0.01386171), 3: (1.344, 0.00848000)},
    "FE": {2: (1.409, 0.01721000), 3: (1.386, 0.01357097)},
    "MG": {2: (1.395, 0.01491700)},
    "V":  {2: (1.476, 0.03198620)},
    "MN": {2: (1.467, 0.02960343)},
    "HG": {2: (1.575, 0.06751391)},
    "CD": {2: (1.506, 0.04090549)},
    "CA": {2: (1.608, 0.08337961)},
    "SN": {2: (1.738, 0.16500296)},
    "SR": {2: (1.753, 0.17618319)},
    "BA": {2: (1.913, 0.31060194)},
    "AL": {3: (1.297, 0.00471279)},
    "IN": {3: (1.461, 0.02808726)},
    "TL": {3: (1.513, 0.04321029)},
    "Y":  {3: (1.602, 0.08034231)},
    "LA": {3: (1.718, 0.15060822)},
    "CE": {3: (1.741, 0.16721338), 4: (1.684, 0.12758274)},
    "PR": {3: (1.733, 0.16134811)},
    "ND": {3: (1.681, 0.12564307)},
    "SM": {3: (1.659, 0.11189491)},
    "EU": {3: (1.666, 0.11617738)},
    "GD": {3: (1.623, 0.09126804)},
    "TB": {3: (1.630, 0.09509276)},
    "DY": {3: (1.609, 0.08389240)},
    "ER": {3: (1.602, 0.08034231)},
    "TM": {3: (1.602, 0.08034231)},
    "LU": {3: (1.588, 0.07351892)},
    "HF": {4: (1.499, 0.03868661)},
    "ZR": {4: (1.519, 0.04525501)},
    "U":  {4: (1.684, 0.12758274)},
    "PU": {4: (1.662, 0.11371963)},
    "TH": {4: (1.708, 0.14364160)},
}

# element -> charge -> (Rmin/2 in Angstrom, epsilon in kcal/mol, C4 in kcal/mol/A^4)
LM_1264_TIP3P = {
    "BE": {2: (1.203, 0.00116124, 186.5)},
    "CU": {2: (1.476, 0.03198620, 290.9)},
    "NI": {2: (1.431, 0.02133669, 212.8)},
    "ZN": {2: (1.455, 0.02662782, 231.6)},
    "CO": {2: (1.458, 0.02735051, 209.7)},
    "CR": {2: (1.431, 0.02133669, 136.8), 3: (1.415, 0.01827024, 258.0)},
    "FE": {2: (1.457, 0.02710805, 163.0), 3: (1.443, 0.02387506, 428.0)},
    "MG": {2: (1.437, 0.02257962, 132.9)},
    "V":  {2: (1.494, 0.03715368, 195.7)},
    "MN": {2: (1.485, 0.03450196, 146.1)},
    "HG": {2: (1.641, 0.10128575, 288.8)},
    "CD": {2: (1.535, 0.05102457, 185.6)},
    "CA": {2: (1.642, 0.10185975, 87.3)},
    "SN": {2: (1.777, 0.19470705, 187.9)},
    "SR": {2: (1.777, 0.19470705, 82.7)},
    "BA": {2: (1.936, 0.33132862, 71.9)},
    "AL": {3: (1.369, 0.01128487, 399.0)},
    "IN": {3: (1.491, 0.03625449, 347.0)},
    "TL": {3: (1.571, 0.06573030, 456.0)},
    "Y":  {3: (1.630, 0.09509276, 216.0)},
    "LA": {3: (1.758, 0.17997960, 152.0)},
    "CE": {3: (1.782, 0.19865859, 230.0), 4: (1.766, 0.18612361, 706.0)},
    "PR": {3: (1.780, 0.19707431, 264.0)},
    "ND": {3: (1.724, 0.15486311, 213.0)},
    "SM": {3: (1.711, 0.14571499, 230.0)},
    "EU": {3: (1.716, 0.14920231, 259.0)},
    "GD": {3: (1.658, 0.11129023, 198.0)},
    "TB": {3: (1.671, 0.11928915, 235.0)},
    "DY": {3: (1.637, 0.09900804, 207.0)},
    "ER": {3: (1.635, 0.09788018, 251.0)},
    "TM": {3: (1.647, 0.10475707, 282.0)},
    "LU": {3: (1.625, 0.09235154, 249.0)},
    "HF": {4: (1.600, 0.07934493, 827.0)},
    "ZR": {4: (1.609, 0.08389240, 761.0)},
    "U":  {4: (1.792, 0.20665151, 1034.0)},
    "PU": {4: (1.752, 0.17542802, 828.0)},
    "TH": {4: (1.770, 0.18922704, 512.0)},
}

SETS = {"iod": IOD_TIP3P, "1264": LM_1264_TIP3P}

WATER_POL = 1.444  # A^3, water oxygen (OW) in lj_1264_pol.dat

POL_FILE = Path(__file__).with_name("data") / "lj_1264_pol.dat"


def lookup(name: str, element: str, charge: float):
    """(sigma in nm, epsilon in kJ/mol) for this ion, or None if not tabulated."""
    table = SETS.get(str(name).lower())
    if table is None:
        raise ValueError(f"unknown metal parameter set {name!r}; available: "
                         f"{', '.join(sorted(SETS))}")
    entry = table.get(element.upper(), {}).get(int(round(charge)))
    if entry is None:
        return None
    rmin_half, eps_kcal = entry[:2]
    return (2.0 * rmin_half / 2.0 ** (1.0 / 6.0) / 10.0, eps_kcal * 4.184)


def metal_parameter_set(inputs: dict, md: dict) -> str:
    """The metal parameter set for a stage after md: its own input, else what md used."""
    which = inputs.get("metal_parameters")
    if which is None:
        which = (md.get("results") or {}).get("metal_parameters")
    if which is None:
        which = (md.get("inputs") or {}).get("metal_parameters", "iod")
    return str(which).lower()


def c4_value(element: str, charge: float):
    """The 12-6-4 ion-water C4 (kcal/mol/A^4), or None if not tabulated."""
    entry = LM_1264_TIP3P.get(element.upper(), {}).get(int(round(charge)))
    return None if entry is None else entry[2]


C4_PAIRS_FILE = Path(__file__).with_name("data") / "c4_pairs.yaml"


def c4_sets(path=C4_PAIRS_FILE) -> dict:
    """The literature pair-C4 sets shipped with mcube (data/c4_pairs.yaml)."""
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def ion_key(element: str, charge: float) -> str:
    """'Ca2+', 'La3+', ...: how c4_pairs names an ion."""
    element = element.strip()
    return f"{element[0].upper()}{element[1:].lower()}{int(round(charge))}+"


def c4_overrides(inputs: dict, md_inputs: dict = None) -> dict:
    """
    Pair-specific C4 for this stage: {ion key: {atom class: C4 kcal/mol/A^4}}.

    `c4_set` names a set in data/c4_pairs.yaml; `c4_pairs` gives pairs explicitly and wins
    over the set, pair by pair. A stage's own keys win over md's, as for
    `polarizabilities`, so QM/MM and umbrella stages see what md used.
    """
    md_inputs = md_inputs or {}
    name = inputs.get("c4_set", md_inputs.get("c4_set"))
    explicit = inputs.get("c4_pairs", md_inputs.get("c4_pairs")) or {}
    pairs: dict = {}
    if name:
        sets = c4_sets()
        if name not in sets:
            raise ValueError(f"unknown c4_set {name!r}; available: {', '.join(sorted(sets))}")
        chosen = sets[name]
        print(f"[INFO] 12-6-4 c4_set {name}: {chosen.get('water_model', '?')} water, fitted to "
              f"{chosen.get('fitted_to', '?')}")
        for ion, classes in (chosen.get("pairs") or {}).items():
            pairs.setdefault(str(ion), {}).update({str(k): float(v) for k, v in classes.items()})
    for ion, classes in explicit.items():
        pairs.setdefault(str(ion), {}).update({str(k): float(v) for k, v in classes.items()})
    return pairs


def read_polarizabilities(path=POL_FILE) -> dict:
    """{Amber atom type: polarisability in A^3} from a two-column file."""
    pol = {}
    for line in open(path):
        parts = line.split()
        if len(parts) >= 2:
            pol[parts[0]] = float(parts[1])
    return pol


def type_polarizability(atom_class: str, table: dict, overrides=None):
    """
    Polarisability of an OpenMM atom class, or None if unknown.

    OpenMM's Amber XMLs keep the Amber type as the class for proteins and
    lipids, but prefix water and ions (`tip3p-O`, `tip3p_standard-Na+`), so those
    are translated first. `overrides` (by class) wins over all.
    """
    if overrides and atom_class in overrides:
        return float(overrides[atom_class])
    name = {"tip3p-O": "OW", "tip3p-H": "HW"}.get(atom_class, atom_class)
    if name.startswith("tip3p_standard-"):
        name = name[len("tip3p_standard-"):]
    return table.get(name)


def atom_classes(topology, forcefields) -> list:
    """The force-field atom class of every atom, in topology order.

    Uses the matching OpenMM's own createSystem does. `_getResidueTemplateMatches`
    is private API, but it is the only way to know which template atom each
    residue atom became, and atom names do not have to agree with the template.
    """
    from openmm.app import ForceField
    ff = ForceField(*forcefields)
    bonded = ff._buildBondedToAtomList(topology)
    classes = [None] * topology.getNumAtoms()
    for residue in topology.residues():
        template, matches = ff._getResidueTemplateMatches(residue, bonded)
        if matches is None:
            raise ValueError(f"No force-field template for residue {residue.name} "
                             f"{residue.chain.id}:{residue.id}")
        for atom, k in zip(residue.atoms(), matches):
            classes[atom.index] = ff._atomTypes[template.atoms[k].type].atomClass
    return classes


def add_c4_force(system, topology, forcefields, ion_indices, ion_c4, qm_atoms=(),
                 polarizabilities=None, ion_keys=None, c4_pairs=None):
    """
    Add the 12-6-4 ion-induced-dipole term as a CustomNonbondedForce. Returns it.

    Per particle: `cfour` (C4, nonzero only for the ions) and `alpha` (alpha /
    alpha_OW). The pair energy -(cfour1*alpha2 + cfour2*alpha1)/r^4 is C4(ion)*alpha_j/alpha_OW
    for an ion and any atom, and the sum of both terms for two ions, as in
    ParmEd. It is restricted to ion x everything by an interaction group, and
    takes the NonbondedForce's cutoff and exceptions.

    `qm_atoms` get cfour = alpha = 0: in QM/MM the QM region is polarised by the
    embedding already, so a QM ion's C4 would count that induction twice.

    `c4_pairs` ({ion key: {atom class: C4}}, see `c4_overrides`) replaces the default C4 for
    those (ion, class) pairs: a table `dc4(islot, cslot)` adds (pair C4 - default C4), where
    `islot` numbers the overridden ions and `cslot` the overridden classes (0 = none).
    """
    import openmm

    existing = next((f for f in system.getForces() if f.getName() == "mcube_1264_C4"), None)
    if existing is not None:
        return existing  # already added to this system (stages may rebuild their theory)
    table = read_polarizabilities()
    classes = atom_classes(topology, forcefields)
    ions = {int(i): float(c) for i, c in zip(ion_indices, ion_c4)}
    qm = {int(a) for a in qm_atoms}

    pol, missing = [], set()
    for index, atom_class in enumerate(classes):
        value = type_polarizability(atom_class, table, polarizabilities)
        if value is None:
            # A C4 ion's own polarisability matters only against another C4 ion
            if index not in ions:
                missing.add(atom_class)
            value = 0.0
        pol.append(value)
    if missing:
        raise ValueError(
            "12-6-4 needs a polarisability for every atom type; none for "
            f"{', '.join(sorted(missing))}. Give them as `polarizabilities: "
            "{class: A^3}` in the stage inputs.")

    nonbonded = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))

    # Pair overrides: one row per overridden ion particle, one column per overridden class
    keys = {int(i): k for i, k in zip(ion_indices, ion_keys or [])}
    c4_pairs = c4_pairs or {}
    present = set(classes)
    over_ions = [i for i in sorted(ions) if keys.get(i) in c4_pairs and i not in qm]
    over_classes = sorted({c for i in over_ions for c in c4_pairs[keys[i]] if c in present})
    applied = {f"{keys[i]}:{c}" for i in over_ions for c in c4_pairs[keys[i]] if c in present}
    unused = {f"{k}:{c}" for k, cl in c4_pairs.items() for c in cl} - applied
    if unused:
        print("[WARNING] c4_pairs not applied (no such ion or atom class in this system): "
              + ", ".join(sorted(unused)))
    islot = {i: n + 1 for n, i in enumerate(over_ions)}
    cslot = {c: n + 1 for n, c in enumerate(over_classes)}
    nx, ny = len(over_ions) + 1, len(over_classes) + 1
    delta = [0.0] * (nx * ny)
    for i in over_ions:
        for c, value in c4_pairs[keys[i]].items():
            if c not in cslot:
                continue
            pol_c = type_polarizability(c, table, polarizabilities) or 0.0
            default = ions[i] * pol_c / WATER_POL
            # kcal/mol/A^4 -> kJ/mol/nm^4
            delta[islot[i] + nx * cslot[c]] = (float(value) - default) * 4.184e-4
            print(f"[INFO] 12-6-4 pair C4 {keys[i]}-{c}: {default:.2f} -> {float(value):.2f} kcal/mol/A^4")

    energy = "-(cfour1*alpha2 + cfour2*alpha1)/r^4"
    if over_ions:
        energy = "-(cfour1*alpha2 + cfour2*alpha1 + dc4(islot1, cslot2) + dc4(islot2, cslot1))/r^4"
    force = openmm.CustomNonbondedForce(energy)
    force.setName("mcube_1264_C4")
    force.addPerParticleParameter("cfour")
    force.addPerParticleParameter("alpha")
    if over_ions:
        force.addPerParticleParameter("islot")
        force.addPerParticleParameter("cslot")
        force.addTabulatedFunction("dc4", openmm.Discrete2DFunction(nx, ny, delta))
    for index in range(system.getNumParticles()):
        if index in qm:
            params = [0.0, 0.0, 0.0, 0.0]
        else:
            # kcal/mol/A^4 -> kJ/mol/nm^4
            params = [ions.get(index, 0.0) * 4.184e-4, pol[index] / WATER_POL,
                      float(islot.get(index, 0)), float(cslot.get(classes[index], 0))]
        force.addParticle(params if over_ions else params[:2])
    # Every exception, not only the ions': OpenMM's CPU and GPU platforms refuse
    # nonbonded forces whose exclusions differ ("All Forces must have identical
    # exclusions"); OpenMM's own Amber 12-6-4 path copies them all the same way
    for k in range(nonbonded.getNumExceptions()):
        i, j = nonbonded.getExceptionParameters(k)[:2]
        force.addExclusion(i, j)
    method = nonbonded.getNonbondedMethod()
    if method == openmm.NonbondedForce.NoCutoff:
        force.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
    elif method == openmm.NonbondedForce.CutoffNonPeriodic:
        force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffNonPeriodic)
    else:
        force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
        # As OpenMM's own Amber 12-6-4 path (AmberPrmtopFile) does
        force.setUseLongRangeCorrection(True)
    force.setCutoffDistance(nonbonded.getCutoffDistance())
    force.addInteractionGroup(sorted(ions), range(system.getNumParticles()))
    system.addForce(force)
    return force


def apply_metal_parameters(system, topology, metals, which: str, forcefields=None,
                           qm_atoms=(), polarizabilities=None, c4_pairs=None) -> list:
    """
    Override the metals' Lennard-Jones terms with a published ion set, and for
    `1264` add the C4 term. Returns the indices changed.

    The one place this happens, for md and every QM/MM stage alike, so a metal
    cannot be parameterised one way in MM and another in QM/MM. The ion is
    identified from the structure's element and the charge the force field
    assigned. Under IOD an untabulated metal is left alone with a warning; under
    1264 it raises, because 12-6-4 terms without the C4, or the reverse, are
    neither model.
    """
    import openmm
    from openmm import unit

    which = str(which).lower()
    nonbonded = next((f for f in system.getForces()
                      if isinstance(f, openmm.NonbondedForce)), None)
    if nonbonded is None:
        print("[WARNING] no NonbondedForce: metal parameters unchanged")
        return []

    wanted = {m.upper() for m in (metals or [])}
    changed, c4s, keys = [], [], []
    for atom in topology.atoms():
        element = (atom.element.symbol if atom.element else "").upper()
        if element not in wanted:
            continue
        charge, sigma, epsilon = nonbonded.getParticleParameters(atom.index)
        q = charge.value_in_unit(unit.elementary_charge)
        entry = lookup(which, element, q)
        if entry is None:
            if which == "1264":
                raise ValueError(f"no 1264 parameters for {element} with charge {q:+.0f} "
                                 "(+2 to +4 ions are tabulated; see ions.py)")
            print(f"[WARNING] no {which} parameters for {element} with charge "
                  f"{charge}: left at the force field's own values")
            continue
        new_sigma, new_epsilon = entry
        nonbonded.setParticleParameters(atom.index, charge,
                                        new_sigma * unit.nanometer,
                                        new_epsilon * unit.kilojoule_per_mole)
        changed.append(atom.index)
        keys.append(ion_key(element, q))
        note = ""
        if which == "1264":
            c4s.append(c4_value(element, q))
            note = f", C4 {c4s[-1]} kcal/mol/A^4"
        print(f"[INFO] {which.upper()} parameters for {element} "
              f"{atom.residue.chain.id}:{atom.residue.id}: sigma "
              f"{sigma.value_in_unit(unit.nanometer):.4f} -> {new_sigma:.4f} nm, "
              f"epsilon {epsilon.value_in_unit(unit.kilojoule_per_mole):.4f} -> "
              f"{new_epsilon:.4f} kJ/mol{note}")

    if c4_pairs and which != "1264":
        raise ValueError("c4_set / c4_pairs only apply to metal_parameters: 1264")
    if which == "1264" and changed:
        if forcefields is None:
            raise ValueError("1264 needs the force-field files, to type every atom")
        add_c4_force(system, topology, forcefields, changed, c4s, qm_atoms, polarizabilities,
                     ion_keys=keys, c4_pairs=c4_pairs)
        if set(changed) & {int(a) for a in qm_atoms}:
            print("[INFO] 12-6-4: C4 is off for QM atoms (the embedding polarises them)")
    return changed
