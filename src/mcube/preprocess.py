#!/usr/bin/env python3
"""
Protein/Ligand preprocessing utilities (pipeline-ready).
"""

from __future__ import annotations

import argparse
import math
import yaml
import os
import shutil
from contextlib import contextmanager
from typing import Dict, Any, Tuple
from pathlib import Path
from datetime import datetime

try:
    import fcntl
except ImportError:  # Windows: no locking (single-user use)
    fcntl = None

from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds
from openmm.app import PDBFile, Topology, ForceField

from ash import (
    Fragment,
    OpenMM_Modeller,
    small_molecule_parameterizer,
    merge_pdb_files,
)

from .membrane import LIPID_XML, check_lipid_type, lipid_forcefields

# ==============================
# Utility
# ==============================


def save_result(
    result: Dict[str, Any],
    filename: str = "config.yaml",
    section: str = None,
) -> str:
    """
    Append pipeline stage outputs to a shared YAML file.

    Assumes:
    - "inputs" already exist in YAML
    - Each stage only writes its own results

    Parameters
    ----------
    result : dict
        Outputs from a pipeline stage
    filename : str
        YAML file path
    section : str
        Stage name ("preprocess", "md", "qmmm")

    Returns
    -------
    str
        Path to YAML file
    """

    path = Path(filename)

    # ==============================
    # 1. Clean result
    # ==============================
    clean_result = {}

    for k, v in result.items():

        if isinstance(v, (str, int, float, bool)) or v is None:
            clean_result[k] = v

        elif isinstance(v, (list, dict)):
            clean_result[k] = v

        elif isinstance(v, set):
            clean_result[k] = sorted(v)

        else:
            continue  # skip non-serializable objects

    # ==============================
    # 2. Stage payload
    # ==============================
    payload = {
        "results": clean_result,
        "timestamp": datetime.now().isoformat(),
    }

    path.parent.mkdir(parents=True, exist_ok=True)

    # Read-modify-write under a lock, so several jobs (e.g. one QM/MM run per
    # snapshot, each with its own --section) can share one config file
    with _locked(path):
        # ==============================
        # 3. Load existing YAML
        # ==============================
        if path.exists():
            with open(path, "r") as f:
                data = yaml.safe_load(f) or {}
        else:
            data = {}

        # ==============================
        # 4. Merge into section
        # ==============================
        if section:
            data.setdefault(section, {})
            data[section].update(payload)
        else:
            data.update(payload)

        # ==============================
        # 5. Atomic write
        # ==============================
        tmp_path = path.with_name(f"{path.name}.{os.getpid()}.tmp")

        with open(tmp_path, "w") as f:
            yaml.safe_dump(
                data,
                f,
                sort_keys=False,
                default_flow_style=False,
            )

        tmp_path.replace(path)

    print(f"[INFO] Updated → {path} [{section}]")

    return str(path)


@contextmanager
def _locked(path: Path):
    """Exclusive lock on `<config>.lock` (POSIX; no-op where fcntl is missing)."""
    if fcntl is None:
        yield
        return
    with open(path.with_name(path.name + ".lock"), "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def load_section(filename: str, section: str) -> Dict[str, Any]:
    """
    Load a section from a YAML results file.

    Parameters
    ----------
    filename : str
        YAML file path
    section : str
        Section to load

    Returns
    -------
    dict
    """

    with open(filename, "r") as f:
        data = yaml.safe_load(f) or {}

    return data.get(section, {})


def atom_names(pdb_file: str) -> list:
    """["chain:resid:atom", ...] in file order: the inverse of a selector lookup."""
    return ["%s:%s:%s" % (line[21], line[22:26].strip(), line[12:16].strip())
            for line in open(pdb_file) if line.startswith(("ATOM", "HETATM"))]


def to_selectors(entries, pdb_file: str) -> list:
    """Rewrite each entry's first two atoms as "chain:resid:atom", leaving any
    trailing numbers (a target distance, a force constant) alone.

    Indices are only meaningful against the structure they were resolved on, and
    solvating again renumbers everything after the solute. Selectors survive that,
    which is why they are what the config uses everywhere else. Entries that are
    already selectors pass through, so configs written before this still work.
    """
    if not entries:
        return []
    names = None
    out = []
    for entry in entries:
        atoms = entry.replace(",", " ").split() if isinstance(entry, str) else list(entry)
        if len(atoms) < 2:
            continue
        pair = []
        for atom in atoms[:2]:
            if isinstance(atom, str) and atom.count(":") == 2:
                pair.append(atom)
            else:
                if names is None:
                    names = atom_names(pdb_file)
                index = int(atom)
                pair.append(names[index] if index < len(names) else str(atom))
        out.append(pair + list(atoms[2:]))
    return out


def check_inherited_restraints(constraints, restraints, pdb_file: str, inputs) -> None:
    """Refuse to constrain a pair that md held in place with a harmonic restraint.

    A restraint says "we are imposing this distance"; a constraint on a detected
    contact says "this distance is structural". Re-running preprocess after a
    restrained md turns the first into the second, and nothing downstream can
    then tell them apart. Both lists are already in the config, so this only
    compares them.

    It sees one config. A run started from another run's structure, with its own
    config, is outside what this can check: that restraint is recorded in the
    other config.
    """
    if not constraints or not restraints:
        return
    imposed = {frozenset(entry[:2]): entry for entry in to_selectors(restraints, pdb_file)}
    allowed = {frozenset(entry[:2]) for entry in
               to_selectors(inputs.get("acknowledge_restrained"), pdb_file)}
    clashes = [imposed[key] for key in
               (frozenset(entry[:2]) for entry in to_selectors(constraints, pdb_file))
               if key in imposed and key not in allowed]
    if not clashes:
        return

    detail = "".join("    %s - %s: restrained to %s A by md\n"
                     % (c[0], c[1], c[2] if len(c) > 2 else "?") for c in clashes)
    raise ValueError(
        "md.results.restraints records distances this run imposed, and preprocess "
        "would now freeze them as detected metal coordination:\n" + detail +
        "\nThat makes an assumption look structural. In the preprocess inputs, either\n"
        "add the residue to `constraint_exclude_residues` to leave it free, or list\n"
        "the pair in `acknowledge_restrained` to assert the coordination is real:\n"
        '    acknowledge_restrained: [["%s", "%s"]]\n' % (clashes[0][0], clashes[0][1]))


def coordination_shell(pdb_file: str, inputs, radius=None, per_residue: bool = True) -> list:
    """The metal coordination sphere: the closest donor of each residue within `cutoff`.

    The one place that rule lives. `find_coordinators` builds the MD constraints
    and the QM region from it, and `md` reports it before and after the run, so a
    contact that appears or disappears during equilibration is visible instead of
    being inherited silently. An MD frame is wrapped into the periodic cell, so
    the minimum image is applied when the file carries an orthorhombic CRYST1.
    """
    metals = {m.upper() for m in (inputs.get("metals") or [])}
    donors = {e.upper() for e in (inputs.get("donor_elements") or ["N", "O", "S"])}
    cutoff = float(radius if radius is not None else (inputs.get("cutoff") or 2.8))

    rows, cell = [], None
    for line in open(pdb_file):
        if line.startswith(("ATOM", "HETATM")):
            element = (line[76:78].strip() or line[12:16].strip()[:1]).upper()
            rows.append((line[21], line[22:26].strip(), line[12:16].strip(),
                         line[17:20].strip(), element,
                         (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
        elif line.startswith("CRYST1"):
            if all(abs(float(line[c:c + 7]) - 90.0) < 1e-3 for c in (33, 40, 47)):
                cell = [float(line[6:15]), float(line[15:24]), float(line[24:33])]

    def separation(a, b):
        v = [q - p for p, q in zip(a, b)]
        if cell:
            v = [x - c * round(x / c) if c else x for x, c in zip(v, cell)]
        return math.sqrt(sum(x * x for x in v))

    found = []
    for mi, metal in enumerate(rows):
        if metal[4] not in metals:
            continue
        closest = {}
        for di, donor in enumerate(rows):
            if di == mi or donor[4] not in donors:
                continue
            r = separation(metal[5], donor[5])
            if r > cutoff:
                continue
            # Per residue for the QM region, which takes whole residues anyway.
            # Per atom for constraints: a carboxylate can rotate its second
            # oxygen in while the first is held, which is how Asp92 OD1 reached
            # 1.75 A with OD2 pinned at 2.80.
            key = (donor[0], donor[1], donor[3]) if per_residue else di
            if key not in closest or r < closest[key][0]:
                closest[key] = (r, di, donor[2])
        for (r, di, aname) in sorted(closest.values()):
            donor = rows[di]
            found.append({
                "metal": "%s:%s:%s" % (metal[0], metal[1], metal[2]),
                "donor": "%s:%s:%s" % (donor[0], donor[1], aname),
                "resname": donor[3],
                "distance": round(r, 2),
                "metal_index": mi,
                "donor_index": di,
            })
    return found


def report_coordination(shell, label: str) -> None:
    print("\n===== Metal coordination: %s =====" % label)
    if not shell:
        print("  none within the cutoff")
    for row in shell:
        print("  %-12s %-14s %-4s %5.2f A"
              % (row["metal"], row["donor"], row["resname"], row["distance"]))


def compare_coordination(before: list, after: list) -> list:
    """Lines describing how the sphere changed, empty if it did not.

    A ligand that appears during MD was not there to be constrained at setup, so
    nothing held it: a non-bonded point-charge metal will pull a neutral donor in
    and there is no Pauli repulsion to stop it. That has happened here (a Gln
    amide oxygen went from 5.54 to 1.79 A), and it is invisible unless something
    looks.
    """
    was = {(r["metal"], r["donor"]): r["distance"] for r in before or []}
    now = {(r["metal"], r["donor"]): r["distance"] for r in after or []}
    lines = []
    for key in sorted(set(now) - set(was)):
        lines.append("  APPEARED  %s - %s at %.2f A" % (key[0], key[1], now[key]))
    for key in sorted(set(was) - set(now)):
        lines.append("  LEFT      %s - %s (was %.2f A)" % (key[0], key[1], was[key]))
    for key in sorted(set(was) & set(now)):
        if abs(now[key] - was[key]) >= 0.2:
            lines.append("  MOVED     %s - %s %.2f -> %.2f A"
                         % (key[0], key[1], was[key], now[key]))
    return lines


def resolve_cores(requested, label: str = "cores") -> int:
    """How many threads a stage asks ASH for: the config value, else the job's.

    Thread *placement* -- OMP_NUM_THREADS, OMP_PROC_BIND, OMP_PLACES -- is set by
    the job script, not by mcube, because the right policy differs per stage and
    cannot be changed once the process has started: see the thread recipes in the
    README and `docs/csr/run.slurm`. This only decides a count.
    """
    if requested is not None:
        return int(requested)
    for var in ("OMP_NUM_THREADS", "SLURM_CPUS_PER_TASK"):
        value = (os.environ.get(var) or "").strip()
        if value.isdigit() and int(value) > 0:
            print(f"[INFO] {label} not set: using {value} from ${var}")
            return int(value)
    print(f"[INFO] {label} not set and neither $OMP_NUM_THREADS nor "
          f"$SLURM_CPUS_PER_TASK is set: using 1")
    return 1


# ==============================
# Atom selectors
# ==============================

def resolve_atoms(value, topology):
    """
    Replace atom selectors in a config value by 0-based atom indices.

    An atom can be given as an index (int) or as "chain:resid:atom", e.g. "B:1:ZN"
    or "A:88:NE2", matched against `topology` (an OpenMM Topology or a PDB file).
    Lists are resolved element by element, so constraint/restraint entries such as
    ["B:1:ZN", "C:1:O2x", 2.2, 20.0] keep their numbers. This lets configs name
    atoms instead of hard-coding indices that change with every system.
    """
    if not _has_selector(value):
        return value
    if isinstance(topology, (str, Path)):
        topology = PDBFile(str(topology)).topology
    lookup = {}
    for atom in topology.atoms():
        res = atom.residue
        lookup.setdefault((res.chain.id, res.id.strip(), atom.name), atom.index)

    def resolve(v):
        if isinstance(v, str):
            parts = [p.strip() for p in v.split(":")]
            if len(parts) != 3:
                raise ValueError(f"Atom selector {v!r} must be 'chain:resid:atom', e.g. 'A:88:NE2'")
            if tuple(parts) not in lookup:
                raise ValueError(f"No atom {v!r} (chain {parts[0]}, residue {parts[1]}, atom {parts[2]}) "
                                 "in the structure")
            return lookup[tuple(parts)]
        if isinstance(v, (list, tuple)):
            return [resolve(x) for x in v]
        return v

    resolved = resolve(value)
    print(f"[INFO] Atom selectors resolved: {value} -> {resolved}")
    return resolved


def _has_selector(value) -> bool:
    if isinstance(value, str):
        return True
    if isinstance(value, (list, tuple)):
        return any(_has_selector(v) for v in value)
    return False


def resolve_atom_inputs(inputs: Dict[str, Any], topology, keys) -> Dict[str, Any]:
    """Copy of `inputs` with the atom selectors under `keys` replaced by indices."""
    out = dict(inputs)
    for k in keys:
        if out.get(k) is not None:
            out[k] = resolve_atoms(out[k], topology)
    return out


# ==============================
# Ligand Preparation
# ==============================

def ligand_resname(inputs: Dict[str, Any]) -> str:
    """Residue name used for the ligand: `ligand_resname` input, else the ligand file stem."""
    if inputs.get("ligand_resname"):
        return inputs["ligand_resname"]
    if inputs.get("ligand"):
        return Path(inputs["ligand"]).stem
    return "LIG"


def ligand_mol_from_pdb(ligand: str, charge: int, smiles: str = None) -> Chem.Mol:
    """
    Read a ligand PDB with explicit hydrogens, taking connectivity only from its CONECT
    records (no distance-based bonding, which adds spurious bonds in strained geometries).
    Bond orders come from the `smiles` template if given, else from RDKit DetermineBondOrders.
    """
    mol = Chem.MolFromPDBFile(ligand, removeHs=False, sanitize=False, proximityBonding=False)
    if mol is None:
        raise ValueError(f"RDKit could not read {ligand}")
    if mol.GetNumBonds() == 0:
        raise ValueError(f"{ligand} has no CONECT records; bonds are needed to parameterise it")

    if smiles:
        template = Chem.AddHs(Chem.MolFromSmiles(smiles))
        try:
            mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
        except ValueError as e:
            raise ValueError(f"ligand_smiles does not match the atoms/bonds of {ligand} "
                             "(check hydrogens and protonation)") from e
    else:
        rdDetermineBonds.DetermineBondOrders(mol, charge=charge)

    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    if Chem.GetFormalCharge(mol) != charge:
        raise ValueError(f"{ligand} has formal charge {Chem.GetFormalCharge(mol)} but charge is {charge}")
    # Printed so the bond-order assignment can be checked
    print("[INFO] Ligand SMILES:", Chem.MolToSmiles(Chem.RemoveHs(mol)))
    return mol


def merge_ligand(protein_pdb: str, ligand: str, model: int = 1, charge: int = 0,
                 resname: str = None, smiles: str = None, merge: bool = True) -> Tuple[str, str]:
    """
    Parameterise the ligand with OpenFF and (if merge) merge it into protein_pdb.
    With merge=False the ligand is assumed to be in protein_pdb already.
    """
    filename = resname or Path(ligand).stem
    if len(filename) > 3:
        print(f"[WARNING] Ligand residue name '{filename}' is longer than 3 characters and will be "
              "truncated in PDB files, so it may not match the generated XML template. "
              "Set `ligand_resname` in the preprocess inputs.")
    if ligand.endswith(".pdbqt"):

        smiles, idx, coords = None, [], []
        model_index = 0

        with open(ligand) as f:
            for line in f:
                if model_index != model:
                    if line.startswith("MODEL"):
                        model_index += 1
                    continue

                if line.startswith("ENDMDL"):
                    break

                if line.startswith("REMARK SMILES ") and "IDX" not in line:
                    smiles = line.split("SMILES ")[1].strip()
                elif line.startswith("REMARK SMILES IDX"):
                    idx += list(map(int, line.split()[3:]))
                elif line.startswith(("ATOM", "HETATM")):
                    coords.append(tuple(map(float, [line[30:38], line[38:46], line[46:54]])))

        if smiles is None:
            raise ValueError("No SMILES found in ligand file")

        mol = Chem.MolFromSmiles(smiles)
        mapping = {idx[i]-1: idx[i+1]-1 for i in range(0, len(idx), 2)}
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(i, coords[mapping[i]])
        mol.AddConformer(conf)
        mol = Chem.AddHs(mol, addCoords=True)
        Chem.MolToMolFile(mol, f"{filename}.mol")
        small_molecule_parameterizer(forcefield_option="OpenFF", molfile=f"{filename}.mol", charge=charge,
                                     resname=filename)

    elif ligand.endswith(".pdb"):
        # Needs explicit hydrogens and CONECT records
        mol = ligand_mol_from_pdb(ligand, charge, smiles)
        Chem.MolToMolFile(mol, f"{filename}.sdf")
        small_molecule_parameterizer(forcefield_option="OpenFF", sdffile=f"{filename}.sdf", charge=charge,
                                     resname=filename)

    elif ligand.endswith(".mol"):
        small_molecule_parameterizer(forcefield_option="OpenFF", molfile=ligand, charge=charge, resname=filename)

    else:
        raise ValueError(f"Not a valid file format for ligand {ligand}: choose .pdbqt, .pdb or .mol")

    ligand_xmlfile = Path(f"openff_{filename}.xml").resolve().as_posix()

    if merge:
        merge_pdb_files(protein_pdb, f"{filename}.pdb")
        protein_pdb = "merged.pdb"

    return protein_pdb, ligand_xmlfile

# ==============================
# Coordination
# ==============================

def find_coordinators(
    pdb_file: str,
    inputs=None,
):
    pdb = PDBFile(pdb_file)
    topology = pdb.topology
    atoms = list(topology.atoms())
    exclude_residues=inputs.get("exclude_residues")
    # Residues whose metal-donor distance is left free (no bond constraint) during MD
    constraint_exclude = inputs.get("constraint_exclude_residues") or {}
    cutoff = inputs.get("cutoff")

    elems = [a.element.symbol if a.element else "" for a in atoms]

    # Chain-aware residue info
    resids = [int(a.residue.id) for a in atoms]
    resnames = [a.residue.name for a in atoms]
    chains = [a.residue.chain.id for a in atoms]

    metal_idx = [i for i, e in enumerate(elems) if e.upper() in inputs.get("metals")]

    if inputs.get("donor_atoms"):
        print("[WARNING] donor_atoms is no longer used; donors are chosen by element (donor_elements)")

    qm_residue_dict = {}  # OpenMM-style: {"A": {80: "HIS"}}
    bondconstraints = []

    print("\n===== Metal Coordination Analysis =====")
    constraint_cutoff = float(inputs.get("constraint_cutoff") or cutoff)
    print(f"Cutoff: {cutoff:.2f} Å (QM region); {constraint_cutoff:.2f} Å (constraints)")

    # Two radii, because they answer different questions. The QM region takes
    # whole residues that coordinate the metal. The constraints have to cover
    # everything that could *become* a contact during MM, since a non-bonded
    # metal will pull a donor in and only an existing constraint stops it;
    # `constraint_cutoff` widens that net without dragging more residues into
    # the QM region. Per atom rather than per residue, so both oxygens of a
    # carboxylate are held.
    shell = coordination_shell(pdb_file, inputs)
    held = coordination_shell(pdb_file, inputs, radius=constraint_cutoff,
                              per_residue=False)

    for m in metal_idx:
        print(f"\nMetal {elems[m]} at {chains[m]}:{resids[m]} (atom index {m})")
        qm_residue_dict.setdefault(chains[m], {})[resids[m]] = resnames[m]

        for row in [r for r in held if r["metal_index"] == m]:
            chain, resid_text, aname = row["donor"].split(":")
            resid = int(resid_text)
            if resid in (constraint_exclude.get(chain) or []):
                print(f"[INFO] No bond constraint for {row['resname']} {chain}:{resid} "
                      f"({aname}) (constraint_exclude_residues)")
            else:
                bondconstraints.append([m, row["donor_index"]])

        print("Coordinating residues:")
        for row in [r for r in shell if r["metal_index"] == m]:
            chain, resid_text, aname = row["donor"].split(":")
            resid, resn = int(resid_text), row["resname"]
            dist, idx = row["distance"], row["donor_index"]

            # Excluded residues keep their bond constraint but are left out of the QM region
            if exclude_residues and chain in exclude_residues:
                if resid in exclude_residues[chain]:
                    continue

            print(
                f"{resn:>3} {chain}:{resid:<4} ({aname:<4}) | "
                f"{dist:5.2f} Å | atom index {idx}"
            )

            # Build QM atoms (all non-backbone atoms of residue)
            for res in topology.residues():
                if (res.chain.id == chain and int(res.id) == resid):
                    qm_residue_dict.setdefault(chain, {})[resid] = resn

    return qm_residue_dict, bondconstraints


# ==============================
# QM Region
# ==============================

BACKBONE_ATOMS = {
    "N", "H", "H1", "H2", "H3",
    "CA", "HA", "HA2", "HA3",
    "C", "O", "OXT",
}

# Includes Amber protonation/disulfide variants so they are cut at CA-CB like
# standard residues instead of being put fully into the QM region.
AMINO_ACIDS = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
    "HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
    "THR","TRP","TYR","VAL",
    "HID","HIE","HIP","HSD","HSE","HSP","CYX","CYM",
    "ASH","GLH","LYN",
}


def build_qm_region(
    pdb_file: str,
    qm_residue_dict: Dict[str, Dict[int, str]],
    backbone_atoms=BACKBONE_ATOMS,
    amino_acids=AMINO_ACIDS,
):
    pdb = PDBFile(pdb_file)
    topology = pdb.topology

    qm_residue_dict = {
        chain: {str(k): v for k, v in residues.items()}
        for chain, residues in qm_residue_dict.items()
    }

    qm_atoms = set()
    boundary_excluded = set()

    print("\n===== QM Region Construction =====")

    for res in topology.residues():
        chain = res.chain.id
        resid = res.id.strip()

        if chain in qm_residue_dict and resid in qm_residue_dict[chain]:
            print(f"[QM] Including residue {res.name} {chain}:{resid}")

            for atom in res.atoms():
                if atom.name not in backbone_atoms:
                    qm_atoms.add(atom.index)

                if res.name not in amino_acids:
                    qm_atoms.add(atom.index)
                    boundary_excluded.add(atom.index)

    print(f"\n[INFO] QM atoms: {len(qm_atoms)}")
    print(f"[INFO] Boundary excluded atoms: {len(boundary_excluded)}")

    return sorted(qm_atoms), sorted(boundary_excluded)


def qm_region_settings(inputs: Dict[str, Any]) -> Dict[str, Any]:
    """backbone_atoms / amino_acids overrides from the preprocess inputs."""
    return {
        "backbone_atoms": set(inputs.get("backbone_atoms") or BACKBONE_ATOMS),
        "amino_acids": set(inputs.get("amino_acids") or AMINO_ACIDS),
    }


def check_residue_numbering(pdb_file: str) -> None:
    """
    ASH/pdbfixer renumber residues 1..N within each chain, and residue_variants,
    extra_qm_residues and the exclusion options are matched against those numbers.
    Require the input to already be numbered that way so the numbers mean the same thing.
    """
    bad = []
    for chain in PDBFile(pdb_file).topology.chains():
        ids = [int(r.id) for r in chain.residues()]
        if ids != list(range(1, len(ids) + 1)):
            bad.append(f"{chain.id} (starts at {ids[0]}, {len(ids)} residues)")
    if bad:
        raise ValueError(
            f"Residues in {pdb_file} are not numbered 1..N within each chain: {', '.join(bad)}. "
            "ASH renumbers them that way, so residue_variants / extra_qm_residues / exclusions would "
            "refer to different residues. Renumber the input (and the config) first.")


def check_his_variants(pdb_file: str, inputs: Dict[str, Any]) -> None:
    """
    List histidines whose tautomer is not set in `residue_variants` (pdbfixer would pick
    one by itself). With `require_his_variants: true` this is an error.
    """
    variants = inputs.get("residue_variants") or {}
    missing = []
    for res in PDBFile(pdb_file).topology.residues():
        if res.name in ("HIS", "HID", "HIE", "HIP") and int(res.id) not in (variants.get(res.chain.id) or {}):
            missing.append(f"{res.chain.id}:{res.id}")
    if missing:
        msg = (f"Histidine tautomer not set in residue_variants for {', '.join(missing)} "
               "(pdbfixer will choose one)")
        if inputs.get("require_his_variants"):
            raise ValueError(msg)
        print(f"[WARNING] {msg}")


def run(args) -> Dict[str, Any]:

    preprocess = load_section(args.config, "preprocess")
    inputs = preprocess.get("inputs")
    if inputs is None:
        raise KeyError(f"No preprocess.inputs section in {args.config}")

    # `default_forcefields` is the documented key (also read by md); `forcefields` is accepted for older configs
    forcefields = list(
        inputs.get("default_forcefields")
        or inputs.get("forcefields")
        or ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]
    )

    membrane = bool(inputs.get("membrane"))
    if membrane:
        if inputs.get("implicit"):
            raise ValueError("`membrane` and `implicit` cannot both be set: a bilayer needs explicit "
                             "lipids and water, while `implicit` adds no solvent at all.")
        # protein.ff14SB + tip3p have no lipid parameters, so addMembrane could not type
        # the lipids without this
        forcefields = lipid_forcefields(forcefields, inputs.get("membrane_lipid_xml") or LIPID_XML)

    resname = ligand_resname(inputs)

    if inputs.get("ligand"):
        print("[STEP] Preparing ligand...")
        pdb_file, ligand_xmlfile = merge_ligand(inputs.get("protein_pdb"), inputs.get("ligand"), inputs.get("model"),
                                                inputs.get("charge"), resname=resname,
                                                smiles=inputs.get("ligand_smiles"),
                                                merge=inputs.get("merge_ligand", True))
    else:
        ligand_xmlfile = inputs.get("ligand_xmlfile")
        pdb_file = inputs.get("protein_pdb")

    if ligand_xmlfile:
        forcefields.append(ligand_xmlfile)

    check_residue_numbering(pdb_file)
    check_his_variants(pdb_file, inputs)

    if not inputs.get("skip_solvate"):
        print("[STEP] Solvating system...")
        print(forcefields)
        forcefield_obj = ForceField(*forcefields)

        if inputs.get("implicit"):
            modeller, fragment = OpenMM_Modeller(
                pdbfile=pdb_file,
                implicit=True,
                forcefield_object=forcefield_obj,
                pH=inputs.get("pH", 7.0),
                residue_variants=inputs.get("residue_variants"),
                use_higher_occupancy=True,
            )
        elif membrane:
            # Explicit lipid bilayer plus water, through OpenMM's Modeller.addMembrane.
            # The membrane normal is +z and the bilayer is centred on `membrane_center_z`:
            # nothing here reorients the protein, so the transmembrane bundle must already
            # be aligned with z with its hydrophobic belt on that plane (e.g. from OPM).
            lipidtype = check_lipid_type(inputs.get("membrane_lipidtype"))
            print(f"[INFO] Building a {lipidtype} bilayer in the xy plane, centred on z = "
                  f"{inputs.get('membrane_center_z', 0.0)} Å. The input structure is assumed to be "
                  "aligned with the membrane normal along z already.")
            modeller, fragment = OpenMM_Modeller(
                pdbfile=pdb_file,
                forcefield_object=forcefield_obj,
                watermodel="tip3p",
                pH=inputs.get("pH", 8.0),
                membrane=True,
                membrane_lipidtype=lipidtype,
                # ASH passes this on as addMembrane's `minimumPadding`: the smallest
                # distance from any pre-existing atom to a box edge
                membrane_padding=inputs.get("membrane_padding", 10.0),
                membraneCenterZ=inputs.get("membrane_center_z", 0.0),
                ionicstrength=inputs.get("ionic_strength", 0.1),
                pos_iontype='Na+',
                neg_iontype='Cl-',
                residue_variants=inputs.get("residue_variants"),
                use_higher_occupancy=True,
            )
        else:
            modeller, fragment = OpenMM_Modeller(
                pdbfile=pdb_file,
                forcefield_object=forcefield_obj,
                watermodel="tip3p",
                pH=inputs.get("pH", 8.0),
                solvent_padding=inputs.get("solvent_padding", 10.0),
                ionicstrength=inputs.get("ionic_strength", 0.1),
                pos_iontype='Na+',
                neg_iontype='Cl-',
                residue_variants=inputs.get("residue_variants"),
                use_higher_occupancy=True,
            )

        pdb_file = Path("finalsystem.pdb").resolve().as_posix()

    # With skip_solvate, pdb_file stays as protein_pdb, or merged.pdb if a ligand was merged in

    if inputs.get("metals"):
        print("[STEP] Detecting coordination constraints...")
        qm_residues, constraints = find_coordinators(pdb_file, inputs)
        # A contact an earlier run restrained into place must not become a rigid
        # constraint just because the detector now measures it there.
        check_inherited_restraints(
            constraints,
            (load_section(args.config, "md").get("results") or {}).get("restraints"),
            pdb_file, inputs)
    else:
        # e.g. protein–protein complexes: no metal site to hold, no metal QM region
        print("[INFO] No metals given: no coordination constraints and no metal QM region")
        qm_residues, constraints = {}, []

    # Residues added to the QM region beyond the metal's coordination sphere
    extra = inputs.get("extra_qm_residues") or {}
    if extra:
        topology = PDBFile(pdb_file).topology
        for chain, resids in extra.items():
            for res in topology.residues():
                if res.chain.id == chain and int(res.id) in resids:
                    qm_residues.setdefault(chain, {})[int(res.id)] = res.name
                    print(f"[INFO] Extra QM residue {res.name} {chain}:{res.id}")
            missing = set(resids) - {int(r) for c, rs in qm_residues.items() if c == chain for r in rs}
            if missing:
                raise ValueError(f"extra_qm_residues: chain {chain} residues {sorted(missing)} not found in {pdb_file}")

    print("[STEP] Building QM region...")
    # Atom indices refer to processed_pdb; OpenMM keeps atom order through MD, so qmmm/neb reuse them
    qm_atoms, boundary_excluded = build_qm_region(pdb_file, qm_residues, **qm_region_settings(inputs))

    result = {
        "ligand_xml": ligand_xmlfile,
        "ligand_resname": resname if inputs.get("ligand") else None,
        "processed_pdb": Path(pdb_file).resolve().as_posix(),
        "qm_residues": qm_residues,
        "qm_atoms": qm_atoms,
        "boundary_excluded": boundary_excluded,
        "constraints": constraints,
    }
    save_result(result, filename=args.config, section="preprocess")

    return result
