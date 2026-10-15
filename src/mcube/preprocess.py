#!/usr/bin/env python3
"""
Protein/Ligand preprocessing utilities (pipeline-ready).
"""

from __future__ import annotations

import argparse
import yaml
import os
import shutil
from typing import Dict, Any, Tuple
from pathlib import Path
from datetime import datetime

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds
from openmm.app import PDBFile, Topology, ForceField

from ash import (
    Fragment,
    OpenMM_Modeller,
    small_molecule_parameterizer,
    merge_pdb_files,
)

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
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(".tmp")

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
    cutoff = inputs.get("cutoff")

    coords = np.array([[p.x, p.y, p.z] for p in pdb.positions]) * 10.0  # nm → Å

    elems = [a.element.symbol if a.element else "" for a in atoms]
    names = [a.name for a in atoms]

    # Chain-aware residue info
    resids = [int(a.residue.id) for a in atoms]
    resnames = [a.residue.name for a in atoms]
    chains = [a.residue.chain.id for a in atoms]

    metal_idx = [i for i, e in enumerate(elems) if e.upper() in inputs.get("metals")]

    donor_elements = {e.upper() for e in inputs.get("donor_elements") or ["N", "O", "S"]}
    if inputs.get("donor_atoms"):
        print("[WARNING] donor_atoms is no longer used; donors are chosen by element (donor_elements)")

    qm_residue_dict = {}  # OpenMM-style: {"A": {80: "HIS"}}
    bondconstraints = []

    print("\n===== Metal Coordination Analysis =====")
    print(f"Cutoff: {cutoff:.2f} Å")

    for m in metal_idx:
        print(f"\nMetal {elems[m]} at {chains[m]}:{resids[m]} (atom index {m})")
        qm_residue_dict.setdefault(chains[m], {})[resids[m]] = resnames[m]

        residue_min_dist = {}

        for i in range(len(atoms)):
            if i == m:
                continue

            # donor atom filter
            if elems[i].upper() not in donor_elements:
                continue

            dist = np.linalg.norm(coords[i] - coords[m])
            if dist > cutoff:
                continue

            key = (chains[i], resids[i], resnames[i])

            if key not in residue_min_dist or dist < residue_min_dist[key][0]:
                residue_min_dist[key] = (dist, i, names[i])

        sorted_res = sorted(residue_min_dist.items(), key=lambda x: x[1][0])

        print("Coordinating residues:")
        for (chain, resid, resn), (dist, idx, aname) in sorted_res:
            
            bondconstraints.append([m, idx])

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

    print("[STEP] Detecting coordination constraints...")
    qm_residues, constraints = find_coordinators(pdb_file, inputs)

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
