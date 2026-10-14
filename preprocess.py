#!/usr/bin/env python3
"""
Protein/Ligand preprocessing utilities (pipeline-ready).
"""

from __future__ import annotations

import argparse
import yaml
import os
import shutil
from typing import Dict, Any
from pathlib import Path
from datetime import datetime

import numpy as np
from rdkit import Chem
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

def merge_ligand(protein_pdb: str, ligand: str, model: int = 1, charge: int = 0) -> Dict[str, Any]:
    filename = Path(ligand).stem
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

        small_molecule_parameterizer(
            forcefield_option="OpenFF",
            molfile=f"{filename}.mol",
            charge=charge,
            resname=filename
        )
        ligand_xmlfile = Path(f"openff_{filename}.xml").resolve().as_posix()

        merge_pdb_files(protein_pdb, f"{filename}.pdb")
        protein_pdb = "merged.pdb"

    elif ligand.endswith(".pdb"):
        filename = Path(ligand).stem
        print("Transforming PDB to .mol file")
        mol = Chem.MolFromPDBFile(ligand, removeHs=False)
        #mol = Chem.AddHs(mol, addCoords=True)
        Chem.MolToMolFile(mol, f"{filename}.mol", includeStereo=True)

        small_molecule_parameterizer(forcefield_option="OpenFF", pdbfile=f"{filename}.mol", charge=charge, resname=filename)
        merge_pdb_files(protein_pdb, f"{filename}.pdb")
        protein_pdb = "merged.pdb"
        ligand_xmlfile = Path(f"openff_{filename}.xml").resolve().as_posix()

    elif ligand.endswith(".mol"):
        filename = Path(ligand).stem
        
        #mol = Chem.MolFromMolFile(ligand, removeHs=False, sanitize=False)
        #Chem.MolToMolFile(mol, f"{filename}.mol")
        small_molecule_parameterizer(forcefield_option="OpenFF", molfile=f"{filename}.mol", charge=charge, resname=filename)
        merge_pdb_files(protein_pdb, f"{filename}.pdb")
        protein_pdb = "merged.pdb"
        ligand_xmlfile = Path(f"openff_{filename}.xml").resolve().as_posix()

    else:
        print("Not a valid file format for ligand: Choose either .pdbqt or .pdb or .mol")
    
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

    qm_atoms = set()
    boundary_excluded = set()
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
            if all(d not in names[i] for d in inputs.get("donor_atoms")):
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

            # Skip excluded residues
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





def run(args) -> Dict[str, Any]:

    preprocess = load_section(args.config, "preprocess")
    inputs = preprocess.get("inputs")

    if inputs.get("forcefields"):
        forcefields = inputs.get("forcefields")
    else:
        forcefields = ["amber14/protein.ff14SB.xml, amber14/tip3p.xml"]

    if inputs.get("ligand"):
        print("[STEP] Preparing ligand...")
        pdb_file, ligand_xmlfile = merge_ligand(inputs.get("protein_pdb"), inputs.get("ligand"), inputs.get("model"), inputs.get("charge"))
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
                residue_variants=inputs.get("residue_variants"),
                use_higher_occupancy=True,
            )
        else:
            modeller, fragment = OpenMM_Modeller(
                pdbfile=pdb_file,
                forcefield_object=forcefield_obj,
                watermodel="tip3p",
                pH=8.0,
                solvent_padding=10.0,
                ionicstrength=0.1,
                pos_iontype='Na+',
                neg_iontype='Cl-',
                residue_variants=inputs.get("residue_variants"),
                use_higher_occupancy=True,
            )
        
        pdb_file = Path("finalsystem.pdb").resolve().as_posix()
            
    else:
        pdb_file = inputs.get("protein_pdb")

    print("[STEP] Detecting coordination constraints...")
    qm_residues, constraints = find_coordinators(pdb_file, inputs)
    
    save_result(
        {
            "ligand_xml": ligand_xmlfile,
            "processed_pdb": Path(pdb_file).resolve().as_posix(),
            "qm_residues": qm_residues,
            "constraints": constraints
        },
        filename=args.config,
        section="preprocess"
    )    

    return
