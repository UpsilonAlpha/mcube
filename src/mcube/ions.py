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

Sources, as given per row in the file itself:
  +2 ions   Li, Roberts, Chakravorty & Merz, J. Chem. Theory Comput. 2013, 9, 2733
  +3/+4     Li, Song & Merz, J. Phys. Chem. B 2015, 119, 883

The 12-6-4 model does better for highly charged ions but needs an r^-4 term
with per-ion C4 values and the partner atom's polarisability, which Amber
supplies through ParmEd rather than in the frcmod file. Not implemented here.
"""

from __future__ import annotations

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

SETS = {"iod": IOD_TIP3P}


def lookup(name: str, element: str, charge: float):
    """(sigma in nm, epsilon in kJ/mol) for this ion, or None if not tabulated."""
    table = SETS.get(str(name).lower())
    if table is None:
        raise ValueError(f"unknown metal parameter set {name!r}; available: "
                         f"{', '.join(sorted(SETS))}")
    entry = table.get(element.upper(), {}).get(int(round(charge)))
    if entry is None:
        return None
    rmin_half, eps_kcal = entry
    return (2.0 * rmin_half / 2.0 ** (1.0 / 6.0) / 10.0, eps_kcal * 4.184)
