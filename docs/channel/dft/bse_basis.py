#!/usr/bin/env python3
"""
Can PySCF get def2-SVP (and its ECP) for every element we need through basis_set_exchange?

PySCF 2.13's bundled def2 files stop at La, and it does not fall back to BSE by itself for
an element they lack, so the basis is fetched from BSE explicitly here and parsed into a
PySCF molecule. Prints the AO count and electron count per ion as a check.

  python bse_basis.py
"""

import basis_set_exchange as bse
from pyscf import gto


def bse_basis(name: str, element: str):
    """(basis, ecp) for one element from BSE, parsed for PySCF (ecp None if all-electron)."""
    text = bse.get_basis(name, elements=[element], fmt="nwchem", header=False)
    basis = gto.basis.parse(text)
    try:
        ecp = gto.basis.parse_ecp(text)
    except Exception:
        ecp = None
    return basis, ecp or None


for symbol, charge in (("La", 3), ("Lu", 3), ("Ca", 2)):
    basis, ecp = bse_basis("def2-SVP", symbol)
    mol = gto.M(atom=f"{symbol} 0 0 0", charge=charge, spin=0,
                basis={symbol: basis}, ecp={symbol: ecp} if ecp else {})
    print(symbol, "AOs", mol.nao, "electrons", mol.nelectron, "ECP core",
          mol.atom_nelec_core(0) if ecp else 0)
