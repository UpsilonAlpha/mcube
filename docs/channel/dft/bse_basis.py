#!/usr/bin/env python3
"""
Can PySCF get def2-SVP (and its ECP) for every element we need through basis_set_exchange?

PySCF 2.13's bundled def2 files stop at La, and it does not fall back to BSE by itself for
an element they lack. This writes one NWChem file from BSE for the elements given and loads
it the way mcube does for `dft_basis_file` / `dft_ecp_file` (pyscf.gto.basis.load and
load_ecp per element), then prints the AO count, electron count and ECP core per ion.

  python bse_basis.py --out def2-svp_HOLuLaCa.nw
"""

import argparse
from pathlib import Path

import basis_set_exchange as bse
from pyscf import gto

p = argparse.ArgumentParser()
p.add_argument("--basis", default="def2-SVP")
p.add_argument("--elements", default="H,O,C,N,S,Ca,La,Lu")
p.add_argument("--out", required=True)
args = p.parse_args()

elements = args.elements.split(",")
text = bse.get_basis(args.basis, elements=elements, fmt="nwchem", header=False)
path = Path(args.out).resolve()
path.write_text(text)
print("wrote", path)

for symbol, charge in (("La", 3), ("Lu", 3), ("Ca", 2)):
    basis = gto.basis.load(str(path), symbol)
    try:
        ecp = gto.basis.load_ecp(str(path), symbol)
    except Exception:
        ecp = None
    mol = gto.M(atom=f"{symbol} 0 0 0", charge=charge, spin=0, basis={symbol: basis},
                ecp={symbol: ecp} if ecp else {})
    print(symbol, "AOs", mol.nao, "electrons", mol.nelectron,
          "ECP core", mol.atom_nelec_core(0))
