#!/usr/bin/env python3
"""
Write one NWChem basis+ECP file from basis_set_exchange for the elements the channel's DFT
needs, and check it loads the way mcube does for `dft_basis_file` / `dft_ecp_file`
(pyscf.gto.basis.load and load_ecp per element): AO count, electrons and ECP core per ion.

PySCF 2.13's bundled def2 files stop at La, and it does not fall back to BSE by itself.
--override puts some elements on another basis family, e.g. the lanthanides on the
4f-in-core lcecp-1 sets (core 46 for La3+, 60 for Lu3+), so both have their 4f shell
treated the same way (decision log C19).

  python bse_basis.py --out def2-svp.nw
  python bse_basis.py --override La,Lu=lcecp-1-SVP --out def2-svp_lcecp1.nw
"""

import argparse
from pathlib import Path

import basis_set_exchange as bse
from pyscf import gto

p = argparse.ArgumentParser()
p.add_argument("--basis", default="def2-SVP")
p.add_argument("--elements", default="H,O,C,N,S,Ca,La,Lu")
p.add_argument("--override", action="append", default=[],
               help="ELEMENTS=BASIS, e.g. La,Lu=lcecp-1-SVP (repeatable)")
p.add_argument("--out", required=True)
args = p.parse_args()

elements = args.elements.split(",")
family = {e: args.basis for e in elements}
for spec in args.override:
    els, name = spec.split("=", 1)
    for e in els.split(","):
        family[e] = name

blocks = []
for name in dict.fromkeys(family.values()):
    group = [e for e in elements if family[e] == name]
    blocks.append(bse.get_basis(name, elements=group, fmt="nwchem", header=False))
path = Path(args.out).resolve()
path.write_text("\n".join(blocks))
print("wrote", path, {e: family[e] for e in elements})

for symbol, charge in (("La", 3), ("Lu", 3), ("Ca", 2)):
    if symbol not in elements:
        continue
    basis = gto.basis.load(str(path), symbol)
    try:
        ecp = gto.basis.load_ecp(str(path), symbol)
    except Exception:
        ecp = None
    mol = gto.M(atom=f"{symbol} 0 0 0", charge=charge, spin=0, basis={symbol: basis},
                ecp={symbol: ecp} if ecp else {})
    print(symbol, family[symbol], "AOs", mol.nao, "electrons", mol.nelectron,
          "ECP core", mol.atom_nelec_core(0))
