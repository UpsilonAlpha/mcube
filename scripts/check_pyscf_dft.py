#!/usr/bin/env python3
"""
Check that PySCF can do DFT/MM with a given functional: functional recognised
(incl. VV10 non-local correlation), point-charge embedding, and an analytic
gradient that matches finite differences for both a QM atom and an MM charge.

Runs on a water molecule with three point charges, so it takes seconds.

Usage:
    python check_pyscf_dft.py [--xc wb97m_v] [--basis def2-svp] [--df]
"""

import argparse
import time

import numpy as np
import pyscf
from pyscf import dft, gto, qmmm
from pyscf.dft import libxc


def energy(mol_coords, mm_coords, mm_charges, args):
    mol = gto.M(atom=[("O", mol_coords[0]), ("H", mol_coords[1]), ("H", mol_coords[2])],
                basis=args.basis, unit="Angstrom", verbose=0)
    mf = dft.RKS(mol)
    mf.xc = args.xc
    if args.df:
        mf = mf.density_fit()
    mf = qmmm.mm_charge(mf, mm_coords, mm_charges, unit="Angstrom")
    mf.conv_tol = 1e-11
    e = mf.kernel()
    return mf, e


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--xc", default="wb97m_v")
    p.add_argument("--basis", default="def2-svp")
    p.add_argument("--df", action="store_true", help="Use density fitting")
    args = p.parse_args()

    print("pyscf", pyscf.__version__, "| libxc", libxc.__version__ if hasattr(libxc, "__version__") else "?")
    print(f"{args.xc}: hybrid={libxc.is_hybrid_xc(args.xc)}, range-separated={libxc.is_hybrid_xc(args.xc) and libxc.rsh_coeff(args.xc)[0] != 0}, "
          f"VV10 non-local={libxc.is_nlc(args.xc)}")
    try:
        import gpu4pyscf  # noqa: F401
        print("gpu4pyscf: available")
    except ImportError:
        print("gpu4pyscf: not installed")

    qm = np.array([[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]])
    mm = np.array([[0.0, 3.0, 0.5], [2.5, -1.0, 0.3], [-2.0, -2.0, -0.4]])
    q = np.array([-0.8, 0.4, 0.4])

    t = time.time()
    mf, e0 = energy(qm, mm, q, args)
    g = mf.nuc_grad_method()
    g_qm = g.kernel()
    dm = mf.make_rdm1()
    g_mm = g.grad_hcore_mm(dm) + g.grad_nuc_mm()  # Hartree/Bohr on MM charges
    print(f"E = {e0:.10f} Eh  ({time.time() - t:.1f} s)")

    h = 1e-3  # Angstrom
    bohr = 0.52917721092
    checks = []
    for label, which, i in (("QM O x", "qm", 0), ("MM q1 y", "mm", 0)):
        ax = 0 if which == "qm" else 1
        ep, em = [], []
        for sign, store in ((1, ep), (-1, em)):
            qq, mmm = qm.copy(), mm.copy()
            (qq if which == "qm" else mmm)[i, ax] += sign * h
            store.append(energy(qq, mmm, q, args)[1])
        fd = (ep[0] - em[0]) / (2 * h) * bohr  # Eh/Bohr
        an = (g_qm if which == "qm" else g_mm)[i, ax]
        checks.append((label, an, fd))
    for label, an, fd in checks:
        ok = abs(an - fd) < 2e-5
        print(f"{label}: analytic {an: .7f}  finite-diff {fd: .7f}  diff {an - fd: .1e}  {'OK' if ok else 'MISMATCH'}")


if __name__ == "__main__":
    main()
