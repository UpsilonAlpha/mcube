#!/usr/bin/env python3
"""
Which basis sets in basis_set_exchange carry an element with a given ECP core size?
Used to find 4f-in-core (large-core) pseudopotentials for Lu3+ (core 60 = [Kr]4d10 4f14).

  python find_ecp.py --element Lu --core 60
"""

import argparse

import basis_set_exchange as bse

p = argparse.ArgumentParser()
p.add_argument("--element", default="Lu")
p.add_argument("--core", type=int, default=60)
args = p.parse_args()

z = bse.lut.element_Z_from_sym(args.element)
found = []
for name, meta in bse.get_metadata().items():
    latest = meta["versions"][meta["latest_version"]]
    if str(z) not in latest["elements"]:
        continue
    try:
        data = bse.get_basis(name, elements=[args.element])
    except Exception:
        continue
    el = data["elements"].get(str(z), {})
    core = el.get("ecp_electrons")
    if core:
        shells = el.get("electron_shells", [])
        ams = sorted({am for s in shells for am in s["angular_momentum"]})
        found.append((core, meta["display_name"], name, ams, len(shells)))
for core, display, name, ams, nsh in sorted(found):
    flag = "  <==" if core == args.core else ""
    print(f"core {core:3d}  {display:40s} [{name}]  l={ams} shells={nsh}{flag}")
