"""
Read-only probe of the ASH environment for the umbrella stage. Prints a report and changes
nothing. Run inside the ASH conda env, with the same PYTHONPATH a job uses:

    python docs/channel/probe_env.py > $L/probe_env.log 2>&1

Answers: which OpenMM and which platforms; which mcube is imported (git checkout or an
installed copy) and whether it has the umbrella stage and the 12-6-4 data; the stock
amber14/tip3p.xml radii of Ca2+, La3+ and Lu3+; and the optional packages.
"""

import hashlib
import importlib
import os
import re
import sys
from pathlib import Path

print("python", sys.version.split()[0], sys.executable)
print("PYTHONPATH", os.environ.get("PYTHONPATH"))

import openmm  # noqa: E402

names = [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())]
print("openmm", openmm.__version__, "platforms", names)

import mcube  # noqa: E402

print("mcube from", Path(mcube.__file__).parent)
for module in ("umbrella", "ions"):
    try:
        importlib.import_module(f"mcube.{module}")
        print(f"mcube.{module}: imports")
    except Exception as exc:  # report, do not stop
        print(f"mcube.{module}: FAILS {type(exc).__name__}: {exc}")
try:
    from mcube import ions
    print("ions sets", sorted(ions.SETS), "pol file", ions.POL_FILE, ions.POL_FILE.exists())
    if ions.POL_FILE.exists():
        print("pol sha256", hashlib.sha256(ions.POL_FILE.read_bytes()).hexdigest())
except Exception as exc:
    print("ions check FAILS", exc)

xml = Path(openmm.__file__).parent / "app" / "data" / "amber14" / "tip3p.xml"
text = xml.read_text()
for ion in ("Ca2+", "La3+", "Lu3+", "Gd3+"):
    m = re.search(r'<Atom [^>]*sigma="([0-9.]+)"[^>]*type="tip3p_standard-%s"' % re.escape(ion), text)
    if m:
        sigma = float(m.group(1))
        print(f"stock {ion}: sigma {sigma:.5f} nm, Rmin/2 {sigma * 2 ** (1 / 6) / 2 * 10:.3f} A")

for package in ("pymbar", "openmmplumed", "basis_set_exchange", "mdtraj", "pyscf"):
    try:
        mod = importlib.import_module(package)
        print(package, getattr(mod, "__version__", "present"))
    except Exception:
        print(package, "ABSENT")
