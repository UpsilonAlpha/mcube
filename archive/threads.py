"""
CPU thread placement for mcube runs.

On SMT nodes (e.g. Bunya's EPYC nodes), SLURM's `--cpus-per-task=N` gives N/2
physical cores plus their hyperthread siblings. Unbound OpenMP threads then float
and stack onto the same CPU; with one thread per physical core and binding, xTB
ran 2-9x faster and without the erratic spikes (docs/running-on-bunya.md).

Binding must only reach the xTB program, not the Python process: an in-process
OpenMP runtime (PySCF) bound with OMP_PROC_BIND pins the main thread to one core,
and threads created later (numpy/BLAS pools) inherit that single-core mask
(measured: a DFT single point ran ~18% CPU per thread on one core). So by default
mcube:
- leaves OMP_PROC_BIND/OMP_PLACES unset in its own process;
- runs xTB through a small wrapper (`xtb_dir()`) that exports OMP_PROC_BIND=close
  and OMP_PLACES=cores for the xTB process only;
- uses the number of physical cores in the process's CPU set when a stage's
  `cores` input is not given, and warns when `cores` asks for more.
If the user sets OMP_PROC_BIND or OMP_PLACES, or MCUBE_THREAD_BINDING=0, mcube
leaves binding to them and runs xTB directly. Stages that run several xTB
processes at once (NEB images, parallel NumFreq) call `unbind_threads()`, since
each would bind its threads to the same first cores.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Optional

_BIND_XTB = "MCUBE_BIND_XTB"  # "1" while mcube's xTB binding is active


def cpu_set() -> list:
    try:
        return sorted(os.sched_getaffinity(0))
    except AttributeError:  # not Linux
        return list(range(os.cpu_count() or 1))


def physical_cores() -> int:
    """Physical cores (distinct package/core ids) among the CPUs this process may use."""
    ids = set()
    for c in cpu_set():
        base = f"/sys/devices/system/cpu/cpu{c}/topology"
        try:
            with open(f"{base}/physical_package_id") as f:
                pkg = f.read().strip()
            with open(f"{base}/core_id") as f:
                core = f.read().strip()
            ids.add((pkg, core))
        except OSError:
            ids.add(("cpu", c))
    return max(1, len(ids))


def configure_threads(verbose: bool = True) -> int:
    """Decide the binding policy (see module docstring). Returns the physical core count."""
    cores = physical_cores()
    if os.environ.get("MCUBE_THREAD_BINDING", "1") == "0":
        os.environ[_BIND_XTB] = "0"
        msg = "disabled (MCUBE_THREAD_BINDING=0)"
    elif "OMP_PROC_BIND" in os.environ or "OMP_PLACES" in os.environ:
        os.environ[_BIND_XTB] = "0"
        msg = (f"left as set by the user: OMP_PROC_BIND={os.environ.get('OMP_PROC_BIND')}, "
               f"OMP_PLACES={os.environ.get('OMP_PLACES')}")
    else:
        os.environ[_BIND_XTB] = "1"
        msg = "xTB runs with OMP_PROC_BIND=close, OMP_PLACES=cores (Python process unbound)"
    if verbose:
        print(f"[INFO] CPUs available: {len(cpu_set())} ({cores} physical cores); thread binding: {msg}")
    return cores


def unbind_threads(reason: str):
    """Stop binding xTB (mcube's default only) before running parallel xTB processes."""
    if os.environ.get(_BIND_XTB) == "1":
        os.environ[_BIND_XTB] = "0"
        print(f"[INFO] xTB thread binding off: {reason}")


def xtb_dir(workdir: str = ".") -> Optional[str]:
    """
    Directory to pass to ASH's xTBTheory(xtbdir=...): a wrapper `xtb` that binds its
    OpenMP threads and then runs the real xTB. None (ASH finds xTB itself) when
    binding is off or xTB is not on PATH.
    """
    if os.environ.get(_BIND_XTB) != "1":
        return None
    real = shutil.which("xtb")
    if real is None or os.name == "nt":
        return None
    real = os.path.realpath(real)
    home = os.path.normpath(os.path.join(os.path.dirname(real), "..", "share", "xtb"))
    wrap_dir = Path(workdir, ".mcube_xtb").resolve()
    wrap_dir.mkdir(exist_ok=True)
    script = wrap_dir / "xtb"
    script.write_text(
        "#!/bin/sh\n"
        "# Written by mcube (threads.py): bind xTB's OpenMP threads, one per physical core\n"
        "export OMP_PROC_BIND=close\n"
        "export OMP_PLACES=cores\n"
        f"export XTBHOME={home}\n"
        f"export XTBPATH={home}${{XTBPATH:+:$XTBPATH}}\n"
        f'exec "{real}" "$@"\n')
    script.chmod(0o755)
    return str(wrap_dir)


def resolve_cores(requested: Optional[int], label: str = "cores") -> int:
    """`requested` if given (warning if it exceeds the physical cores), else the physical core count."""
    physical = physical_cores()
    if requested is None:
        print(f"[INFO] {label} not set: using {physical} (physical cores available)")
        return physical
    requested = int(requested)
    if requested > physical:
        print(f"[WARNING] {label} = {requested} exceeds the {physical} physical cores available: threads will "
              f"share cores (hyperthreads), which was slower for xTB. Consider {label}: {physical}.")
    return requested
