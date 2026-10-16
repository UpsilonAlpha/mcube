"""
CPU thread placement for mcube runs.

On SMT nodes (e.g. Bunya's EPYC nodes), SLURM's `--cpus-per-task=N` gives N/2
physical cores plus their hyperthread siblings. Unbound OpenMP threads then float
and stack onto the same CPU; with one thread per physical core and binding, xTB
ran 2-9x faster and without the erratic spikes (docs/running-on-bunya.md).

So by default mcube:
- sets OMP_PROC_BIND=close and OMP_PLACES=cores before ASH, OpenMM or xTB start
  (unless they are already set, or MCUBE_THREAD_BINDING=0);
- uses the number of physical cores in the process's CPU set when a stage's
  `cores` input is not given, and warns when `cores` asks for more.
Stages that run several worker processes at once (NEB images, parallel NumFreq)
call `unbind_threads()`: each worker's OpenMP runtime would otherwise bind its
threads to the same first cores.
"""

from __future__ import annotations

import os
from typing import Optional

_SET_BY_MCUBE = "MCUBE_SET_THREAD_BINDING"


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
    """Set the OpenMP binding defaults (see module docstring). Returns the physical core count."""
    cores = physical_cores()
    if os.environ.get("MCUBE_THREAD_BINDING", "1") == "0":
        msg = "disabled (MCUBE_THREAD_BINDING=0)"
    elif "OMP_PROC_BIND" in os.environ or "OMP_PLACES" in os.environ:
        msg = (f"left as set: OMP_PROC_BIND={os.environ.get('OMP_PROC_BIND')}, "
               f"OMP_PLACES={os.environ.get('OMP_PLACES')}")
    else:
        os.environ["OMP_PROC_BIND"] = "close"
        os.environ["OMP_PLACES"] = "cores"
        os.environ[_SET_BY_MCUBE] = "1"
        msg = "OMP_PROC_BIND=close, OMP_PLACES=cores"
    if verbose:
        print(f"[INFO] CPUs available: {len(cpu_set())} ({cores} physical cores); thread binding: {msg}")
    return cores


def unbind_threads(reason: str):
    """Undo mcube's binding defaults (not the user's) before running parallel worker processes."""
    if os.environ.pop(_SET_BY_MCUBE, None):
        os.environ["OMP_PROC_BIND"] = "false"
        os.environ.pop("OMP_PLACES", None)
        print(f"[INFO] Thread binding off: {reason}")


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
