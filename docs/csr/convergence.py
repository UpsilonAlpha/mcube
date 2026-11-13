#!/usr/bin/env python
"""Is 10 ns different from 100 ns? Answer it from trajectories we already have.

The question (user, 2026-11-06): before spending ~7 h of H100 on a fresh 100 ns
production run, check whether the structure is still changing between 10 ns and
100 ns. If a 10 ns window already gives the same RMSF profile and the same
thermodynamic state as the full 100 ns, the equilibration we have is enough and
the production run buys nothing.

What it measures, per run:

  RMSD(t)        protein CA, superposed on frame 0. Reported as block means, so
                 a plateau is visible rather than inferred from a single number.
  RMSF           per-residue fluctuation about the window mean, computed on
                 several windows of the SAME trajectory:
                   0-10 ns, 0-100 ns, 10-100 ns, and two independent halves.
                 Two numbers decide convergence - the Pearson r between a
                 window's profile and the full one (does 10 ns rank the mobile
                 residues the same way?) and the ratio of their mean RMSF (does
                 10 ns UNDERSTATE the motion, which is the usual failure?).
  site RMSF      the same, restricted to the residues this project argues about,
                 because a global profile can converge while the site has not.
  T, P, E        from `nvtsim.csv` / `nptsim.csv`: mean, SD and drift (last
                 decile minus first decile) over the same windows. An NVT run
                 has no meaningful pressure; it is read from the NPT phase.

Usage:  python convergence.py RUN_DIR [RUN_DIR ...]

Each RUN_DIR is an `md_*_results` directory. Everything is printed; nothing is
written, so this cannot damage a result or consume quota.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

# The md stage writes its topology as .cif for some runs and .pdb for others,
# and mcube's own MD_TRAJECTORIES constant lists only the .pdb spellings - which
# is why `salvage_results` cannot find a trajectory pair for the wt100/prod100
# runs at all. Try every combination rather than assume one (¶225).
TRAJ_PAIRS = tuple(
    (traj, top)
    for traj in ("NVTtrajectory_imaged.dcd", "NVTtrajectory.dcd")
    for top in ("NVTtrajectory_firstframe.pdb", "NVTtrajectory_firstframe.cif",
                "NVTtrajectory_lastframe.pdb", "NVTtrajectory_lastframe.cif")
)

# The residues this project's argument turns on, in project numbering.
SITE_RESIDUES = (88, 90, 92, 93, 164, 186, 225, 226)


def find_traj(d: Path):
    for traj, top in TRAJ_PAIRS:
        if (d / traj).is_file() and (d / top).is_file():
            return d / traj, d / top
    return None, None


def pearson(a, b):
    n = len(a)
    if n < 2:
        return float("nan")
    ma, mb = sum(a) / n, sum(b) / n
    va = sum((x - ma) ** 2 for x in a)
    vb = sum((x - mb) ** 2 for x in b)
    if va <= 0 or vb <= 0:
        return float("nan")
    cov = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    return cov / math.sqrt(va * vb)


def rmsf(traj, atom_idx):
    """Per-atom RMSF about the mean structure of the frames given, in A."""
    import numpy as np
    xyz = traj.xyz[:, atom_idx, :]              # nm
    mean = xyz.mean(axis=0)
    dev = ((xyz - mean) ** 2).sum(axis=2)       # nm^2 per frame per atom
    return np.sqrt(dev.mean(axis=0)) * 10.0     # A


def csv_stats(path: Path, keys=("Temperature (K)", "Box Volume (nm^3)",
                                "Total Energy (kJ/mole)", "Density (g/mL)")):
    if not path.is_file():
        return None
    rows = list(csv.DictReader(open(path)))
    if not rows:
        return None
    out = {"n": len(rows)}
    header = rows[0].keys()
    for want in keys:
        col = next((h for h in header if h and want.lower() in h.lower()), None)
        if not col:
            continue
        vals = []
        for r in rows:
            try:
                vals.append(float(r[col]))
            except (TypeError, ValueError):
                pass
        if len(vals) < 10:
            continue
        n = len(vals)
        mean = sum(vals) / n
        sd = math.sqrt(sum((v - mean) ** 2 for v in vals) / n)
        dec = max(1, n // 10)
        drift = (sum(vals[-dec:]) / dec) - (sum(vals[:dec]) / dec)
        out[want] = (mean, sd, drift)
    return out


def report(d: Path):
    import mdtraj as mdt

    print("=" * 78)
    print(d.as_posix())
    traj_f, top_f = find_traj(d)
    if not traj_f:
        print("  no NVT trajectory/topology pair here")
        return

    t = mdt.load(traj_f.as_posix(), top=top_f.as_posix())
    ca = t.topology.select("name CA and protein")
    if len(ca) == 0:
        print("  no protein CA atoms")
        return
    print(f"  {t.n_frames} frames, {t.n_atoms} atoms, {len(ca)} CA")

    # --- RMSD vs frame 0, as block means so a plateau shows -----------------
    t.superpose(t, frame=0, atom_indices=ca)
    rms = mdt.rmsd(t, t, frame=0, atom_indices=ca) * 10.0
    nb = min(10, t.n_frames)
    size = max(1, t.n_frames // nb)
    print("  CA RMSD vs frame 0, by decile of the run (A):")
    blocks = []
    for i in range(0, t.n_frames - size + 1, size):
        seg = rms[i:i + size]
        blocks.append(sum(seg) / len(seg))
    print("    " + "  ".join(f"{b:5.2f}" for b in blocks))
    print(f"    first {blocks[0]:.2f} -> last {blocks[-1]:.2f}"
          f"   (change {blocks[-1] - blocks[0]:+.2f})")

    # --- RMSF over windows of the same trajectory ---------------------------
    # Windows are expressed as frame fractions so this works whatever the
    # trajectory length or sampling interval is.
    full = rmsf(t, ca)
    site_sel = t.topology.select(
        "name CA and resid " + " ".join(str(r - 1) for r in SITE_RESIDUES))
    windows = {
        "first 10%": (0.0, 0.10),
        "first 50%": (0.0, 0.50),
        "last  50%": (0.50, 1.0),
        "last  90%": (0.10, 1.0),
    }
    print("  RMSF convergence (CA, A):   mean   vs-full r   mean ratio")
    mean_full = sum(full) / len(full)
    print(f"    full trajectory          {mean_full:6.2f}      -           -")
    for label, (lo, hi) in windows.items():
        i0, i1 = int(lo * t.n_frames), max(int(hi * t.n_frames), 2)
        if i1 - i0 < 2:
            continue
        sub = t[i0:i1]
        w = rmsf(sub, ca)
        mw = sum(w) / len(w)
        print(f"    {label}                {mw:6.2f}   {pearson(list(w), list(full)):8.3f}"
              f"    {mw / mean_full:8.3f}")
    if len(site_sel):
        fs = rmsf(t, site_sel)
        i1 = max(int(0.10 * t.n_frames), 2)
        ws = rmsf(t[0:i1], site_sel)
        print("  site residues (CA RMSF, A): " +
              "  ".join(f"{r}:{v:.2f}" for r, v in zip(SITE_RESIDUES, fs)))
        print("    same, first 10% only:     " +
              "  ".join(f"{r}:{v:.2f}" for r, v in zip(SITE_RESIDUES, ws)))

    # --- thermodynamic state ------------------------------------------------
    for name in ("nvtsim.csv", "nptsim.csv"):
        st = csv_stats(d / name)
        if not st:
            continue
        print(f"  {name}: {st['n']} rows")
        for k, v in st.items():
            if k == "n":
                continue
            mean, sd, drift = v
            print(f"    {k:28s} mean {mean:12.3f}  sd {sd:9.3f}  drift {drift:+.3f}")


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    for a in argv[1:]:
        d = Path(a)
        if not d.is_dir():
            print(f"{a}: not a directory")
            continue
        try:
            report(d)
        except Exception as exc:                     # keep going across runs
            print(f"  FAILED: {type(exc).__name__}: {exc}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
