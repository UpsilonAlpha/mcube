"""How long an ion's protein contacts last in one long umbrella window (decision log C21).

    python residence.py LA/umbrella_res_035 [--cutoff 3.0] [--offset 130]

Reads the window's start.pdb (topology), window.dcd (frames) and cv.csv (cn_protein every CV
step). For every protein oxygen within `cutoff` of the ion in any frame it reports, per residue
(channel numbering = mcube + offset): the fraction of frames in contact, the number of contact
episodes, and their mean and longest length. The censored first/last episodes are flagged. When
cv.csv is complete it also gives z and the coordination numbers by half; a running window (cv.csv
still buffered) is read from its frames alone. Prints a JSON summary and writes residence.json in
the window directory.
"""
import argparse
import json
from pathlib import Path

import numpy as np


def episodes(on):
    """(start, length) of each run of True in a boolean series."""
    out, start = [], None
    for i, v in enumerate(on):
        if v and start is None:
            start = i
        elif not v and start is not None:
            out.append((start, i - start))
            start = None
    if start is not None:
        out.append((start, len(on) - start))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("window")
    ap.add_argument("--cutoff", type=float, default=3.0)
    ap.add_argument("--offset", type=int, default=130)
    ap.add_argument("--frame-ps", type=float, default=25.0,
                    help="ps between frames when cv.csv is not written yet (traj_frequency x timestep)")
    ap.add_argument("--ion", default=None, help="ion residue name (default: LA/LU/CA found in the topology)")
    a = ap.parse_args()
    from mdtraj.formats import DCDTrajectoryFile
    from openmm.app import PDBFile
    from mcube.preprocess import AMINO_ACIDS, pdb_resseq

    win = Path(a.window)
    # OpenMM writes residue numbers past 9999 in its own hex form: read the topology with it
    atoms = list(PDBFile(str(win / "start.pdb")).topology.atoms())
    names = [a.ion] if a.ion else ["LA", "LU", "CA"]
    ion = [at.index for at in atoms if at.residue.name in names and len(list(at.residue.atoms())) == 1]
    if len(ion) != 1:
        raise SystemExit(f"expected one ion, found {len(ion)}")
    ion = ion[0]
    oxy = [at.index for at in atoms
           if at.element is not None and at.element.symbol == "O" and at.residue.name in AMINO_ACIDS]
    with DCDTrajectoryFile(str(win / "window.dcd")) as fh:
        xyz, cell, _ = fh.read()                                         # A
    d = xyz[:, oxy, :] - xyz[:, [ion], :]
    d -= cell[:, None, :] * np.round(d / cell[:, None, :])               # orthorhombic minimum image
    d = np.linalg.norm(d, axis=2)
    n_frames = len(xyz)
    # cv.csv is buffered until the window ends; a running window is read from its frames alone
    cv_file = win / "cv.csv"
    cv = (np.genfromtxt(cv_file, delimiter=",", names=True)
          if cv_file.exists() and cv_file.stat().st_size > 100 else None)
    dt = float(cv["time_ps"][-1]) / n_frames if cv is not None else a.frame_ps

    by_res = {}
    for j, o in enumerate(oxy):
        if not (d[:, j] < a.cutoff).any():
            continue
        r = atoms[o].residue
        key = f"{r.name}{pdb_resseq(r.id) + a.offset}"
        by_res.setdefault(key, np.zeros(n_frames, bool))
        by_res[key] |= d[:, j] < a.cutoff

    residues = {}
    for key, on in sorted(by_res.items(), key=lambda kv: -kv[1].mean()):
        ep = episodes(on)
        lengths = np.array([n for _, n in ep]) * dt
        residues[key] = {
            "fraction": round(float(on.mean()), 3),
            "episodes": len(ep),
            "mean_ps": round(float(lengths.mean()), 1),
            "longest_ps": round(float(lengths.max()), 1),
            "bound_at_start": bool(on[0]),
            "bound_at_end": bool(on[-1]),
        }

    summary = {"window": win.resolve().as_posix(), "frames": n_frames, "frame_ps": round(dt, 2),
               "time_ps": round(n_frames * dt, 1), "running": cv is None}
    if cv is not None:
        half = len(cv) // 2
        summary.update({
            "z_mean": round(float(cv["z_A"].mean()), 2), "z_std": round(float(cv["z_A"].std()), 2),
            "cn_protein_first_half": round(float(cv["cn_protein"][:half].mean()), 2),
            "cn_protein_second_half": round(float(cv["cn_protein"][half:].mean()), 2),
            "cn_water_mean": round(float(cv["cn_water"].mean()), 2),
        })
    summary["residues"] = residues
    (win / "residence.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
