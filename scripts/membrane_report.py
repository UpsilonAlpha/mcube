#!/usr/bin/env python3
"""
Membrane health and partner-pose report for a membrane MD trajectory (any system).

Works on the files `mcube md` writes for a membrane system (`md.results.trajfile` and
`md.results.lastframe`, or the NPT `relaxbox_NPT.dcd` with its last frame).

Always reported, per frame:
- box, gross area per lipid (the protein footprint is included, so it reads high),
  phosphate-plane thickness (P-P), leaflet imbalance, bilayer midplane z;
- lipid phosphates and water oxygens inside the receptor bundle (leak checks).

With a receptor (the first protein chain, or --receptor), also per frame, after
superposing on the receptor's transmembrane core:
- receptor backbone RMSD, plus --track (CA displacement, distance from the bundle axis
  and z), --distance and --angle readouts between named CA atoms;

With a partner (the second protein chain, or --partner: a peptide, binder, ...):
- partner CA RMSD, receptor-partner heavy-atom contacts (< --contact A) and the
  closest approach to --site residues;
- one-row summary (`_partner.tsv`): pose RMSF (in the receptor frame) and internal RMSF,
  CA RMSD and centre-of-mass drift, contact retention, site occupancy and, against a
  reference --footprint, the fraction of the reference contact set it reaches.
  `--write-footprint` on a reference run (a known binder) writes that contact set.

Residue numbers are those of the topology, which keeps the numbering of the input PDB
(mcube inputs are numbered 1..N per chain). Chain letters are those of the processed
PDB, not of the input (see the README: OpenMM relabels chains A, B, C, ... in order).

Usage:
    python membrane_report.py md/NVTtrajectory_lastframe_imaged.pdb md/NVTtrajectory_imaged.dcd \
        --out md/report [--frame-ps 50] [--receptor A] [--partner B] [--core-z 15] \
        [--core-exclude 290-330] [--site 203 345 349] [--footprint fp.txt | --write-footprint fp.txt] \
        [--track tm6=310] [--distance gate=319,322] [--angle kink=315,319,323]

Writes <out>.tsv (per frame), <out>.txt (summary), and with a partner <out>_partner.tsv,
<out>_residues.tsv (partner per-residue RMSF) and <out>_occupancy.tsv (receptor contacts).
"""

import argparse
import sys
from pathlib import Path

import mdtraj
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from mcube.membrane import LIPID_RESNAMES  # noqa: E402


def residue_range(text):
    """'290-330' -> (290, 330); '12' -> (12, 12)."""
    lo, _, hi = text.partition("-")
    return int(lo), int(hi or lo)


def named_residues(text, n):
    """'name=12,40' -> ('name', [12, 40]), with n residues required."""
    name, _, rest = text.partition("=")
    res = [int(x) for x in rest.split(",") if x]
    if not name or len(res) != n:
        raise argparse.ArgumentTypeError(f"expected NAME={','.join(['RESID'] * n)}, got {text!r}")
    return name, res


def protein_chains(top, receptor, partner):
    """(receptor, partner) mdtraj chains: by letter if given, else first and second protein chain."""
    chains = [c for c in top.chains if any(r.is_protein for r in c.residues)]

    def pick(letter, default):
        if letter is None:
            return default
        for c in chains:
            if getattr(c, "chain_id", None) == letter:
                return c
        index = ord(letter.upper()) - ord("A")  # OpenMM writes A, B, C, ... in chain order
        if 0 <= index < top.n_chains and top.chain(index) in chains:
            return top.chain(index)
        raise ValueError(f"no protein chain {letter!r} in the topology")

    rec = pick(receptor, chains[0] if chains else None)
    par = pick(partner, chains[1] if len(chains) > 1 else None)
    return rec, par


def rmsf(xyz):
    """RMSF per atom over frames; xyz is (frames, atoms, 3) in Å."""
    return np.sqrt(((xyz - xyz.mean(0)) ** 2).sum(-1).mean(0))


def rms(xyz, ref, idx):
    return np.sqrt(((xyz[:, idx] - ref[idx]) ** 2).sum(-1).mean(-1))


def angle(a, b, c):
    u, v = a - b, c - b
    cos = (u * v).sum(-1) / np.linalg.norm(u, axis=-1) / np.linalg.norm(v, axis=-1)
    return np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("topology")
    ap.add_argument("traj", nargs="+", help="one or more trajectories; segments are concatenated")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--frame-ps", type=float, default=0.0, help="ps between frames (0: frame number only)")
    ap.add_argument("--skip-frames", type=int, default=0, help="drop this many initial frames")
    ap.add_argument("--receptor", help="receptor chain letter (default: first protein chain)")
    ap.add_argument("--partner", help="partner chain letter (default: second protein chain, if any)")
    ap.add_argument("--core-z", type=float, default=15.0,
                    help="receptor CA within this many Å of the midplane form the superposition core")
    ap.add_argument("--core-exclude", action="append", default=[], type=residue_range, metavar="LO-HI",
                    help="receptor residues kept out of the core (a mobile helix); repeatable")
    ap.add_argument("--site", type=int, nargs="*", default=[], help="receptor residues the partner should reach")
    ap.add_argument("--contact", type=float, default=4.5, help="heavy-atom contact cutoff in Å")
    ap.add_argument("--footprint", help="reference contact set (receptor residue per line)")
    ap.add_argument("--write-footprint", help="write residues touched in >= 50%% of frames")
    ap.add_argument("--track", action="append", default=[], type=lambda s: named_residues(s, 1), metavar="NAME=RES")
    ap.add_argument("--distance", action="append", default=[], type=lambda s: named_residues(s, 2), metavar="NAME=R1,R2")
    ap.add_argument("--angle", action="append", default=[], type=lambda s: named_residues(s, 3), metavar="NAME=R1,R2,R3")
    a = ap.parse_args()

    first = mdtraj.load_frame(a.traj[0], 0, top=a.topology)
    top0 = first.topology
    lipid_p = np.array([x.index for x in top0.atoms
                        if x.residue.name in LIPID_RESNAMES and x.element is not None and x.element.symbol == "P"])
    water_o = np.array([x.index for x in top0.atoms if x.residue.is_water and x.element.symbol == "O"])
    if not len(lipid_p):
        sys.exit(f"No lipid phosphorus found (residue names {sorted(LIPID_RESNAMES)}): not a membrane system?")

    keep = np.unique(np.concatenate([lipid_p, water_o, top0.select("protein")]))
    traj = mdtraj.load(a.traj, top=a.topology, atom_indices=keep)[a.skip_frames:]
    mapping = {int(old): new for new, old in enumerate(keep)}
    top, n = traj.topology, traj.n_frames
    lipid_p = np.array([mapping[int(i)] for i in lipid_p])
    water_o = np.array([mapping[int(i)] for i in water_o])
    mid0 = float(np.median(traj.xyz[0, lipid_p, 2]))  # nm

    rec, par = protein_chains(top, a.receptor, a.partner)

    cols = {}
    xyz = traj.xyz * 10.0  # Å
    ref = xyz[0]
    axis_xy = None
    if rec is not None:
        res = {r.resSeq: r for r in rec.residues}
        ca = lambda s: res[s].atom("CA").index  # noqa: E731
        rec_res = list(rec.residues)
        rec_ca = np.array([r.atom("CA").index for r in rec_res])
        rec_seq = np.array([r.resSeq for r in rec_res])
        rec_bb = np.array([x.index for r in rec_res for x in r.atoms if x.name in ("N", "CA", "C")])
        near = np.abs(xyz[0, rec_ca, 2] - mid0 * 10.0) < a.core_z
        for lo, hi in a.core_exclude:
            near &= ~((rec_seq >= lo) & (rec_seq <= hi))
        core = rec_ca[near]
        traj.superpose(traj, frame=0, atom_indices=core)
        xyz = traj.xyz * 10.0
        ref = xyz[0]
        axis_xy = ref[core, :2].mean(axis=0)
        cols["bb_rmsd_core"] = rms(xyz, ref, core)
        cols["bb_rmsd_receptor"] = rms(xyz, ref, rec_bb)
        for name, (s,) in a.track:
            cols[f"{name}_disp"] = np.linalg.norm(xyz[:, ca(s)] - ref[ca(s)], axis=1)
            cols[f"{name}_radial"] = np.linalg.norm(xyz[:, ca(s), :2] - axis_xy, axis=1)
            cols[f"{name}_z"] = xyz[:, ca(s), 2]
        for name, (s, t) in a.distance:
            cols[f"{name}_d"] = np.linalg.norm(xyz[:, ca(s)] - xyz[:, ca(t)], axis=1)
        for name, (s, t, u) in a.angle:
            cols[f"{name}_deg"] = angle(xyz[:, ca(s)], xyz[:, ca(t)], xyz[:, ca(u)])

    # Membrane health
    lx, ly, lz = traj.unitcell_lengths.T * 10.0
    zp = xyz[:, lipid_p, 2]
    mid = np.median(zp, axis=1, keepdims=True)
    up = zp > mid
    cols["box_x"], cols["box_y"], cols["box_z"] = lx, ly, lz
    cols["apl_gross"] = lx * ly / (len(lipid_p) / 2.0)
    cols["thickness_PP"] = np.array([zp[f][up[f]].mean() - zp[f][~up[f]].mean() for f in range(n)])
    cols["leaflet_imbalance"] = up.sum(1) - (~up).sum(1)
    cols["midplane_z"] = mid[:, 0]
    if axis_xy is not None:
        pin = [((np.linalg.norm(xyz[f, lipid_p, :2] - axis_xy, axis=1) < 8.0)
                & (np.abs(zp[f] - mid[f, 0]) < a.core_z)).sum() for f in range(n)]
        win = [((np.linalg.norm(xyz[f, water_o, :2] - axis_xy, axis=1) < 6.0)
                & (np.abs(xyz[f, water_o, 2] - mid[f, 0]) < 8.0)).sum() for f in range(n)]
        cols["lipid_P_in_bundle"], cols["water_O_in_core"] = np.array(pin), np.array(win)

    # Partner pose
    partner_row = None
    if rec is not None and par is not None:
        p_res = list(par.residues)
        p_ca = np.array([r.atom("CA").index for r in p_res])
        p_heavy = np.array([x.index for r in p_res for x in r.atoms if x.element.symbol != "H"])
        r_atoms, r_start = [], []
        for r in rec_res:
            r_start.append(len(r_atoms))
            r_atoms += [x.index for x in r.atoms if x.element.symbol != "H"]
        r_heavy, r_start = np.array(r_atoms), np.array(r_start)
        site = [int(np.where(rec_seq == s)[0][0]) for s in a.site]
        site_ca = rec_ca[site] if site else np.array([], dtype=int)

        contacts, site_d = [], []
        touched = np.zeros((n, len(rec_res)), dtype=bool)
        for f in range(n):
            d = np.linalg.norm(xyz[f, p_heavy][:, None] - xyz[f, r_heavy][None], axis=-1)
            contacts.append(int((d < a.contact).sum()))
            touched[f] = np.minimum.reduceat(d.min(0), r_start) < a.contact
            if site:
                site_d.append(np.linalg.norm(xyz[f, p_heavy][:, None] - xyz[f, site_ca][None], axis=-1).min())
        contacts = np.array(contacts)
        cols["partner_ca_rmsd"] = rms(xyz, ref, p_ca)
        cols["partner_contacts"] = contacts
        if site:
            cols["partner_site_min_d"] = np.array(site_d)

        k = max(1, n // 5)
        occ = touched.mean(0)
        pose = rmsf(xyz[:, p_ca])
        aligned = mdtraj.Trajectory(traj.xyz.copy(), top)
        aligned.superpose(aligned, frame=0, atom_indices=p_ca)
        internal = rmsf(aligned.xyz[:, p_ca] * 10.0)
        cas = xyz[:, p_ca]
        rmsd = np.sqrt(((cas - cas[0]) ** 2).sum(-1).mean(-1))
        drift = np.linalg.norm(cas.mean(1) - cas[0].mean(0), axis=-1)

        if a.write_footprint:
            Path(a.write_footprint).write_text("\n".join(str(int(s)) for s in rec_seq[occ >= 0.5]) + "\n")
        fp_n, cov_mean, cov_last = 0, float("nan"), float("nan")
        if a.footprint:
            fp = {int(x) for x in Path(a.footprint).read_text().split()}
            sel = np.array([i for i, s in enumerate(rec_seq) if int(s) in fp])
            fp_n = len(sel)
            cov = touched[:, sel].mean(1)
            cov_mean, cov_last = float(cov.mean()), float(cov[-k:].mean())
        d0 = np.linalg.norm(xyz[0, rec_ca][:, None] - xyz[0, p_heavy][None], axis=-1).min(1)
        iface = rec_ca[d0 < 6.0]

        names = ["frames", "partner_len", "pose_rmsf_mean", "pose_rmsf_max", "internal_rmsf_mean",
                 "internal_rmsf_max", "ca_rmsd_mean", "ca_rmsd_last", "com_drift_mean", "com_drift_last",
                 "contacts_first", "contacts_mean", "contacts_retention", "site_d_mean", "site_d_max",
                 "site_occ", "footprint_n", "site_cov_mean", "site_cov_last", "iface_n", "iface_rmsf"]
        nan = float("nan")
        values = [n, len(p_ca), pose.mean(), pose.max(), internal.mean(), internal.max(), rmsd.mean(), rmsd[-1],
                  drift.mean(), drift[-1], contacts[0], contacts.mean(), contacts[-k:].mean() / max(contacts[0], 1),
                  np.mean(site_d) if site else nan, np.max(site_d) if site else nan,
                  occ[site].mean() if site else nan, fp_n, cov_mean, cov_last, len(iface),
                  rmsf(xyz[:, iface]).mean() if len(iface) else nan]
        fmt = lambda v: f"{v:.3f}" if isinstance(v, (float, np.floating)) else str(v)  # noqa: E731
        Path(a.out + "_partner.tsv").write_text("\t".join(names) + "\n" + "\t".join(fmt(v) for v in values) + "\n")
        with open(a.out + "_residues.tsv", "w") as fh:
            fh.write("resSeq\tresname\tpose_rmsf\tinternal_rmsf\n")
            for i, r in enumerate(p_res):
                fh.write(f"{r.resSeq}\t{r.name}\t{pose[i]:.3f}\t{internal[i]:.3f}\n")
        with open(a.out + "_occupancy.tsv", "w") as fh:
            fh.write("resSeq\tresname\toccupancy\n")
            for i, r in enumerate(rec_res):
                if occ[i] > 0:
                    fh.write(f"{r.resSeq}\t{r.name}\t{occ[i]:.3f}\n")
        partner_row = dict(zip(names, values))

    # Outputs
    time = np.arange(n) * a.frame_ps
    names = list(cols)
    with open(a.out + ".tsv", "w") as fh:
        fh.write("frame\ttime_ps\t" + "\t".join(names) + "\n")
        for f in range(n):
            fh.write(f"{f}\t{time[f]:.1f}\t" + "\t".join(f"{cols[k][f]:.3f}" for k in names) + "\n")
    k = max(1, n // 5)
    lines = [f"frames={n} span_ps={time[-1]:.0f} lipids={len(lipid_p)} waters={len(water_o)}",
             f"{'metric':22s} {'first':>9s} {'last':>9s} {'mean':>9s} {'sd':>8s} {'last20%mean':>12s}"]
    for name in names:
        v = np.asarray(cols[name], dtype=float)
        lines.append(f"{name:22s} {v[0]:9.2f} {v[-1]:9.2f} {v.mean():9.2f} {v.std():8.2f} {v[-k:].mean():12.2f}")
    if partner_row:
        lines.append("partner: " + " ".join(f"{k_}={v:.3g}" if isinstance(v, float) else f"{k_}={v}"
                                            for k_, v in partner_row.items()))
    Path(a.out + ".txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
