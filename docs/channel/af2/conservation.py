#!/usr/bin/env python3
"""
How conserved is each residue of the channel among its relatives, from AF2's own MSAs?

Reads the alignments AF2 wrote (msas/uniref90_hits.sto, mgnify_hits.sto, bfd_uniref_hits.a3m),
keeps the query's match columns, weights the sequences (1 / number of neighbours at >= 80 %
identity), and per query position reports:

  coverage   weighted fraction of relatives with a residue there (not a gap): does the
             family even have this segment?
  identity   weighted fraction with the query's own residue
  conservation  1 - Shannon entropy / log(20) over the 20 amino acids at that column
  lumen      the side chain points into the barrel (Cβ nearer the axis than Cα), from the model

and the same averaged over the N-terminal segment, the strands, and the loops, so the
segment can be compared with the barrel it sits in.

  python conservation.py --msas monomer/mature/msas --pdb monomer/mature/ranked_0.pdb --out cons.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"


def read_a3m(path):
    seqs, name = [], None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            seqs.append("")
        elif seqs:
            seqs[-1] += "".join(c for c in line.strip() if not c.islower() and c != ".")
    return seqs


def read_sto(path):
    rows, order = {}, []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("//"):
            continue
        name, seq = line.split()[:2]
        if name not in rows:
            rows[name] = ""
            order.append(name)
        rows[name] += seq
    query = rows[order[0]]
    keep = [i for i, c in enumerate(query) if c not in "-."]
    return ["".join(rows[n][i] for i in keep).upper().replace(".", "-") for n in order]


def weights(msa):
    x = np.array([list(s) for s in msa])
    ident = np.zeros((len(msa), len(msa)))
    for i in range(len(msa)):
        both = (x != "-") & (x[i] != "-")
        same = (x == x[i]) & both
        ident[i] = same.sum(1) / np.maximum(both.sum(1), 1)
    return 1.0 / (ident >= 0.8).sum(1)


def lumen_facing(pdb, first_residue):
    import mdtraj
    traj = mdtraj.load(pdb)
    top, xyz = traj.topology, traj.xyz[0] * 10.0
    dssp = mdtraj.compute_dssp(traj, simplified=True)[0]
    ca = {a.residue.index: a.index for a in top.atoms if a.name == "CA"}
    cb = {a.residue.index: a.index for a in top.atoms if a.name == "CB"}
    runs, start = [], None
    for i, s in enumerate(list(dssp) + ["-"]):
        if s == "E" and start is None:
            start = i
        elif s != "E" and start is not None:
            if i - start >= 5:
                runs.append((start, i - 1))
            start = None
    ref = xyz[ca[runs[0][1]]] - xyz[ca[runs[0][0]]]
    dirs = [(lambda d: d if d @ ref > 0 else -d)(xyz[ca[b]] - xyz[ca[a]]) for a, b in runs]
    axis = np.mean(dirs, axis=0)
    axis /= np.linalg.norm(axis)
    centre = xyz[[ca[i] for a, b in runs for i in range(a, b + 1)]].mean(0)
    radial = lambda p: np.linalg.norm((p - centre) - ((p - centre) @ axis) * axis)
    facing = {}
    for r in ca:
        if r in cb:
            facing[r + first_residue] = bool(radial(xyz[cb[r]]) < radial(xyz[ca[r]]))
    in_strand = {i + first_residue for a, b in runs for i in range(a, b + 1)}
    return facing, in_strand


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--msas", required=True)
    p.add_argument("--pdb", required=True)
    p.add_argument("--query", default="EEGFDTEHIFGFMIGTDVGNPGEREFQSQTTGRFGKGGGTYRALAQEVEIEVVPLPNFRIELGATASLHDITGVPDIDDRRQFNFQGASLDLRYRLLDRERAPFGLTVAAELHGDRIDETSGAKGRMYGTEFTLAFDRELIPNFAIGALNLIYQPEWARFDVAGLSEKSSTIGAAFAGLVRVRPNMLLGGELRYFRQYEGIGLGEFSGQALFVGPTAYFQLSERSRLTLSWSMQAWGRPAGSGANLDLVNFERHQARLVFGVNF")
    p.add_argument("--first-residue", type=int, default=110)
    p.add_argument("--segment", default="110-132", help="N-terminal segment in the lumen")
    p.add_argument("--max-seqs", type=int, default=3000)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    d = Path(args.msas)
    msa = []
    for name in ("uniref90_hits.sto", "mgnify_hits.sto", "bfd_uniref_hits.a3m"):
        f = d / name
        if f.exists():
            seqs = read_sto(f) if f.suffix == ".sto" else read_a3m(f)
            msa += [s for s in seqs[1:] if len(s) == len(args.query)]
    # one copy of each sequence, the query excluded, then a cap for the O(N^2) weighting
    msa = [s for s in dict.fromkeys(msa) if s != args.query][:args.max_seqs]
    w = weights(msa)
    x = np.array([list(s) for s in msa])
    facing, in_strand = lumen_facing(args.pdb, args.first_residue)
    s0, s1 = (int(v) for v in args.segment.split("-"))

    rows = []
    for i, q in enumerate(args.query):
        col = x[:, i]
        present = col != "-"
        cov = float(w[present].sum() / w.sum())
        ident = float(w[col == q].sum() / w.sum())
        freq = np.array([w[col == a].sum() for a in AA])
        cons = None
        if freq.sum() > 0:
            pr = freq / freq.sum()
            ent = -(pr[pr > 0] * np.log(pr[pr > 0])).sum()
            cons = round(float(1 - ent / np.log(20)), 3)
        resnum = i + args.first_residue
        region = "segment" if s0 <= resnum <= s1 else ("strand" if resnum in in_strand else "loop")
        rows.append({"res": f"{q}{resnum}", "region": region, "coverage": round(cov, 3),
                     "identity": round(ident, 3), "conservation": cons,
                     "lumen": facing.get(resnum)})

    def avg(sel, key):
        vals = [r[key] for r in rows if sel(r) and r[key] is not None]
        return round(float(np.mean(vals)), 3) if vals else None

    summary = {}
    for region in ("segment", "strand", "loop"):
        sel = lambda r, g=region: r["region"] == g
        summary[region] = {k: avg(sel, k) for k in ("coverage", "identity", "conservation")}
    summary["strand_lumen_facing"] = {k: avg(lambda r: r["region"] == "strand" and r["lumen"], k)
                                      for k in ("coverage", "identity", "conservation")}
    report = {"sequences": len(msa), "effective_sequences": round(float(w.sum()), 1),
              "summary": summary, "residues": rows}
    Path(args.out).write_text(json.dumps(report, indent=1))
    print(json.dumps({k: v for k, v in report.items() if k != "residues"}, indent=1))
    for r in rows:
        if r["region"] == "segment":
            print(r)


if __name__ == "__main__":
    main()
