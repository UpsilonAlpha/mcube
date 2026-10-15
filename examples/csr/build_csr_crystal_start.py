#!/usr/bin/env python3
"""
Build the CSR Michaelis-complex starting structure from the crystal model.

Decisions (see docs/csr-decision-log.md, D21-D23):
- Protein: heavy atoms of the crystal model (H stripped, altloc A kept,
  crystal/BBL numbering). Hydrogens are re-added by mcube preprocess
  with explicit His tautomers.
- Zn and its waters: taken from PDB 6NC5 (CSR-1 Q121H/K254H, Arg118 kept;
  numbering = crystal - 29) by CA superposition. Zn_alpha (His87/His163 NE2
  + two waters) is used; Zn_beta and adventitious Zn are not. The water not
  bridging to Zn_beta becomes OH- (H pointed at the nearest acceptor); the
  other is kept as a Zn-bound water. (The deposited Zn of the crystal model
  clashes with Asp120/Lys254 and is discarded.)
- Meropenem: pose from the processed model (residue numbers -28),
  superposed on CA atoms, then the smallest rigid-body shift that clears
  contacts and puts C4x 3.0-3.2 A from the OH- O at a 100-115 deg
  O...C4x=O2x angle.

Usage:
    python build_csr_crystal_start.py original.pdb processed.pdb site.pdb outdir \
        [--template-config CSR_ligand/config.yaml]
"""

import argparse
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation

OFFSET = 28          # processed numbering + 28 = crystal numbering
SITE_OFFSET = 29     # 6NC5 numbering + 29 = crystal numbering
HIS_ZN = (116, 192)  # His whose NE2 binds Zn (crystal numbering)
ASP = 120
HIS_VARIANTS = {81: "HIE", 116: "HID", 140: "HID", 178: "HIE", 192: "HID"}  # from crystal riding H
EXTRA_QM = [118, 254]  # Arg118 (=Arg90), Lys254 (=Lys226)


def read_pdb(path):
    atoms, conect, seen = [], [], set()
    for line in Path(path).read_text().splitlines():
        line = line.rstrip("\r")
        if line.startswith(("ATOM", "HETATM")):
            if len(line) < 66:  # truncated/malformed record
                continue
            key = (line[21], line[22:27], line[12:16], line[16])
            if key in seen:  # duplicated record
                continue
            seen.add(key)
            atoms.append(dict(line=line, serial=int(line[6:11]), name=line[12:16].strip(), alt=line[16],
                              res=line[17:20].strip(), chain=line[21], rid=int(line[22:26]),
                              xyz=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                              el=(line[76:78].strip() or line[12:16].strip()[0]).upper()))
        elif line.startswith("CONECT"):
            conect.append([int(line[i:i + 5]) for i in range(6, len(line.rstrip()), 5) if line[i:i + 5].strip()])
    return atoms, conect


def kabsch(mobile, target):
    mc, tc = mobile.mean(0), target.mean(0)
    u, _, vt = np.linalg.svd((mobile - mc).T @ (target - tc))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1, 1, d]) @ u.T
    return rot, tc - mc @ rot.T


def angle(a, b, c):
    u, v = a - b, c - b
    return np.degrees(np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("original")
    p.add_argument("processed")
    p.add_argument("site", help="PDB with the Zn site to transplant (6NC5)")
    p.add_argument("outdir")
    p.add_argument("--template-config")
    args = p.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    orig, _ = read_pdb(args.original)
    prot = [a for a in orig if a["chain"] == "A" and a["res"] not in ("HOH", "ZN") and a["el"] != "H"
            and a["alt"] in (" ", "A")]
    print(f"Protein heavy atoms kept: {len(prot)} (H, waters, metals and altloc B removed)")

    def get(rid, name):
        hits = [a for a in prot if a["rid"] == rid and a["name"] == name]
        assert len(hits) == 1, (rid, name)
        return hits[0]["xyz"]

    nes = [get(r, "NE2") for r in HIS_ZN]
    asp_os = [get(ASP, "OD1"), get(ASP, "OD2")]

    # --- Zn_alpha and its waters from the site structure ---
    site, _ = read_pdb(args.site)
    site_prot = [a for a in site if a["line"].startswith("ATOM") and a["alt"] in (" ", "A")]
    pairs = []
    for a in site_prot:
        if a["name"] == "CA":
            b = [x for x in prot if x["name"] == "CA" and x["rid"] == a["rid"] + SITE_OFFSET]
            if b and b[0]["res"] == a["res"]:
                pairs.append((a["xyz"], b[0]["xyz"]))
    srot, strans = kabsch(np.array([x for x, _ in pairs]), np.array([y for _, y in pairs]))
    tf = lambda x: x @ srot.T + strans
    s_rmsd = np.sqrt(np.mean([np.sum((tf(x) - y) ** 2) for x, y in pairs]))
    print(f"Site structure -> crystal superposition on {len(pairs)} CA: RMSD {s_rmsd:.3f} A")
    sget = lambda rid, name: [a for a in site_prot if a["rid"] == rid - SITE_OFFSET and a["name"] == name][0]["xyz"]
    s_nes = [sget(r, "NE2") for r in HIS_ZN]
    zns = [a for a in site if a["el"] == "ZN"]
    zn_a = min(zns, key=lambda z: sum(np.linalg.norm(z["xyz"] - n) for n in s_nes))
    others_zn = [z for z in zns if z is not zn_a]
    wats = [a for a in site if a["res"] == "HOH" and np.linalg.norm(a["xyz"] - zn_a["xyz"]) < 2.6]
    zn_b = min(others_zn, key=lambda z: np.linalg.norm(z["xyz"] - sget(ASP, "OD2")))
    bridging = [w for w in wats if np.linalg.norm(w["xyz"] - zn_b["xyz"]) < 2.6]
    terminal = [w for w in wats if w not in bridging]
    print(f"Site Zn_alpha {zn_a['chain']}{zn_a['rid']}; waters: " + ", ".join(
        f"HOH{w['rid']} ({'bridging to Zn_beta' if w in bridging else 'terminal'})" for w in wats))
    assert len(terminal) == 1 and len(bridging) == 1, "expected one terminal and one bridging water"
    zn = tf(zn_a["xyz"])
    oh_o = tf(terminal[0]["xyz"])
    wat_o = tf(bridging[0]["xyz"])
    # OH- hydrogen towards the nearest acceptor O (protein or the other water), else away from Zn
    acc = [a["xyz"] for a in prot if a["el"] == "O"] + [wat_o]
    acc = [x for x in acc if 2.4 < np.linalg.norm(x - oh_o) < 3.2]
    target = min(acc, key=lambda x: np.linalg.norm(x - oh_o)) if acc else oh_o + (oh_o - zn)
    oh_h = oh_o + 0.97 * (target - oh_o) / np.linalg.norm(target - oh_o)
    for label, x in (("Zn", zn), ("OH- O", oh_o), ("water O", wat_o)):
        close = sorted((np.linalg.norm(a["xyz"] - x), f"{a['res']}{a['rid']}:{a['name']}") for a in prot)[:4]
        print(f"  {label}: " + ", ".join(f"{n} {d:.2f}" for d, n in close))
    asp_o = min(asp_os, key=lambda o: np.linalg.norm(o - oh_o))

    rng = np.random.default_rng(0)

    # --- Meropenem from the processed model ---
    proc, proc_conect = read_pdb(args.processed)
    ca_pairs = []
    for a in proc:
        if a["name"] == "CA" and a["chain"] == "A":
            b = [x for x in prot if x["name"] == "CA" and x["rid"] == a["rid"] + OFFSET]
            if b and b[0]["res"] == a["res"]:
                ca_pairs.append((a["xyz"], b[0]["xyz"]))
    mob, tgt = np.array([x for x, _ in ca_pairs]), np.array([y for _, y in ca_pairs])
    rot, trans = kabsch(mob, tgt)
    ca_rmsd = np.sqrt(np.mean(np.sum((mob @ rot.T + trans - tgt) ** 2, axis=1)))
    print(f"Processed -> crystal superposition on {len(ca_pairs)} CA: RMSD {ca_rmsd:.3f} A")

    mer = [dict(a) for a in proc if a["res"] == "mer"]
    for a in mer:
        a["xyz"] = a["xyz"] @ rot.T + trans
    mer_xyz0 = np.array([a["xyz"] for a in mer])
    heavy = np.array([a["el"] != "H" for a in mer])
    idx = {a["name"]: i for i, a in enumerate(mer)}
    c4, o2 = idx["C4x"], idx["O2x"]
    env = np.vstack([np.array([a["xyz"] for a in prot]), zn, oh_o, wat_o])
    centroid = mer_xyz0.mean(0)

    def place(x):
        return (mer_xyz0 - centroid) @ Rotation.from_rotvec(x[:3]).as_matrix().T + centroid + x[3:]

    def lig_cost(x):
        xyz = place(x)
        d = np.linalg.norm(xyz[:, None, :] - env[None, :, :], axis=2)
        cut = np.where(heavy, 3.0, 2.2)[:, None]
        c = 10 * np.sum(np.clip(cut - d, 0, None) ** 2)
        r = np.linalg.norm(xyz[c4] - oh_o)
        c += 20 * (max(0, r - 3.2) ** 2 + max(0, 3.0 - r) ** 2)
        a = angle(oh_o, xyz[c4], xyz[o2])
        c += 0.02 * (max(0, a - 115) ** 2 + max(0, 100 - a) ** 2)
        c += 0.05 * np.mean(np.sum((xyz - mer_xyz0) ** 2, axis=1))  # stay close to the original pose
        return c
    best = None
    for k in range(60):
        x0 = np.zeros(6) if k == 0 else np.concatenate([rng.normal(scale=0.3, size=3), rng.normal(scale=1.5, size=3)])
        res = minimize(lig_cost, x0, method="Powell", options={"maxiter": 20000, "xtol": 1e-4, "ftol": 1e-8})
        if best is None or res.fun < best.fun:
            best = res
    mer_xyz = place(best.x)
    for a, x in zip(mer, mer_xyz):
        a["xyz"] = x
    disp = np.linalg.norm(mer_xyz - mer_xyz0, axis=1)
    dmin = np.linalg.norm(mer_xyz[:, None, :] - env[None, :, :], axis=2)
    print(f"Meropenem shift: RMSD {np.sqrt(np.mean(disp ** 2)):.2f} A, max {disp.max():.2f} A, "
          f"rotation {np.degrees(np.linalg.norm(best.x[:3])):.1f} deg")
    print(f"C4x...O(OH) {np.linalg.norm(mer_xyz[c4] - oh_o):.2f} A, O...C4x=O2x {angle(oh_o, mer_xyz[c4], mer_xyz[o2]):.0f} deg")
    print(f"Closest ligand heavy atom - environment: {dmin[heavy].min():.2f} A; "
          f"ligand H - environment: {dmin[~heavy].min():.2f} A")

    # --- Write the start PDB ---
    lines, serial = [], 0
    def rec(kind, name, res, chain, rid, xyz, el):
        nonlocal serial
        serial += 1
        nm = name if len(name) == 4 else f" {name:<3}"
        lines.append(f"{kind:<6}{serial:5d} {nm} {res:>3} {chain}{rid:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {el:>2}")
        return serial
    prev = None
    for a in prot:
        if prev is not None and a["chain"] != prev:
            lines.append("TER")
        prev = a["chain"]
        serial += 1
        l = a["line"]
        lines.append(f"ATOM  {serial:5d}{l[11:16]} {l[17:54]}  1.00  0.00          {a['el']:>2}")
    lines.append("TER")
    rec("HETATM", "ZN", "ZN", "B", 1, zn, "ZN")
    lines.append("TER")
    ser_of = {}
    for a in mer:
        ser_of[a["serial"]] = rec("HETATM", a["name"], "mer", "C", 1, a["xyz"], a["el"])
    o_ser = rec("HETATM", "O01", "OH", "C", 2, oh_o, "O")
    h_ser = rec("HETATM", "H03", "OH", "C", 2, oh_h, "H")
    lines.append("TER")
    rec("HETATM", "O", "HOH", "W", 1, wat_o, "O")  # Zn-bound water from the site structure; H added by preprocess
    lines.append("TER")
    bonds = {}
    for c in proc_conect:
        if c[0] in ser_of:
            bonds.setdefault(ser_of[c[0]], set()).update(ser_of[b] for b in c[1:] if b in ser_of)
    bonds[o_ser], bonds[h_ser] = {h_ser}, {o_ser}
    conect = ["CONECT" + f"{s:5d}" + "".join(f"{b:5d}" for b in sorted(bs)) for s, bs in sorted(bonds.items())]
    (out / "csr_crystal_start.pdb").write_text("\n".join(lines + conect + ["END"]) + "\n")

    # Ligand alone, for parameterisation
    renum = {ser_of[a["serial"]]: i + 1 for i, a in enumerate(mer)}
    lig = [l for l in lines if l[17:20] == "mer"]
    lig = [f"{l[:6]}{renum[int(l[6:11])]:5d}{l[11:]}" for l in lig]
    lig += ["CONECT" + f"{renum[s]:5d}" + "".join(f"{renum[b]:5d}" for b in sorted(bs))
            for s, bs in sorted(bonds.items()) if s in renum]
    (out / "mer_crystal.pdb").write_text("\n".join(lig + ["END"]) + "\n")
    print(f"Wrote {out / 'csr_crystal_start.pdb'} ({serial} atoms) and {out / 'mer_crystal.pdb'} ({len(mer)} atoms)")
    print(f"Zn: His116 {np.linalg.norm(zn - nes[0]):.2f}, His192 {np.linalg.norm(zn - nes[1]):.2f}, "
          f"Asp120 {min(np.linalg.norm(zn - o) for o in asp_os):.2f}, Arg118 NH1 {np.linalg.norm(zn - get(118, 'NH1')):.2f} A; "
          f"OH: Zn-O {np.linalg.norm(oh_o - zn):.2f}, O..Asp {np.linalg.norm(oh_o - asp_o):.2f} A; "
          f"water: Zn-O {np.linalg.norm(wat_o - zn):.2f} A")
    print(f"Attack geometry: OH-...C4x {np.linalg.norm(mer_xyz[c4] - oh_o):.2f} A, water...C4x {np.linalg.norm(mer_xyz[c4] - wat_o):.2f} A")

    if args.template_config:
        import yaml
        c = yaml.safe_load(open(args.template_config))
        here = out.parent.resolve().as_posix()
        i = c["preprocess"]["inputs"]
        i.update(protein_pdb="./csr_crystal_start.pdb", ligand="./mer_crystal.pdb", merge_ligand=False,
                 ligand_resname="MER", residue_variants={"A": HIS_VARIANTS}, require_his_variants=True,
                 extra_qm_residues={"A": EXTRA_QM}, constraint_exclude_residues={"C": [1]}, cutoff=3.0)
        i.pop("donor_atoms", None)
        i["donor_elements"] = ["N", "O", "S"]
        for sec in c.values():
            if isinstance(sec, dict):
                sec.pop("results", None)
                sec.pop("timestamp", None)
                sec.pop("metadata", None)
        q = c["qmmm"]["inputs"]
        q["forcefields"] = ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml", f"{here}/hydroxide.xml",
                            f"{here}/preprocess_results/openff_MER.xml"]
        q["actregion_origin"] = None  # Zn index changes; set after preprocess
        yaml.safe_dump(c, open(out.parent / "config.yaml", "w"), sort_keys=False)
        print(f"Wrote {out.parent / 'config.yaml'}")


if __name__ == "__main__":
    main()
