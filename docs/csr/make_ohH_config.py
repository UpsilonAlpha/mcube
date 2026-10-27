#!/usr/bin/env python3
"""Write the ohH arm's config from the lysH arm's, changing three things.

    python make_ohH_config.py CSR_nac/config.yaml CSR_ohH/config.yaml

The two arms must differ by where one hydrogen sits and in nothing else, so the
config is derived rather than written: every `inputs` block is copied verbatim
and only these three preprocess keys change.

    protein_pdb      -> ./csr_crystal_start_ohH.pdb   (from make_water_nuc.py)
    forcefields      -> ../water_nuc.xml replaces ../hydroxide.xml
    residue_variants -> 226: LYN added to the existing HIS variants

`results` and `timestamp` blocks are dropped: they belong to the lysH run, and
a stale `results.lastframe` is exactly how a stage silently starts from another
run's structure. Sections are not copied either - add them with
`python sections.py site.yaml ohHmd ohHqm` once this has written the file, so
site.yaml stays the single source for what each calculation is.

Refuses rather than guesses if the source is not the config this expects.
"""

from __future__ import annotations

import sys
from pathlib import Path

BASE_STAGES = ("preprocess", "md", "qmmm", "sp", "scan", "neb", "ts", "umbrella")
NEW_PDB = "./csr_crystal_start_ohH.pdb"
OLD_XML, NEW_XML = "../hydroxide.xml", "../water_nuc.xml"


def strip_results(node):
    """Every `results` and `timestamp` key, at any depth."""
    if isinstance(node, dict):
        return {k: strip_results(v) for k, v in node.items()
                if k not in ("results", "timestamp")}
    if isinstance(node, list):
        return [strip_results(v) for v in node]
    return node


def main(argv):
    if len(argv) != 2:
        raise SystemExit(__doc__)
    src, dest = Path(argv[0]), Path(argv[1])
    if not src.is_file():
        raise SystemExit(f"no such file: {src}")
    if dest.exists():
        raise SystemExit(f"{dest} already exists; refusing to overwrite it")

    import yaml

    config = yaml.safe_load(src.read_text())
    pre = (config.get("preprocess") or {}).get("inputs")
    if not pre:
        raise SystemExit(f"{src} has no preprocess.inputs")

    forcefields = list(pre.get("forcefields") or [])
    if OLD_XML not in forcefields:
        raise SystemExit(
            f"preprocess.inputs.forcefields does not list {OLD_XML}, so this is "
            f"not the lysH config this script expects: {forcefields}")

    variants = pre.get("residue_variants") or {}
    chain_a = dict(variants.get("A") or {})
    if 226 in chain_a:
        raise SystemExit(
            f"residue 226 already has a variant ({chain_a[226]}); refusing to "
            "change it")
    if not chain_a:
        raise SystemExit(
            "preprocess.inputs.residue_variants has no chain A entries. The "
            "HIS tautomers must carry over, and require_his_variants would "
            "reject the run without them.")

    # Named sections are dropped, not copied. They are top-level keys like
    # `nacqm` beside the base stages, and the lysH arm's carry `pdb_file`
    # pinned to CSR_nac paths - an explicit pdb_file is honoured, so running
    # one here would silently compute on the lysH structure. Sections come
    # from site.yaml via sections.py instead, which is the single source for
    # what each calculation is.
    sections = [k for k in config if k not in BASE_STAGES]
    out = strip_results({k: v for k, v in config.items() if k in BASE_STAGES})
    pre_out = out["preprocess"]["inputs"]
    pre_out["protein_pdb"] = NEW_PDB
    pre_out["forcefields"] = [NEW_XML if f == OLD_XML else f
                              for f in forcefields]
    chain_a[226] = "LYN"
    pre_out["residue_variants"] = dict(variants, A=chain_a)
    # Record the pH that was previously left to preprocess's own default, so
    # the assumed pH is in the config rather than in the code (it defaults to
    # 8.0 on the explicit-solvent path, which is what the lysH arm used; pH 7
    # and 8 are identical under OpenMM's thresholds, so this changes nothing
    # and documents everything).
    pre_out.setdefault("pH", 8.0)

    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(yaml.safe_dump(out, sort_keys=False, default_flow_style=False))

    print(f"wrote {dest}")
    print(f"  protein_pdb      {pre_out['protein_pdb']}")
    print(f"  forcefields      {pre_out['forcefields']}")
    print(f"  residue_variants {pre_out['residue_variants']}")
    print(f"  pH               {pre_out['pH']}")
    print(f"  stages kept:      {sorted(out)}")
    print(f"  sections dropped: {sorted(sections) or 'none'}")
    print("  next: python sections.py config.yaml ohHmd ohHqm")


if __name__ == "__main__":
    main(sys.argv[1:])
