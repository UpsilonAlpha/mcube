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

    out = strip_results(config)
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
    dropped = sorted(set(config) - set(out)) or "none at the top level"
    print(f"  results/timestamp stripped; stages kept: {sorted(out)}")
    print("  next: python sections.py config.yaml ohHmd ohHqm")


if __name__ == "__main__":
    main(sys.argv[1:])
