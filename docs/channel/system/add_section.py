"""Add an umbrella section (inputs only) to an ion's config.yaml, under mcube's config lock,
so running window jobs can keep writing their results.

    python add_section.py config.yaml umbrella_chbe k=2.5 c4_set=kantakevicius2022_chbe_tip3p

Values are parsed as YAML (2.5 -> float, true -> bool). Windows of the new section seed from
the `umbrella` window at the same z0 (its last frame), so the grid (z_start, z_end, spacing)
stays the umbrella one unless given here. Decision log F12.
"""
import os
import sys
from pathlib import Path

import yaml

from mcube.preprocess import _locked


def main(config, section, *pairs):
    if section == "umbrella":
        raise SystemExit("give a new section name, not umbrella")
    extra = {}
    for p in pairs:
        key, _, value = p.partition("=")
        extra[key] = yaml.safe_load(value)
    path = Path(config)
    with _locked(path):
        data = yaml.safe_load(path.read_text()) or {}
        node = data.setdefault(section, {})
        node.setdefault("inputs", {}).update(extra)
        tmp = path.with_name(f"{path.name}.{os.getpid()}.tmp")
        tmp.write_text(yaml.safe_dump(data, sort_keys=False, default_flow_style=False))
        tmp.replace(path)
    print(f"{config}: {section}.inputs = {node['inputs']}")


if __name__ == "__main__":
    main(*sys.argv[1:])
