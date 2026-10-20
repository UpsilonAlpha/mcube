"""Add calculations from site.yaml to config.yaml, checking every selector first.

    python sections.py site.yaml <name> [<name> ...]
    python sections.py site.yaml --list

Replaces the one-off `mk_*.py` scripts: a new calculation is a few lines of YAML
in site.yaml rather than a new Python file. Each named section is written to
config.yaml as `<stage>_<name>`, which is the section name to pass to mcube.

Nothing is written unless every atom selector in the section resolves in that
section's own `start` structure. That check exists because a scan whose second
coordinate starts beyond hydrogen-bond distance is not a proton transfer, and
two days of compute were lost to exactly that before it was added.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import yaml

from mcube.preprocess import _locked

SELECTOR_KEYS = ("atoms", "atoms2", "opt_constraints", "constraints",
                 "monitor", "hessian_atoms", "partial_hessian_atoms",
                 "actregion_origin")


def selectors_in(inputs):
    """Every 'chain:resid:atom' string anywhere in a section's inputs."""
    found = []

    def walk(value):
        if isinstance(value, str) and value.count(":") == 2:
            found.append(value)
        elif isinstance(value, (list, tuple)):
            for v in value:
                walk(v)
        elif isinstance(value, dict):
            for v in value.values():
                walk(v)

    for key in SELECTOR_KEYS:
        if key in inputs:
            walk(inputs[key])
    return found


def atom_index(pdb_path):
    index = {}
    for i, line in enumerate(l for l in open(pdb_path)
                             if l.startswith(("ATOM", "HETATM"))):
        index["%s:%s:%s" % (line[21], line[22:26].strip(), line[12:16].strip())] = i
    return index


def main(site_file, names):
    spec = yaml.safe_load(open(site_file))
    available = spec.get("sections") or {}
    if not names or names == ["--list"]:
        print("sections defined in %s:" % site_file)
        for name, body in available.items():
            print("  %-14s %-6s %s" % (name, body.get("stage", "?"),
                                       body.get("start", "(start given at run time)")))
        return

    unknown = [n for n in names if n not in available]
    if unknown:
        raise SystemExit("not in %s: %s" % (site_file, ", ".join(unknown)))

    # Check everything before writing anything.
    planned, problems = {}, []
    for name in names:
        body = available[name]
        stage, inputs = body.get("stage", "qmmm"), dict(body.get("inputs") or {})
        section = "%s_%s" % (stage, name)
        start = body.get("start")
        sels = selectors_in(inputs)
        if start and not Path(start).exists():
            problems.append("%s: start structure missing: %s" % (name, start))
        elif start and sels:
            index = atom_index(start)
            missing = [s for s in sels if s not in index]
            if missing:
                problems.append("%s: selectors not in %s: %s"
                                % (name, start, ", ".join(missing)))
        planned[section] = (inputs, start, sels)

    if problems:
        print("NOT WRITING. Problems:")
        for p in problems:
            print("  -", p)
        raise SystemExit(1)

    config = Path("config.yaml")
    with _locked(config):
        data = yaml.safe_load(open(config))
        for section, (inputs, _, _) in planned.items():
            data[section] = {"inputs": inputs}
        yaml.safe_dump(data, open("config.yaml.edit.tmp", "w"),
                       sort_keys=False, default_flow_style=False)
        os.replace("config.yaml.edit.tmp", config)

    for section, (inputs, start, sels) in planned.items():
        print("%-26s start=%s  %d selectors checked"
              % (section, start or "(at run time)", len(sels)))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2:])
