"""Add calculations from site.yaml to config.yaml, checking every selector first.

    python sections.py site.yaml <name> [<name> ...]
    python sections.py site.yaml --list

Replaces the one-off `mk_*.py` scripts: a new calculation is a few lines of YAML
in site.yaml rather than a new Python file. Each named section is written to
config.yaml as `<stage>_<name>`, which is the section name to pass to mcube.

Nothing is written unless every structure the section names exists and every
atom selector in it resolves in that structure. The selector check exists
because a scan whose second coordinate starts beyond hydrogen-bond distance is
not a proton transfer, and two days of compute were lost to exactly that before
it was added.

The structure is `inputs.pdb_file` (for `neb`, `reactant_pdb` and
`product_pdb`), which is also what mcube reads: there is no separate `start` key
any more and no structure on the command line, so one place says which structure
a calculation used and that place is checked here. Paths are read against the
config's directory, as mcube reads them.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import yaml

from mcube.preprocess import _locked

# `restraints` belongs in this list and was missing until D128: an `md` section
# whose only selectors are in its restraint list printed `0 selectors checked`
# and wrote the section anyway. That is the one place a bad selector is most
# expensive - `md.md_restraints` passes the list through `resolve_atoms`, so a
# name that does not resolve surfaces hours into a run rather than here. The
# third element of each row is a distance, not a selector, and `walk` ignores
# it because it only collects strings with exactly two colons.
SELECTOR_KEYS = ("atoms", "atoms2", "opt_constraints", "constraints", "coordination",
                 "monitor", "hessian_atoms", "partial_hessian_atoms",
                 "actregion_origin", "restraints")

# Structure keys, in the order a section's selectors are checked against: for a
# NEB the two endpoints are one system in two geometries, so the reactant's
# atom order is the one the selectors have to resolve in.
STRUCTURE_KEYS = ("pdb_file", "reactant_pdb", "product_pdb")


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
            inputs = body.get("inputs") or {}
            print("  %-14s %-6s %s" % (name, body.get("stage", "?"),
                                       inputs.get("pdb_file")
                                       or inputs.get("reactant_pdb")
                                       or "(from the previous stage)"))
        return

    unknown = [n for n in names if n not in available]
    if unknown:
        raise SystemExit("not in %s: %s" % (site_file, ", ".join(unknown)))

    # The config's directory, which is what mcube resolves these against - not
    # site.yaml's, which is in the repo and holds none of the run's structures.
    config = Path("config.yaml")
    if not config.is_file():
        raise SystemExit("no config.yaml here; run this from the run directory")
    home = config.resolve().parent
    # Check everything before writing anything.
    planned, problems = {}, []
    for name in names:
        body = available[name]
        stage, inputs = body.get("stage", "qmmm"), dict(body.get("inputs") or {})
        section = "%s_%s" % (stage, name)
        if "start" in body:
            problems.append("%s: `start` is gone; put the structure in "
                            "`inputs.pdb_file` (neb: `reactant_pdb` and "
                            "`product_pdb`), where mcube reads it" % name)
            continue
        if stage == "neb" and not all(inputs.get(k) for k in ("reactant_pdb", "product_pdb")):
            problems.append("%s: a neb section needs both `reactant_pdb` and "
                            "`product_pdb` in `inputs`" % name)
            continue
        # Paths are read against the config's directory, as mcube reads them.
        given = [(key, inputs[key], home / inputs[key])
                 for key in STRUCTURE_KEYS if inputs.get(key)]
        for key, value, path in given:
            if not path.exists():
                problems.append("%s: %s missing: %s" % (name, key, value))
        start = given[0][2] if given and given[0][2].exists() else None
        sels = selectors_in(inputs)
        if start and sels:
            index = atom_index(start)
            missing = [s for s in sels if s not in index]
            if missing:
                problems.append("%s: selectors not in %s: %s"
                                % (name, given[0][1], ", ".join(missing)))
        planned[section] = (inputs, given[0][1] if given else None, sels)

    if problems:
        print("NOT WRITING. Problems:")
        for p in problems:
            print("  -", p)
        raise SystemExit(1)

    with _locked(config):
        data = yaml.safe_load(open(config))
        for section, (inputs, _, _) in planned.items():
            data[section] = {"inputs": inputs}
        yaml.safe_dump(data, open("config.yaml.edit.tmp", "w"),
                       sort_keys=False, default_flow_style=False)
        os.replace("config.yaml.edit.tmp", config)

    for section, (inputs, start, sels) in planned.items():
        print("%-26s pdb=%s  %d selectors checked"
              % (section, start or "(from the previous stage)", len(sels)))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2:])
