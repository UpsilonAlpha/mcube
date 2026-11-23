#!/usr/bin/env python3
"""Does every section run the region of the section that produced its structure?

    python region_check.py ../../mcube/docs/csr/site.yaml        # from the run dir
    python region_check.py SITE.yaml --config config.yaml --strict

D207's two-surfaces bug has now been found three times, each time by hand and
each time one stage further on:

    D207   endpoint optimisations held active-region waters rigid while the NEB
           left them free                                  -> rigid_water: false
    D317   aspxtb8/attackmin relaxed 481 atoms at 8 A while step1 freed 665 at
           9 A. Every interior image came out ~100 kcal/mol below both
           endpoints, all of it in the MM term, and the run was unusable
    D320   saddle1/saddle2 inherited `actradius: 8` from the base `qmmm.inputs`
           while refining a 9 A path. That one would have landed on the
           PUBLISHED barrier rather than on an intermediate quantity

They are one defect: **a structure optimised with one movable set, consumed by a
calculation using another.** `frozen_check.py` checks the companion condition -
that the region contains every atom differing between two endpoints - and nothing
checked this one, which is why it kept being found by hand (D318).

This reports it. It is deliberately a separate script rather than a `require`
key: `sections.py` is shared with the channel project, `require` rejects keys it
does not know, and the blast radius of making this fatal has not been measured.
Run it at each stage handoff, as `frozen_check.py` is run.

**What it compares.** For each section with a structure input, the producing
section is read off the path - mcube writes `<stage>_<name>`'s output under
`protdft/<stage>_<name>/` - and these keys are compared as the two runs
EFFECTIVELY resolve them, section value first and the base `qmmm.inputs` second:

    qm_residues        the QM region. D316's defect: 124 atoms against 129
    actradius          the active region's radius. D317's and D320's defect
    actregion_origin   what it is centred on

Resolving against the base is the whole point. `saddle1` set no `actradius` at
all, so reading only its own `inputs` would have called it consistent while it
silently inherited 8 from the base block.

A difference is not automatically wrong - `attackmin` changed `qm_residues` on
purpose and said so in `require.same_as.except`. This script does not try to
read intent; it prints what differs and leaves the judgement where D318 left it.
`--strict` exits 1 on any difference, for a future gate.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import yaml

# The region a run resolves, in `resolve_actatoms`'s own terms.
REGION_KEYS = ("qm_residues", "actradius", "actregion_origin")

# Structure inputs, in the order a section's own region is derived against.
STRUCTURE_KEYS = ("pdb_file", "reactant_pdb", "product_pdb")

# mcube writes a section's output under `protdft/<stage>_<name>/`, which is the
# only link from a consumed file back to the run that produced it.
OUTPUT_DIR = re.compile(r"(?:^|/)protdft/([a-z]+_[A-Za-z0-9_]+)/")


def producing_section(path: str):
    """The `<stage>_<name>` that wrote this file, or None if nothing did."""
    found = OUTPUT_DIR.search(path or "")
    return found.group(1) if found else None


def effective(inputs, base, key):
    """A key as the run resolves it: its own value, else the base block's."""
    if key in (inputs or {}):
        return (inputs or {})[key]
    return (base or {}).get(key)


def declared_except(body):
    """Keys a section's `require.same_as` already says may differ."""
    same = ((body.get("require") or {}).get("same_as")) or {}
    if isinstance(same, str):
        return set()
    return set(same.get("except") or [])


def describe(value):
    if value is None:
        return "(unset)"
    if isinstance(value, dict):
        # A QM region: chain -> {resid: name}. Summarise, since the full dict is
        # unreadable in a table and the point is only whether two differ.
        return "{%s}" % ", ".join(
            "%s:%s" % (chain, ",".join(str(r) for r in sorted(residues)))
            for chain, residues in sorted(value.items()))
    return repr(value)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("site")
    parser.add_argument("--config", default="config.yaml",
                        help="the run's config.yaml, for the base qmmm.inputs "
                             "and for sections already written (default: ./config.yaml)")
    parser.add_argument("--strict", action="store_true",
                        help="exit 1 if any section differs from its producer")
    parser.add_argument("--only", help="check just this section")
    args = parser.parse_args(argv)

    site = yaml.safe_load(open(args.site)) or {}
    sections = site.get("sections") or {}

    config_path = Path(args.config)
    written = {}
    base = {}
    if config_path.is_file():
        written = yaml.safe_load(open(config_path)) or {}
        base = ((written.get("qmmm") or {}).get("inputs")) or {}
    else:
        print("[WARN] no %s, so the base `qmmm.inputs` is unknown and a section "
              "that inherits a key will read as `(unset)` rather than as what it "
              "will actually run. That is exactly D320's defect, so this run "
              "cannot rule it out." % args.config)

    def inputs_of(section_key):
        """A written section's inputs, else the site file's, else None."""
        from_config = (written.get(section_key) or {}).get("inputs")
        if from_config is not None:
            return from_config, "config.yaml"
        for name, body in sections.items():
            if "%s_%s" % (body.get("stage", "qmmm"), name) == section_key:
                return (body.get("inputs") or {}), "site.yaml"
        return None, "nowhere"

    rows, differing, unknown_producer = [], [], []
    for name, body in sorted(sections.items()):
        if args.only and name != args.only:
            continue
        stage = body.get("stage", "qmmm")
        inputs = body.get("inputs") or {}
        allowed = declared_except(body)
        for key in STRUCTURE_KEYS:
            path = inputs.get(key)
            if not path:
                continue
            producer = producing_section(path)
            if producer is None:
                unknown_producer.append((name, key, path))
                continue
            produced_inputs, where = inputs_of(producer)
            if produced_inputs is None:
                unknown_producer.append((name, key, path + "  [no such section]"))
                continue
            diffs = []
            for rk in REGION_KEYS:
                mine = effective(inputs, base, rk)
                theirs = effective(produced_inputs, base, rk)
                if mine != theirs:
                    diffs.append((rk, theirs, mine, rk in allowed))
            rows.append((name, stage, key, producer, where, diffs))
            if any(not declared for _, _, _, declared in diffs):
                differing.append(name)

    print("Region consistency: does each section run the region of the section "
          "that produced its structure?")
    print("Compared as each run resolves them (section value, else the base "
          "`qmmm.inputs`): %s" % ", ".join(REGION_KEYS))
    print()
    for name, stage, key, producer, where, diffs in rows:
        if not diffs:
            print("  [OK]    %-12s %-6s %-14s <- %s" % (name, stage, key, producer))
            continue
        print("  [DIFFERS] %-12s %-6s %-14s <- %s   (producer read from %s)"
              % (name, stage, key, producer, where))
        for rk, theirs, mine, declared in diffs:
            print("        %-18s %s: %s" % (rk, producer, describe(theirs)))
            print("        %-18s %s: %s%s"
                  % ("", name, describe(mine),
                     "   [declared in require.same_as.except]" if declared else ""))

    if unknown_producer:
        print()
        print("  Not checked - the structure was not written by a section "
              "(a crystal, an md_* result, or a hand-built file):")
        for name, key, path in unknown_producer:
            print("    %-12s %-14s %s" % (name, key, path))

    print()
    if differing:
        print("UNDECLARED DIFFERENCES in: %s" % ", ".join(sorted(set(differing))))
        print("A difference is not automatically wrong - `attackmin` changed "
              "`qm_residues` deliberately. But each one is a structure optimised "
              "with one movable set and consumed by a calculation using another, "
              "which is D207/D317/D320's defect, so each needs a reason written "
              "down. Declare an intended one in `require.same_as.except`.")
    else:
        print("No undeclared differences.")
    return 1 if (args.strict and differing) else 0


if __name__ == "__main__":
    sys.exit(main())
