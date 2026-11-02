#!/usr/bin/env python3
"""Make a config's force-field paths relative to the config, for the new rule.

    python migrate_paths.py <run dir>/config.yaml [more ...]     # dry run
    python migrate_paths.py --apply <run dir>/config.yaml [...]

Force-field entries are now read against the config's own directory, so a run
directory is self-contained and a copy of it uses its own files. Entries
written for the old behaviour do not survive that change unaided:

  ../hydroxide.xml                      meant <run dir>/hydroxide.xml only
                                        because the stage ran one level down.
                                        Now it points outside the run; rewritten.

  /scratch/.../<this run>/hydroxide.xml still works, but ties the directory to
                                        its own path, so a copy reaches back
                                        here. Rewritten to a bare name.

  /scratch/.../<another run>/x.xml      REPORTED, NEVER REWRITTEN. Rewriting it
                                        would change which file the run uses,
                                        which is a decision about the science
                                        and not a path fix. Three CSR configs
                                        are in this state, from copying.

  amber14/protein.ff14SB.xml            left alone; OpenMM resolves its own.

Only entries that are force-field files are touched: a list item under a
`forcefields` or `default_forcefields` key ending in `.xml`. Nothing else in
the config is read or rewritten, and a dry run prints every change first.
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

KEYS = ("forcefields", "default_forcefields")
ENTRY = re.compile(r"^(\s*-\s*)(\S+\.xml)\s*$")
REPOINT = False
KEYLINE = re.compile(r"^\s*(\w+)\s*:\s*$")



def same_but_for_a_date(a: Path, b: Path) -> bool:
    """True when two force-field files differ only in a <DateGenerated> line.

    OpenFF writes the generation date into its output, so two runs of the same
    preprocess on the same SMILES produce files that differ in one line and
    nothing else. That is the only difference it is safe to ignore: repointing
    a config at a file whose *parameters* differ changes the force field the
    run used, which for a finished run means the config stops describing it.
    """
    try:
        left, right = a.read_text().splitlines(), b.read_text().splitlines()
    except OSError:
        return False
    if len(left) != len(right):
        return False
    differing = [(x, y) for x, y in zip(left, right) if x != y]
    return all("DateGenerated" in x and "DateGenerated" in y
               for x, y in differing)


def repoint(value: str, home: Path):
    """A foreign absolute path replaced by this run's own copy, if that is safe.

    Returns (new value, why) or (None, why not). Five CSR configs read force
    fields out of another run's directory, inherited by copying. Pointing one
    at its own file is only a no-op when the file is there and holds the same
    parameters; otherwise it is a change to the science wearing the clothes of
    a path fix, so it is refused and reported.
    """
    foreign = Path(value)
    rel = None
    for parent in (foreign.parent, foreign.parent.parent):
        try:
            rel = foreign.relative_to(parent.parent)
            break
        except ValueError:
            continue
    local = home / foreign.name
    nested = home / "preprocess_results" / foreign.name
    own = local if local.is_file() else (nested if nested.is_file() else None)
    if own is None:
        return None, f"no {foreign.name} in this run"
    if own.read_bytes() == foreign.read_bytes():
        return own.relative_to(home).as_posix(), "identical"
    if same_but_for_a_date(own, foreign):
        return own.relative_to(home).as_posix(), "same but for the date"
    return None, "its own copy has DIFFERENT contents"


def classify(value: str, home: Path):
    """(action, new value, why) for one entry."""
    if value.startswith("amber14/") or value.startswith("implicit/"):
        return "keep", value, "OpenMM's own"
    path = Path(value)
    if path.is_absolute():
        resolved = path.resolve()
        if home in resolved.parents:
            return "rewrite", resolved.relative_to(home).as_posix(), "absolute, this run"
        if REPOINT:
            own, why = repoint(value, home)
            if own:
                return "rewrite", own, f"repointed at this run's own ({why})"
            return "foreign", value, f"NOT repointed: {why}"
        return "foreign", value, f"absolute, ANOTHER run ({resolved.parent})"
    if value.startswith(".."):
        # `../hydroxide.xml` meant <run dir>/hydroxide.xml, because the stage
        # ran one level down in md_<section>_results/ or two in protdft/<x>/.
        # So the faithful migration strips the prefix and looks beside the
        # config - NOT at home/".."/value, which is the directory the runs sit
        # in and a different file entirely.
        stripped = re.sub(r"^(\.\./)+", "", value)
        if (home / stripped).is_file():
            return "rewrite", stripped, "../ that meant a file in this run"
        return "stale", value, f"../ and no {stripped} beside the config"
    if (home / value).is_file():
        return "keep", value, "already relative to the config"
    return "missing", value, f"no such file under {home}"


def migrate(config: Path, apply: bool) -> int:
    home = config.resolve().parent
    lines = config.read_text().splitlines(keepends=True)
    in_list, changes, notes = False, [], []
    for index, line in enumerate(lines):
        key = KEYLINE.match(line)
        if key:
            in_list = key.group(1) in KEYS
            continue
        entry = ENTRY.match(line)
        if not entry:
            if line.strip() and not line.lstrip().startswith("-"):
                in_list = False
            continue
        if not in_list:
            continue
        action, value, why = classify(entry.group(2), home)
        if action == "rewrite":
            changes.append((index, entry.group(1) + value + "\n",
                            f"  line {index + 1}: {entry.group(2)} -> {value}"))
        elif action in ("foreign", "stale", "missing"):
            notes.append(f"  line {index + 1}: {entry.group(2)}  [{action}: {why}]")

    print(f"{config}")
    for _, _, text in changes:
        print(text)
    for text in notes:
        print(text)
    if not changes and not notes:
        print("  nothing to do")

    if changes and apply:
        shutil.copy2(config, config.with_suffix(".yaml.bak_paths"))
        for index, text, _ in changes:
            lines[index] = text
        config.write_text("".join(lines))
        print(f"  WROTE {len(changes)} change(s); backup .bak_paths")
    elif changes:
        print(f"  dry run: {len(changes)} change(s) not written")
    return len(changes)


def main(argv):
    global REPOINT
    apply = "--apply" in argv
    REPOINT = "--repoint" in argv
    argv = [a for a in argv if a not in ("--apply", "--repoint")]
    if not argv:
        raise SystemExit(__doc__)
    total = 0
    for name in argv:
        config = Path(name)
        if not config.is_file():
            raise SystemExit(f"no such config: {config}")
        total += migrate(config, apply)
    print(f"\n{total} change(s) {'written' if apply else 'pending'}")


if __name__ == "__main__":
    main(sys.argv[1:])
