#!/usr/bin/env python3
"""argmap - validate a directory of argument-map notes.

Project-agnostic. Standard library only. Run it as:

    python argmap.py [directory]          # default: the current directory

It discovers the thesis from frontmatter (`type: thesis`), builds the edge
graph from the body sections, and runs five checks. Optional per-project
settings live in `argmap.json` beside the notes; everything in it is optional
and every check that depends on it degrades to "skipped" rather than failing.

See SPEC.md for the file format and AGENT.md for how to use it.

Exit status: 0 when nothing is flagged, 1 when something is, 2 on a usage or
structural error (no notes, no thesis).
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

# --------------------------------------------------------------------------
# the vocabulary: four edge types, four node types, three evidence classes
# --------------------------------------------------------------------------

#  heading (lowercased, no trailing colon) -> (relation, direction)
#  "out" = this note asserts the edge; "in" = this note mirrors an edge
#  another note asserts. Both are parsed; the graph is one directed edge set.
HEADINGS = {
    "supports": ("supports", "out"),
    "refutes": ("refutes", "out"),
    "requires": ("requires", "out"),
    "superseded by": ("superseded-by", "out"),
    "superseded-by": ("superseded-by", "out"),
    "supported by": ("supports", "in"),
    "refuted by": ("refutes", "in"),
    "required by": ("requires", "in"),
    "supersedes": ("superseded-by", "in"),
}

RELATIONS = ("supports", "refutes", "requires", "superseded-by")
NODE_TYPES = ("thesis", "claim", "observation", "question")
EVIDENCE_CLASSES = ("experimental", "literature", "ours")

# A node counts as asserted-and-standing for the purposes of check 2.
SUPPORTED = ("supported",)

DEFAULTS = {
    # file listing sources this project has overturned; checked by check 4
    "retracted_sources": None,
    # headings in that file to ignore (substring match, lowercased)
    "retracted_sources_ignore_headings": ["already marked"],
    # how a cited source id looks, e.g. "D\\d+" for a decision log, or
    # "[0-9]{4}-[0-9]{2}-[0-9]{2}" for dated journal entries
    "source_id_pattern": r"[A-Za-z]+\d+",
    # frontmatter key holding the citations
    "source_field": "entries",
    # frontmatter key listing overturned sources whose superseding source the
    # author has read. Check 4 stops flagging a citation once it is listed
    # here, so the list is an assertion: "I read the correction."
    "acknowledged_field": "checked",
    # free-text pointer, printed in the header so a reader knows the journal
    "journal": None,
}


# --------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------

WIKILINK = re.compile(r"\[\[([^\]|#]+?)(?:[|#][^\]]*)?\]\]")
INLINE_LIST = re.compile(r"^\[(.*)\]$")


def parse_frontmatter(text: str) -> tuple[dict, str]:
    """Return (frontmatter, body). Supports the small YAML subset used here:
    `key: scalar`, `key: [a, b]`, and a `-` block list under a bare `key:`."""
    if not text.startswith("---"):
        return {}, text
    lines = text.splitlines()
    end = None
    for i, line in enumerate(lines[1:], start=1):
        if line.strip() == "---":
            end = i
            break
    if end is None:
        return {}, text
    meta: dict = {}
    key = None
    for raw in lines[1:end]:
        line = raw.rstrip()
        if not line.strip():
            continue
        if line.lstrip().startswith("-") and key is not None:
            meta.setdefault(key, [])
            if isinstance(meta[key], list):
                meta[key].append(line.lstrip()[1:].strip().strip("\"'"))
            continue
        if ":" not in line:
            continue
        key, _, value = line.partition(":")
        key = key.strip()
        value = value.strip()
        if not value:
            meta[key] = []
            continue
        m = INLINE_LIST.match(value)
        if m:
            inner = m.group(1).strip()
            meta[key] = [
                v.strip().strip("\"'") for v in inner.split(",") if v.strip()
            ] if inner else []
        else:
            meta[key] = value.strip("\"'")
    return meta, "\n".join(lines[end + 1:])


def parse_edges_directed(body: str) -> list[tuple[str, str, str]]:
    """Return [(relation, direction, target_slug)]."""
    out: list[tuple[str, str, str]] = []
    current: tuple[str, str] | None = None
    for line in body.splitlines():
        stripped = line.strip()
        if stripped.startswith("#"):
            name = stripped.lstrip("#").strip().rstrip(":").lower()
            current = HEADINGS.get(name)
            continue
        if current is None or not stripped.startswith(("-", "*")):
            continue
        for target in WIKILINK.findall(stripped):
            out.append((current[0], current[1], target.strip()))
    return out


class Node:
    def __init__(self, path: Path, meta: dict, body: str):
        self.path = path
        self.slug = path.stem
        self.meta = meta
        self.body = body
        self.type = str(meta.get("type", "")).strip().lower()
        self.status = str(meta.get("status", "")).strip().lower()
        self.evidence = str(meta.get("evidence", "")).strip().lower() or None

    def sources(self, field: str) -> list[str]:
        value = self.meta.get(field, [])
        if isinstance(value, str):
            return [v.strip() for v in value.split(",") if v.strip()]
        return [str(v).strip() for v in value if str(v).strip()]


def load_config(directory: Path) -> dict:
    config = dict(DEFAULTS)
    path = directory / "argmap.json"
    if path.exists():
        try:
            config.update(json.loads(path.read_text(encoding="utf-8")))
        except (OSError, ValueError) as exc:
            print(f"WARNING: could not read {path.name}: {exc}", file=sys.stderr)
    return config


def load_nodes(directory: Path) -> tuple[dict[str, Node], list[str]]:
    nodes: dict[str, Node] = {}
    problems: list[str] = []
    skip = {"readme", "spec", "agent", "index"}
    for path in sorted(directory.glob("*.md")):
        if path.stem.lower() in skip:
            continue
        text = path.read_text(encoding="utf-8")
        meta, body = parse_frontmatter(text)
        if not meta:
            problems.append(f"{path.name}: no frontmatter; not treated as a node")
            continue
        node = Node(path, meta, body)
        if node.type not in NODE_TYPES:
            problems.append(
                f"{path.name}: type {node.type!r} is not one of "
                f"{', '.join(NODE_TYPES)}"
            )
        if node.type == "observation" and node.evidence not in EVIDENCE_CLASSES:
            problems.append(
                f"{path.name}: observation needs evidence: one of "
                f"{', '.join(EVIDENCE_CLASSES)} (got {node.evidence!r})"
            )
        if node.type != "observation" and node.evidence:
            problems.append(
                f"{path.name}: evidence class on a {node.type}; evidence belongs "
                f"on observations"
            )
        nodes[node.slug] = node
    return nodes, problems


# --------------------------------------------------------------------------
# retracted-source index (optional, project-configured)
# --------------------------------------------------------------------------

def load_retracted(path: Path, pattern: str, ignore: list[str]) -> set[str]:
    """Extract source ids from a retractions file.

    Generic over format: it reads ids out of the first cell of any table row
    and out of any `- **ID**` bullet, which is how such an index is normally
    written. Headings matching `ignore` are skipped.
    """
    rx = re.compile(pattern)
    bullet = re.compile(r"^\s*[-*]\s+\**\s*(" + pattern + r")")
    found: set[str] = set()
    skipping = False
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("#"):
            heading = stripped.lstrip("#").strip().lower()
            skipping = any(token in heading for token in ignore)
            continue
        if skipping:
            continue
        if stripped.startswith("|"):
            cells = [c.strip() for c in stripped.strip("|").split("|")]
            if not cells or set(cells[0]) <= set("-: "):
                continue
            found.update(rx.findall(cells[0]))
            continue
        m = bullet.match(stripped)
        if m:
            found.add(m.group(1))
    return found


# --------------------------------------------------------------------------
# the graph
# --------------------------------------------------------------------------

class Graph:
    def __init__(self, nodes: dict[str, Node]):
        self.nodes = nodes
        self.edges: set[tuple[str, str, str]] = set()   # (src, rel, dst)
        self.broken: list[tuple[str, str, str]] = []    # (src, rel, missing)
        for slug, node in nodes.items():
            for rel, direction, target in parse_edges_directed(node.body):
                src, dst = (slug, target) if direction == "out" else (target, slug)
                if target not in nodes:
                    self.broken.append((slug, rel, target))
                    continue
                self.edges.add((src, rel, dst))

    def inbound(self, slug: str, rel: str) -> list[str]:
        return sorted(s for (s, r, d) in self.edges if d == slug and r == rel)

    def outbound(self, slug: str, rel: str) -> list[str]:
        return sorted(d for (s, r, d) in self.edges if s == slug and r == rel)

    def connected_to(self, roots: list[str]) -> set[str]:
        """Undirected reachability: a node participates in the argument if it
        is linked to the thesis by any chain of any of the four edge types."""
        adjacency: dict[str, set[str]] = {s: set() for s in self.nodes}
        for src, _rel, dst in self.edges:
            adjacency.setdefault(src, set()).add(dst)
            adjacency.setdefault(dst, set()).add(src)
        seen = set(roots)
        stack = list(roots)
        while stack:
            for neighbour in adjacency.get(stack.pop(), ()):
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        return seen


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------

class Report:
    def __init__(self) -> None:
        self.findings = 0
        self.skipped: list[str] = []

    def section(self, number: int | str, title: str, items: list[str],
                clean: str = "nothing flagged") -> None:
        print(f"\n[{number}] {title}")
        if not items:
            print(f"    {clean}")
            return
        self.findings += len(items)
        for item in items:
            print(f"    - {item}")

    def skip(self, number: int | str, title: str, why: str) -> None:
        print(f"\n[{number}] {title}")
        print(f"    skipped: {why}")
        self.skipped.append(f"{number} {title}")


def main(argv: list[str]) -> int:
    directory = Path(argv[1] if len(argv) > 1 else ".").resolve()
    if not directory.is_dir():
        print(f"not a directory: {directory}", file=sys.stderr)
        return 2

    config = load_config(directory)
    nodes, problems = load_nodes(directory)
    if not nodes:
        print(f"no argument-map notes found in {directory}", file=sys.stderr)
        return 2

    graph = Graph(nodes)
    theses = [s for s, n in nodes.items() if n.type == "thesis"]

    counts = {t: sum(1 for n in nodes.values() if n.type == t) for t in NODE_TYPES}
    evidence = {
        c: sum(1 for n in nodes.values() if n.evidence == c)
        for c in EVIDENCE_CLASSES
    }

    print(f"argmap: {directory}")
    if config.get("journal"):
        print(f"journal: {config['journal']}  (the graph holds what we believe; "
              f"the journal holds what happened)")
    print("nodes:   " + ", ".join(f"{k} {v}" for k, v in counts.items()))
    print("evidence:" + " " + ", ".join(f"{k} {v}" for k, v in evidence.items()))
    print("edges:   " + ", ".join(
        f"{r} {sum(1 for e in graph.edges if e[1] == r)}" for r in RELATIONS))

    report = Report()

    if problems:
        report.section("0", "Note format", problems)
    else:
        print("\n[0] Note format\n    all notes parse, types and evidence "
              "classes valid")

    if not theses:
        print("\nERROR: no node has `type: thesis`; nothing to check against.",
              file=sys.stderr)
        return 2
    if len(theses) > 1:
        report.section("0b", "More than one thesis",
                       [f"{s} (a map has one root)" for s in sorted(theses)])

    # --- 1 -----------------------------------------------------------------
    items = []
    for slug, node in sorted(nodes.items()):
        if node.type not in ("claim", "thesis") or node.status not in SUPPORTED:
            continue
        supporters = graph.inbound(slug, "supports")
        obs = [s for s in supporters if nodes[s].type == "observation"]
        classes = {nodes[s].evidence for s in obs}
        if obs and classes == {"ours"}:
            items.append(
                f"{slug}: marked supported, and its only supporting "
                f"observations are our own calculations ({', '.join(obs)})"
            )
        elif supporters and not obs:
            items.append(
                f"{slug}: marked supported with no observation supporting it "
                f"directly (supporters: {', '.join(supporters)})"
            )
        elif not supporters:
            items.append(f"{slug}: marked supported with an empty support set")
    report.section(1, "Claims held on our own calculations alone", items)

    # --- 2 -----------------------------------------------------------------
    items = []
    for slug, node in sorted(nodes.items()):
        if node.status not in SUPPORTED:
            continue
        for refuter in graph.inbound(slug, "refutes"):
            if nodes[refuter].status in SUPPORTED:
                ev = nodes[refuter].evidence
                items.append(
                    f"{slug}: still marked supported, but {refuter} "
                    f"({nodes[refuter].type}"
                    f"{', ' + ev if ev else ''}, supported) refutes it"
                )
    report.section(2, "Supported claims with a standing refutation", items)

    # --- 3 -----------------------------------------------------------------
    items = []
    for thesis in sorted(theses):
        for premise in graph.outbound(thesis, "requires"):
            if not graph.inbound(premise, "supports"):
                items.append(
                    f"{premise}: required by {thesis}, support set empty "
                    f"(status {nodes[premise].status or 'unset'})"
                )
    report.section(3, "Unsupported premises of the thesis", items)

    # --- 4 -----------------------------------------------------------------
    source_file = config.get("retracted_sources")
    field = config.get("source_field", "entries")
    if not source_file:
        report.skip(4, "Citations of overturned sources",
                    "no `retracted_sources` configured in argmap.json")
    else:
        path = (directory / source_file).resolve()
        if not path.exists():
            report.skip(4, "Citations of overturned sources",
                        f"{source_file} not found")
        else:
            retracted = load_retracted(
                path,
                config.get("source_id_pattern", DEFAULTS["source_id_pattern"]),
                [t.lower() for t in config.get(
                    "retracted_sources_ignore_headings",
                    DEFAULTS["retracted_sources_ignore_headings"])],
            )
            ack_field = config.get("acknowledged_field", "checked")
            items, notes = [], []
            for slug, node in sorted(nodes.items()):
                acknowledged = set(node.sources(ack_field))
                hits = sorted(
                    s for s in node.sources(field)
                    if s in retracted and s not in acknowledged
                )
                if not hits:
                    continue
                line = f"{slug}: cites {', '.join(hits)}"
                if node.status in SUPPORTED:
                    items.append(line + " - and is marked supported")
                else:
                    notes.append(line)
            report.section(
                4,
                f"Supported nodes citing overturned sources without "
                f"acknowledgement ({len(retracted)} sources listed in "
                f"{source_file})",
                items,
            )
            if notes:
                print(f"    note: {len(notes)} further node(s) cite an "
                      f"overturned source; not marked supported, so not "
                      f"counted:")
                for line in notes:
                    print(f"      . {line}")

    # --- 5 -----------------------------------------------------------------
    items = [f"{src} -> [[{target}]] ({rel}): no such note"
             for src, rel, target in sorted(graph.broken)]
    report.section("5a", "Broken wikilinks", items)

    reachable = graph.connected_to(theses)
    items = [f"{slug} ({nodes[slug].type}): no edge path to the thesis"
             for slug in sorted(set(nodes) - reachable)]
    report.section("5b", "Nodes with no path to the thesis", items)

    print()
    if report.findings:
        print(f"{report.findings} finding(s).", end=" ")
    else:
        print("No findings.", end=" ")
    if report.skipped:
        print(f"Skipped: {'; '.join(report.skipped)}.", end=" ")
    print()
    return 1 if report.findings else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
