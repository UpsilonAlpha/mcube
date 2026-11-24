"""No stage may touch its `result` dict before it is bound: pytest tests/test_result_binding.py

D337: every `ts` run died in under a second with

    File "src/mcube/ts.py", line 119, in run
      result["qm_level"] = qm_level_label(qmmm_inputs)
    UnboundLocalError: cannot access local variable 'result' where it is not
    associated with a value

because two statements - `result["qm_level"] = ...` and
`result.update(actatoms_record(...))` - sat ABOVE the `result: Dict[str, Any] =
{...}` that creates it. The stage had never run, so the bug had never surfaced:
job 33410206 was the first `ts` submission the project ever made.

This is a static check rather than a `ts` test because that is what makes it
cheap and total. It needs no ASH, no OpenMM, no run directory and no structure,
so it runs on any machine; and it checks EVERY function in EVERY stage rather
than the one line that failed, which is the difference between a regression test
and a test for the defect class. A stage that reads `result` before binding it
cannot work, whatever else is true of it, and that is decidable from the source.

Deliberately not caught here: a `result` bound on only one branch of an `if`
(bind-on-some-paths). `ast` sees the assignment and this check passes. Catching
that needs flow analysis, which is a linter's job - `pyflakes`/`ruff F821` do
report the ts.py case, and this file is the cheap always-on part.
"""

from __future__ import annotations

import ast
import io
import pathlib

import pytest

SRC = pathlib.Path(__file__).resolve().parent.parent / "src" / "mcube"

# The name every stage's `run()` accumulates into and hands to `save_result`.
NAME = "result"


def _first_binding(fn: ast.FunctionDef) -> int | None:
    """The earliest line that assigns `NAME`, or None if it is never assigned."""
    lines = [
        node.lineno
        for node in ast.walk(fn)
        if isinstance(node, (ast.Assign, ast.AnnAssign))
        for target in (node.targets if isinstance(node, ast.Assign) else [node.target])
        if isinstance(target, ast.Name) and target.id == NAME
    ]
    return min(lines) if lines else None


def _reads_before(fn: ast.FunctionDef, bind: int) -> list[int]:
    """Lines that LOAD `NAME` before `bind` - a subscript or `.update` counts."""
    return sorted(
        {
            node.lineno
            for node in ast.walk(fn)
            if isinstance(node, ast.Name)
            and node.id == NAME
            and isinstance(node.ctx, ast.Load)
            and node.lineno < bind
        }
    )


def _functions():
    for path in sorted(SRC.glob("*.py")):
        tree = ast.parse(io.open(path, encoding="utf8").read(), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                yield path, node


def test_src_is_where_we_think_it_is():
    """A wrong SRC would make every test below pass by finding nothing."""
    assert (SRC / "ts.py").is_file(), SRC
    assert len(list(SRC.glob("*.py"))) > 5


@pytest.mark.parametrize(
    "path,fn",
    [(p, f) for p, f in _functions() if _first_binding(f) is not None],
    ids=lambda v: v.name if isinstance(v, ast.AST) else pathlib.Path(v).name,
)
def test_result_is_bound_before_it_is_read(path, fn):
    bind = _first_binding(fn)
    early = _reads_before(fn, bind)
    assert not early, (
        "%s:%s() binds `%s` at line %d but reads it at %s. "
        "Move those statements below the assignment - this is D337, and it "
        "raises UnboundLocalError on the first call."
        % (path.name, fn.name, NAME, bind, ", ".join(str(n) for n in early))
    )


def test_the_check_catches_the_shape_it_was_written_for():
    """The D337 shape itself, so a broken checker cannot pass silently."""
    tree = ast.parse(
        "def run(args):\n"
        "    result['qm_level'] = label()\n"
        "    result.update(record())\n"
        "    result = {'pdb_file': 'x'}\n"
        "    return result\n"
    )
    fn = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef))
    bind = _first_binding(fn)
    assert bind == 4
    assert _reads_before(fn, bind) == [2, 3]


def test_the_check_passes_the_corrected_shape():
    tree = ast.parse(
        "def run(args):\n"
        "    result = {'pdb_file': 'x'}\n"
        "    result['qm_level'] = label()\n"
        "    result.update(record())\n"
        "    return result\n"
    )
    fn = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef))
    assert _reads_before(fn, _first_binding(fn)) == []
