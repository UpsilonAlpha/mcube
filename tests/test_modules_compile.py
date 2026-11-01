"""Every module parses: pytest tests/test_modules_compile.py

Several mcube modules import heavy optional dependencies at module scope - ASH,
OpenMM, MDAnalysis - so the tests that exercise them are gated on
`importorskip`. That gating silently hides a *syntax* error: when the dependency
is missing the module is never imported, so the file is never parsed, and the
suite is green. A broken `print("\\n...")` in `postprocess.py` reached Bunya
that way, where it failed at collection in the one environment that has
MDAnalysis (D167).

`compile()` parses without importing, so this runs everywhere and needs nothing
installed. It catches only syntax, which is exactly the class of error the
gating hides.
"""

from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src" / "mcube"
MODULES = sorted(SRC.rglob("*.py"))
HELPERS = sorted((Path(__file__).resolve().parents[1] / "docs" / "csr").glob("*.py"))


def ids(paths):
    return [p.name for p in paths]


@pytest.mark.parametrize("path", MODULES, ids=ids(MODULES))
def test_the_module_parses(path):
    compile(path.read_text(encoding="utf-8"), str(path), "exec")


@pytest.mark.parametrize("path", HELPERS, ids=ids(HELPERS))
def test_the_helper_script_parses(path):
    """The docs/csr scripts are run by hand on Bunya, so nothing else parses
    them until the moment they are needed."""
    compile(path.read_text(encoding="utf-8"), str(path), "exec")


def test_there_is_something_to_check():
    """A glob that silently matches nothing would make this file pass forever."""
    assert len(MODULES) > 5
