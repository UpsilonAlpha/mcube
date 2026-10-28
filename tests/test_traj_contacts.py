"""The nucleophile residue is named differently in the two protonation arms.

`make_water_nuc.py` renames the hydroxide OH to NUW when it protonates it, so
that OpenMM matches it against `water_nuc.xml` instead of against every bulk
water. `traj_contacts.py` has to measure both arms with one probe set, and an
empty selection must never be reported as a missing contact - that is the
mistake D107 had to correct by hand.
"""

import importlib.util
import sys
import types
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "docs" / "csr" / "traj_contacts.py"


def load():
    spec = importlib.util.spec_from_file_location("traj_contacts", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def topology(*resnames):
    """A stand-in with just the `.residues` the resolver looks at."""
    return types.SimpleNamespace(
        residues=[types.SimpleNamespace(name=name) for name in resnames])


def test_the_hydroxide_arm_is_found_by_name():
    assert load().nucleophile_resname(topology("ALA", "OH", "HOH")) == "OH"


def test_the_water_arm_is_found_by_name():
    assert load().nucleophile_resname(topology("ALA", "NUW", "HOH")) == "NUW"


def test_a_topology_with_neither_is_refused():
    with pytest.raises(ValueError, match="no nucleophile residue"):
        load().nucleophile_resname(topology("ALA", "HOH", "ZN"))


def test_a_topology_with_both_is_refused():
    # Two nucleophiles would make every distance a silent minimum over them.
    with pytest.raises(ValueError, match="more than one"):
        load().nucleophile_resname(topology("OH", "NUW"))


def test_bulk_water_is_not_mistaken_for_the_nucleophile():
    with pytest.raises(ValueError, match="no nucleophile residue"):
        load().nucleophile_resname(topology(*["HOH"] * 50))
