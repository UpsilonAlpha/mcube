"""preprocess refuses to turn a restraint into a constraint: pytest tests/test_provenance.py

Needs the ASH environment, because mcube.preprocess imports ASH.
"""

import pytest

pytest.importorskip("ash", reason="mcube.preprocess imports ASH")

from mcube.preprocess import check_inherited_restraints  # noqa: E402

RESTRAINTS = [[4203, 4209, 2.2, 20.0]]          # md pulled Zn-O2x to 2.2 A
DETECTED = [[4203, 4209], [4203, 1349]]         # preprocess now sees both as contacts


def test_refuses_and_says_which_pair_and_how_to_resolve_it():
    with pytest.raises(ValueError) as excinfo:
        check_inherited_restraints(DETECTED, RESTRAINTS, {})
    message = str(excinfo.value)
    assert "4203 - 4209" in message and "2.2" in message
    assert "constraint_exclude_residues" in message and "acknowledge_restrained" in message


def test_pair_order_does_not_matter():
    with pytest.raises(ValueError):
        check_inherited_restraints([[4209, 4203]], RESTRAINTS, {})


def test_acknowledging_the_pair_allows_it():
    check_inherited_restraints(DETECTED, RESTRAINTS, {"acknowledge_restrained": [[4203, 4209]]})


@pytest.mark.parametrize("constraints, restraints", [
    ([[4203, 1349]], RESTRAINTS),   # a contact that was never restrained
    (DETECTED, []),                 # md imposed nothing
    (DETECTED, None),               # md has not run
    (None, RESTRAINTS),             # no metals, so nothing detected
])
def test_leaves_everything_else_alone(constraints, restraints):
    check_inherited_restraints(constraints, restraints, {})
