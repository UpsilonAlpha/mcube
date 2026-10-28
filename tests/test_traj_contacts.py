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


class TestSkipMissing:
    """Older systems carry earlier ligand variants, so a probe can name an atom
    the topology genuinely does not have. That must stay an error by default -
    an empty selection reported as a missing contact is the whole reason
    `atoms()` raises - and be skippable only when asked for explicitly."""

    def probe_topology(self, *extra_resnames):
        """A topology complete enough for probes() to resolve, minus whatever
        the caller leaves out."""
        import types

        def atom(name, res):
            return types.SimpleNamespace(
                name=name, residue=res,
                element=types.SimpleNamespace(symbol=name[0]))

        residues, atoms_ = [], []
        spec = [("LYS", 226, ["NZ"]), ("ARG", 90, ["NE", "NH1", "NH2"]),
                ("ARG", 229, ["NE", "NH1", "NH2"]), ("GLN", 93, ["OE1", "NE2"]),
                ("ASP", 92, ["OD1", "OD2"]), ("OH", 2, ["O01"]),
                ("ZN", 1, ["ZN"]),
                ("mer", 1, ["C4x", "N1x", "O1x", "O2x", "O3x"] + list(extra_resnames))]
        for name, seq, names in spec:
            res = types.SimpleNamespace(name=name, resSeq=seq, atoms=[])
            for n in names:
                a = atom(n, res)
                res.atoms.append(a)
                atoms_.append(a)
            residues.append(res)
        for i, a in enumerate(atoms_):
            a.index = i
        return types.SimpleNamespace(residues=residues, atoms=atoms_)

    def test_a_complete_ligand_needs_no_flag(self):
        mod = load()
        assert mod.probes(self.probe_topology("O4x"))

    def test_a_missing_atom_is_an_error_by_default(self):
        mod = load()
        with pytest.raises(ValueError, match="O4x"):
            mod.probes(self.probe_topology())

    def test_skip_missing_drops_only_that_probe(self, capsys):
        mod = load()
        full = mod.probes(self.probe_topology("O4x"))
        partial = mod.probes(self.probe_topology(), skip_missing=True)
        assert 0 < len(partial) < len(full)
        assert "[SKIP]" in capsys.readouterr().out
        # The probes that name O4x and nothing else are gone; the group probes
        # that merely include it survive, measuring the members that exist.
        assert not [l for l, _, _ in partial if l.endswith("mer O4x")]

    def test_a_partial_group_selection_says_so(self, capsys):
        # "carboxylate O3x/O4x" against a ligand with no O4x measures a
        # different quantity than the label claims, which has to be stated.
        mod = load()
        mod.probes(self.probe_topology(), skip_missing=True)
        out = capsys.readouterr().out
        assert "[NOTE]" in out and "O4x" in out
    def test_the_orientation_probe_is_present(self):
        mod = load()
        labels = [l for l, _, _ in mod.probes(self.probe_topology("O4x"))]
        assert any("which O faces the metal" in l for l in labels)
        assert any("O1x" in l for l in labels)
