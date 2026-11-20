"""Charge/mult resolution for scan, neb and ts: pytest tests/test_charge_resolution.py

Needs the ASH environment, because mcube.qmmm imports ASH; skipped elsewhere.
"""

import pytest

pytest.importorskip("ash", reason="mcube.qmmm imports ASH")

from mcube.qmmm import MULT_DEFAULT, charge_and_mult  # noqa: E402


class TestChargeAndMult:
    """The base `qmmm:` section holds no charge on purpose (¶152), so a stage
    section has to be able to state its own. Job 33244793 died in 7.8 s
    because a scan could not (¶290).
    """

    def test_section_charge_is_used(self):
        assert charge_and_mult("scan_posepull", {"charge": 2}, {"mult": 1}) == (2, 1)

    def test_section_overrides_the_base(self):
        assert charge_and_mult("s", {"charge": 2}, {"charge": 0})[0] == 2

    def test_base_is_the_fallback(self):
        assert charge_and_mult("s", {}, {"charge": 0, "mult": 1}) == (0, 1)

    def test_zero_charge_is_not_treated_as_absent(self):
        """charge 0 is a real charge; `or` would make it fall through, which is
        why this tests against None explicitly."""
        assert charge_and_mult("s", {"charge": 0}, {"charge": 2})[0] == 0

    def test_zero_charge_in_the_base_is_not_treated_as_absent(self):
        assert charge_and_mult("s", {}, {"charge": 0})[0] == 0

    def test_missing_charge_raises_naming_the_section_and_the_key(self):
        with pytest.raises(KeyError) as excinfo:
            charge_and_mult("scan_posepull", {}, {"mult": 1})
        message = str(excinfo.value)
        assert "scan_posepull" in message and "charge" in message

    def test_missing_charge_says_there_is_no_default(self):
        """The message has to rule a default out, or the next reader adds one."""
        with pytest.raises(KeyError) as excinfo:
            charge_and_mult("s", {}, {})
        assert "no default" in str(excinfo.value)

    def test_mult_defaults_rather_than_raising(self):
        """Unlike charge, mult has a sensible default: 1."""
        assert charge_and_mult("s", {"charge": 2}, {}) == (2, MULT_DEFAULT)

    def test_section_mult_overrides_the_base(self):
        assert charge_and_mult("s", {"charge": 0, "mult": 3}, {"mult": 1})[1] == 3

    def test_the_resolved_pair_and_its_source_are_printed(self, capsys):
        """A resolution nothing reports is a resolution nobody can check."""
        charge_and_mult("scan_posepull", {"charge": 2}, {"mult": 1})
        out = capsys.readouterr().out
        assert "charge=2 (from scan_posepull.inputs)" in out
        assert "mult=1 (from qmmm.inputs)" in out

    def test_a_defaulted_mult_names_the_default_as_its_source(self, capsys):
        charge_and_mult("s", {"charge": 2}, {})
        assert "mult=1 (from default)" in capsys.readouterr().out


class TestEveryStageUsesTheOneResolver:
    """scan, neb and ts had three different behaviours for the same config:
    scan passed charge=None to the optimizer, neb raised a bare KeyError, and
    ts did what scan did. ¶291 is that there is now one answer, so this asserts
    the identity rather than trusting that three call sites agree.
    """

    def test_scan_neb_and_ts_share_it(self):
        import mcube.neb
        import mcube.scan
        import mcube.ts
        from mcube.qmmm import charge_and_mult as canonical

        assert mcube.scan.charge_and_mult is canonical
        assert mcube.neb.charge_and_mult is canonical
        assert mcube.ts.charge_and_mult is canonical
