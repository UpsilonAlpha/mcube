"""Scan charge resolution: pytest tests/test_scan_charge.py

Needs the ASH environment, because mcube.scan imports ASH; skipped elsewhere.
"""

import pytest

pytest.importorskip("ash", reason="mcube.scan imports ASH")

from mcube.scan import charge_and_mult  # noqa: E402


class TestChargeAndMult:
    """The base `qmmm:` section holds no charge on purpose (¶152), so a scan
    section has to be able to state its own. Job 33244793 died in 7.8 s
    because it could not (¶290).
    """

    def test_section_charge_is_used(self):
        assert charge_and_mult("scan_posepull", {"charge": 2},
                               {"mult": 1}) == (2, 1)

    def test_section_overrides_the_base(self):
        assert charge_and_mult("s", {"charge": 2}, {"charge": 0})[0] == 2

    def test_base_is_the_fallback(self):
        assert charge_and_mult("s", {}, {"charge": 0, "mult": 1}) == (0, 1)

    def test_zero_charge_is_not_treated_as_absent(self):
        """charge 0 is a real charge; `or` would have turned it into a
        fallthrough, which is the whole reason this is `.get` with a default."""
        assert charge_and_mult("s", {"charge": 0}, {"charge": 2})[0] == 0

    def test_missing_charge_raises_naming_the_section(self):
        with pytest.raises(KeyError) as e:
            charge_and_mult("scan_posepull", {}, {"mult": 1})
        assert "scan_posepull" in str(e.value)

    def test_mult_may_be_absent(self):
        """ASH defaults the multiplicity; only the charge is refused."""
        assert charge_and_mult("s", {"charge": 2}, {}) == (2, None)
