"""QM level resolution: pytest tests/test_qm_levels.py

Needs the ASH environment, because mcube.qmmm imports ASH; it is skipped elsewhere.
"""

import pytest

pytest.importorskip("ash", reason="mcube.qmmm imports ASH")

from mcube.qmmm import level_option, qm_level_inputs  # noqa: E402


class TestLevelOption:
    def test_prefixed_key_wins(self):
        inputs = {"dft_basis": "def2-svp", "sp_basis": "def2-tzvpd"}
        assert level_option(inputs, "dft", "basis") == "def2-svp"

    def test_falls_back_to_sp(self):
        """A config that defines the single-point level need not repeat it."""
        assert level_option({"sp_basis": "def2-tzvpd"}, "dft", "basis") == "def2-tzvpd"

    def test_default_when_neither_is_set(self):
        assert level_option({}, "dft", "basis", "def2-svp") == "def2-svp"

    def test_falsy_values_are_respected(self):
        """densityfit: false must not fall through to the default."""
        assert level_option({"dft_densityfit": False}, "dft", "densityfit", True) is False
        assert level_option({"sp_densityfit": False}, "dft", "densityfit", True) is False

    def test_sp_prefix_reads_its_own_keys(self):
        assert level_option({"sp_functional": "wb97x"}, "sp", "functional") == "wb97x"


class TestQmLevelInputs:
    """A stage section may override the QM level, but nothing else about the QM/MM setup.

    scan/neb/ts take the QM theory from qmmm.inputs, so without this a
    `qm_engine: dft` on one of those sections would be silently ignored and the
    stage would quietly run at xTB.
    """

    base = {"charge": 2, "mult": 1, "cores": 4, "qm_engine": "xtb",
            "sp_functional": "wb97x", "forcefields": ["a.xml"]}

    def test_engine_and_dft_keys_are_overridden(self):
        merged = qm_level_inputs(self.base, {"qm_engine": "dft", "dft_functional": "r2scan"})
        assert merged["qm_engine"] == "dft"
        assert merged["dft_functional"] == "r2scan"

    def test_sp_keys_are_overridden(self):
        merged = qm_level_inputs(self.base, {"sp_basis": "def2-tzvpd"})
        assert merged["sp_basis"] == "def2-tzvpd"

    def test_xtb_keys_are_overridden(self):
        merged = qm_level_inputs(self.base, {"xtb_method": "GFN1"})
        assert merged["xtb_method"] == "GFN1"

    def test_other_keys_are_not_taken_from_the_stage(self):
        """charge, mult, cores and force fields stay with qmmm.inputs."""
        merged = qm_level_inputs(self.base, {"charge": 99, "cores": 64, "atoms": [1, 2],
                                             "forcefields": ["wrong.xml"]})
        assert merged["charge"] == 2
        assert merged["cores"] == 4
        assert merged["forcefields"] == ["a.xml"]
        assert "atoms" not in merged

    def test_qmmm_inputs_are_not_mutated(self):
        qm_level_inputs(self.base, {"qm_engine": "dft"})
        assert self.base["qm_engine"] == "xtb"

    def test_no_stage_inputs(self):
        assert qm_level_inputs(self.base, None) == self.base
        assert qm_level_inputs(self.base, {}) == self.base
