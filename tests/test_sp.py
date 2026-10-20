"""Input layering for the sp stage: pytest tests/test_sp.py

Needs the ASH environment, because mcube.sp imports ASH; it is skipped elsewhere.
"""

import textwrap

import pytest

pytest.importorskip("ash", reason="mcube.sp imports ASH")

from mcube.sp import sp_inputs  # noqa: E402


def write(tmp_path, body):
    path = tmp_path / "config.yaml"
    path.write_text(textwrap.dedent(body))
    return str(path)


class TestSpInputs:
    def test_inherits_the_qmmm_model(self, tmp_path):
        """Force fields, charge and multiplicity should not have to be repeated."""
        config = write(tmp_path, """
            qmmm:
              inputs:
                charge: 2
                mult: 1
                forcefields: [amber14/protein.ff14SB.xml]
            """)
        inputs = sp_inputs(config)
        assert inputs["charge"] == 2
        assert inputs["forcefields"] == ["amber14/protein.ff14SB.xml"]

    def test_sp_overrides_qmmm(self, tmp_path):
        config = write(tmp_path, """
            qmmm:
              inputs: {charge: 2, sp_cores: 4}
            sp:
              inputs: {sp_cores: 64}
            """)
        inputs = sp_inputs(config)
        assert inputs["sp_cores"] == 64
        assert inputs["charge"] == 2

    def test_section_overrides_sp(self, tmp_path):
        """A per-structure section is the most specific layer."""
        config = write(tmp_path, """
            qmmm:
              inputs: {sp_basis: def2-svp, charge: 2}
            sp:
              inputs: {sp_cores: 64}
            sp_product:
              inputs: {sp_basis: def2-tzvp}
            """)
        inputs = sp_inputs(config, "sp_product")
        assert inputs["sp_basis"] == "def2-tzvp"
        assert inputs["sp_cores"] == 64
        assert inputs["charge"] == 2

    def test_works_without_an_sp_section(self, tmp_path):
        """Configs written before this stage existed keep working."""
        config = write(tmp_path, """
            qmmm:
              inputs: {charge: 2, sp_functional: wb97x-v}
            """)
        assert sp_inputs(config)["sp_functional"] == "wb97x-v"

    def test_qmmm_section_name_is_not_applied_twice(self, tmp_path):
        """--section qmmm must not be treated as an extra override layer."""
        config = write(tmp_path, """
            qmmm:
              inputs: {sp_cores: 4}
            sp:
              inputs: {sp_cores: 64}
            """)
        assert sp_inputs(config, "qmmm")["sp_cores"] == 64

    def test_missing_sections_are_not_an_error(self, tmp_path):
        config = write(tmp_path, "preprocess:\n  inputs: {}\n")
        assert sp_inputs(config) == {}
