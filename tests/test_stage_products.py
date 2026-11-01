"""A stage that computed nothing must not exit zero: pytest tests/test_stage_products.py

Job 31186038 asked for 3 ps of QM/MM MD, built the system, printed [DONE] in
14 s and was recorded COMPLETED, with a results record holding only the QM
level and the metal parameters - no energy, no trajectory, no structure (D153).
Slurm's own record said success, so nothing downstream had any reason to doubt
it, and the rate it was submitted to measure stayed unmeasured for a day.

`confirm_produced` is the one place that can catch the whole class, because
every stage's result passes through it. It lives in cli.py, which imports only
argparse/importlib/os, so these run without ASH or OpenMM.
"""

import pytest

from mcube.cli import STAGE_PRODUCTS, confirm_produced


class TestItRaisesOnAResultThatComputedNothing:
    def test_the_qmrate_record_verbatim(self):
        """What 31186038 actually wrote into config.yaml."""
        with pytest.raises(SystemExit) as exc:
            confirm_produced("qmmm", {"qm_level": "GFN2-xTB", "metal_parameters": "iod"})
        assert "computed nothing" in str(exc.value)

    def test_the_message_names_what_was_missing_and_what_was_there(self):
        with pytest.raises(SystemExit) as exc:
            confirm_produced("qmmm", {"qm_level": "GFN2-xTB"})
        message = str(exc.value)
        assert "optimized_pdb" in message and "qm_level" in message

    def test_an_empty_record_says_so(self):
        with pytest.raises(SystemExit) as exc:
            confirm_produced("sp", {})
        assert "an empty record" in str(exc.value)

    def test_a_product_key_present_but_none_does_not_count(self):
        """A key written as a placeholder is not a result."""
        with pytest.raises(SystemExit):
            confirm_produced("sp", {"sp_energy": None})


class TestItPassesRealResults:
    def test_an_energy_is_enough(self):
        confirm_produced("sp", {"sp_energy": -5076.869006116424})

    def test_any_one_of_the_keys_is_enough(self):
        """An MD-only qmmm run has no optimised structure, and is still a run."""
        confirm_produced("qmmm", {"trajfile": "/scratch/QM_MM.dcd"})

    def test_a_zero_energy_is_not_treated_as_missing(self):
        """0.0 is falsy and is a perfectly good relative energy."""
        confirm_produced("ts", {"ts_energy": 0.0})

    def test_the_asph_record_passes(self):
        confirm_produced("qmmm", {"qm_level": "GFN2-xTB", "converged": True,
                                  "opt_energy": -447.44914996398757,
                                  "optimized_pdb": "/scratch/QM_MM_optimized.pdb"})


class TestWhatIsAndIsNotChecked:
    def test_an_unlisted_stage_is_not_checked(self):
        """umbrella's three entry points share no result key; see STAGE_PRODUCTS."""
        confirm_produced("umbrella", {"anything": 1})
        confirm_produced("umbrella", {})

    def test_umbrella_stays_unlisted(self):
        """If someone lists it, they have to decide which of its three shapes to
        demand - this test is here to make that a deliberate change."""
        assert "umbrella" not in STAGE_PRODUCTS

    def test_a_non_dict_result_is_not_checked(self):
        confirm_produced("sp", None)
        confirm_produced("sp", "done")

    def test_every_listed_stage_is_a_real_stage_module(self):
        import importlib.util
        for stage in STAGE_PRODUCTS:
            assert importlib.util.find_spec(f"mcube.{stage}") is not None, stage
