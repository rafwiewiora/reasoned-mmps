from __future__ import annotations

import unittest

from reasoned_mmp.outcomes import compare_measurements


def measurement(relation, value, *, state="PSS-cis"):
    return {
        "endpoint": "pKi",
        "state": state,
        "units": "pLog",
        "assay_context": "human_beta2_HEK293T_competition_binding",
        "relation": relation,
        "value": value,
    }


class OutcomeTests(unittest.TestCase):
    def test_censored_child_is_definitively_worse(self):
        result = compare_measurements(
            measurement("=", 6.4),
            measurement("<", 5.0),
            higher_is_better=True,
            equivalence_margin=0.3,
        )
        self.assertEqual(result["classification"], "worsened")
        self.assertIsNone(result["delta_lower"])
        self.assertEqual(result["delta_upper"], -1.4)
        self.assertTrue(result["censoring_preserved"])

    def test_two_equal_bounds_do_not_become_a_numeric_delta(self):
        result = compare_measurements(
            measurement("<", 5.0),
            measurement("<", 5.0),
            higher_is_better=True,
            equivalence_margin=0.3,
        )
        self.assertEqual(result["classification"], "indeterminate")
        self.assertIsNone(result["delta_lower"])
        self.assertIsNone(result["delta_upper"])

    def test_different_states_are_not_joined(self):
        result = compare_measurements(
            measurement("=", 6.4, state="PSS-cis"),
            measurement("=", 6.4, state="trans"),
            higher_is_better=True,
        )
        self.assertEqual(result["assay_comparability"], "not_comparable")
        self.assertIn("state", result["mismatched_fields"])

    def test_log_scale_compares_fold_changes_and_preserves_bounds(self):
        parent = measurement("=", 10)
        child = measurement("<", 2)
        result = compare_measurements(
            parent,
            child,
            higher_is_better=False,
            equivalence_margin=0.3,
            comparison_scale="log10",
        )
        self.assertEqual(result["classification"], "improved")
        self.assertEqual(result["comparison_scale"], "log10")
        self.assertGreater(result["delta_lower"], 0.69)
        self.assertIsNone(result["delta_upper"])


if __name__ == "__main__":
    unittest.main()
