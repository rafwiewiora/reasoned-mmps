from __future__ import annotations

import json
import unittest
from pathlib import Path

from reasoned_mmp.pipeline import ROOT, build


class PipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = build()
        cls.moves = json.loads((ROOT / "outputs/reasoned_moves.json").read_text())

    def test_manifest_records_core_guardrails(self):
        guardrails = self.manifest["guardrails"]
        self.assertTrue(guardrails["reason_frozen_before_outcome_join"])
        self.assertTrue(guardrails["inferred_comparator_is_not_lineage"])
        self.assertTrue(guardrails["bounds_preserved"])

    def test_direct_reason_outcome_is_not_backfilled_from_proxy(self):
        chlorocyano = next(
            move
            for move in self.moves
            if move["layer_2_extracted_design_intent"]["reason_id"].endswith(
                "chlorocyano_sar"
            )
        )
        outcomes = chlorocyano["layer_4_observed_outcomes"]
        self.assertEqual(
            outcomes["stated_intent_outcome"],
            "indeterminate_direct_endpoint_unavailable",
        )
        self.assertEqual(outcomes["classification_counts"]["worsened"], 3)

    def test_explicit_parent_is_not_recast_as_synthesis_lineage(self):
        azo = next(
            move
            for move in self.moves
            if move["layer_2_extracted_design_intent"]["reason_id"].endswith(
                "azoextension"
            )
        )
        relation = azo["layer_3_inferred_structural_comparison"][
            "author_explicit_relationship"
        ]
        self.assertEqual(relation["parent_chembl_id"], "CHEMBL6189525")
        self.assertFalse(relation["historical_synthesis_lineage_claim"])
        self.assertFalse(relation["valid_primary_mmp_rule"])
        self.assertTrue(relation["valid_single_cut_sensitivity_rule"])


if __name__ == "__main__":
    unittest.main()
