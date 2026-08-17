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
        self.assertEqual(self.manifest["counts"]["papers"], 4)
        self.assertEqual(self.manifest["counts"]["reason_assertions"], 5)
        self.assertEqual(self.manifest["counts"]["outcome_comparisons"], 19)
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

    def test_new_reason_aligned_parents_and_transforms(self):
        by_reason = {
            move["layer_2_extracted_design_intent"]["reason_id"]: move
            for move in self.moves
        }
        expected = {
            "PMC5807869:2:late_stage_oxidation": (
                "CHEMBL4160171",
                "[*:1][H]>>O[*:1]",
                "reason_constrained_mmp_comparator",
            ),
            "PMC4207553:5:dimethylisoxazole_solubility": (
                "CHEMBL1873309",
                "c1ccc([*:1])cc1>>Cc1noc(C)c1[*:1]",
                "reason_constrained_mmp_comparator",
            ),
            "PMC10726475:9:noralkoxy_basicity": (
                "CHEMBL288441",
                "CN1CCN(C[*:1])CC1>>CN1CCN(O[*:1])CC1",
                "retrospective_reason_aligned_mmp_comparator",
            ),
        }
        for reason_id, (parent_id, transform, semantics) in expected.items():
            top = by_reason[reason_id]["layer_3_inferred_structural_comparison"][
                "top_candidates"
            ][0]
            self.assertEqual(top["parent_chembl_id"], parent_id)
            self.assertEqual(top["transformation"], transform)
            self.assertEqual(top["edge_semantics"], semantics)

    def test_retrospective_reason_has_no_prospective_success_label(self):
        move = next(
            move
            for move in self.moves
            if move["layer_2_extracted_design_intent"]["assertion_class"]
            == "retrospective_explanation"
        )
        self.assertEqual(
            move["layer_4_observed_outcomes"]["stated_intent_outcome"],
            "not_applicable_retrospective_explanation",
        )

    def test_static_viewer_assets_are_generated(self):
        data_js = ROOT / "docs/data.js"
        self.assertTrue(data_js.read_text().startswith("window.REASONED_MMP_DATA = "))
        for chembl_id in ("CHEMBL4062397", "CHEMBL3343650", "CHEMBL5435819"):
            svg = ROOT / f"docs/molecules/{chembl_id}.svg"
            self.assertTrue(svg.exists())
            self.assertIn("<svg", svg.read_text())


if __name__ == "__main__":
    unittest.main()
