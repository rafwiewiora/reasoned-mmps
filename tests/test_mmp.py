from __future__ import annotations

import csv
import json
import unittest
from pathlib import Path

from reasoned_mmp.mmp import best_candidate_per_parent, infer_parent_candidates


ROOT = Path(__file__).resolve().parents[1]


def load_fixture():
    with (ROOT / "data/photo_clenbuterol_compounds.csv").open(newline="") as handle:
        compounds = list(csv.DictReader(handle))
    reasons = json.loads(
        (ROOT / "data/photo_clenbuterol_reasons.json").read_text()
    )
    return compounds, reasons


class ParentInferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.compounds, cls.reasons = load_fixture()
        cls.by_id = {row["chembl_id"]: row for row in cls.compounds}

    def test_reason_aligned_cl_to_cn_is_top_candidate(self):
        reason = self.reasons[0]
        child = self.by_id[reason["child_chembl_id"]]
        candidates = best_candidate_per_parent(
            infer_parent_candidates(child, self.compounds, reason)
        )
        top = candidates[0]
        self.assertEqual(top["parent_chembl_id"], "CHEMBL6190593")
        self.assertEqual(top["transformation"], "Cl[*:1]>>N#C[*:1]")
        self.assertEqual(top["scores"]["reason_transform_alignment"], 1.0)
        self.assertEqual(top["edge_semantics"], "reason_constrained_mmp_comparator")
        self.assertFalse(top["historical_parent_claim"])

    def test_hydrogen_change_indexes_secondary_neighbor(self):
        reason = self.reasons[0]
        child = self.by_id[reason["child_chembl_id"]]
        candidates = best_candidate_per_parent(
            infer_parent_candidates(child, self.compounds, reason)
        )
        h_candidate = next(
            row for row in candidates if row["parent_chembl_id"] == "CHEMBL6188786"
        )
        self.assertTrue(h_candidate["hydrogen_change"])
        self.assertEqual(h_candidate["transformation"], "[*:1][H]>>N#C[*:1]")
        self.assertEqual(
            h_candidate["edge_semantics"], "structure_only_mmp_comparator"
        )

    def test_large_explicit_azoextension_is_sensitivity_only(self):
        reason = self.reasons[1]
        child = self.by_id[reason["child_chembl_id"]]
        strict = infer_parent_candidates(child, self.compounds, reason)
        sensitive = infer_parent_candidates(
            child, self.compounds, reason, max_variable_fraction=0.40
        )
        self.assertNotIn(
            "CHEMBL6189525", {row["parent_chembl_id"] for row in strict}
        )
        match = next(
            row for row in sensitive if row["parent_chembl_id"] == "CHEMBL6189525"
        )
        self.assertEqual(match["scores"]["reason_transform_alignment"], 1.0)


if __name__ == "__main__":
    unittest.main()
