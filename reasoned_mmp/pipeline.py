"""Build the frozen photo-clenbuterol reasoned-MMP pilot artifacts."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import rdBase

from .mmp import best_candidate_per_parent, infer_parent_candidates, stable_id
from .outcomes import compare_measurements


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUTPUTS = ROOT / "outputs"


def _read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _read_json(path: Path):
    with path.open() as handle:
        return json.load(handle)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _flatten_candidate(row: dict) -> dict:
    flat = {key: value for key, value in row.items() if key not in {"scores", "alignment_checks"}}
    flat.update(row["scores"])
    flat["alignment_checks"] = json.dumps(row["alignment_checks"], sort_keys=True)
    return flat


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _explicit_parent_edge(
    reason: dict, child: dict, compounds: list[dict]
) -> dict | None:
    explicit_id = reason.get("explicit_parent_chembl_id")
    if not explicit_id:
        return None
    # The primary MMP rule is variable fraction <= 0.30. A sensitivity pass at
    # 0.40 checks large but still one-site series moves such as azoextension.
    sensitivity = infer_parent_candidates(
        child,
        compounds,
        reason,
        max_variable_fraction=0.40,
    )
    matches = [row for row in sensitivity if row["parent_chembl_id"] == explicit_id]
    witness = matches[0] if matches else None
    return {
        "relationship_id": stable_id(
            "rel", [reason["reason_id"], explicit_id, child["chembl_id"]]
        ),
        "parent_chembl_id": explicit_id,
        "child_chembl_id": child["chembl_id"],
        "edge_semantics": "author_explicit_scaffold_parent",
        "paper_uses_parent_term": True,
        "historical_synthesis_lineage_claim": False,
        "valid_primary_mmp_rule": any(
            row["parent_chembl_id"] == explicit_id
            for row in infer_parent_candidates(child, compounds, reason)
        ),
        "valid_single_cut_sensitivity_rule": witness is not None,
        "sensitivity_rule_max_variable_fraction": 0.40,
        "structural_witness": witness,
        "basis": reason.get("explicit_parent_basis"),
    }


def build() -> dict:
    OUTPUTS.mkdir(exist_ok=True)
    compound_path = DATA / "photo_clenbuterol_compounds.csv"
    reason_path = DATA / "photo_clenbuterol_reasons.json"
    measurement_path = DATA / "photo_clenbuterol_measurements.csv"
    spec_path = DATA / "photo_clenbuterol_outcome_specs.json"

    compounds = _read_csv(compound_path)
    reasons = _read_json(reason_path)
    measurements = _read_csv(measurement_path)
    outcome_specs = _read_json(spec_path)
    compounds_by_id = {row["chembl_id"]: row for row in compounds}
    measurements_by_id = {row["measurement_id"]: row for row in measurements}

    outcome_rows: list[dict] = []
    outcomes_by_reason: dict[str, list[dict]] = defaultdict(list)
    for spec in outcome_specs:
        parent = measurements_by_id[spec["parent_measurement_id"]]
        child = measurements_by_id[spec["child_measurement_id"]]
        result = compare_measurements(
            parent,
            child,
            higher_is_better=bool(spec["higher_is_better"]),
            equivalence_margin=float(spec["equivalence_margin"]),
        )
        row = {
            **spec,
            "parent_chembl_id": parent["chembl_id"],
            "child_chembl_id": child["chembl_id"],
            "endpoint": child["endpoint"],
            "state": child["state"],
            "parent_relation": parent["relation"],
            "parent_value": parent["value"],
            "child_relation": child["relation"],
            "child_value": child["value"],
            "units": child["units"],
            **result,
        }
        outcome_rows.append(row)
        outcomes_by_reason[spec["reason_id"]].append(row)

    all_decompositions: list[dict] = []
    best_candidates: list[dict] = []
    reasoned_moves: list[dict] = []
    for reason in reasons:
        child = compounds_by_id[reason["child_chembl_id"]]
        decompositions = infer_parent_candidates(child, compounds, reason)
        best = best_candidate_per_parent(decompositions)
        all_decompositions.extend(decompositions)
        best_candidates.extend(best)
        top = best[:3]
        top_margin = None
        if len(top) >= 2:
            top_margin = round(
                top[0]["scores"]["ranking_score_uncalibrated"]
                - top[1]["scores"]["ranking_score_uncalibrated"],
                6,
            )
        reason_outcomes = outcomes_by_reason.get(reason["reason_id"], [])
        counts = Counter(row["classification"] for row in reason_outcomes)
        direct = [
            row
            for row in reason_outcomes
            if row["relation_to_stated_intent"] == "direct_supporting_endpoint"
        ]
        direct_counts = Counter(row["classification"] for row in direct)
        move = {
            "reasoned_move_id": stable_id("move", reason["reason_id"]),
            "layer_1_author_evidence": reason["evidence"],
            "layer_2_extracted_design_intent": reason,
            "layer_3_inferred_structural_comparison": {
                "candidate_universe": "same_paper",
                "candidate_search_order_planned": [
                    "same_scheme_or_table",
                    "same_named_series",
                    "same_paper",
                    "cited_predecessor",
                    "same_chembl_document",
                    "global_chembl",
                ],
                "implemented_scope": "same_paper",
                "top_candidates": top,
                "top1_top2_margin": top_margin,
                "score_status": "uncalibrated_heuristic",
                "historical_lineage_inferred": False,
                "author_explicit_relationship": _explicit_parent_edge(
                    reason, child, compounds
                ),
            },
            "layer_4_observed_outcomes": {
                "comparisons": reason_outcomes,
                "classification_counts": dict(counts),
                "direct_supporting_endpoint_counts": dict(direct_counts),
                "stated_intent_outcome": (
                    "indeterminate_direct_endpoint_unavailable"
                    if reason["outcome_join_status"]
                    == "direct_stated_property_not_measured_in_pilot"
                    else "evaluate_per_endpoint"
                ),
            },
            "confidence_components": {
                "evidence_entailment": reason["extraction_confidence"],
                "child_resolution": 1.0,
                "parent_inference": (
                    top[0]["scores"]["ranking_score_uncalibrated"] if top else None
                ),
                "assay_comparability": (
                    sum(
                        row["assay_comparability"] == "exact_context_match"
                        for row in reason_outcomes
                    )
                    / len(reason_outcomes)
                    if reason_outcomes
                    else None
                ),
                "outcome_confidence": "not_calibrated",
            },
        }
        reasoned_moves.append(move)

    (OUTPUTS / "reasoned_moves.json").write_text(
        json.dumps(reasoned_moves, indent=2, sort_keys=True) + "\n"
    )
    with (OUTPUTS / "decomposition_witnesses.jsonl").open("w") as handle:
        for row in all_decompositions:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    _write_csv(
        OUTPUTS / "candidate_edges.csv",
        [_flatten_candidate(row) for row in best_candidates],
    )
    _write_csv(OUTPUTS / "outcome_comparisons.csv", outcome_rows)

    manifest = {
        "pilot": "photo_clenbuterol",
        "schema_version": "0.1.0",
        "rdkit_version": rdBase.rdkitVersion,
        "input_sha256": {
            path.name: _sha256(path)
            for path in (compound_path, reason_path, measurement_path, spec_path)
        },
        "counts": {
            "compounds": len(compounds),
            "reason_assertions": len(reasons),
            "mmp_decompositions": len(all_decompositions),
            "unique_reason_parent_candidates": len(best_candidates),
            "outcome_comparisons": len(outcome_rows),
        },
        "guardrails": {
            "reason_frozen_before_outcome_join": True,
            "bounds_preserved": True,
            "inferred_comparator_is_not_lineage": True,
            "primary_max_variable_fraction": 0.30,
            "hydrogen_changes_indexed": True,
        },
    }
    (OUTPUTS / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return manifest
