"""Build the frozen multi-paper reasoned-MMP corpus and evidence viewer."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem.Draw import rdMolDraw2D

from .mmp import best_candidate_per_parent, infer_parent_candidates, stable_id
from .outcomes import compare_measurements


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUTPUTS = ROOT / "outputs"
DOCS = ROOT / "docs"
DATASET_PREFIXES = tuple(
    sorted(
        path.name.removesuffix("_compounds.csv")
        for path in DATA.glob("*_compounds.csv")
    )
)


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


def _write_depiction(compound: dict, path: Path) -> None:
    mol = Chem.MolFromSmiles(compound["canonical_smiles"])
    if mol is None:
        return
    drawer = rdMolDraw2D.MolDraw2DSVG(420, 260)
    options = drawer.drawOptions()
    options.clearBackground = False
    options.padding = 0.08
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    path.write_text(drawer.GetDrawingText())


def _write_viewer_data(
    manifest: dict, moves: list[dict], compounds: list[dict]
) -> None:
    molecule_dir = DOCS / "molecules"
    molecule_dir.mkdir(parents=True, exist_ok=True)
    for compound in compounds:
        _write_depiction(compound, molecule_dir / f"{compound['chembl_id']}.svg")
    payload = {
        "manifest": manifest,
        "moves": moves,
        "compounds": compounds,
    }
    (DOCS / "data.js").write_text(
        "window.REASONED_MMP_DATA = "
        + json.dumps(payload, indent=2, sort_keys=True)
        + ";\n"
    )


def _explicit_parent_edge(
    reason: dict, child: dict, compounds: list[dict]
) -> dict | None:
    explicit_id = reason.get("explicit_parent_chembl_id")
    if not explicit_id:
        return None
    # The primary MMP rule is variable fraction <= 0.30. A 0.50 sensitivity
    # pass is used only for an independently author-stated relationship. This
    # captures small fragment-model pairs and large one-site azoextensions
    # without broadening ordinary inferred-candidate generation.
    sensitivity_limit = 0.50
    sensitivity = infer_parent_candidates(
        child,
        compounds,
        reason,
        max_variable_fraction=sensitivity_limit,
    )
    matches = [row for row in sensitivity if row["parent_chembl_id"] == explicit_id]
    witness = matches[0] if matches else None
    return {
        "relationship_id": stable_id(
            "rel", [reason["reason_id"], explicit_id, child["chembl_id"]]
        ),
        "parent_chembl_id": explicit_id,
        "child_chembl_id": child["chembl_id"],
        "edge_semantics": reason.get(
            "explicit_parent_edge_semantics", "author_explicit_scaffold_parent"
        ),
        "author_relationship_explicit": True,
        "paper_uses_parent_term": bool(reason.get("paper_uses_parent_term", False)),
        "historical_synthesis_lineage_claim": reason.get(
            "explicit_historical_synthesis_lineage", False
        ),
        "valid_primary_mmp_rule": any(
            row["parent_chembl_id"] == explicit_id
            for row in infer_parent_candidates(child, compounds, reason)
        ),
        "valid_single_cut_sensitivity_rule": witness is not None,
        "sensitivity_rule_max_variable_fraction": sensitivity_limit,
        "structural_witness": witness,
        "basis": reason.get("explicit_parent_basis"),
    }


def build() -> dict:
    OUTPUTS.mkdir(exist_ok=True)
    compounds: list[dict] = []
    reasons: list[dict] = []
    measurements: list[dict] = []
    outcome_specs: list[dict] = []
    input_paths: list[Path] = []
    compounds_by_dataset: dict[str, list[dict]] = defaultdict(list)
    for dataset_id in DATASET_PREFIXES:
        compound_path = DATA / f"{dataset_id}_compounds.csv"
        reason_path = DATA / f"{dataset_id}_reasons.json"
        measurement_path = DATA / f"{dataset_id}_measurements.csv"
        spec_path = DATA / f"{dataset_id}_outcome_specs.json"
        input_paths.extend((compound_path, reason_path, measurement_path, spec_path))
        dataset_compounds = _read_csv(compound_path)
        dataset_reasons = _read_json(reason_path)
        dataset_measurements = _read_csv(measurement_path)
        dataset_specs = _read_json(spec_path)
        for row in dataset_compounds:
            row["dataset_id"] = dataset_id
        for row in dataset_reasons:
            row["dataset_id"] = dataset_id
        for row in dataset_measurements:
            row["dataset_id"] = dataset_id
        compounds.extend(dataset_compounds)
        reasons.extend(dataset_reasons)
        measurements.extend(dataset_measurements)
        outcome_specs.extend(dataset_specs)
        compounds_by_dataset[dataset_id].extend(dataset_compounds)
    compounds_by_key = {
        (row["dataset_id"], row["chembl_id"]): row for row in compounds
    }
    measurements_by_id = {row["measurement_id"]: row for row in measurements}
    reason_ids = [reason["reason_id"] for reason in reasons]
    measurement_ids = [measurement["measurement_id"] for measurement in measurements]
    if len(reason_ids) != len(set(reason_ids)):
        raise ValueError("reason_id values must be unique")
    if len(measurement_ids) != len(set(measurement_ids)):
        raise ValueError("measurement_id values must be unique")
    for reason in reasons:
        dataset_ids = {
            compound["chembl_id"] for compound in compounds_by_dataset[reason["dataset_id"]]
        }
        members = reason.get("member_chembl_ids", [reason["child_chembl_id"]])
        if reason["child_chembl_id"] not in members:
            raise ValueError(f"{reason['reason_id']}: anchor child is not a named member")
        missing = sorted(set(members) - dataset_ids)
        if missing:
            raise ValueError(f"{reason['reason_id']}: unresolved named members {missing}")
    known_reasons = set(reason_ids)
    for spec in outcome_specs:
        if spec["reason_id"] not in known_reasons:
            raise ValueError(f"unknown reason_id in outcome spec: {spec['reason_id']}")
        for field in ("parent_measurement_id", "child_measurement_id"):
            if spec[field] not in measurements_by_id:
                raise ValueError(f"unknown {field} in outcome spec: {spec[field]}")

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
            comparison_scale=spec.get("comparison_scale", "linear"),
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
        child = compounds_by_key[(reason["dataset_id"], reason["child_chembl_id"])]
        candidate_compounds = compounds_by_dataset[reason["dataset_id"]]
        decompositions = infer_parent_candidates(child, candidate_compounds, reason)
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
            "entities": {"child": child},
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
                    reason, child, candidate_compounds
                ),
                "additional_author_relationships": reason.get(
                    "related_author_relationships", []
                ),
            },
            "layer_4_observed_outcomes": {
                "comparisons": reason_outcomes,
                "unpaired_facts": reason.get("unpaired_outcome_facts", []),
                "classification_counts": dict(counts),
                "direct_supporting_endpoint_counts": dict(direct_counts),
                "stated_intent_outcome": (
                    "indeterminate_direct_endpoint_unavailable"
                    if reason["outcome_join_status"]
                    == "direct_stated_property_not_measured_in_pilot"
                    else "not_applicable_retrospective_explanation"
                    if reason["outcome_join_status"]
                    == "retrospective_explanation_no_prospective_success_label"
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

    reason_bearing_compounds = {
        chembl_id
        for reason in reasons
        for chembl_id in reason.get("member_chembl_ids", [reason["child_chembl_id"]])
    }
    unique_outcome_pairs = {
        (row["parent_measurement_id"], row["child_measurement_id"])
        for row in outcome_specs
    }
    coverage_by_paper = []
    for document_id in sorted({reason["evidence"]["document_id"] for reason in reasons}):
        paper_reasons = [reason for reason in reasons if reason["evidence"]["document_id"] == document_id]
        paper_members = {
            chembl_id
            for reason in paper_reasons
            for chembl_id in reason.get("member_chembl_ids", [reason["child_chembl_id"]])
        }
        paper_structures = [compound for compound in compounds if compound["source_document_id"] == document_id]
        coverage_by_paper.append({
            "document_id": document_id,
            "source_title": paper_reasons[0].get("source_title", document_id),
            "rationale_episodes": len(paper_reasons),
            "reason_bearing_compounds": len(paper_members),
            "resolved_structures": len(paper_structures),
        })

    manifest = {
        "corpus": "reasoned_mmp_corpus",
        "schema_version": "0.4.0",
        "rdkit_version": rdBase.rdkitVersion,
        "input_sha256": {
            path.name: _sha256(path) for path in input_paths
        },
        "counts": {
            "papers": len(
                {reason["evidence"]["document_id"] for reason in reasons}
            ),
            "rationale_episodes": len(reasons),
            "reason_assertions": len(reasons),
            "reason_bearing_compounds": len(reason_bearing_compounds),
            "resolved_structures": len(compounds),
            "compounds": len(compounds),
            "mmp_decompositions": len(all_decompositions),
            "unique_reason_parent_candidates": len(best_candidates),
            "outcome_comparisons": len(outcome_rows),
            "unique_outcome_pairs": len(unique_outcome_pairs),
            "unpaired_outcome_facts": sum(len(reason.get("unpaired_outcome_facts", [])) for reason in reasons),
        },
        "paper_coverage": coverage_by_paper,
        "guardrails": {
            "reason_frozen_before_outcome_join": True,
            "bounds_preserved": True,
            "inferred_comparator_is_not_lineage": True,
            "primary_max_variable_fraction": 0.30,
            "explicit_relationship_sensitivity_fraction": 0.50,
            "hydrogen_changes_indexed": True,
        },
    }
    (OUTPUTS / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    _write_viewer_data(manifest, reasoned_moves, compounds)
    return manifest
