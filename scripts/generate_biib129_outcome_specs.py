"""Generate the reviewed BIIB129 outcome-pair specification.

The reason records are frozen before this script joins exact-context
measurements.  Keeping the compact declarations here makes duplicated assay
pairs across distinct evidence episodes visible and reproducible.
"""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

ENDPOINTS = {
    "btk_ic50": (False, 0.3, "log10"),
    "log_kinact_ki": (True, 0.3, "linear"),
    "wb_cd69": (False, 0.3, "log10"),
    "tmd8_ic50": (False, 0.3, "log10"),
    "mdr1_ratio": (False, 0.3, "log10"),
    "rlm_cl": (False, 0.3, "log10"),
    "hlm_cl": (False, 0.3, "log10"),
    "kpuu": (True, 0.3, "log10"),
    "gsh_t_half": (True, 0.3, "log10"),
}


PAIRS = [
    {
        "tag": "linker",
        "reason_id": "PMC11129193:series:sp3_linker_optimization",
        "parent": "7",
        "child": "10",
        "endpoints": ["btk_ic50", "log_kinact_ki", "mdr1_ratio"],
        "relations": {
            "btk_ic50": "direct_supporting_endpoint",
            "log_kinact_ki": "direct_supporting_endpoint",
            "mdr1_ratio": "CNS_property_retention",
        },
    },
    {
        "tag": "aminoether",
        "reason_id": "PMC11129193:series:aminoether_metabolic_optimization",
        "parent": "10",
        "child": "13",
        "endpoints": ["log_kinact_ki", "wb_cd69", "rlm_cl", "hlm_cl", "mdr1_ratio", "kpuu", "gsh_t_half"],
        "relations": {
            "rlm_cl": "direct_supporting_endpoint",
            "hlm_cl": "direct_supporting_endpoint",
            "log_kinact_ki": "potency_retention",
            "wb_cd69": "potency_retention",
            "mdr1_ratio": "CNS_property_context",
            "kpuu": "CNS_property_context",
            "gsh_t_half": "warhead_stability_context",
        },
    },
    {
        "tag": "alpha_fluoro",
        "reason_id": "PMC11129193:11:alpha_fluoro_clearance",
        "parent": "10",
        "child": "11",
        "endpoints": ["log_kinact_ki", "wb_cd69", "mdr1_ratio", "rlm_cl", "hlm_cl"],
        "relations": {"rlm_cl": "direct_supporting_endpoint", "hlm_cl": "direct_supporting_endpoint"},
    },
    {
        "tag": "fluoro_relocation",
        "reason_id": "PMC11129193:12:fluoro_inductive_explanation",
        "parent": "11",
        "child": "12",
        "endpoints": ["log_kinact_ki", "wb_cd69", "mdr1_ratio", "rlm_cl", "hlm_cl"],
        "default_relation": "effect_explained_retrospectively",
    },
    {
        "tag": "bicyclic",
        "reason_id": "PMC11129193:series:bicyclic_steric_metabolism",
        "parent": "10",
        "child": "13",
        "endpoints": ["log_kinact_ki", "wb_cd69", "rlm_cl", "hlm_cl", "mdr1_ratio", "kpuu", "gsh_t_half"],
        "relations": {
            "rlm_cl": "direct_supporting_endpoint",
            "hlm_cl": "direct_supporting_endpoint",
            "log_kinact_ki": "potency_retention",
            "wb_cd69": "potency_retention",
            "mdr1_ratio": "CNS_property_context",
            "kpuu": "CNS_property_context",
            "gsh_t_half": "warhead_stability_context",
        },
    },
    {
        "tag": "docking",
        "reason_id": "PMC11129193:series:docking_linker_replacements",
        "parent": "10",
        "child": "17",
        "endpoints": ["log_kinact_ki", "wb_cd69", "mdr1_ratio", "rlm_cl", "hlm_cl"],
        "relations": {
            "log_kinact_ki": "covalent_geometry_context",
            "wb_cd69": "potency_retention",
            "mdr1_ratio": "CNS_property_context",
            "rlm_cl": "direct_supporting_endpoint",
            "hlm_cl": "direct_supporting_endpoint",
        },
    },
    {
        "tag": "four_bond",
        "reason_id": "PMC11129193:20_21:four_bond_geometry_probe",
        "parent": "20",
        "child": "21",
        "endpoints": ["log_kinact_ki", "wb_cd69", "tmd8_ic50", "mdr1_ratio", "rlm_cl", "hlm_cl"],
        "relations": {
            "log_kinact_ki": "direct_supporting_endpoint",
            "wb_cd69": "cellular_translation_context",
            "tmd8_ic50": "cellular_translation_context",
            "mdr1_ratio": "CNS_property_context",
            "rlm_cl": "ADME_context",
            "hlm_cl": "ADME_context",
        },
    },
    {
        "tag": "bridge",
        "reason_id": "PMC11129193:22_23:constrained_cyclobutyl_sar",
        "parent": "22",
        "child": "23",
        "endpoints": ["log_kinact_ki", "wb_cd69", "tmd8_ic50", "mdr1_ratio", "rlm_cl", "hlm_cl", "gsh_t_half"],
        "relations": {
            "log_kinact_ki": "direct_supporting_endpoint",
            "wb_cd69": "direct_supporting_endpoint",
            "tmd8_ic50": "direct_supporting_endpoint",
            "mdr1_ratio": "CNS_property_context",
            "rlm_cl": "ADME_context",
            "hlm_cl": "ADME_context",
            "gsh_t_half": "warhead_stability_context",
        },
    },
    {
        "tag": "deconstruct",
        "reason_id": "PMC11129193:24_25:deconstruct_cellular_disconnect",
        "parent": "23",
        "child": "25",
        "endpoints": ["log_kinact_ki", "wb_cd69", "tmd8_ic50", "mdr1_ratio", "rlm_cl", "hlm_cl", "kpuu", "gsh_t_half"],
        "relations": {
            "wb_cd69": "direct_supporting_endpoint",
            "tmd8_ic50": "direct_supporting_endpoint",
            "log_kinact_ki": "biochemical_potency_context",
            "mdr1_ratio": "CNS_property_context",
            "rlm_cl": "ADME_context",
            "hlm_cl": "ADME_context",
            "kpuu": "CNS_property_context",
            "gsh_t_half": "warhead_stability_context",
        },
    },
    {
        "tag": "methyl_14_26",
        "reason_id": "PMC11129193:26_27:alpha_methyl_transfer",
        "parent": "14",
        "child": "26",
        "endpoints": ["log_kinact_ki", "wb_cd69", "tmd8_ic50", "mdr1_ratio", "rlm_cl", "hlm_cl"],
        "relations": {"wb_cd69": "direct_supporting_endpoint", "tmd8_ic50": "direct_supporting_endpoint"},
    },
    {
        "tag": "methyl_15_27",
        "reason_id": "PMC11129193:26_27:alpha_methyl_transfer",
        "parent": "15",
        "child": "27",
        "endpoints": ["log_kinact_ki", "wb_cd69", "tmd8_ic50", "mdr1_ratio", "rlm_cl", "hlm_cl"],
        "relations": {"wb_cd69": "direct_supporting_endpoint", "tmd8_ic50": "direct_supporting_endpoint"},
    },
    {
        "tag": "conformation",
        "reason_id": "PMC11129193:27:conformational_restriction_explanation",
        "parent": "25",
        "child": "27",
        "endpoints": ["log_kinact_ki", "wb_cd69", "tmd8_ic50", "mdr1_ratio", "rlm_cl", "hlm_cl", "kpuu", "gsh_t_half"],
        "default_relation": "effect_explained_retrospectively",
    },
]


def main() -> None:
    rows = []
    for pair in PAIRS:
        for endpoint in pair["endpoints"]:
            higher_is_better, margin, scale = ENDPOINTS[endpoint]
            relation = pair.get("relations", {}).get(
                endpoint, pair.get("default_relation", "supporting_endpoint")
            )
            rows.append(
                {
                    "comparison_id": f"biib_{pair['parent']}_{pair['child']}_{endpoint}_{pair['tag']}",
                    "reason_id": pair["reason_id"],
                    "parent_measurement_id": f"m_biib_{pair['parent']}_{endpoint}",
                    "child_measurement_id": f"m_biib_{pair['child']}_{endpoint}",
                    "higher_is_better": higher_is_better,
                    "equivalence_margin": margin,
                    "comparison_scale": scale,
                    "relation_to_stated_intent": relation,
                }
            )
    path = ROOT / "data" / "biib129_outcome_specs.json"
    path.write_text(json.dumps(rows, indent=2) + "\n")
    print(f"wrote {len(rows)} outcome links to {path}")


if __name__ == "__main__":
    main()
