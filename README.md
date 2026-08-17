# Reasoned MMPs

[![test](https://github.com/rafwiewiora/reasoned-mmps/actions/workflows/test.yml/badge.svg)](https://github.com/rafwiewiora/reasoned-mmps/actions/workflows/test.yml)

Turn medicinal-chemistry prose into auditable **moves with reasons**:

```text
author evidence → extracted design intent → inferred MMP comparator → observed outcome
```

ChEMBL already captures much of *what was measured*. This project captures
*why a chemist chose a feature*, infers the most defensible matched comparator
afterward, and then asks what happened in genuinely comparable assays.

The central rule is simple: **an inferred parent is an analytic comparator, not
a claim about synthetic or historical lineage**.

## First pilot: photo-clenbuterol

The frozen pilot uses the 2025 paper
[Photo-clenbuterol](https://doi.org/10.1021/acs.jmedchem.5c00792). The authors
selected chloro-cyano compound `12e` using published clenbuterol SAR that had
favored a chloro-cyano pattern over a dichloro pattern in guinea-pig
bronchodilation.

Given only the reason-bearing child and same-paper structures, the ranker finds:

| Rank | Comparator → child | Move | Interpretation |
|---:|---|---|---|
| 1 | `12b → 12e` | `aryl-Cl → aryl-CN` | Fully reason-aligned MMP comparator |
| 2 | `12c → 12e` | `aryl-H → aryl-CN` | Valid structural neighbor, but only partly aligned |

The first move is chemically clean and rationale-aligned, but the observed
transfer is mostly negative: PSS-cis binding falls from `pKi 6.4 ± 0.2` to
`pKi < 5.0`, while PSS-cis antagonism is approximately preserved (`pKb 6.9`
to `6.8`). The paper's stated precedent endpoint was bronchodilation, however,
so the strict outcome for that exact intent remains **indeterminate**. The β2
measurements are explicitly labeled as proxy or adjacent pharmacology.

The same child also has a different relationship: compound `18` is the
author-described scaffold parent for an azoextension. That large one-site move
passes a `0.40` variable-fraction sensitivity rule but not the primary `0.30`
MMP rule. The two relationships remain separate.

## What is implemented

- Child-first rationale records with evidence spans and explicitness classes
- RDKit single-cut fragmentation in both allowable orientations
- Explicit-H indexing, so common `H→R` moves are not lost
- Reason/transform alignment for named added and removed features
- Ranked top-k same-paper comparators with every decomposition witness retained
- Separate author-explicit scaffold relationships
- Censoring-aware assay comparisons that preserve `<` and `>` bounds
- Independent confidence fields for evidence, resolution, inference, assay
  comparability, and outcomes
- Frozen inputs, outputs, hashes, tests, and a reproducible ChEMBL query

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
python -m unittest discover -s tests -v
reasoned-mmps build-pilot
```

The build writes:

- [`outputs/reasoned_moves.json`](outputs/reasoned_moves.json): the four-layer
  derived view
- [`outputs/candidate_edges.csv`](outputs/candidate_edges.csv): top structural
  witness per candidate parent
- [`outputs/decomposition_witnesses.jsonl`](outputs/decomposition_witnesses.jsonl):
  every retained exact context/fragment match
- [`outputs/outcome_comparisons.csv`](outputs/outcome_comparisons.csv):
  censoring-aware assay deltas
- [`outputs/run_manifest.json`](outputs/run_manifest.json): versions, input
  hashes, counts, and guardrails

## Data model

The primary facts are independently versioned layers, not one giant confidence
score:

1. `evidence_span`: source, access quality, line range, and entailment
2. `design_assertion`: child, stated change, property, direction, and class
3. `parent_candidate`: context, transform, provenance tier, component scores
4. `measurement`: raw relation/value/unit and full assay context
5. `outcome_comparison`: comparability, bounded delta, and per-objective result

See [`schemas/reason_assertion.schema.json`](schemas/reason_assertion.schema.json)
and the extraction contract in
[`prompts/extract_child_reasons.md`](prompts/extract_child_reasons.md).

## Paperclip and ChEMBL roles

Paperclip supplies evidence-addressable literature and the ChEMBL SQL surface.
The query used for the pilot is preserved in
[`queries/photo_clenbuterol_chembl.sql`](queries/photo_clenbuterol_chembl.sql).
ChEMBL supplies structures, document linkage, assays, and measurements; it is
not treated as a source of medicinal-chemistry intent.

Reason extraction is human-reviewed in this first pilot. That is deliberate:
an unconstrained model can easily turn retrospective SAR into prospective
intent, invent individual rationales from a series-level statement, or promote
the closest structure into a fictional historical parent.

## Current limits

- One paper and two rationale assertions—not a benchmark yet
- Candidate generation currently stops at the same-paper universe
- Heuristic ranking scores are transparent but uncalibrated
- No tautomer/salt standardization policy beyond deterministic RDKit parsing
- No claim that compound numbering reflects synthesis order
- Outcome success is per objective; there is no universal “potency won” label

Next comes an adjudicated rationale-episode benchmark, then expansion through
cited predecessors, same-document ChEMBL compounds, and global ChEMBL
neighbors. Evaluation will keep reason extraction, parent Recall@k,
reason-transform alignment, assay joining, and end-to-end strict precision
separate.

## Licensing

No open-source license has been selected yet. Public visibility does not grant
reuse rights; add a license once the intended terms are decided.
