"""Exact single-cut MMP enumeration and reason-aware candidate ranking.

The comparison is analytic, not genealogical: an inferred candidate is a
defensible matched comparator unless the paper independently states lineage.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections import defaultdict
from dataclasses import asdict, dataclass
from typing import Iterable

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, rdMMPA


FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


@dataclass(frozen=True)
class FragmentRecord:
    context: str
    variable: str
    variable_heavy_atoms: int
    variable_fraction: float
    hydrogen_change: bool = False


def stable_id(prefix: str, value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return f"{prefix}_{hashlib.sha256(payload.encode()).hexdigest()[:16]}"


def canonical_smiles(smiles: str) -> str:
    """Parse, sanitize, and round-trip an isomeric SMILES deterministically."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _canonical_fragment(fragment: str) -> str:
    mol = Chem.MolFromSmiles(fragment)
    if mol is None:
        raise ValueError(f"Invalid fragment SMILES: {fragment}")
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _heavy_atoms(smiles: str) -> int:
    if smiles == "[*:1][H]":
        return 0
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumHeavyAtoms() if mol else 0


def _hydrogen_contexts(mol: Chem.Mol) -> Iterable[str]:
    """Yield contexts that let an implicit-H molecule match an H→R move."""
    seen: set[str] = set()
    for atom in mol.GetAtoms():
        if atom.GetTotalNumHs() < 1:
            continue
        editable = Chem.RWMol(mol)
        dummy = Chem.Atom(0)
        dummy.SetAtomMapNum(1)
        dummy_idx = editable.AddAtom(dummy)
        editable.AddBond(atom.GetIdx(), dummy_idx, Chem.BondType.SINGLE)
        context_mol = editable.GetMol()
        try:
            Chem.SanitizeMol(context_mol)
        except Exception:
            continue
        context = Chem.MolToSmiles(
            context_mol, canonical=True, isomericSmiles=True
        )
        if context not in seen:
            seen.add(context)
            yield context


def fragment_records(
    smiles: str,
    *,
    max_variable_fraction: float = 0.30,
    max_variable_heavy_atoms: int = 12,
    include_hydrogen_changes: bool = True,
) -> list[FragmentRecord]:
    """Enumerate both allowable single-cut orientations plus implicit-H cuts."""
    standardized = canonical_smiles(smiles)
    mol = Chem.MolFromSmiles(standardized)
    total_heavy = mol.GetNumHeavyAtoms()
    records: set[FragmentRecord] = set()

    for _, chains in rdMMPA.FragmentMol(
        mol, maxCuts=1, resultsAsMols=False
    ):
        parts = chains.split(".")
        if len(parts) != 2:
            continue
        canonical_parts = [_canonical_fragment(part) for part in parts]
        # Index both orientations when the proposed variable satisfies the
        # predeclared size rule. This avoids assuming that the larger component
        # must always be the context.
        for variable, context in (
            (canonical_parts[0], canonical_parts[1]),
            (canonical_parts[1], canonical_parts[0]),
        ):
            variable_heavy = _heavy_atoms(variable)
            fraction = variable_heavy / max(total_heavy, 1)
            if variable_heavy > max_variable_heavy_atoms:
                continue
            if fraction > max_variable_fraction:
                continue
            records.add(
                FragmentRecord(
                    context=context,
                    variable=variable,
                    variable_heavy_atoms=variable_heavy,
                    variable_fraction=round(fraction, 6),
                )
            )

    if include_hydrogen_changes:
        for context in _hydrogen_contexts(mol):
            records.add(
                FragmentRecord(
                    context=context,
                    variable="[*:1][H]",
                    variable_heavy_atoms=0,
                    variable_fraction=0.0,
                    hydrogen_change=True,
                )
            )
    return sorted(records, key=lambda r: (r.context, r.variable))


def _feature_present(fragment: str, feature: str) -> bool:
    normalized = feature.lower().replace("-", "_")
    if normalized in {"hydrogen", "h"}:
        return fragment == "[*:1][H]"
    if fragment == "[*:1][H]":
        return False
    mol = Chem.MolFromSmiles(fragment)
    if mol is None:
        return False
    atoms = {atom.GetSymbol() for atom in mol.GetAtoms()}
    if normalized in {"chloro", "chlorine", "cl"}:
        return "Cl" in atoms
    if normalized in {"fluoro", "fluorine", "f"}:
        return "F" in atoms
    if normalized in {"bromo", "bromine", "br"}:
        return "Br" in atoms
    if normalized in {"iodo", "iodine", "i"}:
        return "I" in atoms
    patterns = {
        "cyano": "C#N",
        "nitrile": "C#N",
        "methyl": "[CH3]",
        "hydroxyl": "[OH]",
        "amino": "[NH2]",
        "phenylazo": "N=Nc1ccccc1",
    }
    pattern = patterns.get(normalized)
    query = Chem.MolFromSmarts(pattern) if pattern else None
    return bool(query and mol.HasSubstructMatch(query))


def feature_alignment(
    parent_fragment: str,
    child_fragment: str,
    required_removed: Iterable[str],
    required_added: Iterable[str],
) -> tuple[float | None, list[dict]]:
    """Score only normalized feature constraints asserted during extraction."""
    checks: list[dict] = []
    for feature in required_removed:
        checks.append(
            {
                "side": "parent",
                "feature": feature,
                "matched": _feature_present(parent_fragment, feature),
            }
        )
    for feature in required_added:
        checks.append(
            {
                "side": "child",
                "feature": feature,
                "matched": _feature_present(child_fragment, feature),
            }
        )
    if not checks:
        return None, []
    return sum(check["matched"] for check in checks) / len(checks), checks


def _series_key(label: str) -> str | None:
    match = re.search(r"\b(\d+)[a-z]?\b", label.lower())
    return match.group(1) if match else None


def _tanimoto(smiles_a: str, smiles_b: str) -> float:
    mol_a = Chem.MolFromSmiles(smiles_a)
    mol_b = Chem.MolFromSmiles(smiles_b)
    return float(
        DataStructs.TanimotoSimilarity(FPGEN.GetFingerprint(mol_a), FPGEN.GetFingerprint(mol_b))
    )


def infer_parent_candidates(
    child: dict,
    compounds: Iterable[dict],
    reason: dict,
    *,
    max_variable_fraction: float = 0.30,
    max_variable_heavy_atoms: int = 12,
) -> list[dict]:
    """Return all reason-child MMP decompositions, ranked without claiming lineage."""
    child_records = fragment_records(
        child["canonical_smiles"],
        max_variable_fraction=max_variable_fraction,
        max_variable_heavy_atoms=max_variable_heavy_atoms,
    )
    child_by_context: dict[str, list[FragmentRecord]] = defaultdict(list)
    for record in child_records:
        child_by_context[record.context].append(record)

    candidates: list[dict] = []
    child_heavy = Chem.MolFromSmiles(child["canonical_smiles"]).GetNumHeavyAtoms()
    removed = reason.get("required_removed_features", [])
    added = reason.get("required_added_features", [])

    for parent in compounds:
        if parent["chembl_id"] == child["chembl_id"]:
            continue
        parent_records = fragment_records(
            parent["canonical_smiles"],
            max_variable_fraction=max_variable_fraction,
            max_variable_heavy_atoms=max_variable_heavy_atoms,
        )
        parent_heavy = Chem.MolFromSmiles(parent["canonical_smiles"]).GetNumHeavyAtoms()
        for parent_record in parent_records:
            for child_record in child_by_context.get(parent_record.context, []):
                if parent_record.variable == child_record.variable:
                    continue
                alignment, alignment_checks = feature_alignment(
                    parent_record.variable, child_record.variable, removed, added
                )
                context_heavy = _heavy_atoms(parent_record.context)
                retained_core_fraction = context_heavy / max(parent_heavy, child_heavy, 1)
                max_variable = max(
                    parent_record.variable_heavy_atoms,
                    child_record.variable_heavy_atoms,
                )
                compactness = 1.0 - min(
                    max_variable / max(max_variable_heavy_atoms, 1), 1.0
                )
                same_series = (
                    _series_key(parent["paper_label"]) is not None
                    and _series_key(parent["paper_label"])
                    == _series_key(child["paper_label"])
                )
                similarity = _tanimoto(
                    parent["canonical_smiles"], child["canonical_smiles"]
                )
                # Feature alignment dominates only when text names a structural
                # change. The individual components remain available for audit.
                score = (
                    0.35 * retained_core_fraction
                    + 0.20 * compactness
                    + 0.10 * float(same_series)
                    + 0.05 * similarity
                    + 0.30 * (alignment if alignment is not None else 0.0)
                )
                edge_semantics = (
                    "reason_constrained_mmp_comparator"
                    if alignment is not None and alignment == 1.0
                    else "structure_only_mmp_comparator"
                )
                transformation = (
                    f"{parent_record.variable}>>{child_record.variable}"
                )
                row = {
                    "candidate_id": stable_id(
                        "cand",
                        [
                            reason["reason_id"],
                            parent["chembl_id"],
                            child["chembl_id"],
                            parent_record.context,
                            transformation,
                        ],
                    ),
                    "reason_id": reason["reason_id"],
                    "parent_chembl_id": parent["chembl_id"],
                    "parent_label": parent["paper_label"],
                    "child_chembl_id": child["chembl_id"],
                    "child_label": child["paper_label"],
                    "edge_semantics": edge_semantics,
                    "historical_parent_claim": False,
                    "context": parent_record.context,
                    "parent_fragment": parent_record.variable,
                    "child_fragment": child_record.variable,
                    "transformation": transformation,
                    "cut_count": 1,
                    "valid_single_cut_mmp": True,
                    "hydrogen_change": (
                        parent_record.hydrogen_change or child_record.hydrogen_change
                    ),
                    "scores": {
                        "reason_transform_alignment": alignment,
                        "retained_core_fraction": round(retained_core_fraction, 6),
                        "edit_compactness": round(compactness, 6),
                        "same_numbered_series": float(same_series),
                        "whole_molecule_tanimoto": round(similarity, 6),
                        "ranking_score_uncalibrated": round(score, 6),
                    },
                    "alignment_checks": alignment_checks,
                    "algorithm": "single_cut_rdkit_mmpa+h_v1",
                }
                candidates.append(row)

    candidates.sort(
        key=lambda row: (
            row["scores"]["ranking_score_uncalibrated"],
            row["scores"]["retained_core_fraction"],
            -max(
                _heavy_atoms(row["parent_fragment"]),
                _heavy_atoms(row["child_fragment"]),
            ),
            row["parent_chembl_id"],
        ),
        reverse=True,
    )
    for rank, row in enumerate(candidates, start=1):
        row["rank"] = rank
    return candidates


def fragment_record_dicts(smiles: str) -> list[dict]:
    return [asdict(record) for record in fragment_records(smiles)]


def best_candidate_per_parent(decompositions: Iterable[dict]) -> list[dict]:
    """Collapse decompositions only after retaining the full witness table."""
    best: dict[str, dict] = {}
    for row in decompositions:
        best.setdefault(row["parent_chembl_id"], row)
    ranked = list(best.values())
    ranked.sort(
        key=lambda row: (
            row["scores"]["ranking_score_uncalibrated"],
            row["scores"]["retained_core_fraction"],
            row["parent_chembl_id"],
        ),
        reverse=True,
    )
    output = []
    for rank, row in enumerate(ranked, start=1):
        copy = dict(row)
        copy["parent_rank"] = rank
        output.append(copy)
    return output
