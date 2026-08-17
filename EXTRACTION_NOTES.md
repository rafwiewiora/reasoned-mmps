# Extraction audit for the five-paper corpus

This pass is comprehensive **within the medicinal-chemistry design narrative of
the five selected papers**. The unit is an evidence-backed decision episode,
not a compound row and not every sentence containing an SAR result.

Included:

- prospective compound and series design rationales;
- explicit feature choices based on prior SAR;
- metabolite-guided and physicochemical-property-guided path choices;
- explicit mechanistic probes made to test a medicinal-chemistry hypothesis;
- retrospective explanations, labeled as retrospective;
- every compound explicitly named as a member of an included series.

Excluded from the reasoned-MMP episode count:

- synthesis reagent or route choices made only for yield, safety, or scale;
- assay-operational choices such as illumination wavelength or dose selection;
- results-only SAR observations with no design or experiment rationale;
- progression/termination decisions that do not define a new structural move;
- implicit rationales inferred only from compound numbering.

## Paper coverage

| Paper | Included episodes | Named reason compounds | Notes |
|---|---:|---:|---|
| Photo-clenbuterol | 3 | 5 | Captures the azoextension series, the 12b geometry hypothesis, and the 12e chlorocyano choice. Selection of 12b as the later key compound and selection of 360 nm for assays are progression/operational choices, not new molecular moves. |
| PF-06815189 | 1 | 1 | Compounds 5 and 6 are oxidation products identified alongside 2; the paper does not give them independent medicinal-chemistry design rationales. Synthetic reagent choices for scale-up are route decisions and are excluded. |
| DOS antimalarial | 12 | 36 | All ChEMBL-resolved analogues 3-40 are connected through explicitly stated appendage or core-series objectives. Nested, more specific episodes preserve oxetane, isoxazole, steric-bulk, lower-basicity, pyridyl, des-methyl, heteroatom, ring-contraction, and matched-R1 reasoning. The decision to abandon the des-urea branch and return to 5 is audited as a campaign decision but is not counted as a new structural move. |
| Bosutinib analogue | 3 | 2 | Separates the prospective design of 9, the 17→18 fragment probe, and the retrospective pKa explanation. These are deliberately not collapsed into one hindsight-contaminated assertion. |
| BIIB129 | 13 | 21 | Captures the CNS-oriented minimum pharmacophore, sp3-linker optimization, Asn484 hypothesis, metabolite-guided amino-ether campaign, explicit fluorine attempt, fluorine-effect explanation, bicyclic and docking-derived linker series, four-bond geometry probe, cyclobutyl bridge SAR, fused-ring deconstruction, α-methyl SAR transfer, and the 25→27 conformational explanation. Unnamed failed warheads and the final 25/27 progression choice are excluded because they do not define a new resolved structural move. |

### BIIB129 source reconciliation

The BIIB129 ChEMBL record reports whole-blood CD69 numeric values with `nM`
as the standardized unit, while the article text and tables report the same
numbers in `μM`. The curated measurement layer preserves the paper's `μM`
units, retains the ChEMBL activity identifiers, and records the discrepancy in
`source_note`. This prevents a silent thousand-fold absolute-unit error while
preserving exact within-assay comparisons.

## Denominators

- `rationale_episodes`: independently evidenced design or explanation episodes;
- `reason_bearing_compounds`: union of the compounds named within those
  episodes, including series members;
- `resolved_structures`: all ChEMBL structures loaded as the same-paper
  candidate universe, whether or not they carry a rationale;
- `unique_outcome_pairs`: unique parent/child measurement pairs;
- `outcome_comparisons`: episode-specific links to those pairs; a pair can
  support both a prospective rationale and a retrospective explanation.

The machine-readable per-paper denominators are regenerated in
[`outputs/run_manifest.json`](outputs/run_manifest.json).
