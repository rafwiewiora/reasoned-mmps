# Child-first medicinal-chemistry rationale extraction

Extract only claims supported by a cited evidence span. The primary record is a
reason-bearing **child compound**; a parent compound is optional.

For each assertion, return:

- the paper label and resolved structure/identifier for the child;
- whether the assertion is prospective intent, an explicit feature choice with
  an SAR basis, series-level intent, implied intent, retrospective explanation,
  generic SAR, or no rationale;
- a conservative paraphrase of the author reason;
- the intended property and direction, preserving `not_stated` when needed;
- any explicitly named added and removed structural features;
- the exact evidence path and line range, access quality, and entailment status;
- an explicit parent only if the source states one.

Do not infer a parent from compound numbering. Do not turn an observed benefit
into prospective intent. Do not assign a series rationale to every member
without marking it series-level. If only an abstract or supplement is present,
record the missing-text status instead of returning a negative extraction.

The output must validate against `schemas/reason_assertion.schema.json`.
