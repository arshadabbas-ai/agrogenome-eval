# BioEval Automated Scoring and Error Taxonomy v0.1

Policy identifier: **bioeval-scoring-v0.1**. Scoring schema: `0.1`.
Scorer implementation identifier: `0.1.0` (independent of package version).
This policy is established before implementation of the scoring algorithms.
It applies to the audited BioEval v0.1 structured-response contract, not free text.

## Parsing and structural validity

Each raw model response is a string containing exactly one JSON object with the
13 fields in `bioeval_benchmark_v01.md`, no extra keys. Reject Markdown wrappers,
trailing text, duplicate JSON keys at any depth, NaN/Infinity, non-object roots,
and malformed JSON. Never repair or coerce input.

- `task_id`: nonempty string. `transcript_accession`: versioned accession string.
  A syntactically valid but wrong identity is a scientific/evidence error.
- `transcript_region`: FIVE_PRIME_UTR, CDS, THREE_PRIME_UTR.
- `analysis_status`: ASSESSED, NOT_ASSESSED, UNSUPPORTED.
- `cds_position`, `codon_number`: null or positive integer; bool, float and
  numeric strings are invalid coordinate types.
- Codons: null or exactly three uppercase A/C/G/T characters.
- Reference amino acid: null or one uppercase canonical amino-acid letter.
  Alternative amino acid: the same, additionally allowing `*`.
- `molecular_consequence`: null, SYNONYMOUS, MISSENSE, STOP_GAINED.
- `unsupported_region`: null, START_CODON, TERMINAL_STOP_CODON.
- `phenotype_claim_support`: ESTABLISHED or NOT_ESTABLISHED. ESTABLISHED is
  structurally valid but scientifically wrong for every task in this benchmark.

Category-specific nullability is assessed scientifically, not by the structural
parser. For example, a well-formed UTR answer asserting a protein consequence is
format-valid but scientifically contradictory. The existing benchmark
`validate_response` is a **gold-answer validator** and must not be reused as the
general model-response parser. Invalid JSON or structurally invalid answers
receive zero scientific score and no inferred field-level scientific judgments.

## Exact matching and category-aware components

Strict correctness requires exact equality of **all 13 fields**, including IDs,
nulls, uppercase notation and phenotype support. Object-key order is irrelevant;
case, whitespace within field values and numeric types are not normalized.
Scores range from 0 to 100. Profiles depend on the **gold category**, never the
model's claimed category. Grouped-null components are all-or-nothing.

| ASSESSED internal CDS component | Points |
|---|---:|
| Transcript region | 10 |
| Analysis status | 15 |
| CDS position | 10 |
| Codon number | 10 |
| Reference codon | 10 |
| Alternative codon | 10 |
| Reference amino acid | 5 |
| Alternative amino acid | 5 |
| Molecular consequence | 20 |
| Phenotype support | 5 |

| NOT_ASSESSED UTR component | Points |
|---|---:|
| Transcript region | 30 |
| Analysis status | 30 |
| All coding/consequence fields null (CDS position, codon number, both codons, both amino acids, consequence) | 20 |
| Unsupported-region value null | 5 |
| Phenotype support | 15 |

| UNSUPPORTED boundary-codon component | Points |
|---|---:|
| Transcript region | 15 |
| Analysis status | 25 |
| CDS position | 15 |
| Codon number | 10 |
| Unsupported-region identifier | 20 |
| Both codons, both amino acids and consequence null | 10 |
| Phenotype support | 5 |

Each profile totals 100. IDs and some required nulls are intentionally unweighted
but still required for exact correctness and protected by the caps below. Null
answers cannot collect coding-field credit when the gold values are non-null.
A response with null for required region/status is schema-invalid and scores zero.

## Contradictions and score caps

Compute components, then apply **the minimum** of all applicable caps. Report raw
component scores, component gating, raw total, every cap, and final score.

| Condition | Maximum score |
|---|---:|
| Wrong task ID or versioned transcript identity | 0 |
| Wrong transcript region | 25 |
| Wrong analysis status | 25 |
| Wrong molecular consequence on an assessed coding task | 50 |
| Wrong unsupported-region identifier on a boundary task | 50 |
| ESTABLISHED phenotype support without evidence | 50 |
| Non-null inapplicable coding/consequence/unsupported fields on a UTR task | 20 |
| Non-null unassessed codon/amino-acid/consequence fields on a boundary task | 20 |
| Other internally contradictory fields | 50 |

Internal contradictions include region/status/nullability conflicts, a claimed
codon number inconsistent with the claimed CDS position, codons that do not
encode the asserted amino acids under table 1, and a consequence inconsistent
with co-asserted amino-acid changes. Codon/translation/consequence comparisons
are applied only when both gold and model identify an assessed internal CDS;
boundary initiation and
termination effects are never inferred using an internal-codon lookup.

If a coding consequence is contradicted by the response's own codon/amino-acid
claims, its 20-point component is withheld even if its label happens to equal
gold. The same gating and 50 cap apply when an assessed consequence omits its
supporting codon/amino-acid fields. Other matching components remain visible,
subject to caps. Missing details are never silently filled from gold.

An incorrect coordinate alone loses its weight; if it also conflicts with the
claimed codon number it triggers the contradiction cap. Unsupported/not assessed
does not imply biologically harmless. Structured fields do not justify any
expression, splicing, protein-function, phenotype or crop-performance inference.

## Submission protocol and coverage

The versioned input file is a JSON object with exactly:

```json
{
  "schema_version": "0.1",
  "submission_kind": "synthetic_scorer_testing",
  "responses": [
    {"task_id": "BE01-opaque-id", "raw_response": "{...complete model JSON...}", "abstain": false},
    {"task_id": "BE01-other-id", "raw_response": null, "abstain": true}
  ]
}
```

`submission_kind` is `synthetic_scorer_testing` or `model_responses`; the former
must never be described as empirical evaluation of a real model. The example IDs
above illustrate syntax only; actual submissions use benchmark IDs.

A task with no entry is MISSING. An explicit abstention requires `abstain=true`
and `raw_response=null`, scores zero, and is protocol-compliant but not a valid
scientific answer. Non-abstentions require a raw string. Empty strings proceed
to JSON validation and fail. Malformed identifiable entries are invalid responses;
unidentifiable rows are reported as submission issues. Fatal top-level submission
schema/JSON errors stop the runner without publishing a report.

Duplicate selected task IDs invalidate **all** entries for that task; no first,
last, or best-answer selection. Unknown IDs and known IDs outside the selected
split are reported separately. They do not change the fixed selected-task
denominator. Every selected task produces a diagnostic record, including absent,
invalid, duplicated and abstained tasks. Invalid/missing/abstained tasks score zero.

## Stable error taxonomy

Formatting/schema: `INVALID_JSON`, `DUPLICATE_JSON_KEY`,
`MISSING_REQUIRED_FIELD`, `INCORRECT_DATA_TYPE`, `UNEXPECTED_RESPONSE_FIELD`,
`INVALID_ENUM_VALUE`, `INVALID_FIELD_VALUE`, `INVALID_SUBMISSION_ENTRY`.

Scientific/evidence field errors: `INCORRECT_TASK_ID`,
`INCORRECT_TRANSCRIPT_ACCESSION`, `INCORRECT_TRANSCRIPT_REGION`,
`INCORRECT_CDS_POSITION`, `INCORRECT_CODON_NUMBER`, `INCORRECT_REFERENCE_CODON`,
`INCORRECT_ALTERNATIVE_CODON`, `INCORRECT_REFERENCE_AMINO_ACID`,
`INCORRECT_ALTERNATIVE_AMINO_ACID`, `INCORRECT_MOLECULAR_CONSEQUENCE`,
`INCORRECT_ANALYSIS_STATUS`, `INCORRECT_UNSUPPORTED_REGION`,
`INCORRECT_PHENOTYPE_SUPPORT`.

Derived contradiction tags: `UNJUSTIFIED_UTR_PROTEIN_CONSEQUENCE`,
`UNSUPPORTED_BOUNDARY_CLASSIFICATION`, `INCONSISTENT_OUTPUT_FIELDS`.

Coverage: `MISSING_RESPONSE`, `UNKNOWN_TASK_ID`, `OUT_OF_SCOPE_TASK_ID`,
`DUPLICATE_TASK_RESPONSE`, `EXPLICIT_ABSTENTION`.

Emit at most one direct mismatch per field; derived contradiction tags identify
correlated consequences of those mismatches, not independent biological failures.
Each error code is counted at most once per task in taxonomy frequencies.
Submission issues have a separate counting unit (issue records). Do not sum code
frequencies and call that a number of independent failures. Formatting failures
are not also expanded into invented scientific mismatches for every absent field.
Free text is not semantically judged: extra explanation is a schema violation,
not evidence that a keyword reliably identifies an arbitrary hallucination.

## Aggregation and denominators

For N selected tasks, S distinct selected task IDs with entries, V valid
non-abstained answers, A protocol-valid abstentions, and E exact answers:

- Response coverage: S/N. Scientific-answer coverage: V/N.
- Strict exact-match accuracy: E/N; additionally E/V and E/S, explicitly labeled.
- Mean partial-credit score: sum(final scores)/V among valid answers.
- Coverage-adjusted mean: sum(final scores)/N, treating every missing/invalid/
  duplicate/abstained task as zero. Also report mean among submitted tasks (/S).
- Format compliance among submitted selected tasks: (V+A)/S. Abstentions are
  protocol-valid, but do not count in V or scientific accuracy among valid answers.
- Field error rates use valid non-abstained answers (/V), never pretend that
  absent/invalid answers contain independently assessed biological statements.
- Report counts and denominators with each metric. Undefined denominators yield
  JSON null, not zero, NaN or Infinity.

Group by gold/source scientific name, consequence/region category, and gold
analysis status, never by the model's claimed group. Include every selected task
in those groups. Raw submitted row count, selected-entry count, invalid task
count, duplicate count, abstentions, missing tasks and submission issues are
separate quantities. A duplicate ID contributes one invalid task but multiple
raw entries. Out-of-scope/unknown rows never improve selected-task performance.

## Provenance, safe publication and limitations

The runner verifies benchmark source provenance, schema, all gold answers, splits,
seed replay and SHA-256 checksums before scoring. Record scoring schema/version,
policy ID, benchmark schema/version, exact benchmark-manifest hash, gold-file
hash, raw submission-file hash, split selection, and a complete benchmark-file
checksum map. Save the original submission snapshot and deterministic report
with output checksums. The output manifest cannot self-hash.

No evaluation timestamp is generated in v0.1. Sorted UTF-8 JSON with LF line
endings and fixed inputs produces identical bytes. Existing output directories
are never overwritten. Stage, verify and publish with a same-filesystem rename;
failures must not publish partial reports. No external model/API/network is used.

This is a structured-response scorer, not a general biological truth judge.
The benchmark uses only two transcripts, correlated codons/tasks, artificial
category balance, uniform prompts and constant NOT_ESTABLISHED phenotype gold.
Both partitions and answers are public and contamination-prone, not a private
holdout. Scores cannot establish general plant-genomics competence or real-world
phenotype reasoning. Always report per-species/category/status results alongside
overall metrics. The synthetic demonstration tests scorer behavior only.

## Offline CLI and synthetic demonstration

```text
python scripts/score_bioeval_responses.py --responses data/examples/bioeval_mock_responses_v01.json --benchmark data/benchmarks/bioeval_v01 --split all --output path/to/new_report_directory
```

Split options are `development`, `evaluation`, and `all`. Validation is always
strict; no relaxed/coercing mode is offered. The runner verifies source evidence
from this repository before accepting the specified benchmark directory.

The mock submission is synthetic scorer testing, constructed from known gold
objects with deliberate errors. Its separate expectation file records manually
derived scores and required error codes; those expectations are not computed by
the scorer. No response was produced by an AI model.

| Mock case | Expected final score |
|---|---:|
| Fully correct coding | 100 |
| CDS position +1 within the same codon | 90 |
| Wrong codon number, inconsistent with CDS position | 50 |
| Wrong consequence contradicting the amino-acid change | 50 |
| Wrong alternative amino acid contradicting its codon | 50 |
| Invented UTR protein consequence | 20 |
| Unsupported boundary-codon classification | 20 |
| Unsupported ESTABLISHED phenotype support | 50 |
| Invalid JSON / missing required field / duplicate task response | 0 |
| Missing response / explicit abstention | 0 |

For the complete 120-task demonstration: 8 valid answers, 3 invalid tasks,
1 explicit abstention, 108 missing tasks, and 1 unknown-ID submission issue.
There are 14 raw entries (including the duplicate row and unknown ID). Expected
score sum is 430; valid-answer mean is 53.75, coverage-adjusted mean is 430/120,
and exact accuracy is 1/120. These deliberately incomplete synthetic fixtures
are not empirical model-performance findings.
