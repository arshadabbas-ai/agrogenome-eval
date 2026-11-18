# BioEval Controlled Pilot Evaluation v0.1

Pilot evaluation schema `0.1`, evaluator `0.1.0`, uses the existing BioEval scorer
and the existing capture schema. It joins verified questions, preserved execution
evidence, primary-response conversion and scientific scoring. The package
version remains `0.1.0`.

## Fixed population and preserved selection

The primary denominator is **14 selected development tasks**. Pilot seed is
20261112. The sorted-ID selection hash is
`aa134cfd42bfd1f9d87387f24cd3caf1d20a7c846ef405fc5df1798f78a2b0cb`.
The pipeline requires the original export and replays its sampling against the
verified benchmark. Another seed, stage, partition or selected-ID population
fails verification.

There are seven maize and seven wheat tasks, with one task per species/category
cell across SYNONYMOUS, MISSENSE, STOP_GAINED, FIVE_PRIME_UTR, THREE_PRIME_UTR,
START_CODON and TERMINAL_STOP_CODON. Every pilot ID belongs to development.
Original task IDs, complete questions, source metadata, RefSeq accessions and
the exact original prompt-export manifest remain bound to the experiment.

Both benchmark splits and all gold answers are public. The pilot shares the
development population and its transcript/codon dependencies. Its exposure and
small size require exploratory interpretation. It provides no hidden holdout or
independent species-level generalization estimate.

## Capture audit and trust boundaries

The Step 21B audit verified metadata, parameter availability, offset-aware
timestamps, question matching, UTF-8 preservation, native-payload separation,
path containment, duplicate rejection, attempt history, refusal/abstention
handling and atomic-style output publication. Model text remains untrusted data.
Origin and arbitrary extraction relationships remain researcher attestations.

Two focused regressions reproduced a derived-package defect. A changed
submission or an invented attempt count could pass the former file verifier
after its output hash was updated. Corrected verification now replays the
preserved metadata, ledger and raw evidence, reconstructs first eligible
responses, and requires exact equality with the saved submission and validation
report. It verifies input selection, source metadata hashes and every referenced
evidence file. This is a compatibility-preserving validation extension. Capture
schema `0.1` and original package filenames are retained.

Hashes detect changes relative to recorded values. They do not authenticate a
claimed model provider. A party able to replace all source evidence and its
attestations could create another internally consistent package. Honest origin
documentation and independent source evidence remain necessary.

Windows prevented native symlink creation during the earlier capture tests.
Those tests retain simulated resolved-path containment checks. Linux CI should
exercise actual symlink resolution and Python 3.12 runtime behavior. This local
step uses Python 3.14.7 and includes Python 3.12 syntax checks.

## Primary attempts and mutually exclusive task outcomes

The capture policy selects the first completed eligible attempt by registered
attempt number. Completed malformed, refused, truncated and empty outputs remain
eligible. Explicit documented abstention is also eligible. Technical failure
before an answer and no-response events remain in the ledger and allow a later
eligible attempt to become primary. Subsequent completions remain secondary.
Selection never consults gold labels or scores.

Every selected task receives exactly one outcome.

| Pilot outcome | Existing scorer status | Treatment |
|---|---|---|
| VALID_SCORED_RESPONSE | VALID | Structurally valid answer, scientifically scored against gold |
| INVALID_SCORED_RESPONSE | INVALID | Original format-invalid text, zero scientific score |
| EXPLICIT_ABSTENTION | ABSTAINED | Documented abstention, zero score |
| MISSING_OR_NO_ELIGIBLE_RESPONSE | MISSING | Absent task record or no eligible response, zero for coverage adjustment |

These four counts sum to fourteen. Each task also records its selected attempt,
observed execution outcome and, for missing responses, whether there was no task
record or recorded execution without an eligible answer. Refusal and truncation
remain observable execution states. The scorer processes their original text
under the same strict format contract. A refusal is never inferred by keyword
matching or automatically converted to an abstention.

Unknown tasks, wrong-partition records, duplicate task/attempt identities,
inconsistent retries, invalid metadata and checksum mismatches block evaluation.
Extraneous records cannot enter another population or improve a denominator.
Accepted packages report zero extraneous records. Attempt counts and secondary
completions remain separate from selected-task outcome counts.

## Scientific scoring and denominators

`evaluate_pilot` calls the published
`evaluate_responses(tasks, gold_answers, submission, task_ids=selected_ids)`
function. It supplies exactly the fourteen fixed IDs. The scientific weights,
minimum applicable caps, exact matching, response validation and stable error
taxonomy remain those of the published scoring policy.

Let N=14, S be selected tasks with a primary submission row, V valid structured
answers, A abstentions, E exact answers and P the sum of final task scores.

| Metric | Definition |
|---|---|
| Response coverage | S/N, including submitted abstentions |
| Scientific-answer coverage | V/N |
| Strict exact accuracy | E/N |
| Exact accuracy among valid answers | E/V |
| Mean partial score among valid answers | P/V |
| Coverage-adjusted mean | P/N |
| Format compliance among submitted tasks | (V+A)/S |

All undefined ratios are JSON null. Field error rates use V as denominator.
Scientific categories, species and analysis-status groups use verified gold/source
metadata. Category counts are two tasks each, species counts seven each, and
analysis-status counts six ASSESSED, four NOT_ASSESSED and four UNSUPPORTED.
Each species/category cell contains one observation.

Captured-task coverage counts selected IDs with execution records, including
failed executions. Attempt count includes repeated events. Neither quantity
substitutes for S or V. Failures without output stay absent from submissions.
No placeholder answer is created. Formatting failures yield no invented
biological field judgments.

Direct field mismatches remain separate from derived contradiction tags.
Taxonomy counts each code at most once per task. Correlated tags cannot be
summed into independent model-error counts. The report retains wrong-region,
coordinate, codon, amino-acid, consequence, boundary, UTR and phenotype-support
diagnostics from the existing scorer. Capture refusal status accompanies those
diagnostics through execution metadata.

**The heuristic policy can retain 90 points for a coordinate-only error and 80
for a coherent wrong coordinate pair.** A high partial score must be reported
alongside exact correctness. Scores do not establish biological truth beyond
the structured benchmark evidence or assess arbitrary free-text explanations.
Category percentages are descriptive and highly unstable. The pipeline produces
no confidence intervals from these fourteen correlated questions.

## Offline execution and deterministic output

```text
python scripts/score_bioeval_pilot.py --prompt-export tmp/bioeval_pilot_v01 --capture-import tmp/bioeval_capture_synthetic_v01/derived --output tmp/bioeval_pilot_preflight_v01
```

The command uses the existing explicitly synthetic capture material. For future
real observations, pass their reviewed import directory. `--benchmark` defaults
to the repository's verified benchmark. The runner verifies the original
reference/CDS/variant evidence, complete benchmark, pilot export and capture
replay before scientific scoring. Existing output files/directories and leaf
symlinks receive overwrite protection. Output must be separate from source
directories. Staging, byte verification, rechecking inputs and same-filesystem
rename protect publication. Failed writes or renames clean up staging.

The output directory contains `report.json`, `experiment.snapshot.json`,
`capture_validation.snapshot.json`, `submission.snapshot.json`,
`pilot_export_manifest.snapshot.json`, `prompts.snapshot.json`,
`capture_import_manifest.snapshot.json` and `manifest.json`.
The capture import remains the preserved source of raw/native evidence. The
report records checksums for all of its input files and every scientific source
file, plus policy/scorer/evaluator versions, selection identity and policy hash.
The output manifest hashes the other seven files without self-hashing.

The report preserves complete experiment metadata, settings availability,
attestation limits and capture warnings. Source execution times are retained.
Evaluation adds no timestamp. Canonical sorted JSON with LF termination permits
byte-identical replay from the same verified inputs. Inputs remain local.

## Independently reconciled synthetic preflight

The Step 21B sample contains thirteen attempt records for eleven captured task
IDs. Ten IDs have eligible primary rows, including one explicit abstention.
Three IDs have no task record, and one recorded task has no eligible response.
Thus `10+3+1=14` selected tasks. There are two additional attempt events and one
secondary eligible completion. The full ledger contains three format-valid text
outputs, but one belongs to a secondary attempt and is excluded from scoring.

The primary task outcomes are **2 valid structured answers, 7 invalid answers,
1 abstention and 4 missing/no-eligible tasks**. They sum to fourteen. One valid
answer exactly matches its reference. The deliberately wrong structured answer
matches region (10), status (15) and phenotype-support (5), giving 30 points.
Its consequence is contradicted by its own amino-acid claims and receives no
consequence credit. Together the valid scores are `100+30=130`.

Expected exact accuracy is `1/14`, valid-answer exact accuracy `1/2`, valid-answer
mean `130/2=65`, coverage-adjusted mean `130/14=9.285714...`, response coverage
`10/14`, and format compliance `(2+1)/10=0.3`. These are synthetic end-to-end
preflight calculations. Real-model performance requires separately documented
actual executions.

## Research limitations

Public answers, shared transcripts/codons, artificial category balance, uniform
templates and constant phenotype-support labels permit memorization and
shortcuts. Model version visibility, UI transformations, hidden instructions,
memory and tools can limit experimental comparability. Source and output files
are loaded into memory. Rechecks detect changes but provide no hostile-writer
filesystem lock or independent model-origin authentication. Future records need
privacy and provenance review before public release. Functional, drought,
yield and broad plant-genome claims require additional evidence.

The [manual execution guide](bioeval_manual_model_pilot_guide_v01.md),
[capture specification](bioeval_response_capture_v01.md),
[empirical protocol](bioeval_empirical_protocol_v01.md) and
[published scoring policy](bioeval_scoring_v01.md) define the surrounding workflow.

## Local Step 21C-1 validation record

Python 3.14.7 completed **631 passing tests**, including the 592-test baseline
and 39 new pilot/audit cases. All 39 new tests also passed with socket networking
blocked. Python 3.12 syntax checks passed for the new pilot code and modified
capture interfaces. Python 3.12 runtime CI and actual Linux symlink verification
remain pending.

The documented CLI produced the ignored local
`tmp/bioeval_pilot_preflight_v01/` package from the original Step 21B synthetic
import. Its fourteen task outcomes, manual 130-point total and fixed-denominator
metrics agree with independent assertions. Tests verify byte-identical replay,
source/derived checksums, exact question preservation, synthetic labeling,
missing/abstained/failed populations, refusal/truncation handling, coordinate
partial credit, capture replay, tampering rejection and failure cleanup.
The original pilot export and Step 21B source/import evidence retain their
saved bytes. Existing tracked scientific files and the preserved root workflow
retain the verified baseline.


## Completed empirical pilot and release boundary

The [de-identified results summary](bioeval_pilot_results_v01.md) records the
completed fourteen-task development pilot. Historical synthetic/preflight counts
in this document retain their original testing context. Published aggregate
results distinguish JSON syntax, schema compliance and scientific grading.

Import and scoring outputs contain original response text, execution metadata
and researcher attestations. Keep these packages in private research storage
until their contents receive explicit publication review. The public summary
omits those records. Full empirical replay requires access to the preserved
private capture package. Repository tests reproduce the software behavior using
explicitly synthetic fixtures generated in temporary directories.
