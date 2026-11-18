# BioEval Empirical AI Model Evaluation Protocol v0.1

Protocol identifier `bioeval-empirical-protocol-v0.1` accompanies exporter
version `0.1.0` and export schema `0.1`. The scientific package remains `0.1.0`.
This protocol defines experimental design and offline question exports.
Execution takes place separately from the software. The completed fourteen-task
development pilot is described in the [public results summary](bioeval_pilot_results_v01.md).
Earlier local validation records below describe their historical implementation stages.

## Research questions and experimental objectives

Measure how an identified model follows the BioEval structured response contract,
maps transcript coordinates, interprets internal coding substitutions, respects
UTR/boundary scope, and represents phenotype evidence. Examine format compliance,
exact correctness, category-sensitive partial credit, coverage and error patterns.
Report model identity, experimental settings and the selected task population
alongside each result. The primary unit is a task response within a declared run.
Two transcripts and shared codons constrain the independence of those units.

A genuine model response must originate from an identifiable, independently
executed model through a recorded manual interface, local inference runtime or
API. Its execution record must support the connection between the actual prompt
and observed raw output. Gold-answer files and synthetic scorer fixtures serve
as evaluation/test materials. Researchers must never relabel their contents as
independently generated model responses. Existing synthetic demonstrations retain
their original designation and cannot establish empirical model performance.

## Experimental materials and information boundaries

| Material | Role and handling |
|---|---|
| Model-facing task | Original `task_id`, `input` and complete `question`, copied from verified `tasks.json` |
| Gold answer | Researcher-side reference from verified `gold_answers.json`, available to the scorer and internal pilot stratification |
| Export manifest | Researcher-side sampling/provenance record, including internal pilot category labels |
| Raw response | Original observed model content, stored before validation or scoring |
| Run metadata | Model identity, access method, settings, execution history and prompt provenance |
| Validation diagnostics | Strict parsing/schema results under the existing response contract |
| Scientific grade | Field comparisons, approved weights, caps and error taxonomy |
| Aggregated results | Run-level and group-level metrics with selected-task denominators and limitations |

Only `prompts.json` contains model-facing records. Supply each record's exact
`question` as the user message in a fresh context. Its `input` contains the same
source facts already represented by that question and supports independent
inspection. The record keeps the original task ID and schema. Keep the manifest,
gold files, examples and scorer results in researcher storage during execution.
Internal category labels belong exclusively to researcher-side audit metadata.
The original question provides response fields and selected consequence/status
labels. Its text omits accepted UTR enum spellings and the required one-letter
amino-acid notation. The completed pilot retains those original questions.
These disclosure gaps limit interpretation of strict schema failures. A revised
question/schema presentation requires a separately versioned experiment.
Preserving original instructions is part of the v0.1 benchmark contract.

## Model selection and run registration

Select models according to the research comparison being proposed, execution
access, identifiable version, available parameter controls and traceability.
Record provider/source and exact visible model name, checkpoint or release
identifier. Record a local checkpoint hash and runtime/library versions when
available. For hosted aliases, retain the visible alias, access time and any
provider response identifiers. Mark hidden or unavailable details explicitly.
Comparisons across access methods must disclose differences in tools, context,
system messages, sampling controls, interface transformations and version drift.

Before execution, register the run ID, research questions, model/version,
selected task IDs, export checksums, stage, prompt configuration, parameter
availability, repetitions, retry policy, stopping rule and primary metrics.
Use one completion per task for the initial pilot. Any later repetitions must
have a declared count and separate replicate IDs. Record all attempted runs,
including failed or abandoned ones. Report post hoc analyses as exploratory.

## Staged evaluation and deterministic pilot

| Stage | Source partition | Planned task population |
|---|---|---:|
| Pilot | Development | 14, one per species-by-category combination |
| Development evaluation | Development | A preregistered subset up to all 80 |
| Public evaluation | Evaluation | 40 |

The source contains maize NM_001305944.1, VALIDATED RefSeq, and wheat
XM_044542390.1, MODEL RefSeq. Each species has seven strata, comprising
SYNONYMOUS, MISSENSE, STOP_GAINED, FIVE_PRIME_UTR, THREE_PRIME_UTR,
START_CODON and TERMINAL_STOP_CODON. Source curation terms retain their distinct
meanings. All fourteen strata are available in the audited development split.
An absent stratum causes an explicit failure so a reduced pilot needs a
separately documented design decision.

Pilot seed is **20261112**, recorded independently from the benchmark generation
seed. Restrict candidates to development IDs, group internally by source
scientific name and gold category, sort strata lexicographically, sort eligible
IDs within each stratum, then call `random.Random(seed).choice` once per stratum.
Finally sort the fourteen selected records by their original opaque IDs.
Record seed, algorithm/version, selected IDs and eligible stratum counts.
Selection uses a local random generator. Replaying identical verified inputs
with this procedure produces the same batch. Test replay on each intended runtime.

The exporter provides complete development/evaluation batches and the 14-task
pilot. A smaller development experiment can predeclare a subset from the
80-record development export. Retain its complete intended ID list in the run
record before collecting responses. Pilot tasks are included in the development
population, so pilot/development scores share tasks and require that disclosure.
The public evaluation IDs are disjoint from development IDs.

Use the pilot to find procedural problems and characterize errors. Record any
prompt or execution changes in a new run specification before subsequent stages.
Keep public evaluation out of prompt-development feedback. Because its questions
and answers are already public, this separation describes experimental conduct
while answer exposure and model-training contamination remain possible.

## Prompt export and reproducible offline use

Run from the repository root after `python -m pip install -e ".[dev]"`.
Choose a new output directory for every export.

```text
python scripts/export_bioeval_model_prompts.py --partition development --pilot --seed 20261112 --output tmp/bioeval_pilot_v01
python scripts/export_bioeval_model_prompts.py --partition development --output tmp/bioeval_development_v01
python scripts/export_bioeval_model_prompts.py --partition evaluation --output tmp/bioeval_evaluation_v01
```

`--benchmark` accepts a benchmark directory and defaults to the repository's
`data/benchmarks/bioeval_v01`. Source evidence is verified against this checkout.
`--pilot` requires `--partition development`. `--seed` controls pilot selection.
Full-partition exports use ascending task-ID order and record a null pilot seed.

The exporter verifies the existing reference/CDS/variant evidence chain,
benchmark schemas, checksums, source identities, independent gold checks and
seed replay. It snapshots and rechecks benchmark bytes around validation and
before publication. Staging and same-filesystem rename publish complete outputs.
Existing files, directories and leaf symlinks receive overwrite protection.
Interrupted staging is cleaned up. The workflow uses local files throughout.

`prompts.json` reuses the existing `{"schema_version":"0.1","tasks":[...]}`
envelope. Every record has exactly the original `task_id`, `input`, `question`.
Question text and source input are preserved. Expected classifications, scores,
reference answers and internal pilot category metadata are excluded from this
payload. Export ordering is independent of the pilot stratum traversal order.

The separate `manifest.json` records protocol/exporter identifiers, partition,
stage, selected IDs, pilot seed and strata, benchmark/source SHA-256 values,
each question's UTF-8 SHA-256, and the exact `prompts.json` byte checksum.
The manifest hashes its payload without a circular self-hash. Export outputs
omit timestamps for byte replay. Actual execution timestamps belong in the run
records. Preserve saved bytes, including line endings, when sharing exports.
Checksums detect byte changes and require trusted provenance to establish identity.

## Capturing genuine experimental responses

For every run and task attempt, retain the following in researcher storage.

| Field group | Required contents |
|---|---|
| Model identity | Provider/source, exact visible model name, available version/checkpoint identifier, identity availability notes |
| Execution | Run/task/replicate/attempt IDs, start/end date-times in RFC 3339 with timezone, access method, execution environment and available request identifiers |
| Software provenance | Package, exporter and scorer versions, repository commit or reviewed source hashes, protocol identifier and scoring-policy identity |
| Settings | Temperature, model seed, token limits, sampling controls, tool configuration and other exposed settings |
| Setting availability | Distinguish explicitly set values, observed defaults, unavailable controls and unknown effective values. Numeric zero is a real value. Use null plus availability/reason for unknown values |
| Messages | Exact user prompt, exact known system/developer messages, roles and ordering, attachments/context and availability of hidden instructions |
| Output | Raw unmodified model text and available native response record, finish reason, refusal, truncation, output failure and any tool interactions |
| History | Retry reason, technical failure, elapsed time and every attempted completion, including discarded secondary attempts |
| Integrity | Export manifest/payload hashes, task ID, original question hash, exact sent-prompt hash, raw-response hash and encoding convention |

Use a fresh context per task. The baseline user message is the original complete
question. Predeclare any system/developer messages, tool permissions, retrieval
or wrappers and retain their exact available text. Record unavailable provider
instructions as unavailable. Keep gold answers outside execution contexts and
prohibit answer retrieval during the controlled experiment. Log deviations and
evaluate them separately. Parameter values and seed support must reflect actual
access, including interface defaults and hidden settings.

Save raw response content before parsing. Preserve line breaks, whitespace,
Markdown fences, refusal text and trailing text. Hash question, sent user text
and raw response as exact UTF-8 bytes without normalization. Record native
transport artifacts separately with their own checksums. For multi-message
requests, also save a deterministic UTF-8 JSON representation of the complete
known message list, documenting serialization and its checksum. Record any
manual copy/interface transformation that affects the observed text.

Predeclare technical retry conditions and bounds. Step 21B uses the first
completed eligible attempt as the primary observation. Failures before an
answer and no-response events remain recorded, and a later eligible completion
can be primary. Completed malformed, refused, truncated and empty outputs stay
eligible. Additional attempts have separate attempt IDs and secondary analyses.
Scientific wrong answers and formatting errors must never trigger selective
replacement of the primary observation. Preserve refusals and truncations as
observed outputs. Record transport failures even when there is no model text.
Credentials and API authorization headers belong outside research artifacts.

## Future conversion to the existing scorer

Future import code will convert observed content into the current submission
envelope with `schema_version="0.1"`, `submission_kind="model_responses"` and
a `responses` list. Each entry uses the associated `task_id`, original
`raw_response` string and `abstain=false`. Retain run metadata in a linked
experiment record because the scorer envelope accepts only its established keys.

Raw refusal/truncation content proceeds through normal strict validation and
can receive a zero score. A transport failure with no observed response remains
missing, with its cause retained in the execution ledger. An explicit protocol
abstention uses `abstain=true` and `raw_response=null` and requires an actual
recorded abstention decision. Never repair JSON, strip fences, invent empty
responses, convert biological labels, or substitute gold content. The capture
ledger retains all attempts. Each scoring submission contains only the declared
primary attempt per selected task, since duplicate task IDs invalidate that task.

After verifying the benchmark and captured inputs, the existing reusable API
`evaluate_responses(tasks, gold_answers, submission, task_ids=selected_ids)`
supports a pilot's fixed fourteen-task denominator or a predeclared development
subset. The existing scoring CLI selects entire development/evaluation splits.
Using `--split development` for only fourteen answers evaluates all eighty
tasks and counts the other sixty-six as missing. Keep these denominators explicit.
Step 21A implements question export. The subsequent capture contract and
offline importer are documented in [response capture v0.1](bioeval_response_capture_v01.md).

## Analysis and reporting requirements

Report selected, submitted, valid, invalid, abstained and missing counts.
Provide strict exact-match accuracy across all selected tasks, accuracy among
valid answers, valid-answer mean partial credit, coverage-adjusted mean score,
format compliance and error frequencies. Include per-species, category and
analysis-status groups, run/replicate identifiers, settings, exposure history
and failure counts. Apply the published scoring weights and minimum applicable
caps. A coordinate-only error can retain 90 points, and a coherent wrong
coordinate pair can retain 80. Distinguish these partial scores from exact answers.
Weights and caps are heuristic choices with no empirical calibration.

Use matched tasks and declared settings for model comparisons. Show raw counts
and avoid treating 120 correlated tasks as independent biological experiments.
Only two transcript clusters provide very limited support for uncertainty
estimation. The fourteen-task pilot has one observation per species/category
cell and cannot deliver precise population performance estimates. Any descriptive
interval must disclose its unit, assumptions and dependence limitations.

Both splits and every gold answer are public. Exposure, memorization and training
contamination make these experiments exploratory. Record known prior exposure
and uncertainty about training data. Shared codons/transcripts, deliberately
balanced categories, uniform templates and constant NOT_ESTABLISHED phenotype
labels permit shortcuts. The public evaluation partition is not a hidden or
contamination-resistant holdout. Generalization across plant genomes, molecular
function, drought tolerance or yield requires independent evidence and additional
experimental designs. This protocol measures controlled structured-response
performance under the declared settings.

Detailed contracts are retained in the [benchmark specification](bioeval_benchmark_v01.md),
[scoring policy](bioeval_scoring_v01.md) and [scoring audit](bioeval_scoring_v01_qa.md).

## Step 21A local validation record

Local Python 3.14.7 passed all 525 tests, comprising the 497-test baseline and
28 new exporter tests. All 28 new tests also passed with socket networking
blocked, including the CLI entry point. Python 3.12 syntax checks passed for
the exporter and its tests. Python 3.12 runtime validation of this new feature
will require subsequent CI execution.

The documented commands produced 14 pilot, 80 development and 40 evaluation
records in ignored local review directories. Tests confirmed one pilot task per
species/category, a pinned independent seed-selection fixture, complete task
preservation, partition isolation, payload/manifest separation, input/output
checksums, byte-identical replay, corruption rejection and failed-publication
cleanup. Empirical model execution and response collection remain subsequent
research activities.
