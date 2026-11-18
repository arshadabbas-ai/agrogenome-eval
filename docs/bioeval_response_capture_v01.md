# BioEval Model-Response Capture and Provenance v0.1

Capture schema `0.1` and implementation `0.1.0` bind locally supplied experiment
records to a verified BioEval prompt export. Package version remains `0.1.0`.
The workflow supports manual interaction, future API observations and future
local inference through the same evidence contract. Actual execution belongs to
the researcher. This implementation validates and preserves supplied artifacts.

## Provenance interpretation and privacy

Researcher-entered model identity and attestations do not independently prove
model origin. SHA-256 checksums demonstrate file integrity relative to recorded
values. Provider authentication would require additional independently verified
evidence. The report uses `researcher_supplied_origin_unverified` for model
captures and `synthetic_capture_system_testing` for labeled fixtures.

Manual copying preserves the bytes supplied to this system. UI rendering,
clipboard transformations and hidden interface instructions limit claims about
the original provider output. Record these limitations explicitly. Future API
captures preserve the native source payload separately from extracted answer
text. Each has its own path and hash. The extraction description and asserted
provider/source remain researcher-supplied provenance. The validator verifies
each saved byte sequence, without authenticating their claimed origin or an
arbitrary extraction procedure.

Keep credentials, authorization headers, access tokens, personal identifiers and
unnecessary private conversations outside capture packages. Review source
payloads for such material before research use. Automatic removal would alter
evidence, so this system performs no redaction. Retain only appropriate source
material and record capture limitations honestly.

## Source capture directory

A researcher prepares `experiment.json`, `records.json` and relative evidence
files. Each JSON root is an object. Duplicate keys, nonfinite numbers and invalid
UTF-8 fail validation. Unknown schema fields fail validation. Paths use portable
relative forward-slash notation within the capture root. Absolute paths,
Windows drives/streams, backslashes, parent traversal and symlink escapes are
rejected. Existing missing files, invalid checksums and inconsistent identities
block conversion. Model content is read as data and is never executed.

`experiment.json` has exactly these fields.

| Field | Contract |
|---|---|
| `capture_schema_version` | `"0.1"` |
| `experiment_id` | Nonempty researcher-assigned run identifier |
| `submission_kind` | `synthetic_scorer_testing` or `model_responses` |
| `benchmark_schema_version` | Exact benchmark schema version, currently `"0.1"` |
| `prompt_export` | Exact binding returned by `VerifiedPromptExport.identity` |
| `model` | `{source, name, version}` with meaningful source and exact visible model name |
| `access_method` | `manual_interface`, `api`, or `local_inference` |
| `execution_started_at` | Actual start timestamp with explicit UTC offset |
| `instructions` | Exactly `system` and `developer` instruction descriptors |
| `generation_parameters` | Availability descriptors, including temperature and seed |
| `tool_policy` | `{mode, details}`, with disabled/permitted/unavailable mode and nonempty details |
| `retry_policy` | `{selection: "first_completed_eligible", max_attempts: positive integer}` |
| `provenance` | `{basis, attestation, capture_limitations}`, all documented |

The prompt binding contains `manifest_sha256`, `protocol_id`, `partition`,
`stage` and `selection_id`. The selection identifier hashes the canonical JSON
array of sorted selected task IDs. The manifest hash covers its exact saved
bytes. This binds all fourteen pilot IDs to the development export, including
the pilot seed and independently replayed sampling selection.

`model.version` is `{availability, value, reason}`. A `recorded` version requires
nonempty text and null reason. An `unavailable` version requires null value and
an explicit reason. Never invent an unavailable version or checkpoint identifier.

Each instruction descriptor is `{availability, text, reason}`. `recorded`
requires exact available text and null reason. `absent` records explicit absence
with null text/reason. `unavailable` records hidden or inaccessible instructions
with null text and a nonempty reason. Record all known system/developer text.

Each generation parameter is `{availability, value, reason}`. `set` and
`observed_default` require the actual finite JSON scalar value and null reason.
`unavailable` requires null and a nonempty reason. Temperature must be a
nonnegative number, and seed a nonnegative integer when available. Boolean
coordinates or sampling seeds cannot stand in for integers. Zero temperature
is distinct from an unavailable temperature. Other exposed scalar settings may
be recorded under their actual names. Software/runtime versions and observed
defaults can be described in provenance limitations and the attestation.

For `model_responses`, provenance `basis` is `researcher_attestation`,
`provider_payload` or `local_execution_record`. Attestation and capture limitations
are mandatory. Synthetic fixtures require `synthetic_fixture` and remain labeled
`synthetic_scorer_testing` through conversion. Metadata completeness supports
reproducibility while origin remains unverified.

## Task records and exact response preservation

`records.json` contains exactly `capture_schema_version: "0.1"` and a `records`
array. Every execution event has exactly these fields.

| Field | Contract |
|---|---|
| `task_id` | Original exported BioEval task ID |
| `attempt` | Positive integer within the registered attempt bound |
| `retry_of` | Null for attempt one, otherwise immediately preceding attempt number |
| `retry_reason` | Null for attempt one, otherwise nonempty reason for the additional attempt |
| `timestamp` | Actual event timestamp, with UTC offset |
| `outcome` | completed/refused/truncated/abstained/no_response/execution_failure |
| `outcome_details` | Nonempty decision/failure description for every outcome except completed, where null is permitted |
| `prompt` | `{path, sha256, record_sha256}` binding exact user-question bytes and original complete task record |
| `response` | `{path, sha256, source_payload}` with nullable answer path/hash |
| `output_error` | Boolean recording an observed output error |
| `tool_use` | `{status, details}` with none/recorded/unavailable status |

Timestamps follow ISO 8601 with seconds and `Z` or an explicit `+HH:MM`/`-HH:MM`
offset. Each must be on or after experiment start. Attempt timestamps for each
task must be chronological. All supplied task attempts are preserved.

The prompt file must contain exactly the exported `question` encoded as UTF-8.
Its SHA-256 covers those bytes. `record_sha256` covers canonical JSON for the
entire original `{task_id, input, question}` record. Canonical JSON uses sorted
keys, indentation of two spaces, finite values and a final LF, encoded as UTF-8.
Task selection, complete question text and input facts must agree with the
original verified benchmark. Questions contain no appended answer material.

The answer file contains the actual text exactly as supplied, including leading
or trailing whitespace, case, Markdown fences, explanations, refusals and empty
output. Preserve UTF-8 bytes and hash them directly. JSON parsing happens only
after capture integrity verification. A malformed model answer remains a valid
capture when its execution and file provenance are consistent. Format diagnostics
are distinct from capture-integrity errors. Scientific scoring happens later.

`source_payload` is null or `{path, sha256, extraction}`. Preserve original native
API bytes here, including their original encoding, with an explicit description
of which answer field was extracted and how. Every API event with answer text
requires its source payload. Native payload and answer bytes are independently
preserved even if their encodings or structures differ. Source payloads may also
document technical failures. Their contents are never executed.

Tool status `none` requires null details. `recorded` or `unavailable` requires a
description. Recorded tool use conflicts with a disabled policy and blocks
conversion. Unknown tool visibility remains documented in the record and run
limitations.

## Outcomes, attempts and scorer mapping

The prespecified default selects the first completed eligible attempt by attempt
number, independently of gold labels, JSON validity or score. Eligible outcomes
are completed, refused, truncated and abstained. Empty completed output is also
eligible. Later attempts remain in the evidence package and cannot replace an
eligible primary observation. The validator requires contiguous attempts from
one, unique `(task_id, attempt)` identities and explicit retry relationships.

Technical failures before an answer and no-response events are retained as
ineligible events. A later eligible attempt may be primary. Record the reason
and registered attempt bound. Attempts after an eligible event are retained for
secondary investigation. A manual prompt/copy error that changes question bytes
is a capture-integrity failure and requires explicit correction of the provenance
or a separately documented run. It cannot silently be treated as a poor answer
and replaced after examining gold labels.

| Captured outcome | Required answer evidence | Submission mapping |
|---|---|---|
| completed | UTF-8 file, including a possible zero-byte file | Original text, abstain false |
| refused | Actual refusal text file | Original text, abstain false |
| truncated | Actual partial text file and outcome details | Original text, abstain false |
| abstained | Null answer path/hash and documented decision | Null raw_response, abstain true |
| no_response | Null answer path/hash and details | Omit unless a later eligible attempt exists |
| execution_failure | Null answer path/hash and failure details | Omit unless a later eligible attempt exists |

No record for a selected task yields `MISSING_TASK_RECORD`. Recorded attempts
without an eligible response yield `NO_RESPONSE_OBTAINED`. Both are reported
as coverage warnings. They preserve the selected population and missing-response
denominator. A duplicate task/attempt, unknown ID, wrong partition, mismatched
prompt, unavailable evidence or inconsistent history is a blocking integrity
error. Multiple attempts for one task are distinguished by attempt number.

Refusal text and explicit abstention remain separate. Invalid JSON, fenced
output and extra explanations are preserved without repair or biological
correction. An empty captured string remains submitted output. Missing output
stays absent. Capturing an arbitrary inner `task_id` in model JSON preserves that
scientific mismatch for the scorer. The associated execution task comes from
the validated outer record and prompt binding.

## Offline importer and output package

Prepare capture records using the actual observed model identity, parameters,
timestamps, messages and copied files. The reusable loader exposes exact record
and export hashes for preparing those records. Manual capture uses local text
files and requires neither API credits nor an external service connection.

```text
python scripts/import_bioeval_model_responses.py --prompt-export tmp/bioeval_pilot_v01 --capture path/to/capture --output path/to/new_import
```

An optional `--benchmark` selects a verified benchmark directory. The default is
the repository benchmark. The importer verifies the original scientific source
chain and benchmark, export hashes, deterministic pilot sampling, experiment
metadata, every task/attempt and source-response bytes. A new output directory
must be separate from the capture/export/benchmark directories. Existing paths
receive overwrite protection. Staging, output verification and same-filesystem
rename publish a complete package. Failures clean up staged files. Capture,
export, benchmark and source bytes are rechecked before publication.

The derived package contains these files.

- `experiment.snapshot.json` preserves the validated experiment metadata bytes.
- `records.snapshot.json` preserves the entire supplied attempt ledger.
- `prompt_export_manifest.snapshot.json` and `prompts.snapshot.json` preserve the original selection and question evidence.
- `capture_validation.json` contains deterministic integrity, coverage and output-format diagnostics, plus selected attempt numbers.
- `submission.json` uses the existing versioned scorer envelope with unchanged submission-kind labeling and primary raw strings.
- `evidence/` preserves every referenced prompt, answer and native payload byte sequence, including secondary attempts.
- `manifest.json` records versions, experiment/selection identity, benchmark hashes, source and output hashes, fixed selected IDs/count, and provenance limitations.

Blocking capture errors produce a deterministic JSON report on standard error,
exit unsuccessfully, and publish no derived package. Library callers receive
the same report through `validate_capture`. Missing tasks and format-invalid
answers can still produce an integrity-valid capture package. The manifest hashes
payloads without hashing itself. Execution timestamps are preserved as supplied.
Conversion generates no additional timestamps and identical captured bytes
produce identical output bytes.

Public library interfaces include `load_verified_prompt_export`,
`validate_experiment_metadata`, `validate_task_record`, `verify_prompt_identity`,
`verify_response_bytes`, `validate_capture` and `convert_to_scoring_submission`.
The source capture root remains the reference for safe path validation.

## Pilot population and Step 21C integration

The pilot remains exactly the fourteen development IDs in its verified export.
The importer preserves selected IDs, selection hash, export manifest, stage,
partition and a `selected_task_count` of fourteen. `submission.json` retains
exactly the existing schema, which has no population or metadata extension fields.
Consequently Step 21C must pass the verified `selected_task_ids` to the existing
`evaluate_responses(..., task_ids=selected_ids)` API for fixed-denominator pilot
scoring. The current `--split development` scoring CLI uses eighty tasks and
counts any absent pilot-external answers as missing. Its coverage-adjusted result
must be described with that eighty-task denominator. Capture validation produces
format diagnostics and provenance, while Step 21C will produce scientific grades.

Both splits and all gold answers are public. Correlated positions/codons, two
transcripts, artificial category balance, uniform prompts and constant phenotype
labels constrain inference. Heuristic score weights retain the documented 90/80
coordinate-error partial credit. Captured model identity remains an attestation
unless independently authenticated. A small pilot supports exploratory procedure
and error analysis. Broad genomic competence, molecular function, drought
tolerance and yield require additional evidence.

## Local Step 21B validation

Python 3.14.7 passed all 592 tests, including the 525-test baseline and 67 new
capture tests. All 67 capture tests also passed with socket networking blocked.
Python 3.12 syntax checks passed. Runtime CI for the new capture implementation
remains pending. Windows restricted native symlink creation during testing, so
the escape regression used a simulated resolved target to test containment
without requesting filesystem privileges.

The ignored `tmp/bioeval_capture_synthetic_v01/` review directory contains a
source package and separately derived import. Its metadata and all fictional
observations are labeled synthetic capture-system testing. Thirteen attempt
records yield ten primary submission rows across fourteen selected tasks,
three missing task records and one recorded task without an eligible response.
One primary row is an explicit abstention. Across all attempts, three answer
texts are format-valid and seven are format-invalid, including secondary
attempts. These are capture-system validation counts.

Tests cover byte preservation, malformed and biologically wrong responses,
refusal/abstention separation, retry selection, prompt/export binding, unsafe
paths, duplicate IDs, corrupted files, output checksums, atomic-style failure
cleanup, source-change detection and deterministic replay. Input files are
loaded into memory. Filesystem byte rechecks provide change detection while
concurrent writers remain outside a process-isolation guarantee. Origin and
arbitrary payload extraction relationships remain researcher attestations.

## Step 21C-1 verification extension

Independent audit reproduced two weaknesses in standalone import verification.
Updating a submission's recorded hash could conceal a change to derived answer
text. Updating the validation-report hash could conceal an invented attempt
count. Verification now replays the preserved metadata, complete ledger and
raw evidence and compares both derived files to the reconstructed results.
`load_verified_capture_import` exposes this replay to the pilot pipeline.
`load_verified_prompt_export(..., imported=True)` and
`validate_capture(..., imported=True)` read the existing snapshot layout.
Schema and output filenames remain compatible with Step 21B packages.


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
