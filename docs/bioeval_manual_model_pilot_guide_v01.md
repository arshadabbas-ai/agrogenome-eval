# Manual BioEval model pilot guide v0.1

This guide describes researcher execution after the local preflight is reviewed.
The implementation uses local files and requires no purchased API credits.
Keep future genuine capture/import/scoring records in a local research-output
directory until privacy, provenance and public-release suitability are reviewed.
The ignored repository `tmp/` directory can hold local working records.

## 1. Register the real experiment

Choose a model according to the research comparison and available access.
Record the provider/source, exact visible name and any available model version
or checkpoint. Mark hidden details unavailable with a reason. Record the actual
access method. Capture the experiment start time when execution actually begins,
including timezone. Register the original fourteen task IDs, export/selection
hashes, settings, tool policy, retry bound and first-completed-eligible policy.

Use a fresh, isolated conversation or session for each question where possible.
Avoid prior benchmark-answer context. Record persistent memory, carried context,
custom instructions, browsing, retrieval and available tools. Disable exposed
tools for a tool-free experiment, recording the actual configuration. Hidden
provider instructions and unavailable controls need explicit availability notes.
An isolated session cannot guarantee absence of training-data exposure.

Record actual temperature and seed only when exposed or explicitly configured.
For unavailable parameters use null plus `availability="unavailable"` and an
honest reason. Preserve numeric zero when it was actually set. Record system and
developer instructions as recorded text, explicitly absent, or unavailable.

## 2. Bind the preserved question export

Use the existing reviewed `tmp/bioeval_pilot_v01` export. It contains exactly
fourteen development tasks, seven per species, with seed 20261112. Retain its
manifest and hashes. The following local inspection retrieves real identities
and bytes without creating responses.

```python
from pathlib import Path
from agrogenome_eval.bioeval_capture import (
    canonical_json, load_verified_prompt_export, sha256,
)

prompts = load_verified_prompt_export(
    Path("tmp/bioeval_pilot_v01"), Path("data/benchmarks/bioeval_v01")
)
print(prompts.identity)
task = prompts.records[sorted(prompts.records)[0]]
question_bytes = task["question"].encode("utf-8")
print(task["task_id"])
print(task["question"])
print(sha256(question_bytes))
print(sha256(canonical_json(task)))
```

Send each complete original question unchanged as the user message. Keep gold
answers, internal strata, manifests, synthetic fixtures and scientific grades
outside the model session. The question supplies response fields and scientific rules, with schema
disclosure limitations described below. Record all additional known messages/context.
This baseline importer requires exact exported user-question bytes. Additional
user wrappers or manual prompt edits require a separately reviewed protocol.

## 3. Preserve observations and execution history

Create a source capture directory with `experiment.json`, `records.json` and
relative prompt/response files according to the capture specification. Use the
actual observed details, never a synthetic template's fictional identity/date.
For each task, preserve the exact UTF-8 question bytes and calculate its record
and question hashes. Save copied output in a UTF-8 file without trimming,
rewriting, fence removal, case changes or JSON repair.

UI copying may transform the source text. Record the copying/export method and
its limitations honestly. Record the actual task timestamp with UTC offset,
attempt number, completion/refusal/truncation outcome, output errors and tool
use. Preserve refusal text. Save an empty observed response as a zero-byte file.
A technical failure or no-response event has null answer path/hash and an
explicit description. A documented abstention has its distinct outcome and
null answer path/hash.

Calculate hashes from saved bytes, for example with this local operation.

```python
from pathlib import Path
from agrogenome_eval.bioeval_capture import sha256

response_path = Path("tmp/my_real_pilot/source/responses/task_response.txt")
print(sha256(response_path.read_bytes()))
```

Use safe relative file paths in the ledger. Preserve native API payloads and
document extraction separately if a future run uses API access. Exclude
credentials, authorization headers and unnecessary private conversation data.

The source ledger contains every attempt. Use contiguous attempt numbers,
chronological timestamps, previous-attempt links and explicit retry reasons.
After an execution failure before any answer, a later eligible completion can
be primary. Once a completed eligible output exists, retain it as primary even
when malformed or scientifically wrong. Additional attempts remain secondary.
Selecting a later completion by inspecting gold answers or scores invalidates
the declared attempt policy. Preserve original evidence when documenting any
capture correction and register a separate run for a changed protocol.

## 4. Validate and import the capture

Choose fresh derived output paths. The source JSON metadata/record fields are
fully documented in [response capture v0.1](bioeval_response_capture_v01.md).

```text
python scripts/import_bioeval_model_responses.py --prompt-export tmp/bioeval_pilot_v01 --capture tmp/my_real_pilot/source --output tmp/my_real_pilot/import
```

Review `capture_validation.json`. Resolve integrity errors through documented
evidence correction. Formatting errors in genuine output stay preserved and
receive the declared scoring treatment. Review missing IDs, failed executions,
unavailable provenance and selected attempt numbers. Confirm `model_responses`
only for actual independently executed observations. The importer preserves
synthetic labeling for test packages and refuses incompatible provenance.

Researcher-entered identities and attestations remain unverified origin claims.
Checksums establish integrity relative to saved values. They do not authenticate
the model provider, establish UI byte identity or validate arbitrary extraction.

## 5. Score the selected fourteen questions

```text
python scripts/score_bioeval_pilot.py --prompt-export tmp/bioeval_pilot_v01 --capture-import tmp/my_real_pilot/import --output tmp/my_real_pilot/scored
```

The pilot scorer replays the imported raw evidence, checks original questions
and invokes the existing scorer on exactly fourteen selected IDs. Inspect
`report.json`, snapshots and `manifest.json`. Confirm the denominator is fourteen
and all task outcomes sum to fourteen. Keep captured-task, primary-response,
attempt and valid-structured-answer counts distinct. Report refusals, truncations,
failures and secondary attempts alongside scientific grades.

The existing general scoring CLI's `--split development` uses eighty tasks.
This pilot-specific command is required for the declared fourteen-task analysis.
Report exact accuracy, valid-answer mean, coverage-adjusted score, coverage,
format compliance, species/category/status results and task diagnostics. A
coordinate error can retain 90 or 80 partial points under the heuristic policy.

## 6. Report the experiment with its limits

Include model identity visibility, execution times, access conditions, settings
availability, full known prompts, prior benchmark exposure, attempt policy,
missing/failure counts, source/output hashes and software/policy versions.
Report descriptive counts before interpreting percentages. Each species/category
cell has one task. Public answers, correlated genes/codons, two transcripts,
uniform prompts and constant phenotype labels constrain interpretation.
Avoid confidence intervals that assume fourteen independent biological samples.
Observed performance supports controlled structured-task analysis under the
recorded conditions. Broad genomic competence, gene function, drought tolerance
and crop performance require independent evidence.


## Completed pilot and schema disclosure

The [public pilot summary](bioeval_pilot_results_v01.md) describes the completed
fourteen-task experiment and its audit. The original questions omit exact UTR
enum spellings and the one-letter amino-acid requirement, although the published
scorer requires those conventions. Preserve this distinction when interpreting
schema failures. A future experiment can supply a complete schema under a new
protocol version. Historical captured responses retain their original bytes.

The summary releases aggregates and methodological findings. Detailed source
captures, execution ledgers, researcher attestations and imported/scored metadata
snapshots remain private pending a separate content and publication decision.
