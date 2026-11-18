# BioEval v0.1 empirical pilot results and release notes

## Scope and release status

This de-identified summary accompanies the proposed Milestone 07 source release.
It describes a completed, audited fourteen-task development-set pilot using the
researcher-attested OpenAI GPT-6 High identity through the manual ChatGPT
interface. Model identity lacks independent authentication. Publication approval
and Linux Python 3.12 validation of this release candidate remain pending.

BioEval v0.1 contains 120 benchmark tasks generated from two versioned reference
transcripts, maize NM_001305944.1 and wheat XM_044542390.1. There are 80 development
and 40 evaluation tasks. The pilot selected seven tasks per species, one per
species/category cell, with sampling seed 20261112. All fourteen selected tasks
have one recorded completed attempt. The full benchmark and its gold answers
are public. Task dependence and possible prior answer exposure constrain
interpretation. [Benchmark specification](bioeval_benchmark_v01.md)

## Results and denominators

N denotes fourteen selected tasks, V the two schema-valid answers, E the two
exact answers, and P the 200-point sum. Scores range from 0 to 100 per task.

| Metric | Count or calculation | Result |
|---|---|---:|
| Selected tasks | N | 14 |
| Recorded attempts | One completed attempt per task | 14 |
| Syntactically valid JSON objects | Unchanged captured text parses as an object | 5 |
| Schema-valid responses | V | 2 |
| Invalid structured responses | Nine JSON syntax failures and three enum failures | 12 |
| JSON syntax failures | INVALID_JSON | 9 |
| Enumeration failures | INVALID_ENUM_VALUE | 3 |
| Strict exact matches | E | 2 |
| Strict exact-match accuracy | E/N = 2/14 | 0.14285714285714285 |
| Exact-match accuracy among valid answers | E/V = 2/2 | 1.0 |
| Final score sum | P | 200 |
| Mean partial credit among valid answers | P/V = 200/2 | 100.0/100 |
| Coverage-adjusted mean score | P/N = 200/14 | 14.285714285714286/100 |
| Response coverage | 14/14 | 1.0 |
| Schema compliance among submitted tasks | 2/14 | 0.14285714285714285 |
| Refusals | Recorded execution outcomes | 0 |
| Explicit abstentions | Recorded execution outcomes | 0 |
| Missing outcomes or responses | Selected-task denominator | 0 |
| Retries, truncations or execution failures | Recorded execution outcomes | 0 |

The twelve invalid answers received zero before biological grading. Their
biological correctness has not been comprehensively measured. Five captured
responses are syntactically valid JSON. Three of those fail the response schema,
leaving two eligible for scientific grading. Both eligible answers match all
thirteen gold fields and score 100. The mean of 100 uses V=2, while the
coverage-adjusted mean uses N=14. Strict exact accuracy is approximately 14.2857%.
[Published scoring policy](bioeval_scoring_v01.md)

| Species | Selected | Schema-valid | Invalid | Exact | Strict exact accuracy | Coverage-adjusted mean /100 |
|---|---:|---:|---:|---:|---:|---:|
| Maize, Zea mays | 7 | 1 | 6 | 1 | 1/7 = 0.14285714285714285 | 14.285714285714286 |
| Wheat, Triticum aestivum | 7 | 1 | 6 | 1 | 1/7 = 0.14285714285714285 | 14.285714285714286 |

| Gold category | Tasks | Exact | Invalid | Strict exact accuracy | Coverage-adjusted mean /100 |
|---|---:|---:|---:|---:|---:|
| FIVE_PRIME_UTR | 2 | 0 | 2 | 0.0 | 0.0 |
| MISSENSE | 2 | 0 | 2 | 0.0 | 0.0 |
| START_CODON | 2 | 1 | 1 | 0.5 | 50.0 |
| STOP_GAINED | 2 | 1 | 1 | 0.5 | 50.0 |
| SYNONYMOUS | 2 | 0 | 2 | 0.0 | 0.0 |
| TERMINAL_STOP_CODON | 2 | 0 | 2 | 0.0 | 0.0 |
| THREE_PRIME_UTR | 2 | 0 | 2 | 0.0 | 0.0 |

| Gold analysis status | Tasks | Exact | Invalid | Strict exact accuracy | Coverage-adjusted mean /100 |
|---|---:|---:|---:|---:|---:|
| ASSESSED | 6 | 1 | 5 | 0.16666666666666666 | 16.666666666666668 |
| NOT_ASSESSED | 4 | 0 | 4 | 0.0 | 0.0 |
| UNSUPPORTED | 4 | 1 | 3 | 0.25 | 25.0 |

These group percentages describe sparse pilot cells. The same fourteen records
supply the species, category and analysis-status views.

## Scientific interpretation and audit findings

One eligible answer correctly classifies an internal maize stop-gained
substitution, and the other correctly identifies an excluded wheat start-codon
case. The scorer emits no field-level scientific mismatch for either eligible
answer. Scientific field-error denominators are two. Invalid answers have null
field matches and receive no inferred coordinate, codon, amino-acid, consequence,
or phenotype-error judgments.

Seven captured responses contain Markdown fences. The published parser rejects
those wrappers, consistent with its pre-existing strict policy. Two further
responses contain literal backslashes and HTML space entities. The origin of
those copying-transformation artifacts remains unresolved between the provider
display and the capture/transfer process. Source checksums preserve the supplied
bytes but cannot establish their equivalence to the provider's original output.

The remaining three invalid responses use UTR spellings outside the accepted
schema enums. Certain output conventions were incompletely specified in the
model-facing questions, including exact UTR enum spellings and required
one-letter amino-acid notation. This limits attribution of schema rejection to
instruction-following or biological reasoning. The original questions and scorer
retain their established definitions. A future protocol can disclose the full
schema prospectively. [Empirical protocol](bioeval_empirical_protocol_v01.md)

The evidence supports low end-to-end strict structured-response success on the
captured pilot. It leaves biological accuracy across twelve rejected responses
incompletely assessed. Two schema-valid successes provide limited evidence about
specific supported tasks. The benchmark is insufficient to establish general AI
performance in plant genomics. Molecular consequences, expression, splicing,
protein function, drought tolerance and yield require their own evidence.
Constant NOT_ESTABLISHED phenotype gold labels also limit evaluation of
positive-versus-negative claim discrimination.

## Provenance and reproducibility boundaries

The researcher supplied a final attestation and execution records. Exact model
checkpoint, temperature, provider generation seed, hidden instructions and
effective tool permissions are unavailable. Memory and personalization conditions
are unknown. A recorded absence of observed tool use does not establish disabled
tools. Identity remains researcher-attested. The task-sampling seed describes
selection and supplies no provider-generation seed information.

The local audit verified all fourteen task mappings, first attempts, prompt and
response hashes, transaction records, imported evidence and scoring snapshots.
It independently recalculated all overall and grouped aggregates and field-error
denominators. Fresh import and scoring replays reproduced every saved output
byte for byte. Capture validation passed with six warnings for explicitly
unavailable provenance fields and no ERROR entries.

These audit statements summarize checks on private evidence. The public summary
contains aggregates and methods. Original responses, detailed execution records,
researcher attestations and full import/scoring snapshots remain private.
Public readers can reproduce task selection and software behavior with synthetic
fixtures. Full empirical replay requires authorized access to the original
private package. The public summary alone cannot independently authenticate the
model or reproduce the private response-to-score chain.

## Public offline workflows

After installing the project from its root, the following checks operate on the
public source and saved reference data. The repository tests construct their own
explicitly synthetic capture fixtures in temporary directories.

```text
python -m pip install -e ".[dev]"
python -m pytest -q
python -m pytest -q tests/test_bioeval_empirical_protocol.py tests/test_bioeval_capture.py tests/test_bioeval_pilot.py tests/test_bioeval_scoring.py
python scripts/export_bioeval_model_prompts.py --partition development --pilot --seed 20261112 --output tmp/release_pilot_questions
python scripts/score_bioeval_responses.py --responses data/examples/bioeval_mock_responses_v01.json --benchmark data/benchmarks/bioeval_v01 --split all --output tmp/release_synthetic_scores
```

Choose new output directories when repeating CLI commands. Task export preserves
the original questions, deterministic order and researcher-only manifest. The
manifest records selection identity, seed, source/benchmark hashes, question
hashes and output checksums. Its internal strata remain outside model sessions.
The mock-scoring command exercises an explicitly synthetic fixture and produces
its separately documented expected results. [Prompt-export protocol](bioeval_empirical_protocol_v01.md)

For a separately collected, reviewed capture, the public commands are below.
`private/bioeval_capture` denotes the researcher's own capture directory containing
experiment metadata, ledger and relative evidence paths. The software issues
no model requests. Record independently observed outcomes before importing.

```text
python scripts/import_bioeval_model_responses.py --prompt-export tmp/release_pilot_questions --capture private/bioeval_capture --output tmp/release_capture_import
python scripts/score_bioeval_pilot.py --prompt-export tmp/release_pilot_questions --capture-import tmp/release_capture_import --output tmp/release_pilot_scores
```

The fourteen-task pilot scorer preserves the selected-task denominator. The
general scorer's development split contains eighty tasks. Output packages include
original response text and metadata, so they need private storage and content
review before any release. [Capture procedure](bioeval_response_capture_v01.md)
[Manual collection guide](bioeval_manual_model_pilot_guide_v01.md)
[Pilot evaluator](bioeval_pilot_evaluation_v01.md)

## Scientific release notes and future work

This release adds offline question export, capture provenance validation,
byte-preserving import, fixed-population pilot scoring, synthetic regression
coverage, methodology documentation and the de-identified pilot summary.
Published reference data, benchmark questions, gold labels, classifier semantics,
and scientific scoring weights/caps retain their established definitions.

Follow-up research should disclose the complete response schema, assess capture
fidelity, expand independent transcript/gene sources, and introduce separately
designed private evaluation material. Semantic evaluation of rejected answers
needs a declared rubric and separate reporting. New experiments should document
model settings, repeated trials, session conditions and attempt budgets before
collection. Mixed positive and negative phenotype evidence would support a more
discriminating claim-support task. Dashboards and future APIs can reuse the
verified scientific modules after these measurement limits are represented in
their interfaces.
