# BioEval scoring v0.1 independent QA

This audit covers local Step 20B outputs on
`feature/bioeval-benchmark-v01`. It is not an empirical AI-model evaluation.
Local `main` and `origin/main` remain
`072f9c6c178266f84f440efa7b9c8b9c8dd186bc`.

## Scoring mathematics and policy findings

The approved policy remains `bioeval-scoring-v0.1`; its weights and caps were
not changed. Each profile totals 100:

- ASSESSED: 10+15+10+10+10+10+5+5+20+5.
- NOT_ASSESSED: 30+30+20+5+15.
- UNSUPPORTED: 15+25+15+10+20+10+5.

For structurally valid answers, final score is the minimum of the nonnegative
component sum and **every** applicable cap. Thus caps cannot increase a score;
all scores remain in [0, 100]. Exact correctness requires all 13 fields,
including identities and expected nulls. Every biological field null is invalid
because region/status/phenotype are required enum strings; null coding details
alone earn only 30 on an otherwise correct ASSESSED answer. Wrong consequence
caps at 50; contradictory amino-acid/codon assertions also gate the consequence
component. UTR and boundary overclassification cap at 20 and do not imply that
unassessed substitutions are biologically harmless.

Important interpretation: a wrong CDS position within the same codon retains
90; an internally coherent wrong CDS-position/codon-number pair retains 80.
These are explicit weighted-partial-credit policy choices, **not** completely
correct molecular answers. Both lose exact-match credit and expose coordinate
error tags. A wrong position inconsistent with the asserted codon number caps
at 50. Readers must not present partial scores as whole-answer accuracy. A
stricter coordinate cap would require an explicitly reviewed policy revision.

Weights and caps are heuristic choices, not empirically calibrated probabilities.
No biological consequence algorithm, gold label, question, or split was changed.

## Independent synthetic demonstration arithmetic

Manual expectations are stored separately from the scorer and never supplied to
it. Checking raw field mismatches and the documented caps gives:

| Deliberate case | Component sum | Final score |
|---|---:|---:|
| Correct coding answer | 100 | 100 |
| CDS position wrong within same codon | 90 | 90 |
| Codon number inconsistent with CDS position | 90 | 50 |
| Wrong consequence contradicting amino-acid change | 80 | 50 |
| Wrong alternative amino acid, consequence gated | 75 | 50 |
| Invented UTR consequence | 80 | 20 |
| Invented boundary consequence | 90 | 20 |
| Unsupported phenotype assertion | 95 | 50 |

The eight valid scores total `100+90+4*50+2*20 = 430`.
Mutually exclusive task counts are `8+3+1+108 = 120`: valid, invalid,
abstained, missing. The duplicated task is one of the three invalid tasks.
Fourteen submitted rows include one extra duplicate row and one unknown-ID row;
only twelve selected task IDs were submitted. Unknown IDs are submission issues,
not additional selected tasks.

| Metric | Numerator / denominator | Value |
|---|---|---:|
| Response coverage | 12 / 120 | 0.1 |
| Scientific-answer coverage | 8 / 120 | 0.0666667 |
| Valid-response rate among submitted selected tasks | 8 / 12 | 0.6666667 |
| Strict accuracy across selected tasks | 1 / 120 | 0.0083333 |
| Strict accuracy among valid answers | 1 / 8 | 0.125 |
| Mean partial score among valid answers | 430 / 8 | 53.75 |
| Coverage-adjusted mean score | 430 / 120 | 3.5833333 |
| Format compliance among submitted selected tasks | (8+1) / 12 | 0.75 |
| Abstention rate across selected tasks | 1 / 120 | 0.0083333 |
| Missing-response rate across selected tasks | 108 / 120 | 0.9 |

All undefined ratios return JSON null. The added rate fields apply identically
to overall and gold-defined species/category/status groups. Field error rates
use valid answers as denominator; task taxonomy codes count once per task,
separately from submission issue records. Direct field mismatches and derived
contradiction tags may be correlated; their sum is not an independent-failure
count. The implementation does not assign causal primary/secondary rankings.

The original `data/examples/derived/bioeval_scoring_v01/` is preserved unchanged.
The audited additive report is in
`data/examples/derived/bioeval_scoring_v01_qa/`. Both remain explicitly labeled
synthetic scorer testing. The original report predates the three added rate
fields; its existing metrics and checksums remain valid.

## Confirmed defects and focused corrections

1. **Missing explicit rate outputs.** Added valid-response rate among submitted
   selected tasks (V/S), abstention rate (A/N), and missing-response rate (M/N),
   with null for zero denominators. Existing metrics retain their meanings.
2. **Benchmark verification/read race.** A regression replaced the benchmark
   manifest immediately after successful validation; the old runner accepted it.
   The runner now snapshots all four files before verification and rejects byte
   changes immediately afterwards, then retains the pre-publication recheck.
   It also snapshots/rechecks scoring-policy bytes and protects existing leaf
   symlinks, including dangling symlinks, from publication attempts.

No incorrect scientific gold labels or incorrect profile arithmetic were found.
The scorer's function names match the requested implementation inventory.

## Validation and integrity

- Full suite: **497 passed, 0 failed**, including all 485 previous tests.
- Scorer suite: 79 tests, including 12 focused audit regression cases.
- Audit coverage includes escaped duplicate keys, empty/multiple JSON roots,
  excessive nesting/integers, cap precedence, coordinate-only partial credit,
  independently grouped demo arithmetic, wholly missing/malformed/abstaining
  submissions, malformed abstentions, task substitution, wrong public split,
  contradictory duplicates, existing-file protection, and the reproduced race.
- Existing tests exercise malformed scientific types, unknown enums, booleans
  as coordinates, nonfinite numbers, corrupted gold/source bytes, deterministic
  replay, and failed-write/rename cleanup.
- Python 3.12 syntax checks pass. Local runtime is Python 3.14.7; Python 3.12
  runtime CI remains pending. Scorer tests run with socket networking blocked.
- Both demonstration checksum manifests verify; current policy and benchmark
  hashes agree with both reports. Replay of the audited report is byte-identical.
- Source-chain verification and seeded benchmark replay pass. Pre-audit exact
  byte fingerprints of all benchmark files, authentic references, previous
  evidence, and existing scientific code remain unchanged.
- Tracked scientific files agree with local main under existing Git text
  normalization. Some pre-existing working copies use CRLF while Git stores LF;
  exact evidence hashes always cover the saved bytes specified by each manifest.
- Manifests hash their payload files without self-hashing; reports record the
  benchmark-manifest hash, so there is no circular hash dependency.
- The unrelated root workflow remains untouched and unstaged. No new dependency,
  package-version change, model call, commit, push, PR, or deployment occurred.

## Remaining limitations and publication assessment

Strict JSON is treated solely as data, without execution or repair. Parser
failures are diagnosed without expanding them into unsupported scientific
judgments. This is not a semantic judge of arbitrary free-text hallucinations.
All input files are loaded in memory; there is no submission-size quota, streaming
parser, process resource isolation, or hostile multi-writer filesystem guarantee.
Repeated byte checks detect intervening changes, but do not provide filesystem
locking or authenticated signatures. The CLI is intended for controlled offline
research, not exposure as a hardened public upload service. File permissions and
OS security controls remain in force; no symlink privileges are requested.

Only two transcripts, correlated/cross-partition codons, artificial category
balance, uniform templates, public gold answers, no private holdout, and constant
NOT_ESTABLISHED phenotype labels constrain validity. Shortcuts and memorization
can produce high scores. No general species/genome competence, crop-performance,
functional prediction, or experimental biological findings follow from them.

The implementation is ready for independent publication review, with the
coordinate partial-credit caveat explicitly retained. Python 3.12 runtime CI is
still required before claiming cross-runtime validation. All work stays local,
uncommitted, and unstaged.
