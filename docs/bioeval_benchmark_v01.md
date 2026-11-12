# BioEval Benchmark Dataset v0.1

BioEval v0.1 is a controlled, public collection of plant molecular sequence
reasoning tasks for researchers developing or evaluating scientific assistants.
It contains no model outputs, model API integration, automated LLM judge,
evaluation scores, or web interface. Package version remains 0.1.0.

## Scientific scope and evidence levels

1. **Reference facts:** versioned transcript sequences, CDS boundaries, reference
   bases, codons, RefSeq curation status, and linked protein accessions from
   locally saved NCBI evidence.
2. **Supported molecular predictions:** synonymous, missense, and stop-gained
   substitutions in internal coding codons, delegated to the existing validated
   classifier through the transcript-aware engine.
3. **Not assessed or unsupported:** UTR protein consequences are not assessed;
   the entire first coding codon and terminal stop codon are unsupported by the
   classifier. Unsupported does not mean biologically harmless. No start-lost
   or stop-lost labels are invented.
4. **Not established:** phenotype, drought tolerance, yield, expression,
   splicing, RNA stability, experimentally established protein function,
   functional equivalence, and one-to-one orthology cannot be inferred here.

All introduced substitutions are computational examples, not observed natural
mutations or verified crop alleles. Coordinates are 1-based inclusive transcript
coordinates, not chromosome/genomic coordinates. CDS starts are not predicted.

| Source | Transcript | CDS | Protein | RefSeq curation |
|---|---|---|---|---|
| Maize (*Zea mays*) | NM_001305944.1 | 156–1259 | NP_001292873.1 | VALIDATED |
| Wheat (*Triticum aestivum*) | XM_044542390.1 | 96–734 | XP_044398325.1 | MODEL |

Sequence-analysis validation status and NCBI RefSeq curation status have different
meanings. The maize record title/CDS-product naming difference and synonyms remain
in the original provenance; the benchmark does not silently rename the gene.

## Dataset and partitions

There are 120 tasks: 60 per transcript/species. Counts below apply **per species**.

| Category | Total | Development | Evaluation |
|---|---:|---:|---:|
| Synonymous | 12 | 8 | 4 |
| Missense | 12 | 8 | 4 |
| Stop-gained | 12 | 8 | 4 |
| 5′ UTR | 6 | 4 | 2 |
| 3′ UTR | 6 | 4 | 2 |
| Start codon, unsupported | 6 | 4 | 2 |
| Terminal stop codon, unsupported | 6 | 4 | 2 |
| Total | 60 | 40 | 20 |

Both species together yield **80 development and 40 evaluation tasks**. IDs and
substitution tuples are unique and disjoint across partitions. All 14 previously
published demonstration substitutions are excluded. Eligible counts after this
exclusion are recorded in the manifest; any shortage blocks generation with an
exact eligible/required/shortage report, rather than duplicating substitutions.

**Both partitions and all answers are public.** The evaluation partition is not
secret, hidden, or contamination-resistant. A genuinely private holdout needs
separate design and storage in a future milestone.

## Files and model-facing task contract

`data/benchmarks/bioeval_v01/` contains:

- `tasks.json`: schema version and tasks only; no gold labels or answer paths.
- `gold_answers.json`: schema version and a task-ID-to-answer mapping.
- `splits.json`: schema version and development/evaluation task-ID lists.
- `manifest.json`: generator version, seed, source identity and checksums,
  eligibility/count reports, generation rules, public disclosure, limitations,
  and output hashes for the other three files. The manifest cannot self-hash.

Every task has exactly `task_id`, `input`, and `question`. IDs are opaque hashes
of substitution tuples; file order is sorted by ID, not scientific category.
The question is constructed from input facts independently of answer serialization.

The strict `input` object contains:

`transcript_accession`, `scientific_name`, `protein_accession`, `refseq_status`,
`transcript_length`, `position`, `reference`, `alternative`, `cds_start`, `cds_end`,
`coordinate_convention`, `coding_context`, and `variant_origin`.

`coding_context` is null outside the CDS. Inside it, the object supplies only
`reference_codon` and `within_codon_offset_1based`. It does not supply the mutated
codon, amino acids, CDS-relative position, codon number, region label, analysis
status, or consequence. Providing a reference triplet is an explicit task-design
choice: this tests molecular interpretation, not sequence retrieval from NCBI.

Questions state the general allowed classifications, unsupported boundaries,
table 1/frame 1, coordinate convention, null rules, and limits of biological
inference. These general rules are identical across categories. No question
contains its correct category, expected result, gold record, or answer-file path.
The allowed label list itself is not a gold answer.

## Strict JSON response and gold-answer contract

Return exactly one JSON object with **all 13 fields below**, no additional fields,
no Markdown fences or explanation. These are reference-answer constraints for
the supplied evidence, not an implemented scoring system for future AI responses.

| Field | Type and valid gold values |
|---|---|
| `task_id` | Nonempty string, exactly matching the task |
| `transcript_accession` | Exact versioned transcript accession from the task |
| `transcript_region` | `FIVE_PRIME_UTR`, `CDS`, or `THREE_PRIME_UTR` |
| `cds_position` | Positive integer in CDS; null in UTR |
| `codon_number` | Positive 1-based integer in CDS; null in UTR |
| `reference_codon` | Uppercase A/C/G/T triplet for assessed internal CDS; otherwise null |
| `alternative_codon` | Uppercase A/C/G/T triplet for assessed internal CDS; otherwise null |
| `reference_amino_acid` | Canonical one-letter amino acid for assessed internal CDS; otherwise null |
| `alternative_amino_acid` | Canonical one-letter amino acid or `*` for assessed internal CDS; otherwise null |
| `molecular_consequence` | `SYNONYMOUS`, `MISSENSE`, or `STOP_GAINED` only when assessed; otherwise null |
| `analysis_status` | `ASSESSED`, `NOT_ASSESSED`, or `UNSUPPORTED` |
| `unsupported_region` | `START_CODON` or `TERMINAL_STOP_CODON` for unsupported boundaries; otherwise null |
| `phenotype_claim_support` | `NOT_ESTABLISHED` for every gold answer in this dataset |

UTR answers use `NOT_ASSESSED` and null for every coding/consequence attribute.
Boundary answers use `UNSUPPORTED`, retain CDS position/codon number, and use null
for codons, amino acids, and molecular consequence, even though a reference
triplet is supplied in the question. Internal CDS answers use `ASSESSED`.
Booleans are not accepted as integers. JSON null is distinct from zero, an empty
string, or the literal string `"null"`. Contract validation and scientific
verification are provided by `validate_response` and `validate_benchmark`.

### Step 20B scoring preparation (not implemented here)

`validate_response` validates **gold/reference answers**, including scientific
constraints. It is not a general model-response parser. For example,
`ESTABLISHED` is an explicitly offered model choice in the prompt, but is a
scientifically wrong phenotype-support answer for these tasks, not necessarily
malformed JSON. A future scorer must separate parse/schema failures from
well-formed but incorrect scientific answers rather than treating every rejected
gold constraint as a syntax error.

Before implementing scores, specify raw JSON handling for missing responses,
duplicate keys, extra fields, non-finite numbers, Markdown wrappers, unknown
labels, booleans used as integers, null versus missing/string-null values, and
amino-acid notation. Specify separately how wrong identities, coordinates,
codons, amino acids, consequences, statuses and phenotype claims enter the error
taxonomy. The structured fields support deterministic comparisons without an
LLM judge, but no scoring, repair, coercion, or evaluation policy is implemented
in this milestone. Report malformed outputs and scientific errors separately.

## Generation, gold verification, and reproduction

Run from the repository root using the installed editable package:

```text
python scripts/generate_bioeval_benchmark.py
python scripts/generate_bioeval_benchmark.py --output path/to/new_directory --seed 20261112
python -m pytest -v
```

Schema version is `0.1`, generator version is `0.1.0`, and the default fixed seed
is **20261112**. Generator version is independent of the package version.

The runner reuses the existing offline accession/organism/GeneID/protein evidence
chain, reference checksums, XML parser, CDS validation, and transcript-aware
classifier. No new biological classification algorithm generates gold labels.
Substitutions are enumerated in ascending transcript position and A/C/G/T
alternative order before sampling. A local `random.Random(seed)` instance
shuffles sorted codon/position groups; round-robin selection maximizes group
diversity before reuse. Every allele is checked against its original sequence.
Each species/category is independently split 2:1 development/evaluation.

All 120 gold answers are independently checked. For all 72 assessed coding
tasks, direct coordinate arithmetic, reference-base checks, single-base codon
substitution, and a separately encoded standard genetic-code lookup verify the
full answer. The audit table follows [NCBI translation table 1](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1).
This independent audit gates output; disagreements are blocking errors, never
automatically corrected labels. Boundary codons do not receive predictions from
an ordinary internal-codon lookup. UTR and boundary statuses are checked separately.

Output JSON uses sorted keys, UTF-8, LF endings, and no generation timestamp.
Identical validated input bytes, seed, and generator version reproduce identical
bytes on supported runtimes. Python 3.14 local validation is performed in this
milestone; Python 3.12 syntax compatibility is checked, but runtime CI is pending
publication. Future changes to algorithms or sampling rules require versioned
regeneration and review, not silent replacement.

An existing output directory is never overwritten. Generation stages and validates
all files before a same-filesystem directory rename. Failed validation or writes
leave no published partial benchmark. SHA-256 values cover exact saved bytes;
Git attributes disable text conversion for benchmark JSON. Source hashes detect
byte changes but are not standalone proof of NCBI identity; versioned accessions,
metadata, source URLs, and sequence comparisons provide that evidence chain.
NCBI attribution and third-party-rights caveats are retained in the source
reference-manifest snapshot; no unconditional public-domain claim is made.

## Interpretation and limitations

This collection tests coordinate mapping, region assignment, codon substitution,
translation, consequence classification, respect for unsupported cases, and
limits of biological inference. It is not a validated measure of broad biological
reasoning or generalization across plant genomes. Accession recognition alone is
insufficient: each question requires interpreting the supplied molecular context.

Tasks share two reference genes and can share positions/codons across partitions.
The start and terminal-stop categories necessarily reuse their respective codons;
even internal tasks may be correlated. The 120 tasks are not 120 independent
biological experiments. Class proportions are deliberately balanced rather than
representing natural variant prevalence. Phenotype support is always not
established, so that field alone is not a balanced phenotype-reasoning benchmark.

Public questions, gold answers, and reproducible sampling permit memorization,
answer leakage, and training-data contamination. No empirical model evaluation,
LLM judge, hidden holdout, statistical validation, or employer-facing performance
claim is provided. Future evaluation must report per-species and per-category
results with uncertainty and task-dependence caveats, not only overall accuracy.

## Independent QA findings for the frozen v0.1 sample

Direct calculations from the original FASTA and saved NCBI XML agree with all
120 gold answers. All 72 assessed reference amino acids also agree with their
version-specific reference proteins. The audit found no wrong scientific labels,
answer-bearing IDs, leaked gold metadata, or substitution overlap in the actual
dataset. Validators additionally enforce tuple-hashed opaque IDs, source metadata,
and the public-answer disclosure to reject those defects in altered datasets.

Measured cross-partition dependence for seed 20261112:

| Source | Shared coding codons (1-based) | Shared transcript positions | Evaluation tasks within 3 nt of a development task |
|---|---|---|---:|
| Maize | 1, 368 | 157, 158, 1257, 1258 | 4 of 20 |
| Wheat | 1, 32, 213 | 97, 733, 734 | 6 of 20 |

The proximity metric includes distance zero. Shared positions use different
alternative alleles; no identical substitution tuple occurs across partitions.
Wheat codon 32 is an internal cross-partition overlap; the other shared codons
are excluded boundary codons. These partitions provide neither gene-level nor
species-level independent generalization.

Difficulty is deliberately constrained. Supplied codon context reduces sequence
retrieval and offset difficulty; its absence is a cue that coding interpretation
does not apply, although the coordinate still determines which UTR. Explicit
classifier rules make unsupported-case recognition rule following rather than
discovery. Every phenotype-support gold value is NOT_ESTABLISHED. A uniform
template and public answers permit shortcut learning and memorization. Necessary
molecular context has not been removed in an attempt to hide those cues.

Recommendations for later, separately reviewed benchmark versions include more
genes/species, gene/codon-grouped splits, independent biological evidence for a
balanced claim-support task, and distinct coordinate-reasoning challenge sets.
Do not silently change this approved 120-task distribution to implement those
design changes. A private holdout also requires separate infrastructure and data.
