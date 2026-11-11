
# AgroGenome-Eval

### An Evidence-Aware Plant Genomics & Scientific AI Evaluation Platform

**Status:** Active Research & Development  
**Version:** 0.1.0  
**Research Focus:** Maize (*Zea mays*) and Wheat (*Triticum aestivum*)

## Overview

AgroGenome-Eval is an open research software project
designed to combine computational plant genomics,
genetic variant interpretation, and scientific AI
evaluation.

The long-term goal is to develop a publicly accessible
research platform supporting plant breeders, molecular
geneticists, bioinformaticians, and AI researchers.

The project emphasizes scientific reproducibility,
biological evidence, transparent evaluation, and
responsible interpretation of genomic information.

## Research Objectives

1. Develop reproducible DNA sequence analysis methods.
2. Analyze protein-coding sequences and genetic variants.
3. Support maize and wheat genomics research.
4. Develop benchmarks for evaluating AI-generated
   biological explanations.
5. Identify scientific reasoning errors and unsupported
   biological claims.
6. Develop an interactive public research application.

## Current Implementation — v0.1.0

Two scientific Python components are implemented: the Genetic Variant
Classifier and the GenomeLab DNA Sequence Analysis Engine. These tools
provide sequence-level calculations; crop-specific research workflows,
AI benchmarks, and the public web platform remain planned work.

### Genetic Variant Classifier

- Single-nucleotide substitutions in complete CDS records
- Missense variant classification
- Synonymous variant classification
- Stop-gained variant classification
- CDS position and codon identification
- Reference and alternative nucleotide validation
- Amino acid consequence prediction
- Automated unit testing using Pytest

The classifier currently supports substitutions in
internal codons of validated coding sequences.

### GenomeLab DNA Sequence Analysis Engine

Import these functions from `agrogenome_eval.sequence_analysis`:

| Function | Purpose and return value |
|---|---|
| `validate_dna(sequence)` | Validate a nonempty DNA string and return uppercase DNA. |
| `count_nucleotides(sequence)` | Return a dictionary of A, T, G, and C counts, including zeros. |
| `sequence_length(sequence)` | Return sequence length in nucleotides as an integer. |
| `gc_percentage(sequence)` | Return `100 * (G + C) / length` as an unrounded float. |
| `complement(sequence)` | Return the uppercase complement in input positional order (3′→5′ for input written 5′→3′). |
| `reverse_complement(sequence)` | Return the uppercase reverse complement in the 5′→3′ convention. |
| `analyze_sequence(sequence)` | Return a JSON-serializable dictionary containing the normalized sequence, length, counts, GC percentage, complement, and reverse complement. |
| `analyze_fasta(path)` | Analyze one FASTA record from a string or path-like file path; add `record_id` and `description` to the sequence-analysis dictionary. |

#### DNA and FASTA input contract

- DNA must be a nonempty string containing only A/C/G/T, in either case.
  Raw sequence strings cannot contain whitespace, ambiguity codes, RNA U, or gaps.
- FASTA input must be one UTF-8 record. The header starts with `>` in column
  one, immediately followed by a nonempty identifier. The identifier is the
  first whitespace-delimited token; `description` retains the full header text
  after `>`, including its spacing.
- Sequence lines are concatenated and normalized to uppercase. LF, CRLF, and
  CR line endings and empty or space/tab-only lines are accepted. Whitespace
  inside nonempty sequence lines is rejected. Comments are unsupported.
- Invalid DNA and malformed, empty, multi-record, or invalid-UTF-8 FASTA inputs
  raise `ValueError`. Non-string DNA raises `TypeError`; missing files raise
  `FileNotFoundError`. Other file-access errors propagate as `OSError` subclasses.

The analysis dictionaries use `sequence`, `length`, `nucleotide_counts`,
`gc_percentage`, `complement`, and `reverse_complement` keys. Values are ordinary
Python strings, integers, floats, and dictionaries, suitable for `json.dumps`.

## Scientific Case Study: Maize

An initial analysis was performed using an NCBI
reference transcript associated with the maize
DREB2A / LOC732788 gene.

Reference information:

- Organism: Zea mays
- Gene locus: LOC732788
- Transcript: NM_001305944.1
- Protein: NP_001292873.1
- Transcript length: 1,507 nucleotides
- Annotated CDS length: 1,104 nucleotides
- Translated protein length: 367 amino acids

The annotated protein translation was successfully
reproduced using Biopython.

Three synthetic mutations were evaluated:

| Variant | Codon Change | Consequence |
|---------|--------------|-------------|
| Missense | ACG to CCG | Threonine to Proline |
| Synonymous | ACG to ACA | Threonine unchanged |
| Stop-gained | TTG to TAG | Premature stop codon |

These mutations are computational examples, not
experimentally confirmed maize variants.

## Automated Testing

The project uses Pytest for scientific software
validation.

Current results:

- 108 passing automated tests: 97 sequence-analysis tests and 11 classifier tests
- Verified locally on Windows with Python 3.14.7 and in GitHub Actions on Linux
  with Python 3.12.16
- Coverage includes independent synthetic reference calculations, complement
  orientation, lowercase normalization, invalid inputs, FASTA parsing and
  provenance, and JSON serialization, alongside the existing classifier tests

Further biological edge-case testing is planned
before public deployment.

## Technology Stack

- Python
- Biopython
- Pytest
- Setuptools
- Google Colab
- Visual Studio Code
- Git and GitHub

## Installation

From the local project directory (`AgroGenome-Eval`):

Create and activate a virtual environment.

Install the package and development dependencies:

    python -m pip install -e ".[dev]"

Run the automated tests:

    python -m pytest -v

## Example Usage

### Genetic variant classification

```python
from agrogenome_eval import classify_variant

cds = "ATGACGTTGTGA"

result = classify_variant(
    cds,
    position=4,
    alternative="C"
)

print(result["Classification"])

# Expected output: MISSENSE
```

### DNA sequence analysis and FASTA provenance

This synthetic example has no assigned crop species or gene identity.
The FASTA file is created in a temporary directory and removed automatically.

```python
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from agrogenome_eval.sequence_analysis import analyze_fasta, analyze_sequence

dna = "ATGGCTACCGTTAAGGCTTACGATCGATGA"
analysis = analyze_sequence(dna)

assert analysis["length"] == 30
assert analysis["nucleotide_counts"] == {"A": 8, "T": 8, "G": 8, "C": 6}
print(f"GC content: {analysis['gc_percentage']:.4f}%")
# GC content: 46.6667%

with TemporaryDirectory() as directory:
    fasta_path = Path(directory) / "synthetic.fasta"
    fasta_path.write_text(
        f">synthetic_01 demonstration sequence\n{dna.lower()}\n",
        encoding="utf-8",
    )
    record = analyze_fasta(fasta_path)

assert record["record_id"] == "synthetic_01"
assert record["description"] == "synthetic_01 demonstration sequence"
assert record["sequence"] == dna
print(json.dumps(record, indent=2))
```

## Planned Research Modules

### 1. GenomeLab extensions
The DNA analysis engine above is implemented. Broader gene-structure,
translation, and integrated variant-interpretation workflows remain planned.

### 2. CropStress Explorer
Maize and wheat stress-response research, with an
initial emphasis on drought and seedling-stage biology.

### 3. CrossCrop Bridge
Evidence-based comparative analysis of selected genes
across maize, wheat, and model plants.

### 4. ClaimCheck AI
Evaluation of scientific claims against supporting
biological evidence and documented uncertainty.

### 5. BioEval Studio
A reproducible benchmark for evaluating AI-generated
genomic calculations and biological interpretations.

## Planned Web Platform

AgroGenome-Eval is intended to evolve into a public
scientific web application.

Planned capabilities include:

- Biological sequence input and analysis
- Genetic variant classification
- Crop stress research tools
- Evidence-aware scientific interpretation
- AI evaluation workflows
- Downloadable research reports

The web application is not yet deployed.

## Scientific Limitations

GenomeLab accepts only unambiguous A/C/G/T DNA and single-record UTF-8 FASTA.
It loads the whole file into memory and targets small individual records,
not streaming analysis of large genomes. Ambiguous nucleotides, RNA, gaps,
FASTA comments, and multiple records are unsupported.

Sequence composition does not establish species, gene function, or phenotype.
GenomeLab reports sequence-derived facts and makes none of these inferences.
The maize/wheat research focus does not assign a crop identity to arbitrary input.

The strict FASTA reader is custom. In the local Windows environment, Application
Control blocks a native extension imported by Biopython `SeqIO`. Biopython is
still used for GC and complementary-strand operations; no security controls
are bypassed. The reader supports only the input contract documented above.

Current predictions describe sequence-level molecular
consequences, not experimentally established phenotypes.

A predicted genetic variant effect must not be
interpreted as proof of improved drought tolerance,
yield, or other complex agricultural traits.

The current classifier does not support indels,
transcript-aware genomic coordinates, alternative
splicing consequences, or all possible variant types.

## Project Author

**Arshad Abbas**

Background: Plant Breeding & Genetics, Bioinformatics,
Data Science, and Applied Artificial Intelligence.

GitHub: https://github.com/arshadabbas-ai

## Development Status

AgroGenome-Eval is under active development.

Scientific modules, benchmark datasets, and the public
web application will be released incrementally after
testing and validation.
