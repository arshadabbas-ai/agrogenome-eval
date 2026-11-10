
# AgroGenome-Eval

### An Evidence-Aware Plant Genomics & Scientific AI Evaluation Platform

**Status:** Active Research & Development  
**Version:** 0.1.0  
**Primary Crops:** Maize (*Zea mays*) and Wheat (*Triticum aestivum*)

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

The first implemented software component is a
Python-based genetic variant classifier.

### Supported Features

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

- 11 automated test cases
- 11 tests passed locally on Windows
- Tests cover basic classification, invalid inputs,
  and lowercase sequence handling

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

## Planned Research Modules

### 1. GenomeLab
DNA sequence analysis, gene structure, translation,
and variant interpretation.

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
