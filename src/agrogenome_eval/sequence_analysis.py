"""Reproducible analysis of unambiguous DNA; no species or phenotype inference.

Raw sequences must contain only A, C, G, T (case insensitive), without
whitespace. FASTA sequence lines are concatenated, with blank lines ignored.
Complement strings retain input order; reverse complements reverse that order.
"""

from os import PathLike
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction


def validate_dna(sequence: str) -> str:
    """Return uppercase DNA or raise TypeError/ValueError for invalid input.

    Ambiguity codes, RNA U, gaps, and whitespace are unsupported.
    """
    if not isinstance(sequence, str):
        raise TypeError("DNA sequence must be a string.")
    if not sequence:
        raise ValueError("DNA sequence must not be empty.")
    if any(base not in "ACGTacgt" for base in sequence):
        raise ValueError("Unsupported DNA characters; only A, C, G, T are allowed.")
    return sequence.upper()


def count_nucleotides(sequence: str) -> dict[str, int]:
    """Count each canonical DNA nucleotide, including zero counts."""
    dna = validate_dna(sequence)
    return {base: dna.count(base) for base in "ATGC"}


def sequence_length(sequence: str) -> int:
    """Return validated DNA length in nucleotides."""
    return len(validate_dna(sequence))


def gc_percentage(sequence: str) -> float:
    """Return 100 * (G + C) / length, without display rounding."""
    return 100.0 * gc_fraction(validate_dna(sequence))


def complement(sequence: str) -> str:
    """Return the complement (3' to 5' for input written 5' to 3')."""
    return str(Seq(validate_dna(sequence)).complement())


def reverse_complement(sequence: str) -> str:
    """Return the uppercase reverse complement in the 5'-to-3' convention."""
    return str(Seq(validate_dna(sequence)).reverse_complement())


def analyze_sequence(sequence: str) -> dict:
    """Return JSON-serializable sequence facts, without biological interpretation."""
    dna = validate_dna(sequence)
    return {
        "sequence": dna,
        "length": sequence_length(dna),
        "nucleotide_counts": count_nucleotides(dna),
        "gc_percentage": gc_percentage(dna),
        "complement": complement(dna),
        "reverse_complement": reverse_complement(dna),
    }


def analyze_fasta(path: str | PathLike[str]) -> dict:
    """Analyze one UTF-8 FASTA record and retain its ID and full description.

    Missing paths raise FileNotFoundError. Invalid encoding, missing headers,
    empty records, multiple records, or invalid DNA raise ValueError. Sequence
    whitespace within a nonempty line is rejected rather than silently removed.
    Headers must start with '>' in column one, immediately followed by an
    identifier. Description spacing is preserved. LF, CRLF, and CR line endings
    and whitespace-only blank lines are accepted; comments are unsupported.
    Other file access errors propagate as OSError subclasses. The complete
    file is held in memory; this reader is intended for small single records.
    """
    fasta_path = Path(path)
    try:
        text = fasta_path.read_text(encoding="utf-8")
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}") from exc
    except UnicodeDecodeError as exc:
        raise ValueError("FASTA file must be valid UTF-8 text.") from exc

    # read_text normalizes CRLF/CR to LF. Do not interpret sequence control
    # characters (e.g. vertical tab) as line boundaries or silently discard them.
    lines = [line for line in text.split("\n") if line.strip(" \t")]
    if not lines or not lines[0].startswith(">"):
        raise ValueError("Malformed FASTA: expected a header starting with '>'.")
    if sum(line.startswith(">") for line in lines) != 1:
        raise ValueError("FASTA must contain exactly one record.")
    if not lines[0][1:].strip():
        raise ValueError("Malformed FASTA: record identifier is missing.")
    description = lines[0][1:]
    if description[0].isspace() or any(
        (ord(char) < 32 and char != "\t") or ord(char) == 127
        for char in description
    ) or ">" in description.split()[0]:
        raise ValueError("Malformed FASTA: invalid record header.")
    sequence = "".join(lines[1:])
    return {
        "record_id": description.split()[0],
        "description": description,
        **analyze_sequence(sequence),
    }
