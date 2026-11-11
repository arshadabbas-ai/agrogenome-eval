"""Tests use synthetic DNA without assigning crop species or phenotypes."""

import json

import pytest

from agrogenome_eval.sequence_analysis import (
    analyze_fasta,
    analyze_sequence,
    complement,
    count_nucleotides,
    gc_percentage,
    reverse_complement,
    sequence_length,
    validate_dna,
)


def test_independent_synthetic_sequence():
    result = analyze_sequence("ATGGCTACCGTTAAGGCTTACGATCGATGA")
    assert result["length"] == 30
    assert result["nucleotide_counts"] == {"A": 8, "T": 8, "G": 8, "C": 6}
    assert result["gc_percentage"] == pytest.approx(100 * 14 / 30)
    assert result["gc_percentage"] == pytest.approx(46.67, abs=0.005)
    assert result["complement"] == "TACCGATGGCAATTCCGAATGCTAGCTACT"
    assert result["reverse_complement"] == "TCATCGATCGTAAGCCTTAACGGTAGCCAT"
    assert json.loads(json.dumps(result)) == result


@pytest.mark.parametrize("sequence", ["ATGC", "atgc", "aTgC"])
def test_normalization(sequence):
    assert validate_dna(sequence) == "ATGC"
    assert analyze_sequence(sequence) == {
        "sequence": "ATGC", "length": 4,
        "nucleotide_counts": {"A": 1, "T": 1, "G": 1, "C": 1},
        "gc_percentage": 50.0, "complement": "TACG",
        "reverse_complement": "GCAT",
    }


@pytest.mark.parametrize("sequence,expected", [("A", 0), ("ttaa", 0), ("gcgc", 100)])
def test_gc_boundaries_and_zero_counts(sequence, expected):
    assert gc_percentage(sequence) == expected
    assert sequence_length(sequence) == len(sequence)
    assert sum(count_nucleotides(sequence).values()) == len(sequence)
    assert set(count_nucleotides(sequence)) == set("ATGC")


def test_complement_orientation():
    assert complement("aagc") == "TTCG"
    assert reverse_complement("aagc") == "GCTT"
    assert reverse_complement(reverse_complement("AAGC")) == "AAGC"


@pytest.mark.parametrize("function", [validate_dna, count_nucleotides, sequence_length,
                                    gc_percentage, complement, reverse_complement,
                                    analyze_sequence])
@pytest.mark.parametrize("sequence", ["", "N", "ACGU", "AC-G", "AC G", "ACG\n", "123", ">id"])
def test_invalid_dna(function, sequence):
    with pytest.raises(ValueError, match="empty|Unsupported"):
        function(sequence)


@pytest.mark.parametrize("sequence", [None, 123, b"ATGC", ["A", "T"]])
def test_non_string_input(sequence):
    with pytest.raises(TypeError, match="string"):
        validate_dna(sequence)


@pytest.mark.parametrize("newline", ["\n", "\r\n"])
def test_fasta_provenance_and_multiline(tmp_path, newline):
    path = tmp_path / "synthetic.fasta"
    path.write_bytes(newline.join([">synthetic_01 laboratory example", "aag", "c", ""]).encode())
    result = analyze_fasta(path)
    assert result == {"record_id": "synthetic_01", "description": "synthetic_01 laboratory example",
                      **analyze_sequence("AAGC")}
    assert json.loads(json.dumps(result)) == result


def test_fasta_blank_lines_and_id_only(tmp_path):
    path = tmp_path / "example.fa"
    path.write_text("\n>example\n\nACGT\n\n", encoding="utf-8")
    assert analyze_fasta(str(path))["record_id"] == "example"
    assert analyze_fasta(path)["description"] == "example"


@pytest.mark.parametrize("text", ["", "ACGT", "comment\n>id\nACGT", ">\nACGT",
                                 ">   \nACGT", ">id\n", ">id\nNNNN",
                                 ">id\nAC GT", ">id\nAC\tGT", " >id\nACGT",
                                 ">one\nACGT\n>two\nACGT", ">one\n>two\nACGT"])
def test_malformed_fasta(tmp_path, text):
    path = tmp_path / "invalid.fa"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError):
        analyze_fasta(path)


def test_missing_fasta(tmp_path):
    with pytest.raises(FileNotFoundError, match="FASTA file not found"):
        analyze_fasta(tmp_path / "missing.fa")


def test_invalid_fasta_encoding(tmp_path):
    path = tmp_path / "invalid.fa"
    path.write_bytes(b">id\nACG\xff")
    with pytest.raises(ValueError, match="UTF-8"):
        analyze_fasta(path)


@pytest.mark.parametrize("text", [">id\nAC\vGT", ">id\nAC\fGT",
                                 ">id\nAC\x85GT", ">id\nAC\u2028GT",
                                 ">id\nAC\n\v\nGT", "> id\nACGT",
                                 ">>id\nACGT", ">id\x00description\nACGT"])
def test_fasta_control_characters_and_invalid_headers(tmp_path, text):
    path = tmp_path / "invalid.fa"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError):
        analyze_fasta(path)


def test_fasta_preserves_description_spacing_and_cr_endings(tmp_path):
    path = tmp_path / "provenance.fa"
    path.write_bytes(b">sample_1  synthetic\tdescription  \rAC\rGT\r")
    result = analyze_fasta(path)
    assert result["record_id"] == "sample_1"
    assert result["description"] == "sample_1  synthetic\tdescription  "
    assert result["sequence"] == "ACGT"


@pytest.mark.parametrize("sequence,counts,gc,comp,rev", [
    ("A", {"A": 1, "T": 0, "G": 0, "C": 0}, 0, "T", "T"),
    ("GGC", {"A": 0, "T": 0, "G": 2, "C": 1}, 100, "CCG", "GCC"),
    ("AGT", {"A": 1, "T": 1, "G": 1, "C": 0}, 100 / 3, "TCA", "ACT"),
])
def test_independent_short_references(sequence, counts, gc, comp, rev):
    result = analyze_sequence(sequence)
    assert result["length"] == sum(counts.values())
    assert result["nucleotide_counts"] == counts
    assert result["gc_percentage"] == pytest.approx(gc)
    assert result["complement"] == comp
    assert result["reverse_complement"] == rev
