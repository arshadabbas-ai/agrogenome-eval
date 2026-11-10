
# AgroGenome-Eval
# Automated scientific tests using pytest

import pytest

from agrogenome_eval import classify_variant


# Synthetic reference CDS
TEST_CDS = "ATGACGTTGTGA"


# Test 1-3: Known mutation consequences
@pytest.mark.parametrize(
    "position, alternative, expected",
    [
        (4, "C", "MISSENSE"),
        (6, "A", "SYNONYMOUS"),
        (8, "A", "STOP_GAINED"),
    ]
)
def test_variant_classification(
    position, alternative, expected
):
    result = classify_variant(
        TEST_CDS, position, alternative
    )

    assert result["Classification"] == expected
    assert result["CDS_Position"] == position


# Test 4-10: Invalid variant inputs
@pytest.mark.parametrize(
    "position, alternative",
    [
        (0, "C"),     # Invalid position
        (13, "C"),    # Beyond sequence
        (4, "N"),     # Invalid nucleotide
        (4, "A"),     # No mutation
        (1, "C"),     # Start codon excluded
        (10, "C"),    # Stop codon excluded
        (4.5, "C"),   # Noninteger position
    ]
)
def test_invalid_variants(position, alternative):

    with pytest.raises(ValueError):
        classify_variant(
            TEST_CDS, position, alternative
        )


# Test 11: Lowercase alternative
def test_lowercase_input():

    result = classify_variant(
        TEST_CDS.lower(), 4, "c"
    )

    assert result["Classification"] == "MISSENSE"
    assert result["Alternative_Base"] == "C"
