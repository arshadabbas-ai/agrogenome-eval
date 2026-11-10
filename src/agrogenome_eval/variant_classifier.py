
# AgroGenome-Eval
# Scientific variant classifier — Version 0.1
# Supports internal single-nucleotide substitutions
# in complete protein-coding sequences.

from Bio.Seq import Seq


def classify_variant(cds, position, alternative, table=1):

    cds = str(cds).upper()
    alternative = str(alternative).upper()

    # Validate the coding sequence
    if len(cds) < 9 or len(cds) % 3 != 0:
        raise ValueError("Invalid CDS length.")

    if any(base not in "ACGT" for base in cds):
        raise ValueError("Unsupported nucleotide detected.")

    # Validate mutation position
    if type(position) is not int:
        raise ValueError("Position must be an integer.")

    if not 1 <= position <= len(cds):
        raise ValueError("Position outside CDS.")

    if len(alternative) != 1 or alternative not in "ACGT":
        raise ValueError("Alternative must be A, C, G or T.")

    # Validate complete reference CDS
    Seq(cds).translate(table=table, cds=True)

    index = position - 1
    codon_start = (index // 3) * 3
    codon_number = codon_start // 3 + 1

    # Handle internal codons only
    if codon_number in (1, len(cds) // 3):
        raise ValueError(
            "Start and terminal stop codons not supported."
        )

    original_base = cds[index]

    if original_base == alternative:
        raise ValueError("No nucleotide change detected.")

    original_codon = cds[codon_start:codon_start + 3]
    offset = index % 3

    modified_codon = (
        original_codon[:offset]
        + alternative
        + original_codon[offset + 1:]
    )

    original_aa = str(
        Seq(original_codon).translate(table=table)
    )

    modified_aa = str(
        Seq(modified_codon).translate(table=table)
    )

    # Predict coding consequence
    if modified_aa == "*":
        classification = "STOP_GAINED"

    elif original_aa == modified_aa:
        classification = "SYNONYMOUS"

    else:
        classification = "MISSENSE"

    return {
        "CDS_Position": position,
        "Reference_Base": original_base,
        "Alternative_Base": alternative,
        "Codon_Number": codon_number,
        "Original_Codon": original_codon,
        "Modified_Codon": modified_codon,
        "Original_AA": original_aa,
        "Modified_AA": modified_aa,
        "Classification": classification
    }
