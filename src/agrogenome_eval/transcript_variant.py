"""Transcript-relative SNV mapping to a validated, annotated CDS.

This wrapper delegates all internal-codon consequence classification to the
existing classifier. UTR and boundary-codon consequences are not assessed.
Unsupported or not assessed does not imply biological harmlessness.
Positions are 1-based inclusive, never genomic coordinates. Accession identity
cannot be inferred from a DNA string: the offline runner verifies XML, FASTA,
and manifest provenance, then supplies a sequence fingerprint and identity.
"""

import hashlib
import json
import re

from .cds_analysis import analyze_cds
from .sequence_analysis import validate_dna
from .variant_classifier import classify_variant


def _allele(value: str, name: str) -> str:
    if not isinstance(value, str) or len(value) != 1 or value not in "ACGTacgt":
        raise ValueError(f"{name} must be one A/C/G/T nucleotide.")
    return value.upper()


def _validate_annotation(transcript: str, annotation: dict) -> dict:
    required = {"transcript_accession", "protein_accession", "start", "end", "strand",
                "codon_start", "translation_table", "reference_translation",
                "coordinate_system", "location"}
    if not isinstance(annotation, dict) or not required <= annotation.keys():
        raise ValueError("Missing structured CDS annotation metadata.")
    for key in ("transcript_accession", "protein_accession"):
        if not isinstance(annotation[key], str) or not re.fullmatch(r"[A-Z]+_[0-9]+\.[0-9]+", annotation[key]):
            raise ValueError(f"Version-specific {key} is required.")
    if annotation["coordinate_system"] != "1-based inclusive":
        raise ValueError("Only 1-based inclusive transcript CDS coordinates are supported.")
    if annotation["location"] != f"{annotation['start']}..{annotation['end']}":
        raise ValueError("Conflicting CDS location and coordinates.")
    if "transcript_length" in annotation and (type(annotation["transcript_length"]) is not int
                                              or annotation["transcript_length"] != len(transcript)):
        raise ValueError("Transcript length differs from annotation provenance.")
    if "transcript_sequence_sha256" in annotation and annotation["transcript_sequence_sha256"] != hashlib.sha256(transcript.encode("ascii")).hexdigest():
        raise ValueError("Transcript sequence fingerprint mismatch.")
    provenance = annotation.get("provenance", {})
    if not isinstance(provenance, dict):
        raise ValueError("Annotation provenance must be a dictionary.")
    for key in ("gene_id", "scientific_name", "refseq_status"):
        if provenance.get(key) is not None and (not isinstance(provenance[key], str) or not provenance[key]):
            raise ValueError(f"Provenance {key} must be a nonempty string or null.")
    try:
        json.dumps(provenance, allow_nan=False)
    except (ValueError, TypeError) as exc:
        raise ValueError("Annotation provenance must be JSON-serializable.") from exc
    for key in ("transcript_accession", "protein_accession"):
        if key in provenance and provenance[key] != annotation[key]:
            raise ValueError(f"Annotation/provenance {key} mismatch.")
    # Compare identities only against explicitly supplied annotation evidence;
    # neither a sequence fingerprint nor composition establishes species.
    if provenance.get("scientific_name") is not None and "organism" in annotation and provenance["scientific_name"] != annotation["organism"]:
        raise ValueError("Annotation/provenance scientific_name mismatch.")
    if provenance.get("gene_id") is not None and "gene_id" in annotation and provenance["gene_id"] != annotation["gene_id"]:
        raise ValueError("Annotation/provenance gene_id mismatch.")
    if "qualifiers" in annotation:
        qualifiers = annotation["qualifiers"]
        if not isinstance(qualifiers, dict):
            raise ValueError("CDS qualifiers must be a dictionary.")
        if set(qualifiers) & {"exception", "transl_except", "pseudo", "pseudogene", "ribosomal_slippage", "trans_splicing"}:
            raise ValueError("Exceptional CDS annotations are unsupported.")
        for qualifier, key in {"codon_start": "codon_start", "transl_table": "translation_table",
                               "protein_id": "protein_accession", "translation": "reference_translation",
                               "gene": "gene"}.items():
            if key not in annotation:
                if qualifier in qualifiers:
                    raise ValueError(f"Missing CDS metadata for qualifier: {qualifier}.")
                continue
            if qualifier in qualifiers and qualifiers[qualifier] != [str(annotation[key])]:
                raise ValueError(f"Conflicting CDS qualifier: {qualifier}.")
        if provenance.get("gene_id") is not None and "db_xref" in qualifiers:
            xrefs = qualifiers["db_xref"]
            if not isinstance(xrefs, list) or any(not isinstance(value, str) for value in xrefs):
                raise ValueError("CDS db_xref qualifiers must be a list of strings.")
            gene_ids = {value for value in xrefs if value.startswith("GeneID:")}
            if gene_ids != {"GeneID:" + provenance["gene_id"]}:
                raise ValueError("CDS GeneID and provenance mismatch.")
    # Reuse complete-CDS validation, including annotated protein agreement.
    try:
        result = analyze_cds(transcript, annotation)
    except (KeyError, TypeError) as exc:
        raise ValueError("Malformed CDS annotation.") from exc
    if not result["annotated_translation_match"]:
        raise ValueError("Reference CDS translation disagrees with annotation.")
    return result


def analyze_transcript_variant(transcript: str, annotation: dict, position: int,
                               reference: str, alternative: str) -> dict:
    """Analyze a single hypothetical or observed substitution without phenotype inference.

    `annotation` follows parse_cds_annotation's dictionary schema. Optional
    provenance and transcript fingerprint fields are checked when supplied.
    Optional `organism`/`gene_id` fields and GeneID qualifiers are cross-checked
    against provenance. Only the runner binds these identities to verified XML;
    a pure function cannot prove an accession/organism from sequence alone.
    No species, curation status, or experimental origin is inferred. Returned
    `analysis_status` is independent of the optional NCBI `refseq_status`.
    """
    dna = validate_dna(transcript)
    if type(position) is not int or not 1 <= position <= len(dna):
        raise ValueError("Transcript position must be an integer within sequence boundaries.")
    reference, alternative = _allele(reference, "Reference allele"), _allele(alternative, "Alternative allele")
    if dna[position - 1] != reference:
        raise ValueError("Reference allele does not match the transcript nucleotide.")
    if reference == alternative:
        raise ValueError("Reference and alternative alleles must differ.")
    coding = _validate_annotation(dna, annotation)
    start, end = annotation["start"], annotation["end"]
    provenance = annotation.get("provenance", {})
    result = {
        "transcript_accession": annotation["transcript_accession"],
        "gene_id": provenance.get("gene_id"), "scientific_name": provenance.get("scientific_name"),
        "protein_accession": annotation["protein_accession"],
        "refseq_status": provenance.get("refseq_status"),
        "transcript_position": position, "reference": reference, "alternative": alternative,
        "coordinate_system": "1-based inclusive transcript", "cds_start": start, "cds_end": end,
        "transcript_region": None, "cds_position": None, "codon_number": None,
        "reference_codon": None, "alternative_codon": None,
        "reference_amino_acid": None, "alternative_amino_acid": None,
        "molecular_consequence": None, "analysis_status": "NOT_ASSESSED",
        "unsupported_region": None, "not_assessed_reason": None,
        "limitations": ["Single-nucleotide substitutions in annotated transcripts only; not genomic coordinates.",
                        "Internal CDS consequences describe sequence-level changes, not function, phenotype, or orthology.",
                        "UTR, start-codon, and terminal-stop-codon protein consequences are not assessed."]}
    if position < start or position > end:
        result["transcript_region"] = "FIVE_PRIME_UTR" if position < start else "THREE_PRIME_UTR"
        result["not_assessed_reason"] = "UTR substitution: no protein consequence or RNA/regulatory effect is assessed."
        return result
    result["transcript_region"] = "CDS"
    cds_position = position - start + 1
    codon_number = (cds_position - 1) // 3 + 1
    result.update(cds_position=cds_position, codon_number=codon_number)
    if codon_number == 1 or codon_number == len(coding["cds_sequence"]) // 3:
        result["unsupported_region"] = "START_CODON" if codon_number == 1 else "TERMINAL_STOP_CODON"
        result["analysis_status"] = "UNSUPPORTED"
        result["not_assessed_reason"] = "The existing classifier excludes this boundary codon; no start-lost or stop-lost prediction is made."
        return result
    consequence = classify_variant(coding["cds_sequence"], cds_position, alternative,
                                   table=annotation["translation_table"])
    result.update(reference_codon=consequence["Original_Codon"], alternative_codon=consequence["Modified_Codon"],
                  reference_amino_acid=consequence["Original_AA"], alternative_amino_acid=consequence["Modified_AA"],
                  molecular_consequence=consequence["Classification"], analysis_status="ASSESSED")
    return result
