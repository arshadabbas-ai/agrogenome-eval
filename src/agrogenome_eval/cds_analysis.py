"""Annotated transcript CDS validation for a deliberately restricted subset.

Supports exactly one complete positive-strand interval, explicit codon_start=1,
and explicit NCBI translation table 1. Coordinates are 1-based and inclusive;
Biopython extraction uses 0-based, half-open coordinates internally. Missing
qualifiers, partial/compound locations, translation exceptions, ambiguous DNA
or proteins, and unsupported reading frames/tables raise ValueError. This is
a GB XML subset reader, not a general GenBank parser or gene predictor.
"""

import re
import xml.etree.ElementTree as ET

from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq
from Bio.SeqFeature import SimpleLocation

from .sequence_analysis import validate_dna


AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def _protein(sequence: str) -> str:
    if not isinstance(sequence, str) or not sequence or set(sequence) - AMINO_ACIDS:
        raise ValueError("Protein must contain only the 20 canonical uppercase amino acids, without stop symbols.")
    return sequence


def parse_protein_fasta(data: bytes, expected_accession: str) -> dict:
    """Read exactly one version-matched protein FASTA, without SeqIO imports."""
    try:
        text = data.decode("utf-8").replace("\r\n", "\n").replace("\r", "\n")
        lines = [line for line in text.split("\n") if line]
    except UnicodeDecodeError as exc:
        raise ValueError("Protein FASTA must be UTF-8.") from exc
    if not lines or not lines[0].startswith(">") or sum(line.startswith(">") for line in lines) != 1:
        raise ValueError("Expected one protein FASTA record.")
    description = lines[0][1:]
    if not description or description[0].isspace() or description.split()[0] != expected_accession:
        raise ValueError("Protein FASTA accession/version mismatch.")
    return {"accession": expected_accession, "description": description,
            "sequence": _protein("".join(lines[1:]))}


def _bounds(transcript: str, annotation: dict) -> tuple[str, int, int]:
    dna = validate_dna(transcript)
    if not isinstance(annotation, dict) or not {"start", "end", "strand", "codon_start", "translation_table"} <= annotation.keys():
        raise ValueError("Missing required CDS coordinate or translation metadata.")
    start, end = annotation["start"], annotation["end"]
    if type(start) is not int or type(end) is not int or not 1 <= start <= end <= len(dna):
        raise ValueError("CDS coordinates must be integers within the transcript.")
    if any(type(annotation[key]) is not int or annotation[key] != 1
           for key in ("strand", "codon_start", "translation_table")):
        raise ValueError("Only positive-strand, codon_start=1, translation table 1 CDS is supported.")
    return dna, start, end


def parse_cds_annotation(xml: bytes, accession: str, transcript: str, *,
                         expected_organism: str | None = None,
                         expected_gene: str | None = None,
                         expected_gene_id: str | None = None) -> dict:
    """Match full transcript sequence/version to GB XML and return CDS facts.

    Optional expected identity fields constrain organism and the selected CDS's
    gene/GeneID. The retrieval workflow supplies all three from the pinned
    reference manifest; a raw sequence alone does not identify an organism.
    """
    dna = validate_dna(transcript)
    if not re.fullmatch(r"[A-Z]+_[0-9]+\.[0-9]+", accession):
        raise ValueError("Version-specific transcript accession is required.")
    try:
        root = ET.fromstring(xml)
    except ET.ParseError as exc:
        raise ValueError("Malformed GB XML annotations.") from exc
    records = root.findall(".//GBSeq")
    if len(records) != 1 or root.findall(".//ERROR"):
        raise ValueError("Expected one error-free GB XML nucleotide record.")
    record = records[0]
    # findtext otherwise silently uses the first of conflicting duplicate tags.
    def unique_children(element, tags):
        for tag in tags:
            if len(element.findall(tag)) > 1:
                raise ValueError(f"Duplicate/conflicting XML field: {tag}.")
    def flag_is_false(flag):
        attribute, text = flag.get("value"), (flag.text or "").strip()
        if attribute and text and attribute.lower() != text.lower():
            raise ValueError("Conflicting XML boolean annotation.")
        value = (attribute or text).lower()
        if value not in {"true", "false"}:
            raise ValueError("Malformed XML boolean annotation.")
        return value == "false"
    unique_children(record, {"GBSeq_accession-version", "GBSeq_moltype", "GBSeq_sequence",
                             "GBSeq_length", "GBSeq_organism", "GBSeq_definition", "GBSeq_feature-table"})
    if record.findtext("GBSeq_accession-version") != accession:
        raise ValueError("Transcript annotation accession/version mismatch.")
    if record.findtext("GBSeq_moltype") != "mRNA":
        raise ValueError("Annotation must describe an mRNA transcript.")
    if expected_organism is not None and record.findtext("GBSeq_organism") != expected_organism:
        raise ValueError("CDS annotation organism mismatch.")
    if record.findtext("GBSeq_sequence", "").upper() != dna or record.findtext("GBSeq_length") != str(len(dna)):
        raise ValueError("Transcript annotation sequence or length mismatch.")
    all_features = record.findall("GBSeq_feature-table/GBFeature")
    for candidate in all_features:
        unique_children(candidate, {"GBFeature_key"})
    features = [f for f in all_features if f.findtext("GBFeature_key") == "CDS"]
    if len(features) != 1:
        raise ValueError("Exactly one annotated CDS is required.")
    feature = features[0]
    unique_children(feature, {"GBFeature_location", "GBFeature_intervals", "GBFeature_quals",
                              "GBFeature_partial5", "GBFeature_partial3", "GBFeature_operator"})
    if feature.find("GBFeature_operator") is not None:
        raise ValueError("Compound CDS feature operators are unsupported.")
    location = feature.findtext("GBFeature_location", "")
    if not re.fullmatch(r"[1-9][0-9]*\.\.[1-9][0-9]*", location):
        raise ValueError("Unsupported partial, compound, reverse or malformed CDS location.")
    for tag in ("GBFeature_partial5", "GBFeature_partial3"):
        flag = feature.find(tag)
        if flag is not None and not flag_is_false(flag):
            raise ValueError("Partial CDS is unsupported.")
    intervals = feature.findall("GBFeature_intervals/GBInterval")
    if len(intervals) != 1:
        raise ValueError("CDS must contain exactly one interval.")
    interval = intervals[0]
    unique_children(interval, {child.tag for child in interval})
    if any(child.tag not in {"GBInterval_from", "GBInterval_to", "GBInterval_accession", "GBInterval_iscomp"}
           for child in interval):
        raise ValueError("Unsupported CDS interval structure.")
    if interval.findtext("GBInterval_accession") != accession or interval.find("GBInterval_point") is not None:
        raise ValueError("CDS interval accession or structure mismatch.")
    reverse = interval.find("GBInterval_iscomp")
    if reverse is not None and not flag_is_false(reverse):
        raise ValueError("Reverse-strand CDS is unsupported.")
    try:
        start = int(interval.findtext("GBInterval_from", ""))
        end = int(interval.findtext("GBInterval_to", ""))
    except ValueError as exc:
        raise ValueError("Missing or malformed CDS interval coordinates.") from exc
    if location != f"{start}..{end}":
        raise ValueError("CDS location string and interval coordinates disagree.")
    qualifiers = {}
    for qualifier in feature.findall("GBFeature_quals/GBQualifier"):
        unique_children(qualifier, {"GBQualifier_name", "GBQualifier_value"})
        qualifiers.setdefault(qualifier.findtext("GBQualifier_name"), []).append(qualifier.findtext("GBQualifier_value", ""))
    if set(qualifiers) & {"exception", "transl_except", "ribosomal_slippage", "pseudo", "pseudogene", "trans_splicing"}:
        raise ValueError("Exceptional or pseudogene CDS translation is unsupported.")
    def required(name):
        values = qualifiers.get(name, [])
        if len(values) != 1 or not values[0]:
            raise ValueError(f"Missing or duplicate CDS qualifier: {name}.")
        return values[0]
    codon_start, table = required("codon_start"), required("transl_table")
    if codon_start != "1" or table != "1":
        raise ValueError("Only explicit codon_start=1 and translation table 1 are supported.")
    protein_id = required("protein_id")
    if not re.fullmatch(r"[A-Z]+_[0-9]+\.[0-9]+", protein_id):
        raise ValueError("Version-specific protein accession is required.")
    if expected_gene is not None and required("gene") != expected_gene:
        raise ValueError("Selected CDS gene association mismatch.")
    if expected_gene_id is not None:
        gene_ids = {value for value in qualifiers.get("db_xref", []) if value.startswith("GeneID:")}
        if gene_ids != {"GeneID:" + expected_gene_id}:
            raise ValueError("Selected CDS GeneID association mismatch.")
    annotation = {"transcript_accession": accession, "start": start, "end": end,
                  "coordinate_system": "1-based inclusive", "location": location,
                  "strand": 1, "codon_start": 1, "translation_table": 1,
                  "protein_accession": protein_id, "reference_translation": _protein(required("translation")),
                  "gene": required("gene"), "product": required("product"),
                  "qualifiers": qualifiers, "record_title": record.findtext("GBSeq_definition", "")}
    _bounds(dna, annotation)
    return annotation


def extract_cds(transcript: str, annotation: dict) -> str:
    """Extract the inclusive annotated CDS; do not infer coding boundaries."""
    dna, start, end = _bounds(transcript, annotation)
    return str(SimpleLocation(start - 1, end, strand=1).extract(Seq(dna)))


def utr_lengths(transcript: str, annotation: dict) -> dict[str, int]:
    """Return transcript flanking lengths for the supported complete CDS subset."""
    dna, start, end = _bounds(transcript, annotation)
    return {"five_prime_utr_length": start - 1, "three_prime_utr_length": len(dna) - end}


def translate_cds(cds: str, translation_table: int = 1) -> str:
    """Validate a complete CDS and translate; terminal stop is not an amino acid."""
    dna = validate_dna(cds)
    if type(translation_table) is not int or translation_table != 1:
        raise ValueError("Only NCBI translation table 1 is supported.")
    try:
        return _protein(str(Seq(dna).translate(table=translation_table, cds=True)))
    except TranslationError as exc:
        raise ValueError(f"Invalid complete CDS: {exc}") from exc


def analyze_cds(transcript: str, annotation: dict, reference_protein: dict | None = None) -> dict:
    """Return coding facts and explicit comparison status, including mismatches.

    A status of validated refers to sequence comparisons, not RefSeq curation.
    Reference protein comparison is not claimed when a protein FASTA is absent.
    """
    cds = extract_cds(transcript, annotation)
    protein = translate_cds(cds, annotation["translation_table"])
    expected = _protein(annotation["reference_translation"])
    discrepancies = []
    annotated_match = protein == expected
    if not annotated_match:
        discrepancies.append("Computed translation differs from annotated CDS translation.")
    protein_match = None
    if reference_protein is not None:
        if reference_protein["accession"] != annotation["protein_accession"]:
            raise ValueError("Reference protein accession/version mismatch.")
        protein_match = protein == _protein(reference_protein["sequence"])
        if not protein_match:
            discrepancies.append("Computed translation differs from version-specific protein FASTA.")
    return {"transcript_accession": annotation["transcript_accession"],
            "annotation": annotation, "transcript_length": len(validate_dna(transcript)),
            "cds_sequence": cds, "cds_length": len(cds), **utr_lengths(transcript, annotation),
            "start_codon": cds[:3], "terminal_stop_codon": cds[-3:],
            "protein_sequence": protein, "protein_length": len(protein),
            "annotated_translation_match": annotated_match, "reference_protein_match": protein_match,
            "validation_status": "mismatch" if discrepancies else ("validated" if protein_match else "annotation_only"),
            "discrepancies": discrepancies,
            "limitations": ["Complete positive-strand single-interval CDS; explicit codon_start=1 and table 1 only.",
                            "Sequence agreement does not establish phenotype, function, or orthology."]}
