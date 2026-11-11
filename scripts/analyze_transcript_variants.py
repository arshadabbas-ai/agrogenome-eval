"""Offline, deterministic examples; no natural-allele or phenotype claims.

Run from repository root using the installed editable package:
  python scripts/analyze_transcript_variants.py --create-examples
  python scripts/analyze_transcript_variants.py
Existing examples and result directories are never overwritten.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory

from agrogenome_eval.cds_analysis import analyze_cds, parse_cds_annotation, parse_protein_fasta
from agrogenome_eval.transcript_variant import analyze_transcript_variant
from agrogenome_eval.variant_classifier import classify_variant
from analyze_coding_sequences import verify_derived_outputs
from fetch_reference_transcripts import verify_saved_reference


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def encode(data: dict) -> bytes:
    return (json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")


def load_verified_inputs(reference: Path, derived: Path) -> tuple[dict, dict]:
    """Bind authentic FASTA, XML, protein, manifest and CDS evidence offline."""
    reference, derived = Path(reference), Path(derived)
    verify_derived_outputs(derived, reference)
    # Use the byte-preserved manifest snapshot as the actual metadata input.
    # The verifier also checks semantic equality with data/reference/manifest.json.
    manifest = json.loads((derived / "reference_manifest.snapshot.json").read_bytes())
    summary = json.loads((derived / "summary.json").read_bytes())
    summary_ids = [r["transcript_accession"] for r in summary["records"]]
    manifest_ids = [r["transcript_accession"] for r in manifest["records"]]
    if len(summary_ids) != len(set(summary_ids)) or len(manifest_ids) != len(set(manifest_ids)):
        raise ValueError("Duplicate transcript accessions in reference/CDS evidence.")
    if set(summary_ids) != set(manifest_ids):
        raise ValueError("CDS evidence and reference manifest accession sets differ.")
    contexts, checksums = {}, {}
    for path in sorted(derived.iterdir()):
        if path.is_file(): checksums["cds_evidence/" + path.name] = digest(path.read_bytes())
    for record in manifest["records"]:
        accession = record["transcript_accession"]
        if accession in contexts:
            raise ValueError("Duplicate reference transcript accession.")
        evidence = next(r for r in summary["records"] if r["transcript_accession"] == accession)
        if evidence["status"] != "validated":
            raise ValueError("Reference CDS/protein evidence is not validated.")
        verify_saved_reference(reference / record["filename"], record)
        saved = json.loads((derived / evidence["analysis_file"]).read_bytes())
        prefix = record["species"] + "_" + accession
        transcript = record["analysis"]["sequence"]
        annotation = parse_cds_annotation((derived / (prefix + ".annotation.xml")).read_bytes(), accession, transcript,
            expected_organism=record["scientific_name"], expected_gene=record["gene"], expected_gene_id=record["gene_id"])
        if annotation["protein_accession"] != record["protein_accession"]:
            raise ValueError("Protein accession conflicts with manifest.")
        protein = parse_protein_fasta((derived / (prefix + ".reference_protein.fasta")).read_bytes(), record["protein_accession"])
        verified = analyze_cds(transcript, annotation, protein)
        if verified["validation_status"] != "validated" or any(saved.get(k) != v for k, v in verified.items()):
            raise ValueError("Reconstructed CDS evidence differs from saved analysis.")
        for key in ("scientific_name", "gene_id", "refseq_status"):
            if saved["provenance"][key] != record[key]:
                raise ValueError(f"Reference provenance mismatch: {key}.")
        annotation = {**annotation, "organism": record["scientific_name"], "gene_id": record["gene_id"],
                      "transcript_length": len(transcript),
                      "transcript_sequence_sha256": digest(transcript.encode("ascii")),
                      "provenance": {"transcript_accession": accession, "protein_accession": record["protein_accession"],
                                     "gene_id": record["gene_id"], "scientific_name": record["scientific_name"],
                                     "refseq_status": record["refseq_status"]}}
        contexts[accession] = {"transcript": transcript, "annotation": annotation,
                               "source_url": record["source_url"], "record_title": record["record_title"],
                               "annotation_limitations": record["limitations"], "species": record["species"]}
        checksums["reference/" + record["filename"]] = digest((reference / record["filename"]).read_bytes())
    return contexts, checksums


def select_hypothetical_examples(contexts: dict) -> dict:
    """Select the first position/ACGT alternative for each available category.

    Selection uses the existing classifier; tests independently verify every
    selected coding change against a separate standard-code lookup.
    """
    variants = []
    for accession, context in contexts.items():
        transcript, annotation = context["transcript"], context["annotation"]
        start, end = annotation["start"], annotation["end"]
        cds = transcript[start - 1:end]
        coding = {}
        for cds_pos in range(4, len(cds) - 2):
            for alt in "ACGT":
                if alt == cds[cds_pos - 1]: continue
                consequence = classify_variant(cds, cds_pos, alt)["Classification"]
                coding.setdefault(consequence, (start + cds_pos - 1, alt))
            if len(coding) == 3: break
        requested = [(name, *coding[name]) for name in ("SYNONYMOUS", "MISSENSE", "STOP_GAINED") if name in coding]
        for name, pos in [("FIVE_PRIME_UTR", 1 if start > 1 else None),
                          ("THREE_PRIME_UTR", end + 1 if end < len(transcript) else None),
                          ("START_CODON", start), ("TERMINAL_STOP_CODON", end)]:
            if pos is not None:
                requested.append((name, pos, next(base for base in "ACGT" if base != transcript[pos - 1])))
        for category, position, alt in requested:
            assessed = category in coding
            boundary = category in {"START_CODON", "TERMINAL_STOP_CODON"}
            variants.append({"id": context["species"] + "_" + category.lower(),
                "transcript_accession": accession, "position": position,
                "reference": transcript[position - 1], "alternative": alt,
                "origin": "computationally_introduced_hypothetical_substitution",
                "expected": {"transcript_region": "CDS" if assessed or boundary else category,
                             "molecular_consequence": category if assessed else None,
                             "analysis_status": "ASSESSED" if assessed else "UNSUPPORTED" if boundary else "NOT_ASSESSED",
                             "unsupported_region": category if boundary else None}})
    return {"schema_version": 1, "description": "Hypothetical computational substitutions, not experimentally observed natural mutations or verified crop alleles.",
            "selection_rule": "First internal position then A,C,G,T alternative for each classification; UTR/boundary positions selected explicitly.",
            "variants": variants}


def save_example_inputs(path: Path, definitions: dict) -> None:
    """Publish complete example bytes with an exclusive atomic hard link.

    A filesystem supporting hard links is required (e.g. NTFS/ext4). Failure
    leaves no partial target; existing targets are never replaced.
    """
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"Preserving existing hypothetical inputs: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=path.parent, prefix=".examples-") as temp:
        staged = Path(temp) / "examples.json"
        data = encode(definitions)
        staged.write_bytes(data)
        if staged.read_bytes() != data:
            raise OSError("Example input write did not preserve all bytes.")
        os.link(staged, path)


def verify_outputs(directory: Path, reference: Path | None = None,
                   derived: Path | None = None, examples: Path | None = None) -> None:
    """Verify complete output coverage, snapshot linkage and optional source bytes.

    Hashes detect byte changes, not NCBI authenticity. Source verification also
    reconstructs the accession/organism/gene/protein evidence chain offline.
    """
    manifest = json.loads((directory / "manifest.json").read_bytes())
    required = {"hypothetical_variants.snapshot.json", "results.json"}
    if set(manifest.get("output_checksums", {})) != required:
        raise ValueError("Incomplete or unexpected variant output checksum coverage.")
    for name, checksum in manifest["output_checksums"].items():
        if digest((directory / name).read_bytes()) != checksum:
            raise ValueError(f"Variant output checksum mismatch: {name}.")
    snapshot = (directory / "hypothetical_variants.snapshot.json").read_bytes()
    if digest(snapshot) != manifest.get("input_checksums", {}).get("hypothetical_variants.json"):
        raise ValueError("Hypothetical input snapshot checksum mismatch.")
    definitions = json.loads(snapshot)
    payload = json.loads((directory / "results.json").read_bytes())
    if [row["example"] for row in payload["results"]] != definitions["variants"]:
        raise ValueError("Results do not preserve all original hypothetical definitions.")
    if any(row["example"].get("origin") != "computationally_introduced_hypothetical_substitution" for row in payload["results"]):
        raise ValueError("Result origin must explicitly remain hypothetical.")
    if reference is not None or derived is not None or examples is not None:
        if reference is None or derived is None or examples is None:
            raise ValueError("Source verification requires reference, derived and examples paths.")
        _, checksums = load_verified_inputs(reference, derived)
        if examples.read_bytes() != snapshot or manifest["input_checksums"] != {**checksums, "hypothetical_variants.json": digest(snapshot)}:
            raise ValueError("Source input checksums differ from the variant output manifest.")


def generate_results(examples: Path, output: Path, reference: Path, derived: Path) -> dict:
    """Verify all pinned evidence, assess explicit inputs, then publish atomically."""
    if output.exists():
        raise FileExistsError(f"Preserving existing variant evidence: {output}")
    contexts, checksums = load_verified_inputs(reference, derived)
    raw = examples.read_bytes()
    definitions = json.loads(raw)
    if not isinstance(definitions, dict) or definitions.get("schema_version") != 1 or not isinstance(definitions.get("variants"), list) or not definitions["variants"]:
        raise ValueError("Expected a nonempty schema-version-1 hypothetical variant list.")
    if not isinstance(definitions.get("description"), str) or not definitions["description"]:
        raise ValueError("Example description is required.")
    results, identifiers = [], set()
    for variant in definitions["variants"]:
        required = {"id", "transcript_accession", "position", "reference", "alternative", "origin", "expected"}
        if not isinstance(variant, dict) or not required <= variant.keys():
            raise ValueError("Missing hypothetical variant definition fields.")
        if not isinstance(variant["id"], str) or not variant["id"] or variant["id"] in identifiers:
            raise ValueError("Variant IDs must be unique nonempty strings.")
        identifiers.add(variant["id"])
        if variant["origin"] != "computationally_introduced_hypothetical_substitution":
            raise ValueError("Example origin must explicitly identify a hypothetical computational substitution.")
        if not isinstance(variant["transcript_accession"], str) or variant["transcript_accession"] not in contexts:
            raise ValueError("Variant transcript accession/version is not in verified reference evidence.")
        context = contexts[variant["transcript_accession"]]
        result = analyze_transcript_variant(context["transcript"], context["annotation"], variant["position"], variant["reference"], variant["alternative"])
        expected = variant["expected"]
        expected_keys = {"transcript_region", "molecular_consequence", "analysis_status", "unsupported_region"}
        if not isinstance(expected, dict) or set(expected) != expected_keys or any(result[k] != value for k, value in expected.items()):
            raise ValueError(f"Example expectation mismatch: {variant['id']}.")
        results.append({"example": variant, "result": result,
                        "source_url": context["source_url"], "record_title": context["record_title"],
                        "annotation_limitations": context["annotation_limitations"]})
    payload = {"schema_version": 1, "description": definitions["description"], "results": results,
               "limitations": ["All examples are computationally introduced, not natural crop alleles.",
                               "No expression, splicing, RNA-stability, phenotype or crop-performance inference."]}
    files = {"hypothetical_variants.snapshot.json": raw, "results.json": encode(payload)}
    manifest = {"schema_version": 1, "input_checksums": {**checksums, "hypothetical_variants.json": digest(raw)},
                "output_checksums": {name: digest(data) for name, data in files.items()},
                "source_roots": {"reference": "data/reference", "cds_evidence": "data/derived/cds_protein"},
                "reproducibility": "No timestamps or network requests; identical verified input bytes produce identical output bytes."}
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent, prefix=".variant-") as temp:
        staging = Path(temp) / "verified"; staging.mkdir()
        for name, data in files.items(): (staging / name).write_bytes(data)
        (staging / "manifest.json").write_bytes(encode(manifest))
        verify_outputs(staging, reference, derived, examples)
        if output.exists(): raise FileExistsError("Variant output appeared during analysis.")
        staging.rename(output)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=Path("data/reference"))
    parser.add_argument("--cds-evidence", type=Path, default=Path("data/derived/cds_protein"))
    parser.add_argument("--examples", type=Path, default=Path("data/examples/hypothetical_variants.json"))
    parser.add_argument("--output", type=Path, default=Path("data/derived/variant_effects"))
    parser.add_argument("--create-examples", action="store_true")
    args = parser.parse_args()
    try:
        if args.create_examples:
            contexts, _ = load_verified_inputs(args.reference, args.cds_evidence)
            definitions = select_hypothetical_examples(contexts)
            save_example_inputs(args.examples, definitions)
            print(f"Created {len(definitions['variants'])} explicit hypothetical definitions.")
        else:
            result = generate_results(args.examples, args.output, args.reference, args.cds_evidence)
            print(f"Analyzed {len(result['results'])} hypothetical substitutions.")
    except (ValueError, OSError) as exc:
        parser.exit(1, f"Variant analysis stopped: {exc}\n")


if __name__ == "__main__":
    main()
