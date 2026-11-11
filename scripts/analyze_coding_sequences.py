"""Derive validated CDS/protein files from pinned transcript references.

Set NCBI_EMAIL in the process environment, then run from the repository root:
    python scripts/analyze_coding_sequences.py
Raw annotation XML and version-specific protein FASTA are retained as evidence.
Outputs go to a new directory; existing references and outputs are not replaced.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory

from agrogenome_eval.cds_analysis import analyze_cds, parse_cds_annotation, parse_protein_fasta
from fetch_reference_transcripts import NCBIClient, verify_metadata, verify_saved_reference


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def verify_derived_outputs(directory: Path, reference: Path) -> None:
    """Verify exact saved input/output bytes and the input-manifest snapshot."""
    summary = json.loads((directory / "summary.json").read_bytes())
    snapshot = (directory / "reference_manifest.snapshot.json").read_bytes()
    if sha256(snapshot) != summary["reference_manifest_sha256"]:
        raise ValueError("Input manifest snapshot checksum mismatch.")
    if json.loads(snapshot) != json.loads((reference / "manifest.json").read_bytes()):
        raise ValueError("Input manifest snapshot content mismatch.")
    for record in summary["records"]:
        if record["status"] != "validated":
            continue
        for name, checksum in record["checksums"].items():
            if sha256((directory / name).read_bytes()) != checksum:
                raise ValueError(f"Derived output checksum mismatch: {name}.")
        analysis = json.loads((directory / record["analysis_file"]).read_bytes())
        for name, checksum in analysis["provenance"]["input_checksums"].items():
            path = (directory / "reference_manifest.snapshot.json" if name == "manifest.json"
                    else directory / name if (directory / name).exists() else reference / name)
            if sha256(path.read_bytes()) != checksum:
                raise ValueError(f"Derived input checksum mismatch: {name}.")
        for name, checksum in analysis["provenance"]["output_checksums"].items():
            if sha256((directory / name).read_bytes()) != checksum:
                raise ValueError(f"Analysis output checksum mismatch: {name}.")


def generate_outputs(output: Path, reference: Path, client: NCBIClient) -> dict:
    """Generate offline-reproducible evidence with per-reference failure reporting."""
    output, reference = Path(output), Path(reference)
    if output.exists():
        raise FileExistsError(f"Preserving existing derived output: {output}")
    manifest_path = reference / "manifest.json"
    manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes)
    summary = {"schema_version": 1, "records": [], "reference_manifest_sha256": sha256(manifest_bytes),
               "limitations": ["Derived from transcript CDS annotations, not genomic coordinates.",
                               "RefSeq curation status is distinct from sequence validation status.",
                               "No functional equivalence, orthology, or phenotype inference."]}
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent, prefix=".cds-") as temporary:
        staging = Path(temporary) / "verified"
        staging.mkdir()
        # Preserve the exact input manifest bytes, whose checkout line endings
        # can differ between platforms, without changing the reference manifest.
        (staging / "reference_manifest.snapshot.json").write_bytes(manifest_bytes)
        for record in manifest["records"]:
            accession = record["transcript_accession"]
            saving = False
            try:
                transcript_path = reference / record["filename"]
                verify_saved_reference(transcript_path, record)
                # This helper verifies metadata organism, locus, GeneID, RefSeq
                # status and protein association in addition to the CDS parser.
                xml = client.fetch(accession, "gb", "xml")
                metadata = verify_metadata(xml, record)
                for key in ("record_title", "gene_description", "refseq_status"):
                    if metadata[key] != record[key]:
                        raise ValueError(f"Reference manifest annotation changed: {key}.")
                transcript = record["analysis"]["sequence"]
                annotation = parse_cds_annotation(xml, accession, transcript,
                    expected_organism=record["scientific_name"], expected_gene=record["gene"],
                    expected_gene_id=record["gene_id"])
                if annotation["protein_accession"] != record["protein_accession"]:
                    raise ValueError("CDS protein differs from reference manifest.")
                protein_bytes = client.fetch(annotation["protein_accession"], "fasta", "text", db="protein")
                reference_protein = parse_protein_fasta(protein_bytes, annotation["protein_accession"])
                result = analyze_cds(transcript, annotation, reference_protein)
                if result["validation_status"] != "validated":
                    raise ValueError("; ".join(result["discrepancies"]))
                prefix = record["species"] + "_" + accession
                files = {
                    prefix + ".annotation.xml": xml,
                    prefix + ".reference_protein.fasta": protein_bytes,
                    prefix + ".cds.fasta": (f">{accession}_CDS {annotation['location']} strand=+ table=1\n{result['cds_sequence']}\n").encode(),
                    prefix + ".protein.fasta": (f">{annotation['protein_accession']} translated_from={accession}\n{result['protein_sequence']}\n").encode(),
                }
                result["provenance"] = {
                    "scientific_name": record["scientific_name"], "gene_id": record["gene_id"],
                    "refseq_status": metadata["refseq_status"], "record_title": metadata["record_title"],
                    "reference_annotation_limitations": record["limitations"],
                    "transcript_source_url": record["source_url"],
                    "protein_source_url": "https://www.ncbi.nlm.nih.gov/protein/" + annotation["protein_accession"],
                    "retrieval_timestamp": datetime.now(timezone.utc).isoformat(),
                    "input_checksums": {record["filename"]: sha256(transcript_path.read_bytes()),
                                        "manifest.json": sha256(manifest_bytes),
                                        prefix + ".annotation.xml": sha256(xml),
                                        prefix + ".reference_protein.fasta": sha256(protein_bytes)},
                    "output_checksums": {name: sha256(data) for name, data in files.items()},
                    "reference_protein_description": reference_protein["description"],
                }
                analysis_name = prefix + ".analysis.json"
                files[analysis_name] = (json.dumps(result, indent=2) + "\n").encode()
                saving = True
                for name, data in files.items():
                    (staging / name).write_bytes(data)
                summary["records"].append({"transcript_accession": accession, "status": "validated",
                                           "analysis_file": analysis_name,
                                           "checksums": {name: sha256(data) for name, data in files.items()}})
            except (ValueError, OSError, RuntimeError) as exc:
                if saving:
                    # A local write failure aborts the entire staged output;
                    # never publish a partly saved successful record.
                    raise
                summary["records"].append({"transcript_accession": accession, "status": "failed", "error": str(exc)})
        (staging / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        verify_derived_outputs(staging, reference)
        if output.exists():
            raise FileExistsError(f"Output appeared during analysis: {output}")
        staging.rename(output)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=Path("data/reference"))
    parser.add_argument("--output", type=Path, default=Path("data/derived/cds_protein"))
    args = parser.parse_args()
    try:
        result = generate_outputs(args.output, args.reference, NCBIClient(os.environ.get("NCBI_EMAIL", "")))
    except (ValueError, OSError, RuntimeError) as exc:
        parser.exit(1, f"CDS analysis stopped: {exc}\n")
    print(json.dumps(result, indent=2))
    if any(record["status"] != "validated" for record in result["records"]):
        parser.exit(1, "One or more references failed validation; see summary.json.\n")


if __name__ == "__main__":
    main()
