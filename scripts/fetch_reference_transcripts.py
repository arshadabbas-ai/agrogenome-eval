"""Fetch and verify two versioned NCBI transcript references.

Run from the repository root after setting NCBI_EMAIL to a real contact email:
    python scripts/fetch_reference_transcripts.py
Uses only the standard library and the installed agrogenome_eval package.
No CDS extraction, translation, orthology, or phenotype inference is performed.
"""

import argparse
from datetime import datetime, timezone
import hashlib
from http.client import IncompleteRead
import json
import os
from pathlib import Path
import re
import socket
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen
import xml.etree.ElementTree as ET

from agrogenome_eval.sequence_analysis import analyze_fasta


BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
CANDIDATES = (
    {"species": "maize", "scientific_name": "Zea mays", "gene_id": "732788",
     "gene": "LOC732788", "transcript_accession": "NM_001305944.1",
     "protein_accession": "NP_001292873.1", "expected_refseq_status": "VALIDATED"},
    {"species": "wheat", "scientific_name": "Triticum aestivum", "gene_id": "542812",
     "gene": "LOC542812", "transcript_accession": "XM_044542390.1",
     "protein_accession": "XP_044398325.1", "expected_refseq_status": "MODEL"},
)


class NCBIClient:
    """Sequential EFetch requests, <=2 attempts/second, timeout and bounded retries."""

    def __init__(self, email: str):
        if not re.fullmatch(r"[^\s@]+@[^\s@]+\.[^\s@]+", email):
            raise ValueError("Set NCBI_EMAIL to a valid contact email before retrieval.")
        self.email = email
        self.last_request = None

    def fetch(self, accession: str, rettype: str, retmode: str, *, db: str = "nuccore") -> bytes:
        if db not in {"nuccore", "protein"}:
            raise ValueError("Only nuccore and protein retrieval is supported.")
        params = {"db": db, "id": accession, "rettype": rettype,
                  "retmode": retmode, "tool": "AgroGenomeEval", "email": self.email}
        url = BASE + "?" + urlencode(params)
        for attempt in range(3):
            if self.last_request is not None:
                time.sleep(max(0, 0.5 - (time.monotonic() - self.last_request)))
            self.last_request = time.monotonic()
            try:
                with urlopen(url, timeout=30) as response:
                    return response.read()
            except HTTPError as exc:
                if exc.code not in (429, 500, 502, 503, 504):
                    raise RuntimeError(f"NCBI HTTP {exc.code} for {accession}") from None
                retry_after = exc.headers.get("Retry-After", "") if exc.headers else ""
                if retry_after.isdigit() and int(retry_after) > 60:
                    raise RuntimeError("NCBI requested a longer wait; retry later.") from None
                delay = max(2 ** attempt, int(retry_after) if retry_after.isdigit() else 0)
            except (URLError, TimeoutError, socket.timeout, IncompleteRead):
                delay = 2 ** attempt
            if attempt < 2:
                time.sleep(delay)
        # Do not expose request URLs containing the contact email in error output.
        raise RuntimeError(f"NCBI retrieval failed after 3 attempts for {accession}") from None


def verify_metadata(xml: bytes, candidate: dict) -> dict:
    """Verify accession, organism, mRNA type, gene linkage and CDS protein metadata."""
    try:
        root = ET.fromstring(xml)
    except ET.ParseError as exc:
        raise ValueError("Malformed NCBI nucleotide metadata XML.") from exc
    records = root.findall(".//GBSeq")
    if len(records) != 1:
        raise ValueError("Expected exactly one NCBI nucleotide metadata record.")
    record = records[0]
    expected = {"GBSeq_accession-version": candidate["transcript_accession"],
                "GBSeq_organism": candidate["scientific_name"], "GBSeq_moltype": "mRNA"}
    for tag, value in expected.items():
        if record.findtext(tag) != value:
            raise ValueError(f"NCBI metadata mismatch: {tag}.")
    title = record.findtext("GBSeq_definition", "")
    if not title or candidate["gene"] not in title:
        raise ValueError("NCBI record title does not confirm the expected locus.")
    gene_confirmed = False
    protein_confirmed = False
    gene_description = None
    gene_synonyms = set()
    for feature in record.findall("GBSeq_feature-table/GBFeature"):
        qualifiers = {}
        for qualifier in feature.findall("GBFeature_quals/GBQualifier"):
            qualifiers.setdefault(qualifier.findtext("GBQualifier_name"), []).append(
                qualifier.findtext("GBQualifier_value", ""))
        if feature.findtext("GBFeature_key") == "gene":
            if (candidate["gene"] in qualifiers.get("gene", []) and
                    "GeneID:" + candidate["gene_id"] in qualifiers.get("db_xref", [])):
                gene_confirmed = True
                for names in qualifiers.get("gene_synonym", []):
                    gene_synonyms.update(name.strip() for name in names.split(";") if name.strip())
        if feature.findtext("GBFeature_key") == "CDS":
            if (candidate["protein_accession"] in qualifiers.get("protein_id", []) and
                    "GeneID:" + candidate["gene_id"] in qualifiers.get("db_xref", [])):
                protein_confirmed = True
                gene_description = next(iter(qualifiers.get("product", [])), None)
    if not gene_confirmed or not protein_confirmed:
        raise ValueError("NCBI gene/protein association mismatch.")
    comment = record.findtext("GBSeq_comment", "")
    status_match = re.search(r"\b(VALIDATED|MODEL|REVIEWED|PROVISIONAL|PREDICTED|INFERRED) REFSEQ\b", comment)
    status = status_match.group(1) if status_match else None
    if status and status != candidate["expected_refseq_status"]:
        raise ValueError("NCBI RefSeq status differs from the expected candidate status.")
    sequence = record.findtext("GBSeq_sequence", "").upper()
    if not sequence or len(sequence) != int(record.findtext("GBSeq_length", "0")):
        raise ValueError("NCBI metadata sequence length mismatch.")
    return {"record_title": title, "gene_description": gene_description,
            "gene_synonyms": sorted(gene_synonyms),
            "gene_description_source": "annotated CDS product qualifier",
            "sequence_type": "mRNA transcript", "refseq_status": status,
            "refseq_status_source": "nuccore RefSeq comment" if status else None,
            "metadata_sequence": sequence}


def verify_saved_reference(path: Path, record: dict) -> None:
    """Reject byte corruption or disagreement between a FASTA and its manifest."""
    if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
        raise ValueError("Saved FASTA checksum mismatch.")
    analysis = analyze_fasta(path)
    if analysis != record["analysis"] or analysis["record_id"] != record["transcript_accession"]:
        raise ValueError("Saved FASTA and manifest analysis mismatch.")


def fetch_dataset(output: Path, client: NCBIClient) -> dict:
    """Validate both records before publishing a new output directory.

    Refuses to overwrite existing output. Stages downloads in a temporary sibling
    directory. Unsupported biological symbols fail validation without alteration.
    """
    from tempfile import TemporaryDirectory

    output = Path(output)
    if output.exists():
        raise FileExistsError(f"Output already exists; preserving it: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    manifest = {"schema_version": 1, "tool": "AgroGenomeEval", "records": [],
                "limitations": ["Transcript references, not genomic DNA or extracted CDS.",
                                "No translation, phenotype or one-to-one orthology inference."]}
    with TemporaryDirectory(dir=output.parent, prefix=".reference-") as temporary:
        staging = Path(temporary) / "verified"
        staging.mkdir()
        for candidate in CANDIDATES:
            accession = candidate["transcript_accession"]
            metadata = verify_metadata(client.fetch(accession, "gb", "xml"), candidate)
            fasta_bytes = client.fetch(accession, "fasta", "text")
            filename = f"{candidate['species']}_{accession}.fasta"
            path = staging / filename
            path.write_bytes(fasta_bytes)
            try:
                analysis = analyze_fasta(path)
            except ValueError as exc:
                raise ValueError(f"{accession}: FASTA unsupported or malformed; sequence was not altered: {exc}") from exc
            if analysis["record_id"] != accession:
                raise ValueError("FASTA accession/version mismatch.")
            if analysis["description"] != accession + " " + metadata["record_title"]:
                raise ValueError("FASTA title differs from NCBI metadata.")
            if analysis["sequence"] != metadata.pop("metadata_sequence"):
                raise ValueError("FASTA sequence differs from NCBI metadata.")
            limitations = []
            product = metadata["gene_description"]
            if product and product.casefold() not in metadata["record_title"].casefold():
                limitations.append("Record title and annotated CDS product wording differ; retained verbatim without reconciling gene names.")
            if metadata["refseq_status"] is None:
                limitations.append("RefSeq status not found in nucleotide metadata; Gene status not independently verified.")
            manifest["records"].append({**candidate, **metadata, "filename": filename,
                "source_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
                "gene_source_url": f"https://www.ncbi.nlm.nih.gov/gene/{candidate['gene_id']}",
                "retrieval_timestamp": datetime.now(timezone.utc).isoformat(),
                "sha256": hashlib.sha256(fasta_bytes).hexdigest(),
                "validation_status": "verified" if not limitations else "verified_with_limitations",
                "limitations": limitations, "analysis": analysis})
        for record in manifest["records"]:
            verify_saved_reference(staging / record["filename"], record)
        (staging / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        # Publish all validated files together rather than leaving a partial
        # dataset if an individual file move fails. Staging is on the same disk.
        if output.exists():
            raise FileExistsError(f"Output appeared during retrieval: {output}")
        staging.rename(output)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("data/reference"))
    args = parser.parse_args()
    try:
        client = NCBIClient(os.environ.get("NCBI_EMAIL", ""))
        manifest = fetch_dataset(args.output, client)
    except (ValueError, RuntimeError, OSError) as exc:
        parser.exit(1, f"Retrieval stopped: {exc}\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
