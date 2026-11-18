"""Validate locally captured observations and publish a scoring submission offline."""

import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from agrogenome_eval.bioeval_capture import (
    CAPTURE_SCHEMA_VERSION, CAPTURE_VERSION, canonical_json, convert_to_scoring_submission,
    load_verified_prompt_export, load_verified_capture_import, safe_capture_path, sha256, validate_capture,
)
from generate_bioeval_benchmark import load_sources, verify_benchmark_files


PROJECT = Path(__file__).resolve().parents[1]
CORE_FILES = {
    "experiment.snapshot.json", "records.snapshot.json", "prompt_export_manifest.snapshot.json",
    "prompts.snapshot.json", "capture_validation.json", "submission.json",
}


class CaptureImportError(ValueError):
    """Invalid capture provenance with a deterministic, machine-readable report."""
    def __init__(self, report: dict):
        super().__init__("Capture provenance validation failed.")
        self.report = report


def verify_import_outputs(output: Path, expected_manifest: dict | None = None, *, prompts=None) -> None:
    """Verify exact output bytes, required snapshots and selected-ID provenance."""
    manifest = json.loads(safe_capture_path(output, "manifest.json").read_bytes())
    if expected_manifest is not None and manifest != expected_manifest:
        raise ValueError("Capture import manifest mismatch.")
    checksums = manifest.get("output_checksums", {})
    if not isinstance(checksums, dict) or not CORE_FILES <= checksums.keys():
        raise ValueError("Incomplete capture-import checksum coverage.")
    if set(checksums) - CORE_FILES != {"evidence/" + name for name in manifest["evidence_checksums"]}:
        raise ValueError("Capture evidence checksum coverage mismatch.")
    for name, checksum in checksums.items():
        if sha256(safe_capture_path(output, name).read_bytes()) != checksum:
            raise ValueError(f"Capture-import checksum mismatch for {name}")
    report = json.loads((output / "capture_validation.json").read_bytes())
    submission = json.loads((output / "submission.json").read_bytes())
    if (not report["valid"] or report["selected_task_ids"] != manifest["selected_task_ids"]
            or report["selected_task_count"] != manifest["selected_task_count"]
            or submission["submission_kind"] != manifest["submission_kind"]):
        raise ValueError("Capture import population or labeling mismatch.")
    for name, checksum in manifest["evidence_checksums"].items():
        if checksums["evidence/" + name] != checksum:
            raise ValueError("Preserved source evidence checksum mismatch.")
    if prompts is None:
        prompts = load_verified_prompt_export(output, PROJECT / "data/benchmarks/bioeval_v01", imported=True)
    load_verified_capture_import(output, prompts)


def run_import(prompt_export: Path, capture: Path, output: Path, *,
               benchmark: Path = PROJECT / "data/benchmarks/bioeval_v01",
               project_root: Path = PROJECT) -> dict:
    """Verify scientific provenance, capture observations, and publish complete outputs."""
    prompt_export, capture, output = Path(prompt_export), Path(capture), Path(output)
    benchmark, project_root = Path(benchmark), Path(project_root)
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Preserving existing output path {output}")
    if any(output.resolve().is_relative_to(root.resolve()) for root in (prompt_export, capture, benchmark)):
        raise ValueError("Derived output must be separate from source capture, export and benchmark directories.")
    contexts, excluded, source_checksums = load_sources(
        project_root / "data/reference", project_root / "data/derived/cds_protein",
        project_root / "data/derived/variant_effects",
        project_root / "data/examples/hypothetical_variants.json",
    )
    prompts = load_verified_prompt_export(prompt_export, benchmark)
    verify_benchmark_files(benchmark, contexts, excluded, source_checksums)
    validation = validate_capture(capture, prompts)
    if not validation.report["valid"]:
        raise CaptureImportError(validation.report)
    submission = convert_to_scoring_submission(validation)
    files = {
        "experiment.snapshot.json": validation.snapshots["experiment.json"],
        "records.snapshot.json": validation.snapshots["records.json"],
        "prompt_export_manifest.snapshot.json": prompts.snapshots["manifest.json"],
        "prompts.snapshot.json": prompts.snapshots["prompts.json"],
        "capture_validation.json": canonical_json(validation.report),
        "submission.json": canonical_json(submission),
    }
    files.update({"evidence/" + name: raw for name, raw in validation.evidence.items()})
    manifest = {
        "capture_schema_version": CAPTURE_SCHEMA_VERSION,
        "capture_version": CAPTURE_VERSION,
        "experiment_id": validation.experiment["experiment_id"],
        "submission_kind": submission["submission_kind"],
        "prompt_export_identity": prompts.identity,
        "selected_task_ids": sorted(prompts.records),
        "selected_task_count": len(prompts.records),
        "partition": prompts.manifest["partition"],
        "stage": prompts.manifest["stage"],
        "benchmark_file_checksums": {name: sha256(raw) for name, raw in prompts.benchmark_snapshots.items()},
        "capture_metadata_checksums": {name: sha256(raw) for name, raw in validation.snapshots.items()},
        "evidence_checksums": {name: sha256(raw) for name, raw in sorted(validation.evidence.items())},
        "output_checksums": {name: sha256(raw) for name, raw in sorted(files.items())},
        "provenance_confidence": validation.report["provenance_confidence"],
        "scoring_population_requirement": "Use exactly selected_task_ids for fixed-denominator scoring.",
        "timestamp_policy": "Preserve supplied execution timestamps; derived output adds no timestamp.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent, prefix=".capture-import-") as temporary:
        staging = Path(temporary) / "verified"
        staging.mkdir()
        for name, raw in files.items():
            path = staging / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(raw)
        (staging / "manifest.json").write_bytes(canonical_json(manifest))
        verify_import_outputs(staging, manifest)
        for name, raw in {**validation.snapshots, **validation.evidence}.items():
            if safe_capture_path(capture, name).read_bytes() != raw:
                raise ValueError("Source capture changed during import.")
        for root, snapshots in [(prompt_export, prompts.snapshots), (benchmark, prompts.benchmark_snapshots)]:
            if any(safe_capture_path(root, name).read_bytes() != raw for name, raw in snapshots.items()):
                raise ValueError("Prompt export or benchmark changed during import.")
        _, _, fresh_checksums = load_sources(
            project_root / "data/reference", project_root / "data/derived/cds_protein",
            project_root / "data/derived/variant_effects",
            project_root / "data/examples/hypothetical_variants.json",
        )
        if fresh_checksums != source_checksums:
            raise ValueError("Scientific source evidence changed during import.")
        if output.exists() or output.is_symlink():
            raise FileExistsError("Output path appeared during publication.")
        staging.rename(output)
    return validation.report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prompt-export", type=Path, required=True)
    parser.add_argument("--capture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--benchmark", type=Path, default=PROJECT / "data/benchmarks/bioeval_v01")
    args = parser.parse_args()
    try:
        report = run_import(args.prompt_export, args.capture, args.output, benchmark=args.benchmark)
    except CaptureImportError as exc:
        parser.exit(1, canonical_json(exc.report).decode("utf-8"))
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(1, f"Capture import stopped: {exc}\n")
    print(json.dumps({"valid": report["valid"], "selected_task_count": report["selected_task_count"],
                      "selected_response_count": len(report["selected_attempts"]),
                      "missing_task_count": len(report["missing_task_ids"])}, indent=2))


if __name__ == "__main__":
    main()
