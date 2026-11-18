"""Score the original fourteen-task pilot from a replay-verified capture import."""

import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from agrogenome_eval.bioeval_capture import (
    canonical_json, load_verified_capture_import, load_verified_prompt_export, safe_capture_path, sha256,
)
from agrogenome_eval.bioeval_pilot import evaluate_pilot, verify_fixed_pilot
from agrogenome_eval.bioeval_scoring import POLICY_ID, SCORER_VERSION, SCORING_SCHEMA_VERSION
from generate_bioeval_benchmark import load_sources, verify_benchmark_files


PROJECT = Path(__file__).resolve().parents[1]
OUTPUT_FILES = {
    "report.json", "experiment.snapshot.json", "capture_validation.snapshot.json", "submission.snapshot.json",
    "pilot_export_manifest.snapshot.json", "prompts.snapshot.json", "capture_import_manifest.snapshot.json",
}


def verify_pilot_outputs(output: Path, expected_provenance: dict | None = None) -> None:
    """Verify all saved pilot outputs and the fixed selected population."""
    manifest = json.loads(safe_capture_path(output, "manifest.json").read_bytes())
    if set(manifest.get("output_checksums", {})) != OUTPUT_FILES:
        raise ValueError("Incomplete pilot output checksum coverage.")
    for name, checksum in manifest["output_checksums"].items():
        if sha256(safe_capture_path(output, name).read_bytes()) != checksum:
            raise ValueError(f"Pilot output checksum mismatch for {name}")
    report = json.loads((output / "report.json").read_bytes())
    if report["provenance"] != manifest["provenance"]:
        raise ValueError("Pilot report provenance mismatch.")
    if expected_provenance is not None and report["provenance"] != expected_provenance:
        raise ValueError("Pilot input provenance mismatch.")
    if report["summary"]["total_selected_tasks"] != 14 or sum(report["pilot_outcome_counts"].values()) != 14:
        raise ValueError("Pilot output denominator must be fourteen.")
    if report["population"]["selected_task_ids"] != manifest["selected_task_ids"]:
        raise ValueError("Pilot output task-ID coverage mismatch.")


def run_pilot(prompt_export: Path, capture_import: Path, output: Path, *,
              benchmark: Path = PROJECT / "data/benchmarks/bioeval_v01",
              project_root: Path = PROJECT) -> dict:
    """Verify source, export and captured evidence before fixed-population scoring."""
    prompt_export, capture_import, output = Path(prompt_export), Path(capture_import), Path(output)
    benchmark, project_root = Path(benchmark), Path(project_root)
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Preserving existing pilot output {output}")
    if any(output.resolve().is_relative_to(root.resolve()) for root in (prompt_export, capture_import, benchmark)):
        raise ValueError("Pilot output must be separate from source evidence directories.")
    contexts, excluded, source_hashes = load_sources(
        project_root / "data/reference", project_root / "data/derived/cds_protein",
        project_root / "data/derived/variant_effects",
        project_root / "data/examples/hypothetical_variants.json",
    )
    prompts = load_verified_prompt_export(prompt_export, benchmark)
    verify_fixed_pilot(prompts)
    verify_benchmark_files(benchmark, contexts, excluded, source_hashes)
    validation, capture_bytes = load_verified_capture_import(capture_import, prompts)
    tasks = json.loads(prompts.benchmark_snapshots["tasks.json"])
    gold = json.loads(prompts.benchmark_snapshots["gold_answers.json"])
    policy_path = project_root / "docs/bioeval_scoring_v01.md"
    policy_bytes = policy_path.read_bytes()
    report = evaluate_pilot(prompts, validation, tasks, gold)
    input_files = {
        **{"benchmark/" + name: sha256(raw) for name, raw in prompts.benchmark_snapshots.items()},
        **{"pilot_export/" + name: sha256(raw) for name, raw in prompts.snapshots.items()},
        **{"capture_import/" + name: sha256(raw) for name, raw in capture_bytes.items()},
        "scoring_policy": sha256(policy_bytes),
    }
    provenance = {
        "pilot_evaluation_schema_version": report["pilot_evaluation_schema_version"],
        "pilot_evaluator_version": report["pilot_evaluator_version"],
        "scoring_schema_version": SCORING_SCHEMA_VERSION, "scorer_version": SCORER_VERSION,
        "scoring_policy_id": POLICY_ID, "scoring_policy_sha256": sha256(policy_bytes),
        "pilot_export_identity": prompts.identity, "input_checksums": dict(sorted(input_files.items())),
        "scientific_source_checksums": source_hashes,
        "timestamp_policy": "Preserve recorded execution times; generate no evaluation timestamp.",
    }
    report["provenance"] = provenance
    files = {
        "report.json": canonical_json(report),
        "experiment.snapshot.json": validation.snapshots["experiment.json"],
        "capture_validation.snapshot.json": capture_bytes["capture_validation.json"],
        "submission.snapshot.json": capture_bytes["submission.json"],
        "pilot_export_manifest.snapshot.json": prompts.snapshots["manifest.json"],
        "prompts.snapshot.json": prompts.snapshots["prompts.json"],
        "capture_import_manifest.snapshot.json": capture_bytes["manifest.json"],
    }
    manifest = {"pilot_evaluation_schema_version": report["pilot_evaluation_schema_version"],
                "provenance": provenance, "selected_task_ids": sorted(prompts.records),
                "output_checksums": {name: sha256(raw) for name, raw in sorted(files.items())}}
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent, prefix=".pilot-scoring-") as temporary:
        staging = Path(temporary) / "verified"
        staging.mkdir()
        for name, raw in files.items():
            (staging / name).write_bytes(raw)
        (staging / "manifest.json").write_bytes(canonical_json(manifest))
        verify_pilot_outputs(staging, provenance)
        for root, snapshots in [(prompt_export, prompts.snapshots), (capture_import, capture_bytes),
                                (benchmark, prompts.benchmark_snapshots)]:
            if any(safe_capture_path(root, name).read_bytes() != raw for name, raw in snapshots.items()):
                raise ValueError("Pilot input evidence changed during scoring.")
        if policy_path.read_bytes() != policy_bytes:
            raise ValueError("Scoring policy changed during pilot evaluation.")
        _, _, fresh_hashes = load_sources(
            project_root / "data/reference", project_root / "data/derived/cds_protein",
            project_root / "data/derived/variant_effects",
            project_root / "data/examples/hypothetical_variants.json",
        )
        if fresh_hashes != source_hashes:
            raise ValueError("Scientific source evidence changed during pilot evaluation.")
        if output.exists() or output.is_symlink():
            raise FileExistsError("Pilot output appeared during publication.")
        staging.rename(output)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prompt-export", type=Path, required=True)
    parser.add_argument("--capture-import", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--benchmark", type=Path, default=PROJECT / "data/benchmarks/bioeval_v01")
    args = parser.parse_args()
    try:
        report = run_pilot(args.prompt_export, args.capture_import, args.output, benchmark=args.benchmark)
    except (ValueError, OSError, KeyError, TypeError) as exc:
        parser.exit(1, f"Pilot scoring stopped: {exc}\n")
    fields = ["total_selected_tasks", "valid_response_count", "invalid_response_count", "abstention_count",
              "missing_response_count", "strict_exact_match_accuracy", "coverage_adjusted_mean_score"]
    print(json.dumps({"submission_kind": report["submission_kind"],
                      "summary": {name: report["summary"][name] for name in fields}}, indent=2))


if __name__ == "__main__":
    main()
