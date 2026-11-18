"""Export verified BioEval questions offline for manual or future model testing."""

import argparse
from collections import defaultdict
import json
from pathlib import Path
import random
from tempfile import TemporaryDirectory

from agrogenome_eval.bioeval_benchmark import TASK_INPUT_FIELDS, answer_category
from analyze_transcript_variants import digest, encode
from generate_bioeval_benchmark import load_sources, verify_benchmark_files


PROJECT = Path(__file__).resolve().parents[1]
PROTOCOL_ID = "bioeval-empirical-protocol-v0.1"
EXPORT_SCHEMA_VERSION = "0.1"
EXPORTER_VERSION = "0.1.0"
DEFAULT_PILOT_SEED = 20261112
BENCHMARK_FILES = ("manifest.json", "tasks.json", "gold_answers.json", "splits.json")


def select_tasks(tasks: dict, gold: dict, splits: dict, partition: str,
                 pilot: bool = False, seed: int = DEFAULT_PILOT_SEED) -> tuple[list[dict], list[dict]]:
    """Select from an already verified benchmark, sorting records by original ID.

    Pilot strata use gold only inside the researcher-side selection process.
    The returned records retain exactly the original model-facing task schema.
    """
    if partition not in {"development", "evaluation"}:
        raise ValueError("Partition must be development or evaluation.")
    if type(pilot) is not bool or type(seed) is not int:
        raise ValueError("Pilot must be boolean and the seed must be an integer.")
    if pilot and partition != "development":
        raise ValueError("Pilot selection requires the development partition.")
    task_map = {task["task_id"]: task for task in tasks["tasks"]}
    selected_ids = sorted(splits[partition])
    strata = []
    if pilot:
        candidates = defaultdict(list)
        for task_id in selected_ids:
            task = task_map[task_id]
            key = (task["input"]["scientific_name"], answer_category(gold["answers"][task_id]))
            candidates[key].append(task_id)
        if len(candidates) != 14:
            raise ValueError("Pilot requires all 14 verified species-by-category strata.")
        generator = random.Random(seed)
        selected_ids = []
        for (species, category), task_ids in sorted(candidates.items()):
            task_id = generator.choice(sorted(task_ids))
            selected_ids.append(task_id)
            strata.append({"scientific_name": species, "category": category,
                           "task_id": task_id, "eligible_development_tasks": len(task_ids)})
    return [task_map[task_id] for task_id in sorted(selected_ids)], strata


def verify_prompt_export(output: Path, expected_manifest: dict | None = None) -> None:
    """Verify exported payload bytes, exact task shape, IDs and question hashes."""
    manifest = json.loads((output / "manifest.json").read_bytes())
    if expected_manifest is not None and manifest != expected_manifest:
        raise ValueError("Prompt-export manifest differs from verified provenance.")
    if (manifest.get("export_schema_version") != EXPORT_SCHEMA_VERSION
            or manifest.get("protocol_id") != PROTOCOL_ID
            or set(manifest.get("output_checksums", {})) != {"prompts.json"}):
        raise ValueError("Invalid prompt-export manifest or checksum coverage.")
    raw = (output / "prompts.json").read_bytes()
    if digest(raw) != manifest["output_checksums"]["prompts.json"]:
        raise ValueError("Prompt-export checksum mismatch.")
    payload = json.loads(raw)
    if set(payload) != {"schema_version", "tasks"} or payload["schema_version"] != "0.1":
        raise ValueError("Export must retain the existing benchmark task envelope.")
    records = payload["tasks"]
    for record in records:
        if (set(record) != {"task_id", "input", "question"}
                or set(record["input"]) != TASK_INPUT_FIELDS
                or not isinstance(record["task_id"], str)
                or not isinstance(record["question"], str)):
            raise ValueError("Invalid model-facing record shape.")
    ids = [record["task_id"] for record in records]
    if (ids != sorted(set(ids)) or ids != manifest["selected_task_ids"]
            or len(ids) != manifest["task_count"]):
        raise ValueError("Prompt-export task coverage or ordering mismatch.")
    hashes = {record["task_id"]: digest(record["question"].encode("utf-8")) for record in records}
    if hashes != manifest["question_utf8_sha256"]:
        raise ValueError("Exported question checksum mismatch.")


def export_prompts(benchmark: Path, partition: str, output: Path, *, pilot: bool = False,
                   seed: int = DEFAULT_PILOT_SEED, project_root: Path = PROJECT) -> dict:
    """Verify sources and benchmark, then atomically publish a question-only export."""
    benchmark, output, project_root = Path(benchmark), Path(output), Path(project_root)
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Preserving existing export path {output}")
    if partition not in {"development", "evaluation"}:
        raise ValueError("Partition must be development or evaluation.")
    if type(pilot) is not bool or type(seed) is not int:
        raise ValueError("Pilot must be boolean and the seed must be an integer.")
    if pilot and partition != "development":
        raise ValueError("Pilot selection requires the development partition.")
    contexts, excluded, source_checksums = load_sources(
        project_root / "data/reference", project_root / "data/derived/cds_protein",
        project_root / "data/derived/variant_effects",
        project_root / "data/examples/hypothetical_variants.json",
    )
    snapshots = {name: (benchmark / name).read_bytes() for name in BENCHMARK_FILES}
    verify_benchmark_files(benchmark, contexts, excluded, source_checksums)
    if any((benchmark / name).read_bytes() != raw for name, raw in snapshots.items()):
        raise ValueError("Benchmark changed during verification.")
    tasks = json.loads(snapshots["tasks.json"])
    gold = json.loads(snapshots["gold_answers.json"])
    splits = json.loads(snapshots["splits.json"])
    benchmark_manifest = json.loads(snapshots["manifest.json"])
    records, strata = select_tasks(tasks, gold, splits, partition, pilot, seed)
    prompts = encode({"schema_version": tasks["schema_version"], "tasks": records})
    manifest = {
        "export_schema_version": EXPORT_SCHEMA_VERSION,
        "exporter_version": EXPORTER_VERSION,
        "protocol_id": PROTOCOL_ID,
        "artifact_role": "researcher_provenance_only",
        "stage": "pilot" if pilot else partition,
        "partition": partition,
        "pilot_sampling_seed": seed if pilot else None,
        "selection_method": "one random choice per sorted species/category stratum" if pilot else "entire requested partition",
        "batch_order": "ascending original opaque task ID",
        "task_count": len(records),
        "selected_task_ids": [record["task_id"] for record in records],
        "pilot_strata": strata,
        "question_utf8_sha256": {
            record["task_id"]: digest(record["question"].encode("utf-8")) for record in records
        },
        "benchmark_schema_version": benchmark_manifest["schema_version"],
        "benchmark_generator_version": benchmark_manifest["generator_version"],
        "benchmark_generation_seed": benchmark_manifest["generation_seed"],
        "benchmark_file_checksums": {name: digest(raw) for name, raw in snapshots.items()},
        "source_checksums": source_checksums,
        "public_benchmark_disclosure": benchmark_manifest["public_benchmark_disclosure"],
        "output_checksums": {"prompts.json": digest(prompts)},
        "timestamp_policy": "Export bytes omit timestamps. Actual model execution times belong in experiment records.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent, prefix=".prompt-export-") as temporary:
        staging = Path(temporary) / "verified"
        staging.mkdir()
        (staging / "prompts.json").write_bytes(prompts)
        (staging / "manifest.json").write_bytes(encode(manifest))
        verify_prompt_export(staging, manifest)
        if any((benchmark / name).read_bytes() != raw for name, raw in snapshots.items()):
            raise ValueError("Benchmark changed during export.")
        if output.exists() or output.is_symlink():
            raise FileExistsError("Export path appeared during publication.")
        staging.rename(output)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, default=PROJECT / "data/benchmarks/bioeval_v01")
    parser.add_argument("--partition", choices=["development", "evaluation"], required=True)
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--seed", type=int, default=DEFAULT_PILOT_SEED)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        manifest = export_prompts(args.benchmark, args.partition, args.output,
                                  pilot=args.pilot, seed=args.seed)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"Prompt export stopped: {exc}\n")
    print(json.dumps({"stage": manifest["stage"], "partition": manifest["partition"],
                      "task_count": manifest["task_count"], "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
