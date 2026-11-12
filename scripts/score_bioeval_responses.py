"""Strict offline BioEval scoring; existing report directories are preserved.

Example:
 python scripts/score_bioeval_responses.py --responses data/examples/bioeval_mock_responses_v01.json --split all --output data/examples/derived/bioeval_scoring_v01
"""

import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from agrogenome_eval.bioeval_scoring import (
    POLICY_ID, SCORER_VERSION, SCORING_SCHEMA_VERSION, evaluate_responses, parse_model_response,
)
from analyze_transcript_variants import digest, encode
from generate_bioeval_benchmark import load_sources, verify_benchmark_files


PROJECT=Path(__file__).resolve().parents[1]


def verify_evaluation_files(output: Path, expected_provenance: dict | None = None) -> None:
    """Verify saved report/submission bytes and their linked provenance."""
    manifest=json.loads((output/"manifest.json").read_bytes())
    if set(manifest.get("output_checksums",{}))!={"report.json","submission.snapshot.json"}:
        raise ValueError("Incomplete scoring output checksum coverage.")
    for name,checksum in manifest["output_checksums"].items():
        if digest((output/name).read_bytes())!=checksum:raise ValueError(f"Scoring output checksum mismatch: {name}.")
    report=json.loads((output/"report.json").read_bytes());provenance=report["provenance"]
    if provenance!=manifest["provenance"] or (expected_provenance is not None and provenance!=expected_provenance):
        raise ValueError("Scoring output provenance mismatch.")
    if digest((output/"submission.snapshot.json").read_bytes())!=provenance["submitted_response_sha256"]:
        raise ValueError("Submission snapshot checksum mismatch.")


def run_evaluation(responses: Path, benchmark: Path, split: str, output: Path,
                    project_root: Path = PROJECT) -> dict:
    """Validate the full benchmark/source chain then evaluate a fixed selected split."""
    if output.exists() or output.is_symlink():raise FileExistsError(f"Preserving existing evaluation output: {output}")
    if split not in {"all","development","evaluation"}:raise ValueError("Unknown benchmark split.")
    contexts,excluded,checksums=load_sources(project_root/"data/reference",project_root/"data/derived/cds_protein",
        project_root/"data/derived/variant_effects",project_root/"data/examples/hypothetical_variants.json")
    benchmark_bytes={name:(benchmark/name).read_bytes() for name in ("manifest.json","tasks.json","gold_answers.json","splits.json")}
    verify_benchmark_files(benchmark,contexts,excluded,checksums)
    if any((benchmark/name).read_bytes()!=data for name,data in benchmark_bytes.items()):
        raise ValueError("Benchmark inputs changed during verification.")
    policy_bytes=(project_root/"docs/bioeval_scoring_v01.md").read_bytes()
    benchmark_manifest=json.loads(benchmark_bytes["manifest.json"])
    tasks=json.loads(benchmark_bytes["tasks.json"]);gold=json.loads(benchmark_bytes["gold_answers.json"])
    splits=json.loads(benchmark_bytes["splits.json"])
    raw=responses.read_bytes()
    try:submission=parse_model_response(raw.decode("utf-8"))
    except UnicodeDecodeError as exc:raise ValueError("Submission must be UTF-8 JSON.") from exc
    selected=None if split=="all" else splits[split]
    report=evaluate_responses(tasks,gold,submission,selected)
    provenance={"scoring_schema_version":SCORING_SCHEMA_VERSION,"scorer_version":SCORER_VERSION,
        "scoring_policy_id":POLICY_ID,"benchmark_version":benchmark_manifest["schema_version"],
        "benchmark_generator_version":benchmark_manifest["generator_version"],
        "benchmark_manifest_sha256":digest(benchmark_bytes["manifest.json"]),
        "gold_answers_sha256":digest(benchmark_bytes["gold_answers.json"]),
        "benchmark_file_checksums":{name:digest(data) for name,data in benchmark_bytes.items()},
        "submitted_response_sha256":digest(raw),"evaluation_split":split,
        "scoring_policy_sha256":digest(policy_bytes),
        "public_benchmark_disclosure":benchmark_manifest["public_benchmark_disclosure"],
        "timestamp_policy":"No timestamp; fixed verified inputs produce byte-deterministic output."}
    report["provenance"]=provenance
    report["evaluation_description"]=("Synthetic scorer testing; not an evaluation of a real AI model."
        if submission["submission_kind"]=="synthetic_scorer_testing" else "Offline evaluation of user-supplied structured responses; no model was called by this runner.")
    files={"report.json":encode(report),"submission.snapshot.json":raw}
    output_manifest={"scoring_schema_version":SCORING_SCHEMA_VERSION,"provenance":provenance,
                     "output_checksums":{name:digest(data) for name,data in files.items()}}
    output.parent.mkdir(parents=True,exist_ok=True)
    with TemporaryDirectory(dir=output.parent,prefix=".scoring-") as temporary:
        staging=Path(temporary)/"verified";staging.mkdir()
        for name,data in files.items():(staging/name).write_bytes(data)
        (staging/"manifest.json").write_bytes(encode(output_manifest))
        verify_evaluation_files(staging,provenance)
        # Detect changed benchmark/submission bytes before publication.
        if responses.read_bytes()!=raw or any((benchmark/name).read_bytes()!=data for name,data in benchmark_bytes.items()):
            raise ValueError("Scoring inputs changed during evaluation.")
        if (project_root/"docs/bioeval_scoring_v01.md").read_bytes()!=policy_bytes:
            raise ValueError("Scoring policy changed during evaluation.")
        if output.exists() or output.is_symlink():raise FileExistsError("Evaluation output appeared during scoring.")
        staging.rename(output)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--responses",type=Path,required=True)
    parser.add_argument("--benchmark",type=Path,default=Path("data/benchmarks/bioeval_v01"))
    parser.add_argument("--split",choices=["all","development","evaluation"],default="all")
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    try:report=run_evaluation(args.responses,args.benchmark,args.split,args.output)
    except (ValueError,OSError) as exc:parser.exit(1,f"Scoring stopped: {exc}\n")
    summary=report["summary"]
    print(json.dumps({key:summary[key] for key in ["total_selected_tasks","valid_response_count","invalid_response_count",
        "missing_response_count","abstention_count","strict_exact_match_accuracy","coverage_adjusted_mean_score"]},indent=2))


if __name__=="__main__":main()
