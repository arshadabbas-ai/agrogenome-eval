"""Offline capture provenance, byte preservation and scorer submission conversion."""

from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path, PureWindowsPath
import random
import re

from .bioeval_benchmark import answer_category
from .bioeval_scoring import ResponseFormatError, parse_model_response, validate_model_response


CAPTURE_SCHEMA_VERSION = "0.1"
CAPTURE_VERSION = "0.1.0"
ELIGIBLE_OUTCOMES = {"completed", "refused", "truncated", "abstained"}
OUTCOMES = ELIGIBLE_OUTCOMES | {"no_response", "execution_failure"}
BENCHMARK_FILES = ("manifest.json", "tasks.json", "gold_answers.json", "splits.json")
IMPORT_CORE_FILES = {
    "experiment.snapshot.json", "records.snapshot.json", "prompt_export_manifest.snapshot.json",
    "prompts.snapshot.json", "capture_validation.json", "submission.json",
}
EXPERIMENT_FIELDS = {
    "capture_schema_version", "experiment_id", "submission_kind", "benchmark_schema_version",
    "prompt_export", "model", "access_method", "execution_started_at", "instructions",
    "generation_parameters", "tool_policy", "retry_policy", "provenance",
}
RECORD_FIELDS = {
    "task_id", "attempt", "retry_of", "retry_reason", "timestamp", "outcome",
    "prompt", "response", "output_error", "tool_use", "outcome_details",
}


def sha256(raw: bytes) -> str:
    """Hash exact saved bytes, without normalization."""
    return hashlib.sha256(raw).hexdigest()


def canonical_json(value: object) -> bytes:
    """Serialize deterministic UTF-8 JSON with sorted keys and LF termination."""
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def _json_object(raw: bytes) -> dict:
    return parse_model_response(raw.decode("utf-8"))


def _text(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _timestamp(value: object) -> datetime:
    if not isinstance(value, str) or not re.fullmatch(
        r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?(?:Z|[+-]\d{2}:\d{2})", value
    ):
        raise ValueError("Timestamp requires ISO 8601 seconds and an explicit UTC offset.")
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.utcoffset() is None:
        raise ValueError("Timestamp requires an explicit UTC offset.")
    return parsed


def safe_capture_path(root: Path, name: str) -> Path:
    """Resolve a portable relative file path and reject traversal or symlink escapes."""
    if (not _text(name) or "\\" in name or ":" in name or PureWindowsPath(name).drive
            or Path(name).is_absolute() or ".." in name.split("/")):
        raise ValueError("Unsafe evidence path.")
    try:
        resolved_root = Path(root).resolve(strict=True)
        resolved = (resolved_root / name).resolve(strict=True)
    except RuntimeError as exc:
        raise ValueError("Cyclic or invalid symlink resolution.") from exc
    if not resolved.is_relative_to(resolved_root) or not resolved.is_file():
        raise ValueError("Evidence path escapes the capture root or is not a file.")
    return resolved


def _read_evidence(root: Path, descriptor: dict, evidence: dict[str, bytes]) -> bytes:
    if (not isinstance(descriptor, dict) or set(descriptor) != {"path", "sha256"}
            or not isinstance(descriptor["sha256"], str)
            or not re.fullmatch(r"[0-9a-f]{64}", descriptor["sha256"])):
        raise ValueError("Evidence needs a relative path and lowercase SHA-256.")
    name = descriptor["path"]
    raw = safe_capture_path(root, name).read_bytes()
    if sha256(raw) != descriptor["sha256"]:
        raise ValueError("Evidence checksum mismatch.")
    if name in evidence and evidence[name] != raw:
        raise ValueError("Evidence changed during validation.")
    evidence[name] = raw
    return raw


@dataclass
class VerifiedPromptExport:
    """Verified original records and immutable byte snapshots for this validation."""
    manifest: dict
    records: dict[str, dict]
    benchmark_ids: set[str]
    partition_ids: set[str]
    snapshots: dict[str, bytes]
    benchmark_snapshots: dict[str, bytes]

    @property
    def identity(self) -> dict:
        """Capture binding that includes the exact selected-ID population."""
        return {
            "manifest_sha256": sha256(self.snapshots["manifest.json"]),
            "protocol_id": self.manifest["protocol_id"],
            "partition": self.manifest["partition"],
            "stage": self.manifest["stage"],
            "selection_id": sha256(canonical_json(sorted(self.records))),
        }


def load_verified_prompt_export(export_root: Path, benchmark: Path, *, imported: bool = False) -> VerifiedPromptExport:
    """Bind export bytes to a trusted benchmark, replaying original pilot selection.

    The importer additionally verifies the benchmark's full scientific source
    chain before calling this reusable byte/schema/selection verifier.
    """
    names = {"manifest.json": "prompt_export_manifest.snapshot.json", "prompts.json": "prompts.snapshot.json"} if imported else {
        "manifest.json": "manifest.json", "prompts.json": "prompts.json",
    }
    snapshots = {name: safe_capture_path(export_root, filename).read_bytes()
                 for name, filename in names.items()}
    manifest = _json_object(snapshots["manifest.json"])
    payload = _json_object(snapshots["prompts.json"])
    benchmark_raw = {name: safe_capture_path(benchmark, name).read_bytes() for name in BENCHMARK_FILES}
    benchmark_manifest = _json_object(benchmark_raw["manifest.json"])
    if manifest.get("benchmark_file_checksums") != {n: sha256(b) for n, b in benchmark_raw.items()}:
        raise ValueError("Prompt export benchmark checksum binding mismatch.")
    if manifest.get("source_checksums") != benchmark_manifest["source_checksums"]:
        raise ValueError("Prompt export source provenance mismatch.")
    for name in BENCHMARK_FILES[1:]:
        if sha256(benchmark_raw[name]) != benchmark_manifest["output_checksums"][name]:
            raise ValueError("Benchmark payload checksum mismatch.")
    if (manifest.get("export_schema_version") != "0.1"
            or manifest.get("exporter_version") != "0.1.0"
            or manifest.get("protocol_id") != "bioeval-empirical-protocol-v0.1"
            or manifest.get("artifact_role") != "researcher_provenance_only"
            or manifest.get("benchmark_schema_version") != benchmark_manifest["schema_version"]
            or manifest.get("benchmark_generator_version") != benchmark_manifest["generator_version"]
            or manifest.get("benchmark_generation_seed") != benchmark_manifest["generation_seed"]
            or manifest.get("public_benchmark_disclosure") != benchmark_manifest["public_benchmark_disclosure"]
            or manifest.get("output_checksums") != {"prompts.json": sha256(snapshots["prompts.json"])}):
        raise ValueError("Prompt export metadata or payload checksum mismatch.")
    tasks = _json_object(benchmark_raw["tasks.json"])
    gold = _json_object(benchmark_raw["gold_answers.json"])
    splits = _json_object(benchmark_raw["splits.json"])
    task_map = {t["task_id"]: t for t in tasks["tasks"]}
    partition = manifest.get("partition")
    stage = manifest.get("stage")
    if partition not in {"development", "evaluation"} or stage not in {"pilot", partition}:
        raise ValueError("Invalid prompt export partition or stage.")
    ids = sorted(splits[partition])
    strata = []
    if stage == "pilot":
        if partition != "development" or type(manifest.get("pilot_sampling_seed")) is not int:
            raise ValueError("Pilot requires development and an integer sampling seed.")
        groups = defaultdict(list)
        for task_id in ids:
            groups[(task_map[task_id]["input"]["scientific_name"], answer_category(gold["answers"][task_id]))].append(task_id)
        if len(groups) != 14:
            raise ValueError("Pilot requires fourteen species/category strata.")
        generator = random.Random(manifest["pilot_sampling_seed"])
        ids = []
        for (species, category), candidates in sorted(groups.items()):
            task_id = generator.choice(sorted(candidates))
            ids.append(task_id)
            strata.append({"scientific_name": species, "category": category, "task_id": task_id,
                           "eligible_development_tasks": len(candidates)})
        ids.sort()
    elif manifest.get("pilot_sampling_seed") is not None:
        raise ValueError("Full-partition exports require a null pilot seed.")
    expected = {"schema_version": "0.1", "tasks": [task_map[i] for i in ids]}
    hashes = {i: sha256(task_map[i]["question"].encode("utf-8")) for i in ids}
    if (payload != expected or manifest.get("selected_task_ids") != ids
            or manifest.get("task_count") != len(ids) or manifest.get("pilot_strata") != strata
            or manifest.get("question_utf8_sha256") != hashes
            or manifest.get("batch_order") != "ascending original opaque task ID"):
        raise ValueError("Prompt records or selected population differ from the original export protocol.")
    return VerifiedPromptExport(manifest, {i: task_map[i] for i in ids}, set(task_map),
                                set(splits[partition]), snapshots, benchmark_raw)


def _issue(code: str, field: str, message: str, severity: str = "ERROR") -> dict:
    return {"code": code, "field": field, "message": message, "severity": severity}


def _keys(value: object, required: set[str], field_name: str) -> None:
    if not isinstance(value, dict) or set(value) != required:
        raise ValueError(f"{field_name} must contain exactly its documented fields.")


def validate_experiment_metadata(experiment: dict, prompts: VerifiedPromptExport) -> list[dict]:
    """Validate required identity/settings, without authenticating model origin."""
    issues = []
    try:
        _keys(experiment, EXPERIMENT_FIELDS, "experiment")
        if experiment["capture_schema_version"] != CAPTURE_SCHEMA_VERSION:
            raise ValueError("Unsupported capture schema version.")
        if not _text(experiment["experiment_id"]):
            raise ValueError("Experiment ID is required.")
        if experiment["submission_kind"] not in ("synthetic_scorer_testing", "model_responses"):
            raise ValueError("Invalid submission kind.")
        if experiment["benchmark_schema_version"] != prompts.manifest["benchmark_schema_version"]:
            raise ValueError("Benchmark version mismatch.")
        if experiment["prompt_export"] != prompts.identity:
            raise ValueError("Experiment prompt-export identity mismatch.")
        _timestamp(experiment["execution_started_at"])
        if experiment["access_method"] not in ("manual_interface", "api", "local_inference"):
            raise ValueError("Invalid execution access method.")
        model = experiment["model"]
        _keys(model, {"source", "name", "version"}, "model")
        for name in ("source", "name"):
            if not _text(model[name]) or model[name].strip().lower() in {"unknown", "n/a", "none", "placeholder"}:
                raise ValueError("Meaningful model source and exact visible name are required.")
        version = model["version"]
        _keys(version, {"availability", "value", "reason"}, "model.version")
        if version["availability"] == "recorded":
            if not _text(version["value"]) or version["reason"] is not None:
                raise ValueError("Recorded model version requires its value and a null reason.")
        elif version["availability"] != "unavailable" or version["value"] is not None or not _text(version["reason"]):
            raise ValueError("Unavailable model version requires null and a reason.")
        instructions = experiment["instructions"]
        _keys(instructions, {"system", "developer"}, "instructions")
        for instruction in instructions.values():
            _keys(instruction, {"availability", "text", "reason"}, "instruction")
            availability = instruction["availability"]
            if availability == "recorded":
                if not _text(instruction["text"]) or instruction["reason"] is not None:
                    raise ValueError("Recorded instructions require exact text.")
            elif availability == "absent":
                if instruction["text"] is not None or instruction["reason"] is not None:
                    raise ValueError("Absent instructions require null text and reason.")
            elif availability != "unavailable" or instruction["text"] is not None or not _text(instruction["reason"]):
                raise ValueError("Unavailable instructions require a reason.")
        parameters = experiment["generation_parameters"]
        if not isinstance(parameters, dict) or not {"temperature", "seed"} <= parameters.keys():
            raise ValueError("Temperature and seed availability must be documented.")
        for name, setting in parameters.items():
            if not isinstance(name, str) or not name:
                raise ValueError("Generation parameter names must be nonempty strings.")
            if name.lower() in {"api_key", "access_token", "authorization", "password", "secret", "token"}:
                raise ValueError("Credentials are excluded from generation parameters.")
            _keys(setting, {"availability", "value", "reason"}, "generation parameter")
            value = setting["value"]
            if setting["availability"] == "unavailable":
                if value is not None or not _text(setting["reason"]):
                    raise ValueError("Unavailable parameters require null and an availability reason.")
                continue
            if setting["availability"] not in ("set", "observed_default") or value is None or setting["reason"] is not None:
                raise ValueError("Known parameters require an explicit value and availability status.")
            if not isinstance(value, (str, bool, int, float)) or isinstance(value, float) and not math.isfinite(value):
                raise ValueError("Generation parameters must be finite JSON scalar values.")
            if name == "temperature" and (type(value) not in (int, float) or value < 0):
                raise ValueError("Temperature requires a nonnegative number.")
            if name == "seed" and (type(value) is not int or value < 0):
                raise ValueError("Seed requires a nonnegative integer.")
        policy = experiment["tool_policy"]
        _keys(policy, {"mode", "details"}, "tool policy")
        if policy["mode"] not in ("disabled", "permitted", "unavailable") or not _text(policy["details"]):
            raise ValueError("Tool-use policy and details are required.")
        retry = experiment["retry_policy"]
        _keys(retry, {"selection", "max_attempts"}, "retry policy")
        if (retry["selection"] != "first_completed_eligible"
                or type(retry["max_attempts"]) is not int or retry["max_attempts"] < 1):
            raise ValueError("Retry policy requires first_completed_eligible and a positive attempt bound.")
        provenance = experiment["provenance"]
        _keys(provenance, {"basis", "attestation", "capture_limitations"}, "provenance")
        if not _text(provenance["attestation"]) or not _text(provenance["capture_limitations"]):
            raise ValueError("Origin attestation and capture limitations are required.")
        if experiment["submission_kind"] == "synthetic_scorer_testing":
            if provenance["basis"] != "synthetic_fixture":
                raise ValueError("Synthetic capture requires the synthetic_fixture provenance basis.")
        elif provenance["basis"] not in ("researcher_attestation", "provider_payload", "local_execution_record"):
            raise ValueError("Model responses require documented execution-source provenance.")
    except (ValueError, TypeError, KeyError) as exc:
        issues.append(_issue("INVALID_EXPERIMENT_METADATA", "experiment", str(exc)))
    return issues


def verify_prompt_identity(record: dict, prompts: VerifiedPromptExport, root: Path,
                           evidence: dict[str, bytes]) -> None:
    """Verify the exact original record identity and the supplied user-prompt bytes."""
    descriptor = record["prompt"]
    _keys(descriptor, {"path", "sha256", "record_sha256"}, "prompt")
    expected = prompts.records[record["task_id"]]
    if descriptor["record_sha256"] != sha256(canonical_json(expected)):
        raise ValueError("Prompt record identity mismatch.")
    raw = _read_evidence(root, {"path": descriptor["path"], "sha256": descriptor["sha256"]}, evidence)
    if raw != expected["question"].encode("utf-8"):
        raise ValueError("Prompt bytes differ from the exported question.")


def verify_response_bytes(response: dict, root: Path, evidence: dict[str, bytes]) -> bytes | None:
    """Read and verify raw answer bytes independently of any native payload."""
    _keys(response, {"path", "sha256", "source_payload"}, "response")
    raw = None
    if response["path"] is not None or response["sha256"] is not None:
        raw = _read_evidence(root, {"path": response["path"], "sha256": response["sha256"]}, evidence)
        raw.decode("utf-8")
    payload = response["source_payload"]
    if payload is not None:
        _keys(payload, {"path", "sha256", "extraction"}, "source payload")
        if not _text(payload["extraction"]):
            raise ValueError("Source-payload extraction or relationship must be documented.")
        _read_evidence(root, {"path": payload["path"], "sha256": payload["sha256"]}, evidence)
    return raw


def validate_task_record(record: dict, experiment: dict, prompts: VerifiedPromptExport,
                         root: Path, evidence: dict[str, bytes]) -> tuple[list[dict], dict | None]:
    """Validate one execution event, preserving format-invalid completed answers."""
    issues = []
    try:
        _keys(record, RECORD_FIELDS, "capture record")
        task_id = record["task_id"]
        if not isinstance(task_id, str) or task_id not in prompts.benchmark_ids:
            return [_issue("UNKNOWN_TASK_ID", "task_id", "Unknown benchmark task ID.")], None
        if task_id not in prompts.partition_ids:
            return [_issue("WRONG_PARTITION", "task_id", "Task belongs to the other partition.")], None
        if task_id not in prompts.records:
            return [_issue("OUTSIDE_EXPORTED_SELECTION", "task_id", "Task is outside this exported batch.")], None
        attempt = record["attempt"]
        if type(attempt) is not int or not 1 <= attempt <= experiment["retry_policy"]["max_attempts"]:
            raise ValueError("Attempt must be a positive integer within the declared bound.")
        timestamp = _timestamp(record["timestamp"])
        if timestamp < _timestamp(experiment["execution_started_at"]):
            raise ValueError("Task timestamp precedes experiment execution.")
        if attempt == 1:
            if record["retry_of"] is not None or record["retry_reason"] is not None:
                raise ValueError("Initial attempt requires null retry metadata.")
        elif type(record["retry_of"]) is not int or record["retry_of"] != attempt - 1 or not _text(record["retry_reason"]):
            raise ValueError("Retry requires the immediately preceding attempt and an explicit reason.")
        if record["outcome"] not in OUTCOMES or type(record["output_error"]) is not bool:
            raise ValueError("Invalid execution outcome or output-error flag.")
        if record["outcome"] != "completed" and not _text(record["outcome_details"]):
            raise ValueError("Non-completed outcomes require an explicit description or decision.")
        if record["outcome_details"] is not None and not _text(record["outcome_details"]):
            raise ValueError("Outcome details must be nonempty text or null.")
        tool = record["tool_use"]
        _keys(tool, {"status", "details"}, "tool use")
        if tool["status"] not in ("none", "recorded", "unavailable"):
            raise ValueError("Invalid tool-use status.")
        if tool["status"] != "none" and not _text(tool["details"]):
            raise ValueError("Recorded or unavailable tool use needs details.")
        if tool["status"] == "none" and tool["details"] is not None:
            raise ValueError("No tool use requires null details.")
        if tool["status"] == "recorded" and experiment["tool_policy"]["mode"] == "disabled":
            raise ValueError("Observed tool use conflicts with the declared policy.")
        verify_prompt_identity(record, prompts, root, evidence)
        raw = verify_response_bytes(record["response"], root, evidence)
        has_answer = record["outcome"] in {"completed", "refused", "truncated"}
        if has_answer != (raw is not None):
            raise ValueError("Execution outcome and answer-file presence are inconsistent.")
        if (has_answer and (experiment["access_method"] == "api"
                           or experiment["provenance"]["basis"] == "provider_payload")
                and record["response"]["source_payload"] is None):
            raise ValueError("API output requires preserved native payload and extraction details.")
        format_errors = []
        format_status = "ABSENT"
        if record["outcome"] == "abstained":
            format_status = "ABSTAINED"
        elif raw is not None:
            try:
                answer = parse_model_response(raw.decode("utf-8"))
                format_errors = validate_model_response(answer)
            except ResponseFormatError as exc:
                format_errors = exc.issues
            format_status = "INVALID" if format_errors else "VALID"
        return issues, {"task_id": task_id, "attempt": attempt, "outcome": record["outcome"],
                        "timestamp": record["timestamp"], "raw_response": raw.decode("utf-8") if raw is not None else None,
                        "model_output_format": format_status, "model_format_errors": format_errors}
    except (ValueError, TypeError, KeyError, OSError, UnicodeError) as exc:
        message = "Capture evidence is missing or inaccessible." if isinstance(exc, OSError) else str(exc)
        issues.append(_issue("INVALID_CAPTURE_RECORD", "record", message))
    return issues, None


@dataclass
class CaptureValidation:
    """Validation report plus preserved evidence used by deterministic conversion."""
    report: dict
    experiment: dict = field(default_factory=dict)
    primary: dict[str, dict] = field(default_factory=dict)
    evidence: dict[str, bytes] = field(default_factory=dict)
    snapshots: dict[str, bytes] = field(default_factory=dict)


def validate_capture(root: Path, prompts: VerifiedPromptExport, *, imported: bool = False) -> CaptureValidation:
    """Validate metadata, all attempts and checksums, reporting missing tasks explicitly."""
    result = CaptureValidation({
        "capture_schema_version": CAPTURE_SCHEMA_VERSION, "capture_version": CAPTURE_VERSION,
        "valid": False, "selected_task_ids": sorted(prompts.records),
        "selected_task_count": len(prompts.records), "issues": [], "attempts": [],
        "selected_attempts": {}, "missing_task_ids": [], "without_response_task_ids": [],
        "provenance_confidence": "researcher_supplied_origin_unverified",
    })
    try:
        names = {"experiment.json": "experiment.snapshot.json", "records.json": "records.snapshot.json"} if imported else {
            "experiment.json": "experiment.json", "records.json": "records.json",
        }
        result.snapshots = {name: safe_capture_path(root, filename).read_bytes()
                            for name, filename in names.items()}
        evidence_root = Path(root) / "evidence" if imported else Path(root)
        if imported and (evidence_root.exists() or evidence_root.is_symlink()):
            if not evidence_root.resolve(strict=True).is_relative_to(Path(root).resolve(strict=True)):
                raise ValueError("Imported evidence directory escapes the package root.")
        result.experiment = _json_object(result.snapshots["experiment.json"])
        result.report["issues"].extend(validate_experiment_metadata(result.experiment, prompts))
        if result.report["issues"]:
            return result
        result.report["submission_kind"] = result.experiment["submission_kind"]
        if result.experiment["submission_kind"] == "synthetic_scorer_testing":
            result.report["provenance_confidence"] = "synthetic_capture_system_testing"
        unavailable = []
        if result.experiment["model"]["version"]["availability"] == "unavailable":
            unavailable.append("model.version")
        unavailable.extend("instructions." + key for key, value in result.experiment["instructions"].items()
                           if value["availability"] == "unavailable")
        unavailable.extend("generation_parameters." + key
                           for key, value in result.experiment["generation_parameters"].items()
                           if value["availability"] == "unavailable")
        if result.experiment["tool_policy"]["mode"] == "unavailable":
            unavailable.append("tool_policy")
        for name in sorted(unavailable):
            result.report["issues"].append(_issue("UNAVAILABLE_PROVENANCE_FIELD", name,
                                                  "Availability limitation explicitly documented.", "WARNING"))
        records = _json_object(result.snapshots["records.json"])
        _keys(records, {"capture_schema_version", "records"}, "records envelope")
        if records["capture_schema_version"] != CAPTURE_SCHEMA_VERSION or not isinstance(records["records"], list):
            raise ValueError("Invalid capture records envelope.")
        groups = defaultdict(list)
        seen = set()
        present = set()
        for index, record in enumerate(records["records"]):
            issues, event = validate_task_record(record, result.experiment, prompts, evidence_root, result.evidence)
            for issue in issues:
                issue["record_index"] = index
            result.report["issues"].extend(issues)
            if event is None:
                continue
            key = (event["task_id"], event["attempt"])
            if key in seen:
                result.report["issues"].append(_issue("DUPLICATE_TASK_ATTEMPT", "records", "Repeated task/attempt identity."))
            seen.add(key)
            present.add(event["task_id"])
            groups[event["task_id"]].append(event)
        for task_id, events in sorted(groups.items()):
            ordered = sorted(events, key=lambda e: e["attempt"])
            if [e["attempt"] for e in ordered] != list(range(1, len(ordered) + 1)):
                result.report["issues"].append(_issue("INCOMPLETE_ATTEMPT_HISTORY", task_id, "Attempt numbers must be contiguous from one."))
            times = [_timestamp(e["timestamp"]) for e in ordered]
            if times != sorted(times):
                result.report["issues"].append(_issue("INCONSISTENT_ATTEMPT_TIME", task_id, "Attempt timestamps must be chronological."))
            eligible = [e for e in ordered if e["outcome"] in ELIGIBLE_OUTCOMES]
            if eligible:
                result.primary[task_id] = eligible[0]
            else:
                result.report["without_response_task_ids"].append(task_id)
            result.report["attempts"].extend({k: v for k, v in e.items() if k != "raw_response"} for e in ordered)
        result.report["missing_task_ids"] = sorted(set(prompts.records) - present)
        for task_id in result.report["missing_task_ids"]:
            result.report["issues"].append(_issue("MISSING_TASK_RECORD", task_id, "No execution record supplied.", "WARNING"))
        for task_id in result.report["without_response_task_ids"]:
            result.report["issues"].append(_issue("NO_RESPONSE_OBTAINED", task_id, "Execution recorded without an eligible response.", "WARNING"))
        result.report["selected_attempts"] = {task_id: e["attempt"] for task_id, e in sorted(result.primary.items())}
        result.report["captured_attempt_count"] = len(result.report["attempts"])
        result.report["primary_response_count"] = len(result.primary)
        result.report["missing_task_count"] = len(result.report["missing_task_ids"])
        result.report["without_response_task_count"] = len(result.report["without_response_task_ids"])
        result.report["format_valid_answer_count"] = sum(e["model_output_format"] == "VALID" for e in result.report["attempts"])
        result.report["format_invalid_answer_count"] = sum(e["model_output_format"] == "INVALID" for e in result.report["attempts"])
        result.report["valid"] = not any(i["severity"] == "ERROR" for i in result.report["issues"])
    except (ValueError, TypeError, KeyError, OSError, UnicodeError) as exc:
        message = "Capture metadata is missing or inaccessible." if isinstance(exc, OSError) else str(exc)
        result.report["issues"].append(_issue("INVALID_CAPTURE_PACKAGE", "capture", message))
    return result


def convert_to_scoring_submission(validation: CaptureValidation) -> dict:
    """Convert only integrity-valid primary observations, retaining raw strings."""
    if not validation.report["valid"]:
        raise ValueError("Capture provenance validation failed; submission conversion is blocked.")
    rows = []
    for task_id, event in sorted(validation.primary.items()):
        rows.append({"task_id": task_id, "raw_response": event["raw_response"],
                     "abstain": event["outcome"] == "abstained"})
    return {"schema_version": "0.1", "submission_kind": validation.experiment["submission_kind"],
            "responses": rows}


def load_verified_capture_import(root: Path, prompts: VerifiedPromptExport) -> tuple[CaptureValidation, dict[str, bytes]]:
    """Replay imported source evidence and reject altered derived answers or reports.

    Hash consistency alone cannot establish which attempt supplied a submission.
    This function reconstructs the declared capture from its preserved snapshots.
    """
    root = Path(root)
    manifest_raw = safe_capture_path(root, "manifest.json").read_bytes()
    manifest = _json_object(manifest_raw)
    checksums = manifest.get("output_checksums")
    evidence = manifest.get("evidence_checksums")
    if not isinstance(checksums, dict) or not isinstance(evidence, dict):
        raise ValueError("Invalid capture-import checksum maps.")
    expected_names = IMPORT_CORE_FILES | {"evidence/" + name for name in evidence}
    if set(checksums) != expected_names:
        raise ValueError("Capture-import checksum coverage mismatch.")
    snapshots = {name: safe_capture_path(root, name).read_bytes() for name in sorted(checksums)}
    for name, raw in snapshots.items():
        if sha256(raw) != checksums[name]:
            raise ValueError(f"Capture-import checksum mismatch for {name}")
    if (snapshots["prompt_export_manifest.snapshot.json"] != prompts.snapshots["manifest.json"]
            or snapshots["prompts.snapshot.json"] != prompts.snapshots["prompts.json"]):
        raise ValueError("Imported prompt snapshots differ from the verified export.")
    expected_metadata = {
        "capture_schema_version": CAPTURE_SCHEMA_VERSION, "capture_version": CAPTURE_VERSION,
        "prompt_export_identity": prompts.identity, "selected_task_ids": sorted(prompts.records),
        "selected_task_count": len(prompts.records), "partition": prompts.manifest["partition"],
        "stage": prompts.manifest["stage"],
        "benchmark_file_checksums": {name: sha256(raw) for name, raw in prompts.benchmark_snapshots.items()},
    }
    if any(manifest.get(key) != value for key, value in expected_metadata.items()):
        raise ValueError("Capture-import population or provenance mismatch.")
    validation = validate_capture(root, prompts, imported=True)
    if not validation.report["valid"]:
        raise ValueError("Imported source capture replay failed provenance validation.")
    if snapshots["submission.json"] != canonical_json(convert_to_scoring_submission(validation)):
        raise ValueError("Imported submission differs from replayed primary response bytes.")
    if snapshots["capture_validation.json"] != canonical_json(validation.report):
        raise ValueError("Imported capture validation report differs from source-attempt replay.")
    if (manifest.get("experiment_id") != validation.experiment["experiment_id"]
            or manifest.get("submission_kind") != validation.experiment["submission_kind"]
            or manifest.get("provenance_confidence") != validation.report["provenance_confidence"]
            or manifest.get("capture_metadata_checksums") != {name: sha256(raw) for name, raw in validation.snapshots.items()}
            or evidence != {name: sha256(raw) for name, raw in validation.evidence.items()}):
        raise ValueError("Imported evidence metadata differs from replayed source provenance.")
    snapshots["manifest.json"] = manifest_raw
    return validation, snapshots
