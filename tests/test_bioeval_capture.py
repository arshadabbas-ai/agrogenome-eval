"""Synthetic capture-system testing only, using fictional execution metadata."""

from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import shutil
import socket
import sys

import pytest

from agrogenome_eval import bioeval_capture as capture
from agrogenome_eval.bioeval_scoring import evaluate_responses


ROOT = Path(__file__).parents[1]
BENCHMARK = ROOT / "data/benchmarks/bioeval_v01"


def script(name):
    directory = str(ROOT / "scripts")
    if directory not in sys.path:
        sys.path.insert(0, directory)
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def exports(tmp_path_factory):
    # CI can run this feature without any ignored local review artifacts.
    root = tmp_path_factory.mktemp("capture-prompts")
    exporter = script("export_bioeval_model_prompts")
    exporter.export_prompts(BENCHMARK, "development", root / "pilot", pilot=True)
    exporter.export_prompts(BENCHMARK, "evaluation", root / "evaluation")
    return root


@pytest.fixture(scope="module")
def prompts(exports):
    return capture.load_verified_prompt_export(exports / "pilot", BENCHMARK)


def metadata(prompts):
    unknown = {"availability": "unavailable", "value": None, "reason": "Synthetic fixture has no model execution."}
    return {
        "capture_schema_version": "0.1", "experiment_id": "synthetic-capture-system-testing",
        "submission_kind": "synthetic_scorer_testing", "benchmark_schema_version": "0.1",
        "prompt_export": prompts.identity,
        "model": {"source": "synthetic test harness", "name": "fictional capture fixture",
                  "version": {"availability": "recorded", "value": "fixture-v1", "reason": None}},
        "access_method": "manual_interface", "execution_started_at": "2000-01-01T00:00:00+00:00",
        "instructions": {name: {"availability": "absent", "text": None, "reason": None}
                         for name in ["system", "developer"]},
        "generation_parameters": {"temperature": deepcopy(unknown), "seed": deepcopy(unknown)},
        "tool_policy": {"mode": "disabled", "details": "Synthetic fixture uses no tools."},
        "retry_policy": {"selection": "first_completed_eligible", "max_attempts": 3},
        "provenance": {"basis": "synthetic_fixture", "attestation": "Synthetic capture-system testing only.",
                       "capture_limitations": "Fictional fixture dates and identity, with no real model observations."},
    }


def synthetic_answer(task_id):
    # Independently read from the first pilot question: wheat position 629,
    # CDS start 96 => CDS 534, codon 178, TAC>TAG gives Y>* under table 1.
    return {"task_id": task_id, "transcript_accession": "XM_044542390.1",
            "transcript_region": "CDS", "cds_position": 534, "codon_number": 178,
            "reference_codon": "TAC", "alternative_codon": "TAG", "reference_amino_acid": "Y",
            "alternative_amino_acid": "*", "molecular_consequence": "STOP_GAINED",
            "analysis_status": "ASSESSED", "unsupported_region": None,
            "phenotype_claim_support": "NOT_ESTABLISHED"}


def add_record(root, prompts, task_id=None, *, raw=b"{malformed\r\n ", outcome="completed", attempt=1):
    task_id = task_id or sorted(prompts.records)[0]
    prompt_bytes = prompts.records[task_id]["question"].encode("utf-8")
    prefix = "files/" + task_id + "-" + str(attempt)
    prompt_path = prefix + ".prompt.txt"
    path = root / prompt_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(prompt_bytes)
    response_path = prefix + ".answer.txt" if raw is not None else None
    if raw is not None:
        (root / response_path).write_bytes(raw)
    return {
        "task_id": task_id, "attempt": attempt, "retry_of": attempt - 1 if attempt > 1 else None,
        "retry_reason": "Synthetic additional-attempt scenario." if attempt > 1 else None,
        "timestamp": f"2000-01-01T00:00:0{attempt}+00:00", "outcome": outcome,
        "outcome_details": "Synthetic documented " + outcome if outcome != "completed" else None,
        "prompt": {"path": prompt_path, "sha256": capture.sha256(prompt_bytes),
                   "record_sha256": capture.sha256(capture.canonical_json(prompts.records[task_id]))},
        "response": {"path": response_path, "sha256": capture.sha256(raw) if raw is not None else None,
                     "source_payload": None},
        "output_error": outcome in {"execution_failure", "truncated"},
        "tool_use": {"status": "none", "details": None},
    }


def save(root, experiment, rows):
    root.mkdir(parents=True, exist_ok=True)
    (root / "experiment.json").write_bytes(capture.canonical_json(experiment))
    (root / "records.json").write_bytes(capture.canonical_json({"capture_schema_version": "0.1", "records": rows}))


@pytest.fixture
def source(tmp_path, prompts):
    root = tmp_path / "capture"
    root.mkdir()
    task_id = sorted(prompts.records)[0]
    row = add_record(root, prompts, raw=capture.canonical_json(synthetic_answer(task_id)))
    experiment = metadata(prompts)
    save(root, experiment, [row])
    return root, experiment, [row]


def assert_invalid(validation, code=None):
    assert not validation.report["valid"]
    if code:
        assert code in {i["code"] for i in validation.report["issues"]}
    with pytest.raises(ValueError, match="blocked"):
        capture.convert_to_scoring_submission(validation)


def test_verified_export_and_independent_hashes(prompts):
    assert capture.sha256(b"abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    assert len(prompts.records) == 14
    assert set(prompts.records) <= prompts.partition_ids
    assert prompts.identity["stage"] == "pilot"
    assert prompts.identity["selection_id"] == capture.sha256(capture.canonical_json(sorted(prompts.records)))


def test_metadata_and_zero_vs_unavailable(prompts):
    experiment = metadata(prompts)
    assert capture.validate_experiment_metadata(experiment, prompts) == []
    experiment["generation_parameters"]["temperature"] = {"availability": "set", "value": 0, "reason": None}
    assert capture.validate_experiment_metadata(experiment, prompts) == []
    experiment["generation_parameters"]["temperature"]["availability"] = "unavailable"
    assert capture.validate_experiment_metadata(experiment, prompts)


@pytest.mark.parametrize("change", ["schema", "identity", "timestamp", "offset", "seed_bool", "unknown_reason", "instruction", "attestation", "secret"])
def test_invalid_metadata_blocks_conversion(source, prompts, change):
    root, experiment, rows = source
    if change == "schema":
        experiment["capture_schema_version"] = "9"
    elif change == "identity":
        experiment["model"]["name"] = " "
    elif change == "timestamp":
        del experiment["execution_started_at"]
    elif change == "offset":
        experiment["execution_started_at"] = "2000-01-01T00:00:00"
    elif change == "seed_bool":
        experiment["generation_parameters"]["seed"] = {"availability": "set", "value": True, "reason": None}
    elif change == "unknown_reason":
        experiment["generation_parameters"]["seed"]["reason"] = None
    elif change == "instruction":
        experiment["instructions"]["system"]["text"] = "Unexpected hidden instruction"
    elif change == "secret":
        experiment["generation_parameters"]["api_key"] = {"availability": "set", "value": "synthetic-forbidden-field", "reason": None}
    else:
        experiment["provenance"]["attestation"] = ""
    save(root, experiment, rows)
    assert_invalid(capture.validate_capture(root, prompts), "INVALID_EXPERIMENT_METADATA")


def test_synthetic_cannot_silently_be_promoted_to_model_responses(prompts):
    experiment = metadata(prompts)
    experiment["submission_kind"] = "model_responses"
    assert capture.validate_experiment_metadata(experiment, prompts)
    experiment["provenance"]["basis"] = "researcher_attestation"
    experiment["model"]["source"] = "unknown"
    assert capture.validate_experiment_metadata(experiment, prompts)


@pytest.mark.parametrize("raw", [b"{", b"```json\n{}\n```\n", b"  Unmodified REFUSAL.\r\n\r\n", b"", b"__import__('os').system('never execute')", b"\xef\xbb\xbf\xc3\xa9\r\n  "])
def test_invalid_answer_bytes_are_preserved(source, prompts, raw):
    root, experiment, _ = source
    row = add_record(root, prompts, raw=raw)
    save(root, experiment, [row])
    result = capture.validate_capture(root, prompts)
    assert result.report["valid"]
    assert result.report["format_invalid_answer_count"] == 1
    submission = capture.convert_to_scoring_submission(result)
    assert submission["submission_kind"] == "synthetic_scorer_testing"
    assert submission["responses"][0]["raw_response"].encode("utf-8") == raw
    assert result.evidence[row["response"]["path"]] == raw


def test_wrong_biology_remains_a_format_valid_observation(source, prompts):
    root, experiment, _ = source
    task_id = sorted(prompts.records)[0]
    answer = synthetic_answer(task_id)
    answer["alternative_amino_acid"] = "A"
    row = add_record(root, prompts, raw=capture.canonical_json(answer))
    save(root, experiment, [row])
    result = capture.validate_capture(root, prompts)
    assert result.report["valid"] and result.report["format_valid_answer_count"] == 1
    assert json.loads(capture.convert_to_scoring_submission(result)["responses"][0]["raw_response"])["alternative_amino_acid"] == "A"


@pytest.mark.parametrize("outcome,raw,rows_expected,abstain", [
    ("refused", b"I decline this synthetic task.\n", 1, False),
    ("truncated", b'{"task_id":', 1, False),
    ("abstained", None, 1, True),
    ("no_response", None, 0, False),
    ("execution_failure", None, 0, False),
])
def test_explicit_outcome_mapping(source, prompts, outcome, raw, rows_expected, abstain):
    root, experiment, _ = source
    row = add_record(root, prompts, raw=raw, outcome=outcome)
    save(root, experiment, [row])
    result = capture.validate_capture(root, prompts)
    assert result.report["valid"]
    rows = capture.convert_to_scoring_submission(result)["responses"]
    assert len(rows) == rows_expected
    if rows:
        assert rows[0]["abstain"] is abstain
        assert rows[0]["raw_response"] == (raw.decode() if raw is not None else None)
    else:
        assert result.report["without_response_task_count"] == 1


def test_retry_after_failure_and_first_invalid_completion_policy(source, prompts):
    root, experiment, _ = source
    first = add_record(root, prompts, outcome="execution_failure", raw=None)
    second = add_record(root, prompts, attempt=2, raw=b"{broken first completion")
    third = add_record(root, prompts, attempt=3, raw=b"later corrected content")
    save(root, experiment, [third, first, second])
    result = capture.validate_capture(root, prompts)
    assert result.report["valid"]
    assert result.report["selected_attempts"] == {first["task_id"]: 2}
    assert capture.convert_to_scoring_submission(result)["responses"][0]["raw_response"] == "{broken first completion"
    assert result.report["captured_attempt_count"] == 3
    assert third["response"]["path"] in result.evidence


@pytest.mark.parametrize("change", ["duplicate", "gap", "retry_link", "retry_reason", "time", "missing_timestamp", "outcome", "abstention_text", "prompt_hash", "record_hash", "response_hash", "edited_bytes", "utf8"])
def test_inconsistent_records_block_conversion(source, prompts, change):
    root, experiment, rows = source
    row = rows[0]
    if change == "duplicate": rows.append(deepcopy(row))
    elif change == "gap": row.update(attempt=2, retry_of=1, retry_reason="Synthetic missing history")
    elif change == "retry_link": row["retry_of"] = 0
    elif change == "retry_reason": row["retry_reason"] = "Unlinked retry"
    elif change == "time": row["timestamp"] = "1999-12-31T23:59:59Z"
    elif change == "missing_timestamp": del row["timestamp"]
    elif change == "outcome": row["outcome"] = "guessed"
    elif change == "abstention_text": row.update(outcome="abstained", outcome_details="Synthetic abstention")
    elif change == "prompt_hash": row["prompt"]["sha256"] = "0" * 64
    elif change == "record_hash": row["prompt"]["record_sha256"] = "0" * 64
    elif change == "response_hash": row["response"]["sha256"] = "0" * 64
    elif change == "edited_bytes": (root / row["response"]["path"]).write_bytes(b"changed")
    elif change == "utf8":
        (root / row["response"]["path"]).write_bytes(b"\xff")
        row["response"]["sha256"] = capture.sha256(b"\xff")
    save(root, experiment, rows)
    assert_invalid(capture.validate_capture(root, prompts))


@pytest.mark.parametrize("kind", ["unknown", "other_partition", "outside_pilot"])
def test_task_population_mismatch(source, prompts, kind):
    root, experiment, rows = source
    if kind == "unknown": rows[0]["task_id"] = "UNKNOWN-SYNTHETIC-ID"
    elif kind == "other_partition": rows[0]["task_id"] = sorted(prompts.benchmark_ids - prompts.partition_ids)[0]
    else: rows[0]["task_id"] = sorted(prompts.partition_ids - set(prompts.records))[0]
    save(root, experiment, rows)
    code = {"unknown": "UNKNOWN_TASK_ID", "other_partition": "WRONG_PARTITION", "outside_pilot": "OUTSIDE_EXPORTED_SELECTION"}[kind]
    assert_invalid(capture.validate_capture(root, prompts), code)


@pytest.mark.parametrize("path", ["../outside.txt", "/absolute.txt", "C:/private.txt", "C:relative.txt", "//host/file", "files\\answer.txt", "files/answer.txt:stream"])
def test_unsafe_paths_are_rejected(source, prompts, path):
    root, experiment, rows = source
    rows[0]["response"]["path"] = path
    save(root, experiment, rows)
    assert_invalid(capture.validate_capture(root, prompts), "INVALID_CAPTURE_RECORD")


def test_symlink_escape_where_supported_with_portable_resolution_fallback(tmp_path, monkeypatch):
    root = tmp_path / "capture"
    root.mkdir()
    outside = tmp_path / "outside.txt"
    outside.write_bytes(b"synthetic outside evidence")
    link = root / "escape.txt"
    try:
        link.symlink_to(outside)
    except OSError:
        # Windows can prohibit symlink creation. Exercise the same resolved-path
        # containment check without requesting extra filesystem privileges.
        original = Path.resolve
        def resolve(path, *args, **kwargs):
            if path == link:
                return outside
            return original(path, *args, **kwargs)
        monkeypatch.setattr(Path, "resolve", resolve)
    with pytest.raises(ValueError, match="escapes"):
        capture.safe_capture_path(root, "escape.txt")


def test_api_payload_and_extracted_text_remain_separate(source, prompts):
    root, experiment, rows = source
    experiment["access_method"] = "api"
    save(root, experiment, rows)
    assert_invalid(capture.validate_capture(root, prompts))
    payload = b'{"synthetic_native_payload": "a separately preserved transport artifact"}\r\n'
    (root / "payload.json").write_bytes(payload)
    rows[0]["response"]["source_payload"] = {"path": "payload.json", "sha256": capture.sha256(payload),
                                               "extraction": "Synthetic independent extraction scenario; origin unverified."}
    save(root, experiment, rows)
    result = capture.validate_capture(root, prompts)
    assert result.report["valid"]
    assert result.evidence["payload.json"] == payload
    assert result.evidence[rows[0]["response"]["path"]] != payload


def test_pilot_denominator_and_current_scoring_schema(source, prompts):
    root, _, _ = source
    result = capture.validate_capture(root, prompts)
    submission = capture.convert_to_scoring_submission(result)
    assert set(submission) == {"schema_version", "submission_kind", "responses"}
    assert set(submission["responses"][0]) == {"task_id", "raw_response", "abstain"}
    tasks = json.loads((BENCHMARK / "tasks.json").read_bytes())
    gold = json.loads((BENCHMARK / "gold_answers.json").read_bytes())
    report = evaluate_responses(tasks, gold, submission, task_ids=result.report["selected_task_ids"])
    assert report["summary"]["total_selected_tasks"] == 14
    assert report["summary"]["missing_response_count"] == 13
    assert report["submission_kind"] == "synthetic_scorer_testing"


@pytest.mark.parametrize("name", ["manifest.json", "prompts.json"])
def test_corrupt_export_is_rejected(exports, tmp_path, name):
    target = tmp_path / "export"
    shutil.copytree(exports / "pilot", target)
    path = target / name
    if name == "manifest.json":
        value = json.loads(path.read_bytes())
        value["pilot_sampling_seed"] += 1
        path.write_bytes(capture.canonical_json(value))
    else:
        path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError):
        capture.load_verified_prompt_export(target, BENCHMARK)


def test_import_byte_replay_and_source_preservation(source, exports, tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Networking forbidden during synthetic capture testing")
    monkeypatch.setattr(socket.socket, "connect", forbidden)
    root, _, rows = source
    importer = script("import_bioeval_model_responses")
    source_bytes = {p.relative_to(root).as_posix(): p.read_bytes() for p in root.rglob("*") if p.is_file()}
    first, second = tmp_path / "first", tmp_path / "second"
    report = importer.run_import(exports / "pilot", root, first)
    importer.run_import(exports / "pilot", root, second)
    importer.verify_import_outputs(first)
    for path in first.rglob("*"):
        if path.is_file():
            assert path.read_bytes() == (second / path.relative_to(first)).read_bytes()
    assert source_bytes == {p.relative_to(root).as_posix(): p.read_bytes() for p in root.rglob("*") if p.is_file()}
    assert (first / "evidence" / rows[0]["response"]["path"]).read_bytes() == source_bytes[rows[0]["response"]["path"]]
    assert report["selected_task_count"] == 14
    assert json.loads((first / "manifest.json").read_bytes())["selected_task_count"] == 14


@pytest.mark.parametrize("kind", ["file", "directory", "inside_source"])
def test_existing_output_and_source_separation(source, exports, tmp_path, kind):
    root, _, _ = source
    output = tmp_path / "output"
    if kind == "file": output.write_bytes(b"preserve")
    elif kind == "directory": output.mkdir()
    else: output = root / "derived"
    with pytest.raises((ValueError, FileExistsError)):
        script("import_bioeval_model_responses").run_import(exports / "pilot", root, output)
    if kind == "file": assert output.read_bytes() == b"preserve"


@pytest.mark.parametrize("failure", ["write", "rename"])
def test_failed_output_cleanup(source, exports, tmp_path, monkeypatch, failure):
    root, _, _ = source
    if failure == "write":
        original = Path.write_bytes
        def fail(path, raw):
            if path.name == "submission.json": raise OSError("synthetic publication failure")
            return original(path, raw)
        monkeypatch.setattr(Path, "write_bytes", fail)
    else:
        def fail(*args): raise OSError("synthetic publication failure")
        monkeypatch.setattr(Path, "rename", fail)
    with pytest.raises(OSError, match="synthetic publication"):
        script("import_bioeval_model_responses").run_import(exports / "pilot", root, tmp_path / "output")
    assert not (tmp_path / "output").exists()
    assert not list(tmp_path.glob(".capture-import-*"))


def test_bad_capture_never_publishes_submission(source, exports, tmp_path):
    root, experiment, rows = source
    rows[0]["response"]["sha256"] = "0" * 64
    save(root, experiment, rows)
    importer = script("import_bioeval_model_responses")
    with pytest.raises(importer.CaptureImportError) as exc:
        importer.run_import(exports / "pilot", root, tmp_path / "output")
    assert not exc.value.report["valid"] and not (tmp_path / "output").exists()


def test_tampered_derived_output_is_rejected(source, exports, tmp_path):
    root, _, _ = source
    importer = script("import_bioeval_model_responses")
    output = tmp_path / "output"
    importer.run_import(exports / "pilot", root, output)
    path = output / "submission.json"
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="checksum"):
        importer.verify_import_outputs(output)


def test_cli_success_and_invalid_report(source, exports, tmp_path, monkeypatch, capsys):
    root, experiment, rows = source
    importer = script("import_bioeval_model_responses")
    output = tmp_path / "output"
    argv = ["import_bioeval_model_responses.py", "--prompt-export", str(exports / "pilot"),
            "--capture", str(root), "--output", str(output)]
    monkeypatch.setattr(sys, "argv", argv)
    importer.main()
    assert json.loads(capsys.readouterr().out)["selected_task_count"] == 14
    rows[0]["prompt"]["sha256"] = "0" * 64
    save(root, experiment, rows)
    argv[-1] = str(tmp_path / "failed")
    with pytest.raises(SystemExit) as exc:
        importer.main()
    assert exc.value.code == 1
    assert json.loads(capsys.readouterr().err)["valid"] is False


def create_synthetic_demo(root, prompts):
    """Construct explicitly fictional capture-system scenarios for local review."""
    root.mkdir(parents=True, exist_ok=True)
    ids = sorted(prompts.records)
    first_answer = capture.canonical_json(synthetic_answer(ids[0]))
    wrong_answer = synthetic_answer(ids[1])
    wrong_answer["alternative_amino_acid"] = "A"
    rows = [
        add_record(root, prompts, ids[0], raw=first_answer),
        add_record(root, prompts, ids[1], raw=capture.canonical_json(wrong_answer)),
        add_record(root, prompts, ids[2], raw=b"{"),
        add_record(root, prompts, ids[3], raw=b"```json\n{}\n```"),
        add_record(root, prompts, ids[4], raw=None, outcome="abstained"),
        add_record(root, prompts, ids[5], raw=b"Synthetic refusal text.\r\n", outcome="refused"),
        add_record(root, prompts, ids[6], raw=b""),
        add_record(root, prompts, ids[7], raw=None, outcome="no_response"),
        add_record(root, prompts, ids[8], raw=None, outcome="execution_failure"),
        add_record(root, prompts, ids[8], raw=b"{", attempt=2),
        add_record(root, prompts, ids[9], raw=b'{"task_id":', outcome="truncated"),
        add_record(root, prompts, ids[10], raw=b"{"),
        add_record(root, prompts, ids[10], raw=first_answer, attempt=2),
    ]
    save(root, metadata(prompts), rows)
    return rows


def test_manually_specified_synthetic_demo(prompts, tmp_path):
    root = tmp_path / "synthetic"
    create_synthetic_demo(root, prompts)
    result = capture.validate_capture(root, prompts)
    report = result.report
    assert report["valid"] and report["provenance_confidence"] == "synthetic_capture_system_testing"
    assert report["selected_task_count"] == 14
    assert report["captured_attempt_count"] == 13
    assert report["primary_response_count"] == 10
    assert report["missing_task_count"] == 3
    assert report["without_response_task_count"] == 1
    assert report["format_valid_answer_count"] == 3
    assert report["format_invalid_answer_count"] == 7
    assert "UNAVAILABLE_PROVENANCE_FIELD" in {i["code"] for i in report["issues"]}
    submission = capture.convert_to_scoring_submission(result)
    assert submission["submission_kind"] == "synthetic_scorer_testing"
    assert len(submission["responses"]) == 10
    assert sum(row["abstain"] for row in submission["responses"]) == 1


def test_export_manifest_identity_mismatch_blocks_capture(source, prompts):
    root, experiment, rows = source
    experiment["prompt_export"]["manifest_sha256"] = "0" * 64
    save(root, experiment, rows)
    assert_invalid(capture.validate_capture(root, prompts), "INVALID_EXPERIMENT_METADATA")


def test_missing_file_diagnostics_are_independent_of_capture_location(source, prompts, tmp_path):
    root, experiment, rows = source
    (root / rows[0]["response"]["path"]).unlink()
    save(root, experiment, rows)
    second = tmp_path / "relocated"
    shutil.copytree(root, second)
    first_report = capture.validate_capture(root, prompts).report
    assert not first_report["valid"]
    assert first_report == capture.validate_capture(second, prompts).report


def test_duplicate_metadata_keys_are_rejected(source, prompts):
    root, _, _ = source
    (root / "experiment.json").write_bytes(b'{"capture_schema_version":"0.1","capture_schema_version":"0.1"}')
    assert_invalid(capture.validate_capture(root, prompts), "INVALID_CAPTURE_PACKAGE")


def test_source_change_during_publication_blocks_import(source, exports, tmp_path, monkeypatch):
    root, _, rows = source
    importer = script("import_bioeval_model_responses")
    original = importer.verify_import_outputs
    def mutate(*args):
        original(*args)
        (root / rows[0]["response"]["path"]).write_bytes(b"edited after validation")
    monkeypatch.setattr(importer, "verify_import_outputs", mutate)
    with pytest.raises(ValueError, match="changed during import"):
        importer.run_import(exports / "pilot", root, tmp_path / "output")
    assert not (tmp_path / "output").exists()
    assert not list(tmp_path.glob(".capture-import-*"))
