"""Independent pilot preflight tests using explicitly synthetic capture fixtures."""

from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import shutil
import socket
import sys

import pytest

from agrogenome_eval import bioeval_capture as capture


ROOT = Path(__file__).parents[1]
BENCHMARK = ROOT / "data/benchmarks/bioeval_v01"


def load_script(name):
    directory = str(ROOT / "scripts")
    if directory not in sys.path:
        sys.path.insert(0, directory)
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def helpers():
    spec = importlib.util.spec_from_file_location("synthetic_capture_helpers", ROOT / "tests/test_bioeval_capture.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def export(tmp_path_factory):
    output = tmp_path_factory.mktemp("pilot-preflight") / "export"
    load_script("export_bioeval_model_prompts").export_prompts(BENCHMARK, "development", output, pilot=True)
    return output


@pytest.fixture(scope="module")
def prompts(export):
    return capture.load_verified_prompt_export(export, BENCHMARK)


@pytest.fixture
def imported(tmp_path, helpers, export, prompts):
    source = tmp_path / "source"
    helpers.create_synthetic_demo(source, prompts)
    output = tmp_path / "import"
    load_script("import_bioeval_model_responses").run_import(export, source, output)
    return source, output


def rewrite_payload_and_hash(output, name, value):
    raw = capture.canonical_json(value)
    (output / name).write_bytes(raw)
    manifest = json.loads((output / "manifest.json").read_bytes())
    manifest["output_checksums"][name] = capture.sha256(raw)
    (output / "manifest.json").write_bytes(capture.canonical_json(manifest))


def test_audit_submission_must_match_preserved_primary_bytes(imported):
    _, output = imported
    submission = json.loads((output / "submission.json").read_bytes())
    submission["responses"][0]["raw_response"] = "edited derived answer"
    rewrite_payload_and_hash(output, "submission.json", submission)
    with pytest.raises(ValueError, match="replay|primary|submission"):
        load_script("import_bioeval_model_responses").verify_import_outputs(output)


def test_audit_validation_report_must_match_attempt_replay(imported):
    _, output = imported
    report = json.loads((output / "capture_validation.json").read_bytes())
    report["captured_attempt_count"] = 999
    rewrite_payload_and_hash(output, "capture_validation.json", report)
    with pytest.raises(ValueError, match="replay|validation"):
        load_script("import_bioeval_model_responses").verify_import_outputs(output)


@pytest.fixture
def perfect_source(tmp_path, helpers, prompts):
    source = tmp_path / "perfect-source"
    source.mkdir()
    gold = json.loads((BENCHMARK / "gold_answers.json").read_bytes())["answers"]
    # Reference copies are explicitly synthetic contract fixtures. Operational
    # capture/import/scoring code never substitutes these values for observations.
    rows = [helpers.add_record(source, prompts, task_id, raw=capture.canonical_json(gold[task_id]))
            for task_id in sorted(prompts.records)]
    experiment = helpers.metadata(prompts)
    helpers.save(source, experiment, rows)
    return source, experiment, rows


def import_source(source, export, tmp_path):
    output = tmp_path / "validated-import"
    load_script("import_bioeval_model_responses").run_import(export, source, output)
    return output


def run_scoring(export, imported, output):
    return load_script("score_bioeval_pilot").run_pilot(export, imported, output)


def test_perfect_synthetic_population_and_composition(perfect_source, export, prompts, tmp_path):
    source, _, _ = perfect_source
    imported = import_source(source, export, tmp_path)
    report = run_scoring(export, imported, tmp_path / "scored")
    summary = report["summary"]
    assert summary["total_selected_tasks"] == 14
    assert summary["valid_response_count"] == summary["exact_response_count"] == 14
    assert summary["score_sum"] == 1400
    assert summary["strict_exact_match_accuracy"] == 1
    assert summary["coverage_adjusted_mean_score"] == 100
    assert report["population"]["selected_task_ids"] == sorted(prompts.records)
    assert report["population"]["counts_by_species"] == {"Zea mays": 7, "Triticum aestivum": 7}
    assert len(report["population"]["counts_by_category"]) == 7
    assert all(n == 2 for n in report["population"]["counts_by_category"].values())
    assert {name: group["total_selected_tasks"] for name, group in summary["per_analysis_status"].items()} == {
        "ASSESSED": 6, "NOT_ASSESSED": 4, "UNSUPPORTED": 4,
    }
    assert report["submission_kind"] == "synthetic_scorer_testing"


@pytest.mark.parametrize("changes,score", [
    ({"cds_position": 533}, 90),
    ({"cds_position": 537, "codon_number": 179}, 80),
    ({"analysis_status": "NOT_ASSESSED"}, 25),
    ({"phenotype_claim_support": "ESTABLISHED"}, 50),
    ({"molecular_consequence": "MISSENSE"}, 50),
    ({"alternative_amino_acid": "A"}, 50),
])
def test_independent_partial_scores(perfect_source, export, tmp_path, changes, score):
    source, _, rows = perfect_source
    path = source / rows[0]["response"]["path"]
    answer = json.loads(path.read_bytes())
    answer.update(changes)
    raw = capture.canonical_json(answer)
    path.write_bytes(raw)
    rows[0]["response"]["sha256"] = capture.sha256(raw)
    experiment = json.loads((source / "experiment.json").read_bytes())
    (source / "records.json").write_bytes(capture.canonical_json({"capture_schema_version": "0.1", "records": rows}))
    imported = import_source(source, export, tmp_path)
    report = run_scoring(export, imported, tmp_path / "scored")
    assert report["tasks"][0]["score"] == score
    assert report["summary"]["score_sum"] == 1300 + score
    assert report["summary"]["strict_exact_match_accuracy"] == 13 / 14
    assert report["summary"]["coverage_adjusted_mean_score"] == (1300 + score) / 14
    assert experiment["submission_kind"] == report["submission_kind"]


def test_synthetic_demo_manual_metrics(imported, export, tmp_path):
    _, package = imported
    report = run_scoring(export, package, tmp_path / "scored")
    s = report["summary"]
    assert report["pilot_outcome_counts"] == {
        "VALID_SCORED_RESPONSE": 2, "INVALID_SCORED_RESPONSE": 7,
        "EXPLICIT_ABSTENTION": 1, "MISSING_OR_NO_ELIGIBLE_RESPONSE": 4,
    }
    assert s["score_sum"] == 100 + (10 + 15 + 5) == 130
    assert s["exact_response_count"] == 1
    assert s["strict_exact_match_accuracy"] == 1 / 14
    assert s["exact_match_accuracy_among_valid"] == 1 / 2
    assert s["mean_partial_credit_score"] == 65
    assert s["coverage_adjusted_mean_score"] == 130 / 14
    assert s["response_coverage"] == 10 / 14
    assert s["response_format_compliance_rate"] == 3 / 10
    execution = report["execution_accounting"]
    assert execution["total_execution_attempts"] == 13
    assert execution["captured_task_count"] == 11
    assert execution["eligible_primary_responses"] == 10
    assert execution["missing_task_record_count"] == 3
    assert execution["recorded_tasks_without_eligible_response"] == 1
    assert execution["retry_attempt_count"] == 2 and execution["secondary_eligible_attempt_count"] == 1
    assert execution["outcome_counts_all_attempts"]["execution_failure"] == 1
    assert "Synthetic end-to-end pilot preflight" in report["evaluation_description"]
    assert len(report["tasks"]) == 14


@pytest.mark.parametrize("mode", ["missing", "malformed", "abstained", "failed"])
def test_empty_valid_populations_preserve_denominators(perfect_source, export, helpers, prompts, tmp_path, mode):
    source, experiment, _ = perfect_source
    rows = []
    if mode != "missing":
        for task_id in sorted(prompts.records):
            outcome = {"malformed": "completed", "abstained": "abstained", "failed": "execution_failure"}[mode]
            rows.append(helpers.add_record(source, prompts, task_id,
                                           raw=b"{" if mode == "malformed" else None, outcome=outcome))
    helpers.save(source, experiment, rows)
    imported = import_source(source, export, tmp_path)
    report = run_scoring(export, imported, tmp_path / "scored")
    s = report["summary"]
    assert s["total_selected_tasks"] == 14 and s["strict_exact_match_accuracy"] == 0
    assert s["coverage_adjusted_mean_score"] == 0 and s["valid_response_count"] == 0
    assert s["exact_match_accuracy_among_valid"] is None and s["mean_partial_credit_score"] is None
    assert sum(report["pilot_outcome_counts"].values()) == 14
    if mode in {"missing", "failed"}:
        assert s["missing_response_count"] == 14 and s["response_format_compliance_rate"] is None
    elif mode == "abstained":
        assert s["abstention_count"] == 14 and s["response_format_compliance_rate"] == 1
    else:
        assert s["invalid_response_count"] == 14 and s["response_format_compliance_rate"] == 0


@pytest.mark.parametrize("outcome,raw", [("refused", b"Synthetic refusal."), ("truncated", b'{"task_id":'), ("completed", b"")])
def test_refusal_truncation_and_empty_output_remain_primary(perfect_source, export, helpers, prompts, tmp_path, outcome, raw):
    source, experiment, rows = perfect_source
    first = helpers.add_record(source, prompts, rows[0]["task_id"], raw=raw, outcome=outcome)
    later = helpers.add_record(source, prompts, rows[0]["task_id"], raw=b"later unscored answer", attempt=2)
    helpers.save(source, experiment, [first, later, *rows[1:]])
    imported = import_source(source, export, tmp_path)
    report = run_scoring(export, imported, tmp_path / "scored")
    first_result = report["tasks"][0]
    assert first_result["primary_attempt"] == 1
    assert first_result["execution_outcome"] == outcome
    assert first_result["pilot_outcome"] == "INVALID_SCORED_RESPONSE"
    assert first_result["score"] == 0
    assert report["summary"]["missing_response_count"] == 0


def test_failed_execution_then_eligible_completion(perfect_source, export, helpers, prompts, tmp_path):
    source, experiment, rows = perfect_source
    original = (source / rows[0]["response"]["path"]).read_bytes()
    first = helpers.add_record(source, prompts, rows[0]["task_id"], raw=None, outcome="execution_failure")
    second = helpers.add_record(source, prompts, rows[0]["task_id"], raw=original, attempt=2)
    helpers.save(source, experiment, [first, second, *rows[1:]])
    imported = import_source(source, export, tmp_path)
    report = run_scoring(export, imported, tmp_path / "scored")
    assert report["summary"]["exact_response_count"] == 14
    assert report["tasks"][0]["primary_attempt"] == 2
    assert report["execution_accounting"]["total_execution_attempts"] == 15


@pytest.mark.parametrize("kind", ["duplicate", "unknown", "wrong_partition"])
def test_extraneous_source_records_block_pipeline(perfect_source, export, helpers, prompts, tmp_path, kind):
    source, experiment, rows = perfect_source
    if kind == "duplicate":
        rows.append(deepcopy(rows[0]))
    elif kind == "unknown":
        rows[0]["task_id"] = "UNKNOWN-SYNTHETIC-TASK"
    else:
        rows[0]["task_id"] = sorted(prompts.benchmark_ids - prompts.partition_ids)[0]
    helpers.save(source, experiment, rows)
    importer = load_script("import_bioeval_model_responses")
    with pytest.raises(importer.CaptureImportError):
        importer.run_import(export, source, tmp_path / "rejected")
    assert not (tmp_path / "rejected").exists()


def test_wrong_inner_task_id_remains_a_scientific_error(perfect_source, export, tmp_path):
    source, _, rows = perfect_source
    path = source / rows[0]["response"]["path"]
    answer = json.loads(path.read_bytes())
    answer["task_id"] = rows[1]["task_id"]
    raw = capture.canonical_json(answer)
    path.write_bytes(raw)
    rows[0]["response"]["sha256"] = capture.sha256(raw)
    (source / "records.json").write_bytes(capture.canonical_json({"capture_schema_version": "0.1", "records": rows}))
    imported = import_source(source, export, tmp_path)
    report = run_scoring(export, imported, tmp_path / "scored")
    assert report["tasks"][0]["score"] == 0
    assert "INCORRECT_TASK_ID" in report["tasks"][0]["error_codes"]


@pytest.mark.parametrize("kind", ["response", "metadata", "duplicate_submission", "wrong_partition_submission"])
def test_tampered_imports_are_rejected(imported, export, prompts, tmp_path, kind):
    _, package = imported
    if kind == "response":
        record = json.loads((package / "records.snapshot.json").read_bytes())["records"][0]
        (package / "evidence" / record["response"]["path"]).write_bytes(b"edited source")
    elif kind == "metadata":
        experiment = json.loads((package / "experiment.snapshot.json").read_bytes())
        experiment["model"]["name"] = ""
        rewrite_payload_and_hash(package, "experiment.snapshot.json", experiment)
    else:
        submission = json.loads((package / "submission.json").read_bytes())
        extra = deepcopy(submission["responses"][0])
        if kind == "wrong_partition_submission":
            extra["task_id"] = sorted(prompts.benchmark_ids - prompts.partition_ids)[0]
        submission["responses"].append(extra)
        rewrite_payload_and_hash(package, "submission.json", submission)
    with pytest.raises(ValueError):
        run_scoring(export, package, tmp_path / "rejected")
    assert not (tmp_path / "rejected").exists()


@pytest.mark.parametrize("kind", ["payload", "seed", "full_development"])
def test_wrong_or_corrupt_export_is_rejected(imported, export, tmp_path, kind):
    _, package = imported
    changed = tmp_path / "changed-export"
    if kind == "full_development":
        load_script("export_bioeval_model_prompts").export_prompts(BENCHMARK, "development", changed)
    elif kind == "seed":
        load_script("export_bioeval_model_prompts").export_prompts(BENCHMARK, "development", changed, pilot=True, seed=7)
    else:
        shutil.copytree(export, changed)
        path = changed / "prompts.json"
        path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError):
        run_scoring(changed, package, tmp_path / "rejected")


def test_no_leakage_and_offline_byte_replay(imported, export, tmp_path, monkeypatch):
    _, package = imported
    def blocked(*args, **kwargs):
        raise AssertionError("Network forbidden during synthetic pilot preflight")
    monkeypatch.setattr(socket.socket, "connect", blocked)
    monkeypatch.setattr(socket, "create_connection", blocked)
    first, second = tmp_path / "first", tmp_path / "second"
    report = run_scoring(export, package, first)
    run_scoring(export, package, second)
    runner = load_script("score_bioeval_pilot")
    runner.verify_pilot_outputs(first, report["provenance"])
    for path in first.iterdir():
        assert path.read_bytes() == (second / path.name).read_bytes()
    original = {t["task_id"]: t for t in json.loads((BENCHMARK / "tasks.json").read_bytes())["tasks"]}
    for record in json.loads((export / "prompts.json").read_bytes())["tasks"]:
        assert record == original[record["task_id"]]
        assert set(record) == {"task_id", "input", "question"}
    assert report["summary"]["total_selected_tasks"] == 14


@pytest.mark.parametrize("kind", ["file", "directory", "inside_import"])
def test_existing_output_protection(imported, export, tmp_path, kind):
    _, package = imported
    output = tmp_path / "output"
    if kind == "file":
        output.write_bytes(b"previous research")
    elif kind == "directory":
        output.mkdir()
    else:
        output = package / "derived"
    with pytest.raises((FileExistsError, ValueError)):
        run_scoring(export, package, output)
    if kind == "file":
        assert output.read_bytes() == b"previous research"


@pytest.mark.parametrize("failure", ["write", "rename"])
def test_publication_failure_cleanup(imported, export, tmp_path, monkeypatch, failure):
    _, package = imported
    if failure == "write":
        original = Path.write_bytes
        def fail(path, raw):
            if path.name == "report.json":
                raise OSError("synthetic pilot publication failure")
            return original(path, raw)
        monkeypatch.setattr(Path, "write_bytes", fail)
    else:
        def fail(*args):
            raise OSError("synthetic pilot publication failure")
        monkeypatch.setattr(Path, "rename", fail)
    with pytest.raises(OSError, match="synthetic pilot"):
        run_scoring(export, package, tmp_path / "output")
    assert not (tmp_path / "output").exists()
    assert not list(tmp_path.glob(".pilot-scoring-*"))


def test_output_checksum_tampering_is_rejected(imported, export, tmp_path):
    _, package = imported
    output = tmp_path / "output"
    run_scoring(export, package, output)
    path = output / "report.json"
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="checksum"):
        load_script("score_bioeval_pilot").verify_pilot_outputs(output)


def test_cli_scores_the_pilot(imported, export, tmp_path, monkeypatch, capsys):
    _, package = imported
    output = tmp_path / "output"
    monkeypatch.setattr(sys, "argv", ["score_bioeval_pilot.py", "--prompt-export", str(export),
                                     "--capture-import", str(package), "--output", str(output)])
    load_script("score_bioeval_pilot").main()
    summary = json.loads(capsys.readouterr().out)
    assert summary["submission_kind"] == "synthetic_scorer_testing"
    assert summary["summary"]["total_selected_tasks"] == 14


def test_audit_retry_times_must_be_chronological(perfect_source, helpers, prompts):
    source, experiment, rows = perfect_source
    later = helpers.add_record(source, prompts, rows[0]["task_id"], raw=b"secondary", attempt=2)
    rows[0]["timestamp"] = "2000-01-01T00:01:00Z"
    helpers.save(source, experiment, [rows[0], later, *rows[1:]])
    validation = capture.validate_capture(source, prompts)
    assert not validation.report["valid"]
    assert "INCONSISTENT_ATTEMPT_TIME" in {issue["code"] for issue in validation.report["issues"]}


def test_imported_evidence_root_cannot_escape_package(imported, prompts, tmp_path, monkeypatch):
    _, package = imported
    evidence = package / "evidence"
    outside = tmp_path / "outside"
    outside.mkdir()
    original = Path.resolve
    def resolve(path, *args, **kwargs):
        if path == evidence:
            return outside
        return original(path, *args, **kwargs)
    monkeypatch.setattr(Path, "resolve", resolve)
    with pytest.raises(ValueError, match="replay"):
        capture.load_verified_capture_import(package, prompts)
