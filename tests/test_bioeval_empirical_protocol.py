"""Offline export checks against the frozen BioEval source records."""

from collections import Counter
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import shutil
import socket
import sys

import pytest

from agrogenome_eval.bioeval_benchmark import RESPONSE_FIELDS, TASK_INPUT_FIELDS


ROOT = Path(__file__).parents[1]
BENCHMARK = ROOT / "data/benchmarks/bioeval_v01"


@pytest.fixture(scope="module")
def exporter():
    scripts = str(ROOT / "scripts")
    if scripts not in sys.path:
        sys.path.insert(0, scripts)
    spec = importlib.util.spec_from_file_location(
        "empirical_prompt_exporter", ROOT / "scripts/export_bioeval_model_prompts.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def benchmark():
    return {name: json.loads((BENCHMARK / (name + ".json")).read_bytes())
            for name in ["tasks", "gold_answers", "splits"]}


@pytest.mark.parametrize("partition,pilot,count", [
    ("development", False, 80), ("evaluation", False, 40), ("development", True, 14),
])
def test_partition_and_pilot_selection(exporter, benchmark, partition, pilot, count):
    records, strata = exporter.select_tasks(
        benchmark["tasks"], benchmark["gold_answers"], benchmark["splits"], partition, pilot
    )
    ids = [record["task_id"] for record in records]
    assert len(ids) == len(set(ids)) == count
    assert ids == sorted(ids)
    assert set(ids) <= set(benchmark["splits"][partition])
    other = "evaluation" if partition == "development" else "development"
    assert not set(ids) & set(benchmark["splits"][other])
    if not pilot:
        assert ids == sorted(benchmark["splits"][partition])
        assert strata == []


def test_pilot_covers_each_species_category_once(exporter, benchmark):
    records, strata = exporter.select_tasks(
        benchmark["tasks"], benchmark["gold_answers"], benchmark["splits"], "development", True
    )
    answers = benchmark["gold_answers"]["answers"]
    coverage = Counter()
    for record in records:
        answer = answers[record["task_id"]]
        category = (answer["molecular_consequence"] or answer["unsupported_region"]
                    or answer["transcript_region"])
        coverage[(record["input"]["scientific_name"], category)] += 1
    categories = {"SYNONYMOUS", "MISSENSE", "STOP_GAINED", "FIVE_PRIME_UTR",
                  "THREE_PRIME_UTR", "START_CODON", "TERMINAL_STOP_CODON"}
    assert coverage == {(species, category): 1 for species in ["Zea mays", "Triticum aestivum"]
                        for category in categories}
    assert len(strata) == 14


def test_selection_seed_is_repeatable_without_global_random_changes(exporter, benchmark):
    args = (benchmark["tasks"], benchmark["gold_answers"], benchmark["splits"], "development", True)
    state = random.getstate()
    first = exporter.select_tasks(*args, seed=20261112)
    # Independently derived from the frozen sorted strata and seed, without
    # calling the exporter to construct these expected task identifiers.
    assert [record["task_id"] for record in first[0]] == [
        "BE01-4e9380e92f86e7db", "BE01-59a97a6dcb591090", "BE01-5d0827fb49185cb7",
        "BE01-6139da68064a293d", "BE01-6142241409dd428f", "BE01-800c18b07bd792aa",
        "BE01-9aafa550770ac853", "BE01-a731ca957135d22a", "BE01-a7d87278e83cba88",
        "BE01-ca3cbb251bcdcf28", "BE01-d1752fdd63a89b7d", "BE01-d95655fcb5a42333",
        "BE01-f493c674898043f2", "BE01-fc388045aa7ba8de",
    ]
    assert first == exporter.select_tasks(*args, seed=20261112)
    assert first != exporter.select_tasks(*args, seed=7)
    assert random.getstate() == state


def test_missing_pilot_stratum_is_an_explicit_failure(exporter, benchmark):
    splits = deepcopy(benchmark["splits"])
    answers = benchmark["gold_answers"]["answers"]
    splits["development"] = [i for i in splits["development"]
                             if answers[i]["molecular_consequence"] != "MISSENSE"]
    with pytest.raises(ValueError, match="14"):
        exporter.select_tasks(benchmark["tasks"], benchmark["gold_answers"], splits, "development", True)


@pytest.mark.parametrize("partition", ["all", "pilot", "invalid"])
def test_invalid_partition_is_rejected(exporter, tmp_path, partition):
    with pytest.raises(ValueError, match="Partition"):
        exporter.export_prompts(BENCHMARK, partition, tmp_path / "export")
    assert not (tmp_path / "export").exists()


@pytest.mark.parametrize("seed", [True, "20261112", 1.5])
def test_invalid_seed_is_rejected(exporter, tmp_path, seed):
    with pytest.raises(ValueError, match="integer"):
        exporter.export_prompts(BENCHMARK, "development", tmp_path / "export", pilot=True, seed=seed)


def test_pilot_cannot_select_evaluation_tasks(exporter, tmp_path):
    with pytest.raises(ValueError, match="development"):
        exporter.export_prompts(BENCHMARK, "evaluation", tmp_path / "export", pilot=True)


def test_export_is_question_preserving_and_gold_free(exporter, benchmark, tmp_path):
    output = tmp_path / "export"
    manifest = exporter.export_prompts(BENCHMARK, "development", output, pilot=True)
    payload = json.loads((output / "prompts.json").read_bytes())
    source = {record["task_id"]: record for record in benchmark["tasks"]["tasks"]}
    assert set(payload) == {"schema_version", "tasks"}
    assert payload["schema_version"] == "0.1"
    assert manifest["artifact_role"] == "researcher_provenance_only"
    assert manifest["pilot_sampling_seed"] == 20261112
    for record in payload["tasks"]:
        assert record == source[record["task_id"]]
        assert set(record) == {"task_id", "input", "question"}
        assert set(record["input"]) == TASK_INPUT_FIELDS
        assert not {"gold", "expected", "category", "score", "analysis_status", "molecular_consequence"} & record["input"].keys()
        assert all(field in record["question"] for field in RESPONSE_FIELDS)
    assert set(path.name for path in output.iterdir()) == {"prompts.json", "manifest.json"}


def test_byte_replay_provenance_and_offline_operation(exporter, tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Network forbidden during prompt export")
    monkeypatch.setattr(socket.socket, "connect", forbidden)
    monkeypatch.setattr(socket, "create_connection", forbidden)
    first = tmp_path / "first"
    second = tmp_path / "second"
    manifest = exporter.export_prompts(BENCHMARK, "development", first, pilot=True)
    exporter.export_prompts(BENCHMARK, "development", second, pilot=True)
    exporter.verify_prompt_export(first, manifest)
    for path in first.iterdir():
        assert path.read_bytes() == (second / path.name).read_bytes()
    assert manifest["output_checksums"]["prompts.json"] == hashlib.sha256((first / "prompts.json").read_bytes()).hexdigest()
    for name, checksum in manifest["benchmark_file_checksums"].items():
        assert checksum == hashlib.sha256((BENCHMARK / name).read_bytes()).hexdigest()
    records = json.loads((first / "prompts.json").read_bytes())["tasks"]
    for record in records:
        assert manifest["question_utf8_sha256"][record["task_id"]] == hashlib.sha256(record["question"].encode("utf-8")).hexdigest()
    assert "timestamp" not in manifest and manifest["benchmark_generation_seed"] == 20261112


@pytest.mark.parametrize("name", ["tasks.json", "gold_answers.json", "splits.json", "manifest.json"])
def test_corrupt_benchmark_is_rejected(exporter, tmp_path, name):
    benchmark = tmp_path / "benchmark"
    shutil.copytree(BENCHMARK, benchmark)
    path = benchmark / name
    if name == "manifest.json":
        value = json.loads(path.read_bytes())
        value["generation_seed"] = 0
        path.write_text(json.dumps(value), encoding="utf-8")
    else:
        path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError):
        exporter.export_prompts(benchmark, "development", tmp_path / "export", pilot=True)
    assert not (tmp_path / "export").exists()


@pytest.mark.parametrize("corrupt", ["payload", "question_hash"])
def test_output_corruption_is_rejected(exporter, tmp_path, corrupt):
    output = tmp_path / "export"
    exporter.export_prompts(BENCHMARK, "evaluation", output)
    if corrupt == "payload":
        path = output / "prompts.json"
        path.write_bytes(path.read_bytes() + b"\n")
    else:
        path = output / "manifest.json"
        value = json.loads(path.read_bytes())
        key = value["selected_task_ids"][0]
        value["question_utf8_sha256"][key] = "0" * 64
        path.write_text(json.dumps(value), encoding="utf-8")
    with pytest.raises(ValueError, match="checksum"):
        exporter.verify_prompt_export(output)


@pytest.mark.parametrize("kind", ["file", "directory"])
def test_existing_output_paths_are_preserved(exporter, tmp_path, kind):
    output = tmp_path / "export"
    if kind == "file":
        output.write_bytes(b"existing evidence")
        marker = output
    else:
        output.mkdir()
        marker = output / "research.txt"
        marker.write_bytes(b"existing evidence")
    with pytest.raises(FileExistsError):
        exporter.export_prompts(BENCHMARK, "development", output, pilot=True)
    assert marker.read_bytes() == b"existing evidence"


@pytest.mark.parametrize("failure", ["write", "rename"])
def test_failed_publication_cleans_staging(exporter, tmp_path, monkeypatch, failure):
    if failure == "write":
        original = Path.write_bytes
        def fail(path, data):
            if path.name == "prompts.json":
                raise OSError("simulated publication failure")
            return original(path, data)
        monkeypatch.setattr(Path, "write_bytes", fail)
    else:
        def fail(*args):
            raise OSError("simulated publication failure")
        monkeypatch.setattr(Path, "rename", fail)
    with pytest.raises(OSError, match="simulated"):
        exporter.export_prompts(BENCHMARK, "development", tmp_path / "export", pilot=True)
    assert list(tmp_path.iterdir()) == []


def test_benchmark_change_after_verification_blocks_export(exporter, tmp_path, monkeypatch):
    benchmark = tmp_path / "benchmark"
    shutil.copytree(BENCHMARK, benchmark)
    original = exporter.verify_benchmark_files
    def mutate(*args):
        original(*args)
        path = benchmark / "manifest.json"
        path.write_bytes(path.read_bytes() + b"\n")
    monkeypatch.setattr(exporter, "verify_benchmark_files", mutate)
    with pytest.raises(ValueError, match="changed during verification"):
        exporter.export_prompts(benchmark, "development", tmp_path / "export", pilot=True)
    assert not (tmp_path / "export").exists()


def test_cli_preserves_current_task_contract(exporter, tmp_path, monkeypatch, capsys):
    output = tmp_path / "cli"
    monkeypatch.setattr(sys, "argv", [
        "export_bioeval_model_prompts.py",
        "--partition", "evaluation", "--output", str(output),
    ])
    exporter.main()
    assert json.loads(capsys.readouterr().out)["task_count"] == 40
    assert json.loads((output / "prompts.json").read_bytes())["schema_version"] == "0.1"


def test_corrupt_source_evidence_blocks_export(exporter, tmp_path):
    project = tmp_path / "project"
    (project / "data").mkdir(parents=True)
    shutil.copytree(ROOT / "data/reference", project / "data/reference")
    shutil.copytree(ROOT / "data/derived", project / "data/derived")
    (project / "data/examples").mkdir()
    shutil.copy2(ROOT / "data/examples/hypothetical_variants.json",
                 project / "data/examples/hypothetical_variants.json")
    source = project / "data/reference/maize_NM_001305944.1.fasta"
    source.write_bytes(source.read_bytes() + b"A\n")
    with pytest.raises(ValueError, match="checksum"):
        exporter.export_prompts(BENCHMARK, "development", tmp_path / "export", project_root=project)
    assert not (tmp_path / "export").exists()
