"""Fixed-population BioEval pilot scoring using the published scientific scorer."""

from collections import Counter

from .bioeval_benchmark import answer_category
from .bioeval_capture import CaptureValidation, VerifiedPromptExport, convert_to_scoring_submission
from .bioeval_scoring import evaluate_responses


PILOT_SCHEMA_VERSION = "0.1"
PILOT_EVALUATOR_VERSION = "0.1.0"
PILOT_SEED = 20261112
PILOT_SELECTION_ID = "aa134cfd42bfd1f9d87387f24cd3caf1d20a7c846ef405fc5df1798f78a2b0cb"
TASK_OUTCOMES = {
    "VALID": "VALID_SCORED_RESPONSE",
    "INVALID": "INVALID_SCORED_RESPONSE",
    "ABSTAINED": "EXPLICIT_ABSTENTION",
    "MISSING": "MISSING_OR_NO_ELIGIBLE_RESPONSE",
}


def verify_fixed_pilot(prompts: VerifiedPromptExport) -> None:
    """Require the exact original fourteen-task development sample and seed."""
    manifest = prompts.manifest
    if (manifest["stage"] != "pilot" or manifest["partition"] != "development"
            or manifest["pilot_sampling_seed"] != PILOT_SEED
            or len(prompts.records) != 14 or prompts.identity["selection_id"] != PILOT_SELECTION_ID
            or not set(prompts.records) <= prompts.partition_ids):
        raise ValueError("Pilot scoring requires the original 14-task development selection and seed 20261112.")


def evaluate_pilot(prompts: VerifiedPromptExport, validation: CaptureValidation,
                   tasks: dict, gold_answers: dict) -> dict:
    """Score fourteen verified IDs, retaining separate execution-attempt accounting.

    Callers must verify scientific inputs and replay the capture-import package
    before supplying these objects. No response is generated or repaired here.
    """
    verify_fixed_pilot(prompts)
    ids = sorted(prompts.records)
    if (not validation.report["valid"] or validation.report["selected_task_ids"] != ids
            or validation.experiment["prompt_export"] != prompts.identity):
        raise ValueError("Capture validation is incompatible with the fixed pilot.")
    task_map = {task["task_id"]: task for task in tasks["tasks"]}
    if any(task_map.get(task_id) != record for task_id, record in prompts.records.items()):
        raise ValueError("Pilot records differ from the verified benchmark task inputs.")
    categories = {"SYNONYMOUS", "MISSENSE", "STOP_GAINED", "FIVE_PRIME_UTR",
                  "THREE_PRIME_UTR", "START_CODON", "TERMINAL_STOP_CODON"}
    composition = Counter((prompts.records[i]["input"]["scientific_name"],
                           answer_category(gold_answers["answers"][i])) for i in ids)
    expected = {(species, category): 1 for species in ("Zea mays", "Triticum aestivum")
                for category in categories}
    if composition != expected:
        raise ValueError("Pilot must retain one task in every species/category cell.")
    submission = convert_to_scoring_submission(validation)
    result = evaluate_responses(tasks, gold_answers, submission, task_ids=ids)
    if result["submission_issues"]:
        raise ValueError("Unexpected extraneous submission rows after capture replay.")
    outcomes = Counter({name: 0 for name in TASK_OUTCOMES.values()})
    for task in result["tasks"]:
        status = task["submission_status"]
        if status not in TASK_OUTCOMES:
            raise ValueError("Unexpected task status after capture conversion.")
        task["pilot_outcome"] = TASK_OUTCOMES[status]
        outcomes[task["pilot_outcome"]] += 1
        event = validation.primary.get(task["task_id"])
        task["primary_attempt"] = event["attempt"] if event else None
        task["execution_outcome"] = event["outcome"] if event else None
        if event is None:
            task["capture_absence"] = ("NO_TASK_RECORD" if task["task_id"] in validation.report["missing_task_ids"]
                                       else "NO_ELIGIBLE_RESPONSE")
        else:
            task["capture_absence"] = None
    if sum(outcomes.values()) != 14 or result["summary"]["total_selected_tasks"] != 14:
        raise ValueError("Pilot outcome accounting must sum to fourteen.")
    attempts = validation.report["attempts"]
    execution_outcomes = Counter(event["outcome"] for event in attempts)
    eligible_count = sum(event["outcome"] in {"completed", "refused", "truncated", "abstained"}
                         for event in attempts)
    result["pilot_evaluation_schema_version"] = PILOT_SCHEMA_VERSION
    result["pilot_evaluator_version"] = PILOT_EVALUATOR_VERSION
    result["population"] = {
        "partition": "development", "stage": "pilot", "selection_id": PILOT_SELECTION_ID,
        "sampling_seed": PILOT_SEED, "selected_task_ids": ids, "selected_task_count": 14,
        "counts_by_species": {"Triticum aestivum": 7, "Zea mays": 7},
        "counts_by_category": {name: 2 for name in sorted(categories)},
    }
    result["pilot_outcome_counts"] = dict(sorted(outcomes.items()))
    result["execution_accounting"] = {
        "captured_task_count": len({event["task_id"] for event in attempts}),
        "captured_task_denominator": 14,
        "total_execution_attempts": len(attempts),
        "retry_attempt_count": sum(event["attempt"] > 1 for event in attempts),
        "eligible_attempt_count": eligible_count,
        "eligible_primary_responses": len(validation.primary),
        "secondary_eligible_attempt_count": eligible_count - len(validation.primary),
        "missing_task_record_count": len(validation.report["missing_task_ids"]),
        "recorded_tasks_without_eligible_response": len(validation.report["without_response_task_ids"]),
        "outcome_counts_all_attempts": dict(sorted(execution_outcomes.items())),
        "extraneous_record_count": 0,
    }
    result["experiment_provenance"] = {
        "experiment": validation.experiment,
        "provenance_confidence": validation.report["provenance_confidence"],
        "primary_attempt_policy": validation.experiment["retry_policy"]["selection"],
        "capture_issues": validation.report["issues"],
    }
    result["evaluation_description"] = (
        "Synthetic end-to-end pilot preflight. All observations are labeled synthetic capture-system testing."
        if submission["submission_kind"] == "synthetic_scorer_testing"
        else "Offline scoring of researcher-supplied captured observations with unverified model-origin attestation."
    )
    result["limitations"].extend([
        "Fourteen development tasks, one per species/category cell; descriptive and highly unstable percentages.",
        "A coordinate-only error may retain 90 points and a coherent wrong coordinate pair 80 under the heuristic policy.",
        "Execution attempts, primary task responses and valid structured answers have distinct counting units.",
        "Researcher metadata and byte checksums do not authenticate model origin or native-payload extraction.",
    ])
    return result
