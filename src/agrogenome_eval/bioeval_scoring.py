"""Deterministic structured-response scoring under bioeval-scoring-v0.1.

No coercion, free-text biological judgment, network requests or model calls.
Gold validation is deliberately distinct from permissive-but-typed model schema
validation: well-formed wrong biological assertions must remain diagnosable.
"""

from collections import Counter, defaultdict
import json
import math
import re

from Bio.Seq import Seq

from .bioeval_benchmark import RESPONSE_FIELDS, answer_category, verify_gold_independently


SCORING_SCHEMA_VERSION = "0.1"
SCORER_VERSION = "0.1.0"
POLICY_ID = "bioeval-scoring-v0.1"
_AA = set("ACDEFGHIKLMNPQRSTVWY")
_DETAILS = ("reference_codon", "alternative_codon", "reference_amino_acid", "alternative_amino_acid", "molecular_consequence")
_CODING = ("cds_position", "codon_number", *_DETAILS)
_ENUMS = {
    "transcript_region": {"FIVE_PRIME_UTR", "CDS", "THREE_PRIME_UTR"},
    "analysis_status": {"ASSESSED", "NOT_ASSESSED", "UNSUPPORTED"},
    "molecular_consequence": {None, "SYNONYMOUS", "MISSENSE", "STOP_GAINED"},
    "unsupported_region": {None, "START_CODON", "TERMINAL_STOP_CODON"},
    "phenotype_claim_support": {"ESTABLISHED", "NOT_ESTABLISHED"},
}
FIELD_CODES = {"task_id": "INCORRECT_TASK_ID", "transcript_accession": "INCORRECT_TRANSCRIPT_ACCESSION",
    "transcript_region": "INCORRECT_TRANSCRIPT_REGION", "cds_position": "INCORRECT_CDS_POSITION",
    "codon_number": "INCORRECT_CODON_NUMBER", "reference_codon": "INCORRECT_REFERENCE_CODON",
    "alternative_codon": "INCORRECT_ALTERNATIVE_CODON", "reference_amino_acid": "INCORRECT_REFERENCE_AMINO_ACID",
    "alternative_amino_acid": "INCORRECT_ALTERNATIVE_AMINO_ACID", "molecular_consequence": "INCORRECT_MOLECULAR_CONSEQUENCE",
    "analysis_status": "INCORRECT_ANALYSIS_STATUS", "unsupported_region": "INCORRECT_UNSUPPORTED_REGION",
    "phenotype_claim_support": "INCORRECT_PHENOTYPE_SUPPORT"}
PROFILES = {
    "ASSESSED": {"transcript_region":10, "analysis_status":15, "cds_position":10,
                 "codon_number":10, "reference_codon":10, "alternative_codon":10,
                 "reference_amino_acid":5, "alternative_amino_acid":5,
                 "molecular_consequence":20, "phenotype_claim_support":5},
    "NOT_ASSESSED": {"transcript_region":30, "analysis_status":30, "null_coding_fields":20,
                     "unsupported_region":5, "phenotype_claim_support":15},
    "UNSUPPORTED": {"transcript_region":15, "analysis_status":25, "cds_position":15,
                    "codon_number":10, "unsupported_region":20, "null_consequence_fields":10,
                    "phenotype_claim_support":5},
}


def _issue(code, field=None, message="", kind="FORMAT", derived=False):
    return {"code":code, "field":field, "message":message, "kind":kind, "derived":derived}


class ResponseFormatError(ValueError):
    """Raw JSON failure with stable machine-readable diagnostics."""
    def __init__(self, issues):
        self.issues = issues
        super().__init__("; ".join(issue["code"]+": "+issue["message"] for issue in issues))


def parse_model_response(raw_response: str) -> dict:
    """Parse exactly one JSON object; reject duplicates and non-finite numbers."""
    if not isinstance(raw_response,str):
        raise ResponseFormatError([_issue("INCORRECT_DATA_TYPE",message="Raw response must be a JSON string.")])
    def pairs(items):
        result={}
        for key,value in items:
            if key in result:
                raise ResponseFormatError([_issue("DUPLICATE_JSON_KEY",key,"Duplicate JSON key.")])
            result[key]=value
        return result
    def constant(value):
        raise ResponseFormatError([_issue("INVALID_JSON",message="Non-finite JSON constant: "+value)])
    def floating(value):
        number=float(value)
        if not math.isfinite(number):constant("overflow")
        return number
    try:
        result=json.loads(raw_response,object_pairs_hook=pairs,parse_constant=constant,parse_float=floating)
    except ResponseFormatError:raise
    except json.JSONDecodeError as exc:
        raise ResponseFormatError([_issue("INVALID_JSON",message=f"JSON decoding failed at line {exc.lineno}, column {exc.colno}.")]) from None
    except (ValueError,RecursionError):
        raise ResponseFormatError([_issue("INVALID_JSON",message="JSON exceeds supported decoding limits.")]) from None
    if not isinstance(result,dict):
        raise ResponseFormatError([_issue("INCORRECT_DATA_TYPE",message="JSON root must be an object.")])
    return result


def validate_model_response(response: dict) -> list[dict]:
    """Return structural schema errors without repairing or grading semantics."""
    if not isinstance(response,dict):
        return [_issue("INCORRECT_DATA_TYPE",message="Response must be an object.")]
    issues=[]
    for field in sorted(set(RESPONSE_FIELDS)-set(response)):
        issues.append(_issue("MISSING_REQUIRED_FIELD",field,"Required response field is missing."))
    for field in sorted(set(response)-set(RESPONSE_FIELDS),key=str):
        issues.append(_issue("UNEXPECTED_RESPONSE_FIELD",str(field),"Additional response fields are not permitted."))
    for field in RESPONSE_FIELDS:
        if field not in response:continue
        value=response[field]
        if field in {"task_id","transcript_accession"}:
            if not isinstance(value,str):issues.append(_issue("INCORRECT_DATA_TYPE",field,"Expected a nonempty string."))
            elif not value or (field=="transcript_accession" and not re.fullmatch(r"[A-Z]+_[0-9]+\.[0-9]+",value)):
                issues.append(_issue("INVALID_FIELD_VALUE",field,"Expected nonempty identity with an accession version where applicable."))
        elif field in {"cds_position","codon_number"}:
            if value is not None and type(value) is not int:
                issues.append(_issue("INCORRECT_DATA_TYPE",field,"Expected null or an integer; booleans/floats/strings are not integers."))
            elif value is not None and value<1:issues.append(_issue("INVALID_FIELD_VALUE",field,"Coordinates must be positive."))
        elif field in _ENUMS:
            if value is not None and not isinstance(value,str):
                issues.append(_issue("INCORRECT_DATA_TYPE",field,"Expected an enum string or permitted null."))
            elif value not in _ENUMS[field]:issues.append(_issue("INVALID_ENUM_VALUE",field,"Value is outside the response enum."))
        elif field in {"reference_codon","alternative_codon"}:
            if value is not None and not isinstance(value,str):issues.append(_issue("INCORRECT_DATA_TYPE",field,"Expected a DNA triplet string or null."))
            elif value is not None and (len(value)!=3 or set(value)-set("ACGT")):
                issues.append(_issue("INVALID_FIELD_VALUE",field,"Codons require exactly three uppercase A/C/G/T bases."))
        else:
            allowed=_AA|({"*"} if field=="alternative_amino_acid" else set())
            if value is not None and not isinstance(value,str):issues.append(_issue("INCORRECT_DATA_TYPE",field,"Expected a one-letter amino-acid string or null."))
            elif value is not None and (len(value)!=1 or value not in allowed):
                issues.append(_issue("INVALID_FIELD_VALUE",field,"Use a canonical uppercase one-letter amino acid; no notation conversion is performed."))
    return issues


def _base_result(task,gold):
    return {"task_id":task["task_id"], "scientific_metadata":{
        "scientific_name":task["input"]["scientific_name"], "transcript_accession":gold["transcript_accession"],
        "category":answer_category(gold), "analysis_status":gold["analysis_status"]},
        "submission_status":"MISSING", "response_entry_count":0, "format_valid":None,
        "strict_exact_match":False, "score":0, "raw_score":0, "components":{},
        "applied_caps":[], "field_matches":{field:None for field in RESPONSE_FIELDS},
        "parsed_response":None, "errors":[], "error_codes":[]}


def _finish(result):
    result["error_codes"]=sorted({issue["code"] for issue in result["errors"]})
    return result


def _invalid(task,gold,issues,count=1,status="INVALID"):
    result=_base_result(task,gold)
    result.update(submission_status=status,response_entry_count=count,format_valid=False,errors=issues)
    return _finish(result)


def _contradictions(response,allow_internal_translation=True):
    reasons=[];gate_consequence=False
    if response["transcript_region"]!="CDS":
        if response["analysis_status"]!="NOT_ASSESSED" or any(response[k] is not None for k in (*_CODING,"unsupported_region")):
            reasons.append("Claimed UTR region conflicts with coding/status/nullability fields.")
    elif response["analysis_status"]=="NOT_ASSESSED":
        reasons.append("Claimed CDS region conflicts with NOT_ASSESSED status.")
    elif response["analysis_status"]=="UNSUPPORTED":
        if response["unsupported_region"] is None or any(response[k] is not None for k in _DETAILS):
            reasons.append("Claimed unsupported CDS conflicts with consequence/unsupported-region fields.")
    else:
        if response["unsupported_region"] is not None:
            reasons.append("An assessed CDS cannot simultaneously claim an unsupported region.")
        if any(response[k] is None for k in _DETAILS):
            reasons.append("Assessed consequence lacks its claimed codon/amino-acid details.");gate_consequence=True
        elif allow_internal_translation:
            for codon,amino in [("reference_codon","reference_amino_acid"),("alternative_codon","alternative_amino_acid")]:
                if str(Seq(response[codon]).translate(table=1))!=response[amino]:
                    reasons.append("Asserted codon and amino acid disagree under table 1.");gate_consequence=True
            old,new,label=response["reference_amino_acid"],response["alternative_amino_acid"],response["molecular_consequence"]
            if ((label=="SYNONYMOUS" and old!=new) or
                (label=="MISSENSE" and (old==new or new=="*")) or
                (label=="STOP_GAINED" and new!="*")):
                reasons.append("Asserted consequence conflicts with asserted amino-acid changes.");gate_consequence=True
    if response["cds_position"] is not None and response["codon_number"] is not None and response["codon_number"]!=(response["cds_position"]-1)//3+1:
        reasons.append("Claimed CDS position and codon number disagree.")
    return reasons,gate_consequence


def score_task_response(task: dict, gold: dict, raw_response: str | None,
                         *, abstain: bool = False) -> dict:
    """Score a response to one verified task; None denotes no submitted response."""
    verify_gold_independently(task,gold)
    result=_base_result(task,gold)
    if type(abstain) is not bool:raise ValueError("abstain must be boolean.")
    if abstain:
        if raw_response is not None:
            return _invalid(task,gold,[_issue("INVALID_SUBMISSION_ENTRY",message="Abstention must have null raw_response.")])
        result.update(submission_status="ABSTAINED",response_entry_count=1,format_valid=True,
                      errors=[_issue("EXPLICIT_ABSTENTION",message="Explicit abstention; score zero.",kind="COVERAGE")])
        return _finish(result)
    if raw_response is None:
        result["errors"]=[_issue("MISSING_RESPONSE",message="No response submitted.",kind="COVERAGE")]
        return _finish(result)
    try:response=parse_model_response(raw_response)
    except ResponseFormatError as exc:return _invalid(task,gold,exc.issues)
    schema_errors=validate_model_response(response)
    if schema_errors:return _invalid(task,gold,schema_errors)
    result.update(submission_status="VALID",response_entry_count=1,format_valid=True,parsed_response=response)
    matches={field:response[field]==gold[field] for field in RESPONSE_FIELDS}
    result["field_matches"]=matches
    result["strict_exact_match"]=all(matches.values())
    result["errors"]=[_issue(FIELD_CODES[field],field,"Field differs from verified gold.",kind="SCIENTIFIC") for field in RESPONSE_FIELDS if not matches[field]]
    status=gold["analysis_status"]
    reasons,gate_consequence=_contradictions(response,allow_internal_translation=status=="ASSESSED")
    for name,weight in PROFILES[status].items():
        if name=="null_coding_fields":earned=all(response[k] is None for k in _CODING)
        elif name=="null_consequence_fields":earned=all(response[k] is None for k in _DETAILS)
        else:earned=matches[name]
        gated=name=="molecular_consequence" and gate_consequence
        result["components"][name]={"weight":weight,"score":weight if earned and not gated else 0,
                                    "matched":earned,"gated_by_contradiction":gated}
    caps=[]
    def cap(code,maximum,reason):caps.append({"code":code,"maximum":maximum,"reason":reason})
    if not matches["task_id"] or not matches["transcript_accession"]:cap("IDENTITY_MISMATCH",0,"Wrong task or transcript identity.")
    if not matches["transcript_region"]:cap("REGION_MISMATCH",25,"Incorrect transcript region.")
    if not matches["analysis_status"]:cap("STATUS_MISMATCH",25,"Incorrect analysis status.")
    if status=="ASSESSED" and not matches["molecular_consequence"]:cap("CONSEQUENCE_MISMATCH",50,"Incorrect assessed molecular consequence.")
    if status=="UNSUPPORTED" and not matches["unsupported_region"]:cap("UNSUPPORTED_REGION_MISMATCH",50,"Incorrect boundary-codon identifier.")
    if not matches["phenotype_claim_support"]:cap("PHENOTYPE_OVERCLAIM",50,"The evidence does not establish phenotype support.")
    if status=="NOT_ASSESSED" and any(response[k] is not None for k in (*_CODING,"unsupported_region")):
        cap("UTR_INAPPLICABLE_FIELDS",20,"UTR coding/consequence fields must remain null.")
        if any(response[k] is not None for k in _DETAILS):
            result["errors"].append(_issue("UNJUSTIFIED_UTR_PROTEIN_CONSEQUENCE",message="Protein/codon consequence asserted for a UTR task.",kind="EVIDENCE",derived=True))
        reasons.append("UTR answer contains inapplicable coding/consequence fields.")
    if status=="UNSUPPORTED" and any(response[k] is not None for k in _DETAILS):
        cap("BOUNDARY_OVERCLASSIFICATION",20,"Excluded boundary-codon consequences are not assessed.")
        result["errors"].append(_issue("UNSUPPORTED_BOUNDARY_CLASSIFICATION",message="Consequence asserted for an excluded boundary codon.",kind="EVIDENCE",derived=True))
    if reasons:
        result["errors"].append(_issue("INCONSISTENT_OUTPUT_FIELDS",message=" ".join(dict.fromkeys(reasons)),kind="EVIDENCE",derived=True))
        cap("OUTPUT_CONTRADICTION",50,"Internally contradictory or category-inapplicable assertions.")
    result["applied_caps"]=caps
    result["raw_score"]=sum(component["score"] for component in result["components"].values())
    result["score"]=min([result["raw_score"],*[item["maximum"] for item in caps]])
    return _finish(result)


def _ratio(numerator,denominator):return numerator/denominator if denominator else None


def _aggregate(results):
    n=len(results);submitted=sum(r["response_entry_count"]>0 for r in results)
    valid=sum(r["submission_status"]=="VALID" for r in results)
    abstained=sum(r["submission_status"]=="ABSTAINED" for r in results)
    invalid=sum(r["submission_status"] in {"INVALID","DUPLICATE"} for r in results)
    exact=sum(r["strict_exact_match"] for r in results);total=sum(r["score"] for r in results)
    frequencies=Counter(code for result in results for code in set(result["error_codes"]))
    field_rates={field:{"errors":sum(r["field_matches"][field] is False for r in results),
                        "denominator_valid_answers":valid,
                        "rate":_ratio(sum(r["field_matches"][field] is False for r in results),valid)} for field in RESPONSE_FIELDS}
    return {"total_selected_tasks":n,"selected_submitted_tasks":submitted,
        "selected_response_entry_count":sum(r["response_entry_count"] for r in results),
        "valid_response_count":valid,"invalid_response_count":invalid,"abstention_count":abstained,
        "missing_response_count":sum(r["submission_status"]=="MISSING" for r in results),
        "duplicate_task_count":sum(r["submission_status"]=="DUPLICATE" for r in results),
        "exact_response_count":exact,"score_sum":total,"response_coverage":_ratio(submitted,n),
        "scientific_answer_coverage":_ratio(valid,n),"strict_exact_match_accuracy":_ratio(exact,n),
        "valid_response_rate_among_submitted":_ratio(valid,submitted),
        "abstention_rate":_ratio(abstained,n),
        "missing_response_rate":_ratio(sum(r["submission_status"]=="MISSING" for r in results),n),
        "exact_match_accuracy_among_valid":_ratio(exact,valid),"exact_match_accuracy_among_submitted":_ratio(exact,submitted),
        "mean_partial_credit_score":_ratio(total,valid),"coverage_adjusted_mean_score":_ratio(total,n),
        "mean_score_among_submitted_tasks":_ratio(total,submitted),"response_format_compliance_rate":_ratio(valid+abstained,submitted),
        "field_error_rates":field_rates,"task_error_code_frequencies":dict(sorted(frequencies.items()))}


def summarize_evaluation(results: list[dict], submission_issues: list[dict] | None = None,
                           submitted_response_count: int | None = None) -> dict:
    """Aggregate fixed task denominators and gold-defined groups; undefined is null."""
    summary=_aggregate(results)
    summary["submitted_response_count"]=submitted_response_count if submitted_response_count is not None else summary["selected_response_entry_count"]
    issues=submission_issues or []
    summary["submission_issue_count"]=len(issues)
    summary["submission_issue_code_frequencies"]=dict(sorted(Counter(i["code"] for i in issues).items()))
    for group,key in [("per_species","scientific_name"),("per_category","category"),("per_analysis_status","analysis_status")]:
        buckets=defaultdict(list)
        for result in results:buckets[result["scientific_metadata"][key] or "UNSPECIFIED"].append(result)
        summary[group]={name:_aggregate(rows) for name,rows in sorted(buckets.items())}
    return summary


def evaluate_responses(tasks: dict, gold_answers: dict, submission: dict,
                         task_ids: list[str] | None = None) -> dict:
    """Evaluate versioned submission rows, rejecting duplicate selected IDs entirely."""
    required={"schema_version","submission_kind","responses"}
    if (not isinstance(submission,dict) or set(submission)!=required or submission["schema_version"]!="0.1"
        or not isinstance(submission["submission_kind"],str)
        or submission["submission_kind"] not in {"synthetic_scorer_testing","model_responses"}
        or not isinstance(submission["responses"],list)):
        raise ValueError("Invalid versioned response-submission envelope.")
    task_map={task["task_id"]:task for task in tasks["tasks"]};gold=gold_answers["answers"]
    if len(task_map)!=len(tasks["tasks"]) or set(task_map)!=set(gold):raise ValueError("Benchmark task/gold coverage mismatch.")
    selected=sorted(task_map if task_ids is None else task_ids)
    if len(selected)!=len(set(selected)) or not set(selected)<=set(task_map):raise ValueError("Invalid selected task IDs.")
    # Corrupt gold is a fatal evaluation-input error, even if all submissions
    # for a task are duplicates or absent; never turn it into a model penalty.
    for task_id in selected:verify_gold_independently(task_map[task_id],gold[task_id])
    groups=defaultdict(list);issues=[]
    for index,entry in enumerate(submission["responses"]):
        if not isinstance(entry,dict) or not isinstance(entry.get("task_id"),str) or not entry["task_id"]:
            issues.append(_issue("INVALID_SUBMISSION_ENTRY",message=f"Row {index} has no usable task ID.",kind="COVERAGE"));continue
        task_id=entry["task_id"]
        if task_id not in task_map:
            issues.append(_issue("UNKNOWN_TASK_ID",task_id,f"Row {index} refers to an unknown benchmark ID.",kind="COVERAGE"));continue
        if task_id not in selected:
            issues.append(_issue("OUT_OF_SCOPE_TASK_ID",task_id,f"Row {index} is outside the selected split.",kind="COVERAGE"));continue
        groups[task_id].append(entry)
    results=[]
    for task_id in selected:
        task,answer=task_map[task_id],gold[task_id];entries=groups[task_id]
        if len(entries)>1:
            results.append(_invalid(task,answer,[_issue("DUPLICATE_TASK_RESPONSE",message="All duplicate entries rejected; no best-answer selection.",kind="COVERAGE")],len(entries),"DUPLICATE"));continue
        if not entries:results.append(score_task_response(task,answer,None));continue
        entry=entries[0]
        if (set(entry)!={"task_id","raw_response","abstain"} or type(entry["abstain"]) is not bool
            or (entry["abstain"] and entry["raw_response"] is not None)
            or (not entry["abstain"] and not isinstance(entry["raw_response"],str))):
            results.append(_invalid(task,answer,[_issue("INVALID_SUBMISSION_ENTRY",message="Expected task_id, raw_response string and abstain=false, or null raw_response and abstain=true.")]));continue
        results.append(score_task_response(task,answer,entry["raw_response"],abstain=entry["abstain"]))
    return {"scoring_schema_version":SCORING_SCHEMA_VERSION,"scorer_version":SCORER_VERSION,
        "scoring_policy_id":POLICY_ID,"submission_kind":submission["submission_kind"],
        "summary":summarize_evaluation(results,issues,len(submission["responses"])),
        "submission_issues":issues,"tasks":results,
        "limitations":["Deterministic structured-response scoring, not free-text biological truth judgment.",
                       "Weights and caps are heuristic policy choices, not empirically calibrated probabilities.",
                       "Coordinate-only mistakes retain weighted partial credit; exact correctness still requires all fields.",
                       "Uniform templates and public answers permit shortcut learning and memorization.",
                       "Two correlated transcripts, artificial category balance and constant NOT_ESTABLISHED phenotype gold.",
                       "Both benchmark partitions and gold answers are public, not a private contamination-resistant holdout.",
                       "Scores do not establish general plant-genomics competence or phenotype reasoning.",
                       "Error frequencies count tags per task; correlated/derived tags are not independent biological failures."]}
