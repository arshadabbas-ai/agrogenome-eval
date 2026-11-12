"""BioEval v0.1: deterministic public tasks with separately verified answers.

Task construction uses input facts only. Gold labels come from the existing
transcript engine; the independent standard-code lookup is an audit gate, not
an alternative production consequence classifier. No network or model APIs.
"""

from collections import Counter, defaultdict
from functools import lru_cache
import hashlib
import random
import re

from .transcript_variant import analyze_transcript_variant
from .variant_classifier import classify_variant


SCHEMA_VERSION = "0.1"
GENERATOR_VERSION = "0.1.0"
DEFAULT_SEED = 20261112
TARGETS = {"SYNONYMOUS": 12, "MISSENSE": 12, "STOP_GAINED": 12,
           "FIVE_PRIME_UTR": 6, "THREE_PRIME_UTR": 6,
           "START_CODON": 6, "TERMINAL_STOP_CODON": 6}
RESPONSE_FIELDS = ("task_id", "transcript_accession", "transcript_region", "cds_position",
                   "codon_number", "reference_codon", "alternative_codon", "reference_amino_acid",
                   "alternative_amino_acid", "molecular_consequence", "analysis_status",
                   "unsupported_region", "phenotype_claim_support")
TASK_INPUT_FIELDS = {"transcript_accession", "scientific_name", "protein_accession", "refseq_status",
                     "transcript_length", "position", "reference", "alternative", "cds_start", "cds_end",
                     "coordinate_convention", "coding_context", "variant_origin"}
_BASES = "TCAG"
_STANDARD_CODE = dict(zip([a+b+c for a in _BASES for b in _BASES for c in _BASES],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def substitution_key(variant: dict) -> tuple:
    """Return the exact transcript-position-reference-alternative identity."""
    required={"transcript_accession","position","reference","alternative"}
    if not isinstance(variant,dict) or not required <= variant.keys():
        raise ValueError("Missing substitution identity fields.")
    if not isinstance(variant["transcript_accession"],str) or not re.fullmatch(r"[A-Z]+_[0-9]+\.[0-9]+",variant["transcript_accession"]):
        raise ValueError("Substitution requires a versioned transcript accession.")
    if type(variant["position"]) is not int or variant["position"] < 1:
        raise ValueError("Substitution position must be a positive integer.")
    for field in ("reference","alternative"):
        value=variant[field]
        if not isinstance(value,str) or len(value)!=1 or value not in "ACGTacgt":
            raise ValueError("Substitution alleles must be single canonical DNA nucleotides.")
    if variant["reference"].upper()==variant["alternative"].upper():
        raise ValueError("Substitution alleles must differ.")
    return (variant["transcript_accession"], variant["position"],
            variant["reference"].upper(), variant["alternative"].upper())


def task_identifier(variant: tuple) -> str:
    """Opaque deterministic ID derived only from the substitution, never its label."""
    return "BE01-" + hashlib.sha256("|".join(map(str,variant)).encode()).hexdigest()[:16]


@lru_cache(maxsize=8)
def _enumerate(accession, transcript, start, end, excluded):
    pools = {category: [] for category in TARGETS}
    excluded = set(excluded)
    cds = transcript[start-1:end]
    for position, reference in enumerate(transcript, 1):
        for alternative in "ACGT":
            if alternative == reference: continue
            key = (accession, position, reference, alternative)
            if key in excluded: continue
            if position < start: category = "FIVE_PRIME_UTR"
            elif position > end: category = "THREE_PRIME_UTR"
            elif position < start+3: category = "START_CODON"
            elif position > end-3: category = "TERMINAL_STOP_CODON"
            else:
                category = classify_variant(cds, position-start+1, alternative)["Classification"]
            pools[category].append(key)
    return tuple((category, tuple(pools[category])) for category in TARGETS)


def enumerate_candidates(context: dict, excluded: set[tuple]) -> dict:
    """Validate the context then enumerate all eligible SNVs in position/ACGT order.

    Returns fresh lists; cached enumeration results are immutable. Enumeration
    uses the unchanged molecular classifier, never the independent audit table.
    """
    transcript, annotation = context["transcript"], context["annotation"]
    if not isinstance(transcript,str) or not transcript:
        raise ValueError("A nonempty validated transcript is required.")
    alternative = next(base for base in "ACGT" if base != transcript[0].upper())
    analyze_transcript_variant(transcript, annotation, 1, transcript[0], alternative)
    return {category: list(pool) for category, pool in _enumerate(
        annotation["transcript_accession"], transcript.upper(), annotation["start"], annotation["end"],
        tuple(sorted(excluded)))}


def construct_task(context: dict, position: int, reference: str, alternative: str, task_id: str) -> dict:
    """Construct a self-contained question from reference facts, not gold labels."""
    transcript, annotation = context["transcript"], context["annotation"]
    # Validate the tuple but do not use its classification to write the question.
    analyze_transcript_variant(transcript, annotation, position, reference, alternative)
    start, end = annotation["start"], annotation["end"]
    coding_context = None
    if start <= position <= end:
        relative = position-start
        offset = relative % 3
        first = start-1 + (relative//3)*3
        coding_context = {"reference_codon": transcript[first:first+3].upper(),
                          "within_codon_offset_1based": offset+1}
    provenance = annotation.get("provenance", {})
    inputs = {"transcript_accession": annotation["transcript_accession"],
        "scientific_name": provenance.get("scientific_name"), "protein_accession": annotation["protein_accession"],
        "refseq_status": provenance.get("refseq_status"), "transcript_length": len(transcript),
        "position": position, "reference": reference.upper(), "alternative": alternative.upper(),
        "cds_start": start, "cds_end": end, "coordinate_convention": "1-based inclusive transcript",
        "coding_context": coding_context, "variant_origin": "computationally_introduced_hypothetical_substitution"}
    context_text = (f"The reference coding triplet is {coding_context['reference_codon']}; the affected base has "
                    f"offset {coding_context['within_codon_offset_1based']} within that triplet (1-based)."
                    if coding_context else "No coding triplet is supplied; coding-codon interpretation may not apply.")
    question = (
        f"Task {task_id}: For {inputs['scientific_name']}, transcript {inputs['transcript_accession']} "
        f"({len(transcript)} nt; associated protein {inputs['protein_accession']}; RefSeq {inputs['refseq_status']}), "
        f"consider a computationally introduced hypothetical {inputs['reference']}>{inputs['alternative']} substitution "
        f"at transcript position {position}. This is not an observed natural allele. The complete positive-strand "
        f"CDS spans {start}..{end}, inclusive. All nucleotide positions are 1-based transcript coordinates, "
        f"not chromosome coordinates; codon_start=1, NCBI translation table=1. {context_text} "
        "Map the position to the annotated transcript region and, where applicable, CDS position and codon number. "
        "The classifier supports internal CDS codons only: SYNONYMOUS, MISSENSE or STOP_GAINED. "
        "It excludes the entire first coding codon and terminal stop codon; do not predict start-lost or stop-lost. "
        "UTR substitutions have no assessed protein consequence. Use ASSESSED, NOT_ASSESSED or UNSUPPORTED "
        "as appropriate; unsupported does not mean harmless. Use null for inapplicable fields. "
        "For excluded boundary codons, provide CDS position/codon number but null codon, amino-acid and molecular-consequence fields. "
        "Does this sequence evidence establish a drought-tolerance or yield claim? Use ESTABLISHED or NOT_ESTABLISHED "
        "for phenotype_claim_support according to the supplied evidence; do not infer expression, splicing or phenotype. "
        "Return exactly one JSON object with these fields: " + ", ".join(RESPONSE_FIELDS) + ".")
    return {"task_id": task_id, "input": inputs, "question": question}


def make_gold_answer(task: dict, context: dict) -> dict:
    """Obtain gold fields from the existing transcript-aware engine."""
    inputs = task["input"]
    result = analyze_transcript_variant(context["transcript"], context["annotation"],
        inputs["position"], inputs["reference"], inputs["alternative"])
    answer = {field: result.get(field) for field in RESPONSE_FIELDS}
    answer["task_id"] = task["task_id"]
    answer["phenotype_claim_support"] = "NOT_ESTABLISHED"
    validate_response(answer)
    return answer


def validate_response(answer: dict) -> None:
    """Enforce the v0.1 gold response contract (not free-text semantic judging)."""
    if not isinstance(answer, dict) or set(answer) != set(RESPONSE_FIELDS):
        raise ValueError("Response must have exactly the versioned contract fields.")
    if not isinstance(answer["task_id"], str) or not answer["task_id"]:
        raise ValueError("Response task_id must be a nonempty string.")
    if not isinstance(answer["transcript_accession"], str) or not re.fullmatch(r"[A-Z]+_[0-9]+\.[0-9]+",answer["transcript_accession"]):
        raise ValueError("Response must retain the version-specific transcript accession.")
    for key in RESPONSE_FIELDS:
        if key not in {"cds_position","codon_number"} and answer[key] is not None and not isinstance(answer[key],str):
            raise ValueError("Response text fields must be strings or null.")
    if answer["phenotype_claim_support"] != "NOT_ESTABLISHED":
        raise ValueError("The supplied benchmark evidence does not establish phenotype claims.")
    region, status = answer["transcript_region"], answer["analysis_status"]
    details = ("reference_codon", "alternative_codon", "reference_amino_acid", "alternative_amino_acid")
    if region in {"FIVE_PRIME_UTR", "THREE_PRIME_UTR"}:
        if status != "NOT_ASSESSED" or any(answer[k] is not None for k in (*details, "cds_position", "codon_number", "molecular_consequence", "unsupported_region")):
            raise ValueError("UTR answers cannot contain protein consequence assessments.")
    elif region == "CDS":
        if any(type(answer[k]) is not int or answer[k] < 1 for k in ("cds_position", "codon_number")):
            raise ValueError("Coding positions and codon numbers must be positive integers.")
        if status == "UNSUPPORTED":
            if answer["unsupported_region"] not in {"START_CODON", "TERMINAL_STOP_CODON"} or any(answer[k] is not None for k in (*details, "molecular_consequence")):
                raise ValueError("Boundary-codon consequences must remain unsupported.")
        elif status == "ASSESSED":
            if answer["molecular_consequence"] not in {"SYNONYMOUS", "MISSENSE", "STOP_GAINED"} or answer["unsupported_region"] is not None:
                raise ValueError("Invalid assessed molecular consequence.")
            for key in ("reference_codon", "alternative_codon"):
                value = answer[key]
                if not isinstance(value, str) or len(value) != 3 or set(value)-set("ACGT"):
                    raise ValueError("Assessed codons must be canonical DNA triplets.")
            for key in ("reference_amino_acid", "alternative_amino_acid"):
                allowed = set("ACDEFGHIKLMNPQRSTVWY" + ("*" if key == "alternative_amino_acid" else ""))
                if not isinstance(answer[key],str) or len(answer[key]) != 1 or answer[key] not in allowed:
                    raise ValueError("Invalid single-letter amino-acid symbol.")
        else: raise ValueError("Invalid CDS analysis status.")
    else: raise ValueError("Invalid transcript region.")


def verify_gold_independently(task: dict, answer: dict) -> None:
    """Block disagreements using direct indexing and an independent table-1 audit.

    This audit never silently changes engine output. Boundary codons are not
    assigned amino-acid consequences using an ordinary internal-codon lookup.
    """
    validate_response(answer)
    inputs = task["input"]
    start, end, position = inputs["cds_start"], inputs["cds_end"], inputs["position"]
    expected = {field: None for field in RESPONSE_FIELDS}
    expected.update(task_id=task["task_id"], transcript_accession=inputs["transcript_accession"],
                    phenotype_claim_support="NOT_ESTABLISHED")
    if position < start or position > end:
        expected.update(transcript_region="FIVE_PRIME_UTR" if position < start else "THREE_PRIME_UTR",
                        analysis_status="NOT_ASSESSED")
    else:
        cds_position = position-start+1
        codon_number = (cds_position-1)//3+1
        expected.update(transcript_region="CDS", cds_position=cds_position, codon_number=codon_number)
        if codon_number in {1, (end-start+1)//3}:
            expected.update(analysis_status="UNSUPPORTED", unsupported_region="START_CODON" if codon_number == 1 else "TERMINAL_STOP_CODON")
        else:
            context = inputs["coding_context"]
            old, offset = context["reference_codon"], (cds_position-1)%3
            if context["within_codon_offset_1based"] != offset+1 or old[offset] != inputs["reference"]:
                raise ValueError("Independent audit: reference base/codon offset disagreement.")
            new = old[:offset]+inputs["alternative"]+old[offset+1:]
            old_aa, new_aa = _STANDARD_CODE[old], _STANDARD_CODE[new]
            label = "STOP_GAINED" if new_aa == "*" else "SYNONYMOUS" if old_aa == new_aa else "MISSENSE"
            expected.update(analysis_status="ASSESSED", reference_codon=old, alternative_codon=new,
                            reference_amino_acid=old_aa, alternative_amino_acid=new_aa, molecular_consequence=label)
    if answer != expected:
        raise ValueError(f"Blocking independent gold-answer disagreement: {task['task_id']}.")


def _varied_sample(pool, count, start, end, rng):
    groups = defaultdict(list)
    for variant in pool:
        position = variant[1]
        group = (position-start)//3 if start <= position <= end else position
        groups[group].append(variant)
    keys = sorted(groups)
    rng.shuffle(keys)
    for key in keys: rng.shuffle(groups[key])
    selected = []
    while len(selected) < count:
        for key in keys:
            if groups[key]: selected.append(groups[key].pop())
            if len(selected) == count: break
    return selected


def answer_category(answer: dict) -> str:
    return answer["molecular_consequence"] or answer["unsupported_region"] or answer["transcript_region"]


def build_benchmark(contexts: dict, excluded: set[tuple], seed: int = DEFAULT_SEED) -> dict:
    """Build 60 balanced tasks per source, with exactly 2/3 development per stratum."""
    if type(seed) is not int or not contexts:
        raise ValueError("An integer seed and nonempty verified contexts are required.")
    rng = random.Random(seed)
    pools, shortages = {}, []
    for accession in sorted(contexts):
        if contexts[accession]["annotation"]["transcript_accession"] != accession:
            raise ValueError("Context key and annotation accession differ.")
        pools[accession] = enumerate_candidates(contexts[accession], excluded)
        for category, target in TARGETS.items():
            actual = len(pools[accession][category])
            if actual < target: shortages.append(f"{accession} {category}: eligible={actual}, required={target}, shortage={target-actual}")
    if shortages:
        raise ValueError("Insufficient eligible substitutions; revise distribution explicitly, do not duplicate: " + "; ".join(shortages))
    tasks, answers, splits = [], {}, {"development": [], "evaluation": []}
    counts, eligible = {}, {}
    for accession in sorted(contexts):
        context = contexts[accession]; annotation = context["annotation"]
        counts[accession] = dict(TARGETS)
        eligible[accession] = {category: len(pool) for category,pool in pools[accession].items()}
        for category, target in TARGETS.items():
            selected = _varied_sample(pools[accession][category], target, annotation["start"], annotation["end"], rng)
            rng.shuffle(selected)
            for index, variant in enumerate(selected):
                _, position, reference, alternative = variant
                task_id = task_identifier(variant)
                if task_id in answers: raise ValueError("Task-ID collision.")
                task = construct_task(context, position, reference, alternative, task_id)
                answer = make_gold_answer(task, context)
                verify_gold_independently(task, answer)
                if answer_category(answer) != category: raise ValueError("Candidate classification disagrees with gold answer.")
                tasks.append(task); answers[task_id] = answer
                splits["development" if index < target*2//3 else "evaluation"].append(task_id)
    result = {"tasks": {"schema_version": SCHEMA_VERSION, "tasks": sorted(tasks, key=lambda t:t["task_id"])},
              "gold_answers": {"schema_version": SCHEMA_VERSION, "answers": dict(sorted(answers.items()))},
              "splits": {"schema_version": SCHEMA_VERSION, **{k: sorted(v) for k,v in splits.items()}},
              "counts_by_accession_category": counts, "eligible_counts": eligible}
    validate_benchmark(result["tasks"], result["gold_answers"], result["splits"], contexts, excluded)
    return result


def validate_benchmark(tasks: dict, gold: dict, splits: dict, contexts: dict, excluded: set[tuple]) -> None:
    """Validate schema, evidence-linked tasks, all gold answers, balance and partitions."""
    for value, keys in [(tasks,{"schema_version","tasks"}), (gold,{"schema_version","answers"}),
                        (splits,{"schema_version","development","evaluation"})]:
        if not isinstance(value,dict) or set(value) != keys or value["schema_version"] != SCHEMA_VERSION:
            raise ValueError("Invalid benchmark envelope schema.")
    if not isinstance(tasks["tasks"], list) or not isinstance(gold["answers"], dict):
        raise ValueError("Invalid benchmark task/answer containers.")
    ids, tuples, counts = set(), set(), Counter()
    lookup = {}
    for task in tasks["tasks"]:
        if not isinstance(task,dict) or set(task) != {"task_id","input","question"} or not isinstance(task["input"],dict) or set(task["input"]) != TASK_INPUT_FIELDS:
            raise ValueError("Invalid task schema or leaked output fields.")
        task_id = task["task_id"]
        if not isinstance(task_id,str) or not task_id or task_id in ids:
            raise ValueError("Task IDs must be unique nonempty strings.")
        ids.add(task_id)
        inputs = task["input"]; key = substitution_key(inputs)
        if task_id != task_identifier(key):
            raise ValueError("Task ID must be the opaque deterministic substitution hash.")
        if key in tuples or key in excluded: raise ValueError("Duplicate or excluded demonstration substitution.")
        tuples.add(key)
        if key[0] not in contexts: raise ValueError("Unknown source accession/version.")
        if task != construct_task(contexts[key[0]], key[1], key[2], key[3], task_id):
            raise ValueError("Task question/context differs from source evidence.")
        if task_id not in gold["answers"]: raise ValueError("Missing gold answer.")
        answer = gold["answers"][task_id]
        verify_gold_independently(task, answer)
        category = answer_category(answer); counts[(key[0],category)] += 1
        lookup[task_id] = (key[0],category)
    if ids != set(gold["answers"]): raise ValueError("Gold-answer coverage differs from task IDs.")
    if [task["task_id"] for task in tasks["tasks"]] != sorted(ids):
        raise ValueError("Task ordering must be by opaque ID, not category.")
    expected_counts = Counter({(accession,category):target for accession in contexts for category,target in TARGETS.items()})
    if counts != expected_counts: raise ValueError("Species/category counts differ from target distribution.")
    if any(not isinstance(splits[name],list) or any(not isinstance(i,str) for i in splits[name]) for name in ("development","evaluation")):
        raise ValueError("Partition IDs must be lists of strings.")
    dev, evaluation = splits["development"], splits["evaluation"]
    if len(dev) != len(set(dev)) or len(evaluation) != len(set(evaluation)) or set(dev)&set(evaluation) or set(dev)|set(evaluation) != ids:
        raise ValueError("Partitions must be unique, disjoint and exhaustive.")
    for name, fraction in (("development",2),("evaluation",1)):
        actual = Counter(lookup[i] for i in splits[name])
        expected = Counter({(a,c):n*fraction//3 for a in contexts for c,n in TARGETS.items()})
        if actual != expected: raise ValueError("Partition stratification is inconsistent.")
