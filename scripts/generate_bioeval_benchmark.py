"""Build BioEval v0.1 entirely offline; never replace existing benchmark outputs.

From the repository root:
    python scripts/generate_bioeval_benchmark.py
"""

import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from agrogenome_eval.bioeval_benchmark import (
    DEFAULT_SEED, GENERATOR_VERSION, SCHEMA_VERSION, build_benchmark,
    substitution_key, validate_benchmark,
)
from analyze_transcript_variants import digest, encode, load_verified_inputs, verify_outputs


PUBLIC_DISCLOSURE = "Both partitions AND all answers are public research materials. Evaluation is not hidden, secret or contamination-resistant."


def source_metadata(contexts: dict) -> list[dict]:
    """Retain identities and terminology from the verified input evidence."""
    return [{"transcript_accession":accession,**context["annotation"]["provenance"],
             "cds_start":context["annotation"]["start"],"cds_end":context["annotation"]["end"],
             "source_url":context["source_url"],"record_title":context["record_title"],
             "annotation_limitations":context["annotation_limitations"]}
            for accession,context in sorted(contexts.items())]


def load_sources(reference: Path, cds: Path, variant_evidence: Path, examples: Path):
    """Verify source bytes and source identities; return contexts and exclusion set."""
    verify_outputs(variant_evidence, reference, cds, examples)
    contexts, checksums = load_verified_inputs(reference, cds)
    expected = {"NM_001305944.1": ("Zea mays","NP_001292873.1",156,1259,"VALIDATED"),
                "XM_044542390.1": ("Triticum aestivum","XP_044398325.1",96,734,"MODEL")}
    if set(contexts) != set(expected): raise ValueError("Benchmark v0.1 requires the two pinned transcript references.")
    for accession, context in contexts.items():
        ann=context["annotation"]; prov=ann["provenance"]
        actual=(prov["scientific_name"],ann["protein_accession"],ann["start"],ann["end"],prov["refseq_status"])
        if actual != expected[accession]: raise ValueError("Pinned benchmark reference identity/annotation mismatch.")
    definitions=json.loads(examples.read_bytes())
    excluded={substitution_key(row) for row in definitions["variants"]}
    if len(excluded) != 14 or len(definitions["variants"]) != 14:
        raise ValueError("Expected exactly 14 distinct demonstration exclusions.")
    checksums={**checksums,"demonstrations/hypothetical_variants.json":digest(examples.read_bytes())}
    for path in sorted(variant_evidence.iterdir()):
        if path.is_file(): checksums["variant_evidence/"+path.name]=digest(path.read_bytes())
    return contexts,excluded,checksums


def verify_benchmark_files(output: Path, contexts: dict, excluded: set,
                           expected_source_checksums: dict) -> None:
    manifest=json.loads((output/"manifest.json").read_bytes())
    required={"tasks.json","gold_answers.json","splits.json"}
    if set(manifest.get("output_checksums",{})) != required:
        raise ValueError("Incomplete benchmark output checksum coverage.")
    for name,checksum in manifest["output_checksums"].items():
        if digest((output/name).read_bytes()) != checksum:
            raise ValueError(f"Benchmark checksum mismatch: {name}.")
    if manifest.get("schema_version") != SCHEMA_VERSION or manifest.get("generator_version") != GENERATOR_VERSION:
        raise ValueError("Benchmark version mismatch.")
    if manifest.get("source_checksums") != expected_source_checksums:
        raise ValueError("Benchmark source checksum mismatch.")
    if manifest.get("sources") != source_metadata(contexts):
        raise ValueError("Benchmark manifest source identity/provenance mismatch.")
    if manifest.get("public_benchmark_disclosure") != PUBLIC_DISCLOSURE:
        raise ValueError("Benchmark public-answer disclosure is inconsistent.")
    tasks=json.loads((output/"tasks.json").read_bytes())
    gold=json.loads((output/"gold_answers.json").read_bytes())
    splits=json.loads((output/"splits.json").read_bytes())
    validate_benchmark(tasks,gold,splits,contexts,excluded)
    independent=manifest.get("independent_gold_validation",{})
    assessed=sum(answer["analysis_status"]=="ASSESSED" for answer in gold["answers"].values())
    if independent.get("all_tasks_verified") is not True or independent.get("assessed_coding_tasks") != assessed:
        raise ValueError("Manifest independent-validation report disagrees with verified answers.")
    if manifest["task_count"] != len(tasks["tasks"]) or manifest["partition_counts"] != {k:len(splits[k]) for k in ("development","evaluation")}:
        raise ValueError("Benchmark manifest counts disagree with files.")
    replay=build_benchmark(contexts,excluded,manifest["generation_seed"])
    if any(actual != replay[key] for actual,key in [(tasks,"tasks"),(gold,"gold_answers"),(splits,"splits")]):
        raise ValueError("Benchmark files do not reproduce the recorded seed/generator rules.")
    species_counts={contexts[a]["species"]:v for a,v in replay["counts_by_accession_category"].items()}
    if (manifest["counts_by_accession_category"] != replay["counts_by_accession_category"]
            or manifest["counts_by_species_category"] != species_counts
            or manifest["eligible_counts_after_demo_exclusion"] != replay["eligible_counts"]):
        raise ValueError("Benchmark manifest category/eligibility counts are inconsistent.")


def generate_benchmark(output: Path, reference: Path, cds: Path, variant_evidence: Path,
                        examples: Path, seed: int = DEFAULT_SEED) -> dict:
    if output.exists(): raise FileExistsError(f"Preserving existing benchmark: {output}")
    contexts,excluded,checksums=load_sources(reference,cds,variant_evidence,examples)
    dataset=build_benchmark(contexts,excluded,seed)
    files={"tasks.json":encode(dataset["tasks"]),"gold_answers.json":encode(dataset["gold_answers"]),
           "splits.json":encode(dataset["splits"])}
    sources=source_metadata(contexts)
    manifest={"schema_version":SCHEMA_VERSION,"generator_version":GENERATOR_VERSION,"generation_seed":seed,
        "task_count":len(dataset["tasks"]["tasks"]),"sources":sources,"source_checksums":checksums,
        "source_roots":{"reference":"data/reference","cds_evidence":"data/derived/cds_protein",
                        "variant_evidence":"data/derived/variant_effects","demonstrations":"data/examples"},
        "counts_by_accession_category":dataset["counts_by_accession_category"],
        "counts_by_species_category":{contexts[a]["species"]:v for a,v in dataset["counts_by_accession_category"].items()},
        "eligible_counts_after_demo_exclusion":dataset["eligible_counts"],
        "partition_counts":{k:len(dataset["splits"][k]) for k in ("development","evaluation")},
        "independent_gold_validation":{"method":"Separate standard-code codon lookup and direct coordinate arithmetic; disagreements block generation.",
                                       "all_tasks_verified":True,"assessed_coding_tasks":72},
        "generation_rules":["Enumerate position ascending and alternatives A,C,G,T; exclude all 14 demo tuples.",
                            "Use a local random.Random(seed) after deterministic enumeration.",
                            "Round-robin shuffled codon/position groups to maximize diversity before reusing groups.",
                            "Stratify 2/3 development and 1/3 evaluation in each species/category.",
                            "Opaque task IDs hash substitution tuples; task order is sorted by ID, not category.",
                            "Gold answers are separate; questions use input facts and the same general rules for every category."],
        "public_benchmark_disclosure":PUBLIC_DISCLOSURE,
        "limitations":["120 correlated hypothetical tasks from only two transcripts; not broad plant-genome generalization.",
                       "Tasks can share genes, positions and codons, including across partitions; boundary categories necessarily share codons.",
                       "Reference codon context and classifier rules are supplied: this tests controlled interpretation, not sequence retrieval.",
                       "Public answers permit leakage/contamination; a genuinely private holdout requires separate future design.",
                       "No experimental natural alleles, phenotype, expression, splicing, crop-performance, function or orthology inference.",
                       "Future evaluation must report per-species/per-category metrics, not only overall accuracy.",
                       "No model APIs, automated LLM judging, web UI or evaluation scores are implemented.",
                       "Manifest hashes cover the other three output files; a manifest cannot contain its own byte hash.",
                       "NCBI source attribution and third-party-rights caveats remain in the verified reference-manifest snapshot."],
        "output_checksums":{name:digest(data) for name,data in files.items()}}
    output.parent.mkdir(parents=True,exist_ok=True)
    with TemporaryDirectory(dir=output.parent,prefix=".benchmark-") as temporary:
        staging=Path(temporary)/"verified";staging.mkdir()
        for name,data in files.items():(staging/name).write_bytes(data)
        (staging/"manifest.json").write_bytes(encode(manifest))
        verify_benchmark_files(staging,contexts,excluded,checksums)
        # Recheck physical source bytes immediately before publication.
        _,_,fresh_checksums=load_sources(reference,cds,variant_evidence,examples)
        if fresh_checksums != checksums: raise ValueError("Source evidence changed during generation.")
        if output.exists(): raise FileExistsError("Benchmark output appeared during generation.")
        staging.rename(output)
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path("data/benchmarks/bioeval_v01"))
    parser.add_argument("--seed",type=int,default=DEFAULT_SEED)
    args=parser.parse_args()
    try:
        result=generate_benchmark(args.output,Path("data/reference"),Path("data/derived/cds_protein"),
                                  Path("data/derived/variant_effects"),Path("data/examples/hypothetical_variants.json"),args.seed)
    except (ValueError,OSError) as exc: parser.exit(1,f"Benchmark generation stopped: {exc}\n")
    print(json.dumps({"tasks":result["task_count"],"partitions":result["partition_counts"],
                      "eligible":result["eligible_counts_after_demo_exclusion"]},indent=2))


if __name__=="__main__":main()
