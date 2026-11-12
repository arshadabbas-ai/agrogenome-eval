"""Independent arithmetic fixtures and offline synthetic scorer integration."""

from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import shutil
import sys
from unittest.mock import patch

import pytest
import agrogenome_eval.bioeval_scoring as scoring

from agrogenome_eval.bioeval_scoring import (
    ResponseFormatError, evaluate_responses, parse_model_response,
    score_task_response, summarize_evaluation, validate_model_response,
)


ROOT=Path(__file__).parents[1]
BENCHMARK=ROOT/'data/benchmarks/bioeval_v01'
MOCK=ROOT/'data/examples/bioeval_mock_responses_v01.json'
EXPECTATIONS=ROOT/'data/examples/bioeval_mock_expectations_v01.json'
DEMO=ROOT/'data/examples/derived/bioeval_scoring_v01'
AUDITED_DEMO=ROOT/'data/examples/derived/bioeval_scoring_v01_qa'


def fixture(kind='coding'):
    # Independent synthetic CDS: CC + ATG ACG TTG TGA + GGG. No crop identity.
    task={'task_id':'synthetic-'+kind,'input':{
        'scientific_name':'synthetic fixture (species unassigned)','transcript_accession':'NM_000001.1',
        'position':6,'reference':'A','alternative':'C','cds_start':3,'cds_end':14,
        'coding_context':{'reference_codon':'ACG','within_codon_offset_1based':1}}}
    gold={'task_id':task['task_id'],'transcript_accession':'NM_000001.1',
        'transcript_region':'CDS','cds_position':4,'codon_number':2,
        'reference_codon':'ACG','alternative_codon':'CCG','reference_amino_acid':'T','alternative_amino_acid':'P',
        'molecular_consequence':'MISSENSE','analysis_status':'ASSESSED','unsupported_region':None,
        'phenotype_claim_support':'NOT_ESTABLISHED'}
    if kind=='utr':
        task['input'].update(position=1,reference='C',alternative='A',coding_context=None)
        gold.update(transcript_region='FIVE_PRIME_UTR',cds_position=None,codon_number=None,
            reference_codon=None,alternative_codon=None,reference_amino_acid=None,alternative_amino_acid=None,
            molecular_consequence=None,analysis_status='NOT_ASSESSED')
    if kind=='boundary':
        task['input'].update(position=3,reference='A',alternative='C',coding_context={'reference_codon':'ATG','within_codon_offset_1based':1})
        gold.update(cds_position=1,codon_number=1,reference_codon=None,alternative_codon=None,
            reference_amino_acid=None,alternative_amino_acid=None,molecular_consequence=None,
            analysis_status='UNSUPPORTED',unsupported_region='START_CODON')
    return task,gold


def scored(kind='coding',**changes):
    task,gold=fixture(kind);answer={**gold,**changes}
    return score_task_response(task,gold,json.dumps(answer))


@pytest.mark.parametrize('kind',['coding','utr','boundary'])
def test_perfect_answers(kind):
    result=scored(kind)
    assert result['score']==100 and result['raw_score']==100
    assert result['strict_exact_match'] and result['format_valid']
    assert result['errors']==[] and result['applied_caps']==[]
    assert sum(c['weight'] for c in result['components'].values())==100


def test_terminal_boundary_is_not_translated_as_internal_codon():
    task,gold=fixture('boundary')
    task['input'].update(position=14,reference='A',alternative='C',coding_context={'reference_codon':'TGA','within_codon_offset_1based':3})
    gold.update(cds_position=12,codon_number=4,unsupported_region='TERMINAL_STOP_CODON')
    assert score_task_response(task,gold,json.dumps(gold))['score']==100


@pytest.mark.parametrize('changes,expected,code',[
    ({'cds_position':5},90,'INCORRECT_CDS_POSITION'),
    ({'codon_number':3},50,'INCORRECT_CODON_NUMBER'),
    ({'molecular_consequence':'SYNONYMOUS'},50,'INCORRECT_MOLECULAR_CONSEQUENCE'),
    ({'alternative_amino_acid':'A'},50,'INCORRECT_ALTERNATIVE_AMINO_ACID'),
    ({'reference_codon':'ATG'},50,'INCORRECT_REFERENCE_CODON'),
    ({'alternative_codon':'ACG'},50,'INCORRECT_ALTERNATIVE_CODON'),
    ({'reference_amino_acid':'A'},50,'INCORRECT_REFERENCE_AMINO_ACID'),
    ({'analysis_status':'NOT_ASSESSED'},25,'INCORRECT_ANALYSIS_STATUS'),
    ({'transcript_region':'THREE_PRIME_UTR'},25,'INCORRECT_TRANSCRIPT_REGION'),
    ({'phenotype_claim_support':'ESTABLISHED'},50,'INCORRECT_PHENOTYPE_SUPPORT'),
    ({'transcript_accession':'NM_000001.2'},0,'INCORRECT_TRANSCRIPT_ACCESSION'),
    ({'task_id':'different-task'},0,'INCORRECT_TASK_ID'),
    ({'unsupported_region':'START_CODON'},50,'INCORRECT_UNSUPPORTED_REGION'),
])
def test_independently_specified_coding_scores(changes,expected,code):
    result=scored(**changes)
    assert result['score']==expected and code in result['error_codes']
    assert result['format_valid'] and not result['strict_exact_match']
    assert result['score']==min([result['raw_score'],*[c['maximum'] for c in result['applied_caps']]])


def test_contradiction_withholds_accidentally_correct_label():
    result=scored(alternative_amino_acid='A')
    # 100 - 5 wrong amino acid - 20 gated consequence = 75, capped at 50.
    assert result['raw_score']==75 and result['score']==50
    assert result['components']['molecular_consequence']['matched']
    assert result['components']['molecular_consequence']['score']==0
    assert result['components']['molecular_consequence']['gated_by_contradiction']


def test_null_only_details_do_not_collect_coding_credit():
    result=scored(cds_position=None,codon_number=None,reference_codon=None,alternative_codon=None,
        reference_amino_acid=None,alternative_amino_acid=None,molecular_consequence=None)
    # Only correct region (10), status (15), phenotype (5).
    assert result['score']==30


def test_fully_wrong_but_structurally_valid_answer():
    result=scored(task_id='wrong',transcript_accession='NM_999999.2',transcript_region='THREE_PRIME_UTR',
        analysis_status='NOT_ASSESSED',cds_position=None,codon_number=None,reference_codon=None,
        alternative_codon=None,reference_amino_acid=None,alternative_amino_acid=None,
        molecular_consequence=None,phenotype_claim_support='ESTABLISHED')
    assert result['format_valid'] and result['score']==0


@pytest.mark.parametrize('kind,expected,code',[('utr',20,'UNJUSTIFIED_UTR_PROTEIN_CONSEQUENCE'),
                                            ('boundary',20,'UNSUPPORTED_BOUNDARY_CLASSIFICATION')])
def test_non_coding_or_excluded_protein_claims(kind,expected,code):
    result=scored(kind,reference_codon='AAA',alternative_codon='AAG',reference_amino_acid='K',
                  alternative_amino_acid='K',molecular_consequence='SYNONYMOUS')
    assert result['score']==expected and code in result['error_codes']
    assert result['format_valid']


def test_boundary_identifier_mismatch_cap():
    result=scored('boundary',unsupported_region='TERMINAL_STOP_CODON')
    assert result['raw_score']==80 and result['score']==50


def test_utr_wrong_status_cannot_receive_null_field_bonus_as_high_score():
    result=scored('utr',analysis_status='UNSUPPORTED')
    assert result['raw_score']==70 and result['score']==25


@pytest.mark.parametrize('raw,code',[
    ('not JSON','INVALID_JSON'), ('```json\n{}\n```','INVALID_JSON'), ('{} trailing','INVALID_JSON'),
    ('{"x":NaN}','INVALID_JSON'), ('{"x":Infinity}','INVALID_JSON'), ('{"x":1e400}','INVALID_JSON'),
    ('{"x":1,"x":2}','DUPLICATE_JSON_KEY'), ('{"outer":{"x":1,"x":2}}','DUPLICATE_JSON_KEY'),
    ('[]','INCORRECT_DATA_TYPE'), ('null','INCORRECT_DATA_TYPE')])
def test_strict_json_parsing(raw,code):
    with pytest.raises(ResponseFormatError) as exc:parse_model_response(raw)
    assert exc.value.issues[0]['code']==code


def test_valid_json_key_order_is_irrelevant():
    task,gold=fixture()
    raw=json.dumps(dict(reversed(list(gold.items()))))
    assert parse_model_response(raw)==gold
    assert score_task_response(task,gold,raw)['strict_exact_match']


@pytest.mark.parametrize('field,value,code',[
    ('cds_position',True,'INCORRECT_DATA_TYPE'),('cds_position',4.0,'INCORRECT_DATA_TYPE'),
    ('cds_position','4','INCORRECT_DATA_TYPE'),('codon_number',0,'INVALID_FIELD_VALUE'),
    ('alternative_amino_acid','Pro','INVALID_FIELD_VALUE'),('reference_amino_acid','*','INVALID_FIELD_VALUE'),
    ('reference_codon','acg','INVALID_FIELD_VALUE'),('analysis_status','SUCCESS','INVALID_ENUM_VALUE'),
    ('phenotype_claim_support',None,'INVALID_ENUM_VALUE')])
def test_schema_errors_score_zero_without_coercion(field,value,code):
    task,gold=fixture();answer={**gold,field:value}
    assert any(e['code']==code for e in validate_model_response(answer))
    result=score_task_response(task,gold,json.dumps(answer))
    assert result['score']==0 and not result['format_valid']
    assert all(value is None for value in result['field_matches'].values())
    assert not any(e['kind']=='SCIENTIFIC' for e in result['errors'])


def test_missing_and_extra_fields_are_explicit():
    task,gold=fixture();answer=deepcopy(gold);del answer['reference_codon'];answer['explanation']='This improves yield'
    result=score_task_response(task,gold,json.dumps(answer))
    assert result['error_codes']==['MISSING_REQUIRED_FIELD','UNEXPECTED_RESPONSE_FIELD']
    assert result['score']==0
    assert 'INCORRECT_PHENOTYPE_SUPPORT' not in result['error_codes']  # No free-text keyword judging.


def test_structured_phenotype_error_is_not_a_format_failure():
    result=scored(phenotype_claim_support='ESTABLISHED')
    assert result['format_valid'] and result['score']==50
    assert 'INVALID_ENUM_VALUE' not in result['error_codes']


def small_benchmark():
    a,ga=fixture();b,gb=fixture('utr')
    b['input']['scientific_name']='second synthetic grouping (not a crop)'
    return {'schema_version':'0.1','tasks':[a,b]}, {'schema_version':'0.1','answers':{a['task_id']:ga,b['task_id']:gb}}


def submission(rows):return {'schema_version':'0.1','submission_kind':'synthetic_scorer_testing','responses':rows}


def entry(task_id,answer):return {'task_id':task_id,'raw_response':json.dumps(answer),'abstain':False}


def test_missing_response_denominators_and_gold_groups():
    tasks,gold=small_benchmark();task_id=tasks['tasks'][0]['task_id']
    report=evaluate_responses(tasks,gold,submission([entry(task_id,gold['answers'][task_id])]))
    summary=report['summary']
    assert summary['total_selected_tasks']==2 and summary['missing_response_count']==1
    assert summary['response_coverage']==0.5 and summary['strict_exact_match_accuracy']==0.5
    assert summary['exact_match_accuracy_among_valid']==1
    assert summary['mean_partial_credit_score']==100 and summary['coverage_adjusted_mean_score']==50
    assert set(summary['per_category'])=={'MISSENSE','FIVE_PRIME_UTR'}
    assert len(summary['per_species'])==2 and set(summary['per_analysis_status'])=={'ASSESSED','NOT_ASSESSED'}
    assert summary['field_error_rates']['cds_position']=={'errors':0,'denominator_valid_answers':1,'rate':0}


def test_explicit_abstention_is_not_missing_or_valid_answer():
    tasks,gold=small_benchmark();task_id=tasks['tasks'][0]['task_id']
    report=evaluate_responses(tasks,gold,submission([{'task_id':task_id,'raw_response':None,'abstain':True}]))
    s=report['summary']
    assert s['abstention_count']==1 and s['missing_response_count']==1 and s['valid_response_count']==0
    assert s['response_format_compliance_rate']==1 and s['scientific_answer_coverage']==0
    assert s['mean_partial_credit_score'] is None and s['coverage_adjusted_mean_score']==0


def test_duplicate_responses_all_rejected():
    tasks,gold=small_benchmark();task_id=tasks['tasks'][0]['task_id'];row=entry(task_id,gold['answers'][task_id])
    result=evaluate_responses(tasks,gold,submission([row,deepcopy(row)]))
    s=result['summary']
    assert s['invalid_response_count']==1 and s['duplicate_task_count']==1
    assert s['submitted_response_count']==2 and s['selected_submitted_tasks']==1
    assert s['coverage_adjusted_mean_score']==0
    assert result['tasks'][0]['submission_status']=='DUPLICATE'


def test_unknown_and_out_of_scope_ids_are_reported():
    tasks,gold=small_benchmark();first,second=[t['task_id'] for t in tasks['tasks']]
    rows=[entry('unknown',gold['answers'][first]),entry(second,gold['answers'][second])]
    result=evaluate_responses(tasks,gold,submission(rows),[first])
    assert {e['code'] for e in result['submission_issues']}=={'UNKNOWN_TASK_ID','OUT_OF_SCOPE_TASK_ID'}
    assert result['summary']['total_selected_tasks']==1 and result['summary']['missing_response_count']==1
    assert result['summary']['response_coverage']==0


def test_malformed_entries_do_not_disappear():
    tasks,gold=small_benchmark();task_id=tasks['tasks'][0]['task_id']
    result=evaluate_responses(tasks,gold,submission([{'task_id':task_id},'unidentifiable']))
    assert result['summary']['invalid_response_count']==1 and result['summary']['submission_issue_count']==1
    assert 'INVALID_SUBMISSION_ENTRY' in result['tasks'][0]['error_codes']


def test_zero_denominators_are_null():
    tasks,gold=small_benchmark();report=evaluate_responses(tasks,gold,submission([]),[])
    s=report['summary']
    for field in ['response_coverage','strict_exact_match_accuracy','mean_partial_credit_score','coverage_adjusted_mean_score','response_format_compliance_rate']:
        assert s[field] is None
    json.dumps(report,allow_nan=False)


def test_error_frequencies_are_per_task_not_per_diagnostic():
    result=scored('utr',reference_codon='AAA',alternative_codon='AAG',reference_amino_acid='K',
        alternative_amino_acid='K',molecular_consequence='SYNONYMOUS')
    s=summarize_evaluation([result])
    assert s['task_error_code_frequencies']['INCONSISTENT_OUTPUT_FIELDS']==1
    assert s['field_error_rates']['molecular_consequence']['rate']==1
    assert any(e['derived'] for e in result['errors'])


def test_groups_use_gold_not_model_claims():
    result=scored(transcript_region='THREE_PRIME_UTR')
    summary=summarize_evaluation([result])
    assert set(summary['per_category'])=={'MISSENSE'}


def runner():
    scripts=str(ROOT/'scripts')
    if scripts not in sys.path:sys.path.insert(0,scripts)
    spec=importlib.util.spec_from_file_location('scoring_runner',ROOT/'scripts/score_bioeval_responses.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def test_demo_against_manually_specified_expectations():
    report=json.loads((DEMO/'report.json').read_bytes())
    expected=json.loads(EXPECTATIONS.read_bytes());lookup={row['task_id']:row for row in report['tasks']}
    for case in expected['cases']:
        actual=lookup[case['task_id']]
        assert actual['score']==case['expected_score'],case['case']
        assert set(case['required_error_codes'])<=set(actual['error_codes'])
    for key,value in expected['expected_summary'].items():assert report['summary'][key]==value
    assert report['summary']['mean_partial_credit_score']==53.75
    assert report['summary']['coverage_adjusted_mean_score']==pytest.approx(430/120)
    assert report['summary']['strict_exact_match_accuracy']==pytest.approx(1/120)
    assert set(report['summary']['per_species'])=={'Zea mays','Triticum aestivum'}
    assert len(report['summary']['per_category'])==7
    assert 'not an evaluation of a real AI model' in report['evaluation_description']


def test_offline_byte_reproducibility_and_provenance(tmp_path):
    module=runner()
    with patch('socket.socket.connect',side_effect=AssertionError('Network forbidden')):
        report=module.run_evaluation(MOCK,BENCHMARK,'all',tmp_path/'report')
    module.verify_evaluation_files(tmp_path/'report',report['provenance'])
    for path in AUDITED_DEMO.iterdir():assert path.read_bytes()==(tmp_path/'report'/path.name).read_bytes()
    assert report['provenance']['evaluation_split']=='all'
    assert report['provenance']['scoring_policy_id']=='bioeval-scoring-v0.1'
    assert 'timestamp' not in report['provenance']


@pytest.mark.parametrize('name',['report.json','submission.snapshot.json'])
def test_corrupted_scoring_output_rejected(tmp_path,name):
    output=tmp_path/'report';shutil.copytree(DEMO,output)
    path=output/name;path.write_bytes(path.read_bytes()+b'\n')
    with pytest.raises(ValueError,match='checksum'):runner().verify_evaluation_files(output)


def test_corrupt_benchmark_rejected_before_scoring(tmp_path):
    benchmark=tmp_path/'benchmark';shutil.copytree(BENCHMARK,benchmark)
    path=benchmark/'gold_answers.json';path.write_bytes(path.read_bytes()+b'\n')
    with pytest.raises(ValueError,match='checksum'):runner().run_evaluation(MOCK,benchmark,'all',tmp_path/'report')
    assert not (tmp_path/'report').exists()


def test_existing_report_is_preserved(tmp_path):
    marker=tmp_path/'research.txt';marker.write_text('preserve')
    with pytest.raises(FileExistsError):runner().run_evaluation(MOCK,BENCHMARK,'all',tmp_path)
    assert marker.read_text()=='preserve'


@pytest.mark.parametrize('failure',['write','rename'])
def test_report_publication_failure_cleanup(tmp_path,monkeypatch,failure):
    if failure=='rename':
        def fail(*args):raise OSError('simulated rename failure')
        monkeypatch.setattr(Path,'rename',fail)
    else:
        original=Path.write_bytes
        def write(path,data):
            if path.name=='report.json':raise OSError('simulated write failure')
            return original(path,data)
        monkeypatch.setattr(Path,'write_bytes',write)
    with pytest.raises(OSError,match='simulated'):runner().run_evaluation(MOCK,BENCHMARK,'all',tmp_path/'report')
    assert not (tmp_path/'report').exists()


def test_all_120_audited_gold_responses_score_exactly_100():
    tasks=json.loads((BENCHMARK/'tasks.json').read_bytes());gold=json.loads((BENCHMARK/'gold_answers.json').read_bytes())
    rows=[entry(task_id,answer) for task_id,answer in gold['answers'].items()]
    report=evaluate_responses(tasks,gold,submission(rows));summary=report['summary']
    assert summary['total_selected_tasks']==120 and summary['exact_response_count']==120
    assert summary['coverage_adjusted_mean_score']==100 and summary['strict_exact_match_accuracy']==1
    assert all(group['strict_exact_match_accuracy']==1 for group in summary['per_category'].values())


def test_invalid_submissions_cannot_improve_overall_denominator():
    tasks,gold=small_benchmark();first,second=[t['task_id'] for t in tasks['tasks']]
    rows=[entry(first,gold['answers'][first]),{'task_id':second,'raw_response':'broken JSON','abstain':False}]
    summary=evaluate_responses(tasks,gold,submission(rows))['summary']
    assert summary['invalid_response_count']==1 and summary['valid_response_count']==1
    assert summary['response_coverage']==1 and summary['response_format_compliance_rate']==0.5
    assert summary['mean_partial_credit_score']==100 and summary['coverage_adjusted_mean_score']==50
    assert summary['strict_exact_match_accuracy']==0.5


def test_corrupt_reference_provenance_stops_runner(tmp_path):
    project=tmp_path/'project';(project/'data/derived').mkdir(parents=True)
    (project/'data/examples').mkdir()
    shutil.copytree(ROOT/'data/reference',project/'data/reference')
    shutil.copytree(ROOT/'data/derived/cds_protein',project/'data/derived/cds_protein')
    shutil.copytree(ROOT/'data/derived/variant_effects',project/'data/derived/variant_effects')
    shutil.copy2(ROOT/'data/examples/hypothetical_variants.json',project/'data/examples/hypothetical_variants.json')
    path=project/'data/reference/maize_NM_001305944.1.fasta';path.write_bytes(path.read_bytes()+b'A\n')
    with pytest.raises(ValueError,match='checksum'):
        runner().run_evaluation(MOCK,BENCHMARK,'all',tmp_path/'report',project)
    assert not (tmp_path/'report').exists()


def test_boundary_overclassification_never_uses_internal_translation(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('Boundary effects cannot use internal-codon translation')
    monkeypatch.setattr(scoring,'Seq',forbidden)
    result=scored('boundary',analysis_status='ASSESSED',unsupported_region=None,
        reference_codon='ATG',alternative_codon='CTG',reference_amino_acid='M',
        alternative_amino_acid='M',molecular_consequence='SYNONYMOUS')
    assert result['score']<=20 and 'UNSUPPORTED_BOUNDARY_CLASSIFICATION' in result['error_codes']


def test_corrupt_gold_blocks_even_duplicate_submission_path():
    tasks,gold=small_benchmark();task_id=tasks['tasks'][0]['task_id']
    row=entry(task_id,gold['answers'][task_id])
    gold['answers'][task_id]['reference_amino_acid']='A'
    with pytest.raises(ValueError,match='Blocking independent'):
        evaluate_responses(tasks,gold,submission([row,deepcopy(row)]))


def test_audit_explicit_rates_and_empty_denominators():
    tasks,gold=small_benchmark();ids=[t['task_id'] for t in tasks['tasks']]
    s=evaluate_responses(tasks,gold,submission([
        {'task_id':ids[0],'raw_response':None,'abstain':True}]))['summary']
    assert s['valid_response_rate_among_submitted']==0
    assert s['abstention_rate']==0.5 and s['missing_response_rate']==0.5
    empty=evaluate_responses(tasks,gold,submission([]),[])['summary']
    for name in ['valid_response_rate_among_submitted','abstention_rate','missing_response_rate']:
        assert empty[name] is None


def test_audit_benchmark_changed_after_verification_is_rejected(tmp_path,monkeypatch):
    module=runner();benchmark=tmp_path/'benchmark';shutil.copytree(BENCHMARK,benchmark)
    original=module.verify_benchmark_files
    def changed_after_verification(*args):
        original(*args)
        path=benchmark/'manifest.json'
        value=json.loads(path.read_bytes());value['generation_seed']=0
        path.write_text(json.dumps(value),encoding='utf-8')
    monkeypatch.setattr(module,'verify_benchmark_files',changed_after_verification)
    with pytest.raises(ValueError,match='changed during verification'):
        module.run_evaluation(MOCK,benchmark,'all',tmp_path/'report')
    assert not (tmp_path/'report').exists()


def test_audit_escaped_duplicate_keys_and_parser_limits():
    for raw,code in [(r'{"task_id":1,"task_\u0069d":2}','DUPLICATE_JSON_KEY'),
                     ('   ','INVALID_JSON'),('{} {}','INVALID_JSON'),
                     ('{"n":'+'9'*5000+'}','INVALID_JSON')]:
        with pytest.raises(ResponseFormatError) as exc:parse_model_response(raw)
        assert exc.value.issues[0]['code']==code
    # Decoder depth thresholds differ by runtime; either a bounded decoding
    # failure or rejection of the non-object root is a safe contract failure.
    with pytest.raises(ResponseFormatError) as exc:
        parse_model_response('['*2000+'0'+']'*2000)
    assert exc.value.issues[0]['code'] in {'INVALID_JSON','INCORRECT_DATA_TYPE'}


def test_audit_manual_demo_arithmetic_and_group_denominators():
    report=json.loads((AUDITED_DEMO/'report.json').read_bytes())
    expected=json.loads(EXPECTATIONS.read_bytes())
    cases={c['task_id']:c for c in expected['cases']}
    # Separately specified arithmetic: 100+90+4*50+2*20 = 430.
    manual_scores={i:c['expected_score'] for i,c in cases.items()}
    assert sum(manual_scores.values())==100+90+4*50+2*20==430
    assert 8+3+1+108==120
    for task in report['tasks']:assert task['score']==manual_scores.get(task['task_id'],0)
    tasks=json.loads((BENCHMARK/'tasks.json').read_bytes())['tasks']
    gold=json.loads((BENCHMARK/'gold_answers.json').read_bytes())['answers']
    for grouping,metadata in [('per_species',lambda t:t['input']['scientific_name']),
                             ('per_category',lambda t:gold[t['task_id']]['molecular_consequence'] or gold[t['task_id']]['unsupported_region'] or gold[t['task_id']]['transcript_region']),
                             ('per_analysis_status',lambda t:gold[t['task_id']]['analysis_status'])]:
        for name,group in report['summary'][grouping].items():
            selected=[t for t in tasks if metadata(t)==name]
            total=sum(manual_scores.get(t['task_id'],0) for t in selected)
            assert group['total_selected_tasks']==len(selected)
            assert group['score_sum']==total
            assert group['coverage_adjusted_mean_score']==total/len(selected)


def test_audit_cap_precedence_and_coordinate_partial_credit():
    result=scored(task_id='wrong',transcript_region='FIVE_PRIME_UTR',
                  analysis_status='NOT_ASSESSED',phenotype_claim_support='ESTABLISHED')
    assert result['score']==0
    assert {c['maximum'] for c in result['applied_caps']}=={0,25,50}
    # Approved policy: coherent but wrong coordinates lose 10+10 points.
    result=scored(cds_position=7,codon_number=3)
    assert result['raw_score']==80 and result['score']==80
    assert not result['strict_exact_match']
    assert {'INCORRECT_CDS_POSITION','INCORRECT_CODON_NUMBER'}<=set(result['error_codes'])


@pytest.mark.parametrize('mode',['missing','malformed','abstained'])
def test_audit_uniform_non_answers_do_not_inflate_accuracy(mode):
    tasks,gold=small_benchmark()
    rows=[] if mode=='missing' else [
        {'task_id':t['task_id'],'raw_response':None if mode=='abstained' else 'invalid',
         'abstain':mode=='abstained'} for t in tasks['tasks']]
    s=evaluate_responses(tasks,gold,submission(rows))['summary']
    assert s['total_selected_tasks']==2 and s['coverage_adjusted_mean_score']==0
    assert s['strict_exact_match_accuracy']==0 and s['exact_match_accuracy_among_valid'] is None
    assert s['valid_response_count']+s['invalid_response_count']+s['abstention_count']+s['missing_response_count']==2


def test_audit_malformed_abstention_and_task_substitution():
    tasks,gold=small_benchmark();first,second=[t['task_id'] for t in tasks['tasks']]
    rows=[{'task_id':first,'raw_response':'{}','abstain':True},entry(second,gold['answers'][first])]
    report=evaluate_responses(tasks,gold,submission(rows))
    lookup={r['task_id']:r for r in report['tasks']}
    assert lookup[first]['submission_status']=='INVALID'
    assert lookup[second]['score']==0 and 'INCORRECT_TASK_ID' in lookup[second]['error_codes']


def test_audit_wrong_split_perfect_responses_cannot_improve_selected_score():
    tasks=json.loads((BENCHMARK/'tasks.json').read_bytes())
    gold=json.loads((BENCHMARK/'gold_answers.json').read_bytes())
    splits=json.loads((BENCHMARK/'splits.json').read_bytes())
    rows=[entry(i,gold['answers'][i]) for i in splits['development']]
    report=evaluate_responses(tasks,gold,submission(rows),splits['evaluation'])
    assert report['summary']['missing_response_count']==40
    assert report['summary']['coverage_adjusted_mean_score']==0
    assert len(report['submission_issues'])==80
    assert {i['code'] for i in report['submission_issues']}=={'OUT_OF_SCOPE_TASK_ID'}


def test_audit_duplicate_contradictory_answers_cannot_select_the_best():
    tasks,gold=small_benchmark();task_id=tasks['tasks'][0]['task_id']
    wrong={**gold['answers'][task_id],'phenotype_claim_support':'ESTABLISHED'}
    rows=[entry(task_id,gold['answers'][task_id]),entry(task_id,wrong)]
    report=evaluate_responses(tasks,gold,submission(rows))
    task=next(t for t in report['tasks'] if t['task_id']==task_id)
    assert task['submission_status']=='DUPLICATE' and task['score']==0
    assert report['summary']['invalid_response_count']==1


def test_audit_existing_file_output_is_preserved(tmp_path):
    output=tmp_path/'file';output.write_bytes(b'previous research')
    with pytest.raises(FileExistsError):runner().run_evaluation(MOCK,BENCHMARK,'all',output)
    assert output.read_bytes()==b'previous research'
