"""Offline benchmark QA with source-linked, independently derived expectations."""

from collections import Counter
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import shutil
import sys
from unittest.mock import patch
import xml.etree.ElementTree as ET

import pytest

import agrogenome_eval.bioeval_benchmark as benchmark


ROOT=Path(__file__).parents[1]
REFERENCE=ROOT/'data/reference'
CDS=ROOT/'data/derived/cds_protein'
VARIANTS=ROOT/'data/derived/variant_effects'
EXAMPLES=ROOT/'data/examples/hypothetical_variants.json'
OUTPUT=ROOT/'data/benchmarks/bioeval_v01'
EXPECTED_COUNTS={'SYNONYMOUS':12,'MISSENSE':12,'STOP_GAINED':12,'FIVE_PRIME_UTR':6,
                 'THREE_PRIME_UTR':6,'START_CODON':6,'TERMINAL_STOP_CODON':6}


def runner():
    scripts=str(ROOT/'scripts')
    if scripts not in sys.path:sys.path.insert(0,scripts)
    spec=importlib.util.spec_from_file_location('benchmark_runner',ROOT/'scripts/generate_bioeval_benchmark.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


@pytest.fixture(scope='module')
def sources():
    with patch('socket.socket.connect',side_effect=AssertionError('Network forbidden')):
        return runner().load_sources(REFERENCE,CDS,VARIANTS,EXAMPLES)


@pytest.fixture(scope='module')
def dataset(sources):
    contexts,excluded,_=sources
    return benchmark.build_benchmark(contexts,excluded)


def category(answer):
    return answer['molecular_consequence'] or answer['unsupported_region'] or answer['transcript_region']


def test_exact_counts_species_and_category(dataset):
    tasks=dataset['tasks']['tasks'];answers=dataset['gold_answers']['answers']
    assert len(tasks)==120 and len(answers)==120
    for accession in ['NM_001305944.1','XM_044542390.1']:
        rows=[t for t in tasks if t['input']['transcript_accession']==accession]
        assert len(rows)==60
        assert Counter(category(answers[t['task_id']]) for t in rows)==EXPECTED_COUNTS


def test_unique_ids_tuples_demo_exclusion_and_coverage(dataset,sources):
    tasks=dataset['tasks']['tasks'];ids=[t['task_id'] for t in tasks]
    keys=[benchmark.substitution_key(t['input']) for t in tasks]
    assert len(set(ids))==120 and len(set(keys))==120
    assert not set(keys)&sources[1]
    assert set(ids)==set(dataset['gold_answers']['answers'])


def test_stratified_disjoint_partitions(dataset):
    answers=dataset['gold_answers']['answers'];splits=dataset['splits']
    assert len(splits['development'])==80 and len(splits['evaluation'])==40
    assert not set(splits['development'])&set(splits['evaluation'])
    assert set(splits['development'])|set(splits['evaluation'])==set(answers)
    for name,numerator in [('development',2),('evaluation',1)]:
        counts=Counter((answers[i]['transcript_accession'],category(answers[i])) for i in splits[name])
        assert counts=={(a,c):n*numerator//3 for a in ['NM_001305944.1','XM_044542390.1'] for c,n in EXPECTED_COUNTS.items()}


@pytest.mark.parametrize('species,accession',[('maize','NM_001305944.1'),('wheat','XM_044542390.1')])
def test_all_source_linked_gold_independently(dataset,species,accession):
    # Read original FASTA/XML directly, bypassing benchmark/variant mappers.
    prefix=species+'_'+accession
    sequence=''.join((REFERENCE/(prefix+'.fasta')).read_text().splitlines()[1:])
    xml=ET.fromstring((CDS/(prefix+'.annotation.xml')).read_bytes())
    feature=next(f for f in xml.findall('.//GBFeature') if f.findtext('GBFeature_key')=='CDS')
    start,end=map(int,feature.findtext('GBFeature_location').split('..'))
    protein=''.join((CDS/(prefix+'.reference_protein.fasta')).read_text().splitlines()[1:])
    bases='TCAG'
    code=dict(zip([a+b+c for a in bases for b in bases for c in bases],
        'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
    answers=dataset['gold_answers']['answers'];seen=0
    for task in dataset['tasks']['tasks']:
        inputs=task['input']
        if inputs['transcript_accession']!=accession:continue
        seen+=1;p=inputs['position'];answer=answers[task['task_id']]
        assert inputs['reference']==sequence[p-1] and inputs['reference']!=inputs['alternative']
        assert (inputs['cds_start'],inputs['cds_end'])==(start,end)
        assert answer['phenotype_claim_support']=='NOT_ESTABLISHED'
        if p<start or p>end:
            assert answer['transcript_region']==('FIVE_PRIME_UTR' if p<start else 'THREE_PRIME_UTR')
            assert answer['analysis_status']=='NOT_ASSESSED' and answer['cds_position'] is None
            assert answer['molecular_consequence'] is None
            continue
        cds_pos=p-start+1;number=(cds_pos-1)//3+1
        assert answer['cds_position']==cds_pos and answer['codon_number']==number
        codon_first=start-1+3*(number-1);old=sequence[codon_first:codon_first+3];offset=(cds_pos-1)%3
        assert inputs['coding_context']=={'reference_codon':old,'within_codon_offset_1based':offset+1}
        if number in [1,(end-start+1)//3]:
            assert answer['analysis_status']=='UNSUPPORTED'
            assert answer['unsupported_region']==('START_CODON' if number==1 else 'TERMINAL_STOP_CODON')
            assert all(answer[k] is None for k in ['reference_codon','alternative_codon','reference_amino_acid','alternative_amino_acid','molecular_consequence'])
            continue
        new=old[:offset]+inputs['alternative']+old[offset+1:]
        aa,new_aa=code[old],code[new]
        assert aa==protein[number-1]
        expected='STOP_GAINED' if new_aa=='*' else 'SYNONYMOUS' if aa==new_aa else 'MISSENSE'
        assert (answer['reference_codon'],answer['alternative_codon'])==(old,new)
        assert (answer['reference_amino_acid'],answer['alternative_amino_acid'])==(aa,new_aa)
        assert answer['molecular_consequence']==expected and answer['analysis_status']=='ASSESSED'
    assert seen==60


def test_task_answer_separation_and_self_containment(dataset):
    for task in dataset['tasks']['tasks']:
        assert set(task)=={'task_id','input','question'}
        inputs=task['input'];question=task['question']
        assert not {'expected','gold','answer','molecular_consequence','analysis_status','transcript_region'}&set(inputs)
        assert 'gold_answers.json' not in question and 'expected classification' not in question.lower()
        assert all(str(inputs[key]) in question for key in ['transcript_accession','scientific_name','position'])
        assert all(label in question for label in ['SYNONYMOUS','MISSENSE','STOP_GAINED','ASSESSED','NOT_ASSESSED','UNSUPPORTED'])
        assert inputs['variant_origin']=='computationally_introduced_hypothetical_substitution'
        assert 'not an observed natural allele' in question


def test_question_construction_does_not_serialize_gold(sources,monkeypatch):
    def forbidden(*args):raise AssertionError('Gold serialization must not construct questions')
    monkeypatch.setattr(benchmark,'make_gold_answer',forbidden)
    context=sources[0]['NM_001305944.1']
    task=benchmark.construct_task(context,160,'C','A','independent-question')
    assert task['input']['reference']=='C' and 'gold' not in task


def test_codon_and_position_diversity(dataset):
    answers=dataset['gold_answers']['answers']
    for accession in ['NM_001305944.1','XM_044542390.1']:
        for name in ['SYNONYMOUS','MISSENSE','STOP_GAINED']:
            rows=[a for a in answers.values() if a['transcript_accession']==accession and a['molecular_consequence']==name]
            assert len({a['codon_number'] for a in rows})==12
        for region in ['FIVE_PRIME_UTR','THREE_PRIME_UTR']:
            positions={t['input']['position'] for t in dataset['tasks']['tasks'] if t['input']['transcript_accession']==accession and answers[t['task_id']]['transcript_region']==region}
            assert len(positions)==6


def synthetic_context():
    transcript='CCATGACGTTGTGAGGG'
    annotation=dict(transcript_accession='NM_000001.1',protein_accession='NP_000001.1',
        start=3,end=14,strand=1,codon_start=1,translation_table=1,reference_translation='MTL',
        coordinate_system='1-based inclusive',location='3..14')
    return {'transcript':transcript,'annotation':annotation,'species':'synthetic'}


@pytest.mark.parametrize('position,alternative,region,status,relative,codon',[
    (2,'A','FIVE_PRIME_UTR','NOT_ASSESSED',None,None),
    (3,'C','CDS','UNSUPPORTED',1,1),(5,'A','CDS','UNSUPPORTED',3,1),
    (6,'C','CDS','ASSESSED',4,2),(12,'A','CDS','UNSUPPORTED',10,4),
    (14,'C','CDS','UNSUPPORTED',12,4),(15,'A','THREE_PRIME_UTR','NOT_ASSESSED',None,None)])
def test_independent_synthetic_boundaries(position,alternative,region,status,relative,codon):
    context=synthetic_context();reference=context['transcript'][position-1]
    task=benchmark.construct_task(context,position,reference,alternative,'synthetic-case')
    answer=benchmark.make_gold_answer(task,context)
    assert (answer['transcript_region'],answer['analysis_status'],answer['cds_position'],answer['codon_number'])==(region,status,relative,codon)
    benchmark.verify_gold_independently(task,answer)


def test_exact_shortage_report_no_fabrication():
    context=synthetic_context()
    with pytest.raises(ValueError,match='STOP_GAINED: eligible=1, required=12, shortage=11'):
        benchmark.build_benchmark({'NM_000001.1':context},set())


@pytest.mark.parametrize('seed',[True,'20261112',None])
def test_invalid_seed(seed,sources):
    with pytest.raises(ValueError,match='integer seed'):benchmark.build_benchmark(sources[0],sources[1],seed)


@pytest.mark.parametrize('field,value',[
    ('phenotype_claim_support','ESTABLISHED'),('cds_position',True),
    ('reference_codon','NNN'),('alternative_amino_acid',[]),
    ('reference_amino_acid','*'),
    ('analysis_status','SUCCESS'),('unsupported_region','START_CODON')])
def test_strict_response_contract(dataset,field,value):
    answer=deepcopy(next(a for a in dataset['gold_answers']['answers'].values() if a['analysis_status']=='ASSESSED'))
    answer[field]=value
    with pytest.raises(ValueError):benchmark.validate_response(answer)


def test_extra_response_fields_are_rejected(dataset):
    answer=deepcopy(next(iter(dataset['gold_answers']['answers'].values())))
    answer['free_text']='unsupported narrative'
    with pytest.raises(ValueError):benchmark.validate_response(answer)


def test_independent_disagreement_blocks_generation(dataset):
    task=next(t for t in dataset['tasks']['tasks'] if dataset['gold_answers']['answers'][t['task_id']]['molecular_consequence']=='MISSENSE')
    answer=deepcopy(dataset['gold_answers']['answers'][task['task_id']]);answer['alternative_amino_acid']='W'
    if answer['reference_amino_acid']=='W':answer['alternative_amino_acid']='A'
    # Force a structurally valid but scientifically wrong value.
    answer['cds_position']+=3
    with pytest.raises(ValueError,match='Blocking independent'):benchmark.verify_gold_independently(task,answer)


@pytest.mark.parametrize('change',['missing_answer','overlap','extra_task_field','wrong_accession','duplicate_tuple'])
def test_dataset_schema_and_consistency_rejection(dataset,sources,change):
    data=deepcopy(dataset)
    task=data['tasks']['tasks'][0]
    if change=='missing_answer':del data['gold_answers']['answers'][task['task_id']]
    elif change=='overlap':data['splits']['evaluation'][0]=data['splits']['development'][0]
    elif change=='extra_task_field':task['expected']='leaked-label'
    elif change=='wrong_accession':task['input']['transcript_accession']='NM_001305944.2'
    else:
        other=data['tasks']['tasks'][1]
        other['input']=deepcopy(task['input'])
    with pytest.raises(ValueError):benchmark.validate_benchmark(data['tasks'],data['gold_answers'],data['splits'],sources[0],sources[1])


def generate(output):
    return runner().generate_benchmark(output,REFERENCE,CDS,VARIANTS,EXAMPLES)


def test_full_offline_byte_reproduction(tmp_path):
    with patch('socket.socket.connect',side_effect=AssertionError('Network forbidden')):
        manifest=generate(tmp_path/'first');generate(tmp_path/'second')
    assert manifest['task_count']==120 and manifest['partition_counts']=={'development':80,'evaluation':40}
    for path in OUTPUT.iterdir():
        assert path.read_bytes()==(tmp_path/'first'/path.name).read_bytes()==(tmp_path/'second'/path.name).read_bytes()


@pytest.mark.parametrize('name',['tasks.json','gold_answers.json','splits.json'])
def test_saved_output_corruption(tmp_path,sources,name):
    output=tmp_path/'benchmark';shutil.copytree(OUTPUT,output)
    path=output/name;path.write_bytes(path.read_bytes()+b'\n')
    with pytest.raises(ValueError,match='checksum'):runner().verify_benchmark_files(output,*sources)


@pytest.mark.parametrize('change',['category_count','seed','source_checksums'])
def test_inconsistent_manifest_rejected(tmp_path,sources,change):
    output=tmp_path/'benchmark';shutil.copytree(OUTPUT,output)
    path=output/'manifest.json';manifest=json.loads(path.read_bytes())
    if change=='category_count':manifest['counts_by_species_category']['maize']['MISSENSE']=999
    elif change=='seed':manifest['generation_seed']+=1
    else:manifest['source_checksums']={}
    from analyze_transcript_variants import encode
    path.write_bytes(encode(manifest))
    with pytest.raises(ValueError):runner().verify_benchmark_files(output,*sources)


def test_source_corruption_blocks_publication(tmp_path):
    reference=tmp_path/'reference';shutil.copytree(REFERENCE,reference)
    path=reference/'maize_NM_001305944.1.fasta';path.write_bytes(path.read_bytes()+b'A\n')
    with pytest.raises(ValueError):runner().generate_benchmark(tmp_path/'benchmark',reference,CDS,VARIANTS,EXAMPLES)
    assert not (tmp_path/'benchmark').exists()


def test_existing_output_preserved(tmp_path):
    marker=tmp_path/'research.txt';marker.write_text('preserve')
    with pytest.raises(FileExistsError):generate(tmp_path)
    assert marker.read_text()=='preserve'


@pytest.mark.parametrize('failure',['write','rename'])
def test_atomic_publication_failure(tmp_path,monkeypatch,failure):
    if failure=='rename':
        def fail(*args):raise OSError('simulated rename failure')
        monkeypatch.setattr(Path,'rename',fail)
    else:
        original=Path.write_bytes
        def fail(path,data):
            if path.name=='gold_answers.json':raise OSError('simulated write failure')
            return original(path,data)
        monkeypatch.setattr(Path,'write_bytes',fail)
    with pytest.raises(OSError,match='simulated'):generate(tmp_path/'benchmark')
    assert not (tmp_path/'benchmark').exists()


def test_answer_bearing_ids_rejected_even_with_consistent_links(dataset,sources):
    import hashlib
    for task in dataset['tasks']['tasks']:
        key=benchmark.substitution_key(task['input'])
        expected='BE01-'+hashlib.sha256('|'.join(map(str,key)).encode()).hexdigest()[:16]
        assert task['task_id']==expected
    data=deepcopy(dataset);task=data['tasks']['tasks'][0]
    old=task['task_id'];new='MISSENSE-ANSWER-LEAK'
    inputs=task['input']
    task.update(benchmark.construct_task(sources[0][inputs['transcript_accession']],
        inputs['position'],inputs['reference'],inputs['alternative'],new))
    answer=data['gold_answers']['answers'].pop(old);answer['task_id']=new
    data['gold_answers']['answers'][new]=answer
    for name in ['development','evaluation']:
        data['splits'][name]=[new if value==old else value for value in data['splits'][name]]
    with pytest.raises(ValueError,match='opaque deterministic'):
        benchmark.validate_benchmark(data['tasks'],data['gold_answers'],data['splits'],sources[0],sources[1])


@pytest.mark.parametrize('change',['species','disclosure','independent_validation'])
def test_manifest_provenance_and_disclosure_audited(tmp_path,sources,change):
    output=tmp_path/'benchmark';shutil.copytree(OUTPUT,output)
    path=output/'manifest.json';manifest=json.loads(path.read_bytes())
    if change=='species':manifest['sources'][0]['scientific_name']='Incorrect organism'
    elif change=='disclosure':manifest['public_benchmark_disclosure']='Secret evaluation holdout'
    else:manifest['independent_gold_validation']['all_tasks_verified']=False
    from analyze_transcript_variants import encode
    path.write_bytes(encode(manifest))
    with pytest.raises(ValueError):runner().verify_benchmark_files(output,*sources)


def test_wrong_codon_offset_blocks_independent_verification(dataset):
    task=deepcopy(next(t for t in dataset['tasks']['tasks'] if dataset['gold_answers']['answers'][t['task_id']]['analysis_status']=='ASSESSED'))
    task['input']['coding_context']['within_codon_offset_1based']=task['input']['coding_context']['within_codon_offset_1based']%3+1
    with pytest.raises(ValueError,match='offset disagreement'):
        benchmark.verify_gold_independently(task,dataset['gold_answers']['answers'][task['task_id']])


def test_gold_accession_must_be_versioned(dataset):
    answer=deepcopy(next(iter(dataset['gold_answers']['answers'].values())))
    answer['transcript_accession']=answer['transcript_accession'].split('.')[0]
    with pytest.raises(ValueError,match='version-specific'):benchmark.validate_response(answer)


def test_missing_gold_response_field_rejected(dataset):
    answer=deepcopy(next(iter(dataset['gold_answers']['answers'].values())))
    del answer['phenotype_claim_support']
    with pytest.raises(ValueError,match='exactly'):benchmark.validate_response(answer)


def test_wrong_amino_acid_rejected_by_independent_code(dataset):
    task=next(t for t in dataset['tasks']['tasks'] if dataset['gold_answers']['answers'][t['task_id']]['analysis_status']=='ASSESSED')
    answer=deepcopy(dataset['gold_answers']['answers'][task['task_id']])
    answer['alternative_amino_acid']=next(aa for aa in 'ACDEFGHIKLMNPQRSTVWY' if aa!=answer['alternative_amino_acid'])
    with pytest.raises(ValueError,match='Blocking independent'):benchmark.verify_gold_independently(task,answer)


def test_answer_metadata_cannot_be_embedded_in_task_inputs(dataset,sources):
    data=deepcopy(dataset);data['tasks']['tasks'][0]['input']['analysis_status']='ASSESSED'
    with pytest.raises(ValueError,match='leaked output fields'):
        benchmark.validate_benchmark(data['tasks'],data['gold_answers'],data['splits'],sources[0],sources[1])


def test_malformed_task_allele_has_explicit_schema_error(dataset,sources):
    data=deepcopy(dataset);data['tasks']['tasks'][0]['input']['reference']=None
    with pytest.raises(ValueError,match='single canonical DNA'):
        benchmark.validate_benchmark(data['tasks'],data['gold_answers'],data['splits'],sources[0],sources[1])
