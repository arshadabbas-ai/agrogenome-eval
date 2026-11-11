"""Independent synthetic and authentic hypothetical SNV tests; no live NCBI."""

from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
from unittest.mock import patch

import pytest

import agrogenome_eval.transcript_variant as engine


TRANSCRIPT = "CCATGACGTTGTGAGGG"
ROOT = Path(__file__).parents[1]
REFERENCE = ROOT / 'data/reference'
EVIDENCE = ROOT / 'data/derived/cds_protein'
EXAMPLES = ROOT / 'data/examples/hypothetical_variants.json'
OUTPUTS = ROOT / 'data/derived/variant_effects'


def annotation():
    return dict(transcript_accession='NM_000001.1', protein_accession='NP_000001.1',
        start=3, end=14, strand=1, codon_start=1, translation_table=1,
        coordinate_system='1-based inclusive', location='3..14', reference_translation='MTL',
        transcript_length=len(TRANSCRIPT),
        transcript_sequence_sha256=hashlib.sha256(TRANSCRIPT.encode()).hexdigest(),
        provenance=dict(transcript_accession='NM_000001.1', protein_accession='NP_000001.1',
                        gene_id='synthetic', scientific_name=None, refseq_status=None))


def analyze(position, alt, ann=None, transcript=TRANSCRIPT, reference=None):
    return engine.analyze_transcript_variant(transcript, annotation() if ann is None else ann,
        position, transcript[position-1] if reference is None else reference, alt)


@pytest.mark.parametrize('pos,alt,cds_pos,codon,old,new,aa,new_aa,classification', [
    (6,'C',4,2,'ACG','CCG','T','P','MISSENSE'),
    (8,'A',6,2,'ACG','ACA','T','T','SYNONYMOUS'),
    (10,'A',8,3,'TTG','TAG','L','*','STOP_GAINED'),
])
def test_independent_synthetic_coding(pos,alt,cds_pos,codon,old,new,aa,new_aa,classification):
    result=analyze(pos,alt)
    assert result['transcript_region']=='CDS' and result['analysis_status']=='ASSESSED'
    assert result['cds_position']==cds_pos and result['codon_number']==codon
    assert (result['reference_codon'],result['alternative_codon'])==(old,new)
    assert (result['reference_amino_acid'],result['alternative_amino_acid'])==(aa,new_aa)
    assert result['molecular_consequence']==classification
    assert json.loads(json.dumps(result))==result


@pytest.mark.parametrize('position,region', [(1,'FIVE_PRIME_UTR'),(2,'FIVE_PRIME_UTR'),
                                         (15,'THREE_PRIME_UTR'),(17,'THREE_PRIME_UTR')])
def test_utr_boundaries(position,region):
    alt=next(b for b in 'ACGT' if b!=TRANSCRIPT[position-1]);r=analyze(position,alt)
    assert r['transcript_region']==region and r['analysis_status']=='NOT_ASSESSED'
    for key in ['cds_position','codon_number','reference_codon','alternative_codon',
                'reference_amino_acid','alternative_amino_acid','molecular_consequence']:
        assert r[key] is None


@pytest.mark.parametrize('pos,cds_pos,codon,region', [(3,1,1,'START_CODON'),(4,2,1,'START_CODON'),
    (5,3,1,'START_CODON'),(12,10,4,'TERMINAL_STOP_CODON'),
    (13,11,4,'TERMINAL_STOP_CODON'),(14,12,4,'TERMINAL_STOP_CODON')])
def test_excluded_boundary_codons(pos,cds_pos,codon,region):
    alt=next(b for b in 'ACGT' if b!=TRANSCRIPT[pos-1]);r=analyze(pos,alt)
    assert r['cds_position']==cds_pos and r['codon_number']==codon
    assert r['unsupported_region']==region and r['analysis_status']=='UNSUPPORTED'
    assert r['molecular_consequence'] is None and r['alternative_amino_acid'] is None


def test_lowercase_normalization_and_missing_optional_provenance():
    ann=annotation();del ann['provenance']
    r=engine.analyze_transcript_variant(TRANSCRIPT.lower(),ann,6,'a','c')
    assert r['reference']=='A' and r['alternative']=='C'
    assert r['molecular_consequence']=='MISSENSE'
    assert r['gene_id'] is None and r['scientific_name'] is None and r['refseq_status'] is None


@pytest.mark.parametrize('position', [0,-1,18,True,False,1.5,'6',None])
def test_invalid_positions(position):
    with pytest.raises(ValueError,match='position'):
        engine.analyze_transcript_variant(TRANSCRIPT,annotation(),position,'A','C')


@pytest.mark.parametrize('allele', ['', 'AA','N','U','-',1,None,' A','\n'])
@pytest.mark.parametrize('field', ['reference','alternative'])
def test_invalid_alleles(allele,field):
    args=dict(transcript=TRANSCRIPT,annotation=annotation(),position=6,reference='A',alternative='C')
    args[field]=allele
    with pytest.raises(ValueError,match='nucleotide'): engine.analyze_transcript_variant(**args)


def test_reference_mismatch():
    with pytest.raises(ValueError,match='does not match'): analyze(6,'C',reference='T')


def test_identical_alleles():
    with pytest.raises(ValueError,match='must differ'): analyze(6,'a')


@pytest.mark.parametrize('transcript', ['', 'ATGN', 'ATG U', TRANSCRIPT[:-1]])
def test_invalid_transcript_or_length(transcript):
    with pytest.raises(ValueError):
        engine.analyze_transcript_variant(transcript,annotation(),1,'C','A')


@pytest.mark.parametrize('key', ['start','end','strand','codon_start','translation_table',
    'reference_translation','protein_accession','transcript_accession','location','coordinate_system'])
def test_missing_annotations(key):
    ann=annotation();del ann[key]
    with pytest.raises(ValueError,match='Missing'): analyze(6,'C',ann)


@pytest.mark.parametrize('ann', [None,[],{},'annotation'])
def test_invalid_annotation_types(ann):
    with pytest.raises(ValueError): engine.analyze_transcript_variant(TRANSCRIPT,ann,6,'A','C')


@pytest.mark.parametrize('key,value', [('transcript_accession','NM_000001'),('protein_accession','NP_000001'),
    ('strand',-1),('codon_start',2),('translation_table',2),('start',0),('end',100),
    ('coordinate_system','0-based'),('location','join(3..8,9..14)'),
    ('reference_translation','MPQ'),('transcript_length',100),('transcript_sequence_sha256','0'*64)])
def test_malformed_or_conflicting_annotation(key,value):
    ann=annotation();ann[key]=value
    with pytest.raises(ValueError): analyze(6,'C',ann)


@pytest.mark.parametrize('key,value', [('transcript_accession','NM_000001.2'),('protein_accession','NP_000001.2')])
def test_provenance_accession_versions(key,value):
    ann=annotation();ann['provenance'][key]=value
    with pytest.raises(ValueError,match='mismatch'): analyze(6,'C',ann)


def test_qualifier_conflict_and_exception():
    ann=annotation();ann['qualifiers']={'transl_table':['2']}
    with pytest.raises(ValueError,match='Conflicting'): analyze(6,'C',ann)
    ann['qualifiers']={'transl_except':['unsupported']}
    with pytest.raises(ValueError,match='Exceptional'): analyze(6,'C',ann)


def test_delegation_to_existing_classifier(monkeypatch):
    original=engine.classify_variant;calls=[]
    def spy(cds,position,alternative,table=1):
        calls.append((cds,position,alternative,table));return original(cds,position,alternative,table)
    monkeypatch.setattr(engine,'classify_variant',spy)
    analyze(6,'C');assert calls==[('ATGACGTTGTGA',4,'C',1)]
    analyze(1,'A');analyze(3,'C');assert len(calls)==1


def runner():
    scripts=str(ROOT/'scripts')
    if scripts not in sys.path:sys.path.insert(0,scripts)
    spec=importlib.util.spec_from_file_location('variant_runner',ROOT/'scripts/analyze_transcript_variants.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


# Independent table-1 codon lookup; the engine does not use this table.
BASES='TCAG'
CODE=dict(zip([a+b+c for a in BASES for b in BASES for c in BASES],
    'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
AUTHENTIC=json.loads(EXAMPLES.read_bytes())['variants']


@pytest.mark.parametrize('example',AUTHENTIC,ids=lambda e:e['id'])
def test_authentic_examples_independently(example):
    module=runner();contexts,_=module.load_verified_inputs(REFERENCE,EVIDENCE)
    context=contexts[example['transcript_accession']];tx=context['transcript'];ann=context['annotation']
    r=engine.analyze_transcript_variant(tx,ann,example['position'],example['reference'],example['alternative'])
    assert tx[example['position']-1]==example['reference']
    for key,value in example['expected'].items():assert r[key]==value
    assert example['origin']=='computationally_introduced_hypothetical_substitution'
    if r['analysis_status']=='ASSESSED':
        p=example['position']-ann['start']+1
        cds=tx[ann['start']-1:ann['end']];i=((p-1)//3)*3
        old=cds[i:i+3];offset=(p-1)%3;new=old[:offset]+example['alternative']+old[offset+1:]
        aa,new_aa=CODE[old],CODE[new]
        expected='STOP_GAINED' if new_aa=='*' else 'SYNONYMOUS' if aa==new_aa else 'MISSENSE'
        assert r['cds_position']==p and r['codon_number']==i//3+1
        assert (r['reference_codon'],r['alternative_codon'])==(old,new)
        assert (r['reference_amino_acid'],r['alternative_amino_acid'])==(aa,new_aa)
        assert r['molecular_consequence']==expected
    else:
        assert r['molecular_consequence'] is None


def test_deterministic_offline_replay(tmp_path):
    module=runner()
    with patch('socket.socket.connect',side_effect=AssertionError('Network forbidden')):
        contexts,_=module.load_verified_inputs(REFERENCE,EVIDENCE)
        assert module.select_hypothetical_examples(contexts)==json.loads(EXAMPLES.read_bytes())
        # Analysis/replay does not require a filesystem capable of hard links.
        # The optional atomic input-creation operation is tested separately.
        module.generate_results(EXAMPLES,tmp_path/'output',REFERENCE,EVIDENCE)
    module.verify_outputs(tmp_path/'output')
    for path in OUTPUTS.iterdir():assert path.read_bytes()==(tmp_path/'output'/path.name).read_bytes()


@pytest.mark.parametrize('target', ['reference','xml','analysis','snapshot'])
def test_corrupted_inputs_rejected(tmp_path,target):
    reference=tmp_path/'reference';derived=tmp_path/'derived'
    shutil.copytree(REFERENCE,reference);shutil.copytree(EVIDENCE,derived)
    choices={'reference':reference/'maize_NM_001305944.1.fasta',
             'xml':derived/'wheat_XM_044542390.1.annotation.xml',
             'analysis':derived/'maize_NM_001305944.1.analysis.json',
             'snapshot':derived/'reference_manifest.snapshot.json'}
    path=choices[target];path.write_bytes(path.read_bytes()+b'\n')
    with pytest.raises(ValueError,match='checksum'):
        runner().generate_results(EXAMPLES,tmp_path/'output',reference,derived)
    assert not (tmp_path/'output').exists()


@pytest.mark.parametrize('name', ['results.json','hypothetical_variants.snapshot.json'])
def test_corrupted_outputs_rejected(tmp_path,name):
    output=tmp_path/'output';shutil.copytree(OUTPUTS,output)
    path=output/name;path.write_bytes(path.read_bytes()+b'\n')
    with pytest.raises(ValueError,match='checksum'):runner().verify_outputs(output)


def test_existing_output_preserved(tmp_path):
    with pytest.raises(FileExistsError):runner().generate_results(EXAMPLES,tmp_path,REFERENCE,EVIDENCE)


@pytest.mark.parametrize('failure',['write','rename'])
def test_atomic_failure(tmp_path,monkeypatch,failure):
    if failure=='rename':
        def fail(*args):raise OSError('simulated rename failure')
        monkeypatch.setattr(Path,'rename',fail)
    else:
        original=Path.write_bytes
        def write(path,data):
            if path.name=='results.json':raise OSError('simulated write failure')
            return original(path,data)
        monkeypatch.setattr(Path,'write_bytes',write)
    with pytest.raises(OSError,match='simulated'):
        runner().generate_results(EXAMPLES,tmp_path/'output',REFERENCE,EVIDENCE)
    assert not (tmp_path/'output').exists()


@pytest.mark.parametrize('change',['version','reference','expectation','origin','duplicate','missing'])
def test_invalid_explicit_examples(tmp_path,change):
    definitions=json.loads(EXAMPLES.read_bytes());variant=definitions['variants'][0]
    if change=='version':variant['transcript_accession']='NM_001305944.2'
    elif change=='reference':variant['reference']='N'
    elif change=='expectation':variant['expected']['molecular_consequence']='MISSENSE'
    elif change=='origin':variant['origin']='natural'
    elif change=='duplicate':definitions['variants'].append(deepcopy(variant))
    else:del variant['position']
    path=tmp_path/'examples.json';path.write_bytes(runner().encode(definitions))
    with pytest.raises(ValueError):runner().generate_results(path,tmp_path/'output',REFERENCE,EVIDENCE)
    assert not (tmp_path/'output').exists()


def test_gene_id_conflicts_with_annotation_qualifier():
    ann=annotation();ann['qualifiers']={'db_xref':['GeneID:123']}
    ann['provenance']['gene_id']='999'
    with pytest.raises(ValueError,match='GeneID'):analyze(6,'C',ann)


def test_species_provenance_conflicts_with_annotation():
    ann=annotation();ann['organism']='Zea mays';ann['provenance']['scientific_name']='Triticum aestivum'
    with pytest.raises(ValueError,match='scientific_name mismatch'):analyze(6,'C',ann)


@pytest.mark.parametrize('key,value',[('gene_id',{}),('scientific_name',[]),('refseq_status',5)])
def test_malformed_provenance_fields(key,value):
    ann=annotation();ann['provenance'][key]=value
    with pytest.raises(ValueError,match='nonempty string'):analyze(6,'C',ann)


def test_duplicate_cds_evidence_accessions(tmp_path):
    derived=tmp_path/'derived';shutil.copytree(EVIDENCE,derived)
    path=derived/'summary.json';summary=json.loads(path.read_bytes())
    summary['records'].append(deepcopy(summary['records'][0]))
    path.write_bytes(runner().encode(summary))
    with pytest.raises(ValueError,match='Duplicate transcript'):
        runner().load_verified_inputs(REFERENCE,derived)


def test_incomplete_output_checksum_coverage(tmp_path):
    output=tmp_path/'output';shutil.copytree(OUTPUTS,output)
    path=output/'manifest.json';manifest=json.loads(path.read_bytes())
    del manifest['output_checksums']['results.json'];path.write_bytes(runner().encode(manifest))
    with pytest.raises(ValueError,match='coverage'):runner().verify_outputs(output)


def test_snapshot_checksum_is_linked_to_input_manifest(tmp_path):
    output=tmp_path/'output';shutil.copytree(OUTPUTS,output)
    path=output/'manifest.json';manifest=json.loads(path.read_bytes())
    manifest['input_checksums']['hypothetical_variants.json']='0'*64
    path.write_bytes(runner().encode(manifest))
    with pytest.raises(ValueError,match='snapshot checksum'):runner().verify_outputs(output)


def test_results_preserve_all_example_definitions(tmp_path):
    output=tmp_path/'output';shutil.copytree(OUTPUTS,output)
    path=output/'results.json';payload=json.loads(path.read_bytes());payload['results'].pop()
    raw=runner().encode(payload);path.write_bytes(raw)
    manifest_path=output/'manifest.json';manifest=json.loads(manifest_path.read_bytes())
    manifest['output_checksums']['results.json']=hashlib.sha256(raw).hexdigest()
    manifest_path.write_bytes(runner().encode(manifest))
    with pytest.raises(ValueError,match='preserve all'):runner().verify_outputs(output)


def test_output_input_checksums_rechecked_against_sources(tmp_path):
    output=tmp_path/'output';shutil.copytree(OUTPUTS,output)
    path=output/'manifest.json';manifest=json.loads(path.read_bytes())
    manifest['input_checksums']['reference/maize_NM_001305944.1.fasta']='0'*64
    path.write_bytes(runner().encode(manifest))
    with pytest.raises(ValueError,match='Source input checksums'):
        runner().verify_outputs(output,REFERENCE,EVIDENCE,EXAMPLES)


def test_existing_example_inputs_never_overwritten(tmp_path):
    path=tmp_path/'examples.json';path.write_bytes(b'existing input')
    with pytest.raises(FileExistsError):runner().save_example_inputs(path,{'variants':[]})
    assert path.read_bytes()==b'existing input'


@pytest.mark.parametrize('failure',['write','link'])
def test_example_input_publication_is_atomic(tmp_path,monkeypatch,failure):
    module=runner();path=tmp_path/'examples.json'
    if failure=='link':
        def fail(*args):raise OSError('simulated link failure')
        monkeypatch.setattr(module.os,'link',fail)
    else:
        original=Path.write_bytes
        def fail(staged,data):
            original(staged,data[:4]);raise OSError('simulated partial write')
        monkeypatch.setattr(Path,'write_bytes',fail)
    with pytest.raises(OSError,match='simulated'):module.save_example_inputs(path,{'variants':[]})
    assert not path.exists()


def test_example_creation_on_current_filesystem(tmp_path):
    import errno
    module=runner();path=tmp_path/'examples.json';definitions=json.loads(EXAMPLES.read_bytes())
    try:
        module.save_example_inputs(path,definitions)
    except OSError as exc:
        unsupported={errno.ENOTSUP,errno.EOPNOTSUPP,errno.ENOSYS,errno.EXDEV,errno.EPERM,errno.EACCES}
        if exc.errno not in unsupported and getattr(exc,'winerror',None) not in {1,50}:
            raise
        assert not path.exists()  # Safe refusal is valid on unsupported filesystems.
    else:
        assert path.read_bytes()==EXAMPLES.read_bytes()


def test_unsupported_hardlinks_preserve_existing_data(tmp_path,monkeypatch):
    import errno
    module=runner();path=tmp_path/'examples.json'
    def unsupported(*args):raise OSError(errno.ENOTSUP,'hard links unavailable')
    monkeypatch.setattr(module.os,'link',unsupported)
    with pytest.raises(OSError):module.save_example_inputs(path,{'variants':[]})
    assert not path.exists()
    path.write_bytes(b'existing scientific input')
    with pytest.raises(FileExistsError):module.save_example_inputs(path,{'variants':[]})
    assert path.read_bytes()==b'existing scientific input'
