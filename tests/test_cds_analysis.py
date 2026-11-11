"""Offline synthetic CDS tests plus pinned authentic annotation evidence."""

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from unittest.mock import Mock
import xml.etree.ElementTree as ET

import pytest

from agrogenome_eval.cds_analysis import (
    analyze_cds, extract_cds, parse_cds_annotation, parse_protein_fasta,
    translate_cds, utr_lengths,
)


ACCESSION = "NM_000001.1"
TRANSCRIPT = "CCATGAAATAAGGG"


def fixture_xml():
    root = ET.Element("GBSet")
    record = ET.SubElement(root, "GBSeq")
    for key, value in {"accession-version": ACCESSION, "moltype": "mRNA",
                       "sequence": TRANSCRIPT.lower(), "length": str(len(TRANSCRIPT)),
                       "definition": "synthetic transcript"}.items():
        ET.SubElement(record, "GBSeq_" + key).text = value
    feature = ET.SubElement(ET.SubElement(record, "GBSeq_feature-table"), "GBFeature")
    ET.SubElement(feature, "GBFeature_key").text = "CDS"
    ET.SubElement(feature, "GBFeature_location").text = "3..11"
    interval = ET.SubElement(ET.SubElement(feature, "GBFeature_intervals"), "GBInterval")
    for key, value in {"from": "3", "to": "11", "accession": ACCESSION}.items():
        ET.SubElement(interval, "GBInterval_" + key).text = value
    qualifiers = ET.SubElement(feature, "GBFeature_quals")
    for name, value in {"codon_start": "1", "transl_table": "1", "protein_id": "NP_000001.1",
                        "translation": "MK", "gene": "synthetic", "product": "synthetic peptide"}.items():
        q = ET.SubElement(qualifiers, "GBQualifier")
        ET.SubElement(q, "GBQualifier_name").text = name
        ET.SubElement(q, "GBQualifier_value").text = value
    return root


def parse(root=None):
    return parse_cds_annotation(ET.tostring(root if root is not None else fixture_xml()), ACCESSION, TRANSCRIPT)


def test_synthetic_extraction_translation_and_json():
    annotation = parse()
    assert extract_cds(TRANSCRIPT.lower(), annotation) == "ATGAAATAA"
    assert utr_lengths(TRANSCRIPT, annotation) == dict(five_prime_utr_length=2, three_prime_utr_length=3)
    protein = parse_protein_fasta(b">NP_000001.1 synthetic\nMK\n", "NP_000001.1")
    result = analyze_cds(TRANSCRIPT, annotation, protein)
    assert result['cds_length'] == 9 and result['protein_length'] == 2
    assert result['protein_sequence'] == 'MK'
    assert result['start_codon'] == 'ATG' and result['terminal_stop_codon'] == 'TAA'
    assert result['annotated_translation_match'] and result['reference_protein_match']
    assert result['validation_status'] == 'validated'
    assert json.loads(json.dumps(result)) == result


def test_annotation_only_not_claimed_as_protein_verified():
    result = analyze_cds(TRANSCRIPT, parse())
    assert result['validation_status'] == 'annotation_only'
    assert result['reference_protein_match'] is None


def test_no_utrs_for_full_length_cds():
    annotation = parse(); annotation.update(start=1, end=9)
    assert utr_lengths('ATGAAATAA', annotation) == dict(five_prime_utr_length=0, three_prime_utr_length=0)


def test_wrong_reference_protein_version():
    with pytest.raises(ValueError, match='accession/version'):
        analyze_cds(TRANSCRIPT, parse(), dict(accession='NP_000001.2', sequence='MK'))


@pytest.mark.parametrize('target', ['annotation', 'protein'])
def test_translation_mismatch_reported(target):
    annotation = parse()
    protein = dict(accession='NP_000001.1', sequence='MK')
    if target == 'annotation': annotation['reference_translation'] = 'MP'
    else: protein['sequence'] = 'MP'
    result = analyze_cds(TRANSCRIPT, annotation, protein)
    assert result['validation_status'] == 'mismatch'
    assert len(result['discrepancies']) == 1
    assert result['protein_sequence'] == 'MK'


@pytest.mark.parametrize('tag,value', [
    ('GBSeq_accession-version', 'NM_000001.2'), ('GBSeq_moltype', 'DNA'),
    ('GBSeq_sequence', 'aaaa'), ('GBSeq_length', '999'),
    ('GBInterval_accession', 'NM_999999.1'), ('GBInterval_from', 'bad'),
    ('GBInterval_to', '12'), ('GBFeature_location', '<3..11'),
    ('GBFeature_location', '3..>11'), ('GBFeature_location', 'join(3..5,6..11)'),
    ('GBFeature_location', 'complement(3..11)'),
])
def test_annotation_mismatch_or_unsupported_location(tag, value):
    root = fixture_xml(); root.find('.//' + tag).text = value
    with pytest.raises(ValueError): parse(root)


@pytest.mark.parametrize('name,value', [('codon_start', '2'), ('codon_start', '0'),
    ('transl_table', '2'), ('transl_table', 'bad'), ('protein_id', 'NP_000001'),
    ('translation', 'MX'), ('translation', 'MK*')])
def test_invalid_qualifiers(name, value):
    root = fixture_xml()
    for q in root.findall('.//GBQualifier'):
        if q.findtext('GBQualifier_name') == name: q.find('GBQualifier_value').text = value
    with pytest.raises(ValueError): parse(root)


@pytest.mark.parametrize('name', ['translation', 'codon_start', 'transl_table', 'protein_id', 'gene', 'product'])
def test_missing_required_qualifiers(name):
    root = fixture_xml(); parent = root.find('.//GBFeature_quals')
    for q in list(parent):
        if q.findtext('GBQualifier_name') == name: parent.remove(q)
    with pytest.raises(ValueError, match='Missing'): parse(root)


@pytest.mark.parametrize('name', ['exception', 'transl_except', 'ribosomal_slippage', 'pseudo', 'pseudogene'])
def test_translation_exceptions_rejected(name):
    root = fixture_xml()
    q = ET.SubElement(root.find('.//GBFeature_quals'), 'GBQualifier')
    ET.SubElement(q, 'GBQualifier_name').text = name
    with pytest.raises(ValueError, match='unsupported'): parse(root)


@pytest.mark.parametrize('tag,parent', [('GBFeature_partial5', 'GBFeature'),
    ('GBFeature_partial3', 'GBFeature'), ('GBInterval_iscomp', 'GBInterval')])
def test_partial_and_reverse_flags(tag, parent):
    root = fixture_xml(); ET.SubElement(root.find('.//' + parent), tag, value='true')
    with pytest.raises(ValueError): parse(root)


@pytest.mark.parametrize('mode', ['no_cds', 'two_cds', 'two_intervals', 'error'])
def test_missing_or_complex_records(mode):
    root = fixture_xml()
    if mode == 'no_cds': root.find('.//GBFeature_key').text = 'gene'
    elif mode == 'two_cds':
        root.find('.//GBSeq_feature-table').append(ET.fromstring(ET.tostring(root.find('.//GBFeature'))))
    elif mode == 'two_intervals':
        root.find('.//GBFeature_intervals').append(ET.fromstring(ET.tostring(root.find('.//GBInterval'))))
    else: ET.SubElement(root, 'ERROR').text = 'Unavailable'
    with pytest.raises(ValueError): parse(root)


@pytest.mark.parametrize('xml', [b'', b'<broken', b'<GBSet/>'])
def test_malformed_xml(xml):
    with pytest.raises(ValueError): parse_cds_annotation(xml, ACCESSION, TRANSCRIPT)


@pytest.mark.parametrize('tag,parent', [('GBFeature_operator', 'GBFeature'), ('GBInterval_interbp', 'GBInterval')])
def test_unsupported_structure_flags(tag, parent):
    root = fixture_xml(); ET.SubElement(root.find('.//' + parent), tag).text = 'join'
    with pytest.raises(ValueError): parse(root)


@pytest.mark.parametrize('cds', ['', 'ATGNNNTAA', 'AAATAA', 'ATGAAA', 'ATGAAATA',
                               'ATGTAATAA', 'CCCAAATAA'])
def test_invalid_complete_cds(cds):
    with pytest.raises(ValueError): translate_cds(cds)


@pytest.mark.parametrize('cds,protein', [('ATGAAATAA', 'MK'), ('ATGCCCTAG', 'MP'),
                                      ('ATGGGGTGA', 'MG'), ('CTGAAATAA', 'MK')])
def test_independent_codon_references(cds, protein):
    # CTG is an accepted alternative initiation codon under NCBI table 1.
    assert translate_cds(cds) == protein


@pytest.mark.parametrize('data', [b'', b'>OTHER.1\nMK\n', b'>NP_000001.1\nMX\n',
    b'>NP_000001.1\nMK*\n', b'>NP_000001.1\nMK\n>second\nMK',
    b'>NP_000001.1\nM\vK', b'>NP_000001.1\nM K'])
def test_invalid_protein_fasta(data):
    with pytest.raises(ValueError): parse_protein_fasta(data, 'NP_000001.1')


@pytest.mark.parametrize('start,end', [(0,11), (3,99), (11,3), (True,11)])
def test_coordinate_bounds(start, end):
    annotation = parse(); annotation.update(start=start, end=end)
    with pytest.raises(ValueError): extract_cds(TRANSCRIPT, annotation)


DERIVED = Path(__file__).parents[1] / 'data/derived/cds_protein'
REFERENCE = Path(__file__).parents[1] / 'data/reference'


@pytest.mark.parametrize('species,accession,start,end,utr5,utr3,aa', [
    ('maize', 'NM_001305944.1', 156, 1259, 155, 248, 367),
    ('wheat', 'XM_044542390.1', 96, 734, 95, 351, 212),
])
def test_authentic_offline_evidence(species, accession, start, end, utr5, utr3, aa):
    from agrogenome_eval.sequence_analysis import analyze_fasta
    prefix = species + '_' + accession
    analysis = json.loads((DERIVED / (prefix + '.analysis.json')).read_text())
    transcript = analyze_fasta(REFERENCE / (prefix + '.fasta'))['sequence']
    annotation = parse_cds_annotation((DERIVED / (prefix + '.annotation.xml')).read_bytes(), accession, transcript)
    protein = parse_protein_fasta((DERIVED / (prefix + '.reference_protein.fasta')).read_bytes(), annotation['protein_accession'])
    result = analyze_cds(transcript, annotation, protein)
    assert all(analysis[k] == v for k,v in result.items())
    assert annotation['start'] == start and annotation['end'] == end
    assert result['cds_sequence'] == transcript[start-1:end]
    assert result['cds_length'] == end-start+1
    assert result['five_prime_utr_length'] == utr5 and result['three_prime_utr_length'] == utr3
    assert result['protein_length'] == aa
    summary = json.loads((DERIVED / 'summary.json').read_text())
    record = next(r for r in summary['records'] if r['transcript_accession'] == accession)
    for name, checksum in record['checksums'].items():
        assert hashlib.sha256((DERIVED / name).read_bytes()).hexdigest() == checksum
    snapshot = (DERIVED / 'reference_manifest.snapshot.json').read_bytes()
    assert hashlib.sha256(snapshot).hexdigest() == summary['reference_manifest_sha256']
    assert json.loads(snapshot) == json.loads((REFERENCE / 'manifest.json').read_text())


def runner():
    scripts = str(Path(__file__).parents[1] / 'scripts')
    if scripts not in sys.path: sys.path.insert(0, scripts)
    spec = importlib.util.spec_from_file_location('coding_runner', Path(scripts) / 'analyze_coding_sequences.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def test_retrieval_failure_has_no_success_files(tmp_path):
    module = runner()
    client = Mock(); client.fetch.side_effect = RuntimeError('offline')
    result = module.generate_outputs(tmp_path / 'outputs', REFERENCE, client)
    assert all(r['status'] == 'failed' for r in result['records'])
    assert not list((tmp_path / 'outputs').glob('*.fasta'))


def test_existing_output_is_preserved(tmp_path):
    with pytest.raises(FileExistsError): runner().generate_outputs(tmp_path, REFERENCE, Mock())


def test_offline_runner_with_saved_responses(tmp_path):
    def fetch(accession, rettype, retmode, **kwargs):
        manifest = json.loads((REFERENCE / 'manifest.json').read_text())
        for r in manifest['records']:
            prefix = r['species'] + '_' + r['transcript_accession']
            if accession == r['transcript_accession']:
                assert rettype == 'gb'; return (DERIVED / (prefix + '.annotation.xml')).read_bytes()
            if accession == r['protein_accession']:
                assert kwargs['db'] == 'protein'
                return (DERIVED / (prefix + '.reference_protein.fasta')).read_bytes()
        raise AssertionError('Unexpected accession')
    result = runner().generate_outputs(tmp_path / 'outputs', REFERENCE, Mock(fetch=Mock(side_effect=fetch)))
    assert all(r['status'] == 'validated' for r in result['records'])
    runner().verify_derived_outputs(tmp_path / 'outputs', REFERENCE)
    for record in result['records']:
        saved = json.loads((DERIVED / record['analysis_file']).read_bytes())
        replayed = json.loads((tmp_path / 'outputs' / record['analysis_file']).read_bytes())
        saved.pop('provenance'); replayed.pop('provenance')
        assert saved == replayed
        for name in record['checksums']:
            if not name.endswith('.analysis.json'):
                assert (tmp_path / 'outputs' / name).read_bytes() == (DERIVED / name).read_bytes()


@pytest.mark.parametrize('parent,tag', [
    ('GBSeq','GBSeq_accession-version'), ('GBSeq','GBSeq_sequence'),
    ('GBSeq','GBSeq_moltype'), ('GBFeature','GBFeature_location'),
    ('GBInterval','GBInterval_from'), ('GBInterval','GBInterval_to'),
    ('GBQualifier','GBQualifier_name'), ('GBQualifier','GBQualifier_value'),
])
def test_duplicate_xml_fields_rejected(parent, tag):
    root = fixture_xml()
    ET.SubElement(root.find('.//' + parent), tag).text = 'conflicting-value'
    with pytest.raises(ValueError, match='Duplicate'): parse(root)


def test_duplicate_translation_qualifier_rejected():
    root = fixture_xml(); parent = root.find('.//GBFeature_quals')
    original = next(q for q in parent if q.findtext('GBQualifier_name') == 'translation')
    parent.append(ET.fromstring(ET.tostring(original)))
    with pytest.raises(ValueError, match='duplicate'): parse(root)


@pytest.mark.parametrize('kwargs', [dict(expected_organism='Wrong organism'),
    dict(expected_gene='Wrong gene'), dict(expected_gene_id='999')])
def test_expected_biological_metadata_rejected(kwargs):
    with pytest.raises(ValueError, match='mismatch'):
        parse_cds_annotation(ET.tostring(fixture_xml()), ACCESSION, TRANSCRIPT, **kwargs)


@pytest.mark.parametrize('attribute,text', [('false','true'), ('invalid',''), ('','')])
def test_conflicting_or_malformed_boolean(attribute, text):
    root = fixture_xml()
    flag = ET.SubElement(root.find('.//GBInterval'), 'GBInterval_iscomp', value=attribute)
    flag.text = text
    with pytest.raises(ValueError, match='XML boolean'): parse(root)


@pytest.mark.parametrize('name', ['strand','codon_start','translation_table','start'])
def test_missing_annotation_mapping_fields(name):
    annotation = parse(); del annotation[name]
    with pytest.raises(ValueError, match='Missing'): extract_cds(TRANSCRIPT, annotation)


@pytest.mark.parametrize('species,accession', [('maize','NM_001305944.1'), ('wheat','XM_044542390.1')])
def test_independent_full_codon_translation(species, accession):
    prefix = species + '_' + accession
    xml = ET.fromstring((DERIVED / (prefix + '.annotation.xml')).read_bytes())
    feature = next(f for f in xml.findall('.//GBFeature') if f.findtext('GBFeature_key') == 'CDS')
    qualifiers = {q.findtext('GBQualifier_name'):q.findtext('GBQualifier_value')
                  for q in feature.findall('GBFeature_quals/GBQualifier')}
    start,end = map(int,feature.findtext('GBFeature_location').split('..'))
    transcript = ''.join((REFERENCE / (prefix + '.fasta')).read_text().splitlines()[1:])
    cds = transcript[start-1:end]
    # Independent table-1 lookup, in T,C,A,G codon order, without engine translation.
    bases = 'TCAG'
    amino = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    table = dict(zip([a+b+c for a in bases for b in bases for c in bases], amino))
    translated = ''.join(table[cds[i:i+3]] for i in range(0,len(cds),3))
    protein = ''.join((DERIVED / (prefix + '.reference_protein.fasta')).read_text().splitlines()[1:])
    assert len(cds)%3 == 0 and translated[-1] == '*' and '*' not in translated[:-1]
    assert translated[:-1] == qualifiers['translation'] == protein


@pytest.mark.parametrize('name', ['maize_NM_001305944.1.annotation.xml',
    'maize_NM_001305944.1.cds.fasta', 'wheat_XM_044542390.1.protein.fasta',
    'wheat_XM_044542390.1.analysis.json', 'reference_manifest.snapshot.json'])
def test_derived_corruption_is_detected(tmp_path, name):
    import shutil
    target = tmp_path / 'evidence'; shutil.copytree(DERIVED, target)
    path = target / name; path.write_bytes(path.read_bytes() + b'\n')
    with pytest.raises(ValueError, match='checksum mismatch'):
        runner().verify_derived_outputs(target, REFERENCE)


def saved_response_client():
    records = json.loads((REFERENCE / 'manifest.json').read_text())['records']
    def fetch(accession, rettype, mode, **kwargs):
        record = next(r for r in records if accession in (r['transcript_accession'], r['protein_accession']))
        prefix = record['species'] + '_' + record['transcript_accession']
        suffix = '.annotation.xml' if rettype == 'gb' else '.reference_protein.fasta'
        return (DERIVED / (prefix + suffix)).read_bytes()
    return Mock(fetch=Mock(side_effect=fetch))


@pytest.mark.parametrize('failure', ['write','rename'])
def test_atomic_publication_failures(tmp_path, monkeypatch, failure):
    output = tmp_path / 'outputs'
    if failure == 'rename':
        def fail(*args): raise OSError('simulated publication failure')
        monkeypatch.setattr(Path, 'rename', fail)
    else:
        original = Path.write_bytes
        def write(path, data):
            if path.name.endswith('.protein.fasta'): raise OSError('simulated write failure')
            return original(path, data)
        monkeypatch.setattr(Path, 'write_bytes', write)
    with pytest.raises(OSError, match='simulated'):
        runner().generate_outputs(output, REFERENCE, saved_response_client())
    assert not output.exists()
