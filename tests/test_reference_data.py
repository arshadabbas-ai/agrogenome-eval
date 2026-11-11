"""Offline synthetic reference fixtures: never presented as authentic NCBI data."""

import hashlib
import importlib.util
import json
from pathlib import Path
from unittest.mock import Mock
from urllib.error import HTTPError, URLError

import pytest

spec = importlib.util.spec_from_file_location(
    "fetch_reference_transcripts", Path(__file__).parents[1] / "scripts/fetch_reference_transcripts.py")
retrieval = importlib.util.module_from_spec(spec)
spec.loader.exec_module(retrieval)


def metadata(candidate, sequence="atgc"):
    return f'''<GBSet><GBSeq>
    <GBSeq_accession-version>{candidate['transcript_accession']}</GBSeq_accession-version>
    <GBSeq_organism>{candidate['scientific_name']}</GBSeq_organism>
    <GBSeq_moltype>mRNA</GBSeq_moltype><GBSeq_length>{len(sequence)}</GBSeq_length>
    <GBSeq_definition>{candidate['gene']} synthetic fixture</GBSeq_definition>
    <GBSeq_comment>{candidate['expected_refseq_status']} REFSEQ: synthetic fixture</GBSeq_comment>
    <GBSeq_sequence>{sequence}</GBSeq_sequence><GBSeq_feature-table>
    <GBFeature><GBFeature_key>gene</GBFeature_key><GBFeature_quals>
    <GBQualifier><GBQualifier_name>gene</GBQualifier_name><GBQualifier_value>{candidate['gene']}</GBQualifier_value></GBQualifier>
    <GBQualifier><GBQualifier_name>db_xref</GBQualifier_name><GBQualifier_value>GeneID:{candidate['gene_id']}</GBQualifier_value></GBQualifier>
    </GBFeature_quals></GBFeature>
    <GBFeature><GBFeature_key>CDS</GBFeature_key><GBFeature_quals>
    <GBQualifier><GBQualifier_name>protein_id</GBQualifier_name><GBQualifier_value>{candidate['protein_accession']}</GBQualifier_value></GBQualifier>
    <GBQualifier><GBQualifier_name>db_xref</GBQualifier_name><GBQualifier_value>GeneID:{candidate['gene_id']}</GBQualifier_value></GBQualifier>
    </GBFeature_quals></GBFeature></GBSeq_feature-table></GBSeq></GBSet>'''.encode()


def fake_client(sequence="atgc"):
    def fetch(accession, rettype, mode):
        candidate = next(c for c in retrieval.CANDIDATES if c['transcript_accession'] == accession)
        if rettype == "gb":
            return metadata(candidate, sequence)
        return f">{accession} {candidate['gene']} synthetic fixture\n{sequence}\n".encode()
    return Mock(fetch=Mock(side_effect=fetch))


def test_offline_dataset_provenance(tmp_path):
    client = fake_client()
    output = tmp_path / "reference"
    result = retrieval.fetch_dataset(output, client)
    assert client.fetch.call_count == 4
    assert json.loads((output / "manifest.json").read_text()) == result
    for record in result['records']:
        assert record['analysis']['length'] == 4
        assert record['analysis']['nucleotide_counts'] == dict(A=1, T=1, G=1, C=1)
        assert record['analysis']['gc_percentage'] == 50
        assert record['sha256'] == hashlib.sha256((output / record['filename']).read_bytes()).hexdigest()
        assert record['validation_status'] == 'verified'
        assert record['retrieval_timestamp'].endswith('+00:00')


@pytest.mark.parametrize('old,new', [
    ('NM_001305944.1', 'NM_001305944.2'), ('Zea mays', 'Triticum aestivum'),
    ('<GBSeq_moltype>mRNA', '<GBSeq_moltype>DNA'), ('GeneID:732788', 'GeneID:999'),
    ('NP_001292873.1', 'NP_000000000.1'), ('VALIDATED REFSEQ', 'MODEL REFSEQ'),
    ('<GBSeq_length>4', '<GBSeq_length>5'),
])
def test_metadata_mismatches(old, new):
    with pytest.raises(ValueError):
        retrieval.verify_metadata(metadata(retrieval.CANDIDATES[0]).replace(old.encode(), new.encode()), retrieval.CANDIDATES[0])


@pytest.mark.parametrize('xml', [b'bad xml', b'<GBSet/>'])
def test_malformed_metadata(xml):
    with pytest.raises(ValueError):
        retrieval.verify_metadata(xml, retrieval.CANDIDATES[0])


def test_missing_status_is_explicit():
    candidate = retrieval.CANDIDATES[0]
    result = retrieval.verify_metadata(metadata(candidate).replace(b'VALIDATED REFSEQ', b'unknown'), candidate)
    assert result['refseq_status'] is None


def test_unsupported_symbols_leave_no_dataset(tmp_path):
    with pytest.raises(ValueError, match='not altered'):
        retrieval.fetch_dataset(tmp_path / 'reference', fake_client('atgn'))
    assert not (tmp_path / 'reference').exists()


def test_existing_output_preserved(tmp_path):
    with pytest.raises(FileExistsError):
        retrieval.fetch_dataset(tmp_path, fake_client())


@pytest.mark.parametrize('email', ['', 'invalid', 'user@host'])
def test_missing_email(email):
    with pytest.raises(ValueError, match='NCBI_EMAIL'):
        retrieval.NCBIClient(email)


def test_retries_are_bounded(monkeypatch):
    opener = Mock(side_effect=URLError('offline'))
    monkeypatch.setattr(retrieval, 'urlopen', opener)
    monkeypatch.setattr(retrieval.time, 'sleep', Mock())
    with pytest.raises(RuntimeError, match='3 attempts'):
        retrieval.NCBIClient('researcher@example.org').fetch('TEST.1', 'fasta', 'text')
    assert opener.call_count == 3


def test_nonretryable_http_error(monkeypatch):
    opener = Mock(side_effect=HTTPError('url', 404, 'missing', {}, None))
    monkeypatch.setattr(retrieval, 'urlopen', opener)
    with pytest.raises(RuntimeError, match='404'):
        retrieval.NCBIClient('researcher@example.org').fetch('TEST.1', 'fasta', 'text')
    assert opener.call_count == 1


def test_request_parameters_and_timeout(monkeypatch):
    from urllib.parse import parse_qs, urlparse
    response = Mock()
    response.__enter__ = Mock(return_value=response)
    response.__exit__ = Mock(return_value=False)
    response.read.return_value = b'>TEST.1\nATGC\n'
    opener = Mock(return_value=response)
    sleeper = Mock()
    monkeypatch.setattr(retrieval, 'urlopen', opener)
    monkeypatch.setattr(retrieval.time, 'sleep', sleeper)
    client = retrieval.NCBIClient('researcher@example.org')
    client.fetch('TEST.1', 'fasta', 'text')
    client.fetch('TEST.1', 'gb', 'xml')
    params = parse_qs(urlparse(opener.call_args.args[0]).query)
    assert params['db'] == ['nuccore']
    assert params['id'] == ['TEST.1']
    assert params['tool'] == ['AgroGenomeEval']
    assert params['email'] == ['researcher@example.org']
    assert opener.call_args.kwargs['timeout'] == 30
    sleeper.assert_called_once()


@pytest.mark.parametrize('accession,length,counts,gc', [
    ('NM_001305944.1', 1507, dict(A=383, T=345, G=412, C=367), 77900 / 1507),
    ('XM_044542390.1', 1085, dict(A=232, T=217, G=321, C=315), 63600 / 1085),
])
def test_downloaded_reference_integrity(accession, length, counts, gc):
    directory = Path(__file__).parents[1] / 'data/reference'
    if not (directory / 'manifest.json').exists():
        pytest.skip('Authentic dataset has not yet been retrieved.')
    manifest = json.loads((directory / 'manifest.json').read_text(encoding='utf-8'))
    record = next(r for r in manifest['records'] if r['transcript_accession'] == accession)
    path = directory / record['filename']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == record['sha256']
    analysis = retrieval.analyze_fasta(path)
    assert analysis == record['analysis']
    assert analysis['record_id'] == accession
    assert analysis['length'] == length
    assert analysis['nucleotide_counts'] == counts
    assert analysis['gc_percentage'] == pytest.approx(gc)


@pytest.mark.parametrize('change', ['accession', 'title', 'sequence', 'empty', 'multiple'])
def test_fasta_response_mismatches(tmp_path, change):
    client = fake_client()
    original = client.fetch.side_effect
    def fetch(accession, rettype, mode):
        data = original(accession, rettype, mode)
        if rettype == 'fasta':
            if change == 'accession': return data.replace(accession.encode(), b'OTHER.1')
            if change == 'title': return data.replace(b'synthetic fixture', b'wrong title')
            if change == 'sequence': return data.replace(b'atgc', b'aaaa')
            if change == 'empty': return b''
            return data + b'>second\nATGC\n'
        return data
    client.fetch.side_effect = fetch
    with pytest.raises(ValueError):
        retrieval.fetch_dataset(tmp_path / 'reference', client)
    assert not (tmp_path / 'reference').exists()


@pytest.mark.parametrize('change', ['checksum', 'analysis', 'file'])
def test_saved_reference_corruption(tmp_path, change):
    output = tmp_path / 'reference'
    record = retrieval.fetch_dataset(output, fake_client())['records'][0]
    path = output / record['filename']
    if change == 'checksum': record['sha256'] = '0' * 64
    elif change == 'analysis': record['analysis']['length'] += 1
    else: path.write_bytes(path.read_bytes() + b'\n')
    with pytest.raises(ValueError, match='mismatch'):
        retrieval.verify_saved_reference(path, record)


def test_incomplete_network_response_retried(monkeypatch):
    from http.client import IncompleteRead
    opener = Mock(side_effect=IncompleteRead(b'partial'))
    monkeypatch.setattr(retrieval, 'urlopen', opener)
    monkeypatch.setattr(retrieval.time, 'sleep', Mock())
    with pytest.raises(RuntimeError, match='3 attempts'):
        retrieval.NCBIClient('researcher@example.org').fetch('TEST.1', 'fasta', 'text')
    assert opener.call_count == 3


@pytest.mark.parametrize('code', [429, 503])
def test_transient_http_errors_retried(monkeypatch, code):
    opener = Mock(side_effect=HTTPError('private-url', code, 'transient', {}, None))
    monkeypatch.setattr(retrieval, 'urlopen', opener)
    monkeypatch.setattr(retrieval.time, 'sleep', Mock())
    with pytest.raises(RuntimeError, match='3 attempts') as exc:
        retrieval.NCBIClient('researcher@example.org').fetch('TEST.1', 'fasta', 'text')
    assert 'private-url' not in str(exc.value)
    assert 'researcher@example.org' not in str(exc.value)
    assert opener.call_count == 3


def test_no_partial_output_on_publish_failure(tmp_path, monkeypatch):
    def fail(*args):
        raise OSError('simulated rename failure')
    monkeypatch.setattr(Path, 'rename', fail)
    with pytest.raises(OSError, match='rename failure'):
        retrieval.fetch_dataset(tmp_path / 'reference', fake_client())
    assert not (tmp_path / 'reference').exists()
