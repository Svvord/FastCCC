"""Fast smoke tests for the one-call API layer (no FastCCC analysis is run here)."""
import json
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from typer.testing import CliRunner

import fastccc
from fastccc import datacheck, resources, results
from fastccc.cli import app


# ── fixtures ─────────────────────────────────────────────────────────────────

def _make_adata(kind: str, n_cells=300, n_genes=50, seed=0) -> anndata.AnnData:
    rng = np.random.default_rng(seed)
    counts = rng.poisson(0.7, size=(n_cells, n_genes)).astype(np.float32)
    obs = pd.DataFrame({
        'cell_type': rng.choice(['T cell', 'B cell', 'Monocyte'], n_cells),
        'condition': rng.choice(['ctrl', 'case'], n_cells),
        'n_counts': counts.sum(1),
    }, index=[f'cell{i}' for i in range(n_cells)])
    var = pd.DataFrame(index=[f'GENE{i}' for i in range(n_genes)])
    if kind == 'raw':
        return anndata.AnnData(X=sp.csr_matrix(counts), obs=obs, var=var)
    if kind == 'log1p':
        X = np.log1p(counts / np.maximum(counts.sum(1, keepdims=True), 1) * 1e4)
        return anndata.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)
    if kind == 'normalized':
        X = counts / np.maximum(counts.sum(1, keepdims=True), 1) * 1e4
        return anndata.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)
    if kind == 'scaled':
        X = (counts - counts.mean(0)) / (counts.std(0) + 1e-6)
        return anndata.AnnData(X=X, obs=obs, var=var)
    raise ValueError(kind)


# ── matrix detection ─────────────────────────────────────────────────────────

@pytest.mark.parametrize('kind,expected', [('raw', 'raw_counts'), ('log1p', 'log1p'), ('normalized', 'normalized'), ('scaled', 'scaled')])
def test_detect_matrix_type(kind, expected):
    info = datacheck.detect_matrix_type(_make_adata(kind).X)
    assert info.kind == expected


def test_detect_dense_matrix():
    a = _make_adata('raw')
    a.X = a.X.toarray()
    assert datacheck.detect_matrix_type(a.X).kind == 'raw_counts'


def test_prepare_for_standard_normalises_raw_without_touching_input():
    a = _make_adata('raw')
    before = a.X.copy()
    out, info, applied = datacheck.prepare_for_standard(a, copy=True)
    assert info.kind == 'raw_counts' and applied == ['normalize_total(target_sum=1e4)', 'log1p']
    assert datacheck.detect_matrix_type(out.X).kind == 'log1p'
    assert (a.X != before).nnz == 0            # caller's object untouched


def test_prepare_for_standard_keeps_log1p():
    a = _make_adata('log1p')
    out, info, applied = datacheck.prepare_for_standard(a)
    assert out is a and applied == []


def test_prepare_for_standard_rejects_scaled():
    with pytest.raises(ValueError, match='negative'):
        datacheck.prepare_for_standard(_make_adata('scaled'))


def test_prepare_for_reference_finds_raw_in_layers_and_raw():
    raw = _make_adata('raw')
    log = _make_adata('log1p')
    log.layers['counts'] = raw.X.copy()
    out, info, applied = datacheck.prepare_for_reference(log)
    assert info.source == 'layers:counts' and sp.isspmatrix_csr(out.X)
    assert datacheck.detect_matrix_type(out.X).kind == 'raw_counts'

    log2 = _make_adata('log1p')
    log2.raw = raw
    out2, info2, _ = datacheck.prepare_for_reference(log2)
    assert info2.source == 'raw' and out2.n_obs == log2.n_obs


def test_prepare_for_reference_accepts_log1p_without_raw():
    a = _make_adata('log1p')
    out, info, applied = datacheck.prepare_for_reference(a)
    assert info.kind == 'log1p' and out is not a and sp.isspmatrix_csr(out.X)
    assert any('rank-based' in x for x in applied)


def test_prepare_for_reference_rejects_scaled():
    with pytest.raises(ValueError, match='negative'):
        datacheck.prepare_for_reference(_make_adata('scaled'))


# ── fractional counts and detection evidence ─────────────────────────────────

def _placenta_like_counts(n_cells=400, n_genes=300, seed=1):
    rng = np.random.default_rng(seed)
    libsize = rng.lognormal(7, 0.8, n_cells)              # realistic library-size spread
    gene_w = rng.pareto(1.2, n_genes) + 0.05
    lam = np.outer(libsize, gene_w / gene_w.sum())
    return rng.poisson(lam).astype(np.float64)


def test_detect_fractional_counts_as_raw():
    rng = np.random.default_rng(0)
    counts = _placenta_like_counts()
    soupx = sp.csr_matrix(counts); soupx.data = np.clip(soupx.data - rng.uniform(0, 0.9, soupx.data.size), 0, None); soupx.eliminate_zeros()
    info = datacheck.detect_matrix_type(soupx)
    assert info.kind == 'raw_counts' and not info.is_integer and info.confidence == 'medium'
    alevin = sp.csr_matrix(counts); idx = rng.choice(alevin.data.size, alevin.data.size // 10, replace=False)
    alevin.data[idx] += 0.5
    assert datacheck.detect_matrix_type(alevin).kind == 'raw_counts'


def test_detect_normalised_variants_by_row_sums():
    import scanpy as sc
    counts = _placenta_like_counts()
    a = anndata.AnnData(X=sp.csr_matrix(counts)); sc.pp.normalize_total(a, target_sum=1e4)
    info = datacheck.detect_matrix_type(a.X)
    assert info.kind == 'normalized' and info.confidence == 'high' and info.rowsum_cv < 0.01
    b = anndata.AnnData(X=sp.csr_matrix(counts)); sc.pp.normalize_total(b, target_sum=1e4); sc.pp.log1p(b)
    info = datacheck.detect_matrix_type(b.X)
    assert info.kind == 'log1p' and info.confidence == 'high' and info.expm1_rowsum_cv < 0.01
    c = anndata.AnnData(X=sp.csr_matrix(counts)); sc.pp.normalize_total(c, target_sum=1e4, exclude_highly_expressed=True)
    assert datacheck.detect_matrix_type(c.X).kind == 'normalized'


def test_prepare_for_standard_input_type_override():
    a = _make_adata('raw')
    out, info, applied = datacheck.prepare_for_standard(a, input_type='log1p')
    assert out is a and applied == [] and info.user_override and info.kind == 'log1p'
    out, info, applied = datacheck.prepare_for_standard(_make_adata('log1p'), input_type='raw')
    assert applied == ['normalize_total(target_sum=1e4)', 'log1p']
    with pytest.raises(ValueError, match='input_type'):
        datacheck.prepare_for_standard(a, input_type='bogus')


def test_celltype_column_scoring_prefers_annotations_over_metadata():
    n = 300
    obs = pd.DataFrame({
        'orig.ident': np.repeat(['s1', 's2', 's3'], 100),
        'nCount_RNA': np.arange(n),
        'sex': np.tile(['F', 'M'], 150),
        'condition': np.repeat(['ctrl', 'case'], 150),
        'cell_ontology_class': np.tile(['T cell', 'B cell', 'NK cell', 'macrophage', 'fibroblast'], 60),
        'seurat_clusters': pd.Categorical(np.tile([str(i) for i in range(10)], 30)),
    })
    key, conf, cands = datacheck.suggest_celltype_key(obs)
    assert key == 'cell_ontology_class' and conf == 'high'
    order = [c['column'] for c in cands]
    assert order.index('seurat_clusters') < order.index('orig.ident')
    assert all(c['score'] < 0 for c in cands if c['column'] in ('orig.ident', 'sex', 'condition'))


def test_inspect_data_suggests_celltype_column():
    info = datacheck.inspect_data(_make_adata('raw'))
    assert info['suggested_celltype_key'] == 'cell_type' and info['suggested_celltype_key_confidence'] == 'high'
    assert info['matrix']['kind'] == 'raw_counts' and info['raw_counts_available']
    assert info['recommended_workflow_input'] == {'standard': 'ok', 'reference': 'ok'}
    assert datacheck.inspect_data(_make_adata('log1p'))['recommended_workflow_input']['reference'].startswith('ok')
    json.dumps(info)   # must be serialisable


# ── resources ────────────────────────────────────────────────────────────────

def test_packaged_databases_resolve():
    names = resources.list_databases()
    assert 'CPDBv5.0.0' in names and 'CPDBv4.1.0' in names
    p = resources.resolve_database('CPDBv5.0.0')
    assert (p / 'interaction_table.csv').is_file()
    # legacy relative path falls back to the packaged copy
    assert resources.resolve_database('./db/CPDBv5.0.0') == p
    assert resources.resolve_database('db/CPDBv5.0.0/') == p
    with pytest.raises(FileNotFoundError):
        resources.resolve_database('NoSuchDB')


def test_reference_home_env(tmp_path, monkeypatch):
    monkeypatch.setenv('FASTCCC_HOME', str(tmp_path))
    assert resources.fastccc_home() == tmp_path
    assert resources.reference_roots()[0] == tmp_path / 'reference'
    assert 'liver' in resources.list_reference_panels(include_downloadable=True)


def test_resolve_reference_custom_dir(tmp_path):
    panel = tmp_path / 'mypanel'
    panel.mkdir()
    (panel / 'config.toml').write_text('reference_name = "mypanel"\nmin_percentile = 0.1\nLRI_database = "CPDBv5.0.0"\n[celltype]\n"T cell" = 10\n')
    info = resources.resolve_reference(str(panel))
    assert info['source'] == 'custom'
    p = results.ReferencePanel.load(panel)
    assert p.name == 'mypanel' and p.n_cells == 10 and p.database == 'CPDBv5.0.0' and not p.is_activated


# ── manifest & result objects ────────────────────────────────────────────────

def test_manifest_roundtrip(tmp_path):
    m = results.write_manifest(tmp_path, {'workflow': 'run', 'task_id': 'abc123', 'path': Path('/x'), 'n': np.int64(3)})
    assert m['fastccc_version'] == fastccc.__version__ and m['path'] == '/x' and m['n'] == 3
    assert results.read_manifest(tmp_path)['task_id'] == 'abc123'
    results.update_manifest(tmp_path, outputs={'report_html': 'r.html'})
    assert results.read_manifest(tmp_path)['outputs']['report_html'] == 'r.html'


def test_runresult_from_minimal_dir(tmp_path):
    tid = 'abc123'
    pd.DataFrame({'A': [0.01, 1.0]}, index=['T cell|B cell', 'B cell|T cell']).to_csv(tmp_path / f'{tid}_Cauchy_pvals.tsv', sep='\t')
    pd.DataFrame({'A': [0.5, 0.2]}, index=['T cell|B cell', 'B cell|T cell']).to_csv(tmp_path / f'{tid}_average_interactions_strength.tsv', sep='\t')
    pd.DataFrame({'sender_celltype': ['T cell'], 'receiver_celltype': ['B cell'], 'LRI_ID': ['A'],
                  'ligand': ['L'], 'receptor': ['R'], 'p-value': [0.01]}).to_csv(tmp_path / f'{tid}_significant_results.tsv', sep='\t', index=False)
    r = results.RunResult.load(tmp_path)          # works without a manifest (legacy directories)
    assert r.task_id == tid and r.celltypes == ['B cell', 'T cell']
    s = r.summary()
    assert s['n_significant'] == 1 and s['top_senders'][0]['celltype'] == 'T cell'
    json.dumps(r.to_dict())
    assert isinstance(results.load_result(tmp_path), results.RunResult)


# ── CLI ──────────────────────────────────────────────────────────────────────

def test_cli_help_and_version():
    runner = CliRunner()
    assert runner.invoke(app, ['--help']).exit_code == 0
    for cmd in ['inspect', 'run', 'ref', 'build-ref', 'panels', 'report']:
        assert runner.invoke(app, [cmd, '--help']).exit_code == 0, cmd
    out = runner.invoke(app, ['--version'])
    assert out.exit_code == 0 and fastccc.__version__ in out.output


def test_cli_inspect(tmp_path):
    a = _make_adata('raw')
    a.write_h5ad(tmp_path / 'x.h5ad')
    out = CliRunner().invoke(app, ['inspect', str(tmp_path / 'x.h5ad'), '--json'])
    assert out.exit_code == 0, out.output
    assert json.loads(out.output.strip().splitlines()[-1])['suggested_celltype_key'] == 'cell_type'


def test_logging_keeps_user_handlers():
    from loguru import logger
    # importing fastccc must not have removed every handler
    assert logger._core.handlers


def test_reference_result_interaction_summary_and_repr(tmp_path):
    rows = [
        # sender|receiver, in_reference, LRI_ID, ligand, receptor, comm_score, ligand_expr_percent, receptor_expr_percent,
        # ligand_expr_percent_ref, receptor_expr_percent_ref, is_significant, is_significant_ref, trend_vs_ref
        ('A|B', True, 'L1_R1', 'L1', 'R1', 5.0, 0.1, 0.5, 0.4, 0.6, False, True, 'Down'),
        ('A|C', True, 'L1_R1', 'L1', 'R1', 2.0, 0.1, 0.4, 0.4, 0.5, False, True, 'Down'),
        ('B|C', True, 'L2_R2', 'L2', 'R2', 9.0, 0.3, 0.3, 0.3, 0.3, False, True, 'Down'),
        ('A|B', True, 'L3_R3', 'L3', 'R3', 7.0, 0.5, 0.5, 0.1, 0.1, True, False, 'Up'),
        ('A|B', True, 'L4_R4', 'L4', 'R4', 1.0, 0.5, 0.5, 0.5, 0.5, False, False, 'Both NS'),
        ('D|D', False, 'L4_R4', 'L4', 'R4', 1.0, 0.5, 0.5, None, None, False, None, None),
    ]
    cols = ['sender|receiver', 'in_reference', 'LRI_ID', 'ligand', 'receptor', 'comm_score', 'ligand_expr_percent',
            'receptor_expr_percent', 'ligand_expr_percent_ref', 'receptor_expr_percent_ref', 'is_significant',
            'is_significant_ref', 'trend_vs_ref']
    pd.DataFrame(rows, columns=cols).to_csv(tmp_path / 'query_infer_results.tsv', sep='\t', index=False)
    panel = tmp_path / 'panel'; panel.mkdir()
    (panel / 'config.toml').write_text('reference_name = "p"\nmin_percentile = 0.1\nLRI_database = "CPDBv5.0.0"\n[celltype]\n"A" = 1\n')
    results.write_manifest(tmp_path, {'workflow': 'compare_to_reference', 'name': 'Q', 'reference': {'path': str(panel), 'name': 'p'},
                                      'celltype_mapping': {'x': 'A'}, 'celltype_coverage': {'n_covered': 3}})
    res = fastccc.load_result(tmp_path)
    assert res.trend_counts == {'Up': 1, 'Both Sig': 0, 'Down': 3, 'Both NS': 1, 'Not in reference': 1}
    assert 'BothNS=1' in repr(res) and 'NotInReference=1' in repr(res)
    assert res.celltype_mapping == {'x': 'A'} and res.celltype_coverage == {'n_covered': 3}
    down = res.interaction_summary('Down')
    assert list(down['ligand']) == ['L1', 'L2'] and down.iloc[0]['n_celltype_pairs'] == 2
    assert down.iloc[0]['top_sender_receiver'] == 'A|B' and down.iloc[0]['max_comm_score'] == 5.0
    assert down.iloc[0]['ligand_pct_ref'] == pytest.approx(0.4) and down.iloc[0]['ligand_pct_query'] == pytest.approx(0.1)
    assert len(res.interaction_summary('Up')) == 1 and res.interaction_summary('Both Sig').empty
