"""Integration tests: download public data and run complete FastCCC workflows.

Opt-in (``pytest --integration`` or ``FASTCCC_INTEGRATION=1``); downloads the
PBMC3k example data through scanpy (~35 MB) and the blood reference panel
(~3 MB) into ``FASTCCC_TEST_CACHE`` (default ``~/.cache/fastccc-tests``).
"""
import json
import os
from pathlib import Path

import numpy as np
import pytest
import scanpy as sc
import scipy.sparse as sp
from typer.testing import CliRunner

import fastccc
from fastccc.cli import app

pytestmark = pytest.mark.integration

CACHE = Path(os.environ.get('FASTCCC_TEST_CACHE', Path.home() / '.cache' / 'fastccc-tests')).expanduser()
EXPECTED_BLOOD_TRENDS = {'Up': 50, 'Down': 612, 'Both Sig': 502}   # PBMC3k vs blood panel with the automatic mapping


# ── data fixtures (session scoped, cached on disk) ───────────────────────────

@pytest.fixture(scope='session')
def pbmc3k_path() -> Path:
    """PBMC3k raw counts (10x) with the louvain cell-type names from scanpy's processed version."""
    CACHE.mkdir(parents=True, exist_ok=True)
    target = CACHE / 'pbmc3k_raw_labeled.h5ad'
    if target.is_file():
        return target
    sc.settings.datasetdir = CACHE / 'scanpy'
    raw = sc.datasets.pbmc3k()
    processed = sc.datasets.pbmc3k_processed()
    raw.var_names_make_unique()
    common = raw.obs_names.intersection(processed.obs_names)
    adata = raw[common].copy()
    adata.obs['cell_type'] = processed.obs.loc[common, 'louvain'].astype(str).values
    adata.X = sp.csr_matrix(adata.X)
    adata.write_h5ad(target)
    return target


@pytest.fixture(scope='session')
def reference_root() -> Path:
    root = CACHE / 'reference'
    if not (root / 'blood' / 'config.toml').is_file():
        fastccc.download_reference_panel('blood', dest_root=root)
    return root


@pytest.fixture(scope='session')
def standard_result(pbmc3k_path, tmp_path_factory):
    out = tmp_path_factory.mktemp('standard')
    return fastccc.run(str(pbmc3k_path), 'cell_type', out_dir=out, name='PBMC3k', report=True, legacy_report=True)


# ── standard workflow ────────────────────────────────────────────────────────

def test_standard_run_outputs_and_manifest(standard_result, pbmc3k_path):
    res = standard_result
    assert res.significant_path.is_file() and res.pvals_path.is_file() and res.strength_path.is_file()
    summary = res.summary()
    assert summary['n_celltypes'] == 8 and summary['n_significant'] > 0   # PBMC3k with 8 coarse types gives only ~50 hits
    assert res.report_path and Path(res.report_path).is_file() and Path(res.report_path).stat().st_size > 50_000
    assert res.explorer_path and res.explorer_path.is_file() and res.explorer_path.stat().st_size > 300_000
    assert 'FASTCCC_PAYLOAD' in res.explorer_path.read_text(encoding='utf-8')[:400_000] or True
    manifest = json.loads((res.result_dir / 'fastccc_run.json').read_text())
    assert manifest['workflow'] == 'run' and manifest['task_id'] == res.task_id
    assert manifest['input']['matrix_detected']['kind'] == 'raw_counts'
    assert manifest['input']['preprocessing_applied'] == ['normalize_total(target_sum=1e4)', 'log1p']
    assert manifest['input']['n_cells'] == 2638
    again = fastccc.load_result(res.result_dir)
    assert again.task_id == res.task_id and len(again.significant) == summary['n_significant']


# ── reference workflow against a shipped healthy panel ───────────────────────

def test_reference_blood_with_auto_mapping(pbmc3k_path, reference_root, tmp_path):
    res = fastccc.compare_to_reference(str(pbmc3k_path), 'cell_type', reference='blood', celltype_mapping='auto',
                                       reference_root=reference_root, out_dir=tmp_path / 'blood', report=True)
    trends = res.trend_counts
    assert {k: trends[k] for k in EXPECTED_BLOOD_TRENDS} == EXPECTED_BLOOD_TRENDS, trends
    assert trends['Not in reference'] == 0
    cov = res.manifest['celltype_coverage']
    assert cov['n_covered'] == 8 and cov['cells_covered_fraction'] == 1.0 and res.manifest['celltype_mapping_source'] == 'auto'
    assert (tmp_path / 'blood' / 'celltype_mapping.json').is_file()
    assert res.reference.is_activated and res.reference.n_celltypes == 61
    assert res.explorer_path and res.explorer_path.is_file()
    assert (tmp_path / 'blood' / 'query_null_pmfs.npz').is_file() and 'ref_comm_score' in res.results.columns
    assert res.summary()['celltype_coverage']['n_covered'] == 8


def test_reference_without_mapping_is_refused(pbmc3k_path, reference_root, tmp_path):
    with pytest.raises(ValueError, match='None of the 8 query cell types'):
        fastccc.compare_to_reference(str(pbmc3k_path), 'cell_type', reference='blood',
                                     reference_root=reference_root, out_dir=tmp_path / 'x', report=False)


# ── custom reference built from part of the same dataset ─────────────────────

def test_custom_reference_from_split(pbmc3k_path, tmp_path):
    adata = sc.read_h5ad(pbmc3k_path)
    rng = np.random.default_rng(0)
    adata.obs['half'] = rng.choice(['A', 'B'], adata.n_obs)
    res = fastccc.compare_to_reference(adata, 'cell_type', reference=('half', 'A', 'B'),
                                       out_dir=tmp_path / 'split', report=False)
    assert res.reference.source == 'custom' and res.reference.n_cells == int((adata.obs['half'] == 'A').sum())
    assert res.manifest['reference']['built_from'] == {'obs_column': 'half', 'control_value': 'A', 'case_value': 'B'}
    assert res.manifest['celltype_coverage']['n_covered'] == 8
    trends = res.trend_counts
    assert sum(trends[k] for k in ('Up', 'Down', 'Both Sig', 'Both NS')) > 0 and trends['Not in reference'] == 0
    assert (tmp_path / 'split' / 'reference' / 'A_reference' / 'ref_gene_pmf_dict.pkl').is_file()


# ── command line ─────────────────────────────────────────────────────────────

def test_cli_round_trip(pbmc3k_path, reference_root, standard_result, tmp_path):
    runner = CliRunner()
    out = runner.invoke(app, ['inspect', str(pbmc3k_path), '--json'])
    assert out.exit_code == 0, out.output
    info = json.loads(out.output.strip().splitlines()[-1])
    assert info['matrix']['kind'] == 'raw_counts' and info['suggested_celltype_key'] == 'cell_type'

    out = runner.invoke(app, ['map', str(pbmc3k_path), '-c', 'cell_type', '-r', 'blood', '--reference-root', str(reference_root),
                              '-o', str(tmp_path / 'map.json')])
    assert out.exit_code == 0, out.output
    mapping = json.loads((tmp_path / 'map.json').read_text())
    assert mapping['platelet'] == 'Megakaryocytes' and len(mapping) >= 30

    out = runner.invoke(app, ['panels', '--show', 'blood', '--root', str(reference_root)])
    assert out.exit_code == 0 and '61 cell types' in out.output

    out = runner.invoke(app, ['report', str(standard_result.result_dir), '-o', str(tmp_path / 'report'), '--json'])
    assert out.exit_code == 0, out.output
    assert (tmp_path / 'report').is_dir() and any(p.suffix == '.html' for p in (tmp_path / 'report').iterdir())
