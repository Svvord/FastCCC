"""Tests for cell-type mapping suggestion, validation and coverage (no FastCCC analysis is run)."""
import json

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from typer.testing import CliRunner

import fastccc
from fastccc import celltype_match as cm
from fastccc.cli import app

BLOOD = ['B-2 B cell', 'CD141-positive myeloid dendritic cell', 'CD16-negative, CD56-bright natural killer cell, human',
         'CD16-positive, CD56-dim natural killer cell, human', 'CD38-negative naive B cell', 'CD4-positive helper T cell',
         'CD4-positive, alpha-beta T cell', 'CD4-positive, alpha-beta cytotoxic T cell', 'CD4-positive, alpha-beta memory T cell',
         'CD8-positive, alpha-beta T cell', 'CD8-positive, alpha-beta cytokine secreting effector T cell',
         'CD8-positive, alpha-beta memory T cell', 'CD8-positive, alpha-beta memory T cell, CD45RO-positive', 'IgA plasmablast',
         'IgG plasmablast', 'T cell', 'T follicular helper cell', 'activated CD4-positive, alpha-beta T cell',
         'activated CD4-positive, alpha-beta T cell, human', 'activated CD8-positive, alpha-beta T cell, human', 'alpha-beta T cell',
         'alveolar macrophage', 'basophil', 'central memory CD4-positive, alpha-beta T cell', 'classical monocyte',
         'common myeloid progenitor', 'conventional dendritic cell', 'dendritic cell', 'double negative thymocyte',
         'effector CD8-positive, alpha-beta T cell', 'effector memory CD4-positive, alpha-beta T cell',
         'effector memory CD8-positive, alpha-beta T cell', 'effector memory CD8-positive, alpha-beta T cell, terminally differentiated',
         'erythrocyte', 'gamma-delta T cell', 'granulocyte', 'group 3 innate lymphoid cell', 'hematopoietic stem cell',
         'immature neutrophil', 'intermediate monocyte', 'lymphocyte', 'macrophage', 'mature NK T cell', 'memory B cell', 'monocyte',
         'mucosal invariant T cell', 'naive B cell', 'naive thymus-derived CD4-positive, alpha-beta T cell',
         'naive thymus-derived CD8-positive, alpha-beta T cell', 'natural killer cell', 'neutrophil', 'non-classical monocyte',
         'plasma cell', 'plasmablast', 'plasmacytoid dendritic cell', 'plasmacytoid dendritic cell, human', 'platelet',
         'regulatory T cell', 'transitional stage B cell', 'type I NK T cell', 'unswitched memory B cell']
BLOOD_COUNTS = {name: 100 + i for i, name in enumerate(BLOOD)}
PBMC = {'CD4 T cells': 1144, 'CD14+ Monocytes': 480, 'B cells': 342, 'CD8 T cells': 316, 'NK cells': 154,
        'FCGR3A+ Monocytes': 150, 'Dendritic cells': 37, 'Megakaryocytes': 15}


# ── normalisation and similarity ─────────────────────────────────────────────

@pytest.mark.parametrize('label,expected', [
    ('CD4 T cells', {'cd4', 't'}),
    ('CD14+ Monocytes', {'cd14', 'classical', 'monocyte'}),
    ('FCGR3A+ Monocytes', {'fcgr3a', 'nonclassical', 'monocyte'}),
    ('NK cells', {'nk'}),
    ('natural killer cell', {'nk'}),
    ('Megakaryocytes', {'platelet'}),
    ('Tregs', {'regulatory', 't', 'cd4'}),
    ('pDC', {'plasmacytoid', 'dendritic'}),
    ('naive thymus-derived CD4-positive, alpha-beta T cell', {'naive', 'cd4', 't'}),
    ('gamma-delta T cell', {'gammadelta', 't'}),
    ('Fibroblasts', {'fibroblast'}),
])
def test_normalise_label(label, expected):
    assert cm.normalise_label(label) == expected


def test_similarity_synonyms_and_conflicts():
    assert cm.label_similarity('platelet', 'Megakaryocytes') == 1.0
    assert cm.label_similarity('natural killer cell', 'NK cells') == 1.0
    assert cm.label_similarity('plasmacytoid dendritic cell', 'pDC') == 1.0
    assert cm.label_similarity('regulatory T cell', 'Tregs') == 1.0
    assert cm.label_similarity('classical monocyte', 'CD14+ Monocytes') >= 0.75
    assert cm.label_similarity('CD4-positive, alpha-beta T cell', 'CD8 T cells') == 0.0      # CD4 vs CD8
    assert cm.label_similarity('gamma-delta T cell', 'CD8 T cells') == 0.0                   # gd vs alpha-beta
    assert cm.label_similarity('B cell', 'T cells') == 0.0                                   # different lineage
    assert cm.label_similarity('classical monocyte', 'FCGR3A+ Monocytes') == 0.0            # classical vs non-classical
    assert cm.label_similarity('regulatory T cell', 'CD4 T cells') > cm.label_similarity('T cell', 'CD4 T cells') > 0


# ── suggestion ───────────────────────────────────────────────────────────────

def test_suggest_mapping_pbmc_vs_blood():
    sug = cm.suggest_mapping_from_labels(BLOOD_COUNTS, PBMC, reference_name='blood')
    assert sug.coverage['n_covered'] == 8 and sug.coverage['cells_covered_fraction'] == 1.0
    assert set(sug.mapping.values()) <= set(PBMC) and set(sug.mapping) <= set(BLOOD)
    assert sug.mapping['platelet'] == 'Megakaryocytes'
    assert sug.mapping['classical monocyte'] == 'CD14+ Monocytes'
    assert sug.mapping['non-classical monocyte'] == 'FCGR3A+ Monocytes'
    assert sug.mapping['regulatory T cell'] == 'CD4 T cells'
    import re
    for ref, q in sug.mapping.items():                              # never cross the CD4/CD8 line
        if re.search(r'\bCD8\b', ref):
            assert q == 'CD8 T cells', ref
        if re.search(r'\bCD4\b', ref):
            assert q == 'CD4 T cells', ref
    ambiguous = {a['reference'] for a in sug.ambiguous}
    assert {'T cell', 'alpha-beta T cell', 'monocyte'} <= ambiguous   # generic types with tied fine-grained queries
    assert 'erythrocyte' in sug.unmatched_reference and not sug.unmatched_query
    tbl = sug.table()
    assert list(tbl['status'].unique()) == ['mapped'] and tbl.loc[tbl.query_celltype == 'CD4 T cells', 'n_reference_types'].item() >= 8
    json.dumps(sug.to_dict())


def test_generic_query_absorbs_ties_and_unmatched_query_reported():
    ref = {'CD4-positive, alpha-beta T cell': 10, 'CD8-positive, alpha-beta T cell': 10, 'T cell': 5, 'hepatocyte': 7}
    sug = cm.suggest_mapping_from_labels(ref, {'T cells': 50, 'Fibroblasts': 20})
    assert sug.mapping == {'CD4-positive, alpha-beta T cell': 'T cells', 'CD8-positive, alpha-beta T cell': 'T cells', 'T cell': 'T cells'}
    assert sug.unmatched_query == ['Fibroblasts'] and sug.unmatched_reference == ['hepatocyte']
    assert sug.coverage['uncovered_cells_fraction'] == pytest.approx(20 / 70, abs=1e-3)


def test_coverage_follows_core_rules():
    ref = {'T cell': 10, 'B cell': 10, 'NK cell': 10}
    qry = {'T cell': 100, 'B cell': 50, 'Myeloid': 50}
    cov = cm.compute_coverage(ref, qry, mapping=None)
    assert cov['n_covered'] == 2 and cov['uncovered_query_types'] == ['Myeloid'] and cov['cells_covered_fraction'] == 0.75
    # a reference type used as a mapping key no longer matches by name
    cov2 = cm.compute_coverage(ref, qry, mapping={'T cell': 'Myeloid'})
    assert set(cov2['uncovered_query_types']) == {'T cell'} and cov2['n_covered'] == 2
    assert cm.compute_coverage(ref, {'X': 1, 'Y': 2})['n_covered'] == 0


def test_validate_mapping_hints(tmp_path):
    with pytest.raises(ValueError) as err:
        cm.validate_celltype_mapping({'classical monocytes': 'CD14+ Monocytes', 'platelet': 'Megakaryocyte'}, BLOOD_COUNTS, PBMC)
    msg = str(err.value)
    assert "did you mean 'classical monocyte'" in msg and "did you mean 'Megakaryocytes'" in msg
    ok = cm.validate_celltype_mapping({'platelet': 'Megakaryocytes'}, BLOOD_COUNTS, PBMC)
    assert ok == {'platelet': 'Megakaryocytes'}
    sug = cm.suggest_mapping_from_labels(BLOOD_COUNTS, PBMC)
    full = sug.to_json(tmp_path / 'full.json', full=True)
    plain = sug.to_json(tmp_path / 'map.json')
    assert cm.load_mapping(full) == cm.load_mapping(plain) == sug.mapping


# ── API and CLI wiring ───────────────────────────────────────────────────────

def _panel_dir(tmp_path, counts, name='mypanel'):
    panel = tmp_path / name
    panel.mkdir()
    body = '\n'.join(f'"{k}" = {v}' for k, v in counts.items())
    (panel / 'config.toml').write_text(f'reference_name = "{name}"\nmin_percentile = 0.1\nLRI_database = "CPDBv5.0.0"\n[celltype]\n{body}\n')
    return panel


def _query_adata(labels, n_genes=30):
    rng = np.random.default_rng(0)
    n = len(labels)
    X = sp.csr_matrix(rng.poisson(1.0, size=(n, n_genes)).astype(np.float32))
    obs = pd.DataFrame({'cell_type': labels, 'donor': rng.choice(['d1', 'd2'], n)}, index=[f'c{i}' for i in range(n)])
    return anndata.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f'G{i}' for i in range(n_genes)]))


def test_api_suggest_and_panel_helpers(tmp_path):
    panel = _panel_dir(tmp_path, BLOOD_COUNTS)
    labels = [k for k, v in PBMC.items() for _ in range(max(1, v // 50))]
    adata = _query_adata(labels)
    sug = fastccc.suggest_celltype_mapping(adata, 'cell_type', reference=str(panel))
    assert sug.reference_name == 'mypanel' and sug.coverage['n_covered'] == 8
    sug2 = fastccc.suggest_celltype_mapping(PBMC, reference=fastccc.ReferencePanel.load(panel))
    assert sug2.mapping == sug.mapping
    p = fastccc.ReferencePanel.load(panel)
    assert p.celltype_counts == BLOOD_COUNTS and p.celltype_table().iloc[0]['celltype'] == BLOOD[-1]
    assert 'cell types' in p.describe()


def test_compare_to_reference_refuses_zero_coverage(tmp_path):
    panel = _panel_dir(tmp_path, {'hepatocyte': 10, 'cholangiocyte': 5})
    adata = _query_adata(['CD4 T cells'] * 20 + ['B cells'] * 10)
    with pytest.raises(ValueError, match='None of the 2 query cell types'):
        fastccc.compare_to_reference(adata, 'cell_type', reference=str(panel), out_dir=tmp_path / 'out', report=False)
    assert not (tmp_path / 'out').exists()
    with pytest.raises(ValueError, match='did you mean'):
        fastccc.compare_to_reference(adata, 'cell_type', reference=str(panel), celltype_mapping={'hepatocytes': 'CD4 T cells'},
                                     out_dir=tmp_path / 'out2', report=False)


def test_build_reference_subset_parsing():
    from fastccc.api import _parse_subset
    assert _parse_subset('condition:control') == ('condition', 'control')
    assert _parse_subset(('condition', 'control')) == ('condition', 'control')
    with pytest.raises(ValueError):
        _parse_subset('nocolon')


def test_cli_map_and_panels_show(tmp_path):
    panel = _panel_dir(tmp_path, BLOOD_COUNTS)
    labels = [k for k, v in PBMC.items() for _ in range(max(1, v // 50))]
    h5 = tmp_path / 'q.h5ad'
    _query_adata(labels).write_h5ad(h5)
    runner = CliRunner()
    out = runner.invoke(app, ['map', str(h5), '-c', 'cell_type', '-r', str(panel), '-o', str(tmp_path / 'map.json')])
    assert out.exit_code == 0, out.output
    assert 'Megakaryocytes' in out.output and 'ambiguous' in out.output and (tmp_path / 'map.json').is_file()
    out = runner.invoke(app, ['map', str(h5), '-c', 'cell_type', '-r', str(panel), '--json'])
    assert out.exit_code == 0 and json.loads(out.output.strip().splitlines()[-1])['coverage']['n_covered'] == 8
    out = runner.invoke(app, ['panels', '--show', str(panel)])
    assert out.exit_code == 0 and 'platelet' in out.output and '61 cell types' in out.output
    assert runner.invoke(app, ['map', '--help']).exit_code == 0
    assert '--auto-map' in runner.invoke(app, ['ref', '--help']).output
    assert '--subset' in runner.invoke(app, ['build-ref', '--help']).output
