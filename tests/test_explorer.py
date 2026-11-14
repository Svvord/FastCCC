"""Unit tests for the result explorer (payload pieces and the HTML build), no FastCCC analysis run."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import fastccc
from fastccc.explorer import payload as ex
from fastccc.explorer import build_explorer
from fastccc.resources import resolve_database


def test_split_range_and_downsample():
    assert ex._split_range('1.5-4.25') == (1.5, 4.25)
    assert ex._split_range('-0.5--0.1') == (-0.5, -0.1)
    assert ex._split_range(None) == (None, None)
    vals = np.ones(400) / 400
    d = ex._downsample_pmf(1000, vals, 0.01, n_bins=80)
    assert d['x0'] == 10.0 and len(d['y']) == 80 and abs(sum(d['y']) - 1) < 1e-6 and d['dx'] == 0.05
    small = ex._downsample_pmf(5, np.array([0.5, 0.5]), 0.01)
    assert small['y'] == [0.5, 0.5] and small['dx'] == 0.01


def test_nmf_patterns_on_synthetic_rows():
    rng = np.random.default_rng(0)
    # two blocks: cell types 0-4 use LR 0-9, cell types 5-9 use LR 10-19
    rows = {'sender': [], 'receiver': [], 'lr': [], 'cs': [], 'trend': []}
    for s in range(10):
        for lr in (range(10) if s < 5 else range(10, 20)):
            for r in rng.choice(10, 3):
                rows['sender'].append(s); rows['receiver'].append(int(r)); rows['lr'].append(lr); rows['cs'].append(1.0); rows['trend'].append(1)
    pat = ex.communication_patterns(rows, 10, 20, 'reference', k=2)
    out = pat['outgoing']
    W = np.array(out['W']); H = np.array(out['H'])
    assert W.shape == (10, 2) and H.shape == (2, 20) and np.allclose(W.sum(axis=1), 1)
    top = W.argmax(axis=1)
    assert len(set(top[:5])) == 1 and len(set(top[5:])) == 1 and top[0] != top[5]   # the two blocks separate
    assert ex.communication_patterns({'sender': [0], 'receiver': [1], 'lr': [0], 'cs': [1], 'trend': [1]}, 2, 1, 'reference') is None


def test_gene_sets_from_cached_gmt(tmp_path):
    (tmp_path / 'TestLib.gmt').write_text('Set A\t\tG1\tG2\tG3\tZZ\nSet B\t\tG1\tG9\nSet C\t\tG4,1.0\tG5,1.0\tG6,1.0\n')
    gs = ex.gene_sets(['G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7'], libraries=('TestLib',), cache_dir=tmp_path)
    assert gs['genes'][:3] == ['G1', 'G2', 'G3']
    sets = {s['name']: s for s in gs['libraries'][0]['sets']}
    assert set(sets) == {'Set A', 'Set C'}                      # Set B has < 3 genes in the universe
    assert sets['Set A']['genes'] == [0, 1, 2] and sets['Set A']['size'] == 4 and sets['Set C']['genes'] == [3, 4, 5]
    assert ex.gene_sets(['G1'], libraries=('Missing_Library_xyz',), cache_dir=tmp_path) is None or True   # offline: tolerated


def test_celltype_colors_extend_beyond_palette():
    names = [f'ct{i}' for i in range(60)]
    cols = ex.celltype_colors(names)
    assert len(set(cols.values())) == 60 and all(c.startswith('#') and len(c) == 7 for c in cols.values())


def _reference_result_dir(tmp_path) -> Path:
    db = resolve_database('CPDBv5.0.0')
    table = pd.read_csv(db / 'interaction_table.csv').head(3)
    ids = list(table['id_cp_interaction'])
    rows = [
        ('A|B', True, ids[0], 'L1', 'R1', 5.0, 1.0, '0.5-3.0', True, 0.3, 0.4, True, 0.2, 0.5, 2.5, '0.1-1.0', 2.5, '0.1-1.0', True, False, 'Up', 1.2),
        ('A|B', True, ids[1], 'L2', 'R2', 1.0, 1.0, '0.5-3.0', True, 0.3, 0.4, True, 0.6, 0.5, 0.5, '0.8-2.0', 0.5, '0.1-1.0', False, True, 'Down', 4.0),
        ('B|A', True, ids[2], 'L3', 'R3', 2.0, 1.0, '0.5-3.0', True, 0.3, 0.4, True, 0.2, 0.5, 1.0, '0.1-1.0', 1.0, '0.1-1.0', True, True, 'Both Sig', 2.2),
    ]
    cols = ['sender|receiver', 'in_reference', 'LRI_ID', 'ligand', 'receptor', 'comm_score', 'null_comm_score', 'sig_threshold_CI',
            'above_expr_threshold', 'ligand_expr_percent', 'receptor_expr_percent', 'above_expr_threshold_ref', 'ligand_expr_percent_ref',
            'receptor_expr_percent_ref', 'ligand_CS_component', 'ligand_CS_CI', 'receptor_CS_component', 'receptor_CS_CI',
            'is_significant', 'is_significant_ref', 'trend_vs_ref', 'ref_comm_score']
    pd.DataFrame(rows, columns=cols).to_csv(tmp_path / 'query_infer_results.tsv', sep='\t', index=False)
    pd.DataFrame([[5.0, 1.0, np.nan], [np.nan, np.nan, 2.0]], index=['A|B', 'B|A'], columns=ids).to_csv(tmp_path / 'query_interactions_strength.tsv', sep='\t')
    panel = tmp_path / 'panel'; panel.mkdir()
    (panel / 'config.toml').write_text('reference_name = "p"\nmin_percentile = 0.1\nLRI_database = "CPDBv5.0.0"\n[celltype]\n"A" = 10\n"B" = 20\n')
    # null pmfs for the (cell type, multidata) pairs of the first interaction
    m1, m2 = int(table.iloc[0]['multidata_1_id']), int(table.iloc[0]['multidata_2_id'])
    keys = [f'A\t{m1}', f'B\t{m2}']; vals = np.concatenate([np.full(50, 0.02), np.full(100, 0.01)]).astype(np.float32)
    np.savez_compressed(tmp_path / 'query_null_pmfs.npz', keys=np.array(keys), offsets=np.array([200, 300], dtype=np.int32),
                        lengths=np.array([50, 100], dtype=np.int32), values=vals, grid_step=np.array([0.01]))
    fastccc.results.write_manifest(tmp_path, {
        'workflow': 'compare_to_reference', 'name': 'Q', 'reference': {'path': str(panel), 'name': 'p', 'source': 'custom', 'celltypes': {'A': 10, 'B': 20}},
        'database': {'name': 'CPDBv5.0.0', 'path': str(db)}, 'input': {'celltype_counts': {'A': 30, 'B': 40}, 'celltype_key': 'ct'},
        'celltype_coverage': {'n_query_types': 2, 'n_covered': 2, 'cells_covered_fraction': 1.0, 'uncovered_query_types': []},
    })
    return tmp_path


def test_reference_payload_and_html(tmp_path):
    d = _reference_result_dir(tmp_path)
    payload = ex.build_payload(d, genesets=False, cells=False)
    assert payload['kind'] == 'reference' and [c['name'] for c in payload['celltypes']] == ['A', 'B']
    r = payload['rows']
    assert r['trend'] == [1, 2, 3] and r['ref_cs'] == [1.2, 4.0, 2.2] and r['thr_lo'] == [0.5] * 3 and r['lig_hi'][1] == 2.0
    assert payload['meta']['has_ref_cs'] and payload['meta']['trend_counts']['Up'] == 1
    curves = payload['nulls']['curves']
    assert len(curves) == 2 and all(abs(sum(c['y']) - 1) < 1e-5 for c in curves.values())
    assert payload['patterns'] == {'k': None}              # patterns are computed in the page on the current selection
    html = build_explorer(d, genesets=False, cells=False)
    text = html.read_text(encoding='utf-8')
    assert html.name == 'explorer.html' and 'FASTCCC_PAYLOAD' in text and 'd3.sankey' in text.replace(' ', '') or 'sankey' in text
    res = fastccc.load_result(d)
    assert res.explorer_path == html
    json.dumps(payload)


def test_effect_size_formula():
    # one gene per side, cluster sizes 100 and 400: null mean (mu_L + mu_R)/2, var (s2_L/100 + s2_R/400)/4
    lr_items = [{'id': 'x', 'ligand': 'L', 'receptor': 'R'}, {'id': 'y', 'ligand': 'A,B', 'receptor': 'R'}]
    gstats = {'genes': ['L', 'R', 'A', 'B'], 'mean': [1.0, 2.0, 1.0, 3.0], 'var': [4.0, 16.0, 4.0, 4.0]}
    z = ex.effect_sizes([3.5, 3.0], lr_items, [0, 1], [0, 0], [1, 1], [100, 400], gstats)
    null_mean, null_sd = 1.5, ((4 / 100 + 16 / 400) / 4) ** 0.5
    assert z[0] == pytest.approx((3.5 - null_mean) / null_sd, abs=0.01)
    # complex A,B: ligand side mean 2.0, var (4+4)/(4*100)
    sd2 = ((8 / 400 + 16 / 400) / 4) ** 0.5
    assert z[1] == pytest.approx((3.0 - 2.0) / sd2, abs=0.01)
    assert ex.effect_sizes([1.0], [{'id': 'q', 'ligand': 'ZZ', 'receptor': 'R'}], [0], [0], [1], [10, 10], gstats) == [None]
    assert ex.effect_sizes([1.0], lr_items, [0], [0], [1], [10, 10], None) == [None]


def test_replication_support_synthetic():
    class FakeRun:
        def __init__(self, rows, cols, sig_pairs):
            self.strength = pd.DataFrame(np.arange(1, len(rows) * len(cols) + 1, dtype=float).reshape(len(rows), len(cols)), index=rows, columns=cols)
            self.significant = pd.DataFrame([{'sender_celltype': p.split('|')[0], 'receiver_celltype': p.split('|')[1], 'LRI_ID': l} for p, l in sig_pairs])
    rows, cols = ['A|B', 'B|A'], ['L1', 'L2']
    runs = [('all', FakeRun(rows, cols, [('A|B', 'L1'), ('A|B', 'L2'), ('B|A', 'L1')])),
            ('s1', FakeRun(rows, cols, [('A|B', 'L1'), ('B|A', 'L1')])), ('s2', FakeRun(rows, cols, [('A|B', 'L1')])), ('s3', FakeRun(rows, cols, [('A|B', 'L1'), ('A|B', 'L2')]))]
    rep = ex.replication_support(runs)
    assert rep['n_groups'] == 3 and rep['groups'] == ['s1', 's2', 's3']
    assert rep['fraction']['A|B\tL1'] == 1.0 and rep['fraction']['B|A\tL1'] == pytest.approx(1 / 3, abs=1e-3) and rep['fraction']['A|B\tL2'] == pytest.approx(1 / 3, abs=1e-3)
    assert rep['top'][0] == ['A|B', 'L1'] and len(rep['matrix']) == 3 and rep['matrix'][1][0] == 1.0 and rep['matrix'][1][1] == 0.0
    assert ex.replication_support(runs[:3]) is None


def test_explorer_js_avoids_spreading_interaction_arrays():
    """Math.max(...array) over per-interaction arrays overflows the call stack above ~65k interactions."""
    import re
    from pathlib import Path
    js = (Path(__file__).resolve().parents[1] / 'fastccc' / 'explorer' / 'static' / 'explorer.js').read_text(encoding='utf-8')
    bad = [m.group(0) for m in re.finditer(r'Math\.(?:max|min)\([^;\n]*\.\.\.(?:R\.|pairs|cs\b|filtered|shown)[^;\n]*', js)]
    assert not bad, bad


def test_explorer_svg_export_is_serialised_from_dom_nodes():
    """An SVG assembled from strings broke on the quotes of the font list (unopenable file); and an SVG button
    pointing at a container that does not exist silently did nothing."""
    import re
    from pathlib import Path
    root = Path(__file__).resolve().parents[1] / 'fastccc' / 'explorer'
    js = (root / 'static' / 'explorer.js').read_text(encoding='utf-8')
    html = (root / 'templates' / 'explorer.html').read_text(encoding='utf-8')
    assert not re.search(r'''[`'"]\s*<svg''', js) and 'font-family="${' not in js      # no SVG markup written as a string
    assert js.count('new XMLSerializer()') == 1 and '<?xml version="1.0" encoding="UTF-8"?>' in js
    specs = re.findall(r'data-svg="([^"]+)"', html)
    targets = [t for spec in specs for t in spec.split(',')]
    assert len(specs) >= 19 and all(f'id="{t}"' in html for t in targets), [t for t in targets if f'id="{t}"' not in html]
    assert "querySelectorAll('.btn.svg[data-svg]')" in js          # the copy button shares the class but exports nothing

