"""Build the data payload of the interactive result explorer.

The explorer is a self-contained HTML page: this module turns a FastCCC result
folder (standard run, split run or reference comparison) into one JSON object
that the page's JavaScript renders and filters.  Everything is derived from the
files FastCCC already writes, plus per-gene expression statistics computed from
the input dataset when it is still available.
"""
from __future__ import annotations

import glob
import math
import os
import re
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from loguru import logger

from ..results import RunResult, MultiRunResult, ReferenceResult, load_result
from ..resources import resolve_database
from ..datacheck import load_anndata, prepare_for_standard

TREND_CODES = {'Up': 1, 'Down': 2, 'Both Sig': 3, 'Both NS': 4}
TREND_NAMES = {0: 'Not in reference', 1: 'Up', 2: 'Down', 3: 'Both Sig', 4: 'Both NS'}

# soft categorical palette in the spirit of the paper figures (cell types are
# always labelled, so hue only needs to be pleasant and distinguishable)
CELLTYPE_PALETTE = [
    '#e07b91', '#6fa8dc', '#f2b134', '#9b8bc4', '#ef8f6a', '#5aa9a0', '#c45b8a', '#a8c36b',
    '#d8a0c0', '#7b9ac9', '#d9b27c', '#b09ec6', '#e9a44b', '#8dbfce', '#9d5c7e', '#c7d8a0',
    '#f29393', '#5c7fa6', '#f5bf85', '#ac6e9b', '#2fa682', '#8298aa', '#d53456', '#fbb040',
    '#4c688e', '#e8bcd6', '#7dc582', '#6d859e', '#f2a450', '#a4d7ee', '#747a9e', '#b11b51',
    '#dbe9be', '#2e5471', '#e9b5d2', '#7493c1', '#d4ae80', '#78acd9', '#9c86b8', '#ed8561',
]


def _round(x, nd=3):
    if x is None:
        return None
    try:
        f = float(x)
    except (TypeError, ValueError):
        return None
    if math.isnan(f) or math.isinf(f):
        return None
    return round(f, nd)


def _pval(x):
    try:
        f = float(x)
    except (TypeError, ValueError):
        return None
    if math.isnan(f):
        return None
    return float(f'{f:.3g}')


def _split_range(s) -> tuple[Optional[float], Optional[float]]:
    """'1.23-4.56' -> (1.23, 4.56); handles negatives via rsplit on the separator."""
    if s is None or (isinstance(s, float) and math.isnan(s)):
        return None, None
    txt = str(s)
    m = re.match(r'^\s*(-?[0-9.eE+-]+?)-(-?[0-9.eE+-]+)\s*$', txt)
    if not m:
        return None, None
    return _round(m.group(1)), _round(m.group(2))


def db_metadata(db_path) -> Optional[dict]:
    """Species, version and provenance of the database (its ``database.json``), for the explorer's About panel."""
    if db_path is None:
        return None
    from ..dbtools import database_info
    m = database_info(db_path)
    keys = ('name', 'species', 'version', 'source', 'citation', 'license', 'url',
            'n_interactions', 'n_genes', 'n_complexes', 'has_classification')
    return {k: m.get(k) for k in keys}


def db_annotations(db_path) -> dict[str, dict]:
    """``{CPI id: {'classification', 'directionality'}}`` from the LRI database."""
    table = pd.read_csv(Path(db_path) / 'interaction_table.csv')
    out = {}
    for _, r in table.iterrows():
        out[str(r['id_cp_interaction'])] = {
            'classification': None if pd.isna(r.get('classification')) else str(r['classification']),
            'directionality': None if pd.isna(r.get('directionality')) else str(r['directionality']),
        }
    return out


def celltype_colors(names: list[str]) -> dict[str, str]:
    """One colour per cell type; beyond the curated palette, evenly spaced soft hues (no repeats)."""
    import colorsys
    out = {}
    for i, n in enumerate(names):
        if i < len(CELLTYPE_PALETTE):
            out[n] = CELLTYPE_PALETTE[i]
        else:
            k = i - len(CELLTYPE_PALETTE)
            h = (k * 0.618033988749895) % 1.0          # golden-angle steps spread hues evenly
            r, g, b = colorsys.hls_to_rgb(h, 0.62 if k % 2 == 0 else 0.5, 0.45)
            out[n] = '#%02x%02x%02x' % (int(r * 255), int(g * 255), int(b * 255))
    return out


def gene_expression_stats(data_path: Optional[str], celltype_key: Optional[str], genes: list[str],
                          max_cells: int = 300_000) -> Optional[dict]:
    """Per cell type: fraction of expressing cells and mean log-normalised expression of ``genes``."""
    if not data_path or not celltype_key or not Path(data_path).is_file():
        return None
    try:
        adata = load_anndata(data_path)
        if celltype_key not in adata.obs.columns:
            return None
        if adata.n_obs > max_cells:
            rng = np.random.default_rng(0)
            adata = adata[np.sort(rng.choice(adata.n_obs, max_cells, replace=False))].copy()
        adata, _, _ = prepare_for_standard(adata, copy=True)
        present = [g for g in dict.fromkeys(genes) if g in adata.var_names]
        if not present:
            return None
        sub = adata[:, present]
        X = sub.X
        labels = sub.obs[celltype_key].astype(str).values
        cts = list(dict.fromkeys(labels))
        pct, mean = [], []
        for ct in cts:
            rows = np.where(labels == ct)[0]
            M = X[rows]
            if hasattr(M, 'toarray'):
                nnz = np.asarray((M > 0).sum(axis=0)).ravel()
                mu = np.asarray(M.mean(axis=0)).ravel()
            else:
                M = np.asarray(M)
                nnz = (M > 0).sum(axis=0)
                mu = M.mean(axis=0)
            pct.append([round(float(v), 3) for v in nnz / max(len(rows), 1)])
            mean.append([round(float(v), 3) for v in mu])
        return {'genes': present, 'celltypes': cts, 'pct': pct, 'mean': mean, 'n_cells_used': int(sub.n_obs)}
    except Exception as err:   # expression stats are optional
        logger.warning(f"Gene expression statistics skipped: {err}")
        return None


def _genes_of(lr_table: list[dict]) -> list[str]:
    genes = []
    for lr in lr_table:
        for part in (lr['ligand'], lr['receptor']):
            genes.extend(g.strip() for g in str(part).split(',') if g.strip())
    return list(dict.fromkeys(genes))


EFFECT_MODULE = 'Mean_Average_Arithmetic'   # the module whose null has closed-form moments


def gene_global_stats(data_path: Optional[str], genes: list[str], split_by: Optional[str] = None,
                      max_cells: int = 300_000) -> Optional[dict]:
    """Mean and variance of the log-normalised expression of ``genes`` over all cells (and per
    ``split_by`` group): the first two moments of FastCCC's permutation null for the mean summary."""
    if not data_path or not Path(data_path).is_file():
        return None
    try:
        adata = load_anndata(data_path)
        if adata.n_obs > max_cells:
            rng = np.random.default_rng(0)
            adata = adata[np.sort(rng.choice(adata.n_obs, max_cells, replace=False))].copy()
        adata, _, _ = prepare_for_standard(adata, copy=True)
        present = [g for g in dict.fromkeys(genes) if g in adata.var_names]
        if not present:
            return None
        X = adata[:, present].X
        X = X.toarray() if hasattr(X, 'toarray') else np.asarray(X)
        def stats(M):
            return {'mean': [round(float(v), 5) for v in M.mean(axis=0)], 'var': [round(float(v), 6) for v in M.var(axis=0)]}
        out = {'genes': present, 'all': stats(X), 'n_cells': int(X.shape[0])}
        if split_by and split_by in adata.obs.columns:
            out['groups'] = {}
            lab = adata.obs[split_by].astype(str).values
            for g in dict.fromkeys(lab):
                out['groups'][g] = stats(X[lab == g])
        return out
    except Exception as err:
        logger.warning(f"Global gene statistics skipped: {err}")
        return None


def _module_strength(result_dir: Path, task_id: str, module: str, pairs: list[str], lris: list[str]) -> Optional[np.ndarray]:
    files = glob.glob(str(result_dir / f'{task_id}_*_{module}_interactions_strength.tsv'))
    if not files:
        return None
    M = pd.read_csv(files[0], sep='\t', index_col=0)
    ri, ci = M.index.get_indexer(pairs), M.columns.get_indexer(lris)
    vals = np.full(len(pairs), np.nan); ok = (ri >= 0) & (ci >= 0)
    vals[ok] = M.values[ri[ok], ci[ok]]
    return vals


def effect_sizes(cs: np.ndarray, lr_items: list[dict], lr_idx: list[int], senders: list[int], receivers: list[int],
                 celltype_counts: list[int], gstats: Optional[dict]) -> list[Optional[float]]:
    """Standardised effect size of the Mean/Average/Arithmetic module: (CS - null mean) / null sd.

    Under FastCCC's null the per-cell-type mean of a gene is the mean of n i.i.d. draws from the
    gene's distribution over all cells (mean mu, variance s2), a complex is the average of its
    subunits and the LR score the arithmetic mean of the two sides."""
    if gstats is None:
        return [None] * len(lr_idx)
    pos = {g: i for i, g in enumerate(gstats['genes'])}
    mu, s2 = gstats['mean'], gstats['var']
    def side(genes: str, n: int):
        ids = [pos.get(g.strip()) for g in genes.split(',')]
        if any(i is None for i in ids) or n <= 0:
            return None
        k = len(ids)
        return sum(mu[i] for i in ids) / k, sum(s2[i] for i in ids) / (k * k * n)
    out = []
    for i, k in enumerate(lr_idx):
        if cs[i] is None or (isinstance(cs[i], float) and math.isnan(cs[i])):
            out.append(None); continue
        L = side(lr_items[k]['ligand'], celltype_counts[senders[i]]); Rr = side(lr_items[k]['receptor'], celltype_counts[receivers[i]])
        if L is None or Rr is None:
            out.append(None); continue
        m = (L[0] + Rr[0]) / 2; v = (L[1] + Rr[1]) / 4
        out.append(round(float((cs[i] - m) / math.sqrt(v)), 2) if v > 0 else None)
    return out


def bootstrap_support(result_dir: Path, sig: pd.DataFrame) -> Optional[list]:
    path = Path(result_dir) / 'bootstrap_support.tsv'
    if not path.is_file():
        return None
    bs = pd.read_csv(path, sep='\t')
    key = lambda df: list(zip(df['sender_celltype'].astype(str), df['receiver_celltype'].astype(str), df['LRI_ID'].astype(str)))
    sup = dict(zip(key(bs), bs['support']))
    return [None if sup.get(k) is None else round(float(sup[k]), 3) for k in key(sig)]


def replication_support(runs: list, min_groups: int = 3, top_k: int = 400) -> Optional[dict]:
    """For split runs with several groups: fraction of groups in which each (pair, LR) is significant,
    plus a samples x interactions score matrix for the most recurrent interactions."""
    groups = [(name, r) for name, r in runs if name != 'all']
    if len(groups) < min_groups:
        return None
    present: dict[tuple, dict] = {}
    for gname, r in groups:
        sig = r.significant
        strength = r.strength
        pairs = sig['sender_celltype'].astype(str) + '|' + sig['receiver_celltype'].astype(str)
        ri, ci = strength.index.get_indexer(pairs), strength.columns.get_indexer(sig['LRI_ID'].astype(str))
        for j, (p_, l_) in enumerate(zip(pairs, sig['LRI_ID'].astype(str))):
            v = float(strength.values[ri[j], ci[j]]) if ri[j] >= 0 and ci[j] >= 0 else None
            present.setdefault((p_, l_), {})[gname] = v
    n = len(groups)
    frac = {k: len(v) / n for k, v in present.items()}
    top = sorted(present, key=lambda k: (-len(present[k]), k))[:top_k]
    matrix = [[(present[k].get(g) if present[k].get(g) is not None else 0.0) for k in top] for g, _ in groups]
    return {'n_groups': n, 'groups': [g for g, _ in groups], 'fraction': {f'{k[0]}\t{k[1]}': round(v, 3) for k, v in frac.items()},
            'top': [[k[0], k[1]] for k in top], 'matrix': [[round(v, 3) for v in row] for row in matrix]}


class _Index:
    def __init__(self):
        self.items: list = []
        self.pos: dict = {}

    def add(self, key, make):
        if key not in self.pos:
            self.pos[key] = len(self.items)
            self.items.append(make())
        return self.pos[key]


def _module_pvalues(result_dir: Path, task_id: str, pairs: list[str], lris: list[str]) -> Optional[dict]:
    """Per-module p-values (16 scoring modules) for the given (pair, LRI) entries of a standard run."""
    files = sorted(p for p in glob.glob(str(result_dir / f'{task_id}_*_pvals.tsv')) if 'Cauchy' not in os.path.basename(p))
    if not files:
        return None
    names, cols = [], []
    for f in files:
        base = os.path.basename(f)
        method = re.sub(rf'^{re.escape(task_id)}_\d{{8}}_\d{{6}}_', '', base).replace('_pvals.tsv', '')
        try:
            M = pd.read_csv(f, sep='\t', index_col=0)
        except Exception as err:
            logger.warning(f"could not read {base}: {err}")
            continue
        ri = M.index.get_indexer(pairs)
        ci = M.columns.get_indexer(lris)
        vals = np.full(len(pairs), np.nan)
        ok = (ri >= 0) & (ci >= 0)
        vals[ok] = M.values[ri[ok], ci[ok]]
        names.append(method)
        cols.append(vals)
    if not names:
        return None
    P = np.vstack(cols).T   # rows x modules
    return {'names': names, 'p': [[_pval(v) for v in row] for row in P],
            'n_sig': [int(np.nansum(row < 0.05)) for row in P]}


# ── standard run ─────────────────────────────────────────────────────────────

def _standard_rows(res: RunResult, ct: _Index, lr: _Index, ann: dict, cond: Optional[int] = None) -> dict:
    sig = res.significant.copy()
    strength = res.strength
    sig['pair'] = sig['sender_celltype'].astype(str) + '|' + sig['receiver_celltype'].astype(str)
    ri = strength.index.get_indexer(sig['pair'])
    ci = strength.columns.get_indexer(sig['LRI_ID'].astype(str))
    cs = np.full(len(sig), np.nan)
    ok = (ri >= 0) & (ci >= 0)
    cs[ok] = strength.values[ri[ok], ci[ok]]
    rows = {'sender': [], 'receiver': [], 'lr': [], 'p': [], 'cs': []}
    if cond is not None:
        rows['cond'] = []
    for i, r in enumerate(sig.itertuples(index=False)):
        s = ct.add(str(r.sender_celltype), lambda n=str(r.sender_celltype): {'name': n})
        t = ct.add(str(r.receiver_celltype), lambda n=str(r.receiver_celltype): {'name': n})
        lid = str(r.LRI_ID)
        k = lr.add(lid, lambda lid=lid, r=r: {'id': lid, 'ligand': str(r.ligand), 'receptor': str(r.receptor),
                                            **ann.get(lid, {'classification': None, 'directionality': None})})
        rows['sender'].append(s); rows['receiver'].append(t); rows['lr'].append(k)
        rows['p'].append(_pval(sig.iloc[i]['p-value']))
        rows['cs'].append(_round(cs[i]))
        if cond is not None:
            rows['cond'].append(cond)
    modules = _module_pvalues(res.result_dir, res.task_id, list(sig['pair']), list(sig['LRI_ID'].astype(str)))
    maa = _module_strength(res.result_dir, res.task_id, EFFECT_MODULE, list(sig['pair']), list(sig['LRI_ID'].astype(str)))
    rows['cs_maa'] = [_round(v) for v in maa] if maa is not None else [None] * len(sig)
    rows['bs'] = bootstrap_support(res.result_dir, sig) or [None] * len(sig)
    return rows, modules


class ExplorerOptions:
    def __init__(self, genesets: bool = True, cells: bool = True, max_cells: int = 4000, patterns: bool = True,
                 n_patterns: Optional[int] = None, nulls: bool = True, include_both_ns: bool = True):
        self.genesets, self.cells, self.max_cells, self.patterns = genesets, cells, max_cells, patterns
        self.n_patterns, self.nulls, self.include_both_ns = n_patterns, nulls, include_both_ns


def _extras(payload: dict, inp: dict, opts: 'ExplorerOptions', result_dir: Optional[Path] = None, db_path=None) -> dict:
    """Attach gene sets, communication patterns, per-cell data and (reference) null curves."""
    genes = _genes_of(payload['lr'])
    species = (payload['meta'].get('database_info') or {}).get('species') or 'human'
    if opts.genesets and species != 'human':
        logger.info(f"The Enrichr gene-set libraries are human; skipped for a {species} database.")
    payload['genesets'] = gene_sets(genes) if opts.genesets and species == 'human' else None
    # communication patterns (NMF) are computed in the page on the current selection; the payload only carries the option
    payload['patterns'] = {'k': opts.n_patterns} if opts.patterns else None
    payload['cells'] = cell_level_data(inp.get('source'), inp.get('celltype_key'), genes, max_cells=opts.max_cells) if opts.cells else None
    if payload['kind'] == 'reference' and opts.nulls and result_dir is not None:
        payload['nulls'] = null_curves(result_dir, payload['rows'], payload['celltypes'], payload['lr'], db_path)
    else:
        payload['nulls'] = None
    payload['calibration'] = calibration_diagnostics(result_dir) if (payload['kind'] == 'reference' and result_dir is not None) else None
    return payload


def calibration_diagnostics(result_dir: Path) -> Optional[dict]:
    """Query-reference comparability: calibration factor k and housekeeping-gene rank agreement."""
    cal = Path(result_dir) / 'query_calibration.json'
    if not cal.is_file():
        return None
    try:
        out = json.loads(cal.read_text(encoding='utf-8'))
        hk = Path(result_dir) / 'query_hk_ranks.tsv'
        if hk.is_file():
            df = pd.read_csv(hk, sep='\t')
            df = df.dropna()
            if len(df) > 1500:
                df = df.sample(1500, random_state=0)
            out['hk'] = {'gene': df['gene'].astype(str).tolist(), 'query': [round(float(v), 3) for v in df['query_hk']],
                         'ref': [round(float(v), 3) for v in df['ref_hk']]}
        return out
    except Exception as err:
        logger.warning(f"calibration diagnostics skipped: {err}")
        return None


def build_standard_payload(res: RunResult, opts: Optional['ExplorerOptions'] = None) -> dict:
    opts = opts or ExplorerOptions()
    ann = db_annotations(res.database_path)
    ct, lr = _Index(), _Index()
    counts = (res.manifest.get('input', {}) or {}).get('celltype_counts') or {}
    for name in counts:
        ct.add(str(name), lambda n=str(name): {'name': n})
    rows, modules = _standard_rows(res, ct, lr, ann)
    names = [c['name'] for c in ct.items]
    colors = celltype_colors(names)
    celltypes = [{'name': n, 'n_cells': int(counts.get(n, 0)), 'color': colors[n]} for n in names]
    inp = res.manifest.get('input', {}) or {}
    genes = gene_expression_stats(inp.get('source'), inp.get('celltype_key'), _genes_of(lr.items))
    gstats = gene_global_stats(inp.get('source'), _genes_of(lr.items))
    rows['z'] = effect_sizes(rows['cs_maa'], lr.items, rows['lr'], rows['sender'], rows['receiver'], [c['n_cells'] for c in celltypes],
                             gstats['all'] | {'genes': gstats['genes']} if gstats else None)
    payload = {
        'kind': 'standard',
        'meta': {
            'name': res.name, 'database': res.database, 'result_dir': str(res.result_dir), 'task_id': res.task_id,
            'database_info': db_metadata(res.database_path),
            'fastccc_version': res.manifest.get('fastccc_version'), 'created': res.manifest.get('created'),
            'n_cells': inp.get('n_cells'), 'n_genes': inp.get('n_genes'), 'celltype_key': inp.get('celltype_key'),
            'matrix': inp.get('matrix_detected'), 'preprocessing': inp.get('preprocessing_applied'),
            'parameters': res.manifest.get('parameters'), 'timing': res.manifest.get('timing_seconds'),
            'n_celltype_pairs_tested': int(res.strength.shape[0]) if res.strength is not None else None,
            'n_lr_tested': int(res.strength.shape[1]) if res.strength is not None else None,
        },
        'celltypes': celltypes, 'lr': lr.items, 'rows': rows, 'modules': modules, 'genes': genes,
        'bootstrap': res.manifest.get('bootstrap'),
    }
    return _extras(payload, inp, opts, db_path=res.database_path)


def build_split_payload(multi: MultiRunResult, opts: Optional['ExplorerOptions'] = None) -> dict:
    opts = opts or ExplorerOptions()
    runs = []
    if multi.all is not None:
        runs.append(('all', multi.all))
    runs.extend((g, r) for g, r in multi.runs.items())
    ann = db_annotations(runs[0][1].database_path)
    ct, lr = _Index(), _Index()
    for _, r in runs:
        for name in ((r.manifest.get('input', {}) or {}).get('celltype_counts') or {}):
            ct.add(str(name), lambda n=str(name): {'name': n})
    rows = {'sender': [], 'receiver': [], 'lr': [], 'p': [], 'cs': [], 'cond': [], 'cs_maa': [], 'bs': []}
    module_rows = []
    for ci, (cond, r) in enumerate(runs):
        part, modules = _standard_rows(r, ct, lr, ann, cond=ci)
        for k in rows:
            rows[k].extend(part[k])
        module_rows.append(modules)
    # modules: concatenate when every run has the same module names
    modules = None
    if all(m is not None for m in module_rows) and len({tuple(m['names']) for m in module_rows}) == 1:
        modules = {'names': module_rows[0]['names'], 'p': sum((m['p'] for m in module_rows), []),
                   'n_sig': sum((m['n_sig'] for m in module_rows), [])}
    names = [c['name'] for c in ct.items]
    colors = celltype_colors(names)
    counts = {}
    for _, r in runs:
        for n, v in ((r.manifest.get('input', {}) or {}).get('celltype_counts') or {}).items():
            counts[n] = max(counts.get(n, 0), int(v))
    celltypes = [{'name': n, 'n_cells': counts.get(n, 0), 'color': colors[n]} for n in names]
    inp0 = runs[0][1].manifest.get('input', {}) or {}
    genes = gene_expression_stats(inp0.get('source'), inp0.get('celltype_key'), _genes_of(lr.items))
    gstats = gene_global_stats(inp0.get('source'), _genes_of(lr.items), split_by=multi.split_by)
    # per-condition cell-type sizes and null moments
    z_all = [None] * len(rows['lr'])
    if gstats:
        for ci, (cond, r) in enumerate(runs):
            cnt = (r.manifest.get('input', {}) or {}).get('celltype_counts') or {}
            sizes = [int(cnt.get(c['name'], 0)) for c in celltypes]
            gs = gstats['all'] if cond == 'all' or cond not in (gstats.get('groups') or {}) else gstats['groups'][cond]
            idx = [i for i in range(len(rows['lr'])) if rows['cond'][i] == ci]
            zz = effect_sizes([rows['cs_maa'][i] for i in idx], lr.items, [rows['lr'][i] for i in idx], [rows['sender'][i] for i in idx],
                              [rows['receiver'][i] for i in idx], sizes, gs | {'genes': gstats['genes']})
            for i, v in zip(idx, zz):
                z_all[i] = v
    rows['z'] = z_all
    rep = replication_support(runs)
    if rep:
        # replication of the full-dataset ('all') rows across the groups
        rows['rep'] = [None] * len(rows['lr'])
        if runs[0][0] == 'all':
            for i in range(len(rows['lr'])):
                if rows['cond'][i] == 0:
                    key = f"{celltypes[rows['sender'][i]]['name']}|{celltypes[rows['receiver'][i]]['name']}\t{lr.items[rows['lr'][i]]['id']}"
                    rows['rep'][i] = rep['fraction'].get(key, 0.0)
    payload = {
        'kind': 'split',
        'meta': {
            'name': multi.manifest.get('name'), 'split_by': multi.split_by, 'database': runs[0][1].database,
            'database_info': db_metadata(runs[0][1].database_path),
            'result_dir': str(multi.result_dir), 'fastccc_version': multi.manifest.get('fastccc_version'),
            'created': multi.manifest.get('created'), 'celltype_key': inp0.get('celltype_key'),
            'conditions': [{'name': c, 'n_cells': (r.manifest.get('input', {}) or {}).get('n_cells'),
                            'n_significant': int(len(r.significant)), 'result_dir': str(r.result_dir)} for c, r in runs],
            'parameters': multi.manifest.get('parameters'),
        },
        'celltypes': celltypes, 'lr': lr.items, 'rows': rows, 'modules': modules, 'genes': genes,
        'replication': ({'n_groups': rep['n_groups'], 'groups': rep['groups'], 'top': rep['top'], 'matrix': rep['matrix']} if rep else None),
    }
    return _extras(payload, inp0, opts, db_path=runs[0][1].database_path)


# ── reference comparison ─────────────────────────────────────────────────────

def build_reference_payload(res: ReferenceResult, opts: Optional['ExplorerOptions'] = None) -> dict:
    opts = opts or ExplorerOptions()
    include_both_ns = opts.include_both_ns
    ann = db_annotations(res.database_path)
    df = res.results.copy()
    if not include_both_ns:
        df = df[df['trend_vs_ref'] != 'Both NS']
    ct, lr = _Index(), _Index()
    inp = res.manifest.get('input', {}) or {}
    counts = inp.get('celltype_counts') or {}
    for name in counts:
        ct.add(str(name), lambda n=str(name): {'name': n})
    sr = df['sender|receiver'].astype(str).str.split('|', n=1, expand=True)
    cols = {
        'sender': [], 'receiver': [], 'lr': [], 'trend': [], 'cs': [], 'ref_cs': [], 'null_cs': [],
        'thr_lo': [], 'thr_hi': [], 'sig': [], 'sig_ref': [], 'expr': [], 'expr_ref': [],
        'lig_pct': [], 'rec_pct': [], 'lig_pct_ref': [], 'rec_pct_ref': [],
        'lig_cs': [], 'lig_lo': [], 'lig_hi': [], 'rec_cs': [], 'rec_lo': [], 'rec_hi': [],
    }
    has_ref_cs = 'ref_comm_score' in df.columns
    for i, r in enumerate(df.itertuples(index=False)):
        rec = df.iloc[i]
        s = ct.add(sr.iloc[i, 0], lambda n=sr.iloc[i, 0]: {'name': n})
        t = ct.add(sr.iloc[i, 1], lambda n=sr.iloc[i, 1]: {'name': n})
        lid = str(rec['LRI_ID'])
        k = lr.add(lid, lambda lid=lid, rec=rec: {'id': lid, 'ligand': str(rec['ligand']), 'receptor': str(rec['receptor']),
                                                **ann.get(lid, {'classification': None, 'directionality': None})})
        trend = rec['trend_vs_ref']
        tl, th = _split_range(rec.get('sig_threshold_CI'))
        ll, lh = _split_range(rec.get('ligand_CS_CI'))
        rl, rh = _split_range(rec.get('receptor_CS_CI'))
        cols['sender'].append(s); cols['receiver'].append(t); cols['lr'].append(k)
        cols['trend'].append(TREND_CODES.get(trend, 0) if isinstance(trend, str) else 0)
        cols['cs'].append(_round(rec.get('comm_score')))
        cols['ref_cs'].append(_round(rec.get('ref_comm_score')) if has_ref_cs else None)
        cols['null_cs'].append(_round(rec.get('null_comm_score')))
        cols['thr_lo'].append(tl); cols['thr_hi'].append(th)
        cols['sig'].append(bool(rec.get('is_significant')) if not pd.isna(rec.get('is_significant')) else None)
        v = rec.get('is_significant_ref')
        cols['sig_ref'].append(None if pd.isna(v) else bool(v))
        cols['expr'].append(bool(rec.get('above_expr_threshold')) if not pd.isna(rec.get('above_expr_threshold')) else None)
        v = rec.get('above_expr_threshold_ref')
        cols['expr_ref'].append(None if pd.isna(v) else bool(v))
        cols['lig_pct'].append(_round(rec.get('ligand_expr_percent'))); cols['rec_pct'].append(_round(rec.get('receptor_expr_percent')))
        cols['lig_pct_ref'].append(_round(rec.get('ligand_expr_percent_ref'))); cols['rec_pct_ref'].append(_round(rec.get('receptor_expr_percent_ref')))
        cols['lig_cs'].append(_round(rec.get('ligand_CS_component'))); cols['lig_lo'].append(ll); cols['lig_hi'].append(lh)
        cols['rec_cs'].append(_round(rec.get('receptor_CS_component'))); cols['rec_lo'].append(rl); cols['rec_hi'].append(rh)
    names = [c['name'] for c in ct.items]
    colors = celltype_colors(names)
    celltypes = [{'name': n, 'n_cells': int(counts.get(n, 0)), 'color': colors[n]} for n in names]
    genes = gene_expression_stats(inp.get('source'), inp.get('celltype_key'), _genes_of(lr.items))
    ref = res.manifest.get('reference', {}) or {}
    payload = {
        'kind': 'reference',
        'meta': {
            'name': res.query_name, 'database': res.database, 'result_dir': str(res.result_dir),
            'database_info': db_metadata(res.database_path),
            'fastccc_version': res.manifest.get('fastccc_version'), 'created': res.manifest.get('created'),
            'n_cells': inp.get('n_cells'), 'n_genes': inp.get('n_genes'), 'celltype_key': inp.get('celltype_key'),
            'matrix': inp.get('matrix_detected'), 'preprocessing': inp.get('preprocessing_applied'),
            'timing': res.manifest.get('timing_seconds'),
            'reference': {
                'name': ref.get('name'), 'display_name': res.manifest.get('reference_display_name'),
                'source': ref.get('source'), 'tissue': ref.get('tissue'), 'path': ref.get('path'),
                'n_cells': ref.get('n_cells'), 'n_celltypes': ref.get('n_celltypes'),
                'celltypes': ref.get('celltypes'), 'built_from': ref.get('built_from'),
                'min_percentile': ref.get('min_percentile'),
            },
            'celltype_mapping': res.manifest.get('celltype_mapping') or {},
            'celltype_mapping_source': res.manifest.get('celltype_mapping_source'),
            'celltype_coverage': res.manifest.get('celltype_coverage'),
            'trend_counts': res.trend_counts, 'has_ref_cs': has_ref_cs,
        },
        'celltypes': celltypes, 'lr': lr.items, 'rows': cols, 'modules': None, 'genes': genes,
    }
    return _extras(payload, inp, opts, result_dir=res.result_dir, db_path=res.database_path)


def build_payload(result_dir, opts: Optional['ExplorerOptions'] = None, **kwargs) -> dict:
    opts = opts or ExplorerOptions(**kwargs)
    res = load_result(result_dir)
    if isinstance(res, ReferenceResult):
        return build_reference_payload(res, opts)
    if isinstance(res, MultiRunResult):
        return build_split_payload(res, opts)
    return build_standard_payload(res, opts)


# ── richer analyses: null curves, gene sets, communication patterns, per-cell data ───────────

import base64
import json
import urllib.request

ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={name}'
DEFAULT_GENESETS = ('MSigDB_Hallmark_2020', 'KEGG_2021_Human', 'GO_Biological_Process_2023')


def lri_multidata(db_path) -> dict[str, tuple[int, int]]:
    """``{CPI id: (ligand multidata id, receptor multidata id)}``."""
    table = pd.read_csv(Path(db_path) / 'interaction_table.csv')
    return {str(r['id_cp_interaction']): (int(r['multidata_1_id']), int(r['multidata_2_id'])) for _, r in table.iterrows()}


def load_null_pmfs(result_dir: Path) -> Optional[dict]:
    path = Path(result_dir) / 'query_null_pmfs.npz'
    if not path.is_file():
        return None
    try:
        z = np.load(path, allow_pickle=False)
        keys, offsets, lengths, values = z['keys'], z['offsets'], z['lengths'], z['values']
        step = float(z['grid_step'][0]) if 'grid_step' in z else 0.01
        out, pos = {}, 0
        for k, off, n in zip(keys, offsets, lengths):
            out[str(k)] = (int(off), values[pos:pos + int(n)])
            pos += int(n)
        return {'step': step, 'pmfs': out}
    except Exception as err:
        logger.warning(f"could not read {path}: {err}")
        return None


def _downsample_pmf(offset: int, vals: np.ndarray, step: float, n_bins: int = 80) -> dict:
    """Sum-pool a fine pmf into at most ``n_bins`` bins; returns x0, dx and bin probabilities."""
    n = len(vals)
    if n <= n_bins:
        return {'x0': round(offset * step, 4), 'dx': step, 'y': [round(float(v), 6) for v in vals]}
    f = int(math.ceil(n / n_bins))
    pad = (-n) % f
    arr = np.concatenate([vals, np.zeros(pad)]) if pad else vals
    pooled = arr.reshape(-1, f).sum(axis=1)
    return {'x0': round(offset * step, 4), 'dx': round(step * f, 6), 'y': [round(float(v), 6) for v in pooled]}


def null_curves(result_dir: Path, rows: dict, celltypes: list[dict], lr_items: list[dict], db_path,
                max_entries: int = 60000) -> Optional[dict]:
    """Null pmf of the ligand side (sender) and receptor side (receiver) for every informative row."""
    nulls = load_null_pmfs(result_dir)
    if nulls is None:
        return None
    mids = lri_multidata(db_path)
    needed = set()
    for i in range(len(rows['sender'])):
        if rows['trend'][i] in (1, 2, 3):
            m = mids.get(lr_items[rows['lr'][i]]['id'])
            if not m:
                continue
            needed.add((rows['sender'][i], m[0])); needed.add((rows['receiver'][i], m[1]))
    if len(needed) > max_entries:
        logger.info(f"null curves limited to {max_entries:,} of {len(needed):,} (cell type, gene) pairs")
        needed = set(sorted(needed)[:max_entries])
    curves = {}
    for ct_idx, mid in needed:
        key = f"{celltypes[ct_idx]['name']}\t{mid}"
        hit = nulls['pmfs'].get(key)
        if hit is None:
            continue
        off, vals = hit
        curves[f"{ct_idx}|{mid}"] = _downsample_pmf(off, vals, nulls['step'])
    return {'curves': curves, 'lr_multidata': {lr['id']: list(mids[lr['id']]) for lr in lr_items if lr['id'] in mids}}


def _fetch_gmt(name: str, cache_dir: Path) -> Optional[dict[str, list[str]]]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f'{name}.gmt'
    if not path.is_file():
        try:
            urllib.request.urlretrieve(ENRICHR_URL.format(name=name), path)
        except Exception as err:
            logger.warning(f"gene set library {name} not available (offline?): {err}")
            return None
    sets = {}
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            genes = [g.split(',')[0] for g in parts[2:] if g]
            sets[parts[0]] = genes
    return sets


def gene_sets(universe: list[str], libraries=DEFAULT_GENESETS, cache_dir: Optional[Path] = None,
              min_overlap: int = 3) -> Optional[dict]:
    """Gene-set libraries pruned to the ligand/receptor universe, for in-page enrichment tests."""
    from ..resources import fastccc_home
    cache_dir = Path(cache_dir) if cache_dir else fastccc_home() / 'genesets'
    uni = list(dict.fromkeys(universe)); pos = {g: i for i, g in enumerate(uni)}
    libs = []
    for name in libraries:
        sets = _fetch_gmt(name, cache_dir)
        if not sets:
            continue
        pruned = []
        for sname, genes in sets.items():
            idx = sorted({pos[g] for g in genes if g in pos})
            if len(idx) >= min_overlap:
                pruned.append({'name': sname, 'genes': idx, 'size': len(genes)})
        if pruned:
            libs.append({'name': name, 'sets': pruned})
    if not libs:
        return None
    return {'genes': uni, 'libraries': libs}


def _nmf(X: np.ndarray, k: int, n_iter: int = 400, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    n, m = X.shape
    W = rng.random((n, k)) + 1e-3
    H = rng.random((k, m)) + 1e-3
    eps = 1e-9
    for _ in range(n_iter):
        H *= (W.T @ X) / (W.T @ W @ H + eps)
        W *= (X @ H.T) / (W @ H @ H.T + eps)
    return W, H


def communication_patterns(rows: dict, n_celltypes: int, n_lr: int, kind: str, k: Optional[int] = None,
                           min_rows: int = 30) -> Optional[dict]:
    """NMF of the cell type x LR matrices (outgoing: senders; incoming: receivers) of significant rows."""
    idx = [i for i in range(len(rows['sender']))
           if (rows['trend'][i] in (1, 3) if kind == 'reference' else (rows.get('cond') is None or rows['cond'][i] == 0))]
    if len(idx) < min_rows or n_celltypes < 3:
        return None
    out = np.zeros((n_celltypes, n_lr)); inc = np.zeros((n_celltypes, n_lr))
    for i in idx:
        w = rows['cs'][i] if rows['cs'][i] is not None else 1.0
        out[rows['sender'][i], rows['lr'][i]] += w
        inc[rows['receiver'][i], rows['lr'][i]] += w
    k = k or int(min(6, max(2, round(math.sqrt(n_celltypes)))))
    result = {'k': k, 'outgoing': None, 'incoming': None}
    for name, M in (('outgoing', out), ('incoming', inc)):
        keep_ct = np.where(M.sum(axis=1) > 0)[0]; keep_lr = np.where(M.sum(axis=0) > 0)[0]
        if len(keep_ct) < 3 or len(keep_lr) < 3:
            continue
        X = M[np.ix_(keep_ct, keep_lr)]
        X = X / (X.max() or 1)
        kk = min(k, len(keep_ct) - 1, len(keep_lr) - 1)
        W, H = _nmf(X, kk)
        # normalise: each cell type's loadings sum to 1; each pattern's LR profile sums to 1
        Wn = W / (W.sum(axis=1, keepdims=True) + 1e-12)
        Hn = H / (H.sum(axis=1, keepdims=True) + 1e-12)
        result[name] = {
            'celltypes': [int(c) for c in keep_ct], 'lr': [int(l) for l in keep_lr],
            'W': [[round(float(v), 4) for v in row] for row in Wn],
            'H': [[round(float(v), 5) for v in row] for row in Hn],
        }
    return result if (result['outgoing'] or result['incoming']) else None


def cell_level_data(data_path: Optional[str], celltype_key: Optional[str], genes: list[str],
                    max_cells: int = 4000, min_per_type: int = 20) -> Optional[dict]:
    """Subsampled per-cell expression (uint8-quantised, base64) and a 2-D embedding for violins and UMAP overlays."""
    if not data_path or not celltype_key or not Path(data_path).is_file():
        return None
    try:
        adata = load_anndata(data_path)
        if celltype_key not in adata.obs.columns:
            return None
        labels = adata.obs[celltype_key].astype(str).values
        rng = np.random.default_rng(0)
        if adata.n_obs > max_cells:   # stratified: proportional, but keep small types visible
            cts, counts = np.unique(labels, return_counts=True)
            quota = {c: max(min_per_type, int(round(max_cells * n / adata.n_obs))) for c, n in zip(cts, counts)}
            pick = []
            for c in cts:
                where = np.where(labels == c)[0]
                pick.extend(where if len(where) <= quota[c] else rng.choice(where, quota[c], replace=False))
            pick = np.sort(np.array(pick))
            adata = adata[pick].copy()
            labels = labels[pick]
        adata, _, _ = prepare_for_standard(adata, copy=True)
        present = [g for g in dict.fromkeys(genes) if g in adata.var_names]
        if not present:
            return None
        X = adata[:, present].X
        X = X.toarray() if hasattr(X, 'toarray') else np.asarray(X)
        gmax = X.max(axis=0); gmax[gmax == 0] = 1.0
        q = np.clip(np.round(X / gmax * 255), 0, 255).astype(np.uint8)
        emb_key = next((k for k in ('X_umap', 'X_tsne', 'X_pca') if k in adata.obsm), None)
        emb = None
        if emb_key is not None:
            E = np.asarray(adata.obsm[emb_key])[:, :2].astype(np.float32)
            emb = {'key': emb_key, 'x': [round(float(v), 2) for v in E[:, 0]], 'y': [round(float(v), 2) for v in E[:, 1]]}
        ct_names = list(dict.fromkeys(labels)); ct_pos = {c: i for i, c in enumerate(ct_names)}
        return {
            'n_cells': int(adata.n_obs), 'genes': present, 'gene_max': [round(float(v), 3) for v in gmax],
            'celltypes': ct_names, 'cell_ct': [ct_pos[c] for c in labels],
            'expr_b64': base64.b64encode(np.ascontiguousarray(q).tobytes()).decode('ascii'),   # row-major cells x genes
            'embedding': emb,
        }
    except Exception as err:
        logger.warning(f"Per-cell data skipped: {err}")
        return None
