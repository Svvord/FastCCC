"""Inspect AnnData inputs: matrix type detection and cell-type column discovery.

FastCCC's standard workflow expects normalised ``log1p`` data, so raw counts
are normalised automatically.  The reference workflows (``build_reference`` /
``infer_query``) rank genes within each cell, which makes them invariant to
per-cell normalisation and log1p; they accept any non-negative matrix.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from loguru import logger

RAW_LAYER_CANDIDATES = ('counts', 'raw_counts', 'raw', 'X_counts', 'umi')

# user-facing names accepted by ``input_type`` / ``--input-type``
INPUT_TYPE_ALIASES = {
    'raw': 'raw_counts', 'counts': 'raw_counts', 'raw_counts': 'raw_counts',
    'normalized': 'normalized', 'normalised': 'normalized',
    'log1p': 'log1p', 'log': 'log1p', 'lognorm': 'log1p',
}

# thresholds for the matrix-type heuristics (see detect_matrix_type)
_INTEGER_FRACTION_MIN = 0.999   # share of non-zero values that are whole numbers -> counts
_LOG1P_MAX = 30.0               # log1p data never exceeds this; counts almost always do
_CONSTANT_ROWSUM_CV = 0.01      # normalize_total makes per-cell totals constant
_LIBRARY_SIZE_CV = 0.2          # genuine library sizes vary far more than this


@dataclass
class MatrixInfo:
    kind: str                    # raw_counts | log1p | normalized | scaled | empty
    is_integer: bool
    min_value: float
    max_value: float
    n_cells_sampled: int
    source: str = 'X'            # X | layers:<name> | raw
    confidence: str = 'high'     # high | medium | low
    integer_fraction: float = 1.0
    rowsum_cv: Optional[float] = None        # CV of per-cell totals of the matrix
    expm1_rowsum_cv: Optional[float] = None  # CV of per-cell totals after expm1 (log1p check)
    reason: str = ''
    user_override: bool = False

    def to_dict(self) -> dict:
        return asdict(self)

    def describe(self) -> str:
        bits = [f"kind={self.kind}", f"confidence={self.confidence}", f"max={self.max_value:.2f}",
                f"integer_fraction={self.integer_fraction:.3f}"]
        if self.rowsum_cv is not None:
            bits.append(f"rowsum_cv={self.rowsum_cv:.3f}")
        if self.expm1_rowsum_cv is not None:
            bits.append(f"expm1_rowsum_cv={self.expm1_rowsum_cv:.3f}")
        return ', '.join(bits)


def _sample_rows(X, n_cells: int, seed: int):
    n = X.shape[0]
    if n > n_cells:
        rng = np.random.default_rng(seed)
        rows = np.sort(rng.choice(n, size=n_cells, replace=False))
        sub = X[rows]
    else:
        sub = X
    if sp.issparse(sub):
        sub = sp.csr_matrix(sub)
        vals = np.asarray(sub.data, dtype=np.float64)
    else:
        sub = np.asarray(sub, dtype=np.float64)
        vals = sub[sub != 0]
    return sub, vals


def _cv(row_sums: np.ndarray) -> float:
    row_sums = row_sums[row_sums > 0]
    if row_sums.size < 2:
        return 0.0
    return float(row_sums.std() / row_sums.mean())


def detect_matrix_type(X, n_cells: int = 2000, seed: int = 0, source: str = 'X') -> MatrixInfo:
    """Classify an expression matrix by sampling up to ``n_cells`` cells.

    Rules, in order:

    1. negative values -> ``scaled`` (z-scored / batch-corrected; unusable);
    2. (almost) all non-zero values are whole numbers -> ``raw_counts``;
    3. values <= 30 and per-cell totals constant after ``expm1`` -> ``log1p``;
    4. per-cell totals constant -> ``normalized`` (``normalize_total`` without log);
    5. otherwise values > 30 with library-size-like variation of per-cell
       totals -> non-integer ``raw_counts`` (SoupX / alevin / RSEM style);
       small variation -> ``normalized`` (e.g. ``exclude_highly_expressed``);
       values <= 30 -> ``log1p``.

    Rules 1-4 are reported with ``confidence='high'``; rule 5 with ``medium`` or
    ``low``.  The evidence is stored on the returned :class:`MatrixInfo`.
    """
    if hasattr(X, 'X') and not hasattr(X, 'shape'):
        X = X.X
    if isinstance(X, anndata.AnnData):
        X = X.X
    sub, vals = _sample_rows(X, n_cells, seed)
    n_sampled = min(X.shape[0], n_cells)
    if vals.size == 0:
        return MatrixInfo('empty', False, 0.0, 0.0, n_sampled, source, 'high', 0.0, reason='no non-zero values')
    vmin, vmax = float(vals.min()), float(vals.max())
    if vmin < 0:
        return MatrixInfo('scaled', False, vmin, vmax, n_sampled, source, 'high', 0.0,
                          reason='negative values present (scaled / centred data)')

    int_frac = float(np.mean(np.isclose(vals, np.rint(vals), atol=1e-6, rtol=0.0)))
    is_int = int_frac >= _INTEGER_FRACTION_MIN
    rowsum_cv = _cv(np.asarray(sub.sum(axis=1), dtype=np.float64).ravel())
    expm1_cv = None
    if vmax <= _LOG1P_MAX:
        e = sub.expm1() if sp.issparse(sub) else np.expm1(sub)
        expm1_cv = _cv(np.asarray(e.sum(axis=1), dtype=np.float64).ravel())

    common = dict(is_integer=is_int, min_value=vmin, max_value=vmax, n_cells_sampled=n_sampled, source=source,
                  integer_fraction=int_frac, rowsum_cv=rowsum_cv, expm1_rowsum_cv=expm1_cv)
    if is_int:
        return MatrixInfo('raw_counts', confidence='high', reason='non-zero values are whole numbers', **common)
    if expm1_cv is not None and expm1_cv < _CONSTANT_ROWSUM_CV:
        return MatrixInfo('log1p', confidence='high',
                          reason='per-cell totals are constant after expm1 (normalize_total + log1p)', **common)
    if rowsum_cv < _CONSTANT_ROWSUM_CV:
        return MatrixInfo('normalized', confidence='high',
                          reason='per-cell totals are constant (normalize_total without log1p)', **common)
    if vmax > _LOG1P_MAX:
        if rowsum_cv >= _LIBRARY_SIZE_CV:
            return MatrixInfo('raw_counts', confidence='medium',
                              reason='non-integer values with library-size-like variation of per-cell totals '
                                     '(ambient-corrected or fractional counts)', **common)
        return MatrixInfo('normalized', confidence='medium',
                          reason='non-integer values with nearly constant per-cell totals '
                                 '(e.g. normalize_total with exclude_highly_expressed)', **common)
    conf = 'medium' if (expm1_cv is not None and expm1_cv < _LIBRARY_SIZE_CV) else 'low'
    return MatrixInfo('log1p', confidence=conf,
                      reason='values <= 30 but per-cell totals vary after expm1 (log1p without normalize_total?)',
                      **common)


def find_raw_counts(adata: anndata.AnnData) -> Optional[tuple[anndata.AnnData, MatrixInfo]]:
    """Return an AnnData whose ``X`` holds raw counts, looking in X, layers and .raw."""
    info = detect_matrix_type(adata.X)
    if info.kind == 'raw_counts':
        return adata, info
    for layer in RAW_LAYER_CANDIDATES:
        if layer in adata.layers:
            linfo = detect_matrix_type(adata.layers[layer], source=f'layers:{layer}')
            if linfo.kind == 'raw_counts':
                out = anndata.AnnData(X=adata.layers[layer], obs=adata.obs.copy(), var=adata.var.copy())
                out.obs_names = adata.obs_names
                out.var_names = adata.var_names
                return out, linfo
    if adata.raw is not None:
        rinfo = detect_matrix_type(adata.raw.X, source='raw')
        if rinfo.kind == 'raw_counts':
            out = adata.raw.to_adata()
            out.obs = adata.obs.copy()
            return out, rinfo
    return None


# ── cell-type column discovery ───────────────────────────────────────────────

_NEGATIVE_NAME_TOKENS = {
    'orig', 'sample', 'samples', 'sampleid', 'donor', 'donors', 'patient', 'patients', 'individual', 'subject',
    'batch', 'batches', 'library', 'libraries', 'lane', 'condition', 'conditions', 'treatment', 'disease',
    'status', 'sex', 'gender', 'age', 'phase', 'tissue', 'organ', 'region', 'site', 'assay', 'technology',
    'chemistry', 'platform', 'dataset', 'study', 'replicate', 'rep', 'timepoint', 'time', 'stage', 'species',
    'organism', 'barcode', 'barcodes', 'well', 'plate', 'id', 'doublet', 'scrublet', 'percent', 'pct', 'qc',
    'outlier', 'group', 'groups', 'project', 'experiment', 'run', 'source', 'protocol', 'suspension',
    'ethnicity', 'smoking', 'bmi', 'mt', 'ribo', 'genes', 'counts', 'umi', 'size', 'score', 'scores',
}
_ANNOTATION_TOKENS = {'annotation', 'annotations', 'annot', 'annotated'}
_TYPE_TOKENS = {'type', 'types', 'subtype', 'subtypes', 'lineage', 'lineages', 'class', 'ontology', 'identity'}
_LABEL_TOKENS = {'label', 'labels'}
_CLUSTER_TOKENS = {'cluster', 'clusters', 'clustering', 'leiden', 'louvain', 'ident'}

_CELL_VOCAB = {
    'cell', 'cells', 't', 'b', 'nk', 'nkt', 'treg', 'tregs', 'th1', 'th2', 'th17', 'tfh', 'cd4', 'cd8', 'ilc',
    'mait', 'macrophage', 'macrophages', 'mono', 'monocyte', 'monocytes', 'dc', 'dcs', 'pdc', 'cdc', 'dendritic',
    'neutrophil', 'neutrophils', 'eosinophil', 'basophil', 'mast', 'plasma', 'plasmablast', 'platelet', 'platelets',
    'erythrocyte', 'erythroid', 'megakaryocyte', 'megakaryocytes', 'hsc', 'hspc', 'progenitor', 'progenitors',
    'stem', 'fibroblast', 'fibroblasts', 'myofibroblast', 'endothelial', 'endothelium', 'endo', 'pericyte',
    'pericytes', 'smooth', 'muscle', 'epithelial', 'epithelium', 'epi', 'keratinocyte', 'hepatocyte',
    'cholangiocyte', 'enterocyte', 'goblet', 'paneth', 'tuft', 'enteroendocrine', 'acinar', 'ductal', 'alpha',
    'beta', 'delta', 'gamma', 'islet', 'neuron', 'neurons', 'neural', 'interneuron', 'astrocyte', 'astrocytes',
    'oligodendrocyte', 'opc', 'microglia', 'schwann', 'glia', 'glial', 'cardiomyocyte', 'podocyte', 'tubule',
    'mesangial', 'osteoblast', 'osteoclast', 'chondrocyte', 'adipocyte', 'melanocyte', 'trophoblast', 'evt',
    'syncytiotrophoblast', 'decidual', 'stromal', 'stroma', 'mesenchymal', 'myeloid', 'lymphoid', 'lymphocyte',
    'lymphocytes', 'leukocyte', 'granulocyte', 'granulocytes', 'naive', 'memory', 'effector', 'exhausted',
    'cycling', 'proliferating', 'tumor', 'tumour', 'cancer', 'malignant', 'epcam', 'basal', 'luminal',
    'secretory', 'ciliated', 'club', 'at1', 'at2', 'alveolar', 'immune', 'myocyte', 'satellite', 'germ',
    'sertoli', 'leydig', 'oocyte', 'retinal', 'photoreceptor', 'bipolar', 'amacrine', 'muller', 'rod', 'cone',
    'cardiac', 'vascular', 'lymphatic', 'ependymal', 'neuroblast', 'radial', 'unknown', 'doublets',
}
_CELL_SUFFIX_RE = re.compile(r'(cyte|cytes|blast|blasts|phage|phages|troph|trophs)$')
_MARKER_RE = re.compile(r'^cd\d+|^(il|ccr|cxcr|klr|fcgr|tcr|ig)[a-z0-9]+$')


def _name_tokens(name: str) -> list[str]:
    spaced = re.sub(r'([a-z0-9])([A-Z])', r'\1 \2', str(name))
    return [t.lower() for t in re.findall(r'[A-Za-z]+|\d+', spaced)]


def _label_looks_like_celltype(label: str) -> bool:
    toks = [t.lower() for t in re.findall(r'[A-Za-z0-9]+', str(label))]
    for t in toks:
        if t in _CELL_VOCAB or _CELL_SUFFIX_RE.search(t) or _MARKER_RE.match(t):
            return True
    return False


def score_celltype_column(name: str, values) -> tuple[float, list[str]]:
    """Score how likely an ``obs`` column holds cell-type labels (higher is better)."""
    toks = _name_tokens(name)
    joined = ''.join(toks)
    score, reasons = 0.0, []
    if any(k in joined for k in ('celltype', 'cellontology', 'cellclass', 'cellannot', 'celllabel', 'cellsubtype', 'cellident')):
        score += 4; reasons.append('name mentions cell type')
    elif set(toks) & _ANNOTATION_TOKENS:
        score += 3; reasons.append('name mentions annotation')
    elif set(toks) & _TYPE_TOKENS:
        score += 2; reasons.append('name mentions type/lineage/class')
    elif set(toks) & _LABEL_TOKENS:
        score += 1.5; reasons.append('name mentions label')
    elif set(toks) & _CLUSTER_TOKENS:
        score += 1; reasons.append('name mentions clusters')
    neg = set(toks) & _NEGATIVE_NAME_TOKENS
    if neg:
        score -= 3; reasons.append(f"name looks like metadata ({', '.join(sorted(neg))})")

    labels = [str(v) for v in pd.unique(pd.Series(values).dropna())][:500]
    n = len(labels)
    if n:
        frac = sum(_label_looks_like_celltype(v) for v in labels) / n
        if frac >= 0.25:
            score += 4 * frac; reasons.append(f'{frac:.0%} of values look like cell types')
        if all(re.fullmatch(r'-?\d+(\.\d+)?', v) for v in labels):
            score -= 1.5; reasons.append('values are numeric cluster ids')
    if 3 <= n <= 150:
        score += 1
    elif n == 2:
        score -= 1; reasons.append('only two values')
    elif n > 150:
        score -= 1; reasons.append('more than 150 values')
    return round(score, 2), reasons


def _candidate_celltype_columns(obs: pd.DataFrame) -> list[dict]:
    """Rank categorical/string columns of ``obs`` by how likely they hold cell types."""
    n = len(obs)
    cands = []
    for col in obs.columns:
        s = obs[col]
        if not (isinstance(s.dtype, pd.CategoricalDtype) or s.dtype == object or pd.api.types.is_string_dtype(s)):
            continue
        nunique = int(s.nunique(dropna=True))
        if nunique < 2 or nunique > min(500, max(2, n // 2)):
            continue
        score, reasons = score_celltype_column(col, s)
        cands.append({'column': col, 'n_unique': nunique, 'score': score, 'reasons': reasons})
    cands.sort(key=lambda d: (-d['score'], d['n_unique'], d['column']))
    return cands


def suggest_celltype_key(obs: pd.DataFrame) -> tuple[Optional[str], str, list[dict]]:
    """Return ``(suggested column or None, confidence, ranked candidates)``."""
    cands = _candidate_celltype_columns(obs)
    if not cands or cands[0]['score'] <= 0:
        return None, 'none', cands
    top = cands[0]['score']
    second = cands[1]['score'] if len(cands) > 1 else float('-inf')
    conf = 'high' if (top >= 3 and top - second >= 1.5) else 'low'
    return cands[0]['column'], conf, cands


def load_anndata(data: str | os.PathLike | anndata.AnnData) -> anndata.AnnData:
    if isinstance(data, anndata.AnnData):
        return data
    path = Path(data).expanduser()
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    return sc.read_h5ad(path)


def inspect_data(data, celltype_key: Optional[str] = None) -> dict:
    """Summarise an input dataset for humans, agents and the UI.

    Returns a JSON-serialisable dict with cell/gene counts, the detected matrix
    type (with evidence and confidence), whether raw counts are recoverable,
    ranked candidate cell-type columns, and a suggested ``celltype_key``.
    """
    adata = load_anndata(data)
    info = detect_matrix_type(adata.X)
    raw = find_raw_counts(adata)
    suggested, suggested_conf, cands = suggest_celltype_key(adata.obs)
    if celltype_key:
        suggested, suggested_conf = celltype_key, 'user'
    var = adata.var_names[:200]
    looks_ensembl = float(np.mean([str(v).startswith('ENS') for v in var])) if len(var) else 0.0

    if info.kind in ('scaled', 'empty'):
        standard = f'unusable ({info.reason})'
    elif info.confidence == 'low':
        standard = f'needs_attention (detected {info.kind} with low confidence; pass input_type to override)'
    else:
        standard = 'ok'
    if raw is not None:
        reference = 'ok'
    elif info.kind in ('scaled', 'empty'):
        reference = f'unusable ({info.reason})'
    else:
        reference = f'ok (adata.X is {info.kind}; the reference workflow is rank-based so this is fine)'

    out = {
        'path': str(data) if not isinstance(data, anndata.AnnData) else None,
        'n_cells': int(adata.n_obs),
        'n_genes': int(adata.n_vars),
        'matrix': info.to_dict(),
        'raw_counts_available': raw is not None,
        'raw_counts_source': raw[1].source if raw else None,
        'gene_id_type': 'ensembl' if looks_ensembl > 0.5 else 'symbol',
        'obs_columns': [
            {'column': c, 'dtype': str(adata.obs[c].dtype), 'n_unique': int(adata.obs[c].nunique(dropna=True))}
            for c in adata.obs.columns
        ],
        'celltype_candidates': cands,
        'suggested_celltype_key': suggested,
        'suggested_celltype_key_confidence': suggested_conf,
        'recommended_workflow_input': {'standard': standard, 'reference': reference},
    }
    if suggested and suggested in adata.obs.columns:
        counts = adata.obs[suggested].astype(str).value_counts()
        out['celltype_counts'] = {str(k): int(v) for k, v in counts.items()}
    return out


# ── preparation for each workflow ────────────────────────────────────────────

def _normalise_input_type(input_type: Optional[str]) -> Optional[str]:
    if input_type is None:
        return None
    key = str(input_type).strip().lower()
    if key not in INPUT_TYPE_ALIASES:
        raise ValueError(f"input_type must be one of raw, normalized, log1p (got '{input_type}')")
    return INPUT_TYPE_ALIASES[key]


def prepare_for_standard(adata: anndata.AnnData, copy: bool = True, input_type: Optional[str] = None,
                         ) -> tuple[anndata.AnnData, MatrixInfo, list[str]]:
    """Return log1p-normalised data for the standard workflow.

    Raw counts (integer or fractional) get ``normalize_total(1e4)`` + ``log1p``;
    normalised-but-not-logged data gets ``log1p``; log1p data is used as is.
    ``input_type`` ('raw' | 'normalized' | 'log1p') overrides the detection.
    """
    info = detect_matrix_type(adata.X)
    override = _normalise_input_type(input_type)
    if override is not None:
        if info.kind == 'scaled':
            raise ValueError("adata.X contains negative values; input_type cannot override that.")
        if override != info.kind:
            logger.info(f"input_type='{input_type}' overrides the detected matrix type ({info.describe()}).")
        info.kind, info.confidence, info.user_override = override, 'high', True
        info.reason = f"set by user (input_type='{input_type}')"
    elif info.confidence == 'low':
        logger.warning(
            f"Matrix type detected with low confidence: {info.describe()}. Treating adata.X as {info.kind}; "
            "pass input_type='raw' | 'normalized' | 'log1p' (CLI: --input-type) if this is wrong."
        )
    elif info.confidence == 'medium':
        logger.info(f"Matrix type detected as {info.kind} ({info.reason}); pass input_type to override.")

    applied: list[str] = []
    if info.kind == 'log1p':
        return adata, info, applied
    if info.kind == 'raw_counts':
        what = 'raw counts' if info.is_integer else 'non-integer counts'
        logger.info(f"Input looks like {what}; applying normalize_total(target_sum=1e4) + log1p for the standard workflow.")
        out = adata.copy() if copy else adata
        sc.pp.normalize_total(out, target_sum=1e4)
        sc.pp.log1p(out)
        applied = ['normalize_total(target_sum=1e4)', 'log1p']
        return out, info, applied
    if info.kind == 'normalized':
        logger.info("Input looks normalised but not log-transformed (max value %.1f); applying log1p." % info.max_value)
        out = adata.copy() if copy else adata
        sc.pp.log1p(out)
        return out, info, ['log1p']
    if info.kind == 'scaled':
        raise ValueError(
            "adata.X contains negative values (scaled/z-scored data). FastCCC needs "
            "normalised log1p expression; pass the unscaled matrix (e.g. adata.raw.to_adata())."
        )
    raise ValueError(f"Cannot use matrix of kind '{info.kind}' for the standard workflow ({info.reason}).")


def prepare_for_reference(adata: anndata.AnnData, copy: bool = True) -> tuple[anndata.AnnData, MatrixInfo, list[str]]:
    """Return a CSR AnnData for reference building / inference.

    The reference workflows rank genes within each cell, so any per-cell
    monotone transform of the counts (``normalize_total``, ``log1p``, CPM,
    sqrt) gives identical results.  Raw counts are preferred when present in
    ``X``, ``layers`` or ``.raw`` (they usually carry the full gene set);
    otherwise ``adata.X`` is used as is.  Only matrices with negative values
    (per-gene scaled / centred data) are rejected.
    """
    applied: list[str] = []
    found = find_raw_counts(adata)
    if found is not None:
        out, info = found
        if info.source != 'X':
            logger.info(f"Using raw counts from {info.source} for the reference workflow.")
            applied.append(f'raw counts taken from {info.source}')
        elif copy:
            out = out.copy()   # rank_preprocess digitises X in place; never touch the caller's object
    else:
        info = detect_matrix_type(adata.X)
        if info.kind == 'scaled':
            raise ValueError(
                "adata.X contains negative values (scaled / centred data), which cannot be ranked within cells. "
                "Pass counts or log-normalised expression instead (e.g. adata.raw.to_adata() or a counts layer)."
            )
        if info.kind == 'empty':
            raise ValueError("adata.X has no non-zero values.")
        logger.info(
            f"No raw count matrix found; using adata.X ({info.kind}). The reference workflow ranks genes within "
            "each cell, so per-cell normalisation and log1p do not change the results."
        )
        out = adata.copy() if copy else adata
        applied.append(f'used adata.X ({info.kind}) directly; the reference workflow is rank-based')
    if not sp.isspmatrix_csr(out.X):
        out.X = sp.csr_matrix(out.X)
        applied.append('converted X to CSR')
    return out, info, applied
