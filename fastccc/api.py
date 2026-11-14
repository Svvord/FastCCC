"""One-call entry points: ``run``, ``build_reference_panel``, ``compare_to_reference``.

These wrap the original workflow functions (which are unchanged) and add:
input inspection and preparation, database/panel resolution, a run manifest,
optional report generation and structured result objects.
"""
from __future__ import annotations

import json
import os
import re
import time
from pathlib import Path
from typing import Optional, Sequence

import anndata
import numpy as np
import pandas as pd
from loguru import logger

from .datacheck import load_anndata, prepare_for_standard, prepare_for_reference, _candidate_celltype_columns
from .celltype_match import (MappingSuggestion, suggest_mapping_from_labels, compute_coverage,
                             validate_celltype_mapping, load_mapping)
from .resources import resolve_database, resolve_reference, database_name
from .results import (
    update_manifest,
    RunResult, MultiRunResult, ReferenceResult, ReferencePanel,
    write_manifest, find_task_ids,
)

DEFAULT_RESULTS_DIR = 'fastccc_results'


# ── helpers ──────────────────────────────────────────────────────────────────

def _slug(text: str) -> str:
    text = re.sub(r'[^\w.\-]+', '_', str(text).strip())
    return text.strip('_') or 'group'


def _dataset_name(data, name: Optional[str], fallback: str) -> str:
    if name:
        return name
    if not isinstance(data, anndata.AnnData):
        return Path(str(data)).stem
    return fallback


def _resolve_labels(adata: anndata.AnnData, celltype_key: Optional[str], celltype_file: Optional[str]) -> pd.Series:
    if celltype_key is not None:
        if celltype_key not in adata.obs.columns:
            cands = [c['column'] for c in _candidate_celltype_columns(adata.obs)]
            raise KeyError(
                f"celltype_key '{celltype_key}' is not a column of adata.obs. "
                f"Candidate columns: {cands or list(adata.obs.columns)[:10]}"
            )
        return adata.obs[celltype_key].astype(str)
    if celltype_file is not None:
        labels = pd.read_csv(celltype_file, sep='\t', index_col=0)
        if labels.shape[1] < 1:
            raise ValueError(f"{celltype_file} must have a cell-id column and a cell-type column.")
        return labels.iloc[:, 0].astype(str)
    cands = [c['column'] for c in _candidate_celltype_columns(adata.obs)]
    raise ValueError(
        "Pass celltype_key (a column of adata.obs) or celltype_file (a TSV). "
        f"Likely cell-type columns: {cands[:5] if cands else 'none detected'}"
    )


def _input_summary(data, adata: anndata.AnnData, labels: pd.Series, celltype_key, celltype_file, matrix_info, applied) -> dict:
    counts = labels.value_counts()
    return {
        'source': None if isinstance(data, anndata.AnnData) else str(Path(str(data)).resolve()),
        'n_cells': int(adata.n_obs), 'n_genes': int(adata.n_vars),
        'celltype_key': celltype_key, 'celltype_file': str(celltype_file) if celltype_file else None,
        'n_celltypes': int(counts.shape[0]),
        'celltype_counts': {str(k): int(v) for k, v in counts.items()},
        'matrix_detected': matrix_info.to_dict(), 'preprocessing_applied': list(applied),
    }


def _new_task_id(result_dir: Path, before: set[str]) -> str:
    after = find_task_ids(result_dir)
    new = [t for t in after if t not in before]
    if len(new) == 1:
        return new[0]
    if not new:
        raise RuntimeError(f"FastCCC finished but no new *_significant_results.tsv appeared in {result_dir}")
    return new[-1]


# ── database fit ─────────────────────────────────────────────────────────────

MIN_DATABASE_GENE_COVERAGE = 0.05


def _database_fit(db, adata, gene_col: str = 'hgnc_symbol') -> dict:
    """Species of the database and the fraction of its genes (``gene_col``: ``hgnc_symbol`` or ``ensembl``) present in ``adata.var_names``."""
    from .dbtools import database_species
    table = pd.read_csv(Path(db) / 'gene_table.csv')
    genes = table[gene_col].dropna().astype(str)
    genes = genes[genes.str.strip() != ''].drop_duplicates()
    found = int(genes.isin(set(map(str, adata.var_names))).sum())
    return {'species': database_species(db), 'gene_id': gene_col, 'gene_coverage': round(found / len(genes), 3) if len(genes) else 0.0,
            'n_database_genes': int(len(genes)), 'n_database_genes_found': found}


def _database_entry(db, adata, gene_col: str = 'hgnc_symbol') -> dict:
    return {'name': database_name(db), 'path': str(db), **_database_fit(db, adata, gene_col)}


def _check_database_fit(db, adata, gene_col: str = 'hgnc_symbol') -> dict:
    """Stop early when the data holds almost none of the database's genes (wrong species or gene identifiers)."""
    fit = _database_fit(db, adata, gene_col)
    name, cov = database_name(db), fit['gene_coverage']
    if cov < MIN_DATABASE_GENE_COVERAGE:
        looks_ensembl = float(np.mean([str(v).startswith('ENS') for v in adata.var_names[:2000]])) > 0.5
        hint = ("var_names look like Ensembl ids: pass convert_type='ensembl' (`--gene-id ensembl`) or rename them to symbols"
                if looks_ensembl and gene_col != 'ensembl' else
                "check the species and that var_names are gene symbols (or Ensembl ids with convert_type='ensembl')")
        raise ValueError(
            f"Only {fit['n_database_genes_found']} of the {fit['n_database_genes']} genes of database '{name}' ({fit['species']}) "
            f"are in the data ({cov:.1%}); {hint}. `fastccc db list --data <file>` shows the coverage of every database.")
    if cov < 0.5:
        logger.warning(f"Only {cov:.0%} of the genes of database '{name}' are in the data; interactions with missing genes are not tested.")
    return fit


# ── standard workflow ────────────────────────────────────────────────────────

def _run_single(
    adata: anndata.AnnData, *, data, name: str, result_dir: Path, db: Path,
    celltype_key, celltype_file, matrix_info, applied, params: dict,
) -> RunResult:
    from .core import Cauchy_combination_of_statistical_analysis_methods

    labels = _resolve_labels(adata, celltype_key, celltype_file)
    before = set(find_task_ids(result_dir))
    t0 = time.perf_counter()
    Cauchy_combination_of_statistical_analysis_methods(
        database_file_path=str(db), celltype_file_path=celltype_file, counts_file_path=adata,
        convert_type=params['convert_type'],
        single_unit_summary_list=list(params['single_unit_summary_list']),
        complex_aggregation_list=list(params['complex_aggregation_list']),
        LR_combination_list=list(params['LR_combination_list']),
        min_percentile=params['min_percentile'], save_path=str(result_dir), meta_key=celltype_key,
        select_list=list(params['select_list']), filter_=params['filter_'], use_DEG=params['use_DEG'],
        small_p_method=params.get('small_p_method', 'acat'),
    )
    elapsed = time.perf_counter() - t0
    task_id = _new_task_id(result_dir, before)
    res = RunResult(result_dir=result_dir.resolve(), task_id=task_id, name=name,
                    database=database_name(db), database_path=db)
    res.manifest = write_manifest(result_dir, {
        'workflow': 'run', 'name': name, 'task_id': task_id,
        'parameters': {k: (list(v) if isinstance(v, (tuple, list)) else v) for k, v in params.items()},
        'database': _database_entry(db, adata, params['convert_type']),
        'input': _input_summary(data, adata, labels, celltype_key, celltype_file, matrix_info, applied),
        'timing_seconds': {'analysis': round(elapsed, 1)},
        'outputs': {
            'significant_results': str(res.significant_path), 'cauchy_pvals': str(res.pvals_path),
            'average_interactions_strength': str(res.strength_path), 'report_html': None,
        },
    })
    logger.success(f"FastCCC run '{name}' done in {elapsed:.1f}s (task_id={task_id}) -> {result_dir}")
    return res


def _bootstrap_support(res: RunResult, adata, *, data, db, celltype_key, celltype_file, matrix_info, applied, params,
                       n: int, fraction: float, seed: int) -> Path:
    """Re-run the analysis on ``n`` stratified cell subsamples and record, for every significant
    interaction of the main run, the fraction of subsamples in which it is significant."""
    import shutil
    import glob as _glob
    labels = _resolve_labels(adata, celltype_key, celltype_file).values
    rng = np.random.default_rng(seed)
    main = res.significant
    keys = list(zip(main['sender_celltype'].astype(str), main['receiver_celltype'].astype(str), main['LRI_ID'].astype(str)))
    hits = {k: 0 for k in keys}
    extra: dict = {}
    bs_root = res.result_dir / 'bootstrap'
    t0 = time.perf_counter()
    for k in range(n):
        keep = []
        for ct in np.unique(labels):
            idx = np.where(labels == ct)[0]
            m = max(1, int(round(len(idx) * fraction)))
            keep.extend(rng.choice(idx, size=min(m, len(idx)), replace=False))
        keep = np.sort(np.array(keep))
        sub = adata[keep].copy()
        bdir = bs_root / f'b{k:02d}'
        logger.info(f"Bootstrap {k + 1}/{n}: {len(keep):,} cells ({fraction:.0%} of every cell type)")
        bres = _run_single(sub, data=data, name=f'{res.name} bootstrap {k + 1}', result_dir=bdir, db=db, celltype_key=celltype_key,
                           celltype_file=celltype_file, matrix_info=matrix_info, applied=applied, params=params)
        sig = bres.significant
        found = set(zip(sig['sender_celltype'].astype(str), sig['receiver_celltype'].astype(str), sig['LRI_ID'].astype(str)))
        for key in found:
            if key in hits:
                hits[key] += 1
            else:
                extra[key] = extra.get(key, 0) + 1
        # keep only the small significant table of each resample; the per-module tables are large
        for f in _glob.glob(str(bdir / '*.tsv')):
            if not f.endswith('_significant_results.tsv'):
                os.remove(f)
    elapsed = time.perf_counter() - t0
    support = pd.DataFrame({
        'sender_celltype': [k[0] for k in keys], 'receiver_celltype': [k[1] for k in keys], 'LRI_ID': [k[2] for k in keys],
        'n_significant': [hits[k] for k in keys], 'n_bootstrap': n, 'support': [hits[k] / n for k in keys],
    })
    path = res.result_dir / 'bootstrap_support.tsv'
    support.to_csv(path, sep='\t', index=False)
    n_stable = int((support['support'] >= 0.8).sum())
    logger.success(f"Bootstrap support written ({n} resamples, {elapsed:.0f}s): {n_stable}/{len(support)} interactions significant in >= 80% of resamples; "
                   f"{len(extra)} interactions appeared only in resamples.")
    res.manifest = update_manifest(res.result_dir, bootstrap={
        'n': n, 'fraction': fraction, 'seed': seed, 'support_file': str(path), 'n_stable_80': n_stable,
        'n_only_in_resamples': len(extra), 'timing_seconds': round(elapsed, 1),
    }) or res.manifest
    return path


def run(
    data,
    celltype_key: Optional[str] = None,
    *,
    celltype_file: Optional[str | os.PathLike] = None,
    database: str | os.PathLike = 'CPDBv5.0.0',
    out_dir: Optional[str | os.PathLike] = None,
    name: Optional[str] = None,
    split_by: Optional[str] = None,
    groups: Optional[Sequence[str]] = None,
    include_all: bool = True,
    report: bool = True,
    report_kwargs: Optional[dict] = None,
    convert_type: str = 'hgnc_symbol',
    single_unit_summary_list: Sequence[str] = ('Mean', 'Median', 'Q3', 'Quantile_0.9'),
    complex_aggregation_list: Sequence[str] = ('Minimum', 'Average'),
    LR_combination_list: Sequence[str] = ('Arithmetic', 'Geometric'),
    min_percentile: float = 0.1,
    use_DEG: bool = False,
    select_list: Sequence[str] = (),
    filter_: bool = False,
    input_type: Optional[str] = None,
    min_cells: int = 20,
    small_p_method: str = 'acat',
    legacy_report: bool = False,
    explorer_kwargs: Optional[dict] = None,
    bootstrap: int = 0,
    bootstrap_fraction: float = 0.8,
    bootstrap_seed: int = 0,
    copy: bool = True,
) -> RunResult | MultiRunResult:
    """Run the standard FastCCC analysis (Cauchy combination of all scoring modules).

    Parameters
    ----------
    data
        Path to an ``.h5ad`` file or an in-memory ``AnnData``. Raw counts
        (integer or fractional) are normalised (``normalize_total(1e4)`` +
        ``log1p``) automatically; log1p data is used as is. The detection and
        its evidence are recorded in the run manifest.
    celltype_key
        Column of ``adata.obs`` holding cell-type labels. Alternatively pass
        ``celltype_file`` (TSV with cell id and cell type columns).
    database
        Packaged LRI database name (see ``fastccc.list_databases()``) or a path.
    out_dir
        Output directory (default ``./fastccc_results/<name>``).
    split_by, groups, include_all
        Run every group of ``adata.obs[split_by]`` separately (plus the full
        dataset when ``include_all``). With exactly two groups and
        ``include_all=True`` the report gains a descriptive comparison tab.
    report
        Build the interactive result explorer (``explorer.html``) after the analysis.
    legacy_report
        Also render the older static HTML report (``fastccc.report.generate_report``).
    bootstrap
        Number of cell-resampling runs (``bootstrap_fraction`` of the cells of every
        cell type, without replacement) used to estimate, for every significant
        interaction, the fraction of resamples in which it stays significant
        (``bootstrap_support.tsv``). 10 is a reasonable default; 0 disables.
    input_type
        ``'raw'``, ``'normalized'`` or ``'log1p'``: state what ``adata.X``
        holds instead of relying on the automatic detection.
    small_p_method
        How extreme module p-values enter the Cauchy combination. ``'acat'``
        (default): p = 1 enters as ``1 - 1/16`` and p below 1e-15 as
        ``1 / (p * pi)``. ``'legacy'``: the previous rule (p at or below 1e-16
        treated as 0 through ``tan``).
    min_cells
        With ``split_by``: groups with fewer cells than this are skipped (with
        a warning) and listed under ``skipped_groups`` in the manifest. A group
        of a handful of cells cannot support the per-cell-type summaries and
        would only add noise to the per-group results and the replication
        fractions.
    copy
        Copy the input before normalising so the caller's object is untouched.

    Returns
    -------
    RunResult, or MultiRunResult when ``split_by`` is given.
    """
    params = {
        'convert_type': convert_type, 'single_unit_summary_list': single_unit_summary_list,
        'complex_aggregation_list': complex_aggregation_list, 'LR_combination_list': LR_combination_list,
        'min_percentile': min_percentile, 'use_DEG': use_DEG, 'select_list': select_list, 'filter_': filter_,
        'input_type': input_type, 'small_p_method': small_p_method,
    }
    if small_p_method not in ('acat', 'legacy'):
        raise ValueError(f"small_p_method must be 'acat' or 'legacy', got {small_p_method!r}")
    report_kwargs = dict(report_kwargs or {})
    db = resolve_database(str(database))
    adata = load_anndata(data)
    _check_database_fit(db, adata, convert_type)
    name = _dataset_name(data, name, 'sample')
    out = Path(out_dir).expanduser() if out_dir else Path(DEFAULT_RESULTS_DIR) / _slug(name)
    adata_p, info, applied = prepare_for_standard(adata, copy=copy, input_type=input_type)
    celltype_file = str(celltype_file) if celltype_file else None

    if split_by is None:
        res = _run_single(adata_p, data=data, name=name, result_dir=out, db=db, celltype_key=celltype_key,
                          celltype_file=celltype_file, matrix_info=info, applied=applied, params=params)
        if bootstrap and bootstrap > 0:
            _bootstrap_support(res, adata_p, data=data, db=db, celltype_key=celltype_key, celltype_file=celltype_file,
                               matrix_info=info, applied=applied, params=params, n=int(bootstrap),
                               fraction=bootstrap_fraction, seed=bootstrap_seed)
        if legacy_report:
            res.report(**report_kwargs)
        if report:
            res.explore(**(explorer_kwargs or {}))
        return res

    # ── split run ────────────────────────────────────────────────────────
    if split_by not in adata_p.obs.columns:
        raise KeyError(f"split_by '{split_by}' is not a column of adata.obs: {list(adata_p.obs.columns)[:15]}")
    labels = adata_p.obs[split_by].astype(str)
    available = sorted(labels.unique())
    groups = [str(g) for g in groups] if groups else available
    missing = [g for g in groups if g not in available]
    if missing:
        raise ValueError(f"Groups {missing} not found in adata.obs['{split_by}']; available: {available}")

    out.mkdir(parents=True, exist_ok=True)
    runs: dict[str, RunResult] = {}
    all_res: Optional[RunResult] = None
    if include_all:
        all_res = _run_single(adata_p, data=data, name=f'{name} (all cells)', result_dir=out / 'all', db=db,
                              celltype_key=celltype_key, celltype_file=celltype_file, matrix_info=info,
                              applied=applied, params=params)
    counts = labels.value_counts()
    skipped = {g: int(counts.get(g, 0)) for g in groups if int(counts.get(g, 0)) < int(min_cells)}
    if skipped:
        shown = ', '.join(f'{g} ({n})' for g, n in list(skipped.items())[:8]) + (' ...' if len(skipped) > 8 else '')
        logger.warning(f"Skipping {len(skipped)} group(s) of '{split_by}' with fewer than {min_cells} cells: {shown} "
                       f"(set min_cells lower to include them).")
    groups = [g for g in groups if g not in skipped]
    if not groups:
        raise ValueError(f"Every group of '{split_by}' has fewer than {min_cells} cells; lower min_cells.")
    for g in groups:
        sub = adata_p[labels == g].copy()
        runs[g] = _run_single(sub, data=data, name=g, result_dir=out / _slug(g), db=db, celltype_key=celltype_key,
                              celltype_file=celltype_file, matrix_info=info, applied=applied, params=params)

    multi = MultiRunResult(result_dir=out.resolve(), split_by=split_by, runs=runs, all=all_res)
    multi.manifest = write_manifest(out, {
        'workflow': 'run_split', 'name': name, 'split_by': split_by,
        'groups': {g: str(r.result_dir) for g, r in runs.items()},
        'min_cells': int(min_cells), 'skipped_groups': skipped,
        'all': str(all_res.result_dir) if all_res else None,
        'database': _database_entry(db, adata, convert_type),
        'parameters': {k: (list(v) if isinstance(v, (tuple, list)) else v) for k, v in params.items()},
        'outputs': {'report_html': None},
    })
    if legacy_report:
        report_split(multi, **report_kwargs)
    if report:
        multi.explore(**(explorer_kwargs or {}))
    return multi


def report_split(multi: MultiRunResult, output_dir: Optional[str | os.PathLike] = None, **kwargs) -> Path:
    """Report for a split run: one 4-tab report when there are two groups plus 'all',
    otherwise one report per run."""
    from .report import generate_report
    from .results import update_manifest

    out = Path(output_dir) if output_dir else multi.result_dir / 'report'
    groups = multi.groups
    if multi.all is not None and len(groups) == 2:
        a, b = groups
        db = str(multi.all.database_path or resolve_database(multi.all.database))
        kwargs.setdefault('sample_name', multi.all.name)
        path = Path(generate_report(
            str(multi.all.result_dir), multi.all.task_id, db, output_dir=str(out),
            cond_a_result_dir=str(multi.runs[a].result_dir), cond_a_task_id=multi.runs[a].task_id, cond_a_name=a,
            cond_b_result_dir=str(multi.runs[b].result_dir), cond_b_task_id=multi.runs[b].task_id, cond_b_name=b,
            **kwargs,
        ))
        multi.report_path = path
        update_manifest(multi.result_dir, outputs={'report_html': str(path)})
        return path
    if len(groups) != 2:
        logger.warning("The condition-comparison tab needs exactly two groups; generating one report per run instead.")
    elif multi.all is None:
        logger.warning("The condition-comparison tab needs include_all=True; generating one report per run instead.")
    paths = []
    for r in ([multi.all] if multi.all else []) + list(multi.runs.values()):
        paths.append(r.report(**kwargs))
    multi.report_path = paths[0] if paths else None
    update_manifest(multi.result_dir, outputs={'report_html': str(multi.report_path) if multi.report_path else None,
                                               'report_html_per_run': [str(p) for p in paths]})
    return multi.report_path


# ── reference panel construction ─────────────────────────────────────────────

def build_reference_panel(
    data,
    celltype_key: Optional[str] = None,
    *,
    name: str,
    out_dir: str | os.PathLike = 'fastccc_reference',
    celltype_file: Optional[str | os.PathLike] = None,
    database: str | os.PathLike = 'CPDBv5.0.0',
    min_percentile: float = 0.1,
    subset: Optional[tuple[str, str] | str] = None,
    copy: bool = True,
) -> ReferencePanel:
    """Build a custom FastCCC reference panel from a dataset.

    Raw counts are preferred when present (``X``, ``layers['counts']`` or
    ``.raw``), but log-normalised data works too: the reference workflow ranks
    genes within each cell, so per-cell normalisation and log1p do not change
    the panel. The panel is written to ``<out_dir>/<name>/`` and can then be
    passed to ``compare_to_reference(reference=<that path>)``.

    ``subset=('condition', 'control')`` (or the string ``'condition:control'``)
    builds the panel from those cells only.
    """
    from .build_reference import build_reference_workflow

    db = resolve_database(str(database))
    adata = load_anndata(data)
    _check_database_fit(db, adata)
    if subset is not None:
        col, value = _parse_subset(subset)
        if col not in adata.obs.columns:
            raise KeyError(f"subset column '{col}' is not a column of adata.obs: {list(adata.obs.columns)[:15]}")
        mask = adata.obs[col].astype(str) == str(value)
        if not mask.any():
            raise ValueError(f"'{value}' not found in adata.obs['{col}']; values: {sorted(set(adata.obs[col].astype(str)))[:20]}")
        logger.info(f"Building the panel from {int(mask.sum()):,} cells with {col} == '{value}'.")
        adata = adata[mask].copy()
        copy = False
    celltype_file = str(celltype_file) if celltype_file else None
    labels = _resolve_labels(adata, celltype_key, celltype_file)
    adata_p, info, applied = prepare_for_reference(adata, copy=copy)
    out = Path(out_dir).expanduser()
    t0 = time.perf_counter()
    build_reference_workflow(
        database_file_path=str(db), reference_counts_file_path=adata_p, celltype_file_path=celltype_file,
        reference_name=name, save_path=str(out), meta_key=celltype_key, min_percentile=min_percentile,
    )
    elapsed = time.perf_counter() - t0
    panel_dir = out / name
    write_manifest(panel_dir, {
        'workflow': 'build_reference', 'name': name,
        'parameters': {'min_percentile': min_percentile, 'subset': list(_parse_subset(subset)) if subset else None},
        'database': _database_entry(db, adata),
        'input': _input_summary(data, adata, labels, celltype_key, celltype_file, info, applied),
        'timing_seconds': {'build': round(elapsed, 1)},
    })
    return ReferencePanel.load(panel_dir, source='custom')


def _parse_subset(subset) -> tuple[str, str]:
    if isinstance(subset, str):
        if ':' not in subset:
            raise ValueError("subset must be 'column:value' or a (column, value) tuple")
        col, value = subset.split(':', 1)
        return col, value
    if len(subset) != 2:
        raise ValueError("subset must be 'column:value' or a (column, value) tuple")
    return str(subset[0]), str(subset[1])


def _label_counts(query, celltype_key, celltype_file) -> dict[str, int]:
    """``{label: n_cells}`` for an AnnData/path, a sequence of labels, or a counts dict."""
    if isinstance(query, dict):
        return {str(k): int(v) for k, v in query.items()}
    if isinstance(query, (str, os.PathLike, anndata.AnnData)) and not isinstance(query, (list, tuple)):
        try:
            adata = load_anndata(query)
        except (FileNotFoundError, OSError):
            adata = None
        if adata is not None:
            labels = _resolve_labels(adata, celltype_key, celltype_file)
            return {str(k): int(v) for k, v in labels.value_counts().items()}
    return {str(k): int(v) for k, v in pd.Series(list(query)).astype(str).value_counts().items()}


def suggest_celltype_mapping(
    query,
    celltype_key: Optional[str] = None,
    *,
    reference,
    celltype_file: Optional[str | os.PathLike] = None,
    reference_root: Optional[str | os.PathLike] = None,
    download: bool = True,
    min_score: float = 0.5,
) -> MappingSuggestion:
    """Suggest a ``{reference cell type: query cell type}`` mapping for ``compare_to_reference``.

    Same argument order as :func:`compare_to_reference`:
    ``suggest_celltype_mapping(data, celltype_key, reference='liver')``.

    Parameters
    ----------
    query
        ``.h5ad`` path / ``AnnData`` (with ``celltype_key``), a list of labels, or a
        ``{label: n_cells}`` dict.
    reference
        Tissue name (``'blood'``), panel directory, or a :class:`ReferencePanel`.
    min_score
        Panel cell types whose best similarity is below this are left unmapped
        (reported in ``unmatched_reference``).

    Returns
    -------
    MappingSuggestion with ``.mapping`` (pass it as ``celltype_mapping``),
    ``.table()``, ``.ambiguous``, ``.unmatched_query`` and ``.coverage``.
    The matching is deterministic (name normalisation, synonyms, weighted
    token overlap); review it before use.
    """
    if isinstance(reference, ReferencePanel):
        panel = reference
    else:
        info = resolve_reference(reference, reference_root=reference_root, download=download)
        panel = ReferencePanel.load(info['path'], source=info['source'], tissue=info['tissue'])
    counts = _label_counts(query, celltype_key, str(celltype_file) if celltype_file else None)
    sug = suggest_mapping_from_labels(panel.celltypes, counts, min_score=min_score, reference_name=panel.name)
    logger.info(f"Cell-type mapping suggestion for panel '{panel.name}': {sug.summary().splitlines()[0]}")
    return sug


# ── reference-based comparison ───────────────────────────────────────────────

def compare_to_reference(
    data,
    celltype_key: Optional[str] = None,
    *,
    reference,
    celltype_file: Optional[str | os.PathLike] = None,
    celltype_mapping: Optional[dict | str | os.PathLike] = None,
    query_relabel: Optional[dict] = None,
    database: str | os.PathLike = 'CPDBv5.0.0',
    out_dir: Optional[str | os.PathLike] = None,
    name: Optional[str] = None,
    reference_name: Optional[str] = None,
    reference_root: Optional[str | os.PathLike] = None,
    report: bool = True,
    report_kwargs: Optional[dict] = None,
    legacy_report: bool = False,
    explorer_kwargs: Optional[dict] = None,
    copy: bool = True,
    download: bool = True,
) -> ReferenceResult:
    """Compare a query dataset against a reference panel (FastCCC's reference-based analysis).

    Parameters
    ----------
    data
        ``.h5ad`` path or ``AnnData``. Raw counts are preferred when present
        (``X``, ``layers['counts']`` or ``.raw``); otherwise ``adata.X`` is used
        as is, which gives identical results because the reference workflow
        ranks genes within each cell. Only scaled data (negative values) is
        rejected.
    reference
        One of: a healthy tissue name such as ``'liver'`` (downloaded on first
        use), a path to a custom panel directory, or a tuple
        ``(obs_column, control_value)`` / ``(obs_column, control_value, case_value)``
        meaning "build a panel from the control cells of this dataset and query
        the case cells against it".
    celltype_mapping
        ``{reference cell type: query cell type}`` dict or JSON path, used when
        the query labels differ from the panel's labels; ``'auto'`` uses
        :func:`suggest_celltype_mapping` and writes the mapping it used to
        ``celltype_mapping.json`` in the output directory. Without a mapping,
        only identically named cell types are compared. The mapping is
        validated and the coverage of the query cell types is checked before
        anything runs: no overlap at all raises an error.
    query_relabel
        ``{query cell type: new label}`` applied to the query labels first,
        e.g. to merge sub-types that the panel does not distinguish.

    Returns
    -------
    ReferenceResult with ``results`` (``trend_vs_ref`` = Up / Down / Both Sig / Both NS).
    """
    from .infer_query import infer_query_workflow

    report_kwargs = dict(report_kwargs or {})
    db = resolve_database(str(database))
    adata = load_anndata(data)
    celltype_file = str(celltype_file) if celltype_file else None
    built_from = None

    if isinstance(reference, (tuple, list)):
        if len(reference) not in (2, 3):
            raise ValueError("reference tuple must be (obs_column, control_value) or (obs_column, control_value, case_value)")
        col, control = str(reference[0]), str(reference[1])
        case = str(reference[2]) if len(reference) == 3 else None
        if col not in adata.obs.columns:
            raise KeyError(f"'{col}' is not a column of adata.obs")
        grp = adata.obs[col].astype(str)
        if control not in set(grp):
            raise ValueError(f"'{control}' not found in adata.obs['{col}']; values: {sorted(set(grp))[:20]}")
        if case is not None and case not in set(grp):
            raise ValueError(f"'{case}' not found in adata.obs['{col}']; values: {sorted(set(grp))[:20]}")
        query_name = name or (case if case else f'non-{control}')
        out = Path(out_dir).expanduser() if out_dir else Path(DEFAULT_RESULTS_DIR) / f'{_slug(query_name)}_vs_{_slug(control)}'
        control_adata = adata[grp == control].copy()
        query = adata[grp == case].copy() if case else adata[grp != control].copy()
        logger.info(f"Building a reference panel from {control_adata.n_obs:,} '{control}' cells; querying {query.n_obs:,} cells.")
        panel = build_reference_panel(control_adata, celltype_key, name=f'{_slug(control)}_reference',
                                      out_dir=out / 'reference', celltype_file=celltype_file, database=db, copy=False)
        built_from = {'obs_column': col, 'control_value': control, 'case_value': case}
        display_ref = reference_name or f'{control} (reference)'
    else:
        info = resolve_reference(reference, reference_root=reference_root, download=download)
        panel = ReferencePanel.load(info['path'], source=info['source'], tissue=info['tissue'])
        query = adata
        query_name = _dataset_name(data, name, 'Query')
        out = Path(out_dir).expanduser() if out_dir else Path(DEFAULT_RESULTS_DIR) / f'{_slug(query_name)}_vs_{_slug(panel.name)}'
        display_ref = reference_name or (f"Healthy {panel.name.replace('_', ' ').title()}" if panel.source == 'healthy_tissue' else panel.name)

    if panel.database and panel.database != database_name(db):
        raise ValueError(
            f"Reference panel '{panel.name}' was built with {panel.database} but database='{database_name(db)}'. "
            "Use the same LRI database for the panel and the query."
        )

    labels = _resolve_labels(query, celltype_key, celltype_file)
    query_p, info_m, applied = prepare_for_reference(query, copy=copy)
    if query_relabel:
        if celltype_key is None:
            raise ValueError("query_relabel needs celltype_key (labels from celltype_file cannot be relabelled)")
        query_p.obs = query_p.obs.copy()
        query_p.obs[celltype_key] = labels.map(lambda v: str(query_relabel.get(v, v))).values
        labels = query_p.obs[celltype_key].astype(str)
        logger.info(f"Relabelled query cell types: {query_relabel}")
    query_counts = {str(k): int(v) for k, v in labels.value_counts().items()}

    mapping_source = 'none'
    suggestion = None
    if isinstance(celltype_mapping, str) and celltype_mapping.lower() == 'auto':
        suggestion = suggest_mapping_from_labels(panel.celltypes, query_counts, reference_name=panel.name)
        mapping = suggestion.mapping
        mapping_source = 'auto'
        for line in suggestion.summary().splitlines():
            logger.info(line)
    elif celltype_mapping is not None:
        mapping = validate_celltype_mapping(load_mapping(celltype_mapping), panel.celltypes, query_counts)
        mapping_source = 'user'
    else:
        mapping = {}
    coverage = compute_coverage(panel.celltypes, query_counts, mapping)
    if coverage['n_covered'] == 0:
        shown = ', '.join(list(panel.celltypes)[:15]) + (' ...' if len(panel.celltypes) > 15 else '')
        raise ValueError(
            f"None of the {coverage['n_query_types']} query cell types matches a cell type of reference panel "
            f"'{panel.name}' ({panel.n_celltypes} types, e.g. {shown}). Pass a celltype_mapping "
            f"{{reference cell type: query cell type}}, use celltype_mapping='auto' "
            f"(CLI: --auto-map), or get a reviewed suggestion with "
            f"`fastccc map <data> -c {celltype_key or '<column>'} -r {panel.tissue or panel.path}`."
        )
    if coverage['uncovered_query_types']:
        logger.warning(
            f"{len(coverage['uncovered_query_types'])} of {coverage['n_query_types']} query cell types "
            f"({coverage['uncovered_cells_fraction']:.1%} of cells) have no counterpart in reference '{panel.name}' and "
            f"will be reported as 'Not in reference': {coverage['uncovered_query_types']}"
        )
    else:
        logger.info(f"All {coverage['n_query_types']} query cell types are covered by reference '{panel.name}'.")
    _check_database_fit(db, adata)
    mapping_arg = mapping or None
    out.mkdir(parents=True, exist_ok=True)
    if mapping:
        with open(out / 'celltype_mapping.json', 'w', encoding='utf-8') as fh:
            json.dump(mapping, fh, indent=2, ensure_ascii=False)

    t0 = time.perf_counter()
    infer_query_workflow(
        database_file_path=str(db), reference_path=str(panel.path), query_counts_file_path=query_p,
        celltype_file_path=celltype_file, save_path=str(out), celltype_mapping_dict=mapping_arg, meta_key=celltype_key,
    )
    elapsed = time.perf_counter() - t0

    res = ReferenceResult(result_dir=out.resolve(), query_name=query_name, reference=panel,
                          database=database_name(db), database_path=db)
    res.manifest = write_manifest(out, {
        'workflow': 'compare_to_reference', 'name': query_name,
        'reference': {**panel.to_dict(), 'built_from': built_from}, 'reference_display_name': display_ref,
        'database': _database_entry(db, adata),
        'input': _input_summary(data, query, labels, celltype_key, celltype_file, info_m, applied),
        'celltype_mapping': mapping, 'celltype_mapping_source': mapping_source,
        'query_relabel': query_relabel, 'celltype_coverage': coverage,
        'celltype_mapping_suggestion': suggestion.to_dict() if suggestion is not None else None,
        'timing_seconds': {'inference': round(elapsed, 1)},
        'outputs': {'results': str(res.results_path), 'interactions_strength': str(res.strength_path),
                    'percents_analysis': str(res.result_dir / 'query_percents_analysis.tsv'), 'report_html': None},
    })
    logger.success(f"Reference comparison '{query_name}' vs '{panel.name}' done in {elapsed:.1f}s -> {out}")
    if legacy_report:
        res.report(**report_kwargs)
    if report:
        res.explore(**(explorer_kwargs or {}))
    return res
