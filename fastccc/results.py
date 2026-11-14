"""Result objects and the per-run manifest (``fastccc_run.json``).

Every workflow writes a manifest next to its output files recording the
FastCCC version, all parameters, an input summary and the output paths.  The
result classes are thin, lazily-loading views over a result directory and can
be rebuilt from the manifest alone with ``.load(result_dir)``.
"""
from __future__ import annotations

import datetime as _dt
import glob
import json
import os
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
from loguru import logger

MANIFEST_NAME = 'fastccc_run.json'
SCHEMA_VERSION = 1


# ── manifest helpers ─────────────────────────────────────────────────────────

def _json_default(obj: Any):
    if isinstance(obj, (Path, os.PathLike)):
        return str(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (set, tuple)):
        return list(obj)
    if isinstance(obj, (_dt.datetime, _dt.date)):
        return obj.isoformat()
    if hasattr(obj, 'to_dict'):
        return obj.to_dict()
    return str(obj)


def write_manifest(result_dir: str | os.PathLike, payload: dict) -> dict:
    """Write ``fastccc_run.json`` into ``result_dir`` and return the completed payload."""
    path = Path(result_dir) / MANIFEST_NAME
    payload = dict(payload)
    payload.setdefault('schema_version', SCHEMA_VERSION)
    from . import __version__
    payload.setdefault('fastccc_version', __version__)
    payload.setdefault('created', _dt.datetime.now().isoformat(timespec='seconds'))
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, indent=2, default=_json_default)
    return json.loads(json.dumps(payload, default=_json_default))


def read_manifest(result_dir: str | os.PathLike) -> Optional[dict]:
    path = Path(result_dir) / MANIFEST_NAME
    if not path.is_file():
        return None
    with open(path, encoding='utf-8') as fh:
        return json.load(fh)


def update_manifest(result_dir: str | os.PathLike, **updates) -> Optional[dict]:
    manifest = read_manifest(result_dir)
    if manifest is None:
        return None
    manifest.update(updates)
    write_manifest(result_dir, manifest)
    return manifest


def find_task_ids(result_dir: str | os.PathLike) -> list[str]:
    hits = sorted(glob.glob(os.path.join(str(result_dir), '*_significant_results.tsv')), key=os.path.getmtime)
    return [os.path.basename(h).replace('_significant_results.tsv', '') for h in hits]


def find_task_id(result_dir: str | os.PathLike) -> str:
    ids = find_task_ids(result_dir)
    if not ids:
        raise FileNotFoundError(f"No *_significant_results.tsv found in {result_dir}")
    if len(ids) > 1:
        logger.warning(f"{len(ids)} FastCCC runs found in {result_dir}; using the newest ({ids[-1]}).")
    return ids[-1]


def _top_counts(series: pd.Series, n: int = 5) -> list[dict]:
    vc = series.value_counts().head(n)
    return [{'celltype': str(k), 'n': int(v)} for k, v in vc.items()]


# ── reference panel ──────────────────────────────────────────────────────────


def _explorer_mixin_explore(self, output=None, **kwargs):
    """Build the interactive result explorer (``explorer.html``) for this result."""
    from .explorer import build_explorer
    path = build_explorer(self.result_dir, output=output, **kwargs)
    self.manifest = read_manifest(self.result_dir) or self.manifest
    self._explorer_path = Path(path)
    return path


def _relocate(stored, *candidates) -> Optional[Path]:
    """Resolve a path recorded in a manifest.

    Manifests store absolute paths, which go stale when a result folder is moved,
    renamed or copied to another machine. The stored path is used when it still
    exists; otherwise the candidate locations (derived from the current result
    folder) are tried in order, and the stored path is returned unchanged when
    none exists so that callers can report it.
    """
    if stored:
        p = Path(stored).expanduser()
        if p.exists():
            return p
    for cand in candidates:
        if cand is not None and Path(cand).exists():
            return Path(cand)
    return Path(stored) if stored else None


def _explorer_path_from_manifest(manifest, result_dir=None) -> Optional[Path]:
    p = ((manifest or {}).get('outputs') or {}).get('explorer_html')
    if not p:
        return None
    return _relocate(p, Path(result_dir) / Path(p).name if result_dir else None)

@dataclass
class ReferencePanel:
    """A built FastCCC reference panel (shipped healthy tissue or custom)."""
    name: str
    path: Path
    source: str                      # healthy_tissue | custom
    database: str = ''
    min_percentile: float = 0.1
    celltypes: dict = field(default_factory=dict)
    tissue: Optional[str] = None

    @classmethod
    def load(cls, path: str | os.PathLike, source: str = 'custom', tissue: Optional[str] = None) -> 'ReferencePanel':
        path = Path(path).expanduser().resolve()
        cfg_path = path / 'config.toml'
        if not cfg_path.is_file():
            raise FileNotFoundError(f"Not a FastCCC reference panel (config.toml missing): {path}")
        with open(cfg_path, 'rb') as fh:
            cfg = tomllib.load(fh)
        return cls(
            name=str(cfg.get('reference_name', path.name)), path=path, source=source,
            database=str(cfg.get('LRI_database', '')), min_percentile=float(cfg.get('min_percentile', 0.1)),
            celltypes={str(k): int(v) for k, v in cfg.get('celltype', {}).items()}, tissue=tissue,
        )

    @property
    def n_cells(self) -> int:
        return int(sum(self.celltypes.values()))

    @property
    def n_celltypes(self) -> int:
        return len(self.celltypes)

    @property
    def is_activated(self) -> bool:
        return (self.path / 'ref_gene_pmf_dict.pkl').is_file()

    @property
    def celltype_counts(self) -> dict:
        """``{cell type: n_cells}`` of the panel (same as ``celltypes``)."""
        return dict(self.celltypes)

    def celltype_table(self) -> pd.DataFrame:
        return pd.DataFrame({'celltype': list(self.celltypes), 'n_cells': list(self.celltypes.values())}
                            ).sort_values('n_cells', ascending=False).reset_index(drop=True)

    def describe(self) -> str:
        lines = [f"Reference panel '{self.name}' ({self.source}{', ' + self.tissue if self.tissue else ''}): "
                 f"{self.n_celltypes} cell types, {self.n_cells:,} cells, database {self.database}, "
                 f"min_percentile {self.min_percentile}, {'activated' if self.is_activated else 'not yet activated'}",
                 f"path: {self.path}"]
        for ct, n in sorted(self.celltypes.items(), key=lambda kv: -kv[1]):
            lines.append(f"  {n:>9,d}  {ct}")
        return '\n'.join(lines)

    def to_dict(self) -> dict:
        return {
            'name': self.name, 'path': str(self.path), 'source': self.source, 'tissue': self.tissue,
            'database': self.database, 'min_percentile': self.min_percentile,
            'n_cells': self.n_cells, 'n_celltypes': self.n_celltypes, 'celltypes': self.celltypes,
        }

    def __repr__(self) -> str:
        return f"ReferencePanel(name={self.name!r}, source={self.source!r}, n_celltypes={self.n_celltypes}, n_cells={self.n_cells:,}, path={str(self.path)!r})"


# ── standard run ─────────────────────────────────────────────────────────────

@dataclass
class RunResult:
    """Outputs of one standard FastCCC run (``fastccc.run``)."""
    result_dir: Path
    task_id: str
    name: str = 'sample'
    database: str = ''
    database_path: Optional[Path] = None
    report_path: Optional[Path] = None
    manifest: dict = field(default_factory=dict)
    _cache: dict = field(default_factory=dict, repr=False)

    # ── construction ────────────────────────────────────────────────────────
    @classmethod
    def load(cls, result_dir: str | os.PathLike, task_id: Optional[str] = None, database: Optional[str] = None) -> 'RunResult':
        result_dir = Path(result_dir).expanduser().resolve()
        manifest = read_manifest(result_dir) or {}
        if manifest and manifest.get('workflow') not in (None, 'run'):
            raise ValueError(f"{result_dir} holds a '{manifest.get('workflow')}' result, not a standard run.")
        task_id = task_id or manifest.get('task_id') or find_task_id(result_dir)
        db_path = database or manifest.get('database', {}).get('path')
        report = _relocate(manifest.get('outputs', {}).get('report_html'), result_dir / 'report' / 'report.html')
        return cls(
            result_dir=result_dir, task_id=task_id, name=manifest.get('name', result_dir.name),
            database=manifest.get('database', {}).get('name', Path(str(db_path)).name if db_path else ''),
            database_path=Path(db_path) if db_path else None,
            report_path=report if report and report.exists() else None, manifest=manifest,
        )

    # ── file paths ──────────────────────────────────────────────────────────
    @property
    def significant_path(self) -> Path:
        return self.result_dir / f'{self.task_id}_significant_results.tsv'

    @property
    def pvals_path(self) -> Path:
        return self.result_dir / f'{self.task_id}_Cauchy_pvals.tsv'

    @property
    def strength_path(self) -> Path:
        return self.result_dir / f'{self.task_id}_average_interactions_strength.tsv'

    # ── lazily loaded tables ────────────────────────────────────────────────
    @property
    def significant(self) -> pd.DataFrame:
        """Significant L-R interactions (long format, p < 0.05)."""
        if 'significant' not in self._cache:
            self._cache['significant'] = pd.read_csv(self.significant_path, sep='\t')
        return self._cache['significant']

    @property
    def pvals(self) -> pd.DataFrame:
        """Cauchy-combined p-values, cell-type pair × LRI."""
        if 'pvals' not in self._cache:
            self._cache['pvals'] = pd.read_csv(self.pvals_path, sep='\t', index_col=0)
        return self._cache['pvals']

    @property
    def strength(self) -> pd.DataFrame:
        """Communication score averaged over scoring modules, cell-type pair × LRI."""
        if 'strength' not in self._cache:
            self._cache['strength'] = pd.read_csv(self.strength_path, sep='\t', index_col=0)
        return self._cache['strength']

    @property
    def celltypes(self) -> list[str]:
        cts = set()
        for pair in self.pvals.index:
            a, b = str(pair).split('|', 1)
            cts.update((a, b))
        return sorted(cts)

    # ── summaries ───────────────────────────────────────────────────────────
    def summary(self) -> dict:
        sig = self.significant
        pairs = sig['sender_celltype'].astype(str) + '|' + sig['receiver_celltype'].astype(str)
        return {
            'n_celltypes': len(self.celltypes),
            'n_celltype_pairs_tested': int(self.pvals.shape[0]),
            'n_lr_tested': int(self.pvals.shape[1]),
            'n_significant': int(len(sig)),
            'n_unique_lr_significant': int(sig['LRI_ID'].nunique()),
            'n_active_celltype_pairs': int(pairs.nunique()),
            'top_senders': _top_counts(sig['sender_celltype']),
            'top_receivers': _top_counts(sig['receiver_celltype']),
            'top_lr_pairs': [
                {'ligand': str(l), 'receptor': str(r), 'n_celltype_pairs': int(n)}
                for (l, r), n in sig.groupby(['ligand', 'receptor']).size().sort_values(ascending=False).head(5).items()
            ],
        }

    def to_dict(self) -> dict:
        d = {
            'workflow': 'run', 'name': self.name, 'task_id': self.task_id,
            'result_dir': str(self.result_dir), 'database': self.database,
            'report_html': str(self.report_path) if self.report_path else None,
            'files': {
                'significant_results': str(self.significant_path),
                'cauchy_pvals': str(self.pvals_path),
                'average_interactions_strength': str(self.strength_path),
            },
        }
        try:
            d['summary'] = self.summary()
        except Exception as exc:   # pragma: no cover - summary is best effort
            d['summary_error'] = str(exc)
        return d

    def report(self, output_dir: Optional[str | os.PathLike] = None, **kwargs) -> Path:
        """(Re)generate the HTML report for this run."""
        from .report import generate_report
        from .resources import resolve_database
        db = resolve_database(str(self.database_path or self.database))
        out = Path(output_dir) if output_dir else self.result_dir / 'report'
        kwargs.setdefault('sample_name', self.name)
        path = Path(generate_report(str(self.result_dir), self.task_id, str(db), output_dir=str(out), **kwargs))
        self.report_path = path
        update_manifest(self.result_dir, outputs={**self.manifest.get('outputs', {}), 'report_html': str(path)})
        return path

    def __repr__(self) -> str:
        try:
            s = self.summary()
            body = f"n_celltypes={s['n_celltypes']}, n_significant={s['n_significant']:,}"
        except Exception:
            body = 'tables not loaded'
        return f"RunResult(name={self.name!r}, task_id={self.task_id!r}, {body}, result_dir={str(self.result_dir)!r})"


@dataclass
class MultiRunResult:
    """Outputs of ``fastccc.run(..., split_by=...)``: one RunResult per group."""
    result_dir: Path
    split_by: str
    runs: dict = field(default_factory=dict)          # group label -> RunResult
    all: Optional[RunResult] = None
    report_path: Optional[Path] = None
    manifest: dict = field(default_factory=dict)

    @classmethod
    def load(cls, result_dir: str | os.PathLike) -> 'MultiRunResult':
        result_dir = Path(result_dir).expanduser().resolve()
        manifest = read_manifest(result_dir)
        if not manifest or manifest.get('workflow') != 'run_split':
            raise ValueError(f"{result_dir} does not hold a split run manifest.")
        runs = {g: RunResult.load(_relocate(p, result_dir / Path(p).name)) for g, p in manifest['groups'].items()}
        all_run = RunResult.load(_relocate(manifest['all'], result_dir / 'all')) if manifest.get('all') else None
        report = _relocate(manifest.get('outputs', {}).get('report_html'), result_dir / 'report' / 'report.html')
        return cls(result_dir=result_dir, split_by=manifest['split_by'], runs=runs, all=all_run,
                   report_path=report if report and report.exists() else None, manifest=manifest)

    @property
    def groups(self) -> list[str]:
        return list(self.runs)

    def __getitem__(self, group: str) -> RunResult:
        return self.runs[group]

    def to_dict(self) -> dict:
        return {
            'workflow': 'run_split', 'split_by': self.split_by, 'result_dir': str(self.result_dir),
            'report_html': str(self.report_path) if self.report_path else None,
            'all': self.all.to_dict() if self.all else None,
            'groups': {g: r.to_dict() for g, r in self.runs.items()},
        }

    def __repr__(self) -> str:
        return f"MultiRunResult(split_by={self.split_by!r}, groups={self.groups}, all={'yes' if self.all else 'no'}, result_dir={str(self.result_dir)!r})"


# ── reference-based inference ────────────────────────────────────────────────

TREND_ORDER = ['Up', 'Both Sig', 'Down', 'Both NS']


@dataclass
class ReferenceResult:
    """Outputs of ``fastccc.compare_to_reference``."""
    result_dir: Path
    query_name: str
    reference: ReferencePanel
    database: str = ''
    database_path: Optional[Path] = None
    report_path: Optional[Path] = None
    manifest: dict = field(default_factory=dict)
    _cache: dict = field(default_factory=dict, repr=False)

    @classmethod
    def load(cls, result_dir: str | os.PathLike) -> 'ReferenceResult':
        result_dir = Path(result_dir).expanduser().resolve()
        manifest = read_manifest(result_dir) or {}
        if manifest and manifest.get('workflow') != 'compare_to_reference':
            raise ValueError(f"{result_dir} holds a '{manifest.get('workflow')}' result, not a reference comparison.")
        if not (result_dir / 'query_infer_results.tsv').is_file():
            raise FileNotFoundError(f"query_infer_results.tsv not found in {result_dir}")
        ref = manifest.get('reference', {})
        ref_path = _relocate(ref.get('path'), result_dir / 'reference' / Path(ref['path']).name if ref.get('path') else None)
        panel = ReferencePanel.load(ref_path, source=ref.get('source', 'custom'), tissue=ref.get('tissue')) \
            if ref_path and ref_path.is_dir() else ReferencePanel(name=ref.get('name', 'reference'), path=ref_path or Path(''), source=ref.get('source', 'custom'))
        db = manifest.get('database', {})
        report = _relocate(manifest.get('outputs', {}).get('report_html'), result_dir / 'report' / 'report.html')
        return cls(result_dir=result_dir, query_name=manifest.get('name', 'Query'), reference=panel,
                   database=db.get('name', ''), database_path=Path(db['path']) if db.get('path') else None,
                   report_path=report if report and report.exists() else None, manifest=manifest)

    @property
    def results_path(self) -> Path:
        return self.result_dir / 'query_infer_results.tsv'

    @property
    def strength_path(self) -> Path:
        return self.result_dir / 'query_interactions_strength.tsv'

    @property
    def results(self) -> pd.DataFrame:
        """Per interaction: query CS, reference null, significance in both, ``trend_vs_ref``."""
        if 'results' not in self._cache:
            self._cache['results'] = pd.read_csv(self.results_path, sep='\t', low_memory=False)
        return self._cache['results']

    @property
    def strength(self) -> pd.DataFrame:
        if 'strength' not in self._cache:
            self._cache['strength'] = pd.read_csv(self.strength_path, sep='\t', index_col=0)
        return self._cache['strength']

    @property
    def trend_counts(self) -> dict:
        tv = self.results['trend_vs_ref']
        counts = {t: int((tv == t).sum()) for t in TREND_ORDER}
        counts['Not in reference'] = int(tv.isna().sum())
        return counts

    def interactions(self, trend: str) -> pd.DataFrame:
        """Rows with a given ``trend_vs_ref`` label (``'Up'``, ``'Down'``, ...)."""
        return self.results[self.results['trend_vs_ref'] == trend].copy()

    def summary(self) -> dict:
        df = self.results
        in_ref = df[df['trend_vs_ref'].notna()]
        up, down = in_ref[in_ref['trend_vs_ref'] == 'Up'], in_ref[in_ref['trend_vs_ref'] == 'Down']
        def _pairs(sub):
            return _top_counts(sub['sender|receiver']) if len(sub) else []
        return {
            'n_tested': int(len(df)),
            'n_in_reference': int(len(in_ref)),
            'n_significant_query': int(df['is_significant'].sum()),
            'trend_counts': self.trend_counts,
            'n_celltype_pairs': int(in_ref['sender|receiver'].nunique()),
            'top_up_celltype_pairs': _pairs(up),
            'top_down_celltype_pairs': _pairs(down),
            'celltype_coverage': self.manifest.get('celltype_coverage'),
            'celltype_mapping_source': self.manifest.get('celltype_mapping_source'),
        }

    def to_dict(self) -> dict:
        d = {
            'workflow': 'compare_to_reference', 'name': self.query_name, 'result_dir': str(self.result_dir),
            'database': self.database, 'reference': self.reference.to_dict(),
            'report_html': str(self.report_path) if self.report_path else None,
            'files': {'results': str(self.results_path), 'interactions_strength': str(self.strength_path),
                      'percents_analysis': str(self.result_dir / 'query_percents_analysis.tsv')},
        }
        try:
            d['summary'] = self.summary()
        except Exception as exc:   # pragma: no cover
            d['summary_error'] = str(exc)
        return d

    def report(self, output_dir: Optional[str | os.PathLike] = None, **kwargs) -> Path:
        from .report import generate_infer_report
        from .resources import resolve_database
        db = resolve_database(str(self.database_path or self.database))
        out = Path(output_dir) if output_dir else self.result_dir / 'report'
        kwargs.setdefault('query_name', self.query_name)
        kwargs.setdefault('reference_name', self.manifest.get('reference_display_name'))
        if self.reference.source == 'healthy_tissue' and self.reference.tissue:
            kwargs.setdefault('reference_tissue', self.reference.tissue)
            kwargs.setdefault('reference_root', str(self.reference.path.parent))
        else:
            kwargs.setdefault('reference_path', str(self.reference.path))
        path = Path(generate_infer_report(str(self.result_dir), str(db), output_dir=str(out), **kwargs))
        self.report_path = path
        update_manifest(self.result_dir, outputs={**self.manifest.get('outputs', {}), 'report_html': str(path)})
        return path

    @property
    def celltype_mapping(self) -> dict:
        """``{reference cell type: query cell type}`` actually used (empty when labels matched by name)."""
        return dict(self.manifest.get('celltype_mapping') or {})

    @property
    def celltype_coverage(self) -> dict:
        return dict(self.manifest.get('celltype_coverage') or {})

    def interaction_summary(self, trend: str = 'Down', n: Optional[int] = 20) -> pd.DataFrame:
        """One row per ligand-receptor pair with the given trend, ranked by how many
        sender-receiver pairs show it (then by the strongest query score).

        Use this instead of sorting rows by ``comm_score``: for ``Down`` the query
        score measures what is left, not what was lost, and the same LR pair often
        repeats across many cell-type pairs.
        """
        df = self.interactions(trend)
        if df.empty:
            return pd.DataFrame(columns=['ligand', 'receptor', 'n_celltype_pairs', 'top_sender_receiver', 'max_comm_score',
                                         'ligand_pct_query', 'ligand_pct_ref', 'receptor_pct_query', 'receptor_pct_ref'])
        df = df.copy()
        df['comm_score'] = pd.to_numeric(df['comm_score'], errors='coerce')
        rows = []
        for (lig, rec), g in df.groupby(['ligand', 'receptor'], sort=False):
            top = g.loc[g['comm_score'].idxmax()] if g['comm_score'].notna().any() else g.iloc[0]
            rows.append({
                'ligand': lig, 'receptor': rec, 'n_celltype_pairs': int(len(g)),
                'top_sender_receiver': top['sender|receiver'], 'max_comm_score': float(g['comm_score'].max()),
                'ligand_pct_query': float(pd.to_numeric(g['ligand_expr_percent'], errors='coerce').mean()),
                'ligand_pct_ref': float(pd.to_numeric(g['ligand_expr_percent_ref'], errors='coerce').mean()),
                'receptor_pct_query': float(pd.to_numeric(g['receptor_expr_percent'], errors='coerce').mean()),
                'receptor_pct_ref': float(pd.to_numeric(g['receptor_expr_percent_ref'], errors='coerce').mean()),
            })
        out = pd.DataFrame(rows).sort_values(['n_celltype_pairs', 'max_comm_score'], ascending=[False, False]).reset_index(drop=True)
        return out.head(n) if n else out

    def __repr__(self) -> str:
        try:
            tc = self.trend_counts
            body = (f"Up={tc['Up']:,}, Down={tc['Down']:,}, BothSig={tc['Both Sig']:,}, BothNS={tc['Both NS']:,}"
                    + (f", NotInReference={tc['Not in reference']:,}" if tc.get('Not in reference') else ''))
        except Exception:
            body = 'tables not loaded'
        return f"ReferenceResult(query={self.query_name!r}, reference={self.reference.name!r}, {body}, result_dir={str(self.result_dir)!r})"


def load_result(result_dir: str | os.PathLike):
    """Load whichever result type a directory holds, based on its manifest."""
    manifest = read_manifest(result_dir)
    wf = (manifest or {}).get('workflow')
    if wf == 'run_split':
        return MultiRunResult.load(result_dir)
    if wf == 'compare_to_reference' or (wf is None and (Path(result_dir) / 'query_infer_results.tsv').is_file()):
        return ReferenceResult.load(result_dir)
    return RunResult.load(result_dir)


for _cls in (RunResult, MultiRunResult, ReferenceResult):
    _cls.explore = _explorer_mixin_explore
    _cls.explorer_path = property(lambda self: getattr(self, '_explorer_path', None)
                                  or _explorer_path_from_manifest(getattr(self, 'manifest', None), getattr(self, 'result_dir', None)))
