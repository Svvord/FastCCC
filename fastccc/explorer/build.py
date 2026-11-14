"""Render the explorer HTML for a FastCCC result folder."""
from __future__ import annotations

import datetime
import importlib.metadata
import json
import os
from pathlib import Path
from typing import Optional

from loguru import logger

from .payload import build_payload
from ..results import update_manifest

HERE = Path(__file__).resolve().parent
KIND_LABELS = {'standard': 'cell-cell communication', 'split': 'per-condition communication', 'reference': 'reference comparison'}


def render_explorer_html(payload: dict) -> str:
    from jinja2 import Environment, FileSystemLoader
    env = Environment(loader=FileSystemLoader(str(HERE / 'templates')), autoescape=False)
    tmpl = env.get_template('explorer.html')
    meta = payload['meta']
    if payload['kind'] == 'reference':
        ref = meta.get('reference') or {}
        title = f"{meta.get('name')} vs {ref.get('display_name') or ref.get('name')}"
        subtitle = f"{meta.get('n_cells') or '?'} cells · reference {ref.get('n_cells') or '?'} cells · FastCCC {meta.get('fastccc_version') or ''}"
    elif payload['kind'] == 'split':
        title = f"{meta.get('name')} by {meta.get('split_by')}"
        subtitle = ' · '.join(f"{c['name']}: {c['n_significant']} significant" for c in meta.get('conditions', []))
    else:
        title = str(meta.get('name'))
        subtitle = f"{meta.get('n_cells') or '?'} cells · {len(payload['celltypes'])} cell types · {len(payload['rows']['sender'])} significant interactions"
    payload_json = json.dumps(payload, separators=(',', ':'), ensure_ascii=False).replace('</', '<\\/')
    return tmpl.render(
        title=f"FastCCC · {title}", kind_label=KIND_LABELS.get(payload['kind'], payload['kind']), subtitle=subtitle,
        version=meta.get('fastccc_version') or importlib.metadata.version('fastccc'),
        created=datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
        d3_js=(HERE / 'static' / 'd3.v7.min.js').read_text(encoding='utf-8'),
        sankey_js=(HERE / 'static' / 'd3-sankey.min.js').read_text(encoding='utf-8'),
        explorer_js=(HERE / 'static' / 'explorer.js').read_text(encoding='utf-8'),
        payload_json=payload_json,
    )


def build_explorer(result_dir: str | os.PathLike, output: Optional[str | os.PathLike] = None,
                   include_both_ns: bool = True, genesets: bool = True, cells: bool = True, max_cells: int = 4000,
                   patterns: bool = True, n_patterns: Optional[int] = None, nulls: bool = True) -> Path:
    """Write ``explorer.html`` (default: inside the result folder) and return its path.

    ``genesets`` downloads Enrichr libraries once (cached under ``FASTCCC_HOME``);
    ``cells`` embeds a stratified subsample of per-cell expression for violins and
    embedding overlays; ``patterns`` runs the NMF communication-pattern analysis;
    ``nulls`` embeds the reference null curves (reference comparisons).
    """
    from .payload import ExplorerOptions
    result_dir = Path(result_dir).expanduser()
    logger.info(f"Building the result explorer for {result_dir} …")
    payload = build_payload(result_dir, ExplorerOptions(genesets=genesets, cells=cells, max_cells=max_cells, patterns=patterns,
                                                        n_patterns=n_patterns, nulls=nulls, include_both_ns=include_both_ns))
    html = render_explorer_html(payload)
    out = Path(output).expanduser() if output else result_dir / 'explorer.html'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding='utf-8')
    size_mb = out.stat().st_size / 1e6
    logger.success(f"Explorer written: {out} ({size_mb:.1f} MB, {len(payload['rows']['sender']):,} interactions)")
    try:
        manifest_dir = result_dir if (result_dir / 'fastccc_run.json').is_file() else None
        if manifest_dir is not None:
            from ..results import read_manifest
            m = read_manifest(manifest_dir) or {}
            update_manifest(manifest_dir, outputs={**(m.get('outputs') or {}), 'explorer_html': str(out)})
    except Exception as err:   # manifest bookkeeping must never fail the build
        logger.debug(f"manifest not updated: {err}")
    return out
