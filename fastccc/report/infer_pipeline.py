"""generate_infer_report() — HTML report for reference-based CCC inference."""

import datetime
import importlib.metadata
import os
import tomllib
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd
from loguru import logger

from .utils import fig_to_base64, save_figure, set_publication_style
from .modules.infer_analysis import (
    build_celltype_reference_profiles,
    InferData,
    plot_trend_distribution,
    plot_ct_pair_trend_breakdown,
    plot_trend_heatmap,
    plot_top_lr_dotplot,
    plot_cs_scatter,
    plot_pathway_breakdown,
)


def _summarize_figure_audit(entries):
    counts = {'generated': 0, 'skipped': 0, 'failed': 0}
    unavailable = []
    for entry in entries:
        counts[entry['status']] += 1
        if entry['status'] != 'generated':
            unavailable.append(entry)
    return {
        'counts': counts,
        'total': len(entries),
        'unavailable': unavailable,
    }


def _find_reference_root(reference_root: Optional[str] = None) -> Path:
    if reference_root is not None:
        root = Path(reference_root).expanduser()
        if root.is_dir():
            return root
        raise FileNotFoundError(f"Reference root does not exist: {root}")

    candidates = [
        Path.cwd() / 'reference',
        Path(__file__).resolve().parents[2] / 'reference',
    ]
    for root in candidates:
        if root.is_dir():
            return root

    raise FileNotFoundError(
        "Could not find a reference panel root. Pass reference_root after "
        "downloading the FastCCC reference panels."
    )


def _load_reference_panel(reference_path: str, source: str, tissue: str = "") -> dict:
    panel_path = Path(reference_path).expanduser()
    config_path = panel_path / 'config.toml'
    if not panel_path.is_dir():
        raise FileNotFoundError(f"Reference panel does not exist: {panel_path}")
    if not config_path.is_file():
        raise FileNotFoundError(
            f"Reference panel is missing config.toml: {panel_path}"
        )

    with config_path.open('rb') as fh:
        config = tomllib.load(fh)

    config_name = config.get('reference_name', panel_path.name)
    source_label = (
        "FastCCC healthy tissue panel"
        if source == 'healthy_tissue' else "Custom reference panel"
    )
    return {
        'path': str(panel_path.resolve()),
        'source': source,
        'source_label': source_label,
        'tissue': tissue,
        'reference_name': config_name,
        'database_name': config.get('LRI_database', ''),
        'min_percentile': config.get('min_percentile', ''),
        'n_celltypes': len(config.get('celltype', {})),
    }


def list_reference_panels(reference_root: Optional[str] = None) -> list[str]:
    """Return healthy tissue panel names available under a reference root."""
    root = _find_reference_root(reference_root)
    return [
        path.name
        for path in sorted(root.iterdir())
        if path.is_dir() and (path / 'config.toml').is_file()
    ]


def _resolve_reference_panel(
    reference_path: Optional[str] = None,
    reference_tissue: Optional[str] = None,
    reference_root: Optional[str] = None,
) -> Optional[dict]:
    if reference_path is None and reference_tissue is None:
        return None
    if reference_path is not None and reference_tissue is not None:
        raise ValueError("Pass reference_path or reference_tissue, not both.")

    if reference_tissue is not None:
        root = _find_reference_root(reference_root)
        available = list_reference_panels(str(root))
        if reference_tissue not in available:
            options = ', '.join(available)
            raise ValueError(
                f"Unknown reference_tissue '{reference_tissue}'. "
                f"Available panels under {root}: {options}"
            )
        return _load_reference_panel(
            str(root / reference_tissue),
            source='healthy_tissue',
            tissue=reference_tissue,
        )

    return _load_reference_panel(reference_path, source='custom')


def _default_reference_name(panel: Optional[dict]) -> str:
    if panel is None:
        return "Reference"

    config_name = str(panel['reference_name']).replace('_', ' ')
    if panel['source'] == 'healthy_tissue':
        return f"Healthy {config_name.title()}"
    return config_name


def _validate_reference_database(panel: dict, database_path: str) -> None:
    panel_db = panel['database_name']
    requested_db = os.path.basename(database_path.rstrip('/\\'))
    if panel_db and requested_db and panel_db != requested_db:
        raise ValueError(
            f"Reference panel '{panel['reference_name']}' was built with "
            f"{panel_db}, but database_path points to {requested_db}."
        )


def generate_infer_report(
    infer_result_dir: str,
    database_path: str,
    output_dir: Optional[str] = None,
    query_name: str = "Query",
    reference_name: Optional[str] = None,
    reference_path: Optional[str] = None,
    reference_tissue: Optional[str] = None,
    reference_root: Optional[str] = None,
    dpi: int = 150,
    save_individual_figures: bool = True,
    top_n_lr: int = 25,
    top_n_celltypes: int = 20,
) -> str:
    """
    Generate an HTML report from reference-based FastCCC inference results.

    Parameters
    ----------
    infer_result_dir
        Directory containing ``query_infer_results.tsv`` and
        ``query_interactions_strength.tsv`` (output of ``infer_query_workflow``).
    database_path
        Path to the LRI database folder (for pathway classification lookup).
    output_dir
        Where to write ``infer_report.html`` and individual figure PNGs.
        Defaults to ``<infer_result_dir>/infer_report/``.
    query_name
        Display name for the query condition (e.g. "PBC").
    reference_name
        Display name for the reference. If omitted with ``reference_tissue``,
        the report uses a "Healthy <tissue>" label.
    reference_path
        Optional custom reference panel path used for panel metadata in the
        report when inference results were generated separately.
    reference_tissue
        Optional healthy tissue panel name under ``reference_root`` (for
        example, ``"liver"``). This is mutually exclusive with
        ``reference_path``.
    reference_root
        Folder containing FastCCC healthy tissue reference panels. If omitted,
        FastCCC checks ``./reference`` and the source checkout reference folder.

    Returns
    -------
    str — absolute path to the generated ``infer_report.html``.
    """
    set_publication_style()
    reference_panel = _resolve_reference_panel(
        reference_path=reference_path,
        reference_tissue=reference_tissue,
        reference_root=reference_root,
    )
    if reference_panel is not None:
        _validate_reference_database(reference_panel, database_path)
    if reference_name is None:
        reference_name = _default_reference_name(reference_panel)

    if output_dir is None:
        output_dir = os.path.join(infer_result_dir, 'infer_report')
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # ── Load data ────────────────────────────────────────────────────────────
    logger.info("Loading inference results…")
    results = pd.read_csv(
        os.path.join(infer_result_dir, 'query_infer_results.tsv'),
        sep='\t',
        low_memory=False,
    )
    strength_path = os.path.join(infer_result_dir, 'query_interactions_strength.tsv')
    strength = pd.read_csv(strength_path, sep='\t', index_col=0) if os.path.exists(strength_path) else pd.DataFrame()

    data = InferData(
        results        = results,
        strength       = strength,
        query_name     = query_name,
        reference_name = reference_name,
        database_path  = database_path,
    )

    logger.info(
        f"Loaded {data.n_tested:,} interactions — "
        f"{data.n_significant:,} significant in {query_name}, "
        f"{data.n_in_reference:,} in reference"
    )

    figure_audit = []

    # ── Figure runner ────────────────────────────────────────────────────────
    def _run(name: str, func, *args, **kwargs):
        logger.info(f"  {name}…")
        audit = {
            'artifact': name,
            'function': func.__name__,
            'status': 'generated',
            'note': '',
        }
        try:
            result = func(*args, **kwargs)
            fig, caption = result if isinstance(result, tuple) else (result, "")
            if fig is None:
                audit['status'] = 'skipped'
                audit['note'] = caption
                figure_audit.append(audit)
                return None, caption
            if save_individual_figures:
                save_figure(fig, out, name, dpi=dpi)
            b64 = fig_to_base64(fig)
            plt.close(fig)
            figure_audit.append(audit)
            return b64, caption
        except Exception as e:
            logger.warning(f"  {name} failed: {e}")
            plt.close('all')
            audit['status'] = 'failed'
            audit['note'] = str(e)
            figure_audit.append(audit)
            return None, f"Figure could not be generated: {e}"

    # ── Generate figures ─────────────────────────────────────────────────────
    logger.info("Generating figures…")
    figs = {}
    figs['trend_dist']      = _run("fig_A_trend_distribution",   plot_trend_distribution,       data)
    figs['ct_breakdown']    = _run("fig_B_ct_pair_breakdown",     plot_ct_pair_trend_breakdown,  data, top_n_celltypes)
    figs['heatmap_up']      = _run("fig_C_heatmap_up",           plot_trend_heatmap,            data, 'Up',   top_n_celltypes)
    figs['heatmap_down']    = _run("fig_D_heatmap_down",         plot_trend_heatmap,            data, 'Down', top_n_celltypes)
    figs['dotplot_up']      = _run("fig_E_dotplot_up",           plot_top_lr_dotplot,           data, 'Up',   top_n_lr, top_n_celltypes)
    figs['dotplot_down']    = _run("fig_F_dotplot_down",         plot_top_lr_dotplot,           data, 'Down', top_n_lr, top_n_celltypes)
    figs['cs_scatter']      = _run("fig_G_cs_scatter",           plot_cs_scatter,               data)
    figs['pathway_breakdown']= _run("fig_H_pathway_breakdown",   plot_pathway_breakdown,        data, top_n_celltypes)
    celltype_reference = build_celltype_reference_profiles(data)

    figure_audit_summary = _summarize_figure_audit(figure_audit)

    # ── Render HTML ──────────────────────────────────────────────────────────
    logger.info("Rendering HTML report…")
    from jinja2 import Environment, FileSystemLoader
    env  = Environment(loader=FileSystemLoader(str(Path(__file__).parent / 'templates')))
    tmpl = env.get_template('infer_report.html')

    html = tmpl.render(
        query_name      = query_name,
        reference_name  = reference_name,
        report_date     = datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
        database_name   = os.path.basename(database_path.rstrip('/\\')),
        version         = importlib.metadata.version("fastccc"),
        reference_panel = reference_panel,
        # stats
        n_tested        = data.n_tested,
        n_significant   = data.n_significant,
        n_in_reference  = data.n_in_reference,
        n_lr_pairs      = data.n_lr_pairs,
        n_ct_pairs      = data.n_ct_pairs,
        trend_counts    = data.trend_counts,
        trend_order     = ['Up', 'Both Sig', 'Down', 'Both NS'],
        trend_colors    = {'Up': '#E64B35', 'Down': '#4DBBD5', 'Both Sig': '#7E6148', 'Both NS': '#dddddd'},
        # figures
        figs            = figs,
        celltype_reference = celltype_reference,
        figure_audit    = figure_audit_summary,
    )

    report_path = out / 'infer_report.html'
    report_path.write_text(html, encoding='utf-8')
    logger.success(f"Report saved to: {report_path}")
    return str(report_path)


def generate_reference_report(
    database_path: str,
    query_counts_file_path,
    infer_result_dir: str,
    celltype_file_path: Optional[str] = None,
    output_dir: Optional[str] = None,
    reference_path: Optional[str] = None,
    reference_tissue: Optional[str] = None,
    reference_root: Optional[str] = None,
    celltype_mapping_dict=None,
    meta_key: Optional[str] = None,
    query_name: str = "Query",
    reference_name: Optional[str] = None,
    dpi: int = 150,
    save_individual_figures: bool = True,
    top_n_lr: int = 25,
    top_n_celltypes: int = 20,
    debug_mode: bool = False,
) -> str:
    """
    Run reference inference and render a reference HTML report.

    Choose a shipped healthy tissue panel with ``reference_tissue`` and
    ``reference_root`` or pass ``reference_path`` for a custom panel.
    """
    if meta_key is None and celltype_file_path is None:
        raise ValueError("Pass meta_key or celltype_file_path for query labels.")

    panel = _resolve_reference_panel(
        reference_path=reference_path,
        reference_tissue=reference_tissue,
        reference_root=reference_root,
    )
    if panel is None:
        raise ValueError("Pass reference_tissue or reference_path.")
    _validate_reference_database(panel, database_path)

    from ..infer_query import infer_query_workflow

    infer_query_workflow(
        database_file_path=database_path,
        reference_path=panel['path'],
        query_counts_file_path=query_counts_file_path,
        celltype_file_path=celltype_file_path,
        save_path=infer_result_dir,
        celltype_mapping_dict=celltype_mapping_dict,
        meta_key=meta_key,
        debug_mode=debug_mode,
    )

    report_reference_path = panel['path'] if panel['source'] == 'custom' else None
    report_reference_tissue = (
        panel['tissue'] if panel['source'] == 'healthy_tissue' else None
    )
    return generate_infer_report(
        infer_result_dir=infer_result_dir,
        database_path=database_path,
        output_dir=output_dir,
        query_name=query_name,
        reference_name=reference_name,
        reference_path=report_reference_path,
        reference_tissue=report_reference_tissue,
        reference_root=reference_root,
        dpi=dpi,
        save_individual_figures=save_individual_figures,
        top_n_lr=top_n_lr,
        top_n_celltypes=top_n_celltypes,
    )
