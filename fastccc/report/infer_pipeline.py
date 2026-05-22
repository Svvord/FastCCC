"""generate_infer_report() — HTML report for reference-based CCC inference."""

import datetime
import importlib.metadata
import os
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd
from loguru import logger

from .utils import fig_to_base64, save_figure, set_publication_style
from .modules.infer_analysis import (
    InferData,
    plot_trend_distribution,
    plot_ct_pair_trend_breakdown,
    plot_trend_heatmap,
    plot_top_lr_dotplot,
    plot_cs_scatter,
    plot_pathway_breakdown,
)


def generate_infer_report(
    infer_result_dir: str,
    database_path: str,
    output_dir: Optional[str] = None,
    query_name: str = "Query",
    reference_name: str = "Reference",
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
        Display name for the reference (e.g. "PSC").

    Returns
    -------
    str — absolute path to the generated ``infer_report.html``.
    """
    set_publication_style()

    if output_dir is None:
        output_dir = os.path.join(infer_result_dir, 'infer_report')
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # ── Load data ────────────────────────────────────────────────────────────
    logger.info("Loading inference results…")
    results = pd.read_csv(os.path.join(infer_result_dir, 'query_infer_results.tsv'), sep='\t')
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

    # ── Figure runner ────────────────────────────────────────────────────────
    def _run(name: str, func, *args, **kwargs):
        logger.info(f"  {name}…")
        try:
            result = func(*args, **kwargs)
            fig, caption = result if isinstance(result, tuple) else (result, "")
            if fig is None:
                return None, caption
            if save_individual_figures:
                save_figure(fig, out, name, dpi=dpi)
            b64 = fig_to_base64(fig)
            plt.close(fig)
            return b64, caption
        except Exception as e:
            logger.warning(f"  {name} failed: {e}")
            plt.close('all')
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
    )

    report_path = out / 'infer_report.html'
    report_path.write_text(html, encoding='utf-8')
    logger.success(f"Report saved to: {report_path}")
    return str(report_path)
