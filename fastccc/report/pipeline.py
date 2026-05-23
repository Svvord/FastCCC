"""Main entry point: generate_report()."""

import datetime
import importlib.metadata
import os
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import pandas as pd
from loguru import logger

from .loader import load_results, CCCData
from .utils import (
    fig_to_base64, get_celltype_colors, save_figure, set_publication_style
)
from .modules.overview import (
    plot_chord_diagram, plot_count_heatmap,
    plot_sender_receiver_bar, plot_strength_heatmap,
)
from .modules.lr_analysis import (
    build_celltype_lr_profiles, build_celltype_pathway_profiles,
    plot_classification_bar, plot_lr_dotplot, plot_pathway_celltype_heatmap,
)
from .modules.enrichment import (
    build_celltype_gene_profiles, plot_ligand_ora, plot_receptor_ora,
    plot_tf_heatmap,
)
from .modules.celltype_profile import (
    plot_io_scatter, plot_interaction_flow, plot_sender_pathway_heatmap,
    plot_receiver_pathway_heatmap, plot_incoming_flow,
)
from .modules.network import (
    plot_network_centrality, plot_sankey_flow,
    plot_autocrine_paracrine, plot_bipartite_lr,
    plot_communication_asymmetry, plot_network_communities,
)
from .modules.advanced import (
    build_celltype_pathway_flow_profiles, plot_pathway_info_flow,
    plot_pathway_lr_multiples,
    plot_lr_specificity, plot_cs_violin,
    plot_lr_cooccurrence, plot_pathway_crosstalk,
)
from .modules.differential import (
    plot_diff_heatmap, plot_diff_volcano, plot_diff_pathway_bar,
    plot_lr_stability,
)
from .modules.pathway_utils import keep_annotated_classifications


def _make_stats(data: CCCData) -> Dict:
    sig = data.significant
    n_sig = len(sig)
    mat   = data.counts_matrix
    annotated_sig = keep_annotated_classifications(sig)
    pathway_counts = (
        annotated_sig['classification'].value_counts()
        if not annotated_sig.empty and 'classification' in annotated_sig.columns
        else pd.Series(dtype=int)
    )
    return {
        'n_celltypes':       len(data.celltypes),
        'n_sig_interactions': n_sig,
        'n_lr_pairs':         int(sig['LRI_ID'].nunique()) if n_sig > 0 else 0,
        'n_ct_pairs':         int((mat > 0).values.sum()),
        'n_pathways':         int(annotated_sig['classification'].nunique()) if not annotated_sig.empty and 'classification' in annotated_sig.columns else 0,
        'top_sender':         mat.sum(axis=1).idxmax() if n_sig > 0 else "",
        'top_sender_count':   int(mat.sum(axis=1).max()) if n_sig > 0 else 0,
        'top_receiver':       mat.sum(axis=0).idxmax() if n_sig > 0 else "",
        'top_receiver_count': int(mat.sum(axis=0).max()) if n_sig > 0 else 0,
        'top_pathway':        pathway_counts.idxmax() if not pathway_counts.empty else "",
        'top_pathway_count':  int(pathway_counts.iloc[0]) if not pathway_counts.empty else 0,
    }


def _summarize_figure_audit(entries: List[Dict]) -> Dict:
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


def generate_report(
    result_dir: str,
    task_id: str,
    database_path: str,
    output_dir: Optional[str] = None,
    sample_name: str = "All Cells",
    pval_threshold: float = 0.05,
    top_n_lr: int = 30,
    top_n_celltypes: int = 20,
    max_chords: int = 180,
    gene_sets: List[str] = ("KEGG_2021_Human", "GO_Biological_Process_2023"),
    dpi: int = 300,
    save_individual_figures: bool = True,
    # Condition A (optional – enables per-condition tabs)
    cond_a_result_dir: Optional[str] = None,
    cond_a_task_id:    Optional[str] = None,
    cond_a_name:       Optional[str] = None,
    # Condition B (optional - requires cond_a; enables condition-comparison tab)
    cond_b_result_dir: Optional[str] = None,
    cond_b_task_id:    Optional[str] = None,
    cond_b_name:       Optional[str] = None,
) -> str:
    """
    Generate a comprehensive post-analysis HTML report from FastCCC results.

    Parameters
    ----------
    result_dir / task_id
        FastCCC results for the full combined dataset (always shown).
    sample_name
        Display name for the full dataset (default "All Cells").
    cond_a_* / cond_b_*
        Optional per-condition FastCCC results.  When provided the report gains
        clickable tabs so the user can switch between the full-dataset view,
        each condition's individual analysis, and (if both are given) a
        condition-comparison tab between the two independently analysed conditions.

    Returns
    -------
    str  –  absolute path to the generated ``report.html``.
    """
    set_publication_style()

    if output_dir is None:
        output_dir = os.path.join(result_dir, f"report_{task_id}")
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # ── Load datasets ────────────────────────────────────────────────────────
    logger.info("Loading FastCCC results…")
    data_all = load_results(result_dir, task_id, database_path, sample_name, pval_threshold)

    data_a, data_b = None, None
    name_a = cond_a_name or "Condition A"
    name_b = cond_b_name or "Condition B"

    if cond_a_result_dir and cond_a_task_id:
        logger.info(f"  Loading condition '{name_a}'…")
        data_a = load_results(cond_a_result_dir, cond_a_task_id, database_path, name_a, pval_threshold)
    if cond_b_result_dir and cond_b_task_id:
        logger.info(f"  Loading condition '{name_b}'…")
        data_b = load_results(cond_b_result_dir, cond_b_task_id, database_path, name_b, pval_threshold)

    has_conditions   = data_a is not None
    has_differential = data_a is not None and data_b is not None

    # ── Shared color palette (sorted cell-type names → consistent across tabs) ─
    all_cts = set(data_all.celltypes)
    if data_a: all_cts |= set(data_a.celltypes)
    if data_b: all_cts |= set(data_b.celltypes)
    colors = get_celltype_colors(list(all_cts))

    figure_audit = []

    # ── Core helper ──────────────────────────────────────────────────────────
    def _run(name: str, func, *args, **kwargs):
        logger.info(f"    {name}…")
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
            logger.warning(f"    {name} failed: {e}")
            plt.close('all')
            audit['status'] = 'failed'
            audit['note'] = str(e)
            figure_audit.append(audit)
            return None, f"Figure could not be generated: {e}"

    # ── Per-condition figure runner ───────────────────────────────────────────
    def _run_cond_figs(data: CCCData, prefix: str) -> Dict:
        """Generate all 21 analysis figures for one condition/dataset."""
        f = {}
        f['count_heatmap']  = _run(f"{prefix}_fig01", plot_count_heatmap,          data)
        f['strength_heatmap']= _run(f"{prefix}_fig02", plot_strength_heatmap,       data)
        f['chord']           = _run(f"{prefix}_fig03", plot_chord_diagram,           data, colors, 'count', max_chords)
        f['sr_bar']          = _run(f"{prefix}_fig04", plot_sender_receiver_bar,     data, colors, top_n_celltypes)
        f['lr_dotplot']      = _run(f"{prefix}_fig05", plot_lr_dotplot,              data, top_n_lr, top_n_celltypes)
        f['celltype_lrs']    = build_celltype_lr_profiles(data)
        f['class_bar']       = _run(f"{prefix}_fig06", plot_classification_bar,      data)
        f['celltype_pathways'] = build_celltype_pathway_profiles(data)
        f['path_ct']         = _run(f"{prefix}_fig07", plot_pathway_celltype_heatmap,data)
        f['lig_ora']         = _run(f"{prefix}_fig08", plot_ligand_ora,              data, gene_sets)
        f['celltype_ligands']= build_celltype_gene_profiles(data, role='ligand')
        f['rec_ora']         = _run(f"{prefix}_fig09", plot_receptor_ora,            data, gene_sets)
        f['celltype_receptors'] = build_celltype_gene_profiles(data, role='receptor')
        f['tf']              = _run(f"{prefix}_fig10", plot_tf_heatmap,              data)
        f['io_scatter']      = _run(f"{prefix}_fig11", plot_io_scatter,              data, colors)
        f['sender_pathway']  = _run(f"{prefix}_fig12", plot_sender_pathway_heatmap,  data)
        f['flow']            = _run(f"{prefix}_fig13", plot_interaction_flow,        data, colors, top_n_celltypes)
        f['network']         = _run(f"{prefix}_fig14", plot_network_centrality,      data, colors)
        f['sankey']          = _run(f"{prefix}_fig15", plot_sankey_flow,             data, colors)
        f['autocrine']       = _run(f"{prefix}_fig16", plot_autocrine_paracrine,     data, colors)
        f['bipartite']       = _run(f"{prefix}_fig17", plot_bipartite_lr,            data, colors)
        f['info_flow']          = _run(f"{prefix}_fig18", plot_pathway_info_flow,          data)
        f['celltype_pathway_flow'] = build_celltype_pathway_flow_profiles(data)
        f['lr_multiples']       = _run(f"{prefix}_fig19", plot_pathway_lr_multiples,       data)
        f['lr_spec']            = _run(f"{prefix}_fig20", plot_lr_specificity,             data)
        f['cs_violin']          = _run(f"{prefix}_fig21", plot_cs_violin,                  data)
        f['receiver_pathway']   = _run(f"{prefix}_fig25", plot_receiver_pathway_heatmap,   data)
        f['incoming_flow']      = _run(f"{prefix}_fig26", plot_incoming_flow,              data, colors)
        f['asymmetry']          = _run(f"{prefix}_fig27", plot_communication_asymmetry,    data, colors)
        f['communities']        = _run(f"{prefix}_fig28", plot_network_communities,        data, colors)
        f['lr_cooccurrence']    = _run(f"{prefix}_fig29", plot_lr_cooccurrence,            data)
        f['pathway_crosstalk']  = _run(f"{prefix}_fig30", plot_pathway_crosstalk,          data)
        return f

    # ── Build tabs ───────────────────────────────────────────────────────────
    tabs = []

    logger.info(f"Generating figures — {sample_name} (full dataset)")
    tabs.append({
        'id':    'all',
        'name':  sample_name,
        'stats': _make_stats(data_all),
        'figs':  _run_cond_figs(data_all, 'all'),
    })

    if data_a is not None:
        logger.info(f"Generating figures — {name_a}")
        tabs.append({
            'id':    'cond_a',
            'name':  name_a,
            'stats': _make_stats(data_a),
            'figs':  _run_cond_figs(data_a, 'cond_a'),
        })

    if data_b is not None:
        logger.info(f"Generating figures — {name_b}")
        tabs.append({
            'id':    'cond_b',
            'name':  name_b,
            'stats': _make_stats(data_b),
            'figs':  _run_cond_figs(data_b, 'cond_b'),
        })

    # ── Differential tab (condition A vs condition B) ────────────────────────
    diff_figs = None
    if has_differential:
        logger.info(f"Generating condition-comparison figures - {name_a} vs {name_b}")
        diff_figs = {
            'heatmap':   _run("diff_fig22", plot_diff_heatmap,     data_a, data_b, name_a, name_b, pval_threshold),
            'volcano':   _run("diff_fig23", plot_diff_volcano,      data_a, data_b, name_a, name_b, pval_threshold),
            'pathway':   _run("diff_fig24", plot_diff_pathway_bar,  data_a, data_b, name_a, name_b, pval_threshold),
            'stability': _run("diff_fig31", plot_lr_stability,      data_a, data_b, name_a, name_b, pval_threshold),
        }

    # ── Supplementary tables (from full dataset) ─────────────────────────────
    sig_all = data_all.significant
    top_interactions = (
        sig_all.sort_values('p-value').head(50).to_dict(orient='records')
        if len(sig_all) > 0 else []
    )
    ct_pair_rows = sorted(
        [{'sender': r, 'receiver': c, 'count': int(data_all.counts_matrix.loc[r, c])}
         for r in data_all.celltypes for c in data_all.celltypes
         if data_all.counts_matrix.loc[r, c] > 0],
        key=lambda x: -x['count'],
    )[:50]

    figure_audit_summary = _summarize_figure_audit(figure_audit)

    # ── Render HTML ──────────────────────────────────────────────────────────
    logger.info("Rendering HTML report…")
    from jinja2 import Environment, FileSystemLoader
    env  = Environment(loader=FileSystemLoader(str(Path(__file__).parent / 'templates')))
    tmpl = env.get_template('report.html')

    html = tmpl.render(
        sample_name    = sample_name,
        task_id        = task_id,
        report_date    = datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
        database_name  = os.path.basename(database_path.rstrip('/\\')),
        database_path  = str(database_path),
        version        = importlib.metadata.version("fastccc"),
        gene_sets_str  = ', '.join(gene_sets),
        pval_threshold = pval_threshold,
        top_n_lr       = top_n_lr,
        top_n_celltypes= top_n_celltypes,
        max_chords     = max_chords,
        dpi            = dpi,

        tabs             = tabs,
        has_conditions   = has_conditions,
        has_differential = has_differential,
        name_a           = name_a,
        name_b           = name_b,
        cond_a_task_id   = cond_a_task_id or "",
        cond_b_task_id   = cond_b_task_id or "",
        diff_figs        = diff_figs,
        figure_audit     = figure_audit_summary,

        top_interactions = top_interactions,
        ct_pair_counts   = ct_pair_rows,
    )

    report_path = out / 'report.html'
    report_path.write_text(html, encoding='utf-8')
    logger.success(f"Report saved to: {report_path}")
    return str(report_path)
