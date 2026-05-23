"""Module 3 – Pathway enrichment: ORA via gseapy, TF prediction heatmap."""

import warnings
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ..loader import CCCData
from ..utils import wrap_labels

try:
    import gseapy as gp
    _GSEAPY_OK = True
except ImportError:
    _GSEAPY_OK = False


# ──────────────────────────────────────────────────────────────────────────────
# Helper: run ORA
# ──────────────────────────────────────────────────────────────────────────────

def _run_ora(
    gene_list: List[str],
    background: List[str],
    gene_sets: List[str],
    organism: str = 'human',
) -> Tuple[Optional[pd.DataFrame], str]:
    if not _GSEAPY_OK:
        return None, "ORA skipped: gseapy is not installed."
    if not gene_list:
        return None, "ORA skipped: no eligible genes were detected."
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            enr = gp.enrichr(
                gene_list=list(gene_list),
                background=list(background),
                gene_sets=gene_sets,
                organism=organism,
                outdir=None,
                verbose=False,
            )
        res = enr.results
        if res is None or res.empty:
            return None, "ORA skipped: Enrichr returned no enrichment terms."
        res = res[res['Adjusted P-value'] < 0.05].copy()
        if res.empty:
            return None, "ORA skipped: no terms passed adjusted p-value < 0.05."
        res['-log10(Padj)'] = -np.log10(res['Adjusted P-value'].clip(1e-10))
        return res.sort_values('-log10(Padj)', ascending=False), ""
    except Exception as e:
        warnings.warn(f"gseapy ORA failed: {e}")
        return None, f"ORA skipped: Enrichr request failed ({e})."


def _plot_ora_bar(
    enr_result: pd.DataFrame,
    title: str,
    top_n: int = 20,
    color: str = '#d9534f',
) -> Tuple[plt.Figure, str]:
    sub = enr_result.head(top_n).copy()
    sub = sub.sort_values('-log10(Padj)', ascending=True)

    fig, ax = plt.subplots(figsize=(9, max(4, len(sub) * 0.38)))
    bars = ax.barh(range(len(sub)), sub['-log10(Padj)'].values, color=color, edgecolor='white', lw=0.4)

    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels(wrap_labels(sub['Term'].tolist(), 50), fontsize=7)
    ax.axvline(x=-np.log10(0.05), color='black', lw=0.8, linestyle='--', alpha=0.6)
    ax.set_xlabel('−log₁₀(Adjusted P-value)', fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)

    for v, bar in zip(sub['-log10(Padj)'].values, bars):
        ax.text(v + 0.05, bar.get_y() + bar.get_height() / 2,
                f'{v:.2f}', va='center', ha='left', fontsize=7)

    fig.tight_layout()
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# Fig 08 – ligand gene ORA
# ──────────────────────────────────────────────────────────────────────────────

def plot_ligand_ora(
    data: CCCData,
    gene_sets: List[str] = ('KEGG_2021_Human', 'GO_Biological_Process_2023'),
    top_n: int = 20,
) -> Tuple[plt.Figure, str]:
    sig = data.significant
    sig_ligands = _explode_genes(sig['ligand'].dropna().unique().tolist())
    background  = _explode_genes(_get_tested_genes(data, role='ligand'))

    enr, ora_note = _run_ora(sig_ligands, background, list(gene_sets))

    if enr is None or enr.empty:
        return None, f"Ligand gene {ora_note}"

    fig = _plot_ora_bar(enr, f'Ligand Gene Enrichment (ORA)\n{data.sample_name}',
                        top_n=top_n, color='#e07b54')

    caption = (
        f"Over-representation analysis (ORA) of the {len(sig_ligands)} unique ligand genes "
        f"detected in significant interactions, tested against {', '.join(gene_sets)} gene sets "
        f"using expression-filtered database ligands as background (n={len(background)}). "
        "The dashed line marks the significance threshold (p_adj = 0.05)."
    )
    return fig, caption


def build_celltype_gene_profiles(
    data: CCCData, role: str = 'ligand', top_n: int = 20
) -> Dict:
    """Build cell-type-specific ligand or receptor gene evidence payload."""
    if role not in {'ligand', 'receptor'}:
        raise ValueError("role must be 'ligand' or 'receptor'")

    sig = data.significant.copy()
    gene_col = 'ligand' if role == 'ligand' else 'receptor'
    ct_col = 'sender_celltype' if role == 'ligand' else 'receiver_celltype'
    partner_col = 'receiver_celltype' if role == 'ligand' else 'sender_celltype'
    required = {gene_col, ct_col, partner_col, 'p-value'}
    if sig.empty or not required.issubset(sig.columns):
        return {
            'profiles': [],
            'default_celltype': '',
            'role': role,
            'caption': f"Cell-type {role} gene explorer: no significant interactions available.",
        }

    evidence_floor = 1e-10
    sig['evidence'] = -np.log10(
        pd.to_numeric(sig['p-value'], errors='coerce').clip(evidence_floor, 1.0)
    )
    sig = sig[np.isfinite(sig['evidence'])].copy()

    rows = []
    for _, row in sig.iterrows():
        for gene in _explode_genes([row[gene_col]]):
            rows.append({
                'celltype': row[ct_col],
                'partner': row[partner_col],
                'gene': gene,
                'evidence': row['evidence'],
                'p_value': row['p-value'],
            })
    if not rows:
        return {
            'profiles': [],
            'default_celltype': '',
            'role': role,
            'caption': f"Cell-type {role} gene explorer: no eligible genes available.",
        }

    gene_df = pd.DataFrame(rows)
    profiles = []
    for celltype in sorted(gene_df['celltype'].dropna().astype(str).unique()):
        sub = gene_df[gene_df['celltype'] == celltype]
        grouped = (
            sub.groupby('gene')
               .agg(
                   evidence=('evidence', 'sum'),
                   interactions=('gene', 'size'),
                   best_p=('p_value', 'min'),
                   partners=('partner', lambda s: sorted(set(s.astype(str)))),
               )
               .sort_values(['evidence', 'interactions'], ascending=False)
               .head(top_n)
        )
        profiles.append({
            'celltype': celltype,
            'n_genes': int(sub['gene'].nunique()),
            'n_interactions': int(len(sub)),
            'total_evidence': round(float(sub['evidence'].sum()), 3),
            'genes': [
                {
                    'gene': gene,
                    'evidence': round(float(row['evidence']), 3),
                    'interactions': int(row['interactions']),
                    'best_p': f"<{evidence_floor:.0e}" if float(row['best_p']) <= evidence_floor else f"{float(row['best_p']):.2e}",
                    'partners': ', '.join(row['partners'][:6]),
                }
                for gene, row in grouped.iterrows()
            ],
        })

    profiles = sorted(
        profiles,
        key=lambda profile: (-profile['total_evidence'], profile['celltype']),
    )
    return {
        'profiles': profiles,
        'default_celltype': profiles[0]['celltype'] if profiles else '',
        'role': role,
        'caption': (
            f"Interactive cell-type {role} gene panel. Genes are ranked within the "
            "selected cell type by cumulative -log10(p-value) over significant "
            "interactions. This is a cell-type evidence view, not a separate ORA test."
        ),
    }


# ──────────────────────────────────────────────────────────────────────────────
# Fig 09 – receptor gene ORA
# ──────────────────────────────────────────────────────────────────────────────

def plot_receptor_ora(
    data: CCCData,
    gene_sets: List[str] = ('KEGG_2021_Human', 'GO_Biological_Process_2023'),
    top_n: int = 20,
) -> Tuple[plt.Figure, str]:
    sig = data.significant
    sig_receptors = _explode_genes(sig['receptor'].dropna().unique().tolist())
    background    = _explode_genes(_get_tested_genes(data, role='receptor'))

    enr, ora_note = _run_ora(sig_receptors, background, list(gene_sets))

    if enr is None or enr.empty:
        return None, f"Receptor gene {ora_note}"

    fig = _plot_ora_bar(enr, f'Receptor Gene Enrichment (ORA)\n{data.sample_name}',
                        top_n=top_n, color='#5b8db8')

    caption = (
        f"Over-representation analysis (ORA) of the {len(sig_receptors)} unique receptor genes "
        f"detected in significant interactions, tested against {', '.join(gene_sets)} gene sets "
        f"using expression-filtered database receptors as background (n={len(background)}). "
        "The dashed line marks the significance threshold (p_adj = 0.05)."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Fig 10 – Receiver cell type × predicted TF heatmap
# ──────────────────────────────────────────────────────────────────────────────

def plot_tf_heatmap(
    data: CCCData, top_n_tf: int = 30, top_n_ct: int = 15
) -> Tuple[Optional[plt.Figure], str]:
    if data.receptor_tf is None:
        return None, "TF heatmap skipped: receptor_to_transcription_factor.csv not available for this database."

    sig = data.significant.copy()
    rtf = data.receptor_tf.copy()

    # Explode complex receptors (e.g. "ITGAV,ITGB3" → individual genes)
    sig['receptor_genes'] = sig['receptor'].apply(lambda x: [g.strip() for g in str(x).split(',')])

    rows = []
    for _, row in sig.iterrows():
        receiver = row['receiver_celltype']
        for rg in row['receptor_genes']:
            matched = rtf[rtf['Receptor'].str.contains(rg, case=False, na=False)]
            for tf in matched['TF'].tolist():
                rows.append({'receiver': receiver, 'TF': tf, 'p_value': row['p-value']})

    if not rows:
        return None, "TF heatmap: no receptor–TF links matched in significant interactions."

    tf_df = pd.DataFrame(rows)

    # Score = -log10(mean p-value) per receiver × TF
    pivot = tf_df.pivot_table(
        index='TF', columns='receiver', values='p_value',
        aggfunc=lambda x: -np.log10(np.mean(x).clip(1e-10)),
    )

    # Select top TFs and cell types
    top_tf = pivot.mean(axis=1).sort_values(ascending=False).head(top_n_tf).index
    top_ct = pivot.mean(axis=0).sort_values(ascending=False).head(top_n_ct).index
    pivot  = pivot.reindex(index=top_tf, columns=top_ct).fillna(0)

    n_r, n_c = pivot.shape
    fig, ax = plt.subplots(figsize=(max(6, n_c * 0.55 + 2), max(5, n_r * 0.35 + 1.5)))

    sns.heatmap(
        pivot, ax=ax, cmap='RdYlBu_r',
        linewidths=0.2, linecolor='#eeeeee',
        cbar_kws={'label': '−log₁₀(mean p-value)', 'shrink': 0.6},
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
    ax.set_xlabel('Receiver cell type', fontsize=9)
    ax.set_ylabel('Transcription factor', fontsize=9)
    ax.set_title(f'Predicted TF Activation by Receiver Cell Type\n{data.sample_name}', pad=10)
    fig.tight_layout()

    caption = (
        f"Heatmap of predicted downstream transcription factor (TF) activation in each receiver "
        "cell type, inferred from significant receptor genes and the receptor→TF mapping in "
        "the LRI database. Colour intensity reflects −log₁₀(mean p-value) of the upstream "
        f"receptor interactions, showing the top {top_n_tf} TFs and top {top_n_ct} receiver "
        "cell types."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

def _load_db_tables(database_path: str):
    """Return (itbl, id2sym, comp_map) for the given database path."""
    import os
    itbl = pd.read_csv(os.path.join(database_path, 'interaction_table.csv'))
    gtbl = pd.read_csv(os.path.join(database_path, 'gene_table.csv'))
    id2sym = gtbl.set_index('protein_id')['hgnc_symbol'].to_dict()
    comp_comp = pd.read_csv(os.path.join(database_path, 'complex_composition_table.csv'))
    comp_map: dict = comp_comp.groupby('complex_multidata_id')['protein_multidata_id'].apply(list).to_dict()
    return itbl, id2sym, comp_map


def _resolve_genes(itbl: pd.DataFrame, id2sym: dict, comp_map: dict, role: str) -> List[str]:
    col = 'multidata_1_id' if role == 'ligand' else 'multidata_2_id'

    def resolve(mid):
        if mid in id2sym:
            return [id2sym[mid]]
        if mid in comp_map:
            return [id2sym[p] for p in comp_map[mid] if p in id2sym]
        return []

    genes = []
    for mid in itbl[col].unique():
        genes.extend(resolve(mid))
    return list(set(genes))


def _get_all_genes(database_path: str, role: str) -> List[str]:
    """Extract all ligand or receptor gene symbols from the database."""
    itbl, id2sym, comp_map = _load_db_tables(database_path)
    return _resolve_genes(itbl, id2sym, comp_map, role)


def _get_tested_genes(data, role: str) -> List[str]:
    """Return gene symbols for LRI interactions actually tested in the analysis.

    Uses data.pvals.columns (the set of tested LRI_IDs) as a filter so the ORA
    background reflects the expression-filtered universe rather than the full DB.
    """
    itbl, id2sym, comp_map = _load_db_tables(data.database_path)
    tested_ids = set(data.pvals.columns)
    itbl_tested = itbl[itbl['id_cp_interaction'].isin(tested_ids)]
    return _resolve_genes(itbl_tested, id2sym, comp_map, role)


def _explode_genes(gene_list: List[str]) -> List[str]:
    """Flatten complex gene names like 'ITGAV,ITGB3' into individual symbols."""
    result = []
    for g in gene_list:
        for sub in str(g).split(','):
            s = sub.strip()
            if s and s != 'nan':
                result.append(s)
    return list(set(result))
