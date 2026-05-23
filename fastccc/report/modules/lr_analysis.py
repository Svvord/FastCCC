"""Module 2 – L-R pair analysis: dotplot, pathway classification, heatmap."""

from typing import Dict, Tuple

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import seaborn as sns

from ..loader import CCCData
from ..utils import wrap_labels
from .pathway_utils import keep_annotated_classifications


# ──────────────────────────────────────────────────────────────────────────────
# Fig 05 – top L-R pairs dot plot
# ──────────────────────────────────────────────────────────────────────────────

def plot_lr_dotplot(
    data: CCCData, top_n: int = 30, max_pairs: int = 20
) -> Tuple[plt.Figure, str]:
    """
    Dot plot: rows = top L-R pairs, columns = cell-type pairs.
    Dot size  = communication strength  (average_interactions_strength).
    Dot color = -log10(p-value).
    """
    sig = data.significant.copy()
    if sig.empty:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, 'No significant interactions found.',
                ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')
        return fig, "Dot plot (no significant interactions)."

    sig['lr_pair']   = sig['ligand'] + ' → ' + sig['receptor']
    sig['ct_pair']   = sig['sender_celltype'] + '\n→ ' + sig['receiver_celltype']

    # Choose top L-R pairs: rank by sum of -log10(p) across all ct_pairs so
    # that both frequent AND highly significant pairs score well.
    lr_score = (
        sig.groupby('lr_pair')['p-value']
           .apply(lambda s: (-np.log10(s.clip(1e-10))).sum())
           .sort_values(ascending=False)
    )
    top_lr = lr_score.head(top_n).index.tolist()
    # Choose top cell-type pairs by total interactions
    top_ct = sig['ct_pair'].value_counts().head(max_pairs).index.tolist()

    sub = sig[sig['lr_pair'].isin(top_lr) & sig['ct_pair'].isin(top_ct)]

    # Build pivot tables for size and color
    pval_piv = sub.pivot_table(
        index='lr_pair', columns='ct_pair', values='p-value', aggfunc='min'
    ).reindex(index=top_lr, columns=top_ct)

    # Merge strength values
    str_df = data.strength.copy()
    rows_list = []
    for _, row in sub.iterrows():
        ct_key = row['sender_celltype'] + '|' + row['receiver_celltype']
        lri    = row['LRI_ID']
        if ct_key in str_df.index and lri in str_df.columns:
            rows_list.append({
                'lr_pair': row['lr_pair'],
                'ct_pair': row['ct_pair'],
                'strength': float(str_df.loc[ct_key, lri]),
            })
    str_pivot = pd.DataFrame(rows_list).pivot_table(
        index='lr_pair', columns='ct_pair', values='strength', aggfunc='mean'
    ).reindex(index=top_lr, columns=top_ct)

    pval_piv  = pval_piv.fillna(1.0)
    str_pivot = str_pivot.fillna(0.0)

    neg_log_pval = -np.log10(pval_piv.values.clip(1e-10, 1))
    size_vals    = str_pivot.values
    size_vals    = (size_vals / (size_vals.max() + 1e-10)) * 200 + 5

    n_lr = len(top_lr)
    n_ct = len(top_ct)
    fig, ax = plt.subplots(figsize=(max(6, n_ct * 0.6 + 2), max(4, n_lr * 0.4 + 1.5)))

    cmap = plt.cm.Reds
    norm = mcolors.Normalize(vmin=0, vmax=neg_log_pval.max() or 1)

    for r, lr in enumerate(top_lr):
        for c, ct in enumerate(top_ct):
            pv = neg_log_pval[r, c]
            sz = size_vals[r, c]
            if pval_piv.iloc[r, c] < 0.05:
                ax.scatter(c, n_lr - 1 - r, s=sz, c=[cmap(norm(pv))],
                           edgecolors='grey', linewidths=0.3, zorder=3)

    ax.set_xticks(range(n_ct))
    ax.set_xticklabels(top_ct, rotation=60, ha='right', fontsize=7)
    ax.set_yticks(range(n_lr))
    ax.set_yticklabels([wrap_labels([lr], 35)[0] for lr in reversed(top_lr)], fontsize=7)
    ax.set_xlim(-0.5, n_ct - 0.5)
    ax.set_ylim(-0.5, n_lr - 0.5)
    ax.grid(True, lw=0.3, color='#dddddd', zorder=0)
    ax.set_axisbelow(True)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label('−log₁₀(p-value)', fontsize=8)

    # Size legend
    for s_label, s_val in [('Low CS', 20), ('Med CS', 80), ('High CS', 200)]:
        ax.scatter([], [], s=s_val, c='grey', alpha=0.6, label=s_label, edgecolors='grey', lw=0.3)
    ax.legend(title='Comm. Score', loc='upper left', bbox_to_anchor=(1.12, 1.0), fontsize=7)

    ax.set_title(f'Top Ligand–Receptor Pairs\n{data.sample_name}', pad=10)
    fig.tight_layout()

    caption = (
        f"Dot plot of the top {top_n} ligand–receptor pairs (ranked by cumulative "
        f"−log₁₀(p-value) across all cell-type pairs) across the top {max_pairs} "
        "cell-type pairs. Dot colour encodes −log₁₀(p-value) and dot size encodes the "
        "communication score (CS)."
    )
    return fig, caption


def _add_interaction_cs(sig: pd.DataFrame, data: CCCData) -> pd.DataFrame:
    sig = sig.copy()
    sig['_ct_key'] = sig['sender_celltype'] + '|' + sig['receiver_celltype']
    str_series = data.strength.stack()
    keys = list(zip(sig['_ct_key'], sig['LRI_ID']))
    sig['cs'] = pd.array([str_series.get(k, np.nan) for k in keys], dtype=float)
    return sig


def _format_pvalue(value: float, floor: float = 1e-10) -> str:
    return f"<{floor:.0e}" if value <= floor else f"{value:.2e}"


def build_celltype_lr_profiles(
    data: CCCData, top_n: int = 15
) -> Dict:
    """Build a cell-type-specific L-R evidence payload for HTML rendering."""
    sig = data.significant.copy()
    required = {
        'sender_celltype', 'receiver_celltype', 'ligand', 'receptor',
        'LRI_ID', 'p-value',
    }
    if sig.empty or not required.issubset(sig.columns):
        return {
            'profiles': [],
            'default_celltype': '',
            'caption': "Cell-type L-R evidence explorer: no significant L-R interactions available.",
        }

    evidence_floor = 1e-10
    sig = _add_interaction_cs(sig, data)
    sig['lr_pair'] = sig['ligand'].astype(str) + ' → ' + sig['receptor'].astype(str)
    if 'classification' not in sig.columns:
        sig['classification'] = 'Unannotated'
    sig['evidence'] = -np.log10(
        pd.to_numeric(sig['p-value'], errors='coerce').clip(evidence_floor, 1.0)
    )
    sig = sig[np.isfinite(sig['evidence'])].copy()
    if sig.empty:
        return {
            'profiles': [],
            'default_celltype': '',
            'caption': "Cell-type L-R evidence explorer: no finite p-value evidence available.",
        }

    def build_rows(sub: pd.DataFrame, celltype: str) -> list[dict]:
        if sub.empty:
            return []
        grouped = (
            sub.groupby(['LRI_ID', 'lr_pair'], dropna=False)
               .agg(
                   evidence=('evidence', 'sum'),
                   best_p=('p-value', 'min'),
                   mean_cs=('cs', 'mean'),
                   interactions=('LRI_ID', 'size'),
                   senders=('sender_celltype', lambda s: sorted(set(s.astype(str)))),
                   receivers=('receiver_celltype', lambda s: sorted(set(s.astype(str)))),
                   pathways=('classification', lambda s: sorted(set(
                       keep_annotated_classifications(
                           pd.DataFrame({'classification': s})
                       )['classification'].tolist()
                   ))),
               )
               .sort_values(['evidence', 'interactions'], ascending=False)
               .head(top_n)
        )
        rows = []
        for (_, lr_pair), row in grouped.iterrows():
            partners = sorted(
                (set(row['senders']) | set(row['receivers'])) - {celltype}
            )
            rows.append({
                'lr_pair': lr_pair,
                'evidence': round(float(row['evidence']), 3),
                'best_p': _format_pvalue(float(row['best_p']), evidence_floor),
                'mean_cs': round(float(row['mean_cs']), 4) if pd.notna(row['mean_cs']) else 0.0,
                'interactions': int(row['interactions']),
                'partners': ', '.join(partners[:6]) if partners else celltype,
                'pathways': ', '.join(row['pathways'][:3]) if row['pathways'] else 'Unannotated',
            })
        return rows

    profiles = []
    for celltype in sorted(data.celltypes):
        outgoing = sig[sig['sender_celltype'] == celltype]
        incoming = sig[sig['receiver_celltype'] == celltype]
        either = sig[
            (sig['sender_celltype'] == celltype)
            | (sig['receiver_celltype'] == celltype)
        ]
        if either.empty:
            continue
        profiles.append({
            'celltype': celltype,
            'n_interactions': int(len(either)),
            'n_outgoing': int(len(outgoing)),
            'n_incoming': int(len(incoming)),
            'total_evidence': round(float(either['evidence'].sum()), 3),
            'modes': {
                'Either role': build_rows(either, celltype),
                'Outgoing': build_rows(outgoing, celltype),
                'Incoming': build_rows(incoming, celltype),
            },
        })

    profiles = sorted(
        profiles,
        key=lambda profile: (-profile['total_evidence'], profile['celltype']),
    )
    return {
        'profiles': profiles,
        'default_celltype': profiles[0]['celltype'] if profiles else '',
        'caption': (
            "Interactive L-R panel. For the selected cell type and communication role, "
            "L-R pairs are ranked by cumulative -log10(p-value) across significant "
            "interactions; p-values are floored at 1e-10 for evidence scoring. "
            "Mean CS is computed over displayed significant instances."
        ),
    }


# ──────────────────────────────────────────────────────────────────────────────
# Fig 06 – pathway classification bar chart
# ──────────────────────────────────────────────────────────────────────────────

def plot_classification_bar(
    data: CCCData, top_n: int = 25
) -> Tuple[plt.Figure, str]:
    sig = data.significant.copy()
    if sig.empty or 'classification' not in sig.columns:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, 'No classification data.', ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')
        return fig, "Classification bar (no data)."

    sig = keep_annotated_classifications(sig)
    if sig.empty:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, 'No annotated pathway classifications.', ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')
        return fig, "Pathway classification bar: no annotated pathway classifications."

    counts = (
        sig.groupby('classification')
           .size()
           .sort_values(ascending=True)
           .tail(top_n)
    )

    pal = sns.color_palette('muted', len(counts))
    fig, ax = plt.subplots(figsize=(8, max(5, len(counts) * 0.38)))
    bars = ax.barh(range(len(counts)), counts.values, color=pal, edgecolor='white', lw=0.4)

    ax.set_yticks(range(len(counts)))
    ax.set_yticklabels(wrap_labels(counts.index.tolist(), 40), fontsize=8)
    ax.set_xlabel('Number of significant interactions', fontsize=9)
    ax.set_title(f'Significant Interactions by Pathway Classification\n{data.sample_name}', fontsize=10)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)

    for v, bar in zip(counts.values, bars):
        ax.text(v + counts.values.max() * 0.01, bar.get_y() + bar.get_height() / 2,
                str(v), va='center', ha='left', fontsize=7)

    fig.tight_layout()

    caption = (
        f"Horizontal bar chart showing the top {top_n} pathway classifications ranked by the "
        "number of significant L-R interactions detected. Classifications are derived from the "
        "curated interaction annotation field in the LRI database; unannotated interactions "
        "are not treated as a pathway class."
    )
    return fig, caption


def build_celltype_pathway_profiles(
    data: CCCData, top_n_pathways: int = 12
) -> Dict:
    """
    Build compact HTML payload for cell-type-specific pathway evidence.

    Pathways are ranked within each cell type by cumulative -log10 p-value over
    significant interactions involving that cell type in either communication
    role. Role counts retain sender and receiver participation separately.
    """
    sig = data.significant.copy()
    required = {'classification', 'sender_celltype', 'receiver_celltype', 'p-value'}
    if sig.empty or not required.issubset(sig.columns):
        return {
            'profiles': [],
            'default_celltype': '',
            'caption': "Cell-type pathway evidence explorer: no pathway annotations available.",
        }

    evidence_floor = 1e-10
    sig = keep_annotated_classifications(sig)
    sig['evidence'] = -np.log10(
        pd.to_numeric(sig['p-value'], errors='coerce').clip(evidence_floor, 1.0)
    )
    sig = sig[np.isfinite(sig['evidence'])].copy()
    if sig.empty:
        return {
            'profiles': [],
            'default_celltype': '',
            'caption': (
                "Cell-type pathway evidence explorer: no annotated pathway "
                "classifications available."
            ),
        }

    profiles = []
    for celltype in sorted(data.celltypes):
        involved_mask = (
            (sig['sender_celltype'] == celltype)
            | (sig['receiver_celltype'] == celltype)
        )
        involved = sig[involved_mask]
        if involved.empty:
            continue

        total = (
            involved.groupby('classification')
                    .agg(
                        evidence=('evidence', 'sum'),
                        interactions=('classification', 'size'),
                        best_p=('p-value', 'min'),
                    )
        )
        outgoing = (
            sig[sig['sender_celltype'] == celltype]
            .groupby('classification')
            .size()
        )
        incoming = (
            sig[sig['receiver_celltype'] == celltype]
            .groupby('classification')
            .size()
        )

        total['outgoing'] = outgoing.reindex(total.index, fill_value=0)
        total['incoming'] = incoming.reindex(total.index, fill_value=0)
        total = total.sort_values(
            ['evidence', 'interactions'], ascending=False
        )
        shown = total.head(top_n_pathways)

        profiles.append({
            'celltype': celltype,
            'total_evidence': round(float(total['evidence'].sum()), 3),
            'n_interactions': int(len(involved)),
            'n_pathways': int(total.shape[0]),
            'pathways': [
                {
                    'pathway': pathway,
                    'evidence': round(float(row['evidence']), 3),
                    'interactions': int(row['interactions']),
                    'outgoing': int(row['outgoing']),
                    'incoming': int(row['incoming']),
                    'best_p': _format_pvalue(float(row['best_p']), evidence_floor),
                }
                for pathway, row in shown.iterrows()
            ],
        })

    profiles = sorted(
        profiles,
        key=lambda profile: (-profile['total_evidence'], profile['celltype']),
    )
    default_celltype = profiles[0]['celltype'] if profiles else ''
    return {
        'profiles': profiles,
        'default_celltype': default_celltype,
        'caption': (
            "Interactive cell-type pathway panel. Within the selected cell type, "
            "pathway classifications are ranked by cumulative -log10(p-value) "
            "across significant interactions in which that cell type participates "
            f"as sender or receiver; p-values are floored at {evidence_floor:.0e} "
            "for this evidence score. Outgoing and incoming counts describe the "
            "cell type's communication role; autocrine interactions contribute to "
            "both role counts. Unannotated interactions are excluded from pathway "
            "ranking rather than shown as a pathway."
        ),
    }


# ──────────────────────────────────────────────────────────────────────────────
# Fig 07 – pathway × sender–receiver heatmap
# ──────────────────────────────────────────────────────────────────────────────

def plot_pathway_celltype_heatmap(
    data: CCCData, top_n_pathways: int = 20, top_n_pairs: int = 20
) -> Tuple[plt.Figure, str]:
    sig = data.significant.copy()
    if sig.empty:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "Pathway heatmap (no data)."

    sig = keep_annotated_classifications(sig)
    if sig.empty:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "Pathway heatmap: no annotated pathway classifications."

    sig['ct_pair'] = sig['sender_celltype'] + ' → ' + sig['receiver_celltype']

    top_paths = sig['classification'].value_counts().head(top_n_pathways).index.tolist()
    top_pairs = sig['ct_pair'].value_counts().head(top_n_pairs).index.tolist()

    sub = sig[sig['classification'].isin(top_paths) & sig['ct_pair'].isin(top_pairs)]
    pivot = sub.pivot_table(index='classification', columns='ct_pair', values='p-value',
                            aggfunc='count').fillna(0)
    pivot = pivot.reindex(index=top_paths, columns=top_pairs, fill_value=0)

    # Remove all-zero rows / cols
    pivot = pivot.loc[pivot.sum(axis=1) > 0, pivot.sum(axis=0) > 0]

    n_r, n_c = pivot.shape
    fig, ax = plt.subplots(figsize=(max(6, n_c * 0.55 + 2), max(4, n_r * 0.4 + 1.5)))

    sns.heatmap(
        pivot, ax=ax, cmap='YlOrRd',
        linewidths=0.3, linecolor='#eeeeee',
        cbar_kws={'label': 'Number of interactions', 'shrink': 0.6},
        xticklabels=True, yticklabels=True,
    )
    ax.set_xticklabels(
        [wrap_labels([x], 22)[0] for x in pivot.columns],
        rotation=45, ha='right', fontsize=7,
    )
    ax.set_yticklabels(
        [wrap_labels([y], 35)[0] for y in pivot.index],
        rotation=0, fontsize=7,
    )
    ax.set_xlabel('Cell-type pair (Sender → Receiver)', fontsize=9)
    ax.set_ylabel('Pathway classification', fontsize=9)
    ax.set_title(
        f'Pathway Activity Across Cell-Type Pairs\n{data.sample_name}', pad=10
    )
    fig.tight_layout()

    caption = (
        f"Heatmap of the top {top_n_pathways} pathway classifications across the top "
        f"{top_n_pairs} sender–receiver cell-type pairs. Cell values indicate the number "
        "of significant ligand–receptor interactions belonging to each annotated pathway "
        "class; unannotated interactions are excluded from this pathway view."
    )
    return fig, caption
