"""Module 2 – L-R pair analysis: dotplot, pathway classification, heatmap."""

from typing import Dict, Tuple

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import seaborn as sns

from ..loader import CCCData
from ..utils import wrap_labels


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

    # Choose top L-R pairs by frequency
    top_lr = sig['lr_pair'].value_counts().head(top_n).index.tolist()
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
        f"Dot plot of the top {top_n} most frequently detected ligand–receptor pairs across "
        f"the top {max_pairs} cell-type pairs. Dot colour encodes −log₁₀(p-value) and dot size "
        "encodes the communication score (CS)."
    )
    return fig, caption


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
        "curated interaction annotation field in the LRI database."
    )
    return fig, caption


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
        "of significant ligand–receptor interactions belonging to each pathway class."
    )
    return fig, caption
