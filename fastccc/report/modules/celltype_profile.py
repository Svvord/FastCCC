"""Module 4 – Cell type communication profiles: I/O scatter, specificity heatmap."""

from typing import Dict, Tuple

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import seaborn as sns

from ..loader import CCCData
from ..utils import wrap_labels


# ──────────────────────────────────────────────────────────────────────────────
# Fig 11 – Outgoing vs incoming interaction count scatter (bubble)
# ──────────────────────────────────────────────────────────────────────────────

def plot_io_scatter(
    data: CCCData, colors_dict: Dict
) -> Tuple[plt.Figure, str]:
    mat = data.counts_matrix
    ct  = data.celltypes

    outgoing = mat.sum(axis=1)  # row sum = sent
    incoming = mat.sum(axis=0)  # col sum = received
    total    = outgoing + incoming

    max_total = total.max() if total.max() > 0 else 1
    sizes = (total / max_total * 500 + 30).values

    fig, ax = plt.subplots(figsize=(7, 6))

    for i, c in enumerate(ct):
        ax.scatter(
            outgoing[c], incoming[c],
            s=sizes[i],
            color=colors_dict.get(c, '#999999'),
            alpha=0.85, edgecolors='white', lw=0.5, zorder=3,
        )
        ax.annotate(
            c, (outgoing[c], incoming[c]),
            textcoords='offset points', xytext=(5, 3),
            fontsize=7, ha='left',
        )

    # Diagonal reference line
    max_val = max(outgoing.max(), incoming.max()) * 1.1
    ax.plot([0, max_val], [0, max_val], 'k--', lw=0.8, alpha=0.4, zorder=1)

    ax.set_xlabel('Number of outgoing interactions (sender)', fontsize=9)
    ax.set_ylabel('Number of incoming interactions (receiver)', fontsize=9)
    ax.set_title(f'Sender vs. Receiver Profile\n{data.sample_name}', pad=10)
    ax.set_xlim(left=-max_val * 0.03)
    ax.set_ylim(bottom=-max_val * 0.03)

    # Quadrant annotations
    ax.text(0.97, 0.03, 'Strong sender\nweak receiver', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=7, color='grey', style='italic')
    ax.text(0.03, 0.97, 'Weak sender\nstrong receiver', transform=ax.transAxes,
            ha='left', va='top', fontsize=7, color='grey', style='italic')

    fig.tight_layout()

    caption = (
        "Scatter plot comparing the total number of significant outgoing (x-axis) and "
        "incoming (y-axis) ligand–receptor interactions per cell type. Bubble size reflects "
        "the total interaction count. The dashed diagonal indicates equal sender and receiver "
        "activity. Cell types above the diagonal are predominantly receivers; those below "
        "are predominantly senders."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Fig 12 – Pathway specificity heatmap per sender cell type
# ──────────────────────────────────────────────────────────────────────────────

def plot_sender_pathway_heatmap(
    data: CCCData, top_n_pathways: int = 20, top_n_ct: int = 15
) -> Tuple[plt.Figure, str]:
    sig = data.significant.copy()
    if sig.empty:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "Sender pathway heatmap (no data)."

    top_ct   = sig['sender_celltype'].value_counts().head(top_n_ct).index.tolist()
    top_path = sig['classification'].value_counts().head(top_n_pathways).index.tolist()

    sub = sig[sig['sender_celltype'].isin(top_ct) & sig['classification'].isin(top_path)]
    pivot = sub.pivot_table(
        index='classification', columns='sender_celltype',
        values='p-value', aggfunc='count',
    ).reindex(index=top_path, columns=top_ct, fill_value=0)
    pivot = pivot.loc[pivot.sum(axis=1) > 0, pivot.sum(axis=0) > 0]

    # Normalise row-wise to show relative enrichment per sender
    row_max = pivot.max(axis=1).replace(0, 1)
    pivot_norm = pivot.div(row_max, axis=0)

    n_r, n_c = pivot_norm.shape
    fig, ax = plt.subplots(figsize=(max(6, n_c * 0.6 + 2.5), max(4, n_r * 0.4 + 1.5)))

    sns.heatmap(
        pivot_norm, ax=ax,
        cmap='Purples',
        linewidths=0.3, linecolor='#eeeeee',
        cbar_kws={'label': 'Row-normalised interaction count', 'shrink': 0.6},
        xticklabels=True, yticklabels=True,
        annot=pivot.values if n_r * n_c <= 200 else False,
        fmt='.0f', annot_kws={'size': 6},
    )
    ax.set_xticklabels(
        [wrap_labels([x], 22)[0] for x in pivot_norm.columns],
        rotation=45, ha='right', fontsize=8,
    )
    ax.set_yticklabels(
        [wrap_labels([y], 40)[0] for y in pivot_norm.index],
        rotation=0, fontsize=8,
    )
    ax.set_xlabel('Sender cell type', fontsize=9)
    ax.set_ylabel('Pathway classification', fontsize=9)
    ax.set_title(f'Sender Cell Type Pathway Specificity\n{data.sample_name}', pad=10)
    fig.tight_layout()

    caption = (
        f"Row-normalised heatmap of the top {top_n_pathways} pathway classifications across "
        f"the top {top_n_ct} sender cell types. Values are normalised per pathway (row) so "
        "that the cell type with the highest interaction count for a given pathway scores 1.0, "
        "highlighting relative pathway specialisation among senders. Raw counts are shown as "
        "annotations where the matrix is small enough."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Fig 13 – Differential interaction network (compare two groups of cell types)
# Useful as a standalone panel showing directed signalling dominance
# ──────────────────────────────────────────────────────────────────────────────

def plot_interaction_flow(
    data: CCCData, colors_dict: Dict, top_n: int = 10
) -> Tuple[plt.Figure, str]:
    """
    Stacked bar showing outgoing interaction counts per sender,
    broken down by top receiver cell types. Gives a 'flow' perspective.
    """
    mat = data.counts_matrix
    ct  = data.celltypes

    # Top senders
    top_senders = mat.sum(axis=1).sort_values(ascending=False).head(top_n).index.tolist()

    # For each sender, breakdown by receiver
    sub_mat = mat.loc[top_senders]

    # Top receivers across those senders
    top_receivers = sub_mat.sum(axis=0).sort_values(ascending=False).head(top_n).index.tolist()
    sub_mat = sub_mat[top_receivers]

    fig, ax = plt.subplots(figsize=(10, max(4, len(top_senders) * 0.5 + 1.5)))

    bottoms = np.zeros(len(top_senders))
    for rec in top_receivers:
        vals = sub_mat[rec].values.astype(float)
        bars = ax.barh(
            range(len(top_senders)), vals, left=bottoms,
            color=colors_dict.get(rec, '#aaaaaa'),
            label=rec, edgecolor='white', lw=0.4,
        )
        bottoms += vals

    ax.set_yticks(range(len(top_senders)))
    ax.set_yticklabels(wrap_labels(top_senders), fontsize=8)
    ax.set_xlabel('Number of significant outgoing interactions', fontsize=9)
    ax.set_title(f'Outgoing Interaction Flow — Top Senders\n{data.sample_name}', pad=10)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    ax.legend(
        title='Receiver', loc='upper right',
        bbox_to_anchor=(1.01, 1), fontsize=7, ncol=1 + len(top_receivers) // 12,
    )
    fig.tight_layout()

    caption = (
        f"Stacked horizontal bar chart showing the outgoing interaction counts for the top "
        f"{top_n} sender cell types, broken down by receiver cell type. Each colour segment "
        "corresponds to a distinct receiver, revealing the directionality and diversity of "
        "cell-cell communication from each major sender."
    )
    return fig, caption
