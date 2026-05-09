"""Module 6 – Advanced L-R analyses: pathway info flow, small multiples, specificity, CS violin."""

from typing import Tuple

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import seaborn as sns

from ..loader import CCCData
from ..utils import wrap_labels


# ──────────────────────────────────────────────────────────────────────────────
# Helper: attach CS values to significant interactions
# ──────────────────────────────────────────────────────────────────────────────

def _attach_cs(sig: pd.DataFrame, data: CCCData) -> pd.DataFrame:
    sig = sig.copy()
    sig['_ct_pair'] = sig['sender_celltype'] + '|' + sig['receiver_celltype']
    str_series = data.strength.stack()
    keys = list(zip(sig['_ct_pair'], sig['LRI_ID']))
    sig['cs'] = pd.array([str_series.get(k, np.nan) for k in keys], dtype=float)
    return sig


# ──────────────────────────────────────────────────────────────────────────────
# Fig 18 – Pathway information flow
# ──────────────────────────────────────────────────────────────────────────────

def plot_pathway_info_flow(
    data: CCCData, top_n: int = 25
) -> Tuple[plt.Figure, str]:
    sig = data.significant.copy()
    if sig.empty or 'classification' not in sig.columns:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "Pathway information flow (no data)."

    sig = _attach_cs(sig, data)

    flow = (
        sig.groupby('classification')['cs']
           .sum()
           .replace(0, np.nan)
           .dropna()
           .sort_values(ascending=True)
           .tail(top_n)
    )
    if flow.empty:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "Pathway information flow (no CS data)."

    pal = sns.color_palette('tab20', len(flow))
    fig, ax = plt.subplots(figsize=(8, max(5, len(flow) * 0.38 + 1.5)))

    bars = ax.barh(range(len(flow)), flow.values, color=pal, edgecolor='white', lw=0.4)
    ax.set_yticks(range(len(flow)))
    ax.set_yticklabels(wrap_labels(flow.index.tolist(), 40), fontsize=8)
    ax.set_xlabel('Total communication score (summed CS)', fontsize=9)
    ax.set_title(f'Pathway Information Flow\n{data.sample_name}', pad=10)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)

    for v, bar in zip(flow.values, bars):
        ax.text(v + flow.values.max() * 0.01, bar.get_y() + bar.get_height() / 2,
                f'{v:.1f}', va='center', ha='left', fontsize=7)

    fig.tight_layout()

    caption = (
        f"Bar chart of the top {top_n} pathway classifications ranked by total communication "
        "score (sum of CS across all significant L-R interactions per pathway), highlighting "
        "which pathways carry the strongest cumulative signalling weight in this sample."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Fig 19 – Per-pathway top L-R pairs (small multiples)
# ──────────────────────────────────────────────────────────────────────────────

def plot_pathway_lr_multiples(
    data: CCCData,
    top_n_pathways: int = 9,
    top_n_lr: int = 8,
) -> Tuple[plt.Figure, str]:
    sig = data.significant.copy()
    if sig.empty or 'classification' not in sig.columns:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "Per-pathway L-R multiples (no data)."

    sig = _attach_cs(sig, data)
    sig['lr_pair'] = sig['ligand'] + ' → ' + sig['receptor']

    top_paths = sig['classification'].value_counts().head(top_n_pathways).index.tolist()

    n_cols = 3
    n_rows = int(np.ceil(len(top_paths) / n_cols))
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(n_cols * 4.5, n_rows * 3.2),
        constrained_layout=True,
    )
    axes_flat = np.array(axes).flatten()

    for idx, path in enumerate(top_paths):
        ax = axes_flat[idx]
        sub = sig[sig['classification'] == path]

        lr_cs = (
            sub.groupby('lr_pair')['cs']
               .sum()
               .sort_values(ascending=True)
               .tail(top_n_lr)
        )

        if lr_cs.empty:
            ax.axis('off')
            continue

        pal = sns.color_palette('husl', len(lr_cs))
        ax.barh(range(len(lr_cs)), lr_cs.values, color=pal, edgecolor='white', lw=0.3)
        ax.set_yticks(range(len(lr_cs)))
        ax.set_yticklabels(wrap_labels(lr_cs.index.tolist(), 28), fontsize=7)
        ax.set_title(wrap_labels([path], 30)[0], fontsize=8, fontweight='bold', pad=4)
        ax.set_xlabel('Total CS', fontsize=7)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='y', length=0, labelsize=7)
        ax.tick_params(axis='x', labelsize=7)

    # Hide unused subplots
    for idx in range(len(top_paths), len(axes_flat)):
        axes_flat[idx].axis('off')

    fig.suptitle(
        f'Top Ligand–Receptor Pairs per Pathway\n{data.sample_name}',
        fontsize=11, fontweight='bold', y=1.01,
    )

    caption = (
        f"Small-multiples bar charts showing the top {top_n_lr} ligand–receptor pairs "
        f"(by total communication score) within each of the top {top_n_pathways} pathway "
        "classifications. Each panel represents one pathway, enabling rapid identification "
        "of the dominant L-R pairs driving a specific signalling programme."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Fig 20 – L-R pair specificity heatmap
# ──────────────────────────────────────────────────────────────────────────────

def plot_lr_specificity(
    data: CCCData,
    top_n_lr: int = 30,
    top_n_pairs: int = 20,
) -> Tuple[plt.Figure, str]:
    sig = data.significant.copy()
    if sig.empty:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "L-R specificity heatmap (no data)."

    sig = _attach_cs(sig, data)
    sig['lr_pair'] = sig['ligand'] + ' → ' + sig['receptor']
    sig['ct_pair'] = sig['sender_celltype'] + ' → ' + sig['receiver_celltype']

    top_lr = sig['lr_pair'].value_counts().head(top_n_lr).index.tolist()
    top_ct = sig['ct_pair'].value_counts().head(top_n_pairs).index.tolist()

    sub = sig[sig['lr_pair'].isin(top_lr) & sig['ct_pair'].isin(top_ct)]
    pivot = sub.pivot_table(
        index='lr_pair', columns='ct_pair', values='cs', aggfunc='sum',
    ).reindex(index=top_lr, columns=top_ct, fill_value=0)

    # Remove all-zero rows
    pivot = pivot.loc[pivot.sum(axis=1) > 0]
    if pivot.empty:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "L-R specificity heatmap (empty pivot)."

    # Row-normalise to show relative specificity
    row_max = pivot.max(axis=1).replace(0, 1)
    pivot_norm = pivot.div(row_max, axis=0)

    n_r, n_c = pivot_norm.shape
    fig, ax = plt.subplots(figsize=(max(6, n_c * 0.55 + 2.5), max(5, n_r * 0.33 + 2)))

    sns.heatmap(
        pivot_norm, ax=ax,
        cmap='RdPu',
        linewidths=0.2, linecolor='#f0f0f0',
        cbar_kws={'label': 'Row-normalised CS', 'shrink': 0.55},
        xticklabels=True, yticklabels=True,
    )
    ax.set_xticklabels(
        [wrap_labels([x], 22)[0] for x in pivot_norm.columns],
        rotation=45, ha='right', fontsize=7,
    )
    ax.set_yticklabels(
        [wrap_labels([y], 38)[0] for y in pivot_norm.index],
        rotation=0, fontsize=7,
    )
    ax.set_xlabel('Cell-type pair (Sender → Receiver)', fontsize=9)
    ax.set_ylabel('L-R pair', fontsize=9)
    ax.set_title(
        f'L-R Pair Specificity across Cell-Type Pairs\n{data.sample_name}', pad=10
    )
    fig.tight_layout()

    caption = (
        f"Row-normalised heatmap of the top {top_n_lr} L-R pairs across the top "
        f"{top_n_pairs} cell-type pairs. Values are normalised per L-R pair so that the "
        "cell-type pair with the highest total CS scores 1.0, revealing whether an L-R "
        "interaction is broadly used (high scores in many columns) or cell-type-specific "
        "(concentrated in one or few columns)."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Fig 21 – CS distribution violin plot
# ──────────────────────────────────────────────────────────────────────────────

def plot_cs_violin(
    data: CCCData, top_n: int = 12
) -> Tuple[plt.Figure, str]:
    sig = data.significant.copy()
    if sig.empty:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "CS violin plot (no data)."

    sig = _attach_cs(sig, data)
    sig = sig.dropna(subset=['cs'])
    if sig.empty:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "CS violin plot (no CS values)."

    # Group by top cell-type pairs
    sig['ct_pair'] = sig['sender_celltype'] + '\n→\n' + sig['receiver_celltype']
    top_pairs = sig['ct_pair'].value_counts().head(top_n).index.tolist()
    sub = sig[sig['ct_pair'].isin(top_pairs)]
    order = sub.groupby('ct_pair')['cs'].median().sort_values(ascending=False).index.tolist()

    fig, ax = plt.subplots(figsize=(max(8, len(order) * 0.9 + 2), 5))

    cs_groups = [sub[sub['ct_pair'] == p]['cs'].values for p in order]
    # Violin
    parts = ax.violinplot(
        cs_groups, positions=range(len(order)),
        showmedians=True, showextrema=False,
        widths=0.7,
    )
    for pc in parts['bodies']:
        pc.set_facecolor('#3C5488')
        pc.set_alpha(0.55)
    parts['cmedians'].set_color('#E64B35')
    parts['cmedians'].set_lw(1.5)

    # Overlay swarm-like jitter for small groups
    for i, vals in enumerate(cs_groups):
        if len(vals) <= 200:
            jitter = np.random.default_rng(seed=i).uniform(-0.12, 0.12, len(vals))
            ax.scatter(i + jitter, vals, s=4, alpha=0.4, color='#666666', zorder=3)

    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(
        [wrap_labels([p.replace('\n→\n', ' → '), ], 26)[0].replace(' → ', '\n→ ') for p in order],
        fontsize=7, ha='center',
    )
    ax.set_ylabel('Communication Score (CS)', fontsize=9)
    ax.set_title(
        f'CS Distribution per Cell-Type Pair\n{data.sample_name}', pad=10
    )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(-0.6, len(order) - 0.4)
    fig.tight_layout()

    caption = (
        f"Violin plots of the communication score (CS) distribution for the top {top_n} "
        "sender→receiver cell-type pairs, ordered by median CS (left = highest). "
        "Red horizontal bars indicate the median. Individual data points are overlaid "
        "(jittered) for pairs with ≤200 interactions, giving a full picture of the "
        "score distribution shape and spread."
    )
    return fig, caption
