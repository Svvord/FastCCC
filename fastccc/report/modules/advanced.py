"""Module 6 – Advanced L-R analyses: pathway info flow, small multiples, specificity, CS violin."""

from typing import Tuple

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import seaborn as sns

from ..loader import CCCData
from ..utils import wrap_labels
from .pathway_utils import keep_annotated_classifications


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

    sig = keep_annotated_classifications(sig)
    if sig.empty:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "Pathway information flow (no annotated pathway classifications)."

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


def build_celltype_pathway_flow_profiles(
    data: CCCData, top_n: int = 12
) -> dict:
    """Build cell-type-specific pathway communication-score payload."""
    sig = data.significant.copy()
    required = {'classification', 'sender_celltype', 'receiver_celltype', 'LRI_ID'}
    if sig.empty or not required.issubset(sig.columns):
        return {
            'profiles': [],
            'default_celltype': '',
            'caption': "Cell-type pathway CS explorer: no pathway annotation data available.",
        }
    sig = keep_annotated_classifications(sig)
    if sig.empty:
        return {
            'profiles': [],
            'default_celltype': '',
            'caption': "Cell-type pathway CS explorer: no annotated pathway classifications available.",
        }
    sig = _attach_cs(sig, data).dropna(subset=['cs'])
    if sig.empty:
        return {
            'profiles': [],
            'default_celltype': '',
            'caption': "Cell-type pathway CS explorer: no communication-score values available.",
        }

    def build_rows(sub: pd.DataFrame) -> list[dict]:
        if sub.empty:
            return []
        grouped = (
            sub.groupby('classification')
               .agg(
                   total_cs=('cs', 'sum'),
                   mean_cs=('cs', 'mean'),
                   interactions=('classification', 'size'),
               )
               .sort_values(['total_cs', 'interactions'], ascending=False)
               .head(top_n)
        )
        return [
            {
                'pathway': pathway,
                'total_cs': round(float(row['total_cs']), 4),
                'mean_cs': round(float(row['mean_cs']), 4),
                'interactions': int(row['interactions']),
            }
            for pathway, row in grouped.iterrows()
        ]

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
            'total_cs': round(float(either['cs'].sum()), 4),
            'n_interactions': int(len(either)),
            'modes': {
                'Either role': build_rows(either),
                'Outgoing': build_rows(outgoing),
                'Incoming': build_rows(incoming),
            },
        })

    profiles = sorted(
        profiles,
        key=lambda profile: (-profile['total_cs'], profile['celltype']),
    )
    return {
        'profiles': profiles,
        'default_celltype': profiles[0]['celltype'] if profiles else '',
        'caption': (
            "Interactive cell-type pathway communication-score panel. Pathways are "
            "ranked by summed CS among significant interactions involving the selected "
            "cell type and role; unannotated interactions are excluded from pathway ranking."
        ),
    }


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

    sig = keep_annotated_classifications(sig)
    if sig.empty:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "Per-pathway L-R multiples (no annotated pathway classifications)."

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
    data: CCCData, top_n: int = 25
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

    top_n_actual = len(order)
    caption = (
        f"Violin plots of the communication score (CS) distribution for the top {top_n_actual} "
        "sender→receiver cell-type pairs (up to {top_n} shown), ordered by median CS "
        "(left = highest). Red horizontal bars indicate the median. Individual data points "
        "are overlaid (jittered) for pairs with ≤200 interactions, giving a full picture "
        "of the score distribution shape and spread."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Fig 29 – L-R pair co-occurrence (Jaccard similarity heatmap)
# ──────────────────────────────────────────────────────────────────────────────

def plot_lr_cooccurrence(
    data: CCCData, top_n_lr: int = 30
) -> Tuple[plt.Figure, str]:
    """Jaccard similarity of active-cell-type-pair sets for top L-R pairs, clustered."""
    sig = data.significant.copy()
    if sig.empty:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "LR co-occurrence (no data)."

    sig['lr_pair'] = sig['ligand'] + ' → ' + sig['receptor']
    sig['ct_pair'] = sig['sender_celltype'] + '|' + sig['receiver_celltype']

    top_lr = sig['lr_pair'].value_counts().head(top_n_lr).index.tolist()
    sub = sig[sig['lr_pair'].isin(top_lr)]

    # Binary matrix: LR pair × cell-type pair
    pivot = (
        sub.groupby(['lr_pair', 'ct_pair'])
           .size()
           .unstack(fill_value=0)
           .reindex(index=top_lr)
           .fillna(0)
           .astype(bool)
    )

    n = len(top_lr)
    arr = pivot.values
    jaccard = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            inter = int((arr[i] & arr[j]).sum())
            union = int((arr[i] | arr[j]).sum())
            v = inter / union if union > 0 else 0.0
            jaccard[i, j] = jaccard[j, i] = v

    jac_df = pd.DataFrame(jaccard, index=top_lr, columns=top_lr)

    try:
        from scipy.cluster.hierarchy import linkage, dendrogram
        from scipy.spatial.distance import squareform
        dist = np.clip(1 - jaccard, 0, None)
        np.fill_diagonal(dist, 0)
        link = linkage(squareform(dist), method='average')
        order = dendrogram(link, no_plot=True)['leaves']
        jac_df = jac_df.iloc[order, order]
    except Exception:
        pass

    sz = max(6, n * 0.35 + 2)
    fig, ax = plt.subplots(figsize=(sz, sz * 0.9))

    sns.heatmap(
        jac_df, ax=ax, cmap='Blues', vmin=0, vmax=1,
        linewidths=0.2, linecolor='#f0f0f0',
        xticklabels=True, yticklabels=True,
        cbar_kws={'label': 'Jaccard similarity', 'shrink': 0.6},
    )
    ax.set_xticklabels(
        [wrap_labels([x], 30)[0] for x in jac_df.columns],
        rotation=60, ha='right', fontsize=6.5,
    )
    ax.set_yticklabels(
        [wrap_labels([y], 30)[0] for y in jac_df.index],
        rotation=0, fontsize=6.5,
    )
    ax.set_title(f'L-R Pair Co-occurrence (Jaccard Similarity)\n{data.sample_name}', pad=10)
    fig.tight_layout()

    caption = (
        f"Pairwise Jaccard similarity of the top {top_n_lr} L-R pairs based on the sets of "
        "cell-type pairs in which each interaction is significant. Pairs clustered together "
        "(high Jaccard, dark blue) tend to co-activate in the same communication contexts, "
        "suggesting shared regulatory programmes or pathway redundancy. Rows and columns are "
        "sorted by hierarchical clustering (average linkage on 1 − Jaccard)."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Fig 30 – Pathway crosstalk heatmap
# ──────────────────────────────────────────────────────────────────────────────

def plot_pathway_crosstalk(
    data: CCCData, top_n: int = 20
) -> Tuple[plt.Figure, str]:
    """Symmetric heatmap: number of cell-type pairs where both pathways are co-active."""
    sig = data.significant.copy()
    if sig.empty or 'classification' not in sig.columns:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "Pathway crosstalk (no data)."

    sig = keep_annotated_classifications(sig)
    if sig.empty:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "Pathway crosstalk (no annotated pathway classifications)."

    sig['ct_pair'] = sig['sender_celltype'] + '|' + sig['receiver_celltype']
    top_paths = sig['classification'].value_counts().head(top_n).index.tolist()
    sub = sig[sig['classification'].isin(top_paths)]

    ct_path = sub.groupby('ct_pair')['classification'].apply(set)

    n = len(top_paths)
    mat = np.zeros((n, n), dtype=int)
    for i, pa in enumerate(top_paths):
        for j, pb in enumerate(top_paths):
            if i > j:
                mat[i, j] = mat[j, i]
                continue
            mat[i, j] = sum(1 for s in ct_path if pa in s and pb in s)

    crosstalk = pd.DataFrame(mat, index=top_paths, columns=top_paths)

    try:
        from scipy.cluster.hierarchy import linkage, dendrogram
        from scipy.spatial.distance import pdist
        link = linkage(pdist(mat.astype(float), metric='euclidean'), method='average')
        order = dendrogram(link, no_plot=True)['leaves']
        crosstalk = crosstalk.iloc[order, order]
    except Exception:
        pass

    sz = max(6, n * 0.45 + 2)
    fig, ax = plt.subplots(figsize=(sz, sz * 0.85))

    mask = np.eye(len(crosstalk), dtype=bool)
    sns.heatmap(
        crosstalk, ax=ax, cmap='YlOrRd', mask=mask,
        linewidths=0.3, linecolor='#eeeeee',
        cbar_kws={'label': 'Shared cell-type pairs', 'shrink': 0.6},
        annot=True, fmt='d', annot_kws={'size': 6},
    )
    ax.set_xticklabels(
        [wrap_labels([x], 28)[0] for x in crosstalk.columns],
        rotation=55, ha='right', fontsize=7,
    )
    ax.set_yticklabels(
        [wrap_labels([y], 28)[0] for y in crosstalk.index],
        rotation=0, fontsize=7,
    )
    ax.set_title(f'Pathway Crosstalk — Co-active Cell-Type Pairs\n{data.sample_name}', pad=10)
    fig.tight_layout()

    caption = (
        f"Symmetric heatmap of pathway crosstalk among the top {top_n} classifications. "
        "Each cell counts the sender–receiver cell-type pairs in which both pathways have "
        "at least one significant L-R interaction simultaneously, revealing which signalling "
        "programmes tend to co-activate in the same communication context. The diagonal is "
        "masked (self-overlap). Rows/columns are sorted by hierarchical clustering."
    )
    return fig, caption
