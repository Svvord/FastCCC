"""Plotting functions for reference-based inference reports."""

from __future__ import annotations
import os
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

TREND_COLORS = {
    'Up':       '#E64B35',
    'Down':     '#4DBBD5',
    'Both Sig': '#7E6148',
    'Both NS':  '#dddddd',
}
TREND_ORDER = ['Up', 'Both Sig', 'Down', 'Both NS']


@dataclass
class InferData:
    """Container for reference-based inference results."""
    results:        pd.DataFrame          # query_infer_results.tsv
    strength:       pd.DataFrame          # query_interactions_strength.tsv
    query_name:     str      = "Query"
    reference_name: str      = "Reference"
    database_path:  str      = ""

    # derived — populated by __post_init__
    results_in_ref: pd.DataFrame = field(init=False, repr=False)

    def __post_init__(self):
        df = self.results.copy()
        df[['sender', 'receiver']] = df['sender|receiver'].str.split('|', expand=True)
        df['ct_pair'] = df['sender'] + ' → ' + df['receiver']
        df['lr_pair'] = df['ligand'] + ' → ' + df['receptor']
        df['comm_score'] = pd.to_numeric(df['comm_score'], errors='coerce')
        self.results = df
        self.results_in_ref = df[df['trend_vs_ref'].notna()].copy()

    # ── summary helpers ──────────────────────────────────────────────────────
    @property
    def trend_counts(self) -> dict:
        return {t: int((self.results_in_ref['trend_vs_ref'] == t).sum()) for t in TREND_ORDER}

    @property
    def n_tested(self) -> int:
        return len(self.results)

    @property
    def n_significant(self) -> int:
        return int(self.results['is_significant'].sum())

    @property
    def n_in_reference(self) -> int:
        return int(self.results['in_reference'].sum())

    @property
    def n_lr_pairs(self) -> int:
        return int(self.results_in_ref['lr_pair'].nunique())

    @property
    def n_ct_pairs(self) -> int:
        return int(self.results_in_ref['ct_pair'].nunique())


# ════════════════════════════════════════════════════════════════════════════
# Fig A — Global trend distribution bar
# ════════════════════════════════════════════════════════════════════════════

def plot_trend_distribution(data: InferData):
    counts = data.trend_counts
    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(
        TREND_ORDER,
        [counts[t] for t in TREND_ORDER],
        color=[TREND_COLORS[t] for t in TREND_ORDER],
        edgecolor='white', lw=0.5,
    )
    for bar, t in zip(bars, TREND_ORDER):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 10,
                f'{counts[t]:,}', ha='center', va='bottom', fontsize=9)
    ax.set_ylabel('Number of L-R interactions', fontsize=10)
    ax.set_title(
        f'{data.query_name} vs {data.reference_name} — Overall Trend Distribution', pad=10
    )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()
    caption = (
        f"<strong>Fig A — Overall Trend Distribution.</strong> "
        f"Each bar shows the number of L-R interactions in each trend category "
        f"(in-reference pairs only, n={data.n_in_reference:,}). "
        f"<em>Up</em>: significant in {data.query_name} but not {data.reference_name} (query-specific). "
        f"<em>Down</em>: not significant in {data.query_name} but significant in {data.reference_name} (reference-specific). "
        f"<em>Both Sig</em>: shared across both. <em>Both NS</em>: background."
    )
    return fig, caption


# ════════════════════════════════════════════════════════════════════════════
# Fig B — Per cell-type pair stacked bar (top 25 most differential pairs)
# ════════════════════════════════════════════════════════════════════════════

def plot_ct_pair_trend_breakdown(data: InferData, top_n: int = 25):
    df = data.results_in_ref
    ct_trend = (
        df.groupby(['ct_pair', 'trend_vs_ref'])
          .size()
          .unstack(fill_value=0)
          .reindex(columns=TREND_ORDER, fill_value=0)
    )
    ct_trend['_diff'] = ct_trend.get('Up', 0) + ct_trend.get('Down', 0)
    ct_trend = ct_trend.sort_values('_diff', ascending=True).drop(columns='_diff').tail(top_n)

    n = len(ct_trend)
    if n == 0:
        return None, "No data."

    fig, ax = plt.subplots(figsize=(10, max(5, n * 0.38 + 1.5)))
    bottoms = np.zeros(n)
    for t in TREND_ORDER:
        if t not in ct_trend.columns:
            continue
        vals = ct_trend[t].values.astype(float)
        ax.barh(range(n), vals, left=bottoms, color=TREND_COLORS[t],
                label=t, edgecolor='white', lw=0.3)
        bottoms += vals
    ax.set_yticks(range(n))
    ax.set_yticklabels(ct_trend.index.tolist(), fontsize=7)
    ax.set_xlabel('Number of L-R interactions', fontsize=9)
    ax.set_title(f'Trend Breakdown per Cell-Type Pair (top {n} most differential)', pad=10)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    ax.legend(loc='lower right', fontsize=8)
    fig.tight_layout()
    caption = (
        f"<strong>Fig B — Per Cell-Type Pair Trend Breakdown.</strong> "
        f"Stacked horizontal bars showing the number of each trend category per "
        f"sender→receiver pair. Top {n} pairs sorted by Up+Down count (most differential on top)."
    )
    return fig, caption


# ════════════════════════════════════════════════════════════════════════════
# Figs C & D — Sender-receiver heatmaps
# ════════════════════════════════════════════════════════════════════════════

def plot_trend_heatmap(data: InferData, trend: str, top_n: int = 20):
    df = data.results_in_ref
    sub = df[df['trend_vs_ref'] == trend]
    if sub.empty:
        return None, f"No {trend} interactions found."

    celltypes = sorted(set(sub['sender'].unique()) | set(sub['receiver'].unique()))
    mat = pd.DataFrame(0, index=celltypes, columns=celltypes, dtype=int)
    for _, row in sub.iterrows():
        mat.loc[row['sender'], row['receiver']] += 1

    totals = mat.sum(axis=1) + mat.sum(axis=0)
    keep = totals.sort_values(ascending=False).head(top_n).index
    mat = mat.reindex(index=keep, columns=keep).fillna(0)
    n = len(keep)

    cmap = 'Reds' if trend == 'Up' else 'Blues'
    label_map = {
        'Up':       f'Up ({data.query_name}-specific)',
        'Down':     f'Down ({data.reference_name}-specific)',
        'Both Sig': 'Both Significant',
        'Both NS':  'Both Non-Significant',
    }
    title_map = {
        'Up':   f'{data.query_name}-specific Interactions (Up)',
        'Down': f'{data.reference_name}-specific Interactions — lost in {data.query_name} (Down)',
    }
    sz = max(5, n * 0.5 + 1.5)
    fig, ax = plt.subplots(figsize=(sz + 1, sz))
    sns.heatmap(
        mat, ax=ax, cmap=cmap,
        linewidths=0.3, linecolor='#eee',
        cbar_kws={'label': f'N {label_map.get(trend, trend)} interactions', 'shrink': 0.65},
        annot=(n <= 15), fmt='d', annot_kws={'size': 6},
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=7)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=7)
    ax.set_xlabel('Receiver', fontsize=9)
    ax.set_ylabel('Sender', fontsize=9)
    ax.set_title(f'{title_map.get(trend, trend)}\n({data.query_name} query vs {data.reference_name} reference)', pad=10)
    fig.tight_layout()

    letter = 'C' if trend == 'Up' else 'D'
    caption = (
        f"<strong>Fig {letter} — Sender-Receiver Count Heatmap ({label_map.get(trend, trend)}).</strong> "
        f"Each cell shows the number of {trend} L-R interactions between a sender (row) and "
        f"receiver (column) cell type. Top {n} cell types shown by total involvement."
    )
    return fig, caption


# ════════════════════════════════════════════════════════════════════════════
# Figs E & F — Top LR pair dotplots
# ════════════════════════════════════════════════════════════════════════════

def plot_top_lr_dotplot(data: InferData, trend: str, top_n_lr: int = 25, top_n_ct: int = 20):
    df = data.results_in_ref
    sub = df[df['trend_vs_ref'] == trend].copy()
    if sub.empty:
        return None, f"No {trend} interactions to plot."

    dot_color = TREND_COLORS.get(trend, '#888')
    lr_score = sub.groupby('lr_pair')['comm_score'].sum().sort_values(ascending=False)
    top_lr = lr_score.head(top_n_lr).index.tolist()
    sub['ct_pair_nl'] = sub['sender'] + '\n→ ' + sub['receiver']
    top_ct = sub['ct_pair_nl'].value_counts().head(top_n_ct).index.tolist()

    sub2 = sub[sub['lr_pair'].isin(top_lr) & sub['ct_pair_nl'].isin(top_ct)]
    if sub2.empty:
        return None, "No data after filtering."

    cs_piv = sub2.pivot_table(
        index='lr_pair', columns='ct_pair_nl', values='comm_score', aggfunc='mean'
    ).reindex(index=top_lr, columns=top_ct, fill_value=0).fillna(0)

    n_lr, n_ct = len(top_lr), len(top_ct)
    fig, ax = plt.subplots(figsize=(max(6, n_ct * 0.65 + 2), max(5, n_lr * 0.4 + 2)))
    max_cs = cs_piv.values.max() or 1
    for r, lr in enumerate(top_lr):
        for c, ct in enumerate(top_ct):
            val = float(cs_piv.loc[lr, ct]) if (lr in cs_piv.index and ct in cs_piv.columns) else 0
            if val > 0:
                ax.scatter(c, n_lr - 1 - r,
                           s=20 + 200 * val / max_cs,
                           c=dot_color, alpha=0.75,
                           edgecolors='white', lw=0.4, zorder=3)

    ax.set_xticks(range(n_ct))
    ax.set_xticklabels(top_ct, rotation=55, ha='right', fontsize=7)
    ax.set_yticks(range(n_lr))
    ax.set_yticklabels(list(reversed(top_lr)), fontsize=7)
    ax.set_xlim(-0.5, n_ct - 0.5)
    ax.set_ylim(-0.5, n_lr - 0.5)
    ax.grid(True, lw=0.3, color='#ddd', zorder=0)
    ax.set_axisbelow(True)
    for lbl, frac in [('Low', 0.1), ('Med', 0.5), ('High', 1.0)]:
        ax.scatter([], [], s=20 + 200 * frac, c=dot_color, alpha=0.75,
                   edgecolors='white', lw=0.4, label=f'{lbl} CS')
    ax.legend(title='Comm. Score', loc='upper left', bbox_to_anchor=(1.02, 1.0), fontsize=7)

    label_map = {'Up': f'{data.query_name}-specific', 'Down': f'{data.reference_name}-specific'}
    letter = 'E' if trend == 'Up' else 'F'
    title = f'Top {label_map.get(trend, trend)} L-R Pairs ({trend})'
    ax.set_title(f'{title}\n({data.query_name} query vs {data.reference_name} reference)', pad=10)
    fig.tight_layout()
    caption = (
        f"<strong>Fig {letter} — Top {trend} L-R Pair Dotplot ({label_map.get(trend, trend)}).</strong> "
        f"Dot presence indicates a significant {trend} interaction between the L-R pair (row) "
        f"and cell-type pair (column). Dot size encodes communication score. "
        f"Top {n_lr} L-R pairs ranked by total communication score across all cell-type pairs."
    )
    return fig, caption


# ════════════════════════════════════════════════════════════════════════════
# Fig G — CS scatter (query CS vs reference significance)
# ════════════════════════════════════════════════════════════════════════════

def plot_cs_scatter(data: InferData):
    df = data.results_in_ref.dropna(subset=['comm_score', 'trend_vs_ref']).copy()
    if df.empty:
        return None, "No data."

    rng = np.random.default_rng(42)
    fig, ax = plt.subplots(figsize=(8, 6))
    for t in TREND_ORDER:
        sub = df[df['trend_vs_ref'] == t]
        if sub.empty:
            continue
        y = sub['is_significant_ref'].astype(float).fillna(-0.1)
        y = y + rng.uniform(-0.05, 0.05, len(sub))
        ax.scatter(sub['comm_score'], y,
                   c=TREND_COLORS[t], s=8, alpha=0.35, lw=0,
                   label=f'{t} (n={len(sub):,})')

    ax.axhline(0.5, color='#999', lw=0.8, ls='--')
    ax.set_xlabel(f'Communication Score ({data.query_name})', fontsize=9)
    ax.set_ylabel(f'Significant in {data.reference_name} [jittered]', fontsize=9)
    ax.set_yticks([-0.1, 0, 1])
    ax.set_yticklabels(['Not in ref', 'Not sig', 'Significant'], fontsize=8)
    ax.set_title(
        f'Query Communication Score vs Reference Significance\n'
        f'({data.query_name} query vs {data.reference_name} reference)',
        pad=10
    )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='upper right', fontsize=8, frameon=True, framealpha=0.9)
    fig.tight_layout()
    caption = (
        f"<strong>Fig G — Communication Score Scatter.</strong> "
        f"Each point is one L-R × cell-type pair. "
        f"X-axis: communication score in {data.query_name}. "
        f"Y-axis (jittered): whether the interaction was significant in the "
        f"{data.reference_name} reference (1 = significant, 0 = not significant). "
        f"Colour encodes trend category."
    )
    return fig, caption


# ════════════════════════════════════════════════════════════════════════════
# Fig H — Pathway breakdown Up vs Down
# ════════════════════════════════════════════════════════════════════════════

def plot_pathway_breakdown(data: InferData, top_n: int = 20):
    itbl_path = os.path.join(data.database_path, 'interaction_table.csv')
    if not os.path.exists(itbl_path):
        return None, "interaction_table.csv not found in database_path."

    itbl = pd.read_csv(itbl_path).set_index('id_cp_interaction')[['classification']]
    df = data.results_in_ref.merge(itbl, left_on='LRI_ID', right_index=True, how='left')
    df['classification'] = df['classification'].fillna('Unknown')

    up_counts   = df[df['trend_vs_ref'] == 'Up'  ]['classification'].value_counts()
    down_counts = df[df['trend_vs_ref'] == 'Down' ]['classification'].value_counts()
    all_paths   = up_counts.index.union(down_counts.index)
    up_counts   = up_counts.reindex(all_paths, fill_value=0)
    down_counts = down_counts.reindex(all_paths, fill_value=0)
    total       = up_counts + down_counts
    top_paths   = total.sort_values(ascending=False).head(top_n).index
    top_sorted  = (up_counts[top_paths] - down_counts[top_paths]).sort_values().index

    y  = np.arange(len(top_sorted))
    bw = 0.38
    fig, ax = plt.subplots(figsize=(9, max(5, len(top_sorted) * 0.38 + 1.5)))
    ax.barh(y + bw / 2, up_counts[top_sorted].values,   height=bw,
            color=TREND_COLORS['Up'],   label=f'Up ({data.query_name}-specific)',
            edgecolor='white', lw=0.3)
    ax.barh(y - bw / 2, down_counts[top_sorted].values, height=bw,
            color=TREND_COLORS['Down'], label=f'Down ({data.reference_name}-specific)',
            edgecolor='white', lw=0.3)
    ax.set_yticks(y)
    ax.set_yticklabels(top_sorted.tolist(), fontsize=7)
    ax.set_xlabel('Number of L-R interactions', fontsize=9)
    ax.set_title(
        f'Pathway Breakdown — Up vs Down (top {len(top_sorted)})\n'
        f'({data.query_name} query vs {data.reference_name} reference)',
        pad=10
    )
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    ax.legend(fontsize=8)
    ax.axvline(0, color='#555', lw=0.8)
    fig.tight_layout()
    caption = (
        f"<strong>Fig H — Pathway Breakdown (Up vs Down).</strong> "
        f"Grouped horizontal bars comparing the number of Up (red, {data.query_name}-specific) "
        f"and Down (blue, {data.reference_name}-specific) interactions per pathway classification. "
        f"Top {len(top_sorted)} pathways by total differential interaction count."
    )
    return fig, caption
