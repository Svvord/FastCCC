"""Module 7 - Condition-comparison summaries between two FastCCC result sets."""

from typing import Tuple

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import seaborn as sns

from ..loader import CCCData
from ..utils import wrap_labels
from .pathway_utils import keep_annotated_classifications


def _cauchy_combine(pa: float, pb: float) -> float:
    """Cauchy combination of two p-values (equal weights).

    More powerful than min(p) * 2 and consistent with FastCCC's default
    combination strategy across methods.
    """
    pa = max(float(pa), 1e-300)
    pb = max(float(pb), 1e-300)
    T = 0.5 * np.tan((0.5 - pa) * np.pi) + 0.5 * np.tan((0.5 - pb) * np.pi)
    return float(np.clip(0.5 - np.arctan(T) / np.pi, 1e-300, 1.0))


# ──────────────────────────────────────────────────────────────────────────────
# Helper: merge two conditions into a long DataFrame
# ──────────────────────────────────────────────────────────────────────────────

def _build_diff_table(
    data_a: CCCData, data_b: CCCData,
    pval_threshold: float = 0.05,
) -> pd.DataFrame:
    """
    For every (ct_pair, LRI_ID) present in either condition, return a row with:
      pval_a, pval_b, cs_a, cs_b, log2fc, sig_a, sig_b
    """
    # Union of all ct_pairs and LRI_IDs
    all_ct = data_a.pvals.index.union(data_b.pvals.index)
    all_lr = data_a.pvals.columns.union(data_b.pvals.columns)

    pvals_a = data_a.pvals.reindex(index=all_ct, columns=all_lr).fillna(1.0)
    pvals_b = data_b.pvals.reindex(index=all_ct, columns=all_lr).fillna(1.0)
    cs_a    = data_a.strength.reindex(index=all_ct, columns=all_lr).fillna(0.0)
    cs_b    = data_b.strength.reindex(index=all_ct, columns=all_lr).fillna(0.0)

    sig_a = pvals_a < pval_threshold
    sig_b = pvals_b < pval_threshold

    # Keep only pairs significant in at least one condition
    keep = sig_a | sig_b
    ct_idx, lr_idx = np.where(keep.values)

    rows = []
    for ci, li in zip(ct_idx, lr_idx):
        ct = all_ct[ci]
        lr = all_lr[li]
        pa = float(pvals_a.iat[ci, li])
        pb = float(pvals_b.iat[ci, li])
        ca = float(cs_a.iat[ci, li])
        cb = float(cs_b.iat[ci, li])
        sender, receiver = ct.split('|')
        rows.append({
            'ct_pair': ct,
            'LRI_ID': lr,
            'sender': sender,
            'receiver': receiver,
            'pval_a': pa, 'pval_b': pb,
            'cs_a': ca,  'cs_b': cb,
            'sig_a': bool(sig_a.iat[ci, li]),
            'sig_b': bool(sig_b.iat[ci, li]),
            'log2fc': np.log2((cb + 1e-6) / (ca + 1e-6)),
            'neg_log_p': -np.log10(_cauchy_combine(pa, pb)),
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # Annotate with pathway classification from data_a (or data_b)
    for src in [data_a, data_b]:
        if 'classification' in src.significant.columns:
            annot = src.significant[['LRI_ID', 'classification', 'ligand', 'receptor']].drop_duplicates('LRI_ID')
            df = df.merge(annot, on='LRI_ID', how='left')
            df['classification'] = df['classification'].fillna('Unannotated')
            break

    return df


# ──────────────────────────────────────────────────────────────────────────────
# Fig 22 - Condition-comparison interaction count heatmap
# ──────────────────────────────────────────────────────────────────────────────

def plot_diff_heatmap(
    data_a: CCCData, data_b: CCCData,
    name_a: str = 'Condition A', name_b: str = 'Condition B',
    pval_threshold: float = 0.05,
) -> Tuple[plt.Figure, str]:
    mat_a = data_a.counts_matrix
    mat_b = data_b.counts_matrix

    all_ct = sorted(set(mat_a.index.tolist()) | set(mat_b.index.tolist()))
    mat_a = mat_a.reindex(index=all_ct, columns=all_ct).fillna(0)
    mat_b = mat_b.reindex(index=all_ct, columns=all_ct).fillna(0)

    diff = mat_b - mat_a  # positive = more interactions in B

    if diff.values.sum() == 0:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, 'No condition-comparison interactions found.',
                ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')
        return fig, "Condition-comparison heatmap (no data)."

    labels = wrap_labels(all_ct)
    n = len(labels)
    sz = max(5, min(13, n * 0.55))

    abs_max = float(np.abs(diff.values).max()) or 1
    cmap = sns.diverging_palette(220, 20, as_cmap=True)

    fig, ax = plt.subplots(figsize=(sz + 2, sz))
    sns.heatmap(
        diff, ax=ax,
        cmap=cmap, center=0, vmin=-abs_max, vmax=abs_max,
        linewidths=0.3, linecolor='#eeeeee',
        annot=(n <= 15), fmt='.0f', annot_kws={'size': 6},
        cbar_kws={'label': f'Δ interactions ({name_b} − {name_a})', 'shrink': 0.65},
        xticklabels=labels, yticklabels=labels, square=True,
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
    ax.set_xlabel(f'Receiver cell type', labelpad=8)
    ax.set_ylabel(f'Sender cell type', labelpad=8)
    ax.set_title(
        f'Condition Comparison of Interaction Counts\n{name_b} - {name_a}',
        pad=10,
    )
    fig.tight_layout()

    caption = (
        f"Heatmap comparing the number of significant interactions per ordered "
        f"sender-receiver cell-type pair ({name_b} minus {name_a}). Red cells indicate "
        f"more FastCCC-positive interactions in {name_b}; blue cells indicate more in "
        f"{name_a}. Values are descriptive count differences after independent per-condition "
        f"analysis (p < {pval_threshold}); they are not p-values for a replicate-aware "
        "between-condition contrast."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Fig 23 - Condition-comparison L-R evidence plot
# ──────────────────────────────────────────────────────────────────────────────

def plot_diff_volcano(
    data_a: CCCData, data_b: CCCData,
    name_a: str = 'Condition A', name_b: str = 'Condition B',
    pval_threshold: float = 0.05,
    top_n_label: int = 20,
) -> Tuple[plt.Figure, str]:
    df = _build_diff_table(data_a, data_b, pval_threshold)

    if df.empty:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, 'No data for condition comparison.',
                ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')
        return fig, "Condition-comparison evidence plot (no data)."

    # ── Classify ──────────────────────────────────────────────────────────────
    def _classify(row):
        if row['sig_a'] and not row['sig_b']:
            return name_a
        if row['sig_b'] and not row['sig_a']:
            return name_b
        return 'Shared'

    df['status'] = df.apply(_classify, axis=1)

    # ── Y-axis cap: distinguish float-underflow zeros from real p-values ──────
    # p=0 stored as 1e-300 → neg_log_p = 300 exactly; treat ≥299 as overflow
    import math
    OVERFLOW_THRESH = 299.0
    finite_y  = df['neg_log_p'].replace([np.inf, -np.inf], np.nan).dropna()
    real_y    = finite_y[finite_y < OVERFLOW_THRESH]
    sig_floor = -np.log10(pval_threshold)

    if len(real_y) >= 10:
        # Cap just above the 99.5th percentile of real (non-underflow) values
        raw_cap = float(np.percentile(real_y, 99.5))
    else:
        raw_cap = sig_floor * 3  # fallback when almost all values are p=0

    # Round up to next multiple of 10, minimum = 2× significance threshold
    y_cap = max(math.ceil(raw_cap / 10) * 10,
                math.ceil(sig_floor * 2 / 10) * 10)

    df['y_plot']    = df['neg_log_p'].clip(upper=y_cap)
    df['is_capped'] = df['neg_log_p'] > y_cap
    n_capped = int(df['is_capped'].sum())

    color_map = {name_a: '#4DBBD5', name_b: '#E64B35', 'Shared': '#aaaaaa'}
    size_map  = {name_a: 18,        name_b: 18,        'Shared': 8}
    alpha_map = {name_a: 0.75,      name_b: 0.75,      'Shared': 0.35}

    fig, ax = plt.subplots(figsize=(9, 7))

    # ── Draw circles (non-capped) and triangles (capped) ─────────────────────
    for status, grp in df.groupby('status'):
        normal  = grp[~grp['is_capped']]
        capped  = grp[grp['is_capped']]
        col     = color_map[status]
        sz      = size_map[status]
        al      = alpha_map[status]
        zo      = 3 if status != 'Shared' else 2

        if not normal.empty:
            ax.scatter(normal['log2fc'], normal['y_plot'],
                       s=sz, c=col, alpha=al, lw=0,
                       label=f'{status} (n={len(grp):,})', zorder=zo)
        elif status in (name_a, name_b):
            # Still add legend entry even if all were capped
            ax.scatter([], [], s=sz, c=col, alpha=al, lw=0,
                       label=f'{status} (n={len(grp):,})')

        if not capped.empty:
            ax.scatter(capped['log2fc'], capped['y_plot'],
                       s=sz * 1.2, c=col, alpha=al, lw=0.5,
                       marker='^', edgecolors=col, zorder=zo + 1)

    # ── Cap line ──────────────────────────────────────────────────────────────
    if n_capped > 0:
        ax.axhline(y_cap, color='#888888', lw=0.8, ls='--', zorder=1)
        ax.text(ax.get_xlim()[1] if ax.get_xlim()[1] != 0 else 1,
                y_cap + y_cap * 0.01,
                f'cap={y_cap:.0f} ({n_capped} points ▲)',
                ha='right', va='bottom', fontsize=7, color='#888888')

    # ── Significance threshold line ───────────────────────────────────────────
    ax.axvline(0, color='#999999', lw=0.8, ls='--', zorder=1)
    p_line = -np.log10(pval_threshold)
    ax.axhline(p_line, color='#999999', lw=0.6, ls=':', zorder=1)

    # ── Label top interactions (adjustText for non-overlapping placement) ────
    # Prefer non-capped points first (spread across y-axis), then fill
    # remaining slots with capped points with extreme log2FC.
    if 'ligand' in df.columns and 'receptor' in df.columns:
        from adjustText import adjust_text

        df['label'] = df['ligand'] + ':' + df['receptor']
        cond_df = df[(df['status'] != 'Shared') & df['ligand'].notna() & df['receptor'].notna()].copy()

        # Combined score: significance × effect size → spreads labels across
        # both axes rather than clustering everything at the top of the plot.
        # Deduplicate by label so the same L-R pair isn't selected multiple
        # times from different cell-type pairs.
        cond_df = cond_df.assign(
            _score=cond_df['y_plot'] * cond_df['log2fc'].abs()
        )
        to_label = (cond_df
                    .sort_values('_score', ascending=False)
                    .drop_duplicates('label')
                    .head(top_n_label))

        texts = []
        for _, row in to_label.iterrows():
            t = ax.text(
                row['log2fc'], row['y_plot'], row['label'],
                fontsize=7, color=color_map[row['status']],
                bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.8),
            )
            texts.append(t)

        adjust_text(
            texts, ax=ax,
            expand=(1.4, 1.8),
            force_text=(0.5, 1.0),
            force_points=(0.3, 0.5),
            arrowprops=dict(arrowstyle='->', color='#aaaaaa', lw=0.6),
        )

    ax.set_xlabel(f'log₂ FC (CS: {name_b} / {name_a})', fontsize=9)
    ylabel = f'−log₁₀(Cauchy p-value)  [capped at {y_cap:.0f}]' if n_capped > 0 \
             else '−log₁₀(Cauchy combined p-value)'
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(f'L-R Condition Comparison - {name_b} vs. {name_a}', pad=10)
    ax.set_ylim(bottom=-y_cap * 0.03, top=y_cap * 1.08)
    ax.legend(loc='upper left', fontsize=8, frameon=True,
              framealpha=0.9, edgecolor='#dddddd')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()

    n_a   = int((df['status'] == name_a).sum())
    n_b   = int((df['status'] == name_b).sum())
    n_sh  = int((df['status'] == 'Shared').sum())

    cap_note = (f" The y-axis is capped at {y_cap:.0f} (99th percentile); "
                f"{n_capped:,} points with higher significance are shown as triangles (▲) at the cap line."
                if n_capped > 0 else "")

    caption = (
        f"Evidence plot comparing ligand-receptor interactions between {name_b} and "
        f"{name_a}. The x-axis shows log2 fold-change in communication score (CS). "
        f"The y-axis shows -log10 of a Cauchy-combined within-condition FastCCC p-value "
        f"used to prioritize interactions supported in either result set; it is not a "
        "formal p-value for a between-condition differential model. "
        f"Blue points (n={n_a:,}) are significant only in {name_a}; "
        f"red points (n={n_b:,}) are significant only in {name_b}; "
        f"grey points (n={n_sh:,}) are shared between conditions.{cap_note} "
        f"The top {top_n_label} condition-specific interactions by comparison evidence are labelled."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Fig 24 - Pathway comparison bar chart
# ──────────────────────────────────────────────────────────────────────────────

def plot_diff_pathway_bar(
    data_a: CCCData, data_b: CCCData,
    name_a: str = 'Condition A', name_b: str = 'Condition B',
    pval_threshold: float = 0.05,
    top_n: int = 20,
) -> Tuple[plt.Figure, str]:
    """Grouped bar: interaction count per pathway, condition A vs B."""

    def _pathway_counts(data):
        sig = data.significant
        if 'classification' not in sig.columns or sig.empty:
            return pd.Series(dtype=float)
        return keep_annotated_classifications(sig)['classification'].value_counts()

    cnt_a = _pathway_counts(data_a)
    cnt_b = _pathway_counts(data_b)

    all_paths = cnt_a.index.union(cnt_b.index)
    cnt_a = cnt_a.reindex(all_paths, fill_value=0)
    cnt_b = cnt_b.reindex(all_paths, fill_value=0)

    # Sort by total, take top N
    total = cnt_a + cnt_b
    top_paths = total.sort_values(ascending=False).head(top_n).index

    cnt_a = cnt_a[top_paths]
    cnt_b = cnt_b[top_paths]
    diff  = (cnt_b - cnt_a).sort_values()

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(14, max(5, top_n * 0.38 + 1.5)))

    # Left: grouped bar
    y = np.arange(len(top_paths))
    bw = 0.38
    top_sorted = total[top_paths].sort_values().index  # bottom-to-top order
    ax_left.barh(y - bw/2, cnt_a[top_sorted].values, height=bw,
                 color='#4DBBD5', label=name_a, edgecolor='white', lw=0.3)
    ax_left.barh(y + bw/2, cnt_b[top_sorted].values, height=bw,
                 color='#E64B35', label=name_b, edgecolor='white', lw=0.3)
    ax_left.set_yticks(y)
    ax_left.set_yticklabels(wrap_labels(top_sorted.tolist(), 38), fontsize=7)
    ax_left.set_xlabel('Number of significant interactions', fontsize=9)
    ax_left.set_title('Pathway Activity by Condition', fontsize=10, fontweight='bold')
    ax_left.spines['left'].set_visible(False)
    ax_left.tick_params(axis='y', length=0)
    ax_left.legend(fontsize=8)

    # Right: delta bar (sorted by diff)
    diff_sorted = diff.sort_values()
    colors = ['#E64B35' if v > 0 else '#4DBBD5' for v in diff_sorted.values]
    ax_right.barh(range(len(diff_sorted)), diff_sorted.values, color=colors,
                  edgecolor='white', lw=0.3)
    ax_right.set_yticks(range(len(diff_sorted)))
    ax_right.set_yticklabels(wrap_labels(diff_sorted.index.tolist(), 38), fontsize=7)
    ax_right.set_xlabel(f'Δ interactions ({name_b} − {name_a})', fontsize=9)
    ax_right.set_title('Pathway Count Difference (Delta)', fontsize=10, fontweight='bold')
    ax_right.axvline(0, color='#555', lw=0.8)
    ax_right.spines['left'].set_visible(False)
    ax_right.tick_params(axis='y', length=0)

    fig.tight_layout()

    caption = (
        f"Descriptive comparison of pathway-level interaction counts between {name_a} (blue) and "
        f"{name_b} (red). Left panel: absolute interaction counts per pathway per condition. "
        f"Right panel: difference (Δ = {name_b} − {name_a}), with red bars indicating "
        f"more FastCCC-positive interactions in {name_b} and blue bars more in {name_a}."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Fig 31 – L-R pair stability ranking
# ──────────────────────────────────────────────────────────────────────────────

def plot_lr_stability(
    data_a: CCCData, data_b: CCCData,
    name_a: str = 'Condition A', name_b: str = 'Condition B',
    pval_threshold: float = 0.05,
    top_n: int = 30,
) -> Tuple[plt.Figure, str]:
    """Stacked bar: L-R pairs ranked by n cell-type pairs significant in both conditions."""
    df = _build_diff_table(data_a, data_b, pval_threshold)
    if df.empty or 'ligand' not in df.columns:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "LR stability (no annotation data available)."

    df['lr_pair'] = df['ligand'].fillna('?') + ' → ' + df['receptor'].fillna('?')

    agg_rows = []
    for lr_pair, grp in df.groupby('lr_pair'):
        agg_rows.append({
            'lr_pair':  lr_pair,
            'n_both':   int((grp['sig_a'] & grp['sig_b']).sum()),
            'n_only_a': int((grp['sig_a'] & ~grp['sig_b']).sum()),
            'n_only_b': int((~grp['sig_a'] & grp['sig_b']).sum()),
        })

    agg = pd.DataFrame(agg_rows)
    agg = agg[agg[['n_both', 'n_only_a', 'n_only_b']].sum(axis=1) > 0]
    agg = (agg.sort_values(['n_both', 'n_only_a', 'n_only_b'], ascending=False)
              .head(top_n)
              .sort_values('n_both', ascending=True))  # bottom-to-top for barh

    if agg.empty:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "LR stability: no data after filtering."

    fig, ax = plt.subplots(figsize=(9, max(5, len(agg) * 0.38 + 1.5)))

    ax.barh(range(len(agg)), agg['n_both'].values,
            color='#3C5488', label=f'Both sig  ({name_a} & {name_b})',
            edgecolor='white', lw=0.3)
    ax.barh(range(len(agg)), agg['n_only_a'].values, left=agg['n_both'].values,
            color='#4DBBD5', label=f'Only {name_a}', edgecolor='white', lw=0.3)
    ax.barh(range(len(agg)), agg['n_only_b'].values,
            left=(agg['n_both'] + agg['n_only_a']).values,
            color='#E64B35', label=f'Only {name_b}', edgecolor='white', lw=0.3)

    ax.set_yticks(range(len(agg)))
    ax.set_yticklabels([wrap_labels([lr], 35)[0] for lr in agg['lr_pair']], fontsize=7)
    ax.set_xlabel('Number of cell-type pairs', fontsize=9)
    ax.set_title(f'L-R Pair Stability Across Conditions\n{name_a} vs. {name_b}', pad=10)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    ax.legend(loc='lower right', fontsize=8, frameon=True)
    fig.tight_layout()

    n_stable = int((agg['n_both'] > 0).sum())
    caption = (
        f"Stacked bar chart ranking the top {top_n} L-R pairs by the number of cell-type "
        f"pairs in which they are significant in both conditions (dark blue = robust). "
        f"Light blue = significant only in {name_a}; red = significant only in {name_b}. "
        f"{n_stable} L-R pairs show stable activity in at least one shared cell-type pair, "
        "representing condition-independent signalling channels."
    )
    return fig, caption
