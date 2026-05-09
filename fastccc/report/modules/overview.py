"""Module 1 – Global CCC overview: heatmaps, chord diagram, sender/receiver bars."""

from typing import Dict, List, Tuple

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.path import Path

from ..loader import CCCData
from ..utils import sequential_cmap, wrap_labels


# ──────────────────────────────────────────────────────────────────────────────
# Fig 01 – interaction count heatmap
# ──────────────────────────────────────────────────────────────────────────────

def plot_count_heatmap(data: CCCData) -> Tuple[plt.Figure, str]:
    mat = data.counts_matrix
    labels = wrap_labels(mat.index.tolist())

    n = len(labels)
    sz = max(5, min(12, n * 0.55))
    fig, ax = plt.subplots(figsize=(sz + 1.5, sz))

    mask = mat.values == 0
    sns.heatmap(
        mat, ax=ax, mask=mask,
        cmap=sequential_cmap(), linewidths=0.3, linecolor='#eeeeee',
        annot=(n <= 15), fmt='d', annot_kws={'size': 7},
        cbar_kws={'label': '# significant interactions', 'shrink': 0.7},
        xticklabels=labels, yticklabels=labels,
        square=True,
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
    ax.set_xlabel('Receiver cell type', labelpad=8)
    ax.set_ylabel('Sender cell type', labelpad=8)
    ax.set_title(f'Number of Significant Interactions\n{data.sample_name}', pad=10)
    fig.tight_layout()

    caption = (
        "Heatmap of the number of statistically significant ligand–receptor interactions "
        "(p < 0.05, Cauchy combination test) between each ordered sender–receiver cell type pair. "
        "Pairs with zero interactions are masked."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Fig 02 – mean interaction strength heatmap
# ──────────────────────────────────────────────────────────────────────────────

def plot_strength_heatmap(data: CCCData) -> Tuple[plt.Figure, str]:
    mat = data.strength_matrix
    labels = wrap_labels(mat.index.tolist())

    n = len(labels)
    sz = max(5, min(12, n * 0.55))
    fig, ax = plt.subplots(figsize=(sz + 1.5, sz))

    mask = mat.values == 0
    sns.heatmap(
        mat, ax=ax, mask=mask,
        cmap='Blues', linewidths=0.3, linecolor='#eeeeee',
        cbar_kws={'label': 'Summed communication score', 'shrink': 0.7},
        xticklabels=labels, yticklabels=labels,
        square=True,
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
    ax.set_xlabel('Receiver cell type', labelpad=8)
    ax.set_ylabel('Sender cell type', labelpad=8)
    ax.set_title(f'Communication Strength\n{data.sample_name}', pad=10)
    fig.tight_layout()

    caption = (
        "Heatmap of the summed communication score (CS) across all significant interactions "
        "for each sender–receiver pair, reflecting both the number and magnitude of active "
        "ligand–receptor interactions."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Fig 03 – chord diagram
# ──────────────────────────────────────────────────────────────────────────────

def _bezier_chord(ax, t1s, t1e, t2s, t2e, color, R=1.0, alpha=0.4):
    """Draw a filled bezier chord between two arc segments."""
    R_i = R * 0.95
    ctrl = (0.0, 0.0)

    # Four attachment points (arc start/end for each cell type)
    p1s = (R_i * np.cos(t1s), R_i * np.sin(t1s))
    p1e = (R_i * np.cos(t1e), R_i * np.sin(t1e))
    p2s = (R_i * np.cos(t2s), R_i * np.sin(t2s))
    p2e = (R_i * np.cos(t2e), R_i * np.sin(t2e))

    # Closed path with non-zero area:
    #   p1s --bezier(ctrl)--> p2e
    #   line p2e -> p2s  (width of chord at arc 2)
    #   p2s --bezier(ctrl)--> p1e
    #   CLOSEPOLY back to p1s  (width of chord at arc 1)
    verts = [p1s, ctrl, p2e, p2s, ctrl, p1e, p1s]
    codes = [
        Path.MOVETO,
        Path.CURVE3, Path.CURVE3,
        Path.LINETO,
        Path.CURVE3, Path.CURVE3,
        Path.CLOSEPOLY,
    ]
    path = Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor=color, edgecolor='none',
                               alpha=alpha, zorder=1)
    ax.add_patch(patch)


def plot_chord_diagram(
    data: CCCData, colors_dict: Dict, weight: str = 'count'
) -> Tuple[plt.Figure, str]:

    mat = data.counts_matrix.values.astype(float) if weight == 'count' else data.strength_matrix.values.astype(float)
    labels = data.celltypes
    n = len(labels)

    totals = mat.sum(axis=1) + mat.sum(axis=0)
    for i in range(n):
        totals[i] -= mat[i, i]

    total_sum = totals.sum()

    sz = max(6, min(10, n * 0.5 + 2))
    fig, ax = plt.subplots(figsize=(sz, sz))

    if total_sum == 0 or n == 0:
        ax.text(0.5, 0.5, 'No significant interactions', ha='center', va='center',
                transform=ax.transAxes, fontsize=12)
        ax.axis('off')
        return fig, "Chord diagram (no data)."

    GAP = max(0.015, 0.04 - n * 0.001) * 2 * np.pi
    available = 2 * np.pi - GAP * n
    arc_sizes = (totals / total_sum) * available

    # Start from top, go clockwise
    starts = np.zeros(n)
    ends   = np.zeros(n)
    pos = np.pi / 2
    for i in range(n):
        starts[i] = pos
        ends[i]   = pos + arc_sizes[i]
        pos = ends[i] + GAP

    row_sums = mat.sum(axis=1)
    col_sums = mat.sum(axis=0)

    # Subdivide each arc: first half = outgoing (per target j),
    # second half = incoming (per source j).
    # sub_out[i][j] = (t_start, t_end) of i's outgoing sub-arc toward j
    # sub_in[i][j]  = (t_start, t_end) of i's incoming sub-arc from j
    sub_out = [{} for _ in range(n)]
    sub_in  = [{} for _ in range(n)]

    for i in range(n):
        out_frac = row_sums[i] / totals[i] if totals[i] > 0 else 0.5
        out_size = arc_sizes[i] * out_frac
        in_size  = arc_sizes[i] * (1 - out_frac)

        # outgoing sub-arcs
        cursor = starts[i]
        for j in range(n):
            if mat[i, j] > 0:
                w = mat[i, j] / row_sums[i] * out_size if row_sums[i] > 0 else 0
                sub_out[i][j] = (cursor, cursor + w)
                cursor += w

        # incoming sub-arcs
        cursor = starts[i] + out_size
        for j in range(n):
            if mat[j, i] > 0:
                w = mat[j, i] / col_sums[i] * in_size if col_sums[i] > 0 else 0
                sub_in[i][j] = (cursor, cursor + w)
                cursor += w

    # Draw chords (smallest weight first so large chords are on top)
    connections = [(mat[i, j], i, j) for i in range(n) for j in range(n) if mat[i, j] > 0]
    connections.sort()
    max_w = max(w for w, _, _ in connections) if connections else 1.0

    for w, i, j in connections:
        if j not in sub_out[i] or i not in sub_in[j]:
            continue
        t1s, t1e = sub_out[i][j]
        t2s, t2e = sub_in[j][i]
        alpha = 0.15 + 0.35 * (w / max_w)
        _bezier_chord(ax, t1s, t1e, t2s, t2e,
                      color=colors_dict[labels[i]], R=1.0, alpha=alpha)

    # Draw outer arcs
    R_OUT, R_IN = 1.0, 0.92
    for i in range(n):
        theta = np.linspace(starts[i], ends[i], 200)
        xo, yo = R_OUT * np.cos(theta), R_OUT * np.sin(theta)
        xi, yi = R_IN  * np.cos(theta), R_IN  * np.sin(theta)
        ax.fill(
            np.concatenate([xo, xi[::-1]]),
            np.concatenate([yo, yi[::-1]]),
            color=colors_dict[labels[i]], zorder=3,
        )

        # Label
        t_mid = (starts[i] + ends[i]) / 2
        lr = 1.1
        lx, ly = lr * np.cos(t_mid), lr * np.sin(t_mid)

        rot = np.degrees(t_mid) % 360
        if 90 < rot < 270:
            rot -= 180

        ha = 'left' if np.cos(t_mid) > 0.05 else ('right' if np.cos(t_mid) < -0.05 else 'center')

        ax.text(lx, ly, labels[i], ha=ha, va='center',
                fontsize=max(5, min(8, 80 // n)),
                rotation=rot, rotation_mode='anchor')

    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(
        f'Cell–Cell Communication Chord Diagram\n{data.sample_name} — weighted by {weight}',
        pad=12, fontsize=10,
    )
    fig.tight_layout()

    caption = (
        f"Chord diagram summarising cell–cell communication weighted by interaction {weight}. "
        "Each arc segment represents a cell type; arc width is proportional to total interactions. "
        "Chords connect communicating pairs and are coloured by the sender cell type."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Fig 04 – top sender / receiver bar chart
# ──────────────────────────────────────────────────────────────────────────────

def plot_sender_receiver_bar(
    data: CCCData, colors_dict: Dict, top_n: int = 15
) -> Tuple[plt.Figure, str]:

    mat = data.counts_matrix
    sender   = mat.sum(axis=1).sort_values(ascending=False).head(top_n)
    receiver = mat.sum(axis=0).sort_values(ascending=False).head(top_n)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, max(4, top_n * 0.35)))

    for ax, series, title in [
        (ax1, sender,   'Top Sender Cell Types'),
        (ax2, receiver, 'Top Receiver Cell Types'),
    ]:
        cts = series.index.tolist()
        vals = series.values
        colors = [colors_dict.get(ct, '#999999') for ct in cts]
        bars = ax.barh(range(len(cts)), vals[::-1], color=colors[::-1], edgecolor='white', lw=0.5)
        ax.set_yticks(range(len(cts)))
        ax.set_yticklabels(wrap_labels(cts[::-1]), fontsize=8)
        ax.set_xlabel('Number of significant interactions', fontsize=9)
        ax.set_title(title, fontsize=10, fontweight='bold')
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis='y', length=0)
        for v, bar in zip(vals[::-1], bars):
            ax.text(v + max(vals) * 0.01, bar.get_y() + bar.get_height() / 2,
                    str(int(v)), va='center', ha='left', fontsize=7)

    fig.suptitle(f'Sender and Receiver Ranking — {data.sample_name}', fontsize=11, fontweight='bold')
    fig.tight_layout()

    caption = (
        f"Horizontal bar charts ranking the top {top_n} sender (left) and receiver (right) cell types "
        "by the total number of significant outgoing and incoming ligand–receptor interactions, respectively."
    )
    return fig, caption
