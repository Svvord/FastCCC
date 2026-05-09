"""Module 5 – Network-level analyses: centrality, Sankey, autocrine/paracrine, bipartite."""

from typing import Dict, Tuple

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import networkx as nx

from ..loader import CCCData
from ..utils import wrap_labels


# ──────────────────────────────────────────────────────────────────────────────
# Fig 14 – Directed communication network with centrality metrics
# ──────────────────────────────────────────────────────────────────────────────

def plot_network_centrality(
    data: CCCData, colors_dict: Dict, top_n: int = 20
) -> Tuple[plt.Figure, str]:
    mat = data.counts_matrix
    ct  = data.celltypes

    G = nx.DiGraph()
    G.add_nodes_from(ct)
    for c1 in ct:
        for c2 in ct:
            w = mat.loc[c1, c2]
            if w > 0:
                G.add_edge(c1, c2, weight=float(w))

    if G.number_of_edges() == 0:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, 'No interactions', ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')
        return fig, "Network centrality (no data)."

    out_deg = dict(G.out_degree(weight='weight'))
    in_deg  = dict(G.in_degree(weight='weight'))
    betw    = nx.betweenness_centrality(G, weight='weight', normalized=True)

    pos = nx.spring_layout(G, weight='weight', seed=42, k=2.5 / max(1, np.sqrt(len(ct))))

    fig, ax = plt.subplots(figsize=(10, 9))

    # Edges
    edges = list(G.edges(data=True))
    max_ew = max(d['weight'] for _, _, d in edges) if edges else 1
    for u, v, d in edges:
        x_start, y_start = pos[u]
        x_end,   y_end   = pos[v]
        lw = 0.4 + 3.0 * d['weight'] / max_ew
        ax.annotate(
            '', xy=(x_end, y_end), xytext=(x_start, y_start),
            arrowprops=dict(
                arrowstyle='->', color='#aaaaaa', lw=lw,
                connectionstyle='arc3,rad=0.12',
            ),
            zorder=1,
        )

    # Nodes
    max_b = max(betw.values()) or 1e-9
    for node in G.nodes():
        x, y = pos[node]
        size_r = 0.015 + 0.055 * betw[node] / max_b
        circle = plt.Circle(
            (x, y), size_r,
            color=colors_dict.get(node, '#999999'),
            ec='white', lw=1.0, zorder=3,
        )
        ax.add_patch(circle)
        ax.text(x, y + size_r + 0.01, wrap_labels([node], 18)[0],
                ha='center', va='bottom', fontsize=6.5, zorder=4)

    # Inset table: top cells by betweenness / out / in
    top_b = sorted(betw.items(), key=lambda x: -x[1])[:5]
    top_o = sorted(out_deg.items(), key=lambda x: -x[1])[:5]
    top_i = sorted(in_deg.items(), key=lambda x: -x[1])[:5]
    lines = (
        ['Betweenness (top 5):'] + [f'  {k[:18]}: {v:.3f}' for k, v in top_b] +
        ['\nOut-strength (top 5):'] + [f'  {k[:18]}: {int(v)}' for k, v in top_o] +
        ['\nIn-strength (top 5):]'] + [f'  {k[:18]}: {int(v)}' for k, v in top_i]
    )
    ax.text(1.01, 0.99, '\n'.join(lines), transform=ax.transAxes,
            fontsize=6.0, va='top', ha='left', family='monospace',
            bbox=dict(facecolor='#f8f8f8', edgecolor='#cccccc', boxstyle='round,pad=0.4'))

    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(
        f'Cell–Cell Communication Network — {data.sample_name}\n'
        '(Node size = betweenness centrality; edge width = interaction count)',
        pad=12,
    )
    fig.tight_layout()

    caption = (
        "Directed communication network. Each node is a cell type; directed edges "
        "indicate significant sender→receiver interactions, with edge width proportional "
        "to total interaction count. Node size scales with betweenness centrality, "
        "identifying hub cell types that mediate indirect signalling between others."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Alluvial helper
# ──────────────────────────────────────────────────────────────────────────────

def _smooth_ribbon(ax, x0, x1, y0_bot, y0_top, y1_bot, y1_top, color, alpha=0.45):
    """Fill a smooth S-shaped ribbon between (x0, y0_bot..y0_top) and (x1, y1_bot..y1_top)."""
    t = np.linspace(0, 1, 120)
    # Hermite / smoothstep interpolation: s(t) = 3t²-2t³
    s = 3 * t**2 - 2 * t**3
    x_curve = x0 + t * (x1 - x0)
    y_lower = y0_bot + s * (y1_bot - y0_bot)
    y_upper = y0_top + s * (y1_top - y0_top)
    ax.fill_between(x_curve, y_lower, y_upper, alpha=alpha, color=color, lw=0, zorder=2)


def _alluvial_3col(
    ax, flows: pd.DataFrame,
    col0: str, col1: str, col2: str,
    colors_col0: Dict,
    gap_frac: float = 0.04,
    block_width: float = 0.06,
):
    """Draw a 3-column alluvial (Sankey) diagram."""
    xs = [0.05, 0.5, 0.95]

    def _build_ypos(col):
        totals = flows.groupby(col)['count'].sum().sort_values(ascending=False)
        total_all = totals.sum()
        gap = gap_frac * total_all
        yp = {}
        y = 0
        for val, wt in totals.items():
            yp[val] = (y, y + wt)
            y += wt + gap
        max_y = y - gap
        # Normalise to [0, 1]
        return {v: (bot / max_y, top / max_y) for v, (bot, top) in yp.items()}, max_y

    yp0, _ = _build_ypos(col0)
    yp1, _ = _build_ypos(col1)
    yp2, _ = _build_ypos(col2)
    all_yp = [yp0, yp1, yp2]
    all_cols = [col0, col1, col2]
    all_xs = xs

    # Draw blocks
    for l, (col, xc, yp) in enumerate(zip(all_cols, all_xs, all_yp)):
        for val, (yb, yt) in yp.items():
            if l == 0:
                color = colors_col0.get(val, '#aaaaaa')
            elif l == 2:
                color = '#aaaaaa'
            else:
                color = '#cccccc'
            ax.fill_betweenx([yb, yt], xc - block_width, xc + block_width,
                              color=color, alpha=0.85, lw=0, zorder=3)
            lbl = wrap_labels([val], 18)[0]
            ha = 'right' if l == 0 else ('left' if l == 2 else 'center')
            xoff = -block_width - 0.01 if l == 0 else (block_width + 0.01 if l == 2 else 0)
            ax.text(xc + xoff, (yb + yt) / 2, lbl, ha=ha, va='center',
                    fontsize=6.5, zorder=5, clip_on=False)

    # Draw ribbons: col0 → col1
    cur_r0 = {v: yb for v, (yb, _) in yp0.items()}
    cur_l1 = {v: yb for v, (yb, _) in yp1.items()}
    flows01 = flows.groupby([col0, col1])['count'].sum().reset_index().sort_values('count', ascending=False)
    for _, row in flows01.iterrows():
        v0, v1, cnt = row[col0], row[col1], row['count']
        if v0 not in yp0 or v1 not in yp1:
            continue
        total_v0 = flows[flows[col0] == v0]['count'].sum()
        total_v1 = flows[flows[col1] == v1]['count'].sum()
        h0 = cnt / total_v0 * (yp0[v0][1] - yp0[v0][0]) if total_v0 > 0 else 0
        h1 = cnt / total_v1 * (yp1[v1][1] - yp1[v1][0]) if total_v1 > 0 else 0
        y0b, y1b = cur_r0[v0], cur_l1[v1]
        _smooth_ribbon(ax, xs[0] + block_width, xs[1] - block_width,
                       y0b, y0b + h0, y1b, y1b + h1,
                       color=colors_col0.get(v0, '#aaaaaa'))
        cur_r0[v0] += h0
        cur_l1[v1] += h1

    # Draw ribbons: col1 → col2
    cur_r1 = {v: yb for v, (yb, _) in yp1.items()}
    cur_l2 = {v: yb for v, (yb, _) in yp2.items()}
    flows12 = flows.groupby([col0, col1, col2])['count'].sum().reset_index()
    flows12_agg = flows12.groupby([col1, col2])['count'].sum().reset_index().sort_values('count', ascending=False)

    # Need sender color for col1→col2 ribbons: use dominant sender for each pathway
    dominant_sender = {}
    if col0 in flows.columns:
        for path in yp1.keys():
            sub = flows[flows[col1] == path]
            if not sub.empty:
                dominant_sender[path] = sub.groupby(col0)['count'].sum().idxmax()

    for _, row in flows12_agg.iterrows():
        v1, v2, cnt = row[col1], row[col2], row['count']
        if v1 not in yp1 or v2 not in yp2:
            continue
        total_v1 = flows12_agg[flows12_agg[col1] == v1]['count'].sum()
        total_v2 = flows[flows[col2] == v2]['count'].sum() if col2 in flows.columns else cnt
        h1 = cnt / total_v1 * (yp1[v1][1] - yp1[v1][0]) if total_v1 > 0 else 0
        h2 = cnt / total_v2 * (yp2[v2][1] - yp2[v2][0]) if total_v2 > 0 else 0
        y1b, y2b = cur_r1[v1], cur_l2[v2]
        color = colors_col0.get(dominant_sender.get(v1, ''), '#aaaaaa')
        _smooth_ribbon(ax, xs[1] + block_width, xs[2] - block_width,
                       y1b, y1b + h1, y2b, y2b + h2, color=color)
        cur_r1[v1] += h1
        cur_l2[v2] += h2

    # Column header labels
    for l, (lbl, xc) in enumerate(zip(['Sender', 'Pathway', 'Receiver'], xs)):
        ax.text(xc, 1.04, lbl, ha='center', va='bottom', fontsize=10,
                fontweight='bold', transform=ax.transAxes)

    ax.set_xlim(-0.12, 1.12)
    ax.set_ylim(-0.02, 1.02)
    ax.axis('off')


# ──────────────────────────────────────────────────────────────────────────────
# Fig 15 – Sankey/Alluvial: Sender → Pathway → Receiver
# ──────────────────────────────────────────────────────────────────────────────

def plot_sankey_flow(
    data: CCCData, colors_dict: Dict,
    top_n_ct: int = 10, top_n_path: int = 12,
) -> Tuple[plt.Figure, str]:
    sig = data.significant.copy()
    if sig.empty or 'classification' not in sig.columns:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "Sankey diagram (no data)."

    top_s = sig['sender_celltype'].value_counts().head(top_n_ct).index.tolist()
    top_r = sig['receiver_celltype'].value_counts().head(top_n_ct).index.tolist()
    top_p = sig['classification'].value_counts().head(top_n_path).index.tolist()

    sub = sig[
        sig['sender_celltype'].isin(top_s) &
        sig['classification'].isin(top_p) &
        sig['receiver_celltype'].isin(top_r)
    ]
    if sub.empty:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "Sankey diagram (no data after filtering)."

    flows = (
        sub.groupby(['sender_celltype', 'classification', 'receiver_celltype'])
           .size()
           .reset_index(name='count')
    )

    fig, ax = plt.subplots(figsize=(13, max(7, min(20, len(top_s) * 0.7 + 3))))
    _alluvial_3col(
        ax, flows,
        col0='sender_celltype', col1='classification', col2='receiver_celltype',
        colors_col0=colors_dict,
    )
    ax.set_title(
        f'Sender → Pathway → Receiver Interaction Flow\n{data.sample_name}',
        pad=18, fontsize=11,
    )
    fig.tight_layout()

    caption = (
        f"Alluvial (Sankey) diagram showing how interactions flow from the top {top_n_ct} "
        f"sender cell types through the top {top_n_path} pathway classifications to the top "
        f"{top_n_ct} receiver cell types. Ribbon width is proportional to the number of "
        "significant interactions; ribbons are coloured by the sender cell type."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Fig 16 – Autocrine vs Paracrine interaction ratio
# ──────────────────────────────────────────────────────────────────────────────

def plot_autocrine_paracrine(
    data: CCCData, colors_dict: Dict
) -> Tuple[plt.Figure, str]:
    mat = data.counts_matrix
    ct  = data.celltypes

    rows = []
    for c in ct:
        autocrine  = int(mat.loc[c, c])
        paracrine  = int(mat.loc[c].sum()) - autocrine
        total      = autocrine + paracrine
        rows.append({'celltype': c, 'autocrine': autocrine,
                     'paracrine': paracrine, 'total': total})

    df = pd.DataFrame(rows).sort_values('total', ascending=True)
    df = df[df['total'] > 0]

    if df.empty:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "Autocrine/paracrine plot (no data)."

    fig, ax = plt.subplots(figsize=(9, max(4, len(df) * 0.38 + 1.5)))

    paracrine_color = '#4DBBD5'
    autocrine_color = '#E64B35'

    ax.barh(range(len(df)), df['paracrine'].values, color=paracrine_color,
            label='Paracrine', edgecolor='white', lw=0.4)
    ax.barh(range(len(df)), df['autocrine'].values, left=df['paracrine'].values,
            color=autocrine_color, label='Autocrine', edgecolor='white', lw=0.4)

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(wrap_labels(df['celltype'].tolist()), fontsize=8)
    ax.set_xlabel('Number of significant outgoing interactions', fontsize=9)
    ax.set_title(f'Autocrine vs. Paracrine Signalling\n{data.sample_name}', pad=10)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)

    # Percentage labels for autocrine fraction
    for i, row in enumerate(df.itertuples()):
        if row.total > 0:
            pct = row.autocrine / row.total * 100
            if pct > 0:
                ax.text(row.paracrine + row.autocrine + df['total'].max() * 0.01,
                        i, f'{pct:.0f}%', va='center', fontsize=6.5, color='#666')

    ax.legend(loc='lower right', fontsize=8)
    fig.tight_layout()

    caption = (
        "Stacked horizontal bar chart quantifying the number of outgoing paracrine "
        "(interactions with other cell types, blue) and autocrine (self-interactions, red) "
        "significant interactions per cell type. The percentage of autocrine interactions "
        "is annotated at the right of each bar."
    )
    return fig, caption


# ──────────────────────────────────────────────────────────────────────────────
# Fig 17 – L-R bipartite network
# ──────────────────────────────────────────────────────────────────────────────

def plot_bipartite_lr(
    data: CCCData, colors_dict: Dict,
    top_n_ligand: int = 18, top_n_receptor: int = 18,
) -> Tuple[plt.Figure, str]:
    sig = data.significant.copy()
    if sig.empty:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "Bipartite L-R network (no data)."

    top_lig = sig['ligand'].value_counts().head(top_n_ligand).index.tolist()
    top_rec = sig['receptor'].value_counts().head(top_n_receptor).index.tolist()
    sub = sig[sig['ligand'].isin(top_lig) & sig['receptor'].isin(top_rec)]
    if sub.empty:
        fig, ax = plt.subplots()
        ax.axis('off')
        return fig, "Bipartite L-R network (no data after filtering)."

    # Node sizes: number of significant interactions
    lig_counts = sub['ligand'].value_counts()
    rec_counts = sub['receptor'].value_counts()

    # Primary sender/receiver color per ligand/receptor
    lig_colors = {
        lig: colors_dict.get(
            sub[sub['ligand'] == lig]['sender_celltype'].value_counts().idxmax(),
            '#4DBBD5'
        ) for lig in top_lig if lig in lig_counts
    }
    rec_colors = {
        rec: colors_dict.get(
            sub[sub['receptor'] == rec]['receiver_celltype'].value_counts().idxmax(),
            '#E64B35'
        ) for rec in top_rec if rec in rec_counts
    }

    # Layout: ligands at x=0, receptors at x=1
    n_l = len(top_lig)
    n_r = len(top_rec)
    lig_y = {lig: 1 - i / max(n_l - 1, 1) for i, lig in enumerate(top_lig)}
    rec_y = {rec: 1 - i / max(n_r - 1, 1) for i, rec in enumerate(top_rec)}

    fig_h = max(6, max(n_l, n_r) * 0.38 + 2)
    fig, ax = plt.subplots(figsize=(10, fig_h))

    # Edges
    lr_edges = sub.groupby(['ligand', 'receptor']).size().reset_index(name='count')
    max_cnt = lr_edges['count'].max() or 1
    for _, erow in lr_edges.iterrows():
        lig, rec, cnt = erow['ligand'], erow['receptor'], erow['count']
        if lig not in lig_y or rec not in rec_y:
            continue
        color = lig_colors.get(lig, '#aaaaaa')
        lw = 0.3 + 2.5 * cnt / max_cnt
        ax.plot([0.2, 0.8], [lig_y[lig], rec_y[rec]],
                color=color, lw=lw, alpha=0.35, zorder=1)

    # Nodes
    max_ls = lig_counts.max() or 1
    max_rs = rec_counts.max() or 1
    for lig in top_lig:
        if lig not in lig_counts:
            continue
        s = 80 + 300 * lig_counts[lig] / max_ls
        ax.scatter([0.2], [lig_y[lig]], s=s, color=lig_colors.get(lig, '#4DBBD5'),
                   edgecolors='white', lw=0.8, zorder=3)
        ax.text(0.18, lig_y[lig], lig, ha='right', va='center', fontsize=7)

    for rec in top_rec:
        if rec not in rec_counts:
            continue
        s = 80 + 300 * rec_counts[rec] / max_rs
        ax.scatter([0.8], [rec_y[rec]], s=s, color=rec_colors.get(rec, '#E64B35'),
                   edgecolors='white', lw=0.8, zorder=3)
        ax.text(0.82, rec_y[rec], rec, ha='left', va='center', fontsize=7)

    ax.text(0.2, 1.04, 'Ligand', ha='center', fontsize=11, fontweight='bold', transform=ax.transAxes)
    ax.text(0.8, 1.04, 'Receptor', ha='center', fontsize=11, fontweight='bold', transform=ax.transAxes)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.05, 1.1)
    ax.axis('off')
    ax.set_title(
        f'Ligand–Receptor Bipartite Network\n{data.sample_name}',
        pad=14,
    )
    fig.tight_layout()

    caption = (
        f"Bipartite network connecting the top {top_n_ligand} ligands (left) and top "
        f"{top_n_receptor} receptors (right) by significant L-R pairings. Node size "
        "reflects interaction frequency. Ligand nodes are coloured by their dominant "
        "sender cell type; edges are weighted by pairwise interaction count."
    )
    return fig, caption
