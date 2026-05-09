"""Shared plot settings, color palettes, and I/O helpers."""

import io
import base64
from pathlib import Path
from typing import Dict, List

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def set_publication_style():
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 9,
        'axes.titlesize': 10,
        'axes.titleweight': 'bold',
        'axes.labelsize': 9,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
        'legend.frameon': False,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
    })


# Curated 40-color palette for single-cell publication figures.
# Sources: ggsci NPG/AAAS/JCO/Lancet palettes + ColorBrewer Dark2/Set1 +
#          hand-picked extensions. Ordered by visual distinctiveness.
_SC_PALETTE = [
    '#E64B35',  # NPG red
    '#4DBBD5',  # NPG teal
    '#00A087',  # NPG emerald
    '#3C5488',  # NPG navy
    '#F39B7F',  # NPG salmon
    '#8491B4',  # NPG periwinkle
    '#91D1C2',  # NPG mint
    '#7E6148',  # NPG brown
    '#3B4992',  # AAAS blue
    '#EE0000',  # AAAS red
    '#008B45',  # AAAS green
    '#631879',  # AAAS purple
    '#008280',  # AAAS teal
    '#BB0021',  # AAAS crimson
    '#5F559B',  # AAAS violet
    '#A20056',  # AAAS magenta
    '#0073C2',  # JCO blue
    '#EFC000',  # JCO gold
    '#CD534C',  # JCO brick
    '#003C67',  # JCO dark navy
    '#A73030',  # JCO dark red
    '#1B9E77',  # CB Dark2 teal
    '#D95F02',  # CB Dark2 orange
    '#7570B3',  # CB Dark2 purple
    '#E7298A',  # CB Dark2 pink
    '#66A61E',  # CB Dark2 olive
    '#E6AB02',  # CB Dark2 amber
    '#A6761D',  # CB Dark2 sienna
    '#2166AC',  # CB RdBu blue
    '#B2182B',  # CB RdBu red
    '#4DAF4A',  # CB Set1 green
    '#984EA3',  # CB Set1 purple
    '#FF7F00',  # CB Set1 orange
    '#A65628',  # CB Set1 brown
    '#2D6A4F',  # forest green
    '#6D4C41',  # coffee brown
    '#455A64',  # blue-grey
    '#AD1457',  # deep pink
    '#1565C0',  # deep blue
    '#558B2F',  # dark lime green
]


def get_celltype_colors(celltypes: List[str]) -> Dict[str, tuple]:
    n = len(celltypes)
    palette = _SC_PALETTE[:n] if n <= len(_SC_PALETTE) else (
        _SC_PALETTE + sns.color_palette('husl', n - len(_SC_PALETTE))
    )
    return {ct: palette[i] for i, ct in enumerate(sorted(celltypes))}


def save_figure(fig: plt.Figure, output_dir: Path, name: str, dpi: int = 300):
    fig_dir = output_dir / 'figures'
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / f'{name}.png', dpi=dpi, bbox_inches='tight')
    try:
        fig.savefig(fig_dir / f'{name}.svg', bbox_inches='tight')
    except Exception:
        pass


def fig_to_base64(fig: plt.Figure, dpi: int = 150) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight')
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('utf-8')


def wrap_labels(labels: List[str], max_len: int = 25) -> List[str]:
    result = []
    for lb in labels:
        if len(lb) > max_len:
            lb = lb[:max_len - 2] + '..'
        result.append(lb)
    return result


def diverging_cmap():
    return sns.diverging_palette(220, 20, as_cmap=True)


def sequential_cmap():
    return 'YlOrRd'
