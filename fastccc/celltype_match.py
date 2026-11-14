"""Match query cell-type labels to the cell types of a reference panel.

The reference workflow needs a ``{reference cell type: query cell type}``
mapping whenever the two label sets differ (shipped panels use Cell Ontology
names such as ``"CD4-positive, alpha-beta T cell"`` while users write
``"CD4 T cells"``).  This module builds that mapping deterministically from
name normalisation, a curated synonym table and weighted token overlap, and
validates mappings / reports coverage before a comparison is run.
"""
from __future__ import annotations

import difflib
import json
import os
import re
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Iterable, Mapping, Optional, Sequence

import pandas as pd

# ── vocabulary ───────────────────────────────────────────────────────────────

# tokens that carry no identity information
_STOPWORDS = {
    'cell', 'cells', 'of', 'the', 'a', 'an', 'and', 'human', 'derived', 'thymus', 'terminally',
    'differentiated', 'stage', 'lineage', 'like', 'subset', 'population', 'cluster', 'positive',
}

# lineage tokens decide what a cell *is*; they weigh more than modifiers
_LINEAGE = {
    't', 'b', 'nk', 'monocyte', 'macrophage', 'dendritic', 'neutrophil', 'eosinophil', 'basophil', 'mast',
    'plasma', 'plasmablast', 'platelet', 'erythrocyte', 'erythroblast', 'erythroid', 'granulocyte',
    'thymocyte', 'lymphocyte', 'leukocyte', 'ilc', 'stem', 'progenitor', 'precursor',
    'fibroblast', 'myofibroblast', 'endothelial', 'epithelial', 'pericyte', 'smooth', 'muscle', 'mesenchymal',
    'stromal', 'adipocyte', 'chondrocyte', 'osteoblast', 'osteoclast', 'keratinocyte', 'melanocyte',
    'hepatocyte', 'cholangiocyte', 'kupffer', 'langerhans', 'microglia', 'astrocyte', 'oligodendrocyte',
    'opc', 'neuron', 'interneuron', 'schwann', 'glial', 'ependymal', 'cardiomyocyte', 'podocyte', 'mesangial',
    'goblet', 'enterocyte', 'enteroendocrine', 'paneth', 'tuft', 'ciliated', 'club', 'basal', 'luminal',
    'acinar', 'ductal', 'alveolar', 'ionocyte', 'trophoblast', 'decidual', 'megakaryocyte', 'photoreceptor',
    'bipolar', 'amacrine', 'horizontal', 'ganglion', 'muller', 'rod', 'cone', 'germ', 'sertoli', 'leydig',
    'oocyte', 'myocyte', 'satellite', 'tubule', 'intercalated', 'principal', 'secretory', 'transit',
}

# abbreviations / synonyms -> canonical token list (applied to whole tokens)
_ALIASES: dict[str, tuple[str, ...]] = {
    'nk': ('nk',), 'nks': ('nk',), 'natural': ('natural',), 'killer': ('killer',),
    'nkt': ('nk', 't'), 'dc': ('dendritic',), 'dcs': ('dendritic',), 'pdc': ('plasmacytoid', 'dendritic'),
    'pdcs': ('plasmacytoid', 'dendritic'), 'cdc': ('conventional', 'dendritic'), 'cdcs': ('conventional', 'dendritic'),
    'cdc1': ('conventional', 'dendritic', '1'), 'cdc2': ('conventional', 'dendritic', '2'),
    'mdc': ('myeloid', 'dendritic'), 'modc': ('monocyte', 'dendritic'),
    'mono': ('monocyte',), 'monos': ('monocyte',), 'mac': ('macrophage',), 'macs': ('macrophage',),
    'mph': ('macrophage',), 'mphi': ('macrophage',), 'macro': ('macrophage',),
    'treg': ('regulatory', 't'), 'tregs': ('regulatory', 't'), 'tfh': ('follicular', 'helper', 't'),
    'th': ('helper', 't'), 'th1': ('helper', 't', '1'), 'th2': ('helper', 't', '2'), 'th17': ('helper', 't', '17'),
    'mait': ('mucosal', 'invariant', 't'), 'gd': ('gamma', 'delta'), 'gdt': ('gamma', 'delta', 't'),
    'γδ': ('gamma', 'delta'), 'ctl': ('cytotoxic', 't'), 'cd4t': ('cd4', 't'), 'cd8t': ('cd8', 't'),
    'tc': ('t',), 'tcell': ('t',), 'tcells': ('t',), 'bcell': ('b',), 'bcells': ('b',), 'nkcell': ('nk',),
    'nkcells': ('nk',), 'tn': ('naive', 't'), 'tcm': ('central', 'memory', 't'), 'tem': ('effector', 'memory', 't'),
    'temra': ('effector', 'memory', 't', 'cd45ra'), 'trm': ('resident', 'memory', 't'),
    'hsc': ('hematopoietic', 'stem'), 'hscs': ('hematopoietic', 'stem'), 'hspc': ('hematopoietic', 'stem', 'progenitor'),
    'cmp': ('common', 'myeloid', 'progenitor'), 'gmp': ('granulocyte', 'monocyte', 'progenitor'),
    'clp': ('common', 'lymphoid', 'progenitor'), 'mep': ('megakaryocyte', 'erythroid', 'progenitor'),
    'ilc': ('ilc',), 'ilcs': ('ilc',), 'ilc1': ('ilc', '1'), 'ilc2': ('ilc', '2'), 'ilc3': ('ilc', '3'),
    'innate': ('innate',), 'mk': ('platelet',), 'megakaryocyte': ('platelet',), 'megakaryocytes': ('platelet',),
    'platelets': ('platelet',), 'rbc': ('erythrocyte',), 'rbcs': ('erythrocyte',), 'ery': ('erythroid',),
    'neut': ('neutrophil',), 'neuts': ('neutrophil',), 'eos': ('eosinophil',), 'baso': ('basophil',),
    'pb': ('plasmablast',), 'pc': ('plasma',), 'fib': ('fibroblast',), 'fibro': ('fibroblast',), 'fb': ('fibroblast',),
    'fbs': ('fibroblast',), 'myofib': ('myofibroblast',), 'caf': ('fibroblast',), 'cafs': ('fibroblast',),
    'endo': ('endothelial',), 'ec': ('endothelial',), 'ecs': ('endothelial',), 'endothelium': ('endothelial',),
    'lec': ('lymphatic', 'endothelial',), 'epi': ('epithelial',), 'epithelium': ('epithelial',),
    'smc': ('smooth', 'muscle'), 'smcs': ('smooth', 'muscle'), 'vsmc': ('smooth', 'muscle'), 'peri': ('pericyte',),
    'cm': ('cardiomyocyte',), 'cms': ('cardiomyocyte',), 'hep': ('hepatocyte',), 'chol': ('cholangiocyte',),
    'msc': ('mesenchymal',), 'stroma': ('stromal',), 'astro': ('astrocyte',), 'oligo': ('oligodendrocyte',),
    'ol': ('oligodendrocyte',), 'opc': ('opc',), 'opcs': ('opc',), 'micro': ('microglia',), 'mg': ('microglia',),
    'kera': ('keratinocyte',), 'mel': ('melanocyte',), 'adipo': ('adipocyte',), 'chondro': ('chondrocyte',),
    'sct': ('trophoblast',), 'vct': ('trophoblast',), 'evt': ('trophoblast',), 'mono_derived': ('monocyte',),
    'non': ('non',), 'nonclassical': ('nonclassical',), 'classical': ('classical',),
    'prolif': ('proliferating',), 'cycling': ('proliferating',), 'proliferative': ('proliferating',),
    'naïve': ('naive',), 'mem': ('memory',), 'eff': ('effector',), 'act': ('activated',),
}

# attribute groups: values on both sides must agree; (query-only, reference-only) handled in _compatible
_ATTRIBUTE_GROUPS = (
    {'cd4', 'cd8'},
    {'classical', 'nonclassical'},
    {'conventional', 'plasmacytoid'},
    {'alphabeta', 'gammadelta'},
)
# tokens that imply an attribute when the group is not already specified
_IMPLIED = {
    'regulatory': 'cd4', 'helper': 'cd4', 'follicular': 'cd4', 'cytotoxic': 'cd8',
    'cd14': 'classical', 'cd16': 'nonclassical', 'fcgr3a': 'nonclassical', 'gammadelta': 'gammadelta',
}
_MARKER_RE = re.compile(r'^(cd\d+[a-z]*|[a-z]{2,8}\d[a-z0-9]*|ig[agmde])$')


def normalise_label(label: str) -> set[str]:
    """Turn a cell-type label into a canonical token set."""
    s = str(label).strip().lower()
    s = s.replace('γδ', ' gammadelta ').replace('αβ', ' alphabeta ').replace('φ', 'phi')
    s = re.sub(r'\bnon[\s\-]?classical\b', 'nonclassical', s)
    s = re.sub(r'\bgamma[\s\-]?delta\b', 'gammadelta', s)
    s = re.sub(r'\balpha[\s\-]?beta\b', 'alphabeta', s)
    s = re.sub(r'([a-z0-9])\+', r'\1 positive ', s)           # CD14+  -> cd14 positive
    s = re.sub(r'([a-z0-9])[\-−]\s', r'\1 negative ', s)      # CD16-  -> cd16 negative
    s = re.sub(r'([a-z0-9])[\-−]$', r'\1 negative ', s)
    s = re.sub(r'[^a-z0-9γδ]+', ' ', s)
    raw = s.split()
    toks: list[str] = []
    i = 0
    while i < len(raw):
        t = raw[i]
        nxt = raw[i + 1] if i + 1 < len(raw) else ''
        if t == 'natural' and nxt == 'killer':
            toks.append('nk'); i += 2; continue
        if nxt in ('positive', 'negative', 'bright', 'dim', 'high', 'low', 'hi', 'lo') and (_MARKER_RE.match(t) or t.startswith('cd')):
            toks.append(t if nxt == 'positive' else f'{t}{nxt}'); i += 2; continue
        toks.extend(_ALIASES.get(t, (t,)))
        i += 1
    out: set[str] = set()
    for t in toks:
        if t in _STOPWORDS or t == 'alphabeta':
            continue
        t = _singular(t)
        out.add(t)
    for t in list(out):                                       # implied attributes
        imp = _IMPLIED.get(t)
        if imp:
            group = next((g for g in _ATTRIBUTE_GROUPS if imp in g), None)
            if group is None or not (out & group):
                out.add(imp)
    if 'nonclassical' in out:
        out.discard('non'); out.discard('classical')
    return out


def _singular(t: str) -> str:
    if t in _ALIASES or len(t) <= 3 or t.isdigit():
        return t
    for suf, rep in (('cytes', 'cyte'), ('blasts', 'blast'), ('phages', 'phage'), ('ils', 'il'), ('ocytes', 'ocyte')):
        if t.endswith(suf):
            return t[: -len(suf)] + rep
    if t.endswith('s') and not t.endswith(('ss', 'us', 'is', 'os')):
        return t[:-1]
    return t


def _weight(tok: str) -> float:
    return 3.0 if tok in _LINEAGE else 1.0


def _compatible(ref: set[str], qry: set[str]) -> bool:
    for group in _ATTRIBUTE_GROUPS:
        r, q = ref & group, qry & group
        if r and q and r != q:
            return False
    # a gamma-delta T reference must not feed a CD4/CD8 (alpha-beta) query
    if ('gammadelta' in ref) != ('gammadelta' in qry) and ({'cd4', 'cd8'} & (ref | qry)):
        return False
    return True


def label_similarity(reference_label: str, query_label: str) -> float:
    """Weighted Jaccard similarity of the canonical token sets (0 = no shared lineage)."""
    r, q = normalise_label(reference_label), normalise_label(query_label)
    return _token_similarity(r, q)


def _token_similarity(r: set[str], q: set[str]) -> float:
    if not r or not q or not _compatible(r, q):
        return 0.0
    shared = r & q
    if not (shared & _LINEAGE) and not ((r & _LINEAGE) == set() and (q & _LINEAGE) == set() and shared):
        return 0.0
    num = sum(_weight(t) for t in shared)
    den = sum(_weight(t) for t in r | q)
    return round(num / den, 3) if den else 0.0


# ── suggestion ───────────────────────────────────────────────────────────────

@dataclass
class MappingSuggestion:
    """Result of :func:`suggest_celltype_mapping`."""
    mapping: dict[str, str]                 # reference cell type -> query cell type (ready for compare_to_reference)
    entries: list[dict]                     # one per mapped reference type: reference, query, score, confidence
    ambiguous: list[dict]                   # reference types with tied candidates (not mapped)
    unmatched_reference: list[str]
    unmatched_query: list[str]
    query_counts: dict[str, int]
    reference_counts: dict[str, int]
    coverage: dict
    reference_name: str = ''
    min_score: float = 0.5

    def to_dict(self) -> dict:
        return asdict(self)

    def to_json(self, path: str | os.PathLike, full: bool = False) -> Path:
        """Write the mapping (or, with ``full=True``, the whole suggestion) as JSON."""
        path = Path(path).expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, 'w', encoding='utf-8') as fh:
            json.dump(self.to_dict() if full else self.mapping, fh, indent=2, ensure_ascii=False)
        return path

    def table(self) -> pd.DataFrame:
        """One row per query cell type: mapped reference types, cells, status."""
        rows = []
        by_query: dict[str, list[dict]] = {}
        for e in self.entries:
            by_query.setdefault(e['query'], []).append(e)
        for q, n in sorted(self.query_counts.items(), key=lambda kv: -kv[1]):
            ents = by_query.get(q, [])
            exact = q in self.reference_counts and q not in self.mapping
            refs = [e['reference'] for e in ents] + ([q] if exact else [])
            conf = min((e['confidence'] for e in ents), key=_conf_rank) if ents else ('high' if exact else '')
            rows.append({
                'query_celltype': q, 'n_cells': n, 'n_reference_types': len(refs),
                'reference_celltypes': '; '.join(refs), 'confidence': conf,
                'status': 'mapped' if refs else 'NOT IN REFERENCE',
            })
        return pd.DataFrame(rows)

    def format(self, max_refs: int = 12) -> str:
        """Human-readable listing: one block per query cell type with its mapped reference types."""
        by_query: dict[str, list[dict]] = {}
        for e in self.entries:
            by_query.setdefault(e['query'], []).append(e)
        lines = []
        for q, n in sorted(self.query_counts.items(), key=lambda kv: -kv[1]):
            ents = sorted(by_query.get(q, []), key=lambda e: -e['score'])
            exact = q in self.reference_counts and q not in self.mapping
            if not ents and not exact:
                lines.append(f"{q} ({n:,} cells)  <- NOT IN REFERENCE")
                continue
            total = len(ents) + (1 if exact else 0)
            lines.append(f"{q} ({n:,} cells)  <- {total} reference cell type{'s' if total != 1 else ''}")
            if exact:
                lines.append(f"    {q}  [same name]")
            for e in ents[:max_refs]:
                lines.append(f"    {e['reference']}  [{e['confidence']}, {e['score']:.2f}]")
            if len(ents) > max_refs:
                lines.append(f"    ... {len(ents) - max_refs} more")
        if self.ambiguous:
            lines.append("")
            lines.append("ambiguous reference cell types (not mapped; add them to the JSON yourself if one query type is right):")
            for a in self.ambiguous:
                cands = ' / '.join(c['query'] for c in a['candidates'])
                lines.append(f"    {a['reference']}  ~  {cands}")
        return '\n'.join(lines)

    def summary(self) -> str:
        cov = self.coverage
        lines = [f"{len(self.mapping)} reference cell types mapped onto {cov['n_covered']}/{cov['n_query_types']} "
                 f"query cell types ({cov['cells_covered_fraction']:.1%} of query cells)"]
        if self.unmatched_query:
            lines.append(f"query cell types without a reference counterpart: {', '.join(self.unmatched_query)}")
        if self.ambiguous:
            amb = '; '.join(f"{a['reference']} ~ {' / '.join(c['query'] for c in a['candidates'])}" for a in self.ambiguous)
            lines.append(f"ambiguous (not mapped, pick manually): {amb}")
        return '\n'.join(lines)

    def __repr__(self) -> str:
        return f"MappingSuggestion({self.summary().splitlines()[0]})"


def _conf_rank(c: str) -> int:
    return {'high': 0, 'medium': 1, 'low': 2}.get(c, 3)


def _confidence(score: float) -> str:
    return 'high' if score >= 0.9 else ('medium' if score >= 0.65 else 'low')


def _as_counts(labels) -> dict[str, int]:
    if isinstance(labels, Mapping):
        return {str(k): int(v) for k, v in labels.items()}
    s = pd.Series(list(labels)).astype(str)
    return {str(k): int(v) for k, v in s.value_counts().items()}


def suggest_mapping_from_labels(
    reference_counts: Mapping[str, int], query_counts: Mapping[str, int],
    min_score: float = 0.5, reference_name: str = '',
) -> MappingSuggestion:
    """Core matcher working on ``{label: n_cells}`` dicts for both sides."""
    ref_counts = _as_counts(reference_counts)
    qry_counts = _as_counts(query_counts)
    ref_tok = {r: normalise_label(r) for r in ref_counts}
    qry_tok = {q: normalise_label(q) for q in qry_counts}

    mapping: dict[str, str] = {}
    entries, ambiguous, unmatched_ref = [], [], []
    for r, rt in ref_tok.items():
        if r in qry_counts:                       # identical label: the core uses it directly
            continue
        scored = sorted(((_token_similarity(rt, qt), q) for q, qt in qry_tok.items()), key=lambda x: (-x[0], x[1]))
        best = scored[0][0] if scored else 0.0
        if best < min_score:
            unmatched_ref.append(r)
            continue
        top = [sq for sq in scored if sq[0] == best]
        if len(top) > 1:
            generic = [q for _, q in top if qry_tok[q] <= rt]     # a coarser query label absorbs the tie
            if len(generic) == 1:
                top = [(best, generic[0])]
        if len(top) > 1:
            ambiguous.append({'reference': r, 'score': best,
                              'candidates': [{'query': q, 'score': s} for s, q in top]})
            continue
        q = top[0][1]
        mapping[r] = q
        shared = sorted(rt & qry_tok[q])
        entries.append({'reference': r, 'query': q, 'score': best, 'confidence': _confidence(best),
                        'reason': f"shared: {', '.join(shared)}"})

    coverage = compute_coverage(ref_counts, qry_counts, mapping)
    return MappingSuggestion(
        mapping=mapping, entries=entries, ambiguous=ambiguous, unmatched_reference=unmatched_ref,
        unmatched_query=coverage['uncovered_query_types'], query_counts=qry_counts, reference_counts=ref_counts,
        coverage=coverage, reference_name=reference_name, min_score=min_score,
    )


def compute_coverage(reference_counts: Mapping[str, int], query_counts: Mapping[str, int],
                     mapping: Optional[Mapping[str, str]] = None) -> dict:
    """Which query cell types the reference will cover, following the core's matching rule.

    A query type is covered when a mapping entry points to it, or when a
    reference type with the identical name exists and is not itself used as a
    mapping key.
    """
    mapping = dict(mapping or {})
    qry = _as_counts(query_counts)
    covered = set(mapping.values()) & set(qry)
    covered |= {q for q in qry if q in reference_counts and q not in mapping}
    uncovered = [q for q in qry if q not in covered]
    total = sum(qry.values()) or 1
    return {
        'n_query_types': len(qry), 'n_covered': len(covered),
        'cells_covered_fraction': round(sum(qry[q] for q in covered) / total, 4),
        'uncovered_query_types': sorted(uncovered, key=lambda q: -qry[q]),
        'uncovered_cells_fraction': round(sum(qry[q] for q in uncovered) / total, 4),
        'n_reference_types_used': len(mapping) + len({q for q in covered if q in reference_counts and q not in mapping}),
    }


def validate_celltype_mapping(mapping: Mapping[str, str], reference_counts: Mapping[str, int],
                              query_counts: Mapping[str, int]) -> dict[str, str]:
    """Check a user mapping against both label sets; raise ValueError with hints on mistakes."""
    ref, qry = set(_as_counts(reference_counts)), set(_as_counts(query_counts))
    problems = []
    for k, v in mapping.items():
        if str(k) not in ref:
            hint = difflib.get_close_matches(str(k), ref, n=1, cutoff=0.6)
            problems.append(f"reference cell type '{k}' is not in the panel" + (f" (did you mean '{hint[0]}'?)" if hint else ''))
        if str(v) not in qry:
            hint = difflib.get_close_matches(str(v), qry, n=1, cutoff=0.6)
            problems.append(f"query cell type '{v}' is not in the data" + (f" (did you mean '{hint[0]}'?)" if hint else ''))
    if problems:
        raise ValueError("celltype_mapping problems (keys are reference cell types, values are query cell types):\n  - "
                         + "\n  - ".join(problems))
    return {str(k): str(v) for k, v in mapping.items()}


def load_mapping(mapping: Mapping[str, str] | str | os.PathLike) -> dict[str, str]:
    if isinstance(mapping, Mapping):
        return {str(k): str(v) for k, v in mapping.items()}
    path = Path(mapping).expanduser()
    if not path.is_file():
        raise FileNotFoundError(f"celltype mapping file not found: {path}")
    with open(path, encoding='utf-8') as fh:
        data = json.load(fh)
    if isinstance(data, dict) and 'mapping' in data and isinstance(data['mapping'], dict):
        data = data['mapping']                       # a full MappingSuggestion dump
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a JSON object {{reference cell type: query cell type}}")
    return {str(k): str(v) for k, v in data.items()}
