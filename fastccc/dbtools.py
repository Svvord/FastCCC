"""Ligand–receptor databases in FastCCC's table format: build, validate, describe, import.

A database is a directory holding five CSV tables (the layout of CellPhoneDB, which the
scoring code reads) and a ``database.json`` with its metadata::

    interaction_table.csv          id_interaction, id_cp_interaction, multidata_1_id, multidata_2_id,
                                   [source], [directionality], [classification]
    gene_table.csv                 id_gene, ensembl, gene_name, hgnc_symbol, protein_name, protein_id
    protein_table.csv              id_protein, protein_name, protein_multidata_id
    complex_table.csv              id_complex, complex_multidata_id, [name]
    complex_composition_table.csv  id_complex_composition, complex_multidata_id, protein_multidata_id, total_protein
    database.json                  name, species, source, version, citation, license, url, built,
                                   n_interactions, n_genes, n_complexes, has_classification

``multidata_1_id`` is the ligand side, ``multidata_2_id`` the receptor side; each is either a
protein's ``protein_multidata_id`` or a complex's ``complex_multidata_id``. The gene table maps
one gene (``hgnc_symbol`` holds the gene symbol of the database's species; ``ensembl`` is filled
when a mapping is available) to each protein. ``classification`` is optional and feeds the
explorer's pathway-family panel.

:func:`build_database` turns a plain table of ligand–receptor pairs (subunits joined by ``_``
or ``,``) into this layout; :func:`validate_database` checks a directory; :func:`database_info`
reads its metadata; :func:`import_database` installs a directory or a pair table under
``$FASTCCC_HOME/databases`` so that ``--db <name>`` finds it.
"""
from __future__ import annotations

import datetime as _dt
import json
import os
import re
import shutil
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

from .resources import PACKAGE_DIR, PACKAGED_DB_DIR, user_database_dir, _is_database_dir

TABLES = {
    'interaction_table.csv': ['id_interaction', 'id_cp_interaction', 'multidata_1_id', 'multidata_2_id'],
    'gene_table.csv': ['id_gene', 'ensembl', 'gene_name', 'hgnc_symbol', 'protein_name', 'protein_id'],
    'protein_table.csv': ['id_protein', 'protein_name', 'protein_multidata_id'],
    'complex_table.csv': ['id_complex', 'complex_multidata_id'],
    'complex_composition_table.csv': ['id_complex_composition', 'complex_multidata_id', 'protein_multidata_id', 'total_protein'],
}
METADATA_FILE = 'database.json'
METADATA_KEYS = ('name', 'species', 'version', 'source', 'citation', 'license', 'url', 'notes', 'gene_id', 'built',
                 'n_interactions', 'n_genes', 'n_complexes', 'has_classification')
SPECIES = ('human', 'mouse', 'rat', 'zebrafish', 'other')
SUBUNIT_SEP = re.compile(r'[_,+&]\s*')
ENSEMBL_MAPS = {'human': PACKAGE_DIR / 'data' / 'db' / 'ensembl_hgnc_Homo_sapiens.GRCh38.112.csv'}
NAME_RE = re.compile(r'^[A-Za-z0-9][A-Za-z0-9._-]*$')


# ── helpers ──────────────────────────────────────────────────────────────────

def split_subunits(name) -> tuple[str, ...]:
    """``'ITGA1_ITGB1'`` / ``'ITGA1,ITGB1'`` → ``('ITGA1', 'ITGB1')``; a ``COMPLEX:`` prefix is dropped."""
    s = str(name).strip()
    if s.upper().startswith('COMPLEX:'):
        s = s[8:]
    return tuple(p.strip() for p in SUBUNIT_SEP.split(s) if p.strip())


def load_ensembl_map(species: str = 'human') -> dict[str, str]:
    """``{symbol: Ensembl gene id}`` for the species, from the mapping shipped with the package (empty if none)."""
    path = ENSEMBL_MAPS.get(species)
    if path is None or not path.is_file():
        return {}
    m = pd.read_csv(path)
    cols = {c.lower(): c for c in m.columns}
    ens, sym = cols.get('ensembl'), cols.get('hgnc') or cols.get('symbol')
    m = m[[ens, sym]].dropna().drop_duplicates(subset=[sym], keep='first')
    return dict(zip(m[sym].astype(str), m[ens].astype(str)))


def _has_classification(inter: pd.DataFrame) -> bool:
    return bool('classification' in inter.columns and inter['classification'].notna().any())


def _check_name(name: str) -> None:
    if not NAME_RE.match(name or ''):
        raise ValueError(f"database name {name!r} must start with a letter or digit and contain only letters, digits, '.', '_' or '-'")


# ── build ────────────────────────────────────────────────────────────────────

def build_database(pairs: pd.DataFrame, name: str, out_dir: str | os.PathLike, *, species: str = 'human',
                   source: Optional[str] = None, version: Optional[str] = None, citation: Optional[str] = None,
                   license: Optional[str] = None, url: Optional[str] = None, notes: Optional[str] = None,
                   ensembl_map: Optional[dict[str, str]] = None, overwrite: bool = False) -> Path:
    """Write the database directory ``out_dir/name`` from a table of pairs and return it.

    ``pairs`` needs columns ``ligand`` and ``receptor`` (gene symbols; complex subunits joined
    by ``_`` or ``,``). Optional columns ``id`` (unique interaction id), ``classification``
    (pathway family), ``directionality`` and ``source`` (evidence) are carried into the
    interaction table. Rows with an empty side are dropped and duplicate pairs keep their first row.
    """
    if species not in SPECIES:
        raise ValueError(f"species must be one of {SPECIES}, got {species!r}")
    _check_name(name)
    out = Path(out_dir).expanduser() / name
    if out.exists() and any(out.iterdir()) and not overwrite:
        raise FileExistsError(f"{out} exists; pass overwrite=True to replace it")
    df = pairs.copy()
    for col in ('ligand', 'receptor'):
        if col not in df.columns:
            raise ValueError(f"pairs needs a '{col}' column (got {list(df.columns)})")
    df['lig'] = df['ligand'].map(split_subunits)
    df['rec'] = df['receptor'].map(split_subunits)
    df = df[(df['lig'].map(len) > 0) & (df['rec'].map(len) > 0)]
    df = df.drop_duplicates(subset=['lig', 'rec'], keep='first').reset_index(drop=True)
    if df.empty:
        raise ValueError('no usable ligand–receptor pairs')

    genes: list[str] = list(dict.fromkeys(g for side in ('lig', 'rec') for tup in df[side] for g in tup))
    gene_ix = {g: i for i, g in enumerate(genes)}
    complexes: list[tuple[str, ...]] = list(dict.fromkeys(t for side in ('lig', 'rec') for t in df[side] if len(t) > 1))
    complex_ix = {t: len(genes) + i for i, t in enumerate(complexes)}
    emap = ensembl_map if ensembl_map is not None else load_ensembl_map(species)

    gene_table = pd.DataFrame({
        'id_gene': range(len(genes)), 'ensembl': [emap.get(g, '') for g in genes], 'gene_name': genes,
        'hgnc_symbol': genes, 'protein_name': genes, 'protein_id': range(len(genes)),
    })
    protein_table = pd.DataFrame({'id_protein': range(len(genes)), 'protein_name': genes,
                                  'protein_multidata_id': range(len(genes))})
    complex_table = pd.DataFrame({'id_complex': range(len(complexes)),
                                  'complex_multidata_id': [complex_ix[t] for t in complexes],
                                  'name': ['_'.join(t) for t in complexes]})
    comp_rows = [(complex_ix[t], gene_ix[g], len(t)) for t in complexes for g in t]
    composition = pd.DataFrame(comp_rows, columns=['complex_multidata_id', 'protein_multidata_id', 'total_protein'])
    composition.insert(0, 'id_complex_composition', range(len(composition)))

    def mid(t):
        return complex_ix[t] if len(t) > 1 else gene_ix[t[0]]
    if 'id' in df.columns:
        ids = df['id'].astype(str).tolist()
    else:
        ids = ['_'.join(l) + '__' + '_'.join(r) for l, r in zip(df['lig'], df['rec'])]
    if len(set(ids)) != len(ids):
        raise ValueError('interaction ids are not unique')
    inter = pd.DataFrame({'id_interaction': range(len(df)), 'id_cp_interaction': ids,
                          'multidata_1_id': [mid(t) for t in df['lig']], 'multidata_2_id': [mid(t) for t in df['rec']]})
    for col in ('source', 'directionality', 'classification'):
        if col in df.columns:
            inter[col] = df[col].where(df[col].notna(), None).values

    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    inter.to_csv(out / 'interaction_table.csv', index=False)
    gene_table.to_csv(out / 'gene_table.csv', index=False)
    protein_table.to_csv(out / 'protein_table.csv', index=False)
    complex_table.to_csv(out / 'complex_table.csv', index=False)
    composition.to_csv(out / 'complex_composition_table.csv', index=False)
    write_metadata(out, name=name, species=species, source=source, version=version, citation=citation,
                   license=license, url=url, notes=notes)
    return out


def pairs_from_table(path: str | os.PathLike, ligand_col: str = 'ligand', receptor_col: str = 'receptor',
                     id_col: Optional[str] = None, classification_col: Optional[str] = None,
                     directionality_col: Optional[str] = None, source_col: Optional[str] = None,
                     sep: Optional[str] = None) -> pd.DataFrame:
    """Read a CSV/TSV of pairs and rename the chosen columns to the ones :func:`build_database` expects."""
    path = Path(path)
    if sep is None:
        sep = '\t' if path.suffix.lower() in ('.tsv', '.txt') else ','
    t = pd.read_csv(path, sep=sep)
    missing = [c for c in (ligand_col, receptor_col) if c not in t.columns]
    if missing:
        raise ValueError(f"columns {missing} not in {path.name}; available: {list(t.columns)[:20]}")
    out = pd.DataFrame({'ligand': t[ligand_col].astype(str), 'receptor': t[receptor_col].astype(str)})
    for new, col in (('id', id_col), ('classification', classification_col),
                     ('directionality', directionality_col), ('source', source_col)):
        if col:
            if col not in t.columns:
                raise ValueError(f"column {col!r} not in {path.name}")
            out[new] = t[col].values
    return out


def pairs_from_liana(csv_path: str | os.PathLike, resource: str) -> pd.DataFrame:
    """Pairs of one resource of LIANA's ``omni_resource.csv`` (columns source_genesymbol, target_genesymbol, resource)."""
    t = pd.read_csv(csv_path, index_col=0)
    sub = t[t['resource'] == resource]
    if sub.empty:
        raise ValueError(f"resource {resource!r} not in {csv_path}; available: {sorted(t['resource'].unique())}")
    return pd.DataFrame({'ligand': sub['source_genesymbol'].astype(str).values,
                         'receptor': sub['target_genesymbol'].astype(str).values})


# ── metadata ─────────────────────────────────────────────────────────────────

def database_info(path: str | os.PathLike) -> dict:
    """Metadata of a database directory: ``database.json`` when present, else counts computed from the tables."""
    d = Path(path).expanduser()
    meta: dict = {'name': d.name, 'species': 'human'}
    f = d / METADATA_FILE
    if f.is_file():
        try:
            meta.update(json.loads(f.read_text(encoding='utf-8')))
        except json.JSONDecodeError:
            pass
    if 'n_interactions' not in meta and (d / 'interaction_table.csv').is_file():
        inter = pd.read_csv(d / 'interaction_table.csv')
        meta['n_interactions'] = int(len(inter))
        meta['has_classification'] = _has_classification(inter)
        if (d / 'gene_table.csv').is_file():
            meta['n_genes'] = int(pd.read_csv(d / 'gene_table.csv')['hgnc_symbol'].nunique())
        if (d / 'complex_table.csv').is_file():
            meta['n_complexes'] = int(len(pd.read_csv(d / 'complex_table.csv')))
    return meta


def database_species(path: str | os.PathLike) -> str:
    return str(database_info(path).get('species') or 'human')


def write_metadata(path: str | os.PathLike, **fields) -> dict:
    """Write or update ``database.json`` of a database directory; counts are recomputed from the tables."""
    d = Path(path).expanduser()
    meta = {k: v for k, v in database_info(d).items() if k in METADATA_KEYS}
    meta.update({k: v for k, v in fields.items() if v is not None})
    if meta.get('species') not in SPECIES:
        raise ValueError(f"species must be one of {SPECIES}, got {meta.get('species')!r}")
    inter = pd.read_csv(d / 'interaction_table.csv')
    genes = pd.read_csv(d / 'gene_table.csv')
    comp = pd.read_csv(d / 'complex_table.csv')
    meta.update({'name': meta.get('name') or d.name, 'n_interactions': int(len(inter)),
                 'n_genes': int(genes['hgnc_symbol'].nunique()), 'n_complexes': int(len(comp)),
                 'has_classification': _has_classification(inter)})
    meta.setdefault('gene_id', 'symbol')
    meta.setdefault('built', _dt.date.today().isoformat())
    ordered = {k: meta.get(k) for k in METADATA_KEYS}
    ordered.update({k: v for k, v in meta.items() if k not in ordered})
    (d / METADATA_FILE).write_text(json.dumps(ordered, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    return ordered


def database_catalog(include_user: bool = True) -> list[dict]:
    """Metadata of every available database (packaged, then user-imported), with ``origin`` and ``path``."""
    rows = []
    roots = [(PACKAGED_DB_DIR, 'packaged')] + ([(user_database_dir(), 'user')] if include_user else [])
    for root, origin in roots:
        if not root.is_dir():
            continue
        for p in sorted(root.iterdir()):
            if _is_database_dir(p):
                rows.append({**database_info(p), 'origin': origin, 'path': str(p.resolve())})
    return rows


# ── validate ─────────────────────────────────────────────────────────────────

def validate_database(path: str | os.PathLike) -> dict:
    """Check a database directory.

    Returns ``{'ok', 'errors', 'warnings', 'summary'}``; ``ok`` is False when any error is found.
    """
    d = Path(path).expanduser()
    errors: list[str] = []
    warnings: list[str] = []
    summary: dict = {}
    if not d.is_dir():
        return {'ok': False, 'errors': [f'{d} is not a directory'], 'warnings': warnings, 'summary': summary}
    tables = {}
    for fname, cols in TABLES.items():
        f = d / fname
        if not f.is_file():
            errors.append(f'missing {fname}')
            continue
        t = pd.read_csv(f)
        missing = [c for c in cols if c not in t.columns]
        if missing:
            errors.append(f'{fname}: missing columns {missing}')
        tables[fname] = t
    if errors:
        return {'ok': False, 'errors': errors, 'warnings': warnings, 'summary': summary}
    inter, genes, prots = tables['interaction_table.csv'], tables['gene_table.csv'], tables['protein_table.csv']
    comp, compo = tables['complex_table.csv'], tables['complex_composition_table.csv']
    protein_mids = set(prots['protein_multidata_id'])
    complex_mids = set(comp['complex_multidata_id'])
    both = protein_mids & complex_mids
    if both:
        errors.append(f'{len(both)} multidata ids are both a protein and a complex')
    used = set(inter['multidata_1_id']) | set(inter['multidata_2_id'])
    unknown = used - protein_mids - complex_mids
    if unknown:
        errors.append(f'{len(unknown)} interaction multidata ids are neither a protein nor a complex (e.g. {sorted(unknown)[:5]})')
    orphan = set(compo['complex_multidata_id']) - complex_mids
    if orphan:
        errors.append(f'{len(orphan)} complex compositions refer to unknown complexes')
    empty = complex_mids - set(compo['complex_multidata_id'])
    if empty:
        errors.append(f'{len(empty)} complexes have no subunits')
    bad_sub = set(compo['protein_multidata_id']) - protein_mids
    if bad_sub:
        errors.append(f'{len(bad_sub)} complex subunits are not proteins')
    linked = set(genes['protein_id'])
    unmapped = set(prots['id_protein']) - linked
    if unmapped:
        errors.append(f'{len(unmapped)} proteins have no gene (e.g. {sorted(unmapped)[:5]})')
    sym = genes['hgnc_symbol']
    if sym.isna().any() or (sym.astype(str).str.strip() == '').any():
        errors.append('gene_table has empty hgnc_symbol entries')
    if inter['id_cp_interaction'].isna().any():
        errors.append('interaction_table has empty interaction ids')
    ndup = int(inter['id_cp_interaction'].duplicated().sum())
    if ndup:
        errors.append(f'{ndup} duplicated interaction ids')
    dup_pairs = int(inter.duplicated(subset=['multidata_1_id', 'multidata_2_id']).sum())
    if dup_pairs:
        warnings.append(f'{dup_pairs} duplicated ligand–receptor pairs (only the first is scored)')
    unused = (protein_mids | complex_mids) - used - set(compo['protein_multidata_id'])
    if unused:
        warnings.append(f'{len(unused)} proteins or complexes take part in no interaction')
    if not (d / METADATA_FILE).is_file():
        warnings.append(f'no {METADATA_FILE}: species and provenance unknown (treated as human)')
    meta = database_info(d)
    summary = {'name': meta.get('name'), 'species': meta.get('species'), 'n_interactions': int(len(inter)),
               'n_genes': int(sym.nunique()), 'n_complexes': int(len(comp)), 'has_classification': _has_classification(inter),
               'ensembl_filled': float(genes['ensembl'].notna().mean()) if len(genes) else 0.0}
    return {'ok': not errors, 'errors': errors, 'warnings': warnings, 'summary': summary}


# ── import ───────────────────────────────────────────────────────────────────

def import_database(path: str | os.PathLike, name: Optional[str] = None, dest_root: Optional[str | os.PathLike] = None,
                    *, overwrite: bool = False, pairs_kwargs: Optional[dict] = None, **fields) -> Path:
    """Install a database under ``$FASTCCC_HOME/databases/<name>`` and return its directory.

    ``path`` is either a database directory (copied after validation) or a CSV/TSV of
    ligand–receptor pairs (converted with :func:`build_database`; ``pairs_kwargs`` chooses its
    columns, see :func:`pairs_from_table`). ``fields`` (species, source, version, citation,
    license, url, notes) go into ``database.json``.
    """
    src = Path(path).expanduser()
    dest_root = Path(dest_root).expanduser() if dest_root else user_database_dir()
    if src.is_dir():
        name = name or src.name
        _check_name(name)
        report = validate_database(src)
        if not report['ok']:
            raise ValueError(f"{src} is not a valid database: " + '; '.join(report['errors']))
        dest = dest_root / name
        if dest.exists():
            if not overwrite:
                raise FileExistsError(f"{dest} exists; pass overwrite=True to replace it")
            shutil.rmtree(dest)
        dest_root.mkdir(parents=True, exist_ok=True)
        shutil.copytree(src, dest)
        write_metadata(dest, name=name, **fields)
        return dest.resolve()
    if src.is_file():
        if not name:
            raise ValueError('a name is required when importing a table of pairs')
        pairs = pairs_from_table(src, **(pairs_kwargs or {}))
        return build_database(pairs, name, dest_root, overwrite=overwrite, **fields).resolve()
    raise FileNotFoundError(f'{src} is neither a database directory nor a table of pairs')


def gene_coverage(path: str | os.PathLike, var_names: Iterable[str]) -> float:
    """Fraction of the database's genes present in ``var_names`` (quick fit check for a dataset).

    Symbols are compared with ``hgnc_symbol``; when most of ``var_names`` look like Ensembl ids
    (``ENS...``) the database's ``ensembl`` column is used instead.
    """
    names = set(map(str, var_names))
    col = 'ensembl' if names and sum(n.startswith('ENS') for n in names) > len(names) / 2 else 'hgnc_symbol'
    genes = pd.read_csv(Path(path) / 'gene_table.csv')[col].dropna().astype(str)
    genes = genes[genes.str.strip() != ''].drop_duplicates()
    return float(genes.isin(names).mean()) if len(genes) else 0.0
