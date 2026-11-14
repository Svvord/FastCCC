"""Locate packaged LRI databases and healthy-tissue reference panels.

Databases ship inside the wheel (``fastccc/data/db``); user-imported databases live in
``$FASTCCC_HOME/databases`` (see :mod:`fastccc.dbtools`).  Reference panels are
downloaded on first use into ``$FASTCCC_HOME/reference`` (default
``~/.fastccc/reference``); a source checkout's ``reference/`` folder is also
recognised so developers keep working offline.
"""
from __future__ import annotations

import os
import urllib.error
import urllib.request
from pathlib import Path
from typing import Optional

from loguru import logger

PACKAGE_DIR = Path(__file__).resolve().parent
PACKAGED_DB_DIR = PACKAGE_DIR / 'data' / 'db'
DEFAULT_DATABASE = 'CPDBv5.0.0'

REFERENCE_TISSUES = [
    'adipose', 'blood', 'bone_marrow', 'breast', 'cerebral_cortex', 'colon',
    'eye', 'heart', 'kidney', 'liver', 'lung', 'lymph_node', 'midbrain',
    'prostate_gland', 'skin', 'small_intestine', 'spleen', 'stomach', 'thymus',
]
REFERENCE_PANEL_FILES = [
    'config.toml', 'basic_info_dict.pkl', 'complex_table.pkl',
    'interactions.pkl', 'ref_hk.txt', 'ref_mean_counts.pkl', 'ref_percents.pkl',
]
DEFAULT_REFERENCE_URL = 'https://raw.githubusercontent.com/Svvord/FastCCC/main/reference'


# ── home directory ───────────────────────────────────────────────────────────

def fastccc_home() -> Path:
    """Writable cache directory (``$FASTCCC_HOME`` or ``~/.fastccc``)."""
    return Path(os.environ.get('FASTCCC_HOME', Path.home() / '.fastccc')).expanduser()


# ── LRI databases ────────────────────────────────────────────────────────────

def _is_database_dir(path: Path) -> bool:
    return path.is_dir() and (path / 'interaction_table.csv').is_file()


def user_database_dir() -> Path:
    """Directory of user-imported databases (``$FASTCCC_HOME/databases``)."""
    return fastccc_home() / 'databases'


def database_roots() -> list[Path]:
    """Directories searched for a database name: packaged first, then the user directory."""
    return [PACKAGED_DB_DIR, user_database_dir()]


def list_databases(include_user: bool = True) -> list[str]:
    """Names of the LRI databases shipped with FastCCC (and, by default, those imported by the user)."""
    roots = database_roots() if include_user else [PACKAGED_DB_DIR]
    names: list[str] = []
    for root in roots:
        if root.is_dir():
            names.extend(p.name for p in root.iterdir() if _is_database_dir(p))
    return sorted(dict.fromkeys(names))


def find_database(name: str) -> Optional[Path]:
    """Directory of a packaged or user-imported database, or ``None``."""
    for root in database_roots():
        if _is_database_dir(root / name):
            return (root / name).resolve()
    return None


def resolve_database(database: Optional[str] = None) -> Path:
    """Turn a database name (``'CPDBv5.0.0'``, ``'liana_consensus'``) or directory path into a path.

    Names are looked up among the packaged databases and then in
    ``$FASTCCC_HOME/databases``. A path that does not exist but whose last
    component is a known name falls back to that database, so legacy
    ``'./db/CPDBv5.0.0'`` arguments keep working after the data moved into the package.
    """
    if database is None:
        database = DEFAULT_DATABASE
    raw = str(database)
    candidate = Path(raw).expanduser()
    if _is_database_dir(candidate):
        return candidate.resolve()

    name = candidate.name or candidate.parent.name   # handles trailing slash
    found = find_database(name)
    if found is not None:
        if name != raw:
            logger.info(f"Database path '{raw}' not found; using database '{name}' from {found.parent}.")
        return found

    available = ', '.join(list_databases())
    raise FileNotFoundError(
        f"Unknown LRI database '{raw}'. Pass a database directory or one of these names: {available}. "
        f"`fastccc db list` shows them with species and size; `fastccc db import` adds your own."
    )


def legacy_database_path(database: str | os.PathLike) -> str:
    """Database path for the original entry points: an existing directory is used as given,
    a name or a former checkout path such as ``./db/CPDBv4.1.0/`` resolves to the packaged copy."""
    return str(resolve_database(str(database)))


def database_name(database: str | os.PathLike) -> str:
    return Path(str(database).rstrip('/\\')).name


# ── reference panels ─────────────────────────────────────────────────────────

def _is_panel_dir(path: Path) -> bool:
    return path.is_dir() and (path / 'config.toml').is_file()


def source_checkout_reference_dir() -> Optional[Path]:
    """``reference/`` next to the package when running from a git checkout."""
    candidate = PACKAGE_DIR.parent / 'reference'
    return candidate if candidate.is_dir() else None


def reference_roots(reference_root: Optional[str | os.PathLike] = None) -> list[Path]:
    """Directories searched for healthy tissue panels, in priority order."""
    roots: list[Path] = []
    if reference_root is not None:
        roots.append(Path(reference_root).expanduser())
    roots.append(fastccc_home() / 'reference')
    checkout = source_checkout_reference_dir()
    if checkout is not None:
        roots.append(checkout)
    cwd_ref = Path.cwd() / 'reference'
    if cwd_ref.is_dir() and cwd_ref not in roots:
        roots.append(cwd_ref)
    seen, ordered = set(), []
    for r in roots:
        key = str(r.resolve()) if r.exists() else str(r)
        if key not in seen:
            seen.add(key)
            ordered.append(r)
    return ordered


def list_reference_panels(
    reference_root: Optional[str | os.PathLike] = None,
    include_downloadable: bool = False,
) -> list[str]:
    """Healthy tissue panels available locally (optionally plus downloadable ones)."""
    names: set[str] = set()
    for root in reference_roots(reference_root):
        if root.is_dir():
            names.update(p.name for p in root.iterdir() if _is_panel_dir(p))
    if include_downloadable:
        names.update(REFERENCE_TISSUES)
    return sorted(names)


def find_reference_panel(
    tissue: str, reference_root: Optional[str | os.PathLike] = None,
) -> Optional[Path]:
    for root in reference_roots(reference_root):
        panel = root / tissue
        if _is_panel_dir(panel):
            return panel.resolve()
    return None


def download_reference_panel(
    tissue: str,
    dest_root: Optional[str | os.PathLike] = None,
    base_url: Optional[str] = None,
    overwrite: bool = False,
) -> Path:
    """Fetch one healthy tissue panel (≈5 MB) from GitHub into ``dest_root``."""
    if tissue not in REFERENCE_TISSUES:
        raise ValueError(
            f"Unknown tissue '{tissue}'. Downloadable panels: {', '.join(REFERENCE_TISSUES)}"
        )
    base_url = (base_url or os.environ.get('FASTCCC_REFERENCE_URL', DEFAULT_REFERENCE_URL)).rstrip('/')
    dest = Path(dest_root).expanduser() if dest_root else fastccc_home() / 'reference'
    panel_dir = dest / tissue
    if _is_panel_dir(panel_dir) and not overwrite:
        return panel_dir.resolve()
    panel_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Downloading FastCCC reference panel '{tissue}' to {panel_dir} …")
    for fname in REFERENCE_PANEL_FILES:
        url = f"{base_url}/{tissue}/{fname}"
        target = panel_dir / fname
        tmp = target.with_suffix(target.suffix + '.part')
        try:
            with urllib.request.urlopen(url, timeout=120) as resp, open(tmp, 'wb') as fh:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            tmp.replace(target)
        except urllib.error.URLError as exc:
            tmp.unlink(missing_ok=True)
            raise ConnectionError(
                f"Could not download {url}: {exc}. If this machine is offline, copy the "
                f"'reference/{tissue}' folder from the FastCCC repository into {dest}."
            ) from exc
    logger.success(f"Reference panel '{tissue}' is ready.")
    return panel_dir.resolve()


def resolve_reference(
    reference: str | os.PathLike,
    reference_root: Optional[str | os.PathLike] = None,
    download: bool = True,
) -> dict:
    """Resolve a tissue name or panel directory.

    Returns ``{'path', 'source', 'tissue', 'root'}`` where ``source`` is
    ``'healthy_tissue'`` for a shipped panel or ``'custom'`` for a user panel.
    """
    raw = str(reference)
    candidate = Path(raw).expanduser()
    if _is_panel_dir(candidate):
        path = candidate.resolve()
        for root in reference_roots(reference_root):
            if root.exists() and path.parent == root.resolve() and path.name in REFERENCE_TISSUES:
                return {'path': path, 'source': 'healthy_tissue', 'tissue': path.name, 'root': root.resolve()}
        return {'path': path, 'source': 'custom', 'tissue': None, 'root': path.parent}

    found = find_reference_panel(raw, reference_root)
    if found is not None:
        return {'path': found, 'source': 'healthy_tissue', 'tissue': raw, 'root': found.parent}

    if raw in REFERENCE_TISSUES:
        if not download:
            raise FileNotFoundError(
                f"Reference panel '{raw}' is not available locally and download=False."
            )
        dest = Path(reference_root).expanduser() if reference_root else None
        path = download_reference_panel(raw, dest_root=dest)
        return {'path': path, 'source': 'healthy_tissue', 'tissue': raw, 'root': path.parent}

    raise FileNotFoundError(
        f"'{raw}' is neither a reference panel directory nor a known tissue. "
        f"Tissues: {', '.join(REFERENCE_TISSUES)}."
    )
