"""``fastccc`` command-line interface (typer)."""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import List, Optional

import typer

app = typer.Typer(
    help="FastCCC: permutation-free, reference-based cell-cell communication analysis.",
    no_args_is_help=True, add_completion=False, rich_markup_mode=None,
)


def _print_json(obj) -> None:
    typer.echo(json.dumps(obj, default=str))


def _call(func, *args, **kwargs):
    """Run an API function; turn expected input errors into a short message and exit code 1."""
    try:
        return func(*args, **kwargs)
    except (ValueError, KeyError, FileNotFoundError) as err:
        msg = err.args[0] if err.args else str(err)
        typer.echo(f"error: {msg}", err=True)
        raise typer.Exit(code=1)


def _version_callback(value: bool):
    if value:
        from . import __version__
        typer.echo(f"fastccc {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Optional[bool] = typer.Option(None, '--version', '-V', callback=_version_callback, is_eager=True,
                                           help="Show the FastCCC version and exit."),
):
    """FastCCC command line. Run ``fastccc <command> --help`` for details."""
    from . import log_to
    log_to(sys.stderr)          # progress messages on stderr; stdout carries results and --json


@app.command()
def inspect(
    data: Path = typer.Argument(..., help="Input .h5ad file."),
    celltype: Optional[str] = typer.Option(None, '--celltype', '-c', help="Cell-type column to summarise."),
    json_out: bool = typer.Option(False, '--json', help="Print machine-readable JSON."),
):
    """Summarise an .h5ad: matrix type, raw-count availability, candidate cell-type columns."""
    from .datacheck import inspect_data
    info = _call(inspect_data, str(data), celltype)
    if json_out:
        _print_json(info)
        return
    m = info['matrix']
    typer.echo(f"{data}: {info['n_cells']:,} cells x {info['n_genes']:,} genes")
    typer.echo(f"adata.X looks like: {m['kind']} (confidence {m['confidence']}: {m['reason']})")
    evidence = f"  max={m['max_value']:.2f}, integer fraction={m['integer_fraction']:.3f}, per-cell total CV={m['rowsum_cv']:.3f}"
    if m.get('expm1_rowsum_cv') is not None:
        evidence += f", after expm1={m['expm1_rowsum_cv']:.3f}"
    typer.echo(evidence)
    typer.echo(f"raw counts available: {info['raw_counts_available']}" + (f" (from {info['raw_counts_source']})" if info['raw_counts_available'] else ''))
    typer.echo(f"gene ids: {info['gene_id_type']}")
    typer.echo("candidate cell-type columns (score, unique values):")
    for c in info['celltype_candidates'][:8]:
        why = '; '.join(c['reasons'][:2])
        typer.echo(f"  {c['column']:<30s} {c['score']:>5.1f} {c['n_unique']:>5d}  {why}")
    if info['suggested_celltype_key']:
        typer.echo(f"suggested celltype_key: {info['suggested_celltype_key']} (confidence {info['suggested_celltype_key_confidence']})")
        counts = info.get('celltype_counts', {})
        if counts:
            shown = ', '.join(f"{k} ({v})" for k, v in list(counts.items())[:10])
            typer.echo(f"  {shown}{' ...' if len(counts) > 10 else ''}")
    else:
        typer.echo("suggested celltype_key: none found; pass -c <column> explicitly")
    typer.echo(f"standard workflow input: {info['recommended_workflow_input']['standard']}")
    typer.echo(f"reference workflow input: {info['recommended_workflow_input']['reference']}")


@app.command()
def run(
    data: Path = typer.Argument(..., help="Input .h5ad file (raw counts or log1p)."),
    celltype: Optional[str] = typer.Option(None, '--celltype', '-c', help="Cell-type column in adata.obs."),
    celltype_file: Optional[Path] = typer.Option(None, '--celltype-file', help="TSV with cell id and cell type columns."),
    db: str = typer.Option('CPDBv5.0.0', '--db', help="Packaged LRI database name or path."),
    out: Optional[Path] = typer.Option(None, '--out', '-o', help="Output directory."),
    name: Optional[str] = typer.Option(None, '--name', help="Sample display name."),
    split_by: Optional[str] = typer.Option(None, '--split-by', help="obs column; run each group separately."),
    group: Optional[List[str]] = typer.Option(None, '--group', '-g', help="Group value to include (repeat). Default: all."),
    no_all: bool = typer.Option(False, '--no-all', help="With --split-by: skip the full-dataset run."),
    min_cells: int = typer.Option(20, '--min-cells', help="With --split-by: skip groups with fewer cells than this (recorded in the manifest)."),
    gene_id: str = typer.Option('symbol', '--gene-id', help="What adata.var_names are: 'symbol' (HGNC/MGI gene symbols) or 'ensembl' (Ensembl gene ids; human databases carry them)."),
    small_p: str = typer.Option('acat', '--small-p', help="How extreme module p-values enter the Cauchy combination: acat (default: p = 1 truncated to 1 - 1/16, p < 1e-15 as 1/(pi*p)) or legacy (the published rule: p <= 1e-16 treated as 0 through tan)."),
    no_report: bool = typer.Option(False, '--no-report', help="Skip building the interactive result explorer."),
    legacy_report: bool = typer.Option(False, '--legacy-report', help="Also render the older static HTML report."),
    no_cells: bool = typer.Option(False, '--no-cells', help="Explorer: do not embed per-cell expression (smaller file)."),
    bootstrap: int = typer.Option(0, '--bootstrap', help="Re-run on N stratified cell subsamples and record each interaction's bootstrap support (0 = off; 10 is reasonable)."),
    bootstrap_fraction: float = typer.Option(0.8, '--bootstrap-fraction', help="Fraction of every cell type kept in each bootstrap subsample."),
    input_type: Optional[str] = typer.Option(None, '--input-type', help="What adata.X holds: raw | normalized | log1p (overrides auto-detection)."),
    summary: Optional[List[str]] = typer.Option(None, '--summary', help="Expression summaries to score (repeat): Mean, Median, Q2, Q3 or Quantile_<q> such as Quantile_0.8. Default: Mean, Median, Q3, Quantile_0.9."),
    aggregation: Optional[List[str]] = typer.Option(None, '--aggregation', help="Complex aggregations to score (repeat): Minimum, Average. Default: both."),
    combination: Optional[List[str]] = typer.Option(None, '--combination', help="Ligand-receptor combinations to score (repeat): Arithmetic, Geometric. Default: both."),
    dpi: int = typer.Option(150, '--dpi', help="Report figure resolution."),
    json_out: bool = typer.Option(False, '--json', help="Print a JSON summary as the last line."),
):
    """Standard FastCCC analysis of one dataset (optionally per group with --split-by)."""
    from .api import run as _run
    modules = {}
    import re
    def _is_summary(v):
        m = re.fullmatch(r'Quantile_(0?\.\d+)', v)
        return v in ('Mean', 'Median', 'Q2', 'Q3') or (m is not None and 0 < float(m.group(1)) < 1)
    for key, values, allowed in (('single_unit_summary_list', summary, ('Mean', 'Median', 'Q2', 'Q3', 'Quantile_<q>')),
                                 ('complex_aggregation_list', aggregation, ('Minimum', 'Average')),
                                 ('LR_combination_list', combination, ('Arithmetic', 'Geometric'))):
        if values:
            bad = [v for v in values if not (_is_summary(v) if key == 'single_unit_summary_list' else v in allowed)]
            if bad:
                flag = {'single_unit_summary_list': '--summary', 'complex_aggregation_list': '--aggregation', 'LR_combination_list': '--combination'}[key]
                raise typer.BadParameter(f"{', '.join(bad)} not in {', '.join(allowed)}", param_hint=flag)
            modules[key] = list(dict.fromkeys(values))
    if gene_id not in ('symbol', 'ensembl'):
        typer.echo(f"error: --gene-id must be 'symbol' or 'ensembl', got {gene_id!r}", err=True)
        raise typer.Exit(code=1)
    res = _call(_run, str(data), celltype, celltype_file=celltype_file, database=db, out_dir=out, name=name, convert_type='ensembl' if gene_id == 'ensembl' else 'hgnc_symbol',
                split_by=split_by, groups=group or None, include_all=not no_all, min_cells=min_cells, small_p_method=small_p, report=not no_report,
                legacy_report=legacy_report, explorer_kwargs={'cells': not no_cells},
                bootstrap=bootstrap, bootstrap_fraction=bootstrap_fraction,
                input_type=input_type, report_kwargs={'dpi': dpi}, **modules)
    _finish(res, json_out)


@app.command()
def ref(
    data: Path = typer.Argument(..., help="Query .h5ad file (counts preferred; normalised or log1p data give the same result)."),
    reference: str = typer.Option(..., '--reference', '-r', help="Tissue name (e.g. liver), panel directory, or COLUMN:CONTROL[:CASE] to build the panel from this dataset's control cells."),
    celltype: Optional[str] = typer.Option(None, '--celltype', '-c', help="Cell-type column in adata.obs."),
    celltype_file: Optional[Path] = typer.Option(None, '--celltype-file', help="TSV with cell id and cell type columns."),
    mapping: Optional[Path] = typer.Option(None, '--mapping', '-m', help="JSON {reference cell type: query cell type} (see `fastccc map`)."),
    auto_map: bool = typer.Option(False, '--auto-map', help="Use the automatic cell-type mapping suggestion (written to celltype_mapping.json in the output)."),
    db: str = typer.Option('CPDBv5.0.0', '--db', help="Packaged LRI database name or path."),
    out: Optional[Path] = typer.Option(None, '--out', '-o', help="Output directory."),
    name: Optional[str] = typer.Option(None, '--name', help="Query display name."),
    reference_name: Optional[str] = typer.Option(None, '--reference-name', help="Reference display name."),
    reference_root: Optional[Path] = typer.Option(None, '--reference-root', help="Folder holding downloaded tissue panels."),
    no_report: bool = typer.Option(False, '--no-report', help="Skip building the interactive result explorer."),
    legacy_report: bool = typer.Option(False, '--legacy-report', help="Also render the older static HTML report."),
    no_cells: bool = typer.Option(False, '--no-cells', help="Explorer: do not embed per-cell expression (smaller file)."),
    dpi: int = typer.Option(150, '--dpi', help="Report figure resolution."),
    json_out: bool = typer.Option(False, '--json', help="Print a JSON summary as the last line."),
):
    """Reference-based comparison: which communications are Up or Down versus a reference."""
    from .api import compare_to_reference
    ref_arg = reference
    if ':' in reference and not Path(reference).exists():
        ref_arg = tuple(reference.split(':'))
    if mapping and auto_map:
        raise typer.BadParameter("use either --mapping or --auto-map, not both")
    res = _call(compare_to_reference, str(data), celltype, reference=ref_arg, celltype_file=celltype_file,
                celltype_mapping='auto' if auto_map else (str(mapping) if mapping else None),
                database=db, out_dir=out,
                name=name, reference_name=reference_name, reference_root=reference_root,
                report=not no_report, legacy_report=legacy_report, explorer_kwargs={'cells': not no_cells},
                report_kwargs={'dpi': dpi})
    _finish(res, json_out)


@app.command('build-ref')
def build_ref(
    data: Path = typer.Argument(..., help="Reference .h5ad file (counts preferred; log-normalised works too)."),
    name: str = typer.Option(..., '--name', '-n', help="Panel name (becomes the folder name)."),
    subset: Optional[str] = typer.Option(None, '--subset', help="COLUMN:VALUE; build the panel from these cells only (e.g. condition:control)."),
    celltype: Optional[str] = typer.Option(None, '--celltype', '-c', help="Cell-type column in adata.obs."),
    celltype_file: Optional[Path] = typer.Option(None, '--celltype-file', help="TSV with cell id and cell type columns."),
    out: Path = typer.Option(Path('fastccc_reference'), '--out', '-o', help="Folder under which <name>/ is created."),
    db: str = typer.Option('CPDBv5.0.0', '--db', help="Packaged LRI database name or path."),
    min_percentile: float = typer.Option(0.1, '--min-percentile'),
    json_out: bool = typer.Option(False, '--json', help="Print machine-readable JSON."),
):
    """Build a custom reference panel from a dataset (optionally from a subset of its cells)."""
    from .api import build_reference_panel
    panel = _call(build_reference_panel, str(data), celltype, name=name, out_dir=out, celltype_file=celltype_file,
                  database=db, min_percentile=min_percentile, subset=subset)
    if json_out:
        _print_json(panel.to_dict())
    else:
        typer.echo(f"Reference panel '{panel.name}' built: {panel.n_celltypes} cell types, {panel.n_cells:,} cells -> {panel.path}")


@app.command()
def panels(
    show: Optional[str] = typer.Option(None, '--show', '-s', help="Tissue name or panel directory: print its cell types and cell counts (a tissue panel is downloaded first if needed)."),
    download: Optional[str] = typer.Option(None, '--download', '-d', help="Tissue to download, or 'all'."),
    root: Optional[Path] = typer.Option(None, '--root', help="Reference root to list / download into."),
    json_out: bool = typer.Option(False, '--json', help="Print machine-readable JSON."),
):
    """List reference panels (local and downloadable), show one panel's cell types (downloads it if needed), or download panels."""
    from . import resources
    from .results import ReferencePanel
    if show:
        info = _call(resources.resolve_reference, show, reference_root=root, download=True)
        panel = ReferencePanel.load(info['path'], source=info['source'], tissue=info['tissue'])
        _print_json(panel.to_dict()) if json_out else typer.echo(panel.describe())
        return
    if download:
        tissues = resources.REFERENCE_TISSUES if download == 'all' else [download]
        for t in tissues:
            resources.download_reference_panel(t, dest_root=root)
    local = resources.list_reference_panels(root)
    rows = []
    for t in resources.REFERENCE_TISSUES:
        path = resources.find_reference_panel(t, root)
        if path is not None:
            p = ReferencePanel.load(path, source='healthy_tissue', tissue=t)
            rows.append({'tissue': t, 'status': 'local', 'path': str(path), 'n_celltypes': p.n_celltypes,
                         'n_cells': p.n_cells, 'database': p.database, 'activated': p.is_activated})
        else:
            rows.append({'tissue': t, 'status': 'downloadable', 'path': None})
    for t in local:
        if t not in resources.REFERENCE_TISSUES:
            path = resources.find_reference_panel(t, root)
            p = ReferencePanel.load(path)
            rows.append({'tissue': t, 'status': 'local (custom)', 'path': str(path), 'n_celltypes': p.n_celltypes,
                         'n_cells': p.n_cells, 'database': p.database, 'activated': p.is_activated})
    if json_out:
        _print_json({'reference_home': str(resources.fastccc_home() / 'reference'), 'panels': rows})
        return
    typer.echo(f"reference cache: {resources.fastccc_home() / 'reference'}  (override with FASTCCC_HOME)")
    for r in rows:
        extra = f"  {r['n_celltypes']:>3d} cell types, {r['n_cells']:>9,d} cells, {r['database']}" if r['status'].startswith('local') else ''
        typer.echo(f"  {r['tissue']:<18s} {r['status']:<15s}{extra}")


db_app = typer.Typer(help="Ligand-receptor databases: list, show, validate, import.", no_args_is_help=True, rich_markup_mode=None)
app.add_typer(db_app, name='db')


@db_app.command('list')
def db_list(
    data: Optional[Path] = typer.Option(None, '--data', help="An .h5ad file: also report the fraction of each database's genes present in it (low values mean a species or gene-id mismatch)."),
    json_out: bool = typer.Option(False, '--json', help="Print machine-readable JSON."),
):
    """List the available databases (packaged and imported) with species, size and pathway annotation."""
    from . import resources
    from .dbtools import database_catalog, gene_coverage
    rows = database_catalog()
    if data is not None:
        import anndata
        var_names = _call(lambda: anndata.read_h5ad(str(data), backed='r').var_names)
        for r in rows:
            r['gene_coverage'] = round(gene_coverage(r['path'], var_names), 3)
    if json_out:
        _print_json({'user_database_dir': str(resources.user_database_dir()), 'databases': rows})
        return
    typer.echo(f"user databases: {resources.user_database_dir()}  (override with FASTCCC_HOME)")
    head = f"  {'name':<24s} {'species':<8s} {'interactions':>12s} {'genes':>6s} {'complexes':>9s}  {'pathways':<8s} {'origin':<8s}"
    if data is not None:
        head += f" {'in data':>7s}"
    typer.echo(head)
    for r in rows:
        line = (f"  {r['name']:<24s} {str(r.get('species') or '?'):<8s} {r.get('n_interactions', 0):>12,d} {r.get('n_genes', 0):>6,d} "
                f"{r.get('n_complexes', 0):>9,d}  {('yes' if r.get('has_classification') else 'no'):<8s} {r['origin']:<8s}")
        if data is not None:
            line += f" {r['gene_coverage']:>6.0%}"
        typer.echo(line)
    typer.echo(f"default: {resources.DEFAULT_DATABASE}; use another with --db <name>. `fastccc db show <name>` prints provenance and citation.")


@db_app.command('show')
def db_show(
    name: str = typer.Argument(..., help="Database name or directory."),
    json_out: bool = typer.Option(False, '--json', help="Print machine-readable JSON."),
):
    """Show one database: path, species, size, provenance, citation and license."""
    from .dbtools import database_info
    from .resources import resolve_database
    path = _call(resolve_database, name)
    meta = {**database_info(path), 'path': str(path)}
    if json_out:
        _print_json(meta)
        return
    for k in ('name', 'path', 'species', 'version', 'n_interactions', 'n_genes', 'n_complexes', 'has_classification',
              'source', 'url', 'citation', 'license', 'notes', 'built'):
        if meta.get(k) not in (None, ''):
            typer.echo(f"{k:<18s} {meta[k]}")


@db_app.command('validate')
def db_validate(
    path: Path = typer.Argument(..., help="A database directory (five CSV tables)."),
    json_out: bool = typer.Option(False, '--json', help="Print machine-readable JSON."),
):
    """Check a database directory for the tables, columns and cross-references FastCCC needs."""
    from .dbtools import validate_database
    rep = validate_database(path)
    if json_out:
        _print_json(rep)
    else:
        s = rep['summary']
        if s:
            typer.echo(f"{path}: {s['n_interactions']:,} interactions, {s['n_genes']:,} genes, {s['n_complexes']:,} complexes, "
                       f"species {s['species']}, pathway classification: {'yes' if s['has_classification'] else 'no'}")
        for e in rep['errors']:
            typer.echo(f"  error: {e}")
        for w in rep['warnings']:
            typer.echo(f"  warning: {w}")
        typer.echo('OK' if rep['ok'] else 'INVALID')
    if not rep['ok']:
        raise typer.Exit(code=1)


@db_app.command('import')
def db_import(
    source: Path = typer.Argument(..., help="A database directory in FastCCC/CellPhoneDB format, or a CSV/TSV with one ligand-receptor pair per row (complex subunits joined by '_' or ',')."),
    name: Optional[str] = typer.Option(None, '--name', '-n', help="Name used with --db (default: the directory name)."),
    species: str = typer.Option('human', '--species', help="human, mouse, rat, zebrafish or other."),
    ligand_col: str = typer.Option('ligand', '--ligand-col', help="Pair table: ligand column."),
    receptor_col: str = typer.Option('receptor', '--receptor-col', help="Pair table: receptor column."),
    id_col: Optional[str] = typer.Option(None, '--id-col', help="Pair table: unique interaction id column (default: built from the gene names)."),
    classification_col: Optional[str] = typer.Option(None, '--classification-col', help="Pair table: pathway/family column (shown in the explorer)."),
    source_col: Optional[str] = typer.Option(None, '--source-col', help="Pair table: evidence/source column."),
    provenance: Optional[str] = typer.Option(None, '--source', help="Provenance of the database (e.g. 'CellChatDB v2')."),
    version: Optional[str] = typer.Option(None, '--version', help="Version of the source database."),
    citation: Optional[str] = typer.Option(None, '--citation', help="Reference to cite."),
    license_: Optional[str] = typer.Option(None, '--license', help="License of the source database."),
    url: Optional[str] = typer.Option(None, '--url', help="Where the source lives."),
    dest: Optional[Path] = typer.Option(None, '--dest', help="Install into this directory instead of the FastCCC home."),
    overwrite: bool = typer.Option(False, '--overwrite', help="Replace an existing database of the same name."),
    json_out: bool = typer.Option(False, '--json', help="Print machine-readable JSON."),
):
    """Install a database so that `--db <name>` finds it (converted from a pair table if needed)."""
    from .dbtools import import_database, database_info
    pairs_kwargs = dict(ligand_col=ligand_col, receptor_col=receptor_col, id_col=id_col,
                        classification_col=classification_col, source_col=source_col)
    try:
        path = import_database(source, name=name, dest_root=dest, overwrite=overwrite, pairs_kwargs=pairs_kwargs,
                               species=species, source=provenance, version=version, citation=citation, license=license_, url=url)
    except (ValueError, FileNotFoundError, FileExistsError) as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1)
    meta = database_info(path)
    if json_out:
        _print_json({**meta, 'path': str(path)})
    else:
        typer.echo(f"database '{meta['name']}' installed: {meta['n_interactions']:,} interactions, {meta['n_genes']:,} genes, "
                   f"{meta['n_complexes']:,} complexes ({meta['species']}) -> {path}")
        typer.echo(f"use it with: fastccc run <data.h5ad> -c <celltype> --db {meta['name']}")


@app.command('map')
def map_celltypes(
    data: Path = typer.Argument(..., help="Query .h5ad file."),
    reference: str = typer.Option(..., '--reference', '-r', help="Tissue name (e.g. blood) or panel directory."),
    celltype: Optional[str] = typer.Option(None, '--celltype', '-c', help="Cell-type column in adata.obs."),
    celltype_file: Optional[Path] = typer.Option(None, '--celltype-file', help="TSV with cell id and cell type columns."),
    out: Optional[Path] = typer.Option(None, '--out', '-o', help="Write the mapping JSON here (edit it, then pass to `fastccc ref -m`)."),
    min_score: float = typer.Option(0.5, '--min-score', help="Minimum similarity (0-1) for a reference type to be mapped."),
    reference_root: Optional[Path] = typer.Option(None, '--reference-root', help="Folder holding downloaded tissue panels."),
    json_out: bool = typer.Option(False, '--json', help="Print the full suggestion as JSON."),
):
    """Suggest how the query cell types map onto a reference panel's cell types (review before use)."""
    from .api import suggest_celltype_mapping
    sug = _call(suggest_celltype_mapping, str(data), celltype, reference=reference, celltype_file=celltype_file,
                reference_root=reference_root, min_score=min_score)
    if out:
        sug.to_json(out)
    if json_out:
        _print_json({**sug.to_dict(), 'written_to': str(out) if out else None})
        return
    typer.echo(f"query cell types -> reference panel '{sug.reference_name}' cell types")
    typer.echo(sug.format())
    typer.echo('')
    typer.echo(sug.summary().splitlines()[0])
    if sug.unmatched_reference:
        typer.echo(f"reference cell types left unused ({len(sug.unmatched_reference)}, best similarity below --min-score {min_score}):")
        for name in sug.unmatched_reference:
            typer.echo(f"    {name}")
    if out:
        typer.echo(f"mapping written to {out}  (edit if needed, then: fastccc ref {data} -c {celltype or '<column>'} -r {reference} -m {out})")


skill_app = typer.Typer(help="Show or install the FastCCC agent skill (for Claude Code, Codex, Cursor, ...).",
                        rich_markup_mode=None, no_args_is_help=True)
app.add_typer(skill_app, name='skill')


@skill_app.command('show')
def skill_show(
    path_only: bool = typer.Option(False, '--path', help="Print the packaged SKILL.md path instead of its content."),
):
    """Print the packaged SKILL.md (what an agent reads to operate FastCCC)."""
    from .skill import skill_path, skill_text
    typer.echo(str(skill_path()) if path_only else skill_text())


@skill_app.command('install')
def skill_install(
    agent: Optional[List[str]] = typer.Option(None, '--agent', '-a', help="claude, codex, cursor, opencode, gemini, copilot or all (repeatable). Default: claude and codex."),
    global_: bool = typer.Option(False, '--global', '-g', help="Install for the user (e.g. ~/.claude/skills) instead of the current project."),
    project_dir: Optional[Path] = typer.Option(None, '--dir', help="Project directory for a project-level install (default: current directory)."),
    force: bool = typer.Option(False, '--force', help="Overwrite an existing installation."),
    json_out: bool = typer.Option(False, '--json'),
):
    """Copy the skill into the agents' skill folders so `/fastccc` (Claude Code) or `$fastccc` (Codex) works."""
    from .skill import install_skill, DEFAULT_AGENTS
    records = _call(install_skill, agent or DEFAULT_AGENTS, global_=global_, project_dir=project_dir, force=force)
    if json_out:
        _print_json(records)
        return
    for r in records:
        note = {'installed': 'installed', 'updated': 'updated', 'exists': 'already present (use --force to overwrite)'}[r['status']]
        typer.echo(f"{', '.join(r['agents']):<22s} {r['path']}  [{note}]")
    typer.echo("Alternative without pip: npx skills add Svvord/FastCCC")


@app.command()
def explore(
    result_dir: Path = typer.Argument(..., help="A FastCCC result folder (standard run, --split-by run, or reference comparison)."),
    out: Optional[Path] = typer.Option(None, '--out', '-o', help="Output HTML file (default: <result_dir>/explorer.html)."),
    no_both_ns: bool = typer.Option(False, '--no-both-ns', help="Reference comparisons: leave out 'Both NS' rows to shrink the file."),
    no_cells: bool = typer.Option(False, '--no-cells', help="Do not embed per-cell expression (violins, embedding overlay)."),
    max_cells: int = typer.Option(4000, '--max-cells', help="Cells kept for the per-cell panels (stratified subsample)."),
    no_genesets: bool = typer.Option(False, '--no-genesets', help="Skip the pathway enrichment libraries (no download)."),
    patterns: Optional[int] = typer.Option(None, '--patterns', help="Number of communication patterns (default: automatic)."),
    json_out: bool = typer.Option(False, '--json'),
):
    """Build the interactive result explorer: one self-contained HTML page with filters, table, evidence and figures."""
    from .explorer import build_explorer
    path = _call(build_explorer, result_dir, output=out, include_both_ns=not no_both_ns, cells=not no_cells, max_cells=max_cells,
                 genesets=not no_genesets, n_patterns=patterns)
    if json_out:
        _print_json({'explorer_html': str(path)})
    else:
        typer.echo(f"explorer: {path}")


@app.command()
def report(
    result_dir: Path = typer.Argument(..., help="A FastCCC result directory (standard, split or reference)."),
    out: Optional[Path] = typer.Option(None, '--out', '-o', help="Report output directory."),
    dpi: int = typer.Option(150, '--dpi', help="Report figure resolution."),
    json_out: bool = typer.Option(False, '--json', help="Print machine-readable JSON."),
):
    """(Re)generate the HTML report for an existing result directory."""
    from .results import load_result, MultiRunResult
    res = load_result(result_dir)
    if isinstance(res, MultiRunResult):
        from .api import report_split
        report_split(res, output_dir=out, dpi=dpi)
    else:
        res.report(output_dir=out, dpi=dpi)
    _finish(res, json_out)


def _finish(res, json_out: bool) -> None:
    if json_out:
        _print_json(res.to_dict())
        return
    typer.echo(repr(res))
    cov = (getattr(res, 'manifest', None) or {}).get('celltype_coverage')
    if cov:
        src = (res.manifest or {}).get('celltype_mapping_source', 'none')
        line = (f"cell types: {cov['n_covered']}/{cov['n_query_types']} query cell types covered by the reference "
                f"({cov['cells_covered_fraction']:.1%} of cells); mapping: {src}")
        if cov['uncovered_query_types']:
            line += f"; not in reference: {', '.join(cov['uncovered_query_types'])}"
        typer.echo(line)
    ep = getattr(res, 'explorer_path', None)
    if ep:
        typer.echo(f"explorer: {ep}   (open in a browser)")
    rp = getattr(res, 'report_path', None)
    if rp:
        typer.echo(f"report: {rp}")


if __name__ == '__main__':   # pragma: no cover
    app()
