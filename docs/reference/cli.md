---
title: Command line
layout: default
parent: Reference
nav_order: 1
description: Every fastccc command and option, generated from the built-in help.
---

# Command line reference

Generated from `fastccc <command> --help`; every command also accepts `--json` for a machine-readable summary.

```text
Usage: fastccc [OPTIONS] COMMAND [ARGS]...

  FastCCC: permutation-free, reference-based cell-cell communication analysis.

Options:
  -V, --version  Show the FastCCC version and exit.
  --help         Show this message and exit.

Commands:
  inspect    Summarise an .h5ad: matrix type, raw-count availability,...
  run        Standard FastCCC analysis of one dataset (optionally per group...
  ref        Reference-based comparison: which communications are Up or...
  build-ref  Build a custom reference panel from a dataset (optionally from...
  panels     List reference panels (local and downloadable), show one...
  map        Suggest how the query cell types map onto a reference panel's...
  explore    Build the interactive result explorer: one self-contained HTML...
  report     (Re)generate the HTML report for an existing result directory.
  db         Ligand-receptor databases: list, show, validate, import.
  skill      Show or install the FastCCC agent skill (for Claude Code,...
```

## `fastccc inspect`

```text
Usage: fastccc inspect [OPTIONS] {data}

  Summarise an .h5ad: matrix type, raw-count availability, candidate cell-type
  columns.

Arguments:
  data  Input .h5ad file.  [required]

Options:
  -c, --celltype <str>  Cell-type column to summarise.
  --json                Print machine-readable JSON.
  --help                Show this message and exit.
```

## `fastccc run`

```text
Usage: fastccc run [OPTIONS] {data}

  Standard FastCCC analysis of one dataset (optionally per group with --split-
  by).

Arguments:
  data  Input .h5ad file (raw counts or log1p).  [required]

Options:
  -c, --celltype <str>          Cell-type column in adata.obs.
  --celltype-file <path>        TSV with cell id and cell type columns.
  --db <str>                    Packaged LRI database name or path.  [default:
                                CPDBv5.0.0]
  -o, --out <path>              Output directory.
  --name <str>                  Sample display name.
  --split-by <str>              obs column; run each group separately.
  -g, --group <str>             Group value to include (repeat). Default: all.
  --no-all                      With --split-by: skip the full-dataset run.
  --min-cells <int>             With --split-by: skip groups with fewer cells
                                than this (recorded in the manifest).  [default:
                                20]
  --gene-id <str>               What adata.var_names are: 'symbol' (HGNC/MGI
                                gene symbols) or 'ensembl' (Ensembl gene ids;
                                human databases carry them).  [default: symbol]
  --small-p <str>               How extreme module p-values enter the Cauchy
                                combination: acat (default: p = 1 truncated to 1
                                - 1/16, p < 1e-15 as 1/(pi*p)) or legacy (the
                                published rule: p <= 1e-16 treated as 0 through
                                tan).  [default: acat]
  --no-report                   Skip building the interactive result explorer.
  --legacy-report               Also render the older static HTML report.
  --no-cells                    Explorer: do not embed per-cell expression
                                (smaller file).
  --bootstrap <int>             Re-run on N stratified cell subsamples and
                                record each interaction's bootstrap support (0 =
                                off; 10 is reasonable).  [default: 0]
  --bootstrap-fraction <float>  Fraction of every cell type kept in each
                                bootstrap subsample.  [default: 0.8]
  --input-type <str>            What adata.X holds: raw | normalized | log1p
                                (overrides auto-detection).
  --summary <str>               Expression summaries to score (repeat): Mean,
                                Median, Q2, Q3 or Quantile_<q> such as
                                Quantile_0.8. Default: Mean, Median, Q3,
                                Quantile_0.9.
  --aggregation <str>           Complex aggregations to score (repeat): Minimum,
                                Average. Default: both.
  --combination <str>           Ligand-receptor combinations to score (repeat):
                                Arithmetic, Geometric. Default: both.
  --dpi <int>                   Report figure resolution.  [default: 150]
  --json                        Print a JSON summary as the last line.
  --help                        Show this message and exit.
```

## `fastccc ref`

```text
Usage: fastccc ref [OPTIONS] {data}

  Reference-based comparison: which communications are Up or Down versus a
  reference.

Arguments:
  data  Query .h5ad file (counts preferred; normalised or log1p data give the
        same result).  [required]

Options:
  -r, --reference <str>    Tissue name (e.g. liver), panel directory, or
                           COLUMN:CONTROL[:CASE] to build the panel from this
                           dataset's control cells.  [required]
  -c, --celltype <str>     Cell-type column in adata.obs.
  --celltype-file <path>   TSV with cell id and cell type columns.
  -m, --mapping <path>     JSON {reference cell type: query cell type} (see
                           `fastccc map`).
  --auto-map               Use the automatic cell-type mapping suggestion
                           (written to celltype_mapping.json in the output).
  --db <str>               Packaged LRI database name or path.  [default:
                           CPDBv5.0.0]
  -o, --out <path>         Output directory.
  --name <str>             Query display name.
  --reference-name <str>   Reference display name.
  --reference-root <path>  Folder holding downloaded tissue panels.
  --no-report              Skip building the interactive result explorer.
  --legacy-report          Also render the older static HTML report.
  --no-cells               Explorer: do not embed per-cell expression (smaller
                           file).
  --dpi <int>              Report figure resolution.  [default: 150]
  --json                   Print a JSON summary as the last line.
  --help                   Show this message and exit.
```

## `fastccc map`

```text
Usage: fastccc map [OPTIONS] {data}

  Suggest how the query cell types map onto a reference panel's cell types
  (review before use).

Arguments:
  data  Query .h5ad file.  [required]

Options:
  -r, --reference <str>    Tissue name (e.g. blood) or panel directory.
                           [required]
  -c, --celltype <str>     Cell-type column in adata.obs.
  --celltype-file <path>   TSV with cell id and cell type columns.
  -o, --out <path>         Write the mapping JSON here (edit it, then pass to
                           `fastccc ref -m`).
  --min-score <float>      Minimum similarity (0-1) for a reference type to be
                           mapped.  [default: 0.5]
  --reference-root <path>  Folder holding downloaded tissue panels.
  --json                   Print the full suggestion as JSON.
  --help                   Show this message and exit.
```

## `fastccc build-ref`

```text
Usage: fastccc build-ref [OPTIONS] {data}

  Build a custom reference panel from a dataset (optionally from a subset of its
  cells).

Arguments:
  data  Reference .h5ad file (counts preferred; log-normalised works too).
        [required]

Options:
  -n, --name <str>          Panel name (becomes the folder name).  [required]
  --subset <str>            COLUMN:VALUE; build the panel from these cells only
                            (e.g. condition:control).
  -c, --celltype <str>      Cell-type column in adata.obs.
  --celltype-file <path>    TSV with cell id and cell type columns.
  -o, --out <path>          Folder under which <name>/ is created.  [default:
                            fastccc_reference]
  --db <str>                Packaged LRI database name or path.  [default:
                            CPDBv5.0.0]
  --min-percentile <float>  [default: 0.1]
  --json                    Print machine-readable JSON.
  --help                    Show this message and exit.
```

## `fastccc panels`

```text
Usage: fastccc panels [OPTIONS]

  List reference panels (local and downloadable), show one panel's cell types
  (downloads it if needed), or download panels.

Options:
  -s, --show <str>      Tissue name or panel directory: print its cell types and
                        cell counts (a tissue panel is downloaded first if
                        needed).
  -d, --download <str>  Tissue to download, or 'all'.
  --root <path>         Reference root to list / download into.
  --json                Print machine-readable JSON.
  --help                Show this message and exit.
```

## `fastccc explore`

```text
Usage: fastccc explore [OPTIONS] {result_dir}

  Build the interactive result explorer: one self-contained HTML page with
  filters, table, evidence and figures.

Arguments:
  result_dir  A FastCCC result folder (standard run, --split-by run, or
              reference comparison).  [required]

Options:
  -o, --out <path>   Output HTML file (default: <result_dir>/explorer.html).
  --no-both-ns       Reference comparisons: leave out 'Both NS' rows to shrink
                     the file.
  --no-cells         Do not embed per-cell expression (violins, embedding
                     overlay).
  --max-cells <int>  Cells kept for the per-cell panels (stratified subsample).
                     [default: 4000]
  --no-genesets      Skip the pathway enrichment libraries (no download).
  --patterns <int>   Number of communication patterns (default: automatic).
  --json
  --help             Show this message and exit.
```

## `fastccc report`

```text
Usage: fastccc report [OPTIONS] {result_dir}

  (Re)generate the HTML report for an existing result directory.

Arguments:
  result_dir  A FastCCC result directory (standard, split or reference).
              [required]

Options:
  -o, --out <path>  Report output directory.
  --dpi <int>       Report figure resolution.  [default: 150]
  --json            Print machine-readable JSON.
  --help            Show this message and exit.
```

## `fastccc db list`

```text
Usage: fastccc db list [OPTIONS]

  List the available databases (packaged and imported) with species, size and
  pathway annotation.

Options:
  --data <path>  An .h5ad file: also report the fraction of each database's
                 genes present in it (low values mean a species or gene-id
                 mismatch).
  --json         Print machine-readable JSON.
  --help         Show this message and exit.
```

## `fastccc db show`

```text
Usage: fastccc db show [OPTIONS] {name}

  Show one database: path, species, size, provenance, citation and license.

Arguments:
  name  Database name or directory.  [required]

Options:
  --json  Print machine-readable JSON.
  --help  Show this message and exit.
```

## `fastccc db validate`

```text
Usage: fastccc db validate [OPTIONS] {path}

  Check a database directory for the tables, columns and cross-references
  FastCCC needs.

Arguments:
  path  A database directory (five CSV tables).  [required]

Options:
  --json  Print machine-readable JSON.
  --help  Show this message and exit.
```

## `fastccc db import`

```text
Usage: fastccc db import [OPTIONS] {source}

  Install a database so that `--db <name>` finds it (converted from a pair table
  if needed).

Arguments:
  source  A database directory in FastCCC/CellPhoneDB format, or a CSV/TSV with
          one ligand-receptor pair per row (complex subunits joined by '_' or
          ',').  [required]

Options:
  -n, --name <str>            Name used with --db (default: the directory name).
  --species <str>             human, mouse, rat, zebrafish or other.  [default:
                              human]
  --ligand-col <str>          Pair table: ligand column.  [default: ligand]
  --receptor-col <str>        Pair table: receptor column.  [default: receptor]
  --id-col <str>              Pair table: unique interaction id column (default:
                              built from the gene names).
  --classification-col <str>  Pair table: pathway/family column (shown in the
                              explorer).
  --source-col <str>          Pair table: evidence/source column.
  --source <str>              Provenance of the database (e.g. 'CellChatDB v2').
  --version <str>             Version of the source database.
  --citation <str>            Reference to cite.
  --license <str>             License of the source database.
  --url <str>                 Where the source lives.
  --dest <path>               Install into this directory instead of the FastCCC
                              home.
  --overwrite                 Replace an existing database of the same name.
  --json                      Print machine-readable JSON.
  --help                      Show this message and exit.
```

## `fastccc skill install`

```text
Usage: fastccc skill install [OPTIONS]

  Copy the skill into the agents' skill folders so `/fastccc` (Claude Code) or
  `$fastccc` (Codex) works.

Options:
  -a, --agent <str>  claude, codex, cursor, opencode, gemini, copilot or all
                     (repeatable). Default: claude and codex.
  -g, --global       Install for the user (e.g. ~/.claude/skills) instead of the
                     current project.
  --dir <path>       Project directory for a project-level install (default:
                     current directory).
  --force            Overwrite an existing installation.
  --json
  --help             Show this message and exit.
```

## `fastccc skill show`

```text
Usage: fastccc skill show [OPTIONS]

  Print the packaged SKILL.md (what an agent reads to operate FastCCC).

Options:
  --path  Print the packaged SKILL.md path instead of its content.
  --help  Show this message and exit.
```
