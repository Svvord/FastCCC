---
title: Python API
layout: default
parent: Reference
nav_order: 2
description: The FastCCC Python API - run, compare_to_reference, build_reference_panel, suggest_celltype_mapping, inspect, load_result and the result objects.
---

# Python API

```python
import fastccc
```

Every function accepts a path to an `.h5ad` file or an in-memory `AnnData`. Keyword arguments after `*` are keyword-only.

## `fastccc.inspect(data, celltype_key=None) -> dict`

Matrix type with evidence and confidence, raw-count availability, scored cell-type column candidates, the suggested column and its counts, gene id type. JSON-serialisable; the same information as `fastccc inspect --json`.

## `fastccc.run(data, celltype_key=None, *, ...) -> RunResult | MultiRunResult`

| Argument | Default | Meaning |
|---|---|---|
| `celltype_key` | required (or `celltype_file`) | column of `adata.obs` with the cell types |
| `celltype_file` | `None` | TSV with cell id and cell type instead of a column |
| `database` | `'CPDBv5.0.0'` | packaged database name or path |
| `out_dir` | `fastccc_results/<name>` | output folder |
| `name` | file stem | display name |
| `split_by`, `groups`, `include_all` | `None`, all, `True` | per-group runs (see [Conditions and samples]({{ site.baseurl }}/guides/conditions.html)) |
| `min_cells` | `20` | with `split_by`: skip groups with fewer cells (listed as `skipped_groups` in the manifest) |
| `small_p_method` | `'acat'` | how extreme module p-values enter the Cauchy combination (`'acat'`: p = 1 truncated to 1 − 1/16, 1/(πp) below 1e-15; `'legacy'`: the published rule) |
| `report` | `True` | build `explorer.html` |
| `legacy_report` | `False` | also render the older static HTML report |
| `explorer_kwargs` | `None` | options for `explore()` (`cells`, `max_cells`, `genesets`, `patterns`, `n_patterns`) |
| `input_type` | detected | `'raw'`, `'normalized'` or `'log1p'` |
| `bootstrap`, `bootstrap_fraction`, `bootstrap_seed` | `0`, `0.8`, `0` | bootstrap support |
| `min_percentile` | `0.1` | expression filter |
| `single_unit_summary_list`, `complex_aggregation_list`, `LR_combination_list` | all 4, 2, 2 | which scoring modules to run |
| `use_DEG`, `select_list`, `filter_`, `convert_type` | as in the original function | |
| `copy` | `True` | never modify the caller's `AnnData` |

## `fastccc.compare_to_reference(data, celltype_key=None, *, reference, ...) -> ReferenceResult`

| Argument | Default | Meaning |
|---|---|---|
| `reference` | required | a tissue name (`'liver'`), a panel directory, or `(column, control)` / `(column, control, case)` to build the panel from a control group of the same dataset |
| `celltype_mapping` | `None` | `{reference cell type: query cell type}` dict, a JSON path, or `'auto'`; validated before the run |
| `query_relabel` | `None` | `{query label: new label}` applied first (merge sub-types) |
| `reference_root`, `download` | `None`, `True` | where panels live; whether to download a missing tissue panel |
| `out_dir`, `name`, `reference_name` | derived | output folder and display names |
| `report`, `legacy_report`, `explorer_kwargs` | as for `run` | |

Raises `ValueError` when no query cell type is covered by the panel, when the mapping has unknown names (with "did you mean" hints) or when the panel was built with another database.

## `fastccc.build_reference_panel(data, celltype_key=None, *, name, out_dir='fastccc_reference', subset=None, ...) -> ReferencePanel`

Builds `<out_dir>/<name>/`; `subset=('condition', 'control')` uses those cells only.

## `fastccc.suggest_celltype_mapping(query, celltype_key=None, *, reference, min_score=0.5, ...) -> MappingSuggestion`

`query` is a file, an `AnnData`, a list of labels or a `{label: n_cells}` dict; `reference` as above or a `ReferencePanel`. The result has `.mapping` (pass it as `celltype_mapping`), `.entries` (score and confidence per mapped panel type), `.ambiguous`, `.unmatched_query`, `.unmatched_reference`, `.coverage`, `.format()`, `.table()`, `.to_json(path)`.

## Panels and resources

```python
fastccc.list_reference_panels(include_downloadable=True)
fastccc.download_reference_panel("liver")
fastccc.resolve_reference("liver")              # {'path', 'source', 'tissue', 'root'}
fastccc.ReferencePanel.load(path)               # .celltypes, .n_cells, .n_celltypes, .database, .describe()
fastccc.list_databases(); fastccc.resolve_database("CPDBv5.0.0")   # packaged and imported databases
from fastccc import dbtools
dbtools.database_catalog()                      # metadata of every database (species, size, citation)
dbtools.build_database(pairs_df, "mydb", out_dir, species="mouse")   # ligand/receptor columns, subunits joined by "_"
dbtools.import_database("pairs.csv", name="mydb", species="mouse", pairs_kwargs={"classification_col": "pathway"})
dbtools.validate_database(path)                 # {"ok", "errors", "warnings", "summary"}
fastccc.fastccc_home()                          # ~/.fastccc or $FASTCCC_HOME
```

## Result objects

`fastccc.load_result(folder)` returns the right object for any result folder.

**RunResult** (one dataset): `result_dir`, `task_id`, `name`, `database`, `manifest`; lazy DataFrames `significant`, `pvals`, `strength`; `celltypes`; `summary()`; `explore(**options) -> Path`; `explorer_path`; `report() -> Path`; `report_path`; `to_dict()`.

**MultiRunResult** (`split_by`): `groups`, `runs` (dict), `all` (RunResult or `None`), `multi['treated']`, `explore()`, `explorer_path`, `to_dict()`.

**ReferenceResult**: `results` (DataFrame of `query_infer_results.tsv`), `strength`, `trend_counts`, `interactions(trend)`, `interaction_summary(trend, n=20)` (per ligand–receptor pair: number of cell-type pairs, strongest pair, expression percentages query vs reference), `summary()`, `celltype_mapping`, `celltype_coverage`, `reference` (ReferencePanel), `explore()`, `explorer_path`, `report()`.

## Explorer

```python
from fastccc.explorer import build_explorer
build_explorer("results/sample", output=None, include_both_ns=True, genesets=True, cells=True, max_cells=4000, patterns=True, n_patterns=None, nulls=True)
```

## Lower-level functions

The functions the one-call API wraps are unchanged and documented under [Legacy API]({{ site.baseurl }}/legacy): `fastccc.Cauchy_combination_of_statistical_analysis_methods`, `fastccc.build_reference.build_reference_workflow`, `fastccc.infer_query.infer_query_workflow`, `fastccc.report.generate_report`, `fastccc.report.generate_reference_report`. They produce the same tables; the new entrypoints add input checks, mapping tools and the run manifest.

## Distribution toolkit

```python
from fastccc.distrib import order_statistic_pmf, distribution_from_values, Distribution, get_pvalue_from_pmf
order_statistic_pmf([0, 1, 2], [0.5, 0.3, 0.2], n=10, k=5)   # exact pmf of the 5th of 10 draws, ties included
d = Distribution(samples=x); (d ** 30) / 30                   # mean of 30 i.i.d. cells; d1 + d2, d1 * d2, d1 & d2
get_pvalue_from_pmf(value, d)                                 # exact upper-tail probability
```

See the [Distribution toolkit]({{ site.baseurl }}/methods/distributions.html).
