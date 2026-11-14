# FastCCC workflows: complete examples

Use after `SKILL.md`. Every command has a `--json` flag that prints a machine-readable
summary as its last line, and every output folder contains `fastccc_run.json`.

## 0. Inspect a dataset

```bash
fastccc inspect data.h5ad            # human readable
fastccc inspect data.h5ad --json     # for parsing
```

Key fields: `matrix.kind` and `matrix.confidence`, `raw_counts_available`,
`celltype_candidates` (column, n_unique, score, reasons), `suggested_celltype_key`,
`suggested_celltype_key_confidence`, `celltype_counts`, `gene_id_type`.

```python
import fastccc
info = fastccc.inspect('data.h5ad')
info['suggested_celltype_key'], info['matrix']['kind']
```

## 1. One dataset

User says: "run cell-cell communication on this h5ad", "which cell types talk to each other".

```bash
fastccc run data.h5ad -c cell_type -o results/sample --name "Sample 1"
```

```python
res = fastccc.run('data.h5ad', 'cell_type', out_dir='results/sample', name='Sample 1')
res.summary()          # n_celltypes, n_significant, top senders / receivers / LR pairs
res.significant.head() # sender_celltype, receiver_celltype, ligand, receptor, p-value
res.report_path        # HTML report (written automatically; report=False to skip)
```

Options: `--celltype-file labels.tsv` (cell id, cell type) instead of `-c`;
`--input-type raw|normalized|log1p` to override the matrix detection; `--no-report`;
`--bootstrap 10` to re-run on 10 stratified cell subsamples and record each
interaction's bootstrap support (a robustness measure; adds ~10x the run time);
`--dpi 300`; `--db <name or path>` for another CellPhoneDB-format database.

Runtime: about 1 min per 10k cells with 20-40 cell types on 4 cores; the per-method
intermediate tables take several GB for large datasets.

## 2. Several conditions in one dataset

User says: "compare treated vs control", "run it per condition", "disease vs healthy
(descriptive)". For a statistical Up/Down call against the control, use workflow 4 instead.

```bash
fastccc run cohort.h5ad -c cell_type --split-by condition -o results/by_condition
fastccc run cohort.h5ad -c cell_type --split-by condition -g treated -g control --no-all -o results/pair
```

- Every group gets its own folder (`results/by_condition/<group>/`) plus `all/` for the
  full dataset unless `--no-all`.
- With exactly two groups and the full-dataset run, the report gains a
  condition-comparison tab (shared, gained and lost interactions).
- With three or more groups (samples, donors) the explorer adds, for every
  interaction of the full dataset, the fraction of samples in which it is
  significant on its own (replication), plus a samples x interactions matrix:
  use `--split-by sample` when the user asks whether findings replicate across
  donors.

```python
multi = fastccc.run('cohort.h5ad', 'cell_type', split_by='condition', out_dir='results/by_condition')
multi.groups            # ['control', 'treated']
multi['treated'].summary()
multi.all.summary()
multi.report_path
```

## 3. Compare against a built-in healthy tissue panel

User says: "which interactions are up or down compared with normal liver", "use the
healthy reference", "what is abnormal in this tumour sample".

```bash
fastccc panels                                   # which panels exist / are downloaded
fastccc panels --show liver                      # the panel's cell types and cell counts (downloads it if needed)
fastccc map tumour.h5ad -c cell_type -r liver -o liver_map.json
# review liver_map.json (see SKILL.md), then
fastccc ref tumour.h5ad -c cell_type -r liver -m liver_map.json -o results/tumour_vs_liver --name "Tumour"
```

```python
sug = fastccc.suggest_celltype_mapping('tumour.h5ad', 'cell_type', reference='liver')
print(sug.format())                 # per query cell type: pooled panel types, scores
sug.ambiguous, sug.unmatched_query, sug.unmatched_reference   # decide these, then edit sug.mapping
mapping = dict(sug.mapping); mapping['T cell'] = 'T cells'      # example decision
q = fastccc.compare_to_reference('tumour.h5ad', 'cell_type', reference='liver',
                                 celltype_mapping=mapping, out_dir='results/tumour_vs_liver')
q.trend_counts                      # {'Up': ..., 'Down': ..., 'Both Sig': ..., 'Both NS': ..., 'Not in reference': ...}
q.interaction_summary('Down')       # per LR pair: n sender-receiver pairs, strongest pair, expression query vs reference
q.interaction_summary('Up')
q.interactions('Down')              # the raw rows
q.celltype_mapping                  # mapping actually used
q.summary()['celltype_coverage']
```

Panels are downloaded on first use into `~/.fastccc/reference` (change with the
`FASTCCC_HOME` environment variable or `--reference-root`). A panel is "activated"
(null distributions computed, ~30 s) the first time it is used.

Mapping format (`-m`): JSON object, keys = panel cell types, values = query cell
types. Several panel types may map to one query type; they are pooled weighted by
their cell numbers. Panel types absent from the JSON are used only when a query type
has exactly the same name.

## 4. Control group as the reference (same dataset)

User says: "use the controls as reference", "disease vs control, which communications
are gained or lost", "treated vs untreated with a proper null".

```bash
fastccc ref cohort.h5ad -c cell_type -r condition:control:disease -o results/disease_vs_control
fastccc ref cohort.h5ad -c cell_type -r condition:control -o results/all_vs_control   # query = every non-control cell
```

The control cells become a custom panel under `results/disease_vs_control/reference/`,
and the disease cells are the query. Labels match by name, so no mapping is needed.

```python
q = fastccc.compare_to_reference('cohort.h5ad', 'cell_type', reference=('condition', 'control', 'disease'),
                                 out_dir='results/disease_vs_control')
```

## 5. Build a reusable reference panel

User says: "build a reference from this atlas", "make a panel from the healthy donors".

```bash
fastccc build-ref atlas.h5ad -c cell_type -n healthy_pancreas -o panels
fastccc build-ref atlas.h5ad -c cell_type -n healthy_pancreas -o panels --subset status:healthy
fastccc ref query.h5ad -c cell_type -r panels/healthy_pancreas -m map.json -o results/query_vs_pancreas
```

```python
panel = fastccc.build_reference_panel('atlas.h5ad', 'cell_type', name='healthy_pancreas',
                                      out_dir='panels', subset=('status', 'healthy'))
panel.celltypes, panel.n_cells, panel.path
```

Whole-transcriptome input (all genes, not only highly variable genes) gives better
panels. Counts, normalised or log1p data all work (rank-based).

## 6. Explorer and reports from existing results

```bash
fastccc explore results/sample                 # interactive explorer.html (default output of run / ref)
fastccc explore results/tumour_vs_liver --no-cells --max-cells 2000
fastccc report results/sample                  # older static HTML report (also: run/ref --legacy-report)
```

```python
fastccc.load_result('results/sample').report(dpi=300)
```

## 7. Interpreting outputs

Standard run (`<task>` is a 6-character id stored in `fastccc_run.json`):

| File | Content |
|---|---|
| `<task>_significant_results.tsv` | sender_celltype, receiver_celltype, LRI_ID, ligand, receptor, p-value (p < 0.05) |
| `<task>_Cauchy_pvals.tsv` | p-values, rows `sender\|receiver`, columns LR pairs |
| `<task>_average_interactions_strength.tsv` | communication scores, same shape |
| `<task>_<timestamp>_<method>_*.tsv` | 16 per-method intermediates; ignore |

Reference comparison:

Rows: every sender|receiver x LR pair expressed in the query (ligand and receptor
each in at least `min_percentile` = 10% of the sender / receiver cells), plus the
reference-significant pairs the query does not express (`above_expr_threshold` False,
always `Down`).

| Column of `query_infer_results.tsv` | Meaning |
|---|---|
| `sender\|receiver`, `ligand`, `receptor`, `LRI_ID` | the interaction |
| `in_reference` | both cell types exist in the reference (after mapping) |
| `comm_score` | communication score in the query |
| `null_comm_score`, `sig_threshold_CI` | null expectation and the significance band derived from the reference |
| `above_expr_threshold`, `above_expr_threshold_ref` | passes the expression filter in the query / reference |
| `ligand_expr_percent`, `receptor_expr_percent` (+ `_ref`) | fraction of expressing cells in the query / reference |
| `ligand_CS_component`, `receptor_CS_component` (+ `_CI`) | ligand and receptor parts of the score, with the reference interval |
| `is_significant`, `is_significant_ref` | significant in the query / in the reference |
| `trend_vs_ref` | `Up`, `Down`, `Both Sig`, `Both NS`; empty when the cell types are not in the reference |

Report the Up and Down counts and the top LR pairs from
`q.interaction_summary('Up')` / `q.interaction_summary('Down')` (per LR pair: number
of sender-receiver pairs, strongest pair, expression percentages in query versus
reference), and which cell types were not covered. `Down` means the reference has the
communication and the query does not.

## 8. Legacy entrypoints (still available)

Older scripts use `fastccc.Cauchy_combination_of_statistical_analysis_methods(...)`,
`fastccc.build_reference.build_reference_workflow(...)`,
`fastccc.infer_query.infer_query_workflow(...)`, `fastccc.report.generate_report(...)`
and `fastccc.report.generate_reference_report(...)`. They are unchanged and produce the
same tables; the entrypoints above wrap them and add input checks, the mapping tools
and `fastccc_run.json`. Prefer the new ones for new work.
