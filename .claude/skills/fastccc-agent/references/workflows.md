# FastCCC Workflows

Use this file after `SKILL.md` triggers.

## 1. Single dataset CCC analysis

Use when the user says things like:

- "帮我跑这个 h5ad 的细胞通讯"
- "Run FastCCC on this dataset"
- "Analyze CCC and give me a report"

Preferred flow:

1. Inspect the `.h5ad` metadata and candidate cell type columns.
2. Confirm whether the matrix is already normalized `log1p` if that is not obvious.
3. Run `fastccc.Cauchy_combination_of_statistical_analysis_methods(...)`.
4. Detect the `task_id` from `*_significant_results.tsv`.
5. Generate a report with `fastccc.report.generate_report(...)` unless the user only wants raw tables.

Minimum parameters:

- `database_file_path`
- `counts_file_path`
- `meta_key` or `celltype_file_path`
- `save_path`

## 2. Multi-condition comparison

Use when the user wants A/B comparison, disease vs control, treated vs untreated, or timepoint comparison.

Preferred flow:

1. Inspect the dataset and identify `condition` and `cell type` columns.
2. Decide which runs to produce:
   - **Full four-tab report** (most common): run FastCCC for all cells + condition A + condition B, then call `generate_report` with all three result dirs.
   - **Two-tab report** (no differential): run FastCCC for all cells + condition A only; pass only `cond_a_*` to `generate_report` and omit `cond_b_*`.
3. Run each condition in a separate result directory. Use `ProcessPoolExecutor` for parallel runs (see `examples/report_multi_condition.py`).
4. Cache already-completed runs by checking for `*_significant_results.tsv` before re-running.
5. Generate one report with `generate_report(...)`.

Use the project example as the shape of the workflow:

- `examples/report_multi_condition.py`

Needed user choices:

- which column defines the condition
- exact values for condition A and B
- whether to also run the full combined dataset as the primary tab

## 3. Reference-based inference

Use when the user says the query should be compared with a tissue reference panel.

Preferred flow:

1. Verify the reference path exists and that its `config.toml` LRI database matches the one you will use for the query.
2. Confirm the query data is raw counts (not log1p), as FastCCC rank-preprocesses internally.
3. If query cell type names differ from those in the reference, prepare a `celltype_mapping_dict` (JSON or dict, mapping reference cell type name → query cell type name).
4. Run `fastccc.infer_query.infer_query_workflow(...)`.
5. Summarize the output files listed below.

Important:

- These workflows expect raw-count input — do not use normalized data.
- Do not swap in the standard normalized workflow unless the user explicitly wants that instead.

Expected outputs (all in `save_path`):

- `query_infer_results.tsv` — main result table; key columns: `sender|receiver`, `ligand`, `receptor`, `comm_score`, `is_significant`, `is_significant_ref`, `trend_vs_ref` (`Up` / `Down` / `Both Sig` / `Both NS`)
- `query_interactions_strength.tsv` — full CS matrix (celltype_pair × LRI)
- `query_percents_analysis.tsv` — expression percentage filter matrix

There is no `task_id` for reference-based inference results.

## 4. Reference building

Use when the user wants to create a new panel from a large reference dataset.

Preferred flow:

1. Confirm the reference data is raw-count based.
2. Verify cell type labels.
3. Run `fastccc.build_reference.build_reference_workflow(...)`.
4. Return the created reference directory and key files.

Expected outputs (all in `save_path/<reference_name>/`):

- `config.toml` — reference metadata (name, LRI database, min_percentile, cell type counts)
- `ref_gene_pmf_dict.pkl` — per-gene null distributions (used at inference time)
- `ref_mean_counts.pkl` — per-cell-type mean expression
- `ref_percents.pkl` — per-cell-type expression percentages
- `complex_table.pkl` — complex composition table filtered to present genes
- `interactions.pkl` — LRI table filtered to present genes
- `ref_hk.txt` — housekeeping gene rank means (used for calibration factor k)

Note: `basic_info_dict.pkl` is only saved when `for_uploading=True` (used for panel upload, not normal local use).

## 5. Existing result interpretation

Use when the user already has FastCCC outputs.

First determine which type of results they have — the output formats are different:

### Standard FastCCC results (from workflows 1 or 2)

Locate the following files using the `task_id`:

- `{task_id}_significant_results.tsv` — significant LR pairs (filtered by p < 0.05)
- `{task_id}_Cauchy_pvals.tsv` — full p-value matrix (authoritative for Cauchy combination, the default)
- `{task_id}_average_interactions_strength.tsv` — mean CS across all method variants (authoritative strength file)

The individual method files (`{task_id}_{timestamp}_{method_key}_pvals.tsv` etc.) are intermediate results; use the Cauchy and average files for interpretation.

If a report does not exist, offer or generate one with `generate_report(...)`.

### Reference-based inference results (from workflow 3)

Locate:

- `query_infer_results.tsv` — main result; interpret via `trend_vs_ref` column (`Up` = stronger than reference, `Down` = weaker, `Both Sig` = significant in both, `Both NS` = non-significant in both)
- `query_interactions_strength.tsv` — full CS matrix

No task_id exists; no `generate_report` is available for these results.

### Useful interpretation points (both types)

- strongest sender-receiver cell type pairs by interaction count or CS
- top pathways or classifications
- interaction count density across cell type pairs
- differences between conditions if two runs are available

## 6. Failure handling

If the run fails, classify the problem before asking the user:

- missing file path
- missing `obs` column
- wrong matrix preprocessing (log1p used where raw counts needed, or vice versa)
- incompatible database or reference pairing (check `config.toml` LRI database field)
- empty or nearly empty overlap after gene mapping
- cell type name mismatch between query and reference (fix with `celltype_mapping_dict`)

Only ask the user for the missing fact that actually blocks the next step.
