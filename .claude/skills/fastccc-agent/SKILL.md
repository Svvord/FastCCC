---
name: fastccc-agent
description: Use when the user wants to run FastCCC through natural language instead of manually wiring Python calls, paths, and parameters. This skill helps an agent choose the right FastCCC workflow, collect the minimum required inputs, run analysis or reporting steps, and explain outputs for single-dataset analysis, multi-condition comparison, reference building, reference-based inference, and result interpretation.
---

# FastCCC Agent

Use this skill when the user wants an agent to operate FastCCC directly from plain-language requests.

## What this skill does

- Maps user intent to one FastCCC workflow
- Collects only the required inputs
- Chooses safe defaults when the user does not care
- Runs the relevant FastCCC entrypoints
- Explains where outputs were written and what they mean

Read [references/workflows.md](references/workflows.md) before executing a workflow.

## Workflow selection

Choose exactly one primary workflow first:

- Single dataset CCC analysis: user has one `.h5ad` or `AnnData` and wants CCC results
- Multi-condition comparison: user wants per-condition runs plus one differential report
- Reference-based inference: user wants to compare a query dataset to a tissue reference panel
- Reference building: user wants to create a new FastCCC reference panel from raw-count reference data
- Result explanation: user already has FastCCC outputs and wants interpretation or a report

If the request is ambiguous, ask one short blocking question.

## Required inputs by workflow

Single dataset CCC analysis:
- dataset path or in-memory `AnnData`
- cell type column name, or a metadata TSV

Multi-condition comparison:
- dataset path
- condition column name
- exact values for condition A and condition B in that column
- cell type column name
- whether to include only two conditions or also run the full combined dataset

Reference-based inference:
- query dataset path (raw counts — FastCCC rank-preprocesses internally)
- reference panel path
- cell type column name, or a metadata TSV
- `celltype_mapping_dict` if query cell type names differ from those in the reference (JSON file or dict mapping reference name → query name)

Reference building:
- raw-count reference dataset path
- reference name
- output directory
- cell type column name, or a metadata TSV

Result explanation:
- result directory
- task id (for standard FastCCC outputs) — or leave empty for reference-based inference outputs
- database path if report or annotation is needed
- clarify whether results came from standard FastCCC or reference-based inference, as output formats differ

## Safe defaults

Use these defaults unless the user specifies otherwise:

- Database: `db/CPDBv5.0.0`
- Gene identifier type: `hgnc_symbol`
- Standard analysis entrypoint: `fastccc.Cauchy_combination_of_statistical_analysis_methods(...)`
- Cell type key preference order: `cell_type`, `celltype`, `CellType`, `annotation`, `cell_label` — inspect `adata.obs.columns` and pick the best match; FastCCC does not auto-detect this
- Report generation: run `fastccc.report.generate_report(...)` after successful analysis
- Gene sets for reports: `KEGG_2021_Human`, `GO_Biological_Process_2023`
- Report DPI: 300 (use 150 for faster test runs)
- Report p-value threshold: 0.05
- Report top N LR pairs shown: 30
- Report top N cell types shown: 20

## Input safety rules

- Do not assume raw counts and normalized `log1p` counts are interchangeable.
- For standard FastCCC analysis, prefer normalized `log1p` data.
- For reference building and reference-based inference, use raw-count data because FastCCC rank-preprocesses those workflows internally.
- If the data scale is unclear and using the wrong assumption would change results materially, ask one short question before running.
- Do not silently guess the wrong cell type or condition column if inspection shows multiple plausible choices.

## Execution rules

- Inspect the dataset metadata first when possible.
- Prefer using the project's packaged functions rather than rewriting analysis logic.
- Create result directories clearly, usually under `results/` or a user-provided output directory.
- After standard FastCCC analysis, identify the produced `task_id` from `*_significant_results.tsv`.
- For reference-based inference, there is no `task_id`; look for `query_infer_results.tsv` instead.
- If the user asked for interpretation, summarize:
  - output files created
  - number of cell types
  - number of significant interactions
  - obvious dominant sender, receiver, or pathway patterns if a report or result table exists

## Response style

When using this skill, the agent should:

- state the chosen workflow in one sentence
- state any assumptions that matter
- ask only the minimum blocking questions
- return exact output paths
- suggest the next useful action, such as generating a report or interpreting the top interactions
