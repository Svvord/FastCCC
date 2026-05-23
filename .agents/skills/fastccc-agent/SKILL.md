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
- Multi-condition comparison: user wants per-condition runs plus one condition-comparison report
- Reference-based inference: user wants to compare a query dataset to a healthy tissue panel or custom reference panel
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
- reference panel choice: either `reference_tissue` under a reference root, or a custom `reference_path`
- cell type column name, or a metadata TSV
- `celltype_mapping_dict` if query cell type names differ from those in the reference (JSON file or dict mapping reference name → query name)

Reference-based inference from two groups in one dataset:
- dataset path with a condition/group column
- reference/control group value used to build a custom reference panel
- query/case group value compared against that custom reference
- cell type column name, or a metadata TSV

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
- Report chord cap: `max_chords=180`
- Reference panel selection: prefer `fastccc.report.list_reference_panels(reference_root)` and `fastccc.report.generate_reference_report(...)` when the user asks for a healthy tissue reference report

## Environment prerequisites

Before running FastCCC commands, verify the execution environment rather than assuming it exists:

- Python can import the local `fastccc` package from the current checkout or installed environment.
- Core dependencies needed for the selected workflow are available, especially `scanpy`, `anndata`, `numpy`, `pandas`, `scipy`, `matplotlib`, and `loguru`.
- The LRI database path exists; use `db/CPDBv5.0.0` by default when present.
- Input `.h5ad` or metadata files are readable from the current machine.
- Standard CCC analysis uses normalized log1p data or applies normalization before running.
- Reference building and reference-based inference use raw count data.
- Healthy tissue reference workflows need a downloaded reference root containing tissue panel folders with `config.toml`.
- Custom reference workflows need an already built FastCCC reference panel directory unless the task is to build that panel first.
- Enrichment/ORA panels may need network access for Enrichr through `gseapy`; if unavailable, continue report generation and surface the skipped ORA reason from the report audit.

If the user only asks for a plan, command template, or result explanation, do not require a runnable FastCCC environment. If the user asks to run analysis, check the environment and report the first missing prerequisite that blocks execution.

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
- To run reference inference and report generation together, use `fastccc.report.generate_reference_report(...)`.
- To render a report from existing reference inference outputs, use `fastccc.report.generate_infer_report(...)`.
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
