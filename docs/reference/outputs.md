---
title: Output files
layout: default
parent: Reference
nav_order: 3
description: Every file and column FastCCC writes - run manifest, significant results, p-value and score matrices, reference comparison table, bootstrap support, calibration diagnostics.
---

# Output files

## Every result folder

| File | Content |
|---|---|
| `fastccc_run.json` | the run manifest: `workflow` (`run`, `run_split`, `compare_to_reference`, `build_reference`), `name`, inputs (file, cell numbers, cell-type column and counts, detected matrix type and preprocessing applied), parameters, database (name, path, species, fraction of its genes present in the data), timing, output paths; reference comparisons add the panel, the mapping used, its source and the coverage; bootstrap runs add `bootstrap` |
| `explorer.html` | the interactive result page |

## Standard analysis (`run`)

`<task>` is the six-character task id in the manifest.

| File | Rows × columns | Content |
|---|---|---|
| `<task>_significant_results.tsv` | one row per significant interaction | `sender_celltype`, `receiver_celltype`, `LRI_ID` (CellPhoneDB interaction id), `ligand`, `receptor` (gene symbols; complex subunits comma-separated), `p-value` (Cauchy-combined) |
| `<task>_Cauchy_pvals.tsv` | `sender\|receiver` × `LRI_ID` | combined p-values for every tested pair |
| `<task>_Cauchy_stats.tsv` | same | the Cauchy statistics |
| `<task>_average_interactions_strength.tsv` | same | communication score averaged over the modules |
| `<task>_<timestamp>_<module>_pvals.tsv`, `..._interactions_strength.tsv`, `..._percents_analysis.tsv` | same | per-module p-values, scores and expression-filter matrices; `<module>` is `<summary>_<aggregation>_<combination>`, for example `Mean_Minimum_Arithmetic` |
| `bootstrap_support.tsv` (with `--bootstrap`) | one row per significant interaction | `n_significant`, `n_bootstrap`, `support` |
| `bootstrap/bNN/` | | the significant table of each resample |

Split runs put each group in its own sub-folder (`all/`, `<group>/`) with the files above; the manifest at the top lists them.

## Reference comparison (`ref`)

`query_infer_results.tsv` has one row per sender|receiver × ligand–receptor pair that passes the expression filter in the query, plus the reference-significant pairs the query does not express (these are always `Down`, with `above_expr_threshold = False`).

| Column | Meaning |
|---|---|
| `sender\|receiver` | sender and receiver cell types |
| `in_reference` | both cell types exist in the reference (after mapping) |
| `LRI_ID`, `ligand`, `receptor` | the interaction |
| `comm_score` | communication score in the query |
| `null_comm_score` | expected score under the reference-derived null |
| `sig_threshold_CI` | the significance band (`low-high`) derived from the reference null and the calibration factor |
| `above_expr_threshold`, `above_expr_threshold_ref` | passes the expression filter in the query / in the reference |
| `ligand_expr_percent`, `receptor_expr_percent`, `..._ref` | fraction of expressing cells in the query / the reference |
| `ligand_CS_component`, `ligand_CS_CI`, `receptor_CS_component`, `receptor_CS_CI` | the ligand-side and receptor-side summaries of the query, with their reference intervals |
| `is_significant`, `is_significant_ref` | significant in the query / in the reference |
| `trend_vs_ref` | `Up`, `Down`, `Both Sig`, `Both NS`; empty when `in_reference` is false |
| `ref_comm_score` | communication score of the same interaction in the reference |

Also written: `query_interactions_strength.tsv` and `query_percents_analysis.tsv` (score and filter matrices), `celltype_mapping.json` (mapping used), `query_calibration.json` and `query_hk_ranks.tsv` (calibration factor *k*, housekeeping-gene ranks in query and reference), `query_null_pmfs.npz` (compact null distributions, used by the explorer's evidence panel). A control-group comparison also writes the panel it built under `reference/`.

## Reference panel (`build-ref`)

`<out_dir>/<name>/`: `config.toml` (name, database, `min_percentile`, cell types and counts), `ref_mean_counts.pkl`, `ref_percents.pkl`, `complex_table.pkl`, `interactions.pkl`, `ref_hk.txt`, `fastccc_run.json`, and after first use `ref_gene_pmf_dict.pkl` (the activated null distributions; large, regenerated if deleted).
