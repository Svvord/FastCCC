# Changelog

## [1.1.0] - 2026-11-14

### Added
- **One-call Python API and a command line.** `fastccc.run`, `fastccc.compare_to_reference`,
  `fastccc.build_reference_panel`, `fastccc.inspect` and `fastccc.load_result`, and the commands
  `fastccc inspect | run | ref | build-ref | map | panels | explore | db | skill`. Every output folder
  carries a `fastccc_run.json` manifest, and every command has `--json` output.
- **Interactive result explorer.** Each run writes `explorer.html`, a self-contained page with filters and a
  Highlights preset, sender × receiver heatmap, chord diagram, dot plot, signalling roles, flows, pathway
  families and enrichment, communication patterns, per-cell expression, per-interaction evidence, a figure
  builder, SVG and CSV export, and a methods paragraph with the run's numbers.
- **Reference comparison in one command.** Healthy tissue panels are downloaded on first use
  (`fastccc panels`); `fastccc map` and `--auto-map` map your cell-type labels to a panel; a control group of
  the same dataset can be the reference (`-r condition:control:disease`); `fastccc build-ref` turns any
  dataset into a reusable panel. The explorer adds query-versus-reference, drivers, disturbance and
  comparability panels.
- **22 ligand–receptor databases** in the package: CellPhoneDB v4 and v5, CellChatDB, NicheNet v1 and v2,
  the LIANA+ consensus and its 15 source resources, and a mouse consensus. `fastccc db list | show |
  validate | import` manages them and installs your own.
- **Input inspection.** The matrix type is detected (raw counts, normalised, log-normalised), raw counts
  are normalised automatically, and a cell-type column is suggested. Ensembl gene ids are supported with
  `--gene-id ensembl`.
- **Robustness measures.** Bootstrap support over cell resamples (`--bootstrap N`), replication across
  samples for `--split-by` runs, and standardised effect sizes.
- **Choice of scoring modules** with `--summary`, `--aggregation` and `--combination`; the CellPhoneDB
  score is `--summary Mean --aggregation Minimum --combination Arithmetic`.
- **Agent skill** `fastccc` for Claude Code, Codex and other coding agents (`npx skills add Svvord/FastCCC`
  or `fastccc skill install`).
- **Distribution toolkit** usable on its own: `fastccc.distrib.order_statistic_pmf` and
  `distribution_from_values`.
- Test suite, continuous integration and `CITATION.cff`.

### Changed
- The null distribution of the quantile summaries (median, upper quartile, percentiles) is evaluated in
  closed form: exact for any number of cells and about ten times faster.
- The Cauchy combination of the scoring modules follows the ACAT conventions for p-values of exactly 0
  or 1, and module p-values are computed as tail probabilities. Results therefore differ from 1.0.x;
  `--small-p legacy` (`small_p_method='legacy'`) applies the previous combination rule, and
  `fastccc==1.0.1` reproduces the published numbers exactly.
- The reference workflows accept raw counts, normalised or log-normalised input and give the same result
  for all three.
- Result folders can be moved or renamed; recorded paths are resolved relative to the current folder.
- The command line writes progress messages to stderr, so stdout carries only results.
- The original entry points (`fastccc.statistical_analysis_method`,
  `fastccc.Cauchy_combination_of_statistical_analysis_methods`, `build_reference_workflow`,
  `infer_query_workflow`) remain available and resolve the former `./db/...` paths to the packaged databases.

## [1.0.1] - 2026-05-22

- Patch release.

## [1.0.0] - 2026-05-22

- HTML report generation and agent skill.
