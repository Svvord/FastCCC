---
name: fastccc
description: "Run FastCCC cell-cell communication (CCC) analysis on single-cell RNA-seq data: ligand-receptor interactions between cell types, permutation-free so it finishes in minutes, two-condition comparison, and reference-based comparison against built-in healthy tissue panels (19 tissues) that reports which communications are Up or Down versus normal. Use when the user mentions cell-cell communication, ligand-receptor signalling, CellPhoneDB / CellChat / LIANA style analysis, 'which interactions are increased or decreased compared with healthy tissue', or has an .h5ad and wants to know how cell types talk to each other."
---

# FastCCC

FastCCC infers cell-cell communication from single-cell RNA-seq: which ligand-receptor
(LR) pairs connect which sender and receiver cell types. Its null distribution is
analytic (no permutations), so 100k cells take minutes, and it can compare a dataset
against a reference, either one of 19 built-in healthy tissue panels or a control group,
labelling every LR interaction as `Up`, `Down`, `Both Sig` or `Both NS` versus that
reference. Read [references/workflows.md](references/workflows.md) for complete
examples before running anything.

## Install and check

```bash
pip install fastccc          # Python 3.11-3.13; the LR database ships inside the package
fastccc --version
fastccc panels               # built-in healthy panels; downloaded on first use (~3 MB each)
```

If `fastccc` is not on PATH, use `python -m fastccc.cli`. Nothing else needs cloning.

## Pick the workflow

| The user wants | Command | Python |
|---|---|---|
| Look at a dataset first | `fastccc inspect data.h5ad` | `fastccc.inspect('data.h5ad')` |
| CCC of one dataset | `fastccc run data.h5ad -c cell_type -o out/sample` | `fastccc.run('data.h5ad', 'cell_type', out_dir='out/sample')` |
| Each condition separately, plus a comparison tab when there are exactly two | `fastccc run data.h5ad -c cell_type --split-by condition -o out/by_condition` | `fastccc.run(..., split_by='condition')` |
| Which interactions are up or down versus healthy tissue | `fastccc map data.h5ad -c cell_type -r liver -o map.json` then `fastccc ref data.h5ad -c cell_type -r liver -m map.json -o out/vs_liver` | `fastccc.suggest_celltype_mapping(data, 'cell_type', reference='liver')`, `fastccc.compare_to_reference(data, 'cell_type', reference='liver', celltype_mapping=...)` |
| Disease vs control inside one dataset, control as the reference | `fastccc ref data.h5ad -c cell_type -r condition:control:disease -o out/disease_vs_control` | `fastccc.compare_to_reference(..., reference=('condition', 'control', 'disease'))` |
| Build a reusable reference panel | `fastccc build-ref atlas.h5ad -c cell_type -n my_panel -o panels` (`--subset condition:control` to use part of it) | `fastccc.build_reference_panel(...)` |
| Report for results that already exist | `fastccc report out/sample` | `fastccc.load_result('out/sample').report()` |

Built-in panels: adipose, blood, bone_marrow, breast, cerebral_cortex, colon, eye, heart,
kidney, liver, lung, lymph_node, midbrain, prostate_gland, skin, small_intestine, spleen,
stomach, thymus. `fastccc panels --show liver` lists a panel's cell types (and downloads
the panel, ~3 MB, into `FASTCCC_HOME`, default `~/.fastccc`, if it is not there yet).

Choose one primary workflow. If the request is ambiguous between "one dataset" and
"compare against something", ask one short question. Lead with the reference
comparison when the user talks about disease, treatment, perturbation or "normal".

## Before running

1. `fastccc inspect data.h5ad --json` tells you: the matrix type of `adata.X`
   (raw_counts / normalized / log1p / scaled, with confidence), whether raw counts are
   available, candidate cell-type columns with scores, and a suggested `celltype_key`.
2. Use the suggested column only when its confidence is `high` and it matches what the
   user said; otherwise ask which column holds the cell types. Never guess silently.
3. The matrix type is handled for you: `run` normalises raw counts
   (`normalize_total(1e4)` + `log1p`) and uses log1p data as is; `ref` and `build-ref`
   rank genes within each cell, so counts, normalised and log1p input give identical
   results. Only scaled data (negative values) is rejected. If the detection is reported
   with low confidence and the user knows what the matrix is, pass
   `--input-type raw|normalized|log1p` to `run`.
4. Gene names must be gene symbols matching the database (HGNC for the human databases,
   MGI for `liana_mouseconsensus`). `fastccc db list --data <file>` shows how many of each
   database's genes the data contains; a low fraction means a species or gene-id mismatch.

## Reference comparison: map the cell types first

Built-in panels use Cell Ontology names (`CD4-positive, alpha-beta T cell`); user
labels rarely match them. Without a mapping FastCCC only compares identically named
cell types and refuses to run when nothing matches.

1. `fastccc map data.h5ad -c cell_type -r liver -o map.json` prints, for every query
   cell type, the panel cell types it would be pooled from, with a score and
   confidence, plus `ambiguous` panel types (tied between query types), query types
   with no counterpart, and the panel types left unused (best similarity below
   `--min-score`, default 0.5; lower it to see weaker candidates).
2. Review it with your biological knowledge before using it: check `low` entries,
   decide the ambiguous ones (for example a panel's generic `T cell` can only go to one
   query type, or to none), and make sure nothing crosses lineages (CD4 to CD8,
   classical to non-classical monocytes). Edit `map.json` (keys are panel cell types,
   values are query cell types; several keys may share one value, which pools them
   weighted by cell numbers). Tell the user what you changed.
3. Run `fastccc ref ... -m map.json`. `--auto-map` applies the suggestion unreviewed;
   use it only when the user asks for a quick look and say so in the answer.
4. The run prints the coverage ("7/8 query cell types covered, 99% of cells"). Query
   cell types without a counterpart appear as `Not in reference` in the results; name
   them in the answer.

When the user's query labels are finer than the panel's (two query types for one panel
type), merge them with `query_relabel={'CD4 naive': 'CD4 T', 'CD4 memory': 'CD4 T'}`
in Python, or ask the user which one should receive the panel type.

## Reading results

Every output folder has `fastccc_run.json` (inputs, parameters, timing, output paths,
coverage, mapping used). `fastccc ... --json` prints the same summary as the last line.

Every run also writes `explorer.html`: a self-contained interactive page (filters,
sender × receiver heatmap, chord, dot plot, signalling roles, pathway families,
communication patterns, pathway enrichment, per-cell violins and embedding, an
evidence panel per interaction, shortlist with CSV/SVG export). Point the user to it
first; `fastccc explore <result_dir>` rebuilds it (`--no-cells` for a smaller file,
`--legacy-report` on `run`/`ref` adds the older static report).

Standard run (`out/sample/`):
- `<task>_significant_results.tsv`: one row per significant sender, receiver, LR pair (p < 0.05)
- `<task>_Cauchy_pvals.tsv`: p-value matrix, cell-type pairs x LR pairs (the authoritative values)
- `<task>_average_interactions_strength.tsv`: communication score matrix
- `report/`: HTML report (overview, LR heatmaps, chord diagrams, pathway enrichment)
- Python: `r = fastccc.load_result('out/sample')`; `r.summary()`, `r.significant`, `r.pvals`, `r.strength`

Reference comparison (`out/vs_liver/`):
- `query_infer_results.tsv`: one row per sender|receiver x LR pair that is expressed in
  the query (ligand and receptor in at least 10% of the sender / receiver cells), plus
  the reference-significant interactions the query does not express at all (these are
  `Down` with `above_expr_threshold` False). `trend_vs_ref` is `Up` (significant in
  the query, not in the reference), `Down` (significant in the reference, not in the
  query), `Both Sig`, `Both NS`, or empty when the cell types are not in the
  reference; `comm_score`, `is_significant`, `is_significant_ref` and the
  `*_expr_percent` / `*_expr_percent_ref` columns give the detail (full column list in
  references/workflows.md)
- `query_interactions_strength.tsv`, `query_percents_analysis.tsv`
- `celltype_mapping.json`: the mapping actually used; `report/`: HTML report
- Python: `q = fastccc.load_result('out/vs_liver')`; `q.trend_counts`,
  `q.interaction_summary('Down')`, `q.interactions('Down')`, `q.summary()`,
  `q.celltype_mapping`, `q.reference.celltypes`

What to tell the user: number of cell types and significant interactions; the dominant
sender and receiver cell types; for a reference comparison, the Up and Down counts and
the top Up and Down ligand-receptor pairs; which query cell types were not in the
reference; the exact output paths. Rank Up and Down by ligand-receptor pair with
`q.interaction_summary('Up')` / `q.interaction_summary('Down')` (how many
sender-receiver pairs show the change, the strongest pair, expression in query versus
reference). Do not rank `Down` rows by `comm_score`: that is the query's score, which
measures what is left rather than what was lost, and one LR pair repeats across many
cell-type pairs. "Down" is as informative as "Up": communications present in the
reference that the query does not show. Do not interpret `Both NS` rows.

## Rules

- Use the entrypoints above; do not re-implement scoring, null distributions, Cauchy
  combination or the mapping logic.
- Defaults that are fine unless the user says otherwise: database `CPDBv5.0.0`,
  `min_percentile 0.1`, p < 0.05, report DPI 150 (300 for publication figures).
- 22 databases ship with the package (`fastccc db list`): CellPhoneDB v4/v5, CellChatDB,
  NicheNet, the LIANA+ consensus and its 15 sources, `liana_mouseconsensus` for mouse data.
  Only CPDBv5.0.0 carries a pathway classification. A user's own pair table is installed
  with `fastccc db import pairs.csv --name <name> --species <sp>` and then used via `--db <name>`.
  Mouse data: use `--db liana_mouseconsensus` (or an imported mouse database); the built-in
  reference panels are human and cannot be used with mouse queries.
- Keep the user's object untouched: the API copies before normalising.
- Standard runs write 16 intermediate per-method tables (several GB for large
  datasets); point `-o` at a disk with space.
- The standard-run report calls Enrichr through the network for pathway enrichment;
  without network it still builds and notes the skipped panel. Reference reports use
  the database's own pathway annotation and need no network.
- If the user only wants a plan or an explanation of existing outputs, do not run anything.

## Common errors

- `None of the N query cell types matches a cell type of reference panel ...`: run
  `fastccc map` and pass `-m`.
- `celltype_mapping problems ... did you mean '...'`: fix the spelling in the JSON.
- `celltype_key '...' is not a column of adata.obs`: use a column from `fastccc inspect`.
- `adata.X contains negative values`: scaled data; pass the unscaled matrix
  (`adata.raw.to_adata()` or a counts layer).
- `Reference panel ... was built with X but database='Y'`: panel and query must use the
  same LR database.
- `Disk quota exceeded`: standard runs need several GB; use another `-o`.

## Response style

State the chosen workflow and the cell-type column in one sentence, state the mapping
decisions you made, give exact output paths, summarise the biology in a few lines, and
offer the next step (report, a different panel, the two-condition comparison).
