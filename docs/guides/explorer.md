---
title: Result explorer
layout: default
parent: Guides
nav_order: 6
description: Every panel of the FastCCC result explorer explained - filters and the Highlights preset, heatmap, chord, dot plot, roles, pathways, patterns, enrichment, evidence, stability, exports and the figure builder.
---

# Result explorer

`explorer.html` is written at the end of every `run` and `ref` (and by `fastccc explore <result_dir>` for any existing result). It is one file: data, D3 and the page code are inside, nothing is loaded from the internet, and it can be sent as an attachment. Every chart has an **SVG** button; the table exports CSV.

<p align="center">
  <img src="{{ site.baseurl }}/images/explorer/overview.png" width="860" alt="Explorer overview">
</p>

## Filters and presets

The left panel filters everything on the page at once: trend (reference comparisons), condition (split runs), gene text, sender and receiver cell types, thresholds, scoring modules, signal type, pathway classification. In reference comparisons the *Both NS* class (significant in neither dataset) is special: ticking it lists those pairs in the table and shows them as the background of the query-vs-reference scatter, but they never count as interactions in any other panel. *Not in reference* (significant in the query, cell types not covered by the panel) appears when such pairs exist and counts as an interaction without an Up/Down call. Three preset buttons show which state the filters are in; the active one is highlighted:

- **Default**: the run's own criteria, every significant interaction (p < 0.05, ligand and receptor expressed in at least 10% of the cells; for reference comparisons every Up, Down and Both Sig interaction).
- **Highlights**: Default plus dataset-adaptive thresholds for specificity, expression support, module agreement and, when present, bootstrap support. They are quantiles of the current dataset, start strict and relax step by step until at least 20 interactions remain; the values chosen are printed under the buttons. For reference comparisons it also restricts to Up and Down.
- **Custom**: your last manual settings. Changing any filter switches to Custom automatically; clicking Custom returns to those settings after a visit to Default or Highlights. The *reset* link forgets them and returns to Default.

The thresholds are sliders whose ranges come from the data:

| Slider | What it requires |
|---|---|
| Communication score | minimum score |
| Specificity | the ligand–receptor pair is significant in at most this many cell-type pairs |
| Expressing cells | ligand in the sender and receptor in the receiver expressed in at least this fraction of cells; starts at the run's 10% (reference comparisons: in the query or in the reference, starts at 0) |
| p-value | Cauchy p below this (standard runs) |
| Effect size | at least this many null standard deviations above expectation (standard runs) |
| Bootstrap support, Replication | when the run has them |

**Module support** (standard runs) is a stricter filter inside the set of interactions called significant by the 16-module Cauchy combination; it never re-tests. Choose which modules count, as 4 × 2 × 2 factors: expression summary (Mean, Median, Q75, Q90), complex aggregation (Min, Avg) and ligand–receptor combination (Arith, Geo); a module counts when all three of its levels are ticked. Below it, set how many of the selected modules must call the interaction significant on their own. By default all sixteen are ticked and the number is 0 (off), so the run's own call stands; changing the selection switches the number to 1, and setting it back to 0 switches the filter off again. To test a different module combination rather than filter within the published one, choose the modules when running FastCCC (`--summary`, `--aggregation`, `--combination`, or the `*_list` arguments of `fastccc.run`). The summaries have a biological reading: *Mean* is the average over all cells and can be carried by a few very high cells; *Median* is significant only when a typical cell expresses the gene, a broad signal; *Q75* and *Q90* (upper quartile, 90th percentile) are significant when the top quarter or top tenth of the cells express it, a sub-population signal. *Min* requires every subunit of a complex to pass where *Avg* averages them; *Arith* lets one strongly expressed side carry the pair where *Geo* needs both. Ticking only Q90 and requiring all four of its modules, for example, lists interactions carried by sub-populations; the *Signal type* filter below is the automatic version of the same reading.

## Panels

**Key findings.** A few sentences computed by fixed rules from the current selection: counts, strongest sender and receiver, most recurrent and most specific interactions, largest signalling family, Up/Down counts and driver shares, most disturbed cell type, stability and replication summaries, top enriched gene set. Numbers only, no interpretation.

**Sender × receiver.** Interactions per cell-type pair; for reference comparisons the colour is Up minus Down. Click a cell to focus the whole page on that pair; click again to release.

**Communication flow.** A directed chord diagram; ribbons run from sender to receiver, coloured by sender or (reference comparisons) by trend. Beyond 36 cell types the least connected ones are pooled into *other*.

**Ligand–receptor pairs across cell-type pairs.** The paper's dot plot: rows are the most recurrent ligand–receptor pairs of the selection, columns are cell-type pairs grouped by sender, dot size is the communication score (reference) or −log10 p (standard), colour is trend or condition. Choose how many pairs and columns to show.

**Signalling roles.** Outgoing versus incoming interactions (or summed score) per cell type; click a dot to see that cell type as sender.

**Flow: sender → ligand–receptor → receiver.** A Sankey of the top pairs of the selection.

**Pathway families by cell type.** Significant interactions per cell type (as sender, as receiver) grouped by family. *Group by* offers the CellPhoneDB classification (163 named families plus *Unclassified* for the 136 database pairs without one) or a gene-set library (MSigDB Hallmark, KEGG, GO Biological Process): with a library, a ligand–receptor pair counts in every set that contains its ligand or receptor, and the busiest 30 sets of the selection are shown. *Rows* either lists the top families overall, as a ranked list, or (the default for reference comparisons) the *cell-type specific* view: for every cell type (separately as sender and as receiver) the three families most over-represented among its interactions (hypergeometric test of the family's share in that cell type against its share in all selected interactions, at least three interactions) are collected, colours show the log2 enrichment of each family in each cell type relative to the average cell type (red above, blue below), and rows and columns are ordered by hierarchical clustering so shared families appear as bands and cell-type-specific ones as blocks. The two blocks have their own rows; with many cell types they are stacked vertically. *Rank by* orders families by number of interactions or by enrichment of the selection's ligand–receptor pairs against all pairs of the result (hypergeometric test, BH-adjusted q shown next to the name; only informative once the selection is narrower than the whole result). *Show top* sets the row limit. Clicking a row filters the whole page to that family or set; clicking it again clears the filter. A database without a classification column falls back to the gene-set libraries.

**Effect size vs significance** (standard runs). The standardised effect of the mean-based module against −log10 p, with the p-value floor marked: interactions at the floor are ranked only by effect size.

**Stability** (runs with `--bootstrap`). Histogram of bootstrap support.

**Samples × interactions** (split runs with three or more groups). Which interactions replicate across samples, with a PCA of the samples.

**Signal type from the 16 scoring modules** (standard runs). Per sender, how many interactions are *broad*, *sub-population*, *mean-driven* or *weak*, plus *complex-limited* and *one-sided* counts; the **module grid** shows the 16 p-values of the selected interaction.

**Disturbance per cell type** (reference comparisons). Two shares per cell type, counting an interaction for its sender and its receiver: *gained* = Up / (Up + Both Sig), the share of the cell type's query communications that are new versus the reference; *lost* = Down / (Down + Both Sig), the share of its reference communications that are missing in the query. The denominators follow the sender, receiver, gene and classification filters but not the trend ticks.

**Query–reference comparability** (reference comparisons). See [Reference comparison]({{ site.baseurl }}/guides/reference.html).

**Query vs reference.** Every interaction of the selection, coloured by trend. The default axes are *fold over expectation*: log2 of the score over the null expectation (the score the pair would have if sender and receiver were random cells of the reference), for the reference (x) and the query (y). The origin is "no cell-type specificity", the diagonal "same as the reference", so the four quadrants read as specific in both (Both Sig, along the diagonal), only in the query (Up), only in the reference (Down) or neither; outlined points have a query score above its significance band. *Combined score* shows the raw scores; there the diagonal is not the decision line, because the Up/Down call is made side by side: an interaction is Up when its ligand or receptor summary in the query exceeds the reference's interval (reference value × (1 ± 1.96/k), k the calibration factor) and its score is above its null, Down when a side falls below. Switch the axes to *ligand side* or *receptor side* to see that interval as a band through the origin: Up points sit above it on at least one side, Down points below. Both NS pairs appear here as the grey background when ticked.

**What drives the change.** For each cell type, stacked bars of its Up and Down interactions split by the side that moved: the ligand summary in the sender left its reference interval, the receptor summary in the receiver did, both did, or neither (the combined score alone left its band).

**Difference between conditions** (split runs). A − B interactions per cell-type pair; click a cell for the gained and lost pairs.

**Communication patterns.** Computed in the page on the current selection (non-negative matrix factorisation, recomputed a moment after the filters change; the number of patterns is automatic or chosen in the panel). Non-negative matrix factorisation of the cell type × ligand–receptor matrix, for outgoing and for incoming signalling, drawn as cell types → patterns → defining pairs.

**Pathway enrichment of the selection.** Hypergeometric test of the selected ligand and receptor genes against all genes of the database (Hallmark, KEGG and GO libraries from Enrichr, pruned to the database), Benjamini–Hochberg adjusted, recomputed in the page whenever the selection changes. It asks which signalling families are over-represented among the selected interactions.

**Cells.** Violins of the selected gene per cell type and the gene on the dataset's embedding (`obsm['X_umap']`, `X_tsne` or `X_pca`), from a stratified subsample of up to 4,000 cells (`--max-cells`, `--no-cells`).

**Interactions table.** Sortable; click a row for its evidence; tick rows to build a **shortlist** (export as CSV, download the dot plot of the shortlist, or show only the shortlist).

**Evidence.** For the selected interaction: p-value or trend, scores, specificity, the module grid (standard) or the reference null distributions with the significance band, the query statistic and the reference's own score (reference comparisons), and the fraction of cells expressing the ligand and the receptor in every cell type.

**Figure builder.** Tick panels and download one SVG with lettered panels, laid out in two columns, ready for Inkscape or Illustrator.

**Methods text.** A paragraph with the run's numbers (version, database, cells, parameters, panel, mapping, calibration, bootstrap, enrichment libraries) to paste into a manuscript; the *Run details* section below it lists the manifest, the panel and the mapping that was used.

## Rebuilding and options

```bash
fastccc explore results/sample                        # default: everything
fastccc explore results/sample --no-cells             # smaller file
fastccc explore results/vs_liver --no-both-ns         # drop Both NS rows (reference comparisons)
fastccc explore results/sample --no-genesets --patterns 4 --max-cells 2000 -o report.html
```

File size grows with the number of rows and cell types: about 2 MB for a small dataset, 25 MB for a 42-cell-type reference comparison with all rows and per-cell data.

## Sharing a view

The page keeps its state in the URL hash: `explorer.html#highlights`, `#row=123` (select and pin an interaction), `#pair=3,5` (focus a sender–receiver pair) open the page in that state.
