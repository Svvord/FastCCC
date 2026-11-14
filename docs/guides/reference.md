---
title: Reference comparison
layout: default
parent: Guides
nav_order: 4
description: Compare a dataset against a healthy tissue reference panel or a control group - cell-type mapping, Up and Down calls, comparability checks, building custom panels.
---

# Reference comparison

The reference-based analysis asks a different question from the standard one: not "which interactions are significant in my data" but "which interactions differ from a reference". The reference is a panel: per-gene expression distributions of a large dataset, stored compactly. FastCCC ships panels for [19 healthy human tissues]({{ site.baseurl }}/reference/panels.html), and any dataset or control group can become one.

## Against a built-in healthy panel

```bash
fastccc panels                                        # the 19 panels; downloaded on first use
fastccc panels --show liver                           # the panel's cell types and cell numbers
fastccc map   data.h5ad -c cell_type -r liver -o map.json
fastccc ref   data.h5ad -c cell_type -r liver -m map.json -o results/vs_liver --name "Tumour"
```

### Cell-type mapping

Panels use Cell Ontology names (`CD4-positive, alpha-beta T cell`, `classical monocyte`), your labels rarely do. Without a mapping, only identically named cell types are compared, and `ref` refuses to run when nothing matches.

`fastccc map` proposes a mapping from name normalisation, a synonym table (NK, Treg, pDC, CD14+ monocytes, FCGR3A+/CD16+ monocytes, megakaryocytes/platelets, ...) and token overlap, never crossing lineages (CD4 vs CD8, classical vs non-classical). It prints, per query cell type, the panel types that would be pooled into it with a score and confidence, the **ambiguous** panel types (tied between query types, left unmapped), the query types with no counterpart, and the panel types left unused. Review it: delete what is biologically wrong (a pDC pooled into a myeloid DC cluster), decide the ambiguous ones, then pass the JSON with `-m`. Keys are panel cell types, values are query cell types; several keys may share a value, in which case they are pooled weighted by cell numbers.

`--auto-map` applies the suggestion without review (also `celltype_mapping="auto"` in Python); the mapping actually used is always written to `celltype_mapping.json` and recorded in the manifest. Before anything runs, the coverage is checked: an error if no query cell type is covered, a warning listing uncovered types otherwise (they appear as *Not in reference* in the results).

If your labels are finer than the panel's, merge them first: `compare_to_reference(..., query_relabel={"CD4 naive": "CD4 T", "CD4 memory": "CD4 T"})`.

### Reading Up and Down

Every interaction (sender, receiver, ligand–receptor pair) gets a `trend_vs_ref`:

| Trend | Meaning |
|---|---|
| **Up** | significant in the query, not in the reference |
| **Down** | significant in the reference, not in the query |
| Both Sig | significant in both |
| Both NS | significant in neither |
| (empty) | the cell types are not in the reference |

Down is as informative as Up: communications present in healthy tissue that the query has lost. The explorer and the result object give you what you need to read both:

- **Comparability.** The *Comparability* panel reports the calibration factor *k* from housekeeping genes and the median ratio of expressing-cell fractions query/reference, so the depth of your data is read next to the panel's.
- **Ranking.** `q.interaction_summary("Down")` ranks the lost interactions per ligand–receptor pair (number of sender–receiver pairs, strongest pair, expression percentages query versus reference), and the explorer's dot plot ranks by the number of cell-type pairs.

The *Drivers* panel says whether a change comes from the ligand side (the sender's ligand summary left its reference interval), the receptor side, or both; the *Disturbance* panel gives, per cell type, the share of its query communications that are new (Up / (Up + Both Sig)) and the share of its reference communications that are missing (Down / (Down + Both Sig)); the *Query vs reference score* scatter shows all interactions at once; and the *Evidence* panel draws, for any interaction, the reference-derived null distributions of the ligand side, the receptor side and their combination, the significance band, the query statistic and the reference's own score.

## Against a control group of the same dataset

```bash
fastccc ref cohort.h5ad -c cell_type -r condition:control:disease -o results/disease_vs_control
fastccc ref cohort.h5ad -c cell_type -r condition:control -o results/rest_vs_control     # query = all non-control cells
```

The control cells become a panel (written to `results/.../reference/`), the disease cells are the query, labels match by name. This is the statistically grounded version of a two-condition comparison.

## Building a reusable panel

```bash
fastccc build-ref atlas.h5ad -c cell_type -n healthy_pancreas -o panels
fastccc build-ref atlas.h5ad -c cell_type -n healthy_pancreas -o panels --subset status:healthy
fastccc ref query.h5ad -c cell_type -r panels/healthy_pancreas -m map.json -o results/query_vs_pancreas
```

A panel folder holds `config.toml` (name, database, cell types and counts), per-cell-type mean ranks and expression fractions, housekeeping-gene ranks and, after first use, the null distributions. Counts, normalised or log-normalised input all give the same panel (the workflow is rank-based within cells); whole-transcriptome input is better than a highly-variable-gene subset. A panel and a query must use the same database.

## Outputs

| File | Content |
|---|---|
| `query_infer_results.tsv` | one row per sender|receiver × ligand–receptor pair expressed in the query, plus reference-significant pairs the query does not express (always Down); columns in the [output reference]({{ site.baseurl }}/reference/outputs.html) |
| `query_interactions_strength.tsv`, `query_percents_analysis.tsv` | score and expression-filter matrices |
| `celltype_mapping.json`, `query_calibration.json`, `query_hk_ranks.tsv`, `query_null_pmfs.npz` | mapping used, calibration diagnostics, null distributions for the explorer |
| `fastccc_run.json`, `explorer.html` | manifest and the interactive page |

## Python

```python
sug = fastccc.suggest_celltype_mapping("data.h5ad", "cell_type", reference="liver")
print(sug.format()); sug.ambiguous; sug.unmatched_query
mapping = dict(sug.mapping); mapping["T cell"] = "T cells"
q = fastccc.compare_to_reference("data.h5ad", "cell_type", reference="liver", celltype_mapping=mapping, out_dir="results/vs_liver")
q.trend_counts; q.interaction_summary("Up"); q.interactions("Down"); q.celltype_coverage; q.celltype_mapping; q.reference.celltypes
panel = fastccc.build_reference_panel("atlas.h5ad", "cell_type", name="healthy_pancreas", out_dir="panels", subset=("status", "healthy"))
```
