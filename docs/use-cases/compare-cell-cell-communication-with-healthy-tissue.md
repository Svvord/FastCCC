---
title: Compare cell–cell communication with healthy tissue
layout: default
parent: Use cases
nav_order: 1
description: How to find which cell-cell communications are gained, lost or unchanged in a disease or tumour single-cell RNA-seq sample compared with healthy tissue, a control group or your own reference - FastCCC's 19 healthy human tissue reference panels, and one command to build your own.
---

# How to compare cell–cell communication in a disease sample with healthy tissue

**Short answer.** Run `fastccc ref` against one of FastCCC's 19 healthy human tissue reference panels. You do not need control samples of your own. Every ligand–receptor interaction between every pair of cell types is labelled **Up** (significant in your sample, not in healthy tissue), **Down** (significant in healthy tissue, not in your sample) or unchanged.

```bash
pip install fastccc
fastccc ref tumour.h5ad -c cell_type -r liver --auto-map -o results/vs_liver
```

The input is an `.h5ad` file with a cell-type column. Raw counts, normalised and log-normalised matrices give the same result, because the reference workflow ranks genes within cells. The output folder holds `query_infer_results.tsv` (one row per sender, receiver and ligand–receptor pair with its `trend_vs_ref`) and `explorer.html`, an interactive page that opens in any browser.

The same in Python:

```python
import fastccc
q = fastccc.compare_to_reference("tumour.h5ad", "cell_type", reference="liver",
                                 celltype_mapping="auto", out_dir="results/vs_liver")
q.trend_counts                  # how many Up, Down, Both Sig, Both NS
q.interaction_summary("Up")     # what your sample gained
q.interaction_summary("Down")   # what your sample lost
```

## The point is the comparison

A cell–cell communication analysis of one dataset tells you which interactions are significant in that dataset. It cannot tell you which of them are new, and it says nothing about the interactions that are missing. A comparison against a reference answers all three questions:

| | Analysis of one dataset | FastCCC comparison against a reference |
|---|---|---|
| What is significant in my sample? | yes | yes |
| Which of these are new? | – | **Up**: significant in your sample, not in the reference |
| Which are unchanged? | – | **Both Sig**: significant in your sample and in the reference |
| What has my sample lost? | – | **Down**: significant in the reference, not in your sample |

The interactions a sample has lost never appear in an analysis of that sample alone; with a reference they are listed, ranked and shown with their evidence.

## Three kinds of reference

**The built-in healthy tissue panels.** They were built from healthy data in CELLxGENE Census: 19 tissues, about 16 million cells and more than 450 cell types ([Hou, Ma & Zhou, *Nat Commun* 2025](https://doi.org/10.1038/s41467-025-66272-z)). The tissues are adipose, blood, bone marrow, breast, cerebral cortex, colon, eye, heart, kidney, liver, lung, lymph node, midbrain, prostate gland, skin, small intestine, spleen, stomach and thymus. `fastccc panels` lists them and `fastccc panels --show liver` prints a panel's cell types and cell numbers ([all panels]({{ site.baseurl }}/reference/panels.html)). A panel stores compact per-cell-type summaries of the atlas, not the cells themselves, so the atlas is never downloaded.

**A control group of your own dataset.** No mapping is needed, the labels are the same:

```bash
fastccc ref cohort.h5ad -c cell_type -r condition:control:disease -o results/disease_vs_control
```

**Your own reference, built with one command.** Any annotated dataset becomes a reusable panel: a published atlas, your lab's healthy cohort, another tissue, another species.

```bash
fastccc build-ref atlas.h5ad -c cell_type -n my_reference -o panels
fastccc ref query.h5ad -c cell_type -r panels/my_reference -o results/vs_my_reference
```

`--subset status:healthy` builds the panel from the healthy cells of a mixed dataset, and counts, normalised or log-normalised input all give the same panel. For mouse data add `--db liana_mouseconsensus` to both commands. Once built, the panel serves every later sample, so a whole project is compared against the same baseline.

## What FastCCC gives you to read the result

- **Cell-type mapping.** Panels use Cell Ontology names. `fastccc map tumour.h5ad -c cell_type -r liver -o map.json` proposes how your labels map onto them, lists what is ambiguous, and writes a file you can review and pass with `-m map.json`; `--auto-map` applies the suggestion directly.
- **Evidence for every interaction.** The *Evidence* panel draws the reference-derived null distributions and where your sample falls; the *Drivers* panel says whether a change comes from the ligand side, the receptor side or both.
- **Ranked summaries.** `q.interaction_summary("Up")` and `q.interaction_summary("Down")` rank the gained and the lost interactions by the number of cell-type pairs they affect.
- **Comparability check.** The explorer's *Comparability* panel compares the sequencing depth of your sample with the panel's, from housekeeping genes and detection rates.

The [reference comparison guide]({{ site.baseurl }}/guides/reference.html) covers each in detail.

Prefer not to run commands yourself? With the [agent skill]({{ site.baseurl }}/guides/agent.html), ask your coding agent: "Compare tumour.h5ad with the healthy liver reference panel and tell me which communications are up, down or unchanged."
