---
title: Home
layout: home
nav_order: 1
description: FastCCC infers cell-cell communication from single-cell RNA-seq without permutations and compares it against healthy tissue reference panels, reporting which ligand-receptor interactions are up or down versus normal.
---

<div style="margin: 0 auto; text-align: center;">
  <img src="{{ site.baseurl }}/images/logo.png" width="420" alt="FastCCC" />
</div>

# Cell–cell communication that tells you what changed

FastCCC infers ligand–receptor communication between cell types in single-cell RNA-seq data **without permutations**, so a million cells take minutes, and compares your dataset against **built-in healthy tissue reference panels** (19 tissues, 16 million cells) to report which communications are **Up** or **Down** versus normal. Every run ends in a self-contained interactive page with publication-ready figures.

[Quick start]({{ site.baseurl }}/quickstart.html){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 }
[Reference comparison]({{ site.baseurl }}/guides/reference.html){: .btn .btn-green .fs-5 .mb-4 .mb-md-0 .mr-2 }
[GitHub](https://github.com/Svvord/FastCCC){: .btn .btn-blue .fs-5 .mb-4 .mb-md-0 .mr-2 }
[Paper](https://doi.org/10.1038/s41467-025-66272-z){: .btn .fs-5 .mb-4 .mb-md-0 }

## Three ways to use it

**Command line**

```bash
pip install fastccc
# one dataset
fastccc run data.h5ad -c cell_type -o results/sample
# the same dataset versus healthy liver
fastccc ref data.h5ad -c cell_type -r liver --auto-map -o results/vs_liver
```

**Python**

```python
import fastccc
res = fastccc.run("data.h5ad", "cell_type", out_dir="results/sample")
q = fastccc.compare_to_reference("data.h5ad", "cell_type",
                                 reference="liver", celltype_mapping="auto")
q.trend_counts, q.interaction_summary("Down")
```

**A coding agent** (Claude Code, Codex, Cursor, ...): `npx skills add Svvord/FastCCC`, then ask in plain language. See [Agent skill]({{ site.baseurl }}/guides/agent.html).

## What a result looks like

Every run ends in one self-contained HTML file. Below: PBMC3k compared with the healthy blood panel, a standard PBMC3k run with bootstrap, and a placenta dataset split by sample. Click a tile to open it at full size.

<style>
.fc-hero { margin: 0.6rem 0 0.9rem; }
.fc-hero img, .fc-gallery img { display: block; width: 100%; border: 1px solid #e3e1db; border-radius: 8px; background: #f7f6f2; }
.fc-hero figcaption, .fc-gallery figcaption { font-size: 0.8rem; line-height: 1.4; color: #5f5d57; margin-top: 0.4rem; }
.fc-gallery { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 16px 16px; margin: 0 0 1.4rem; }
.fc-gallery figure, .fc-hero figure { margin: 0; display: flex; flex-direction: column; }
.fc-gallery figure a { display: block; aspect-ratio: 16 / 10; }
.fc-gallery figure a img { width: 100%; height: 100%; object-fit: cover; }
.fc-gallery .tall { grid-row: span 2; }
.fc-gallery .tall a { flex: 1; aspect-ratio: auto; min-height: 0; }
.fc-gallery .tall a img { object-fit: contain; }
.fc-gallery figcaption b { color: #27251f; }
@media (max-width: 640px) { .fc-gallery { grid-template-columns: 1fr; } .fc-gallery .tall { grid-row: auto; } }
</style>

<div class="fc-hero">
<figure>
  <a href="{{ site.baseurl }}/images/explorer/overview.png" target="_blank"><img src="{{ site.baseurl }}/images/explorer/overview.png" alt="Result explorer: filter sidebar, summary tiles, key findings, sender-receiver heatmap and communication chord"></a>
  <figcaption><b>The explorer.</b> Filters and presets on the left (Default, dataset-adaptive <i>Highlights</i>, Custom); summary tiles, key findings written from the current selection, the sender × receiver heatmap and the communication chord. Every chart exports SVG.</figcaption>
</figure>
</div>

<div class="fc-gallery">
<figure>
  <a href="{{ site.baseurl }}/images/explorer/dotplot.png" target="_blank"><img src="{{ site.baseurl }}/images/explorer/dotplot.png" alt="Ligand-receptor pairs across cell-type pairs, grouped by sender"></a>
  <figcaption><b>Ligand–receptor pairs across cell-type pairs</b>, grouped by sender and coloured Up / Down / unchanged versus the reference.</figcaption>
</figure>
<figure class="tall">
  <a href="{{ site.baseurl }}/images/explorer/evidence.png" target="_blank"><img src="{{ site.baseurl }}/images/explorer/evidence.png" alt="Evidence panel: reference null distributions, significance band and the query value of one interaction"></a>
  <figcaption><b>Evidence for one interaction.</b> The reference-derived nulls of the ligand side, the receptor side and their combined score, the significance band, and where the query falls.</figcaption>
</figure>
<figure>
  <a href="{{ site.baseurl }}/images/explorer/patterns.png" target="_blank"><img src="{{ site.baseurl }}/images/explorer/patterns.png" alt="Communication patterns: NMF of cell types by ligand-receptor pairs"></a>
  <figcaption><b>Communication patterns.</b> Non-negative matrix factorisation of cell type × ligand–receptor pair: which cell types share which signals.</figcaption>
</figure>
<figure class="tall">
  <a href="{{ site.baseurl }}/images/explorer/samples.png" target="_blank"><img src="{{ site.baseurl }}/images/explorer/samples.png" alt="Samples by interactions matrix of a dataset split by sample"></a>
  <figcaption><b>Replication across samples.</b> Which interactions recur in which samples of a split run, samples ordered by their communication profile.</figcaption>
</figure>
<figure>
  <a href="{{ site.baseurl }}/images/explorer/queryref.png" target="_blank"><img src="{{ site.baseurl }}/images/explorer/queryref.png" alt="Query versus reference scatter of every interaction, coloured by trend"></a>
  <figcaption><b>Query vs reference.</b> Every interaction's fold over the null in the query and in the reference, coloured by trend; outlined points sit above their significance band.</figcaption>
</figure>
<figure>
  <a href="{{ site.baseurl }}/images/explorer/comparability.png" target="_blank"><img src="{{ site.baseurl }}/images/explorer/comparability.png" alt="Query-reference comparability: housekeeping-gene calibration and detection-rate shift"></a>
  <figcaption><b>Is my data comparable to the reference?</b> Housekeeping-gene calibration and the detection-rate shift, with a plain-language reading.</figcaption>
</figure>
<figure>
  <a href="{{ site.baseurl }}/images/explorer/stability.png" target="_blank"><img src="{{ site.baseurl }}/images/explorer/stability.png" alt="Bootstrap support of the selected interactions"></a>
  <figcaption><b>Which hits are solid?</b> Bootstrap support over cell resamples shows which interactions survive resampling the cells.</figcaption>
</figure>
<figure>
  <a href="{{ site.baseurl }}/images/explorer/pathways.png" target="_blank"><img src="{{ site.baseurl }}/images/explorer/pathways.png" alt="Pathway families by cell type, as sender and as receiver"></a>
  <figcaption><b>Pathway families by cell type.</b> Enrichment of signalling families per cell type as sender and as receiver, rows and columns clustered.</figcaption>
</figure>
</div>

The [result explorer]({{ site.baseurl }}/guides/explorer.html) also has a sortable table with a shortlist and CSV export, signalling roles, sender → ligand–receptor → receiver flows, in-page pathway enrichment, per-cell violins and embedding overlays, effect sizes, a figure builder that composes panels into one lettered SVG, and a methods paragraph with the run's numbers ready to paste.

## What is different about FastCCC

1. **Analytic p-values.** The null distribution of every communication score is computed instead of permuted: means by FFT convolution, quantile summaries by the exact distribution of order statistics of a discrete variable (ties included, a case the textbook formula gets wrong). Exact p-values, run time independent of the number of permutations, and datasets of millions of cells (the 1M-cell AIDA and 2M-cell HLCA atlases in the paper) run to completion. The machinery is usable on its own: [Distribution toolkit]({{ site.baseurl }}/methods/distributions.html).
2. **Sixteen scoring modules.** Four expression summaries × two complex aggregations × two ligand–receptor combinations, combined with the Cauchy test. The modules also separate *broad* signals (a typical cell expresses the ligand) from *sub-population* signals (a fraction of the cells carries it).
3. **Reference-based comparison.** Healthy tissue panels for 19 tissues, or a control group of your own dataset, serve as the null: every interaction is labelled Up, Down, Both Sig or Both NS, with the reference's own distributions shown as evidence and a comparability check of your data against the panel.
4. **Robustness you can see.** Bootstrap support over cell resamples, replication across samples, and standardised effect sizes that rank interactions whose p-values are at the floor.

## Where to go next

- [Installation]({{ site.baseurl }}/install.html) and [Quick start]({{ site.baseurl }}/quickstart.html)
- Use cases: [Compare with healthy tissue]({{ site.baseurl }}/use-cases/compare-cell-cell-communication-with-healthy-tissue.html) · [A million cells or more]({{ site.baseurl }}/use-cases/cell-cell-communication-million-cells.html) · [FastCCC compared with other tools]({{ site.baseurl }}/use-cases/fastccc-vs-cellphonedb-cellchat-liana.html)
- Guides: [Inputs]({{ site.baseurl }}/guides/inputs.html) · [One dataset]({{ site.baseurl }}/guides/standard.html) · [Conditions and samples]({{ site.baseurl }}/guides/conditions.html) · [Reference comparison]({{ site.baseurl }}/guides/reference.html) · [Databases]({{ site.baseurl }}/guides/databases.html) · [Result explorer]({{ site.baseurl }}/guides/explorer.html) · [Agent skill]({{ site.baseurl }}/guides/agent.html)
- Methods: [How FastCCC works]({{ site.baseurl }}/methods/how-it-works.html) · [Distribution toolkit]({{ site.baseurl }}/methods/distributions.html)
- Reference: [Command line]({{ site.baseurl }}/reference/cli.html) · [Python API]({{ site.baseurl }}/reference/api.html) · [Output files]({{ site.baseurl }}/reference/outputs.html) · [Reference panels]({{ site.baseurl }}/reference/panels.html) · [Databases]({{ site.baseurl }}/reference/databases.html) · [FAQ]({{ site.baseurl }}/reference/faq.html)

## Citation

Hou, S., Ma, W. & Zhou, X. FastCCC: a permutation-free framework for scalable, robust, and reference-based cell-cell communication analysis in single cell transcriptomics studies. *Nature Communications* 16, 11428 (2025). [doi:10.1038/s41467-025-66272-z](https://doi.org/10.1038/s41467-025-66272-z)
