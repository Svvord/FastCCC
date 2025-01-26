---
title: Code snippets
layout: default
parent: Basic usage
nav_order: 3
---
<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.js">
</script>

# Code Snippets

## Introduction
This section provides multiple templates, starting from the folder and data storage structure, followed by data download, processing, and the various ways to use FastCCC for different conditions. It aims to guide users in understanding how to apply our methods and help them reproduce our results effectively.

## Code snippets for using FastCCC with a single $$CS$$ statistic

### Example1: using CPDBTD dataset

Data is available for download [here](https://github.com/ventolab/CellphoneDB/blob/v5.0.1/notebooks/data_tutorial.zip).

#### Directory Structure
```
┌─ ...
├─ (your files)
├─ FastCCC (download from github)
    ├─ db
        ├─ CPDBv4.1.0 (LRI database provided.)
        └─ ...
    ├─ ...
    ├─ code_snippets.py (run your code here.)
    └─ ...
├─ data
    ├─ examples (save your downloads here.)
        ├─ metadata.tsv
        ├─ normalised_log_counts.h5ad
        └─ ...
    └─ ...
├─ (your files)
└─ ...
```

#### templates
```python
# copy following templates to code_snippets.py
import fastccc

LRI_db_file_path = './db/CPDBv4.1.0/'
meta_file_path = '../data/examples/metadata.tsv'
counts_file_path = '../data/examples/normalised_log_counts.h5ad'
convert_type = 'hgnc_symbol'

interactions_strength, pvals, percents_analysis = fastccc.statistical_analysis_method(
    LRI_db_file_path, 
    meta_file_path, 
    counts_file_path,
    convert_type
)
```

<blockquote class="new-title"> <p>output</p>
<div class="highlight"><pre class="highlight"><code><span class="sr">2025-01-25 21:33:20</span> | INFO     | Task id is 41bc60.
<span class="sr">2025-01-25 21:33:20</span> | INFO     | Results will be saved to "./results/".
<span class="sr">2025-01-25 21:33:20</span> | INFO     | Directory already exists: ./results/
<span class="sr">2025-01-25 21:33:21</span> | <span class="sb">SUCCESS </span> | <span class="sb">Data preprocessing done.</span>
<span class="sr">2025-01-25 21:33:21</span> | INFO     | Running:
-> Mean for single-unit summary function.
-> Minimum for multi-unit complex aggregation.
-> Arithmetic for L-R combination to compute the CS.
-> Percentile is 0.1.
<span class="sr">2025-01-25 21:33:45</span> | <span class="sb">SUCCESS </span> | <span class="sb">FastCCC calculation done.</span>
</code></pre></div>
</blockquote>

### Example2: using $$90^{th}$$ percentile as single-unit summary

#### templates
```python
# copy following templates to code_snippets.py
import fastccc

LRI_db_file_path = './db/CPDBv4.1.0/'
meta_file_path = '../data/examples/metadata.tsv'
counts_file_path = '../data/examples/normalised_log_counts.h5ad'
convert_type = 'hgnc_symbol'

interactions_strength, pvals, percents_analysis = fastccc.statistical_analysis_method(
    LRI_db_file_path, 
    meta_file_path,
    counts_file_path,
    convert_type,
    single_unit_summary = 'Quantile_0.9',
)
```

<blockquote class="new-title"> <p>output</p>
<div class="highlight"><pre class="highlight"><code><span class="sr">2025-01-26 09:01:53</span> | INFO     | Task id is 504f09.
<span class="sr">2025-01-26 09:01:53</span> | INFO     | Results will be saved to "./results/".
<span class="sr">2025-01-26 09:01:53</span> | INFO     | Directory already exists: ./results/
<span class="sr">2025-01-26 09:01:54</span> | <span class="sb">SUCCESS </span> | <span class="sb">Data preprocessing done.</span>
<span class="sr">2025-01-26 09:01:54</span> | INFO     | Running:
-> Quantile_0.9 for single-unit summary function.
-> Minimum for multi-unit complex aggregation.
-> Arithmetic for L-R combination to compute the CS.
-> Percentile is 0.1.
<span class="sr">2025-01-26 09:02:14</span> | <span class="sb">SUCCESS </span> | <span class="sb">FastCCC calculation done.</span>
</code></pre></div>
</blockquote>


### Example3: using geometric average as $$h(\cdot)$$

#### templates
```python
# copy following templates to code_snippets.py
import fastccc

LRI_db_file_path = './db/CPDBv4.1.0/'
meta_file_path = '../data/examples/metadata.tsv'
counts_file_path = '../data/examples/normalised_log_counts.h5ad'
convert_type = 'hgnc_symbol'

interactions_strength, pvals, percents_analysis = fastccc.statistical_analysis_method(
    LRI_db_file_path, 
    meta_file_path,
    counts_file_path,
    convert_type,
    single_unit_summary = 'Quantile_0.9',
    LR_combination = 'Geometric',
)
```

<blockquote class="new-title"> <p>output</p>
<div class="highlight"><pre class="highlight"><code><span class="sr">2025-01-26 09:07:09</span> | INFO     | Task id is 1de6e1.
<span class="sr">2025-01-26 09:07:09</span> | INFO     | Results will be saved to "./results/".
<span class="sr">2025-01-26 09:07:09</span> | INFO     | Directory already exists: ./results/
<span class="sr">2025-01-26 09:07:10</span> | <span class="sb">SUCCESS </span> | <span class="sb">Data preprocessing done.</span>
<span class="sr">2025-01-26 09:07:10</span> | INFO     | Running:
-> Quantile_0.9 for single-unit summary function.
-> Minimum for multi-unit complex aggregation.
-> Geometric for L-R combination to compute the CS.
-> Percentile is 0.1.
<span class="sr">2025-01-26 09:07:42</span> | <span class="sb">SUCCESS </span> | <span class="sb">FastCCC calculation done.</span>
</code></pre></div>
</blockquote>


## Code snippets for using FastCCC with multiple $$CS$$ statistics

### Example1: using default FastCCC's multiple analysis

#### templates
```python
# copy following templates to code_snippets.py
import fastccc

LRI_db_file_path = './db/CPDBv4.1.0/'
meta_file_path = '../data/examples/metadata.tsv'
counts_file_path = '../data/examples/normalised_log_counts.h5ad'
convert_type = 'hgnc_symbol'

fastccc.Cauchy_combination_of_statistical_analysis_methods(
    LRI_db_file_path, 
    meta_file_path,
    counts_file_path,
    convert_type
)
```

<blockquote class="new-title"> <p>output</p>
<div class="highlight"><pre class="highlight"><code><span class="sr">2025-01-26 12:16:33</span> | INFO     | Task id is 70b408.
<span class="sr">2025-01-26 12:16:33</span> | INFO     | Results will be saved to "./results/".
<span class="sr">2025-01-26 12:16:33</span> | INFO     | Directory already exists: ./results/
<span class="sr">2025-01-26 12:16:34</span> | <span class="sb">SUCCESS </span> | <span class="sb">Data preprocessing done.</span>
<span class="sr">2025-01-26 12:16:42</span> | INFO     | Running:
-> Mean for single-unit summary function..
-> Minimum for multi-unit complex aggregation.
-> Arithmetic for L-R combination to compute the CS.
-> Percentile is 0.1.
<span class="sr">2025-01-26 12:17:01</span> | <span class="sb">SUCCESS </span> | <span class="sb">CS scoring module calculation done.</span>
<span class="sr">2025-01-26 12:17:01</span> | INFO     | Running:
-> Mean for single-unit summary function..
-> Minimum for multi-unit complex aggregation.
-> Geometric for L-R combination to compute the CS.
-> Percentile is 0.1.
<span class="sr">2025-01-26 12:17:25</span> | <span class="sb">SUCCESS </span> | <span class="sb">CS scoring module calculation done.</span>
<span class="sr">2025-01-26 12:17:26</span> | INFO     | Running:
-> Mean for single-unit summary function..
-> Average for multi-unit complex aggregation.
-> Arithmetic for L-R combination to compute the CS.
-> Percentile is 0.1.
<span class="sr">2025-01-26 12:17:46</span> | <span class="sb">SUCCESS </span> | <span class="sb">CS scoring module calculation done.</span>
<span class="sr">2025-01-26 12:17:46</span> | INFO     | Running:
-> Mean for single-unit summary function..
-> Average for multi-unit complex aggregation.
...
<span class="sr">2025-01-26 12:21:36</span> | <span class="sb">SUCCESS </span> | <span class="sb">CS scoring module calculation done.</span>
<span class="sr">2025-01-26 12:21:36</span> | <span class="sb">SUCCESS </span> | <span class="sb">All CCC branches calculation done.</span>
<span class="sr">2025-01-26 12:21:36</span> | INFO     | Task ID for combining is :70b408
<span class="sr">2025-01-26 12:21:36</span> | INFO     | There are 16 pval files.
<span class="sr">2025-01-26 12:22:13</span> | <span class="sb">SUCCESS </span> | <span class="sb">Cauthy combination done.</span>
<span class="sr">2025-01-26 12:22:13</span> | INFO     | Integrating 16 interactions_strength files.
<span class="sr">2025-01-26 12:22:18</span> | <span class="sb">SUCCESS </span> | <span class="sb">Average CS across all scoring methods calculation done.</span>
</code></pre></div>
</blockquote>

### Example2: using $$2\times1\times2=4$$ scoring modules

#### templates
```python
# copy following templates to code_snippets.py
import fastccc

LRI_db_file_path = './db/CPDBv4.1.0/'
meta_file_path = '../data/examples/metadata.tsv'
counts_file_path = '../data/examples/normalised_log_counts.h5ad'
convert_type = 'hgnc_symbol'

fastccc.Cauchy_combination_of_statistical_analysis_methods(
    LRI_db_file_path, 
    meta_file_path,
    counts_file_path,
    convert_type,
    single_unit_summary_list = ['Mean', 'Quantile_0.9'],
    complex_aggregation_list = ['Minimum'],
    LR_combination_list = ['Arithmetic', 'Geometric'],
    save_path = '../../results/temp/',
    use_DEG = True
)
```

<blockquote class="new-title"> <p>output</p>
<div class="highlight"><pre class="highlight"><code><span class="sr">2025-01-26 09:38:39</span> | INFO     | Task id is c68f2f.
<span class="sr">2025-01-26 09:38:39</span> | INFO     | Directory already exists: ../../results/temp/
<span class="sr">2025-01-26 09:38:40</span> | <span class="sb">SUCCESS </span> | <span class="sb">Data preprocessing done.</span>
<span class="sr">2025-01-26 09:38:48</span> | INFO     | Running:
-> Mean for single-unit summary function..
-> Minimum for multi-unit complex aggregation.
-> Arithmetic for L-R combination to compute the CS.
-> Percentile is 0.1.
<span class="sr">2025-01-26 09:39:08</span> | <span class="sb">SUCCESS </span> | <span class="sb">CS scoring module calculation done.</span>
<span class="sr">2025-01-26 09:39:08</span> | INFO     | Running:
-> Mean for single-unit summary function..
-> Minimum for multi-unit complex aggregation.
-> Geometric for L-R combination to compute the CS.
-> Percentile is 0.1.
<span class="sr">2025-01-26 09:39:31</span> | <span class="sb">SUCCESS </span> | <span class="sb">CS scoring module calculation done.</span>
<span class="sr">2025-01-26 09:39:36</span> | INFO     | Running:
-> Quantile_0.9 for single-unit summary function..
-> Minimum for multi-unit complex aggregation.
-> Arithmetic for L-R combination to compute the CS.
-> Percentile is 0.1.
<span class="sr">2025-01-26 09:39:52</span> | <span class="sb">SUCCESS </span> | <span class="sb">CS scoring module calculation done.</span>
<span class="sr">2025-01-26 09:39:52</span> | INFO     | Running:
-> Quantile_0.9 for single-unit summary function..
-> Minimum for multi-unit complex aggregation.
-> Geometric for L-R combination to compute the CS.
-> Percentile is 0.1.
<span class="sr">2025-01-26 09:40:38</span> | <span class="sb">SUCCESS </span> | <span class="sb">CS scoring module calculation done.</span>
<span class="sr">2025-01-26 09:40:38</span> | <span class="sb">SUCCESS </span> | <span class="sb">All scoring modules calculation done.</span>
<span class="sr">2025-01-26 09:40:38</span> | INFO     | Task ID for combining is :c68f2f
<span class="sr">2025-01-26 09:40:38</span> | INFO     | There are 4 pval files.
<span class="sr">2025-01-26 09:40:44</span> | <span class="sb">SUCCESS </span> | <span class="sb">Cauthy combination done.</span>
<span class="sr">2025-01-26 09:40:58</span> | <span class="sb">SUCCESS </span> | <span class="sb">DEGs selection done.</span>
<span class="sr">2025-01-26 09:40:58</span> | INFO     | Integrating 4 interactions_strength files.
<span class="sr">2025-01-26 09:41:01</span> | <span class="sb">SUCCESS </span> | <span class="sb">Average CS across all scoring methods calculation done.</span>
</code></pre></div>
</blockquote>