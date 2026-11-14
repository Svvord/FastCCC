---
title: Legacy API
layout: default
nav_order: 7
has_children: true
permalink: legacy
description: Documentation of the original FastCCC functions (Cauchy_combination_of_statistical_analysis_methods, build_reference_workflow, infer_query_workflow, the static HTML report). Still supported; new work should use the one-call API.
---

# Legacy API

These pages document the original entrypoints, which remain available and produce the same tables as the one-call API (`fastccc.run`, `fastccc.compare_to_reference`, ...). The one-call API adds input inspection, cell-type mapping, the run manifest and the result explorer; prefer it for new work (see [Quick start]({{ site.baseurl }}/quickstart.html) and the [Python API]({{ site.baseurl }}/reference/api.html)).

Paths in these pages such as `./db/CPDBv5.0.0` still work: the database now ships inside the package and such names fall back to the packaged copy.

- [Basic usage]({{ site.baseurl }}/usage): the standard analysis with `Cauchy_combination_of_statistical_analysis_methods`.
- [Human CCC reference]({{ site.baseurl }}/cccref): building panels and inferring with `build_reference_workflow` / `infer_query_workflow`.
- [Toolkits]({{ site.baseurl }}/toolkits.html).
