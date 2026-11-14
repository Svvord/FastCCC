---
title: Agent skill
layout: default
parent: Guides
nav_order: 7
description: Use FastCCC from Claude Code, Codex, Cursor and other coding agents through the shipped agent skill.
---

# Agent skill

FastCCC ships a skill: a set of instructions that teaches a coding agent which FastCCC command answers which request, how to check the input, how to review the cell-type mapping before a reference comparison, and how to read the results. With it, the agent runs the analysis and summarises the findings from a plain-language request.

## Install

```bash
npx skills add Svvord/FastCCC        # into the agents it detects (Claude Code, Codex, Cursor, OpenCode, Gemini CLI, Copilot, ...)
fastccc skill install                # after pip install fastccc: into ./.claude/skills and ./.agents/skills
fastccc skill install -g -a claude   # for the user instead of the project; -a all for every agent
```

`fastccc skill show` prints the skill; the files are also in the repository under `.claude/skills/fastccc`.

## Ask

In Claude Code the skill is `/fastccc`, in Codex `$fastccc`; other agents pick it up from the request:

```text
/fastccc Compare ./data/tumour.h5ad with the healthy liver reference panel. The cell type column is cell_type. Which cell-cell communications are up or down compared with normal liver?
```

```text
/fastccc Run cell-cell communication analysis on ./data/sample.h5ad (cell types in cell_type) and tell me which interactions are specific and robust.
```

```text
/fastccc In ./data/cohort.h5ad use condition == control as the reference and tell me which communications the disease cells gained or lost.
```

```text
/fastccc Run FastCCC per sample in ./data/cohort.h5ad (column sample) and tell me which interactions replicate across samples.
```

## What the agent does

1. Checks the environment (`fastccc --version`) and the file (`fastccc inspect`), and asks about the cell-type column only when the suggestion is not confident.
2. Picks one workflow: one dataset, per condition or sample, versus a healthy panel, versus a control group, build a panel, or explain existing results.
3. For a healthy panel, runs `fastccc map`, reviews the suggestion with its own biological knowledge (ambiguous and low-confidence entries, lineage crossings), edits the JSON and tells you what it changed, then runs `fastccc ref -m`.
4. Reads `fastccc_run.json` and the result objects, reports counts, the dominant senders and receivers, the top Up and Down ligand–receptor pairs ranked by how many cell-type pairs show them, which query cell types were not in the reference, the exact output paths, and points you to `explorer.html`.

Every command has `--json` for machine-readable output, and every result folder has `fastccc_run.json`, so agents without the skill can still work from the command line help.
