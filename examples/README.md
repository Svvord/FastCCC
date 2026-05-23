# FastCCC Examples

End-to-end scripts for running FastCCC and generating automated HTML reports.

| Script | Description |
|--------|-------------|
| [`report_single_condition.py`](report_single_condition.py) | Run FastCCC on one dataset and generate a report |
| [`report_multi_condition.py`](report_multi_condition.py) | Run FastCCC on a full dataset + two conditions in parallel, generate a report with a condition-comparison tab |

## Quick start

```bash
# single sample
python report_single_condition.py

# two conditions with condition comparison
python report_multi_condition.py
```

Edit the configuration block at the top of each script (file paths, condition column name, condition labels) before running.

## What the report looks like

**Executive Summary and run parameters**
![Overview](report_overview.png)

**CS distribution violin plot and supplementary interaction table**
![CS Violin and Table](report_cs_violin_and_table.png)

**L-R comparison evidence plot (two-condition comparison)**
![L-R Comparison Evidence Plot](report_differential_volcano.png)

## What the report contains

Each script produces a self-contained `report.html` that opens in any browser.  
The multi-condition report has four interactive tabs:

- **All Cells** — full dataset figures plus cell-type evidence explorers
- **Condition A** — condition-specific figures plus cell-type evidence explorers
- **Condition B** — condition-specific figures plus cell-type evidence explorers
- **Condition A vs Condition B** — interaction-count comparison, L-R evidence plot, pathway count comparison

For the full parameter reference and more details, see the
[Report Tutorial](https://svvord.github.io/FastCCC/usage/report.html).
