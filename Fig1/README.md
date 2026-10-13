# Figure 1 demo

## Requirements

R 4.3+ with `readr`, `dplyr`, `ggplot2`, `cowplot`, `scales` and `vegan`. The tested environment used R 4.6.0. Scripts do not install packages or access the network.

## Inputs

| Path | Required fields / role |
| --- | --- |
| `inputs/metadata.csv` | `Sample`, `Clean_Reads_GBp`, `City`, `Sample_Type`, `Sample_Date`; 96 selected, anonymized records for b/c |
| `inputs/precomputed_setting_curves.csv` | Stored setting-level accumulation curves for d |
| `inputs/precomputed_city_curves.csv` | Stored anonymous city-level accumulation curves for d |
| `e_technical_demo/inputs/abundance.csv` and `metadata.csv` | Separate 72-sample technical PCoA branch for e |

The b/c/d inputs are selected processed records, whereas e uses a separate abundance subset with synthetic setting labels. These demonstrations do not represent a single cohort or reproduce the complete manuscript figure.

## Run

```sh
Rscript path/to/Fig1/run_demo.R --output path/to/new_results
```

The output directory must not already exist. Input paths are resolved relative to the script.

## Outputs and checks

The output directory contains b/c, d and e PDF/PNG files, numerical summaries, `coverage.csv`, `session_info.txt`, and the e technical branch outputs. Panel d plots stored curves; it does not recompute accumulation from reads. Setting-based e statistics are software checks only.

```sh
Rscript path/to/Fig1/tests/test_demo.R
Rscript path/to/Fig1/tests/test_technical.R
```
