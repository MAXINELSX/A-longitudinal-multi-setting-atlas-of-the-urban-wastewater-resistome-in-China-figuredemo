# Figure 3 demo

## Requirements

R 4.3+ with `data.table`, `dplyr`, `ggplot2`, `ggalluvial`, `multcompView`, `readr` and `scales`. Scripts do not install packages or access the network.

## Inputs

| Path | Role |
| --- | --- |
| `input/` | Six-table legacy bundle for panel d and input audit |
| `panel_a_demo/input/` | Relative abundance and anonymous metadata for a |
| `panel_b_demo/input/` | Precomputed fractions for b |
| `panel_c_demo/input/` | Simulated ARG-carrying-contig counts and matching summaries for c |
| `panel_e_demo/input/` | Selected species × setting × ARG-type cells for e |

Panel c uses simulated contig counts; panel d uses synthetic site/date metadata. Other panel inputs are selected processed examples with anonymous identifiers. Panel e's side bars summarize selected class carriage, not the manuscript's full ARG-burden estimates. No complete Figure 3 or panel f reproduction is claimed.

## Run

```sh
Rscript path/to/Fig3/run_demo.R --output path/to/new_results
```

Use a new output directory. `--input DIR` optionally replaces only the legacy six-table bundle; panel a/b/c/e inputs remain script-relative.

## Outputs and checks

The output directory contains `panel_a_demo/`, `panel_b_demo/`, `panel_c_demo/`, `fig3d_SYNTHETIC_metadata/`, `panel_e_demo/` and `input_audit/`, including plots, numerical tables and status records.

```sh
Rscript path/to/Fig3/verify_outputs.R path/to/new_results
Rscript path/to/Fig3/test_fig3_audit.R
Rscript path/to/Fig3/test_fig3_demo.R
Rscript path/to/Fig3/test_fig3_clean.R
```
