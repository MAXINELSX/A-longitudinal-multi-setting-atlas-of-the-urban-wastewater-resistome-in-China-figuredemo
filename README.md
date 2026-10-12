# Analysis demo

Companion code for *A longitudinal multi-setting atlas of the urban wastewater resistome in China*.

This package contains small example inputs, scripts, tests and reference outputs. It demonstrates selected downstream workflows, not complete reproduction of the manuscript. Figure 3c and Figure 4e use small simulated datasets; their values are not research estimates. Other examples use selected supplied values and, where documented, synthetic grouping metadata. The demonstrations must not be interpreted as new biological findings.

All methodological explanations, data provenance and interpretation notes are kept in README files. Code comments and explanatory figure captions are omitted. Normal axes, legends, statistical labels, numeric outputs, validation errors and model diagnostics are retained.

Source repository: [manuscript code](https://github.com/MAXINELSX/A-longitudinal-multi-setting-atlas-of-the-urban-wastewater-resistome-in-China).

## Coverage

| Module | Included branches |
| --- | --- |
| Fig1 | b/c metadata displays; d precomputed curve redraw; separate e technical PCoA example |
| Fig2 | a-h using a small abundance cohort and synthetic setting/site/month assignments |
| Fig3 | a/b/e processed-table displays; c simulated example; d burden functions with synthetic metadata |
| Fig4 | c cluster display; d selected pair-metric display; e simulated summary example |
| Fig5 | b-g, using a separate classification fixture and the Fig2-compatible abundance subset |
| Fig6 | a/e only, using a reduced hospital predictor table |

Successful execution means that these included branches ran, not that every panel or upstream analysis in the paper was reproduced. Module README files specify input schemas and the exact calculations.

## Requirements

Tested on Windows 11 x64 with R 4.6.0. Other operating systems have not been tested. No GPU or unusual hardware is required. R Cairo graphics support is required. Fonts may differ between computers.

Install R and these packages:

```r
install.packages(c(
  "cowplot", "data.table", "dplyr", "forcats", "ggalluvial", "ggnewscale",
  "ggplot2", "lme4", "lmerTest", "multcompView", "openxlsx", "patchwork",
  "purrr", "randomForest", "readr", "readxl", "rfPermute", "rstatix",
  "scales", "stringr", "tibble", "tidyr", "vegan", "VennDiagram"
))
install.packages("performance")
```

`performance` is optional for Fig6 diagnostics. Exact tested package versions are recorded in `reference_output/package_versions.csv`. These installation commands do not freeze dependency versions. A clean-install duration has not been measured; most packages were already present on the test computer. Python 3.12 is optional for the independent Fig2 validator; the main runner requires only R.

If using a separate R library, set `R_LIBS_USER` before starting R. `DEMO_R_LIB` can point to an additional installed library.

## Run

Unzip the complete package, open a terminal in its root and run:

```sh
Rscript --vanilla check_dependencies.R
Rscript --vanilla run_all.R
Rscript --vanilla verify_outputs.R --output generated_output
```

Or select another new output directory:

```sh
Rscript --vanilla run_all.R --output my_run
Rscript --vanilla verify_outputs.R --output my_run
```

Windows PowerShell, if R is not on PATH:

```powershell
& 'C:/Program Files/R/R-4.6.0/bin/Rscript.exe' --vanilla './run_all.R' --output './my_run'
```

The root output directory must be new or empty. Inputs are resolved relative to each script, so execution from another working directory is supported. Quote paths containing spaces. The scripts do not download data, install packages or upload files. A failed module stops the root run with a nonzero exit status.

## Outputs and tests

`reference_output/` contains a tested run: PDF/PNG figures, numerical tables, fit diagnostics, package versions, logs and module runtimes. `run_summary.csv` records each module's exit status and time. Allow several minutes on a desktop computer; actual times are recorded for the tested version. Raw-sequence processing is not part of this demo.

`verify_outputs.R` checks saved numerical results, input/output dimensions, included branches and nonempty figures. It also calls the Fig3-Fig6 module validators. Optional checks:

```sh
python Fig2/validate_demo.py --output generated_output/Fig2
Rscript --vanilla Fig1/tests/test_demo.R
Rscript --vanilla Fig1/tests/test_technical.R
Rscript --vanilla Fig3/test_fig3_audit.R
Rscript --vanilla Fig3/test_fig3_demo.R
Rscript --vanilla Fig4/test_demo.R
```

Compare numeric results rather than PDF/PNG byte hashes, which can change with graphics metadata and fonts. `PACKAGE_SHA256.csv` records the exact distributed files. Serialized R workspaces/models are excluded from the distributed reference.

## Data and calculations

- Fig2 and Fig5 share anonymous sample identifiers and compatible synthetic metadata. Their selected abundance values are retained; other modules' similar-looking aliases must not be assumed to identify the same samples.
- Fig3c uses simulated contig counts and fractions, with setting summaries calculated from the same simulated rows. Fig4e uses a small simulated count system with matching context/coverage totals and a separately identifiable classified subset. These examples do not use the manuscript's full totals.
- Fig3a/b/e start with processed tables. They do not reconstruct absent raw measurements, complete registries or upstream workflows. Fig3d retains the documented input row ambiguities and uses synthetic site/date metadata. See its README for handling.
- Fig2 retains its documented optimizer retry after a numerical failure; the model formula and supplied values are unchanged. Singular fits and warnings remain in diagnostics.
- Fig6 fits only the included a/e examples with reduced predictors and a small permutation budget. Because they all the same.

Private source filenames, location/sample identities, source paths and reverse mappings are excluded. This removes identifiers but does not guarantee that retained numerical subsets could never be linked to other records. Original user files and earlier packages are unchanged.

## Reuse and licence

To use other data, follow each module's schema and supply genuine metadata for real analysis. Some assertions intentionally target the bundled example dimensions and must be updated deliberately. Do not use the synthetic groups or simulated counts as manuscript data.

No new software or data licence has been assigned. Authors must confirm redistribution permission and select a licence before public release. Clean-install testing and review by an unfamiliar colleague remain separate tasks. This package has not been pushed to GitHub.
