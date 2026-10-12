# Figure 5 b-g software demonstration

This module exercises the existing plotting and downstream statistical functions. It does not generate new biological priorities or reproduce the manuscript's full-cohort estimates.

## Inputs

- `input/abundance.csv`: 7,344 unchanged values, comprising 102 supplied subtypes for the same 72 distinct source samples used in Fig2. These are the first 24 sample columns in each of three source worksheet strata, from the supplied filtered workbook.
- `input/sample_metadata.csv`: sample/city aliases plus explicitly synthetic settings, sites and dates. All abundance observations originate from one source setting. The four setting labels, 24 sites and January-March 2000 dates only exercise the expected data structure. They are not real cross-setting or longitudinal observations.
- `input/classification.csv`, `classification_subset.csv`, `core_membership.csv` and `classification_provenance.csv`: a separate fixture drawn from existing author-provided classification outputs. The first eight records in each existing level were retained, giving 32 records and eight existing Tier I records. Their labels, mobility summaries and core memberships were retained. This artificially balanced selection does not represent overall level frequencies. It is not inferred from the 72-sample abundance subset.

Fig2 and Fig5 share sample aliases. Matching aliases in other modules must not be treated as the same individuals or samples. Original filenames, sample identifiers, city names, local source paths and reverse mappings are not distributed. Source selection and adaptations are documented in this README.

## Run and verify

```sh
Rscript --vanilla run_demo.R --output demo_output
Rscript --vanilla test_demo.R demo_output
```

Required R packages: `ggplot2`, `dplyr`, `tidyr`, `tibble`, `patchwork`, `ggalluvial`, `VennDiagram`, `vegan`, `lme4`, `multcompView`. Install them before running. Optional `DEMO_R_LIB` selects an additional local package library; scripts do not install dependencies or use the network.

## What runs

| Panels | Executed workflow |
| --- | --- |
| b/c | Original alluvial, membership Venn and class-donut functions, using the existing 32-record/8-Tier-I classification fixture |
| d | Bray-Curtis PCoA, unrestricted sample-level PERMANOVA and PERMDISP; 99 permutations, seed 123 |
| e | Original Gaussian mixed model for total abundance, with six setting contrasts and BH adjustment |
| f/g | Original class composition and monthly summaries for all 102 supplied subtypes |

Synthetic settings/sites/months mean that d-g estimates, P values, letters and trajectories have no scientific interpretation. The classification fixture in b/c has real existing memberships but is a selected, artificially balanced example. Successful software execution does not validate either fixture as a representative study sample.

Outputs include panel figures, plotted data, sample summaries, model diagnostics, contrasts, permutation statistics, `panel_status.csv`, package versions and measured runtime. The checker independently reconciles sample sums, finite contrasts and P values, permutation counts and six principal PDFs. The public reference outputs omit serialized R objects; the original functions can generate those objects during a new local run.

Input routing and neutral figure titles are adapted from the source identified by Git blob `a31cdb9879281f1f0b34884fd9eec9b1f42f1e4f`. Explanations and comments are collected in this README; no `DEMO_SCOPE.txt` or blank-policy text file is emitted. Scientific labels, numeric annotations, significance letters, warning messages and model diagnostics remain. The statistical calculations are retained. Full study reproduction, clinical interpretation and the original classification construction are outside this demo.

For panel b, each subtype contributes once per setting in which it is core. Panel c uses exact core-membership regions, with H = Hospital, W = WWTP, C = Community and M = Wet market. Donut areas are not scaled by subtype count; their centre labels give region totals.

Panel f uses arithmetic means across individual samples within each city and assigned setting. Blank city-setting or month entries indicate no usable abundance data; crosses indicate sampled groups for which all class means are zero. The demo dates are January–March 2000. Panel g clips positive bubble-size abundance at the global 5th/95th percentiles and then uses log10 size scaling. Crosses indicate observed zeros and blanks indicate no usable abundance data. Colour is `(monthly mean - row mean) / row mean`, with a symmetric 5th/95th-percentile cap; row means include sampled zeros. These calculations do not give the assigned dates a real temporal interpretation.

In the optional workbook-reading helper, blank subtype abundance cells and numeric zeros denote nondetection; missing sample metadata are not imputed. Output status codes distinguish successful fixture workflows and singular mixed-model fits. The diagnostics tables retain singularity, convergence messages and any captured warnings.
