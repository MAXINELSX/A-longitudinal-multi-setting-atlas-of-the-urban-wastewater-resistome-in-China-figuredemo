# Figure 2 workflow demo

This package runs panels a–h on 72 distinct sample columns and all 2,646 subtype rows (190,512 unchanged abundance values). It demonstrates code execution. Every assigned setting, site and month is synthetic, so the resulting P values, letters and curves are not scientific comparisons or temporal estimates.

All study samples originate from one source category. `Hospital`, `WWTP`, `Community` and `Wet market` are retained solely as the four labels expected by the supplied analysis functions. `City01`–`City03` are anonymous source strata. `Site01`–`Site24` and dates in January–March 2000 are generated labels. Each assigned group has 18 distinct samples. Each assigned site has three distinct observations. The input values are neither repeated nor imputed.

## Run

```sh
Rscript run_demo.R
Rscript run_demo.R --output output_check
```

The script resolves inputs relative to its own location. No absolute source paths, private mappings or original sample identifiers are needed. Install the packages listed in `required` near the top of `run_demo.R` in a normal R library. `ggnewscale` enables the original panel f circular figure. If it is unavailable, the original f tables/tree coordinates still run and the supplied alternative circle/UpSet renderer is used with this difference recorded. Optional `DEMO_R_LIB` selects an existing additional R library; the script does not install packages.

```r
install.packages('ggnewscale')
```

## Inputs and adaptation

`input/subtype_abundance.csv` contains the first 24 sample columns from each of the first three source strata, with every subtype retained. `input/sample_metadata.csv` declares the anonymous city strata and all synthetic assignments. Numeric zeros are retained; missing, nonnumeric, negative or nonfinite values stop execution.

Panel d also uses `input/external_benchmark.csv`: all 35 examples pasted by the user, with numeric values retained unchanged. Region, country and sample labels are anonymous. The supplied pipeline removes the two `This Study` rows and replaces them with the 18 assigned `WWTP` demo totals, leaving 33 external records from four region groups plus the 18 demo observations (51 analysis rows). Only observed groups are plotted; missing regions are not fabricated. The assigned `WWTP` subset supplies the study comparison solely for workflow demonstration. These are not actual wastewater treatment observations.

The Kruskal–Wallis/Dunn comparison is exploratory and does not model repeated observations at physical sites or external study effects. Exact external duplicate rows are dropped by default, retaining their first occurrence. Audit tables record the duplicate policy and counts; `EXPLORATORY_KW_DUNN` is the corresponding status code.

The inputs comprise selected abundance columns and author-provided benchmark examples, prepared using the selection rules above. `input_checksums.csv` records their SHA-256 hashes. Source files are not modified and private reverse mappings are excluded.

`upstream/analysis_functions.R` is an adapted copy of the supplied Figure 2 code. The model formulas, six pairwise Wald contrasts, BH adjustment, prevalence and abundance thresholds, marker definitions, and core calculations are retained. The demo changes input routing, identifiers and location labels, the city display list, CSV export, figure captions and raster resolution. One numerical adaptation retries the negative-binomial fit with `Nelder_Mead` if the original `bobyqa` optimizer fails with a non-positive-definite update; the data, formula, family and Laplace approximation are unchanged, and both attempts are recorded. The full-cohort loader and CLI are omitted; a portable validated input adapter replaces them. Model environments are not exported to RDS. A neutral color for an unprefixed subtype label supports the optional alternative f renderer without changing its classification or counts.

## Outputs and interpretation

`output/panel_coverage.csv` records execution status. `output/input_validation.csv` checks sample counts, completeness and the balanced full-rank model design. Model directories contain fixed effects, six contrasts, variance components, residual plots and convergence/singularity diagnostics. A singular or warned fit remains explicitly reported; execution success does not imply a scientifically meaningful fit to synthetic groups. `output/warnings.csv`, package versions and the run record document the actual execution.

Figures contain normal titles, axes and legends without explanatory captions. Outputs include b/c model displays, d anonymous benchmark comparisons, e core prevalence/scatter displays, f core-sharing displays, g marker detection and h profiles across three assigned months, plus panel a composition. The supplied panel h spline is descriptive and its interpolation is retained as code behavior, not inferred time dynamics.

Core membership requires within-setting prevalence of at least 70% and mean abundance of at least 1e-5 copies/cell, with zeros retained. Detection means abundance greater than zero; the denominator includes all samples in a group. The circular tree is a Ward.D2 descriptive clustering with four clusters, not a phylogeny. Its outer bars summarize type means of log10 subtype abundance, where each subtype is averaged only over settings in which it is core. The alternative circular renderer omits the tree; its tracks run Wet market, Community, WWTP, Hospital from inside out, with grey indicating non-core membership. UpSet intersections are exclusive combinations. These definitions are recorded here rather than in separate annotation files.

The delivered run used the original f renderer. The b fit required the documented optimizer retry and captured an initialization convergence warning; its final convergence-message field is empty. The c fit has a singular random-effect estimate. Both diagnostics are retained. They reinforce that this is a workflow demonstration, with no claim that the assigned groups support inference.

Optional independent reconciliation uses only Python's standard library:

```sh
python validate_demo.py
python validate_demo.py --output output_check
```
