# Fig6 shared function demo — only a and e actually run

This package runs fresh RF and mixed-effects fits, then uses the current public figure functions to create representative a/e plots. It is a small, pseudonymized hospital-sewage subset demonstration. It does not reproduce the manuscript's all-setting network estimates.

The public a/e panels are the **network** Total ARG RF/LMM panels; hospital-specific models are c/d/g/h. The supplied source tables contain hospital-sewage observations and 30-day context summaries. Consequently, this demo exercises the a/e function workflow on the available hospital subset with the adaptations below. It does not substitute hospital model results for an exact network reproduction.

## Run

Use R with `randomForest`, `rfPermute`, `lme4`, `lmerTest` and `ggplot2`. `performance` is optional for R². Packages must already be installed; the demo does not install packages or use the network.

```sh
Rscript run_demo.R --output demo_output
Rscript tests/verify_outputs.R --output demo_output
```

The script resolves inputs relative to its own location, so it can be run from a different working directory. Output paths containing spaces are supported when quoted. If packages are in a separate library, configure `R_LIBS_USER` before starting R.

Running the module creates a/e PNG and PDF files, fitted coefficient/importance tables, diagnostics, package versions, elapsed time, and the verification report in the chosen output directory. These are newly fitted outputs, not copied manuscript result tables.

## Inputs and provenance

- `input/samples.tsv`: 65 rows from 367 available source records, covering four city groups and six hospitals; sample, city and hospital IDs use `Sample...`, `City...`, `Hospital...` aliases. No patient-level records are included.
- The sample includes at most four city groups, two hospitals per city and 12 evenly spaced chronological samples per hospital. Selection uses metadata, not outcomes or significance.
- Total abundance is recovered from the supplied logged aggregate outcome as `max(10^logged_total - 0.0980514610309349, 0)`, allowing only tiny floating-point negatives within tolerance `1e-10` before clipping. All 18 selected exposure values are preserved.
- The source categories are sample-level aggregate outcomes, sample-level 30-day context summaries, and a supplied aggregate-outcome pseudocount table. Original source filenames, local paths, city/hospital names, original IDs and reverse mappings are omitted from delivery.
- `source/panel_helpers.R` retains the verified public function calculations with a generic local filename. Its source Git blob is `88c4b268223771de2887c3207c2e8888b3606103`, from source tree SHA `3fb17277a1dbc107a7b8adf9e29fbad28eb09b1c` (a Git tree identifier, not a commit identifier). Code comments and explanation-only metadata have been moved to README documentation; the runtime adapter is in `run_demo.R`.

## Explicit demo adaptations

1. Only the supplied 1–30-day exposure summaries are available. They remain `*_lag1`, and LMM plots show `Lag30`. No lag0, lag60 or lag90 data are created. The daily-to-window aggregation helper cannot be tested without daily source records.
2. RF uses those same 30-day predictors instead of the publication's contemporaneous inputs. The constant Hospital setting is excluded. City and month remain predictors.
3. RF uses the 63 complete cases, excluding two weather-incomplete records. This is an explicit demo change from source RF kNN imputation; it avoids adding an imputation dependency. LMM retains the unmodified source modelwise complete-case handling on all 65 input rows.
4. RF uses 200 trees, 99 response permutations, seed 123 and one core. The source permutation helper is called directly to bypass only the publication-scale 1000-tree/1000-permutation guard. Nominal permutation P values and the source signed importance normalization are retained. The reduced permutation count has coarse P-value resolution.
5. The unmodified source LMM function fits each of four exposure classes using all its original member predictors, `factor(season) + (1 | city) + (1 | site_nm)`, REML, Wald intervals and Satterthwaite tests. It min–max scales the outcome/exposures before modelwise complete cases. BH correction is within outcome × exposure class across the available single window, so its family is smaller than the publication's multi-window family.
6. Figure titles are neutral (`Figure 6a` and `Figure 6e`, with the Total ARG outcome). The hospital-subset and 30-day limitations are documented only in this README. The precipitation display label is “Precipitation proxy”; the supplied numeric field is retained without asserting that it is a 30-day cumulative total.

The reference run returned four successful LMM fits and 18 coefficients. Two fits (meteorological factors and digital prescriptions) were singular; this is retained in `lmm_diagnostics.csv` and `run_metadata.csv`. Passing the software checks means that models returned nonempty finite results and the requested plots/tables were generated. It does not validate the scientific estimates or establish that every random effect is supported by this subset. Operating-system timezone text is omitted from shared session metadata to avoid including a location name.

## Panel coverage

| Panels | Public function branch | This package |
|---|---|---|
| a | Total ARG network RF | Actually fitted and rendered using the explicit hospital/30-day adapter |
| e | Total ARG network LMM | Actually fitted and rendered, lag1 only |
| b, f | Tier I network RF/LMM | Not tested; the required separate Tier I outcome is not included |
| c, d | Hospital RF | Not tested; separate multi-table model inputs are required |
| g, h | Hospital LMM | Not tested; separate fitting/result-selection schemas are required |

All public a–h helper functions are included as a shared reference. The small input is only wired to the representative a/e workflow; the other panels must not be described as verified. `panel_status.csv` records this distinction explicitly.

`panel_status.csv` contains only `Panel` and compact `Status` fields: `ACTUALLY_RUN_ADAPTED_DEMO` for a/e and `NOT_TESTED` for the other panels. The scope explanations in the table above are not repeated as output sentence columns. `rf_imputation = complete_case` in `run_metadata.csv` identifies the complete-case adapter. Model formulas, numerical diagnostics, singularity flags and captured warnings are retained. In the optional hospital-quarter helper, coincident injection and combined estimates are not independent evidence; this interpretation belongs here rather than in a repeated `route_note` output field.

Verification requires an actual RF fit, 20 finite RF predictor importances with permutation P values, four successful LMM fits, 18 finite coefficients, eight displayed LMM predictors, the correct 30-day window/formula/BH rule, valid class shares, and four nonempty plot files. Empty result tables fail.
