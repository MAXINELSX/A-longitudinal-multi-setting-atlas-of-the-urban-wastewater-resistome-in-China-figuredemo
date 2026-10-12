# Fig1 metadata demonstration

This offline module runs the current source's Fig1 b/c workflow on 96 historical metadata records, redraws d from stored historical accumulation summaries, and runs e on a separate 72-sample abundance subset with synthetic setting labels. These are separate demonstrations, not a reproduction of a single study population or the complete paper figure.

Run from any working directory:

```sh
Rscript path/to/Fig1/run_demo.R --output path/to/new_results
Rscript path/to/Fig1/tests/test_demo.R
Rscript path/to/Fig1/tests/test_technical.R
```

Use R 4.3 or later with `readr`, `dplyr`, `ggplot2`, `cowplot`, `scales`, and `vegan`. The runner checks dependencies and never installs packages or uses the network. It preserves existing result files by requiring a new output destination. R 4.6.0 was used for the supplied run.

`inputs/metadata.csv` contains only the five fields the original metadata loader requires. The first three eligible source city strata were chosen deterministically. Within each city/setting, eight evenly spaced records were selected after ordering by date and source sample identifier. This gives eight records per city/setting, with no resampling or generated numerical values. Sample and city identifiers were replaced with generic labels; source filenames, paths and reverse mappings are not distributed. Sequencing depth values and dates are unchanged. The inputs originate from the user-provided historical figure archive. The reviewed source-code snapshot is Git blob `2e14e7d0280e23df9034de02796976dbbe3c9cb3`.

`source/plot_functions.R` retains the calculations of the original `theme_fig1`, `normalize_setting`, `load_metadata`, `plot_panels_bc`, `plot_panel_d`, and `save_pdf_safely` functions from the reviewed source snapshot. The city palette and city lookup were adapted for generic metadata labels. The runner supplies relative input paths, neutral figure titles, compact output layouts, validation and CSV exports. Explanations and comments are collected in the READMEs. These are metadata, display and packaging changes; the original counting, logarithm and GBp-to-TBp calculations remain unchanged.

The output includes separate PDF/PNG figures for metadata, precomputed curves and technical PCoA, plus numeric summaries, a coverage table, and R/package versions. Metadata summaries are checked against `expected/summary.csv` and `expected/site_totals.csv`; the references were calculated separately from the source-selected decimal values. Tests verify input validation, exact subset dimensions, group alignment and an independent Bray-Curtis formula check.

Panel d uses `inputs/precomputed_setting_curves.csv` and `inputs/precomputed_city_curves.csv`: all 132 stored rows from four setting summaries and 718 rows from three anonymous city summaries. Column names and display labels were adapted; numerical means, sample counts and available setting SD values were preserved, including small fluctuations in stored means. This executes the original plotting function only. It does not reconstruct an abundance matrix or recompute random accumulation, and the curves represent historical source groups rather than the 96 selected metadata rows. Source endpoints and row counts are checked against `expected/precomputed_curves.csv`.

Panel e lives in `e_technical_demo`, with its own copies of 72 abundance samples and 2,646 subtype rows from a separate source. It executes the original PCoA and permutation-test functions. City labels identify anonymous source strata; setting/site/date fields are synthetic. Assigned-setting statistics are only software checks, not evidence of environmental differences. It produces a separate results subdirectory and figure. This limitation is documented here and in the branch README; plots and statistical tables contain only scientific labels, results and compact status codes. Generic sample/city labels must not be used to join the 72-sample and 96-sample branches.

The historical archive contains no demonstrated sample-by-subtype matrix matching its metadata or saved PCoA coordinates. Its saved workspace was inspected as data without running saved functions. None of these subset counts or technical results should be used as study-wide findings.

For verification, compare numeric CSVs with the expected references at a relative tolerance of 1e-12. PNG/PDF byte equality is not required across platforms because fonts and graphics libraries can differ. The plotting jitter uses the source's fixed seed.

`coverage.csv` uses `ran` for the b/c metadata calculations, `precomputed_plot_only` for the stored d curves, and `separate_technical_branch` for the independent e analysis. `Setting_status = SYNTHETIC_DEMO_ONLY` identifies assigned settings in the technical scores. The runner emits no separate scope file or interpretation sentence columns; genuine validation errors, permutation results and R/package versions remain available.
