# Independent Fig1 e technical example

This branch executes the original source's Bray–Curtis PCoA, PERMANOVA and PERMDISP functions on 72 source-derived abundance samples and all 2,646 subtype rows. It has its own input copies and works without any other module.

```sh
Rscript path/to/e_technical_demo/run_demo.R --output path/to/new_results
```

The main Fig1 runner also runs this branch into a separate `e_technical_demo` results directory. Required packages are `vegan`, `dplyr`, `ggplot2`, `scales`, and `readr`; the main module lists its additional plotting requirements. No packages are installed by the scripts.

Abundance values are unchanged. The first 24 sample columns in each of the first three worksheet strata were retained from the separate user-provided subtype-abundance workbook, along with all 2,646 subtype rows. Generic city labels identify source worksheet strata. Setting, site and date fields are synthetic demonstration assignments; all observations came from one source category. PCoA coordinates therefore demonstrate calculation on actual abundance values, while assigned-setting PERMANOVA/PERMDISP results are software checks and cannot support environmental comparisons.

These 72 samples are separate from the 96 samples used for the main metadata panels and from the precomputed accumulation curves. Reused generic labels do not identify the same observations across branches. The figures are intentionally saved separately.

The calculations of the original `align_metadata`, `compute_pcoa`, `plot_panel_e` and `theme_fig1` functions are preserved with an anonymous city palette. The reviewed source-code snapshot is Git blob `2e14e7d0280e23df9034de02796976dbbe3c9cb3`. Metadata columns are adapted to the original function interface. The demonstration uses 99 permutations and seed 123; source defaults allow other values. The figure uses the neutral title `Figure 1e`; scope and interpretation explanations appear only in READMEs. `Setting_status = SYNTHETIC_DEMO_ONLY` remains a compact machine-readable flag in the score table. No separate scope text file or interpretation sentence columns are emitted. No private identifiers, source filenames, local paths or reverse mapping are included.
