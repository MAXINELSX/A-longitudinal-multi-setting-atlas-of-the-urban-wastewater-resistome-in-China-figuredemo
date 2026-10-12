# Panel b: processed-fraction demo

This compact, runnable example retains the first four ARG classes in source order and all their observed category combinations: **41 of 296 source rows**. Selection did not use fraction size, significance or biological interpretation. The selected fractions sum to **0.273137 (27.3137%) of the original global denominator**.

The input is a processed fraction table, not integer counts. Its exact source column names and numerical values are preserved. The adapted alluvial renderer uses those fractional weights directly. It does not renormalize the subset, reconstruct counts or add missing category combinations. Individual ARG class labels identify the four retained classes; the separate source display-label column is preserved in the exported input.

The source is the processed global-fraction table from [the supplied public repository](https://github.com/MAXINELSX/A-longitudinal-multi-setting-atlas-of-the-urban-wastewater-resistome-in-China). This is a partial plotting demonstration. The original count-based panel analysis was not run because compatible count input is unavailable.

## Run

Requires R and the `ggplot2` and `ggalluvial` packages. Paths below are relative examples; the script resolves its default input relative to its own location, so it works from any current directory when given the script path.

```sh
Rscript --vanilla run_demo.R --output results
Rscript --vanilla test_outputs.R --output results
Rscript --vanilla tests/test_adapter.R
```

With no arguments, output goes to `reference_output/`. An alternate `--input DIR` is supported for a table with the same schema and the same four demo classes. If packages are in an extra library, set `DEMO_R_LIB` to that library directory before running R.

## Outputs and verification

`fraction_flow.png` and `fraction_flow.pdf` show the fraction-weighted paths. `fraction_plot_input.tsv` retains all 41 input rows; `fraction_totals.tsv` reports their unscaled sum. `run_status.csv` distinguishes the successful adapted plot from the original analysis that was not run. Package versions are exported; separate explanation files are not generated.

Tests check row preservation, the global fraction sum, invalid and duplicate inputs, the plotted flow mass, and saved numerical output. Source-to-subset reconciliation was also performed privately without distributing original identifiers or local source paths.
