# Fig3e precomputed cell-percentage demo

This example displays 150 supplied summary cells: six species by five ARG classes across the overall table and four setting categories. It preserves the historical plotting measurement, `pct_mags`, directly. The accompanying `n_mag` and `n_arg_carrying_mags` values are copied from the source cells, and their agreement with the supplied percentages is checked. No missing integer counts are synthesized or recovered by rounding.

This is a partial technical display. Aggregate counts and percentages do not establish the underlying complete MAG membership or the union-denominator workflow used in the current analysis. This demo does not deduplicate or reconstruct MAGs, rerun the complete Figure 3e pipeline, estimate uncertainty, or support manuscript inference. Overall cells come from their supplied overall table; they are not obtained by averaging setting percentages.

The six species and five ARG classes are selected lexically from the supplied display catalogues. The catalogue membership was preselected upstream; no new ranking by biological importance, magnitude or significance is introduced. All supplied overall and four setting cells for this subset are retained. Numerical input strings are unchanged. The setting labels are broad categories rather than specific institutions; no sample, city or institution identifiers occur. Original filenames and paths are excluded.

## Run and verify

Requires R and the `ggplot2` and `scales` packages. No package is installed automatically. From this directory:

```sh
Rscript --vanilla run_demo.R --output reference_output
Rscript --vanilla test_adapter.R
Rscript --vanilla test_outputs.R --output reference_output
```

Absolute script paths also work from any other current directory. Inputs resolve relative to the script; an explicit output directory resolves relative to the caller. `DEMO_R_LIB` optionally adds a package-library directory. Omit `--output` to use the script's own `reference_output` directory. On Windows, an inherited unsupported locale can be avoided by setting `LC_ALL=C` before launching R.

`displayed_cells.tsv` contains the plot input, including supplied counts and percentages. `panel_e_demo.pdf` and `.png` are the figures. `diagnostics.tsv`, `versions.tsv` and `run_summary.txt` record warning class/count codes, versions and runtime. Full R warnings remain visible in the console if they occur. No separate scope or warning-text files are generated. The inputs come from supplied overall and setting-specific precomputed cell tables and species/class display catalogues; the historical plotting code informed the renderer.

The compact ggplot2 facets retain the historical direct percentage-to-size and percentage-to-color mappings on a 0–100 scale. Layout, lexical display order and label style are adapted; the result is not pixel-identical to the historical figure. The zero-percent point uses the historical minimum point size. Source code was inspected, not executed. Tests verify retained cell values, ratios, output dimensions and the absence of inference fields.
