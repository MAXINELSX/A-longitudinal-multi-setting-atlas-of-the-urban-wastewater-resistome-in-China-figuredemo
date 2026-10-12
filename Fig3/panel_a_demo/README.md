# Fig3a historical relative-abundance demo

This small runnable example reproduces the historical central bubble-matrix transformation: take the arithmetic mean of supplied sample relative-abundance fractions within city-by-setting groups. The 16 supplied samples span two anonymous cities and four setting categories, with two samples per group. Six species produce 48 displayed group means.

This is a partial technical display, not the current count-pooled Figure 3a calculation. The supplied fractions retain their original normalization; the selected six species are not renormalized to sum to one. No raw counts or denominators are reconstructed. The historical global count and ARG-composition sidebars, full-data ranking and inferential claims are outside this example. The source scripts were inspected as data and were not executed.

The extraction selected the first two city codes lexically, then the first two matching sample identifiers in each city and setting. Species are the first six names lexically. This deterministic choice does not use abundance values, significance or health relevance. Original numerical input text is unchanged. Sample/city labels are pseudonyms; ordinary setting categories and biological names are retained. Private reverse mappings and original filenames are excluded.

## Run and verify

Requires R and the `ggplot2` and `scales` packages. No package is installed automatically. From this directory:

```sh
Rscript --vanilla run_demo.R --output reference_output
Rscript --vanilla test_adapter.R
Rscript --vanilla test_outputs.R --output reference_output
```

Absolute script paths also work from any other current directory. Inputs resolve relative to the script; an explicit output directory resolves relative to the caller. `DEMO_R_LIB` optionally adds a package-library directory. Omit `--output` to use the script's own `reference_output` directory. On Windows, an inherited unsupported locale can be avoided by setting `LC_ALL=C` before launching R.

`group_means.tsv` is the numerical plot input; `panel_a_demo.pdf` and `.png` are the figures. `diagnostics.tsv`, `versions.tsv` and `run_summary.txt` record warning class/count codes, versions and runtime. Full R warnings remain visible in the console if they occur. No separate scope or warning-text files are generated. Input values come from the supplied relative-abundance matrix and grouping metadata; the historical plotting code informed the transformation.

The renderer uses ggplot2 instead of the historical heatmap dependency, retains circle area proportional to the mean and the historical color saturation at 0.05, and uses a compact layout. It is not a pixel-identical figure recreation. All species use the same label style; no biological prioritization is added. The tests independently recompute all means from the supplied subset and check saved plot dimensions, labels and warnings.
