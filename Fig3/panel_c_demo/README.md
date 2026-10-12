# Panel c: simulated contig demonstration

Both figures use the same deliberately simulated table of 12 samples, with three samples per setting. These are illustrative numbers, not observed research counts or a subset of the full study. No full-study sample sizes are used.

`input/contig_observations.tsv` contains generic sample IDs and integer ARG-carrying contig counts. Each sample has 40–100 total contigs. Mobile equals plasmid plus virus, and total equals mobile plus chromosome. Each fraction is the corresponding count divided by total ARG-carrying contigs. The table has 840 total and 328 mobile contigs across the 12 illustrative samples.

`input/setting_summary.tsv` contains four summaries of those same sample-level mobile contig fractions. `N_samples` is three for every setting. Mean and median are 0.4 in every setting; the example only varies dispersion. The runner recomputes sample count, mean, median, R type-7 quartiles, minimum, maximum and sample standard deviation, then rejects a mismatching bundled summary. These are unweighted summaries of sample fractions, not pooled contig ratios. The summary plot shows mean, median and Q1–Q3; quartiles are descriptive intervals, not confidence intervals.

The violin plot displays all 12 points. With three points per setting, density shapes are illustrative. The original annotation-level mixed model is not run: this demonstration uses contig counts and contains no site or date metadata. No significance letters or inferential tests are produced. Setting colours are retained from the existing figure.

## Run

Requires R and `ggplot2`. The script resolves its default input relative to its own location and works from any current directory when given the script path.

```sh
Rscript --vanilla run_demo.R --output results
Rscript --vanilla test_outputs.R --output results
Rscript --vanilla tests/test_adapter.R
```

With no arguments, output goes to `reference_output/`. An alternate `--input DIR` can contain the same two schemas. If packages are in an extra library, set `DEMO_R_LIB` to that library directory before running R.

## Outputs and verification

`setting_summary.png/.pdf` and `contig_violin.png/.pdf` display the same sample table at summary and sample levels. The exported plot-input tables and `contig_descriptive_statistics.tsv` can be reconciled directly. Status codes and package versions are exported. No separate explanation files are written.

Tests independently recompute all summary statistics, verify the 12-sample total, count identities, 40–100 contig range, saved numerical values, plot units and the absence of model results or explanatory output products. An alternate input directory must contain matching sample and summary tables with the same schemas.
