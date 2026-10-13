# Figure 4 demo

Tested: Windows, R 4.6.0, ggplot2 4.0.3, patchwork 1.3.2, Cairo graphics. Arial is used for plots. Identifiers are anonymous; Figure 4d display metadata and Figure 4e counts are simulated.

```r
install.packages(c("ggplot2", "patchwork"))
```

```sh
Rscript --vanilla run_demo.R
Rscript --vanilla verify_saved_outputs.R
Rscript --vanilla test_demo.R
```

Output: `reference_output/` with Figure 4c, 4d, and 4e PNG/PDF files and TSV checks. A tested run took approximately 34 seconds; runtime varies by computer. To choose an output directory, run `Rscript --vanilla run_demo.R --output PATH`.
