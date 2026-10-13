# Fig2 demo

Requirements: Windows, R 4.6.0, Python 3.11+ for optional validation; no special hardware.

R packages (tested versions): readxl 1.5.0, openxlsx 4.2.8.1, ggplot2 4.0.3, dplyr 1.2.1, tidyr 1.3.2, stringr 1.6.0, purrr 1.2.2, scales 1.4.0, tibble 3.3.1, patchwork 1.3.2, lme4 2.0.1, multcompView 0.1.11, rstatix 0.7.3, readr 2.2.0, forcats 1.0.1, png 0.1.9, ggnewscale 0.5.2.

## Install

```r
install.packages(c("readxl","openxlsx","ggplot2","dplyr","tidyr","stringr","purrr","scales","tibble","patchwork","lme4","multcompView","rstatix","readr","forcats","png","ggnewscale"))
```

Typical CRAN-binary installation: approximately 10–20 minutes, depending on network speed.

## Run

```sh
Rscript run_demo.R --output output
python validate_demo.py --output output
```

Expected R output: `PASS`. Demo run time on the tested Windows system: approximately 2 minutes.

Main output: `output/Fig2_demo_composite.png` and `output/Fig2_demo_composite.pdf`. Panel files and tables are in `output/panel_*` and `output/panels_*`. Run status is in `output/panel_coverage.csv`.

Inputs: `input/subtype_abundance.csv`, `input/sample_metadata.csv`, `input/external_benchmark.csv`. Sample, site, city and date identifiers in the demo metadata are synthetic.
