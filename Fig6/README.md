# Figure 6 demo

## Requirements

R 4.6.0; `randomForest`, `rfPermute`, `lme4`, `lmerTest`, `ggplot2`. `performance` is optional. No GPU.

## Run

```sh
Rscript --vanilla run_demo.R --output demo_output
Rscript --vanilla tests/verify_outputs.R --output demo_output
```

## Input

`input/samples.tsv` contains 65 pseudonymized hospital-wastewater records with total ARG abundance and 30-day predictor summaries. `Sample`, `City` and `Hospital` identifiers are aliases. The code derives synthetic lag0, lag60 and lag90 predictor columns from these 30-day values solely to demonstrate the manuscript's four-window plotting layout; those derived windows are not measurements or manuscript estimates. RF uses the supplied 30-day values and complete cases.

## Output

Panels a and e are fitted and rendered as PNG/PDF. The run also writes RF importance, 72 LMM coefficient records from 16 fitted class-by-window models, diagnostics, plot tables, metadata and `test_results.csv`. Panels b–d and f–h are not executed. The reduced RF uses 200 trees and 99 permutations. Numerical values and significance in this demo are not manuscript results.
