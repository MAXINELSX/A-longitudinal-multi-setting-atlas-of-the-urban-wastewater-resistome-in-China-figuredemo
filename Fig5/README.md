# Figure 5 demo

Requirements: R 4.6.0; `ggplot2`, `dplyr`, `tidyr`, `tibble`, `patchwork`, `ggalluvial`, `VennDiagram`, `cowplot`, `vegan`, `lme4`, `multcompView`; Cairo graphics.

Inputs:

- `input/abundance.csv`: `Sample`, `Subtype`, `Abundance`.
- `input/sample_metadata.csv`: `Sample`, `City`, `Setting`, `PhysicalSite`, `Sample_Date`, `SamplingMonth`, `Metadata_status`, `City_status`, `Source_setting`.
- `input/classification.csv` and `input/classification_subset.csv`: subtype, ARG type, tier, mobility, host breadth and setting membership.
- `input/core_membership.csv` and `input/classification_provenance.csv`: subtype membership and contig-support keys.

Partial demo only: the abundance and classification tables are selected subsets; setting, site and date metadata are synthetic.

```sh
Rscript --vanilla run_demo.R --output demo_output
Rscript --vanilla test_demo.R demo_output
```

Outputs: `demo_output/bc`, `d`, `e` and `fg` contain Figure 5b-g PDFs/PNGs and corresponding tables; `panel_status.csv`, `package_versions.csv` and `run_metadata.csv` contain verification and runtime information. A tested run took approximately one minute; runtime varies by computer. Set `DEMO_R_LIB` if packages are in an additional R library.
