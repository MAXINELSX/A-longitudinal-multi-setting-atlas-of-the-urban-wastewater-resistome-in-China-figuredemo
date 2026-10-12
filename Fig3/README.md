# Figure 3 small downstream demos

All five panel labels have a runnable small demonstration. These are software examples, not a full reproduction of the manuscript figure or its upstream analysis. Panel c uses a deliberately simulated 12-sample contig table and summaries computed from that same table. Other demo numerical inputs are retained from the supplied processed examples; they are not validated research observations.

The inputs use generic filenames. Samples, cities and bins have generic identifiers. Known inconsistencies in the legacy inputs are retained for the audit. Private reverse mappings are excluded. Explanations and provenance are confined to these README files; output tables contain numerical results, categorical labels and compact diagnostic fields.

## Included branches

| Panel | Demonstration | Boundary |
| --- | --- | --- |
| a | Historical relative-abundance heatmap from a small selected table; arithmetic means within anonymous city/setting groups | Relative values are not pooled integer counts. This does not reproduce a count-pooled current panel. |
| b | Small fraction-weighted flow display from precomputed setting/context/ARG-type proportions | No raw annotation counts are available in this input; no counts are reconstructed and no selected subset is rescaled to 100%. |
| c | A violin plot and setting summaries from the same 12 simulated samples, three per setting | Counts are ARG-carrying contigs. The four summaries are recomputed and checked against the bundled summary table. No annotation-level mixed model is run. |
| d | Retained burden aggregation, Kruskal-Wallis/Dunn/BH functions and display with synthetic site/date metadata | Demo inference only; source rows do not establish distinct genomes. |
| e | Small display of supplied species/setting/ARG-type cell counts and percentages | Supplied MAG denominators are retained per cell. Full MAG registry/union and unavailable burden inputs are not reconstructed. |

See `panel_a_demo`, `panel_b_demo`, `panel_c_demo` and `panel_e_demo` for their input schemas, selection rules, provenance and checks. Their aliases must not be joined merely because they look similar. Ordinary setting categories are retained; real locations and sample identities are not.

## Retained burden input and legacy audit

The six-table legacy bundle in `input/` remains for audit. `input_audit/panel_readiness.csv` records compact codes for that bundle's limitations: missing species counts, missing classified counts, contig/annotation unit mismatch, missing site/date metadata with duplicate IDs, and missing MAG registry. These concern the original complete workflows, not the runnable display branches.

The burden table contains 216 rows. Forty-five lack species assignments and are excluded from the species-level demonstration. Three sample/bin combinations recur with conflicting values. All rows are retained in the audit, and each row receives a unique surrogate identifier; that identifier does not establish a distinct biological genome. The 171 species-assigned rows enter the supplied 3d aggregation and Kruskal–Wallis/Dunn/BH functions where estimable. Fifteen species are displayed using the supplied frequency ordering.

Panel d site identities are generated from generic city/setting combinations, with successive synthetic monthly dates. The `synthetic_metadata` boolean records this status. The resulting tests and letters demonstrate software behavior and are not scientific evidence. `INCLUDED` and `NO_SPECIES` in the row map identify inclusion in the species analysis. The legacy contig table has 18 samples and is audited only; panel c uses the separate simulated 12-row table.

## Run and verify

```r
install.packages(c('data.table', 'dplyr', 'ggplot2', 'ggalluvial', 'multcompView', 'readr', 'scales'))
```

```sh
Rscript run_demo.R
Rscript run_demo.R --output output_check
Rscript verify_outputs.R reference_output
Rscript test_fig3_audit.R
Rscript test_fig3_demo.R
Rscript test_fig3_clean.R
```

The default output is `reference_output/`. `--input DIR` changes only the legacy six-table bundle used for d and its audit; new display branches use inputs relative to their own submodule directories. Every new panel runner accepts `--output DIR`. Script-relative input resolution supports execution from another working directory. No packages are installed automatically. Use a new output directory to avoid mixing different versions.

The output includes input audits, generic row identifiers, synthetic d metadata and the five panel demonstrations. Audit/adapter tests check this example's known row counts and defects, not arbitrary replacement data. `verify_outputs.R` also invokes each new panel's `test_outputs.R --output DIR` for saved numerical checks. Run times, package versions and session details are exported; timezone identifiers are omitted. Panels a and e export warning class/count diagnostics, with full R warnings remaining visible in the console if a warning occurs. No separate scope, explanation or warning-text files are generated. Plots carry only simple figure titles, axes, legends, categorical labels and statistical results.

## Code adaptations

`upstream/analysis_functions.R` is an adapted copy of previously supplied figure-analysis code. City lookup/display literals and input filenames are generic. The wrapper handles missing/duplicated input records explicitly; statistical aggregation and testing functions are retained. The source copy is not byte-identical to the original. The legacy tables were supplied processed examples; panels a and e retain selected supplied example values, while panel b retains selected precomputed global fractions from the public source identified in its README. Panel c was replaced with a fresh, small simulated example for this version.

New display adapters are documented separately and do not recompute omitted upstream stages. The derived identifier table maps already-anonymous rows to surrogate demo rows, not back to original identities. Identifier removal does not guarantee that retained numerical values cannot be linked to other data. No new code or data licence is conferred by this package.
