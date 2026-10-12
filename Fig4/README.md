# Figure 4 demo

This module renders Figure 4c, a four-pair subset for Figure 4d, and the three-component Figure 4e format. Figure 4e summary-count inputs are simulated; the Boolean evidence marks are retained technical fixtures. The simulated counts are not manuscript values or full-study totals and must not be interpreted as biological findings, prevalence estimates, or host or allele rankings. Figure 4c and Figure 4d retain the numerical values from the previous demo. No genomic matching, clustering, read mapping, allele comparison, candidate scoring, or biological inference is performed.

| Input | Rows | Origin and purpose |
| --- | ---: | --- |
| `input/cluster_counts.tsv` | 3 | Existing anonymous cluster-count example for 4c |
| `input/pair_demo.tsv` | 4 | Existing anonymous processed pair subset for 4d |
| `input/e_carbapenemase_overlap.tsv` | 12 | Boolean fixture for 4e evidence marks |
| `input/e_mobile_genetic_context.tsv` | 3 | Simulated context counts |
| `input/taxon_counts.tsv` | 7 | Simulated taxonomy counts |
| `input/e_classification_coverage.tsv` | 2 | Simulated classification coverage |

All explanations and data provenance are recorded here. Plots contain only panel titles, axes, legends, category labels, and numerical labels. Run metadata use machine-readable fields. Original filenames, source paths, location fields, and reverse identifier mappings are excluded. The plotting code is adapted from the previously supplied figure functions; the category format and ordering are retained. No new code or data license is conferred by this package.

## Run and verify

Requirements: R, `ggplot2`, `patchwork`, and Cairo graphics support. The scripts do not install packages. Install dependencies in an R library if needed:

Verified on Windows with R 4.6.0, ggplot2 4.0.3, patchwork 1.3.2, Cairo support, and the C locale. A standalone run into an output path containing spaces completed in 30.210 seconds; elapsed time varies by system. All three generated PNG plots were visually checked.

```r
install.packages(c("ggplot2", "patchwork"))
```

```sh
Rscript --vanilla run_demo.R
Rscript --vanilla run_demo.R --output "output check"
Rscript --vanilla verify_saved_outputs.R --output "output check"
Rscript --vanilla test_demo.R
```

Inputs are found relative to `run_demo.R`, so calling the script by absolute path from another working directory works. Quote paths containing spaces. The default output is `reference_output/`; positional input and output directories are also accepted. `DEMO_R_LIB` can add an existing R package library. Arial is requested by the supplied plotting functions; available fonts can affect appearance on other systems.

`test_demo.R` checks retained 4c/4d values, rejects invalid 4d flags and counts, checks simulated 4e counts and percentages, rejects inconsistent context and coverage totals, checks the absence of plot subtitles and captions, renders into a temporary path containing spaces, and verifies saved results. `verify_saved_outputs.R` checks an existing run without rendering. Seventeen nonempty output files are expected: six PDF/PNG plots, eight data tables, and three run metadata files. No serialized R objects are saved. `runtime.txt` records elapsed seconds; `run_manifest.tsv` records row counts and mode codes; `sessionInfo.txt` records R and package versions while omitting timezone identifiers and local library paths.

## Figure 4c and 4d

Figure 4c uses three existing clusters: clinical counts 1, 1, 3; wastewater counts 1, 1, 1; and wastewater percentages 50, 50, 25. It retains the count-and-percentage plotting format with a simple title and axis label.

Figure 4d retains four rows at previously selected fixed positions across the processed source order. The selection was not based on biological metrics, and the rows are not a representative sample. Pair, clinical, and MAG identifiers use consistent aliases. The original input numeric strings and Boolean flags are unchanged.

The partial 4d plot uses `fixed_snp_per_mbp_pair`, `callable_fraction_pair`, `clinical_ARG_count`, `mag_ARG_count`, and `shared_ARG_count`; `ANI_min` remains in the exported table. The supplied pair-level callable fraction is retained. The three supplied flags are `direct_remap_supported_bool`, `instrain_partial_support_bool`, and `shared_carbapenemase_present_bool`. No specimen information or high-priority ARG annotation is inferred. The shared-carbapenemase Boolean retains its own meaning.

`d_partial/Fig4d_displayed_pairs.tsv` has four rows and 13 columns: three anonymous identifiers, `setting`, six numeric fields, and three Boolean fields. Figure 4a/b and the complete manuscript 4d are outside this example.

## Figure 4e simulated values

The 4e example represents 120 simulated contigs. Genetic context counts are 48 plasmid associated, 60 chromosome associated, and 12 phage associated, giving 40%, 50%, and 10%. Classification coverage is 80 classified and 40 unclassified, giving exact percentages of 66.666...% and 33.333...%, stored in the input as 66.7% and 33.3%. Context and coverage therefore use the same total of 120.

The seven taxonomy categories sum to the classified count of 80: `Other identified taxa` contains 20 and each of the six named categories contains 10. Their percentages are 25% and 12.5%, respectively. These arbitrary demonstration counts exercise the existing categories and plotting format without asserting a biological distribution. The category ordering is fixed and is not a ranking derived from these counts.

| Summary | Denominator | Export |
| --- | ---: | --- |
| Genetic context | 120 | `e/Fig4e_context_recalculated.tsv` |
| Classified taxonomy | 80 | `e/Fig4e_taxa_recalculated.tsv` |
| Classification coverage | 120 | `e/Fig4e_classification_coverage.tsv` |

The 12-row Boolean evidence fixture retains the previous Yes/No values as a technical example: each of three allele labels has four evidence categories. The existing renderer displays KPC-2 and OXA-48, giving eight displayed marks; NDM-5 remains in the full exported fixture. No evidence is reassessed, and these entries do not establish actual presence, sharing, or host assignment. Coverage is tabulated rather than plotted as a fourth component.

`e/Fig4e_evidence_all.tsv` and `e/Fig4e_evidence_displayed.tsv` contain `Evidence`, `Allele`, `Supported`, and `Support_score`. Category exports preserve input columns and add `Percent_recalculated`; the recalculation uses each summary's count total. `e/Fig4e_denominators.tsv` records the three denominators. Validation requires the context total to equal the coverage total and the taxonomy total to equal the classified coverage count. Bar limits and label positions scale to the small simulated counts.
