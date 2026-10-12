verify_outputs <- function(output_dir) {
  tab <- function(name) read.csv(file.path(output_dir, name), check.names = FALSE, stringsAsFactors = FALSE)
  rf <- tab("rf_importance.csv")
  lmm <- tab("lmm_coefficients.csv")
  diagnostics <- tab("lmm_diagnostics.csv")
  status <- tab("panel_status.csv")
  metadata <- tab("run_metadata.csv")
  meta <- setNames(metadata$Value, metadata$Key)
  a <- tab("plot_data/Fig6a_predictors.csv")
  cl <- tab("plot_data/Fig6a_classes.csv")
  e <- tab("plot_data/Fig6e_display.csv")
  checks <- c(
    rf_model_actually_fitted = grepl("rfPermute", meta[["rf_model_class"]]),
    rf_nonempty_finite_importance = nrow(rf) == 20L && all(is.finite(rf[["%IncMSE"]])),
    rf_actual_permutation_probabilities = all(is.finite(rf[["%IncMSE.pval"]])) && all(rf[["%IncMSE.pval"]] >= 0 & rf[["%IncMSE.pval"]] <= 1),
    rf_positive_signed_total = sum(rf[["%IncMSE"]]) > 0,
    rf_source_signed_percentage_rule = isTRUE(all.equal(a$Relative_unclipped, a$Raw_importance/sum(a$Raw_importance)*100)),
    rf_class_shares_sum_100 = abs(sum(cl$Percent)-100) < 1e-8,
    lmm_four_successful_models = nrow(diagnostics) == 4L && all(diagnostics$Status == "ok"),
    lmm_all_eighteen_exposures_estimated = nrow(lmm) == 18L && all(is.finite(lmm$coef)) && all(is.finite(lmm$p)),
    lmm_source_formula_retained = all(grepl("factor(season)", diagnostics$Formula, fixed = TRUE)) && all(grepl("(1 | city) + (1 | site_nm)", diagnostics$Formula, fixed = TRUE)),
    lmm_only_real_30day_window = all(lmm$lag == "lag1") && all(e$Lag == "Lag30"),
    lmm_nonempty_display = nrow(e) == 8L && all(e$CI_low <= e$Estimate & e$Estimate <= e$CI_high),
    lmm_BH_within_outcome_class = all(vapply(split(lmm, interaction(lmm$outcome, lmm$class)), function(z) isTRUE(all.equal(z$fdr, p.adjust(z$p, "BH"))), logical(1))),
    compact_panel_status = identical(names(status), c("Panel", "Status")),
    unrun_panels_explicit = all(status$Status[status$Panel %in% c("b","c","d","f","g","h")] == "NOT_TESTED"),
    no_unrun_panel_images = !any(file.exists(file.path(output_dir, paste0("Fig6", c("b","c","d","f","g","h"), ".png"))))
  )
  image_files <- file.path(output_dir, c("Fig6a.pdf", "Fig6a.png", "Fig6e.pdf", "Fig6e.png"))
  checks <- c(checks, four_nonempty_plot_files = all(file.exists(image_files)) && all(file.info(image_files)$size > 1000))
  report <- data.frame(Check = names(checks), Result = ifelse(checks, "PASS", "FAIL"))
  write.csv(report, file.path(output_dir, "test_results.csv"), row.names = FALSE)
  if (!all(checks)) stop("Verification failed: ", paste(names(checks)[!checks], collapse = ", "))
  invisible(report)
}

test_script <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(test_script) && basename(sub("^--file=", "", test_script[1])) == "verify_outputs.R") {
  test_args <- commandArgs(TRUE)
  if (length(test_args) != 2L || test_args[1] != "--output") stop("Usage: Rscript tests/verify_outputs.R --output DIR")
  verify_outputs(test_args[2])
  message("PASS: Fig6 a/e")
}
